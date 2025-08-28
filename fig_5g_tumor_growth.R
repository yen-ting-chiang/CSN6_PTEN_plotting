# 可選：設定工作資料夾
setwd("C:/Users/danny/Documents/R_project/CSN6_PTEN_plotting/fig_5g_tumor_growth")

# =========================================================
# Nature journal style line plots for Tumor Growth assay
# Day 0 一律設為 0 並繪製；其餘天照 input
# 輸出 SD/SEM、逐日 Welch (Holm) 檢定 CSV、圖檔 .tiff/.pdf/.jpg
# =========================================================

# ---- 0) 基本設定 ----
input_path      <- "fig_5g_tumor_growth_cleaned.xlsx"
sheet_name      <- 1
figure_tag      <- "Fig_5G_TumorGrowth_line"
out_dir         <- "figure_outputs"

fig_width_in  <- 3.35
fig_height_in <- 3.0

# 造型參數（與 MTT 版一致）
errorbar_color_mode <- "group"   # "group" 或 "black"
errorbar_linewidth  <- 0.10
mean_point_size     <- 2.8
mean_point_stroke   <- 0.30
y_label <- "Tumor growth (a.u.)"

# ---- 1) 套件 ----
pkgs <- c("readxl","tidyverse","rstatix","ggpubr")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 2) 讀取與轉長表 ----
stopifnot(file.exists(input_path))
raw <- readxl::read_excel(path = input_path, sheet = sheet_name)

nm    <- names(raw); nm_lc <- tolower(nm)
day_cols <- nm[grepl("^\\s*day\\s*\\d+\\s*$", nm_lc)]
if (length(day_cols) == 0) stop("找不到任何 'day 0/1/2...' 欄位。")

# 嘗試抓 Group；同時加上每列的樣本 ID（.RID）以便補 D0
group_candidates <- nm[nm_lc %in% c("group","condition","treatment","genotype","line","cellline")]
raw <- raw %>% dplyr::mutate(.RID = dplyr::row_number())
if (length(group_candidates) >= 1) {
  group_col <- group_candidates[1]
  long_df <- raw %>%
    dplyr::select(all_of(c(group_col, ".RID", day_cols))) %>%
    tidyr::pivot_longer(cols = all_of(day_cols),
                        names_to = "Day_raw", values_to = "Value") %>%
    dplyr::rename(Group = !!group_col)
} else {
  long_df <- raw %>%
    dplyr::select(all_of(c(".RID", day_cols))) %>%
    tidyr::pivot_longer(cols = all_of(day_cols),
                        names_to = "Day_raw", values_to = "Value") %>%
    dplyr::mutate(Group = "Group")
}

# 轉成數字日、標籤
long_df <- long_df %>%
  dplyr::mutate(
    Day_num   = readr::parse_number(tolower(Day_raw)),
    Day_label = paste0("D", Day_num),
    Group     = as.character(Group)
  )

# ====== 關鍵修正：確保「每個樣本的 Day 0 都存在且值=0」 ======
has_d0 <- any(long_df$Day_num == 0)
if (!has_d0) {
  day0_tbl <- long_df %>%
    dplyr::distinct(Group, .RID) %>%
    dplyr::mutate(Day_raw = "day 0", Day_num = 0, Day_label = "D0", Value = 0)
  long_df <- dplyr::bind_rows(day0_tbl, long_df)
} else {
  long_df <- long_df %>% dplyr::mutate(Value = ifelse(Day_num == 0, 0, Value))
}
# 去掉非 D0 的 NA；D0 一律保留（即使原本 NA）
long_df <- long_df %>% dplyr::filter(Day_num == 0 | !is.na(Value))

# Day 順序（含 D0）
day_levels <- long_df %>%
  dplyr::distinct(Day_label, Day_num) %>%
  dplyr::arrange(Day_num) %>% dplyr::pull(Day_label)
long_df <- long_df %>% dplyr::mutate(Day_label = factor(Day_label, levels = day_levels))

# ---- 3) 固定 Group 顯示順序 ----
desired_order <- c("CMV5","PTEN","PTEN+CSN6","CSN6")
long_df <- long_df %>%
  dplyr::mutate(
    Group = dplyr::recode(
      Group,
      "CMV"         = "CMV5",
      "PTEN + CSN6" = "PTEN+CSN6",
      "CSN6 + PTEN" = "PTEN+CSN6",
      .default = Group
    )
  )
levs <- desired_order[desired_order %in% unique(long_df$Group)]
long_df <- long_df %>% dplyr::mutate(Group = factor(Group, levels = levs))

# ---- 4) 色盤 ----
pal_groups_named <- c(
  "CMV5"      = "#E41A1C",
  "PTEN"      = "#FFB300",
  "PTEN+CSN6" = "#A6D854",
  "CSN6"      = "#4DBBD5"
)
pal_groups <- pal_groups_named[levels(long_df$Group)]

# ---- 5) 摘要與檢定 ----
summary_by_group_day <- function(df) {
  df %>%
    dplyr::group_by(Group, Day_label, Day_num) %>%
    dplyr::summarise(
      n    = dplyr::n(),
      mean = mean(Value, na.rm = TRUE),
      sd   = sd(Value, na.rm = TRUE),
      sem  = sd / sqrt(n),
      .groups = "drop"
    )
}

pairwise_tests_by_day <- function(df, mode_tag, out_dir_path, figure_tag_prefix) {
  days <- df %>% dplyr::distinct(Day_label) %>% dplyr::pull(Day_label)
  res_list <- lapply(as.character(days), function(dlv) {
    subdf <- df %>% dplyr::filter(Day_label == dlv)
    ok <- subdf %>% dplyr::count(Group) %>% dplyr::pull(n)
    if (dplyr::n_distinct(subdf$Group) < 2 || any(ok < 2)) return(NULL)
    tryCatch(
      rstatix::pairwise_t_test(
        data = subdf, formula = Value ~ Group,
        p.adjust.method = "holm", pool.sd = FALSE,
        paired = FALSE, detailed = TRUE
      ) %>% dplyr::mutate(Day_label = dlv),
      error = function(e) NULL
    )
  })
  res <- dplyr::bind_rows(res_list)
  if (nrow(res) > 0) {
    export <- res %>%
      dplyr::transmute(
        Day_label,
        comparison  = paste(group1, "vs", group2),
        group1, group2, n1, n2, statistic, df,
        p_value     = p,
        p_adj_holm  = p.adj,
        signif_adj  = p.adj.signif,
        conf_low    = conf.low, conf_high = conf.high,
        method
      )
    if (!dir.exists(out_dir_path)) dir.create(out_dir_path, recursive = TRUE)
    readr::write_csv(
      export,
      file.path(out_dir_path, paste0(figure_tag_prefix, "_", mode_tag, "_pairwise_ttest_pvalues.csv"))
    )
  }
  invisible(res)
}

# ---- 6) 存檔 ----
save_all_formats <- function(plot, tag, w = fig_width_in, h = fig_height_in, dir = out_dir){
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  ggsave(file.path(dir, paste0(tag, ".pdf")),
         plot = plot, width = w, height = h, units = "in", device = cairo_pdf)
  ggsave(file.path(dir, paste0(tag, ".tiff")),
         plot = plot, width = w, height = h, units = "in",
         dpi = 600, device = "tiff", compression = "lzw")
  ggsave(file.path(dir, paste0(tag, ".jpg")),
         plot = plot, width = w, height = h, units = "in",
         dpi = 600, device = "jpeg")
}

# ---- 7) 作圖（SD/SEM；誤差棒可用 group/black）----
make_lineplot <- function(sum_df, err = c("SD","SEM"),
                          y_lab = y_label,
                          palette = pal_groups,
                          tag_suffix = "D0is0_SD_or_SEM",
                          errorbar_color = errorbar_color_mode,
                          errbar_lwd = errorbar_linewidth,
                          mean_pt_size = mean_point_size,
                          mean_pt_stroke = mean_point_stroke) {
  err <- match.arg(err)
  err_col <- if (err == "SD") "sd" else "sem"
  dfp <- sum_df %>% dplyr::mutate(err = .data[[err_col]])
  
  grp_lv <- if (is.factor(dfp$Group)) levels(dfp$Group) else unique(dfp$Group)
  base_shapes <- c(15,18,17,16,8,7,3,4,23,24)
  shape_map <- setNames(rep(base_shapes, length.out = length(grp_lv)), grp_lv)
  
  ymax <- max(dfp$mean + dfp$err, na.rm = TRUE)
  base_top <- ifelse(is.finite(ymax) && ymax > 0, ymax, 1)
  ylim_top <- base_top * 1.12
  
  p <- ggplot(dfp, aes(x = Day_label, y = mean, group = Group, color = Group)) +
    geom_line(linewidth = 0.6)
  
  if (identical(errorbar_color, "black")) {
    p <- p + geom_errorbar(aes(ymin = mean - err, ymax = mean + err),
                           width = 0.15, linewidth = errbar_lwd, color = "black")
  } else {
    p <- p + geom_errorbar(aes(ymin = mean - err, ymax = mean + err),
                           width = 0.15, linewidth = errbar_lwd)
  }
  
  p <- p +
    geom_point(aes(shape = Group), size = mean_pt_size, stroke = mean_pt_stroke) +
    scale_color_manual(values = palette) +
    scale_shape_manual(values = shape_map) +
    guides(color = guide_legend(title = NULL, nrow = 1),
           shape = guide_legend(title = NULL, nrow = 1)) +
    labs(x = NULL, y = y_lab) +
    coord_cartesian(ylim = c(0, ylim_top)) +
    theme_classic(base_size = 8, base_family = "sans") +
    theme(
      axis.text.x  = element_text(size = 8, color = "black"),
      axis.text.y  = element_text(size = 8, color = "black"),
      axis.title.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
      axis.line    = element_line(linewidth = 0.4, color = "black"),
      axis.ticks   = element_line(linewidth = 0.4, color = "black"),
      legend.position = "top",
      legend.key.height = grid::unit(0.3, "lines"),
      plot.margin  = margin(t = 2, r = 4, b = 2, l = 4)
    )
  attr(p, "tag_suffix") <- tag_suffix
  p
}

# ---- 8) 摘要 → 畫圖 → 存檔 → 檢定 ----
sum_df <- summary_by_group_day(long_df)

p_SD  <- make_lineplot(sum_df, err = "SD",  tag_suffix = "D0is0_SD")
p_SEM <- make_lineplot(sum_df, err = "SEM", tag_suffix = "D0is0_SEM")

save_all_formats(p_SD,  paste0(figure_tag, "_", attr(p_SD,  "tag_suffix")))
save_all_formats(p_SEM, paste0(figure_tag, "_", attr(p_SEM, "tag_suffix")))

pairwise_tests_by_day(long_df, mode_tag = "D0is0",
                      out_dir_path = out_dir, figure_tag_prefix = figure_tag)

message("Done. Files saved to: ", normalizePath(out_dir))
print(p_SD); print(p_SEM)
