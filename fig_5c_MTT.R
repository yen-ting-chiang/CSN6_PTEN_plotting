# 可選：設定工作資料夾
setwd("C:/Users/danny/Documents/R_project/CSN6_PTEN_plotting/fig_5c_MTT")

# =========================================================
# Nature journal style line plots for MTT proliferation assay
# - 讀取 fig_5c_MTT_cleaned.xlsx（欄：day 0 ~ day 6，另有 Group/Condition/Treatment 等）
# - 4 種繪圖方式 × (SD/SEM) 誤差棒
# - Welch's t-test（逐日、Holm 校正）輸出 CSV
# - 輸出：.tiff (LZW, 600 dpi), .pdf (向量), .jpg (600 dpi)
# =========================================================

# ---- 0) 基本設定 ----
input_path      <- "fig_5c_MTT_cleaned.xlsx"  # 依需要改路徑
sheet_name      <- 1
figure_tag      <- "Fig_5C_MTT_line"
out_dir         <- "figure_outputs"
norm_as_percent <- FALSE    # TRUE: baseline=100%；FALSE: baseline=1（倍數）

# 尺寸：Nature 單欄寬 ~3.35 in；高度加高
fig_width_in  <- 3.35
fig_height_in <- 3.0   # 加高

# 造型參數（可調）
errorbar_color_mode <- "group"  # "group" 或 "black"；預設各組顏色
errorbar_linewidth  <- 0.10     # 誤差棒線寬（更細）
mean_point_size     <- 2.8      # 平均值點大小（更大）
mean_point_stroke   <- 0.30     # 平均值點外框

# ---- 1) 套件安裝 / 載入 ----
pkgs <- c("readxl","tidyverse","rstatix","ggpubr")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 2) 讀取 + 依 day 0 ~ day 6 欄位轉長表 ----
stopifnot(file.exists(input_path))
raw <- readxl::read_excel(path = input_path, sheet = sheet_name)

nm    <- names(raw)
nm_lc <- tolower(nm)
day_cols <- nm[grepl("^\\s*day\\s*\\d+\\s*$", nm_lc)]
if (length(day_cols) == 0) stop("找不到任何 'day 0' ~ 'day N' 這類欄位，請確認欄名。")

# 嘗試找 Group 欄（常見命名），若沒有就把所有列視為同一組
group_candidates <- nm[nm_lc %in% c("group","condition","treatment","genotype","line","cellline")]
if (length(group_candidates) >= 1) {
  group_col <- group_candidates[1]
  long_df <- raw %>%
    dplyr::select(all_of(c(group_col, day_cols))) %>%
    tidyr::pivot_longer(cols = all_of(day_cols),
                        names_to = "Day_raw", values_to = "Value") %>%
    dplyr::rename(Group = !!group_col)
} else {
  long_df <- raw %>%
    dplyr::select(all_of(day_cols)) %>%
    tidyr::pivot_longer(cols = everything(),
                        names_to = "Day_raw", values_to = "Value") %>%
    dplyr::mutate(Group = "Group")
}

# "day 0" -> Day_num=0；x 軸標籤用 D0, D1...
long_df <- long_df %>%
  dplyr::mutate(
    Day_num   = readr::parse_number(tolower(Day_raw)),
    Day_label = paste0("D", Day_num),
    Group     = as.character(Group)
  ) %>%
  dplyr::filter(!is.na(Value))

day_levels <- long_df %>%
  dplyr::distinct(Day_label, Day_num) %>%
  dplyr::arrange(Day_num) %>%
  dplyr::pull(Day_label)
long_df <- long_df %>%
  dplyr::mutate(Day_label = factor(Day_label, levels = day_levels))

# ---- 3.0) 固定 Group 名稱與圖例順序（CMV5 → PTEN → PTEN+CSN6 → CSN6）----
desired_order <- c("CMV5","PTEN","PTEN+CSN6","CSN6")

# 名稱正規化（把常見寫法歸一）
long_df <- long_df %>%
  dplyr::mutate(
    Group = dplyr::recode(
      Group,
      "CMV"          = "CMV5",
      "PTEN + CSN6"  = "PTEN+CSN6",
      "CSN6 + PTEN"  = "PTEN+CSN6",
      .default = Group
    )
  )

# 只保留資料中出現的群組，但依 desired_order 固定順序
levs <- desired_order[desired_order %in% unique(long_df$Group)]
long_df <- long_df %>% dplyr::mutate(Group = factor(Group, levels = levs))

# ---- 3) 色盤（依固定四組命名；會依 factor levels 排序）----
pal_groups_named <- c(
  "CMV5"      = "#E41A1C",  # vivid red
  "PTEN"      = "#FFB300",  # vivid orange
  "PTEN+CSN6" = "#A6D854",  # light green
  "CSN6"      = "#4DBBD5"   # sky blue
)
pal_groups <- pal_groups_named[levels(long_df$Group)]

# ---- 4) 摘要與統計函式 ----
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

# 顯式傳入 out_dir/figure_tag，避免 promise under evaluation 問題
pairwise_tests_by_day <- function(df, mode_tag,
                                  out_dir_path,
                                  figure_tag_prefix) {
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

# ---- 5) Normalization 函式 ----
normalize_to_baseline <- function(df, baseline_day_num, to_percent = FALSE) {
  base_tbl <- df %>%
    dplyr::filter(Day_num == baseline_day_num) %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(baseline = mean(Value, na.rm = TRUE), .groups = "drop")
  df2 <- df %>%
    dplyr::left_join(base_tbl, by = "Group") %>%
    dplyr::mutate(
      Value = Value / baseline,
      Value = if (to_percent) Value * 100 else Value
    ) %>%
    dplyr::filter(is.finite(Value))
  df2
}

# ---- 6) 輸出多格式 ----
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

# ---- 6.5) 作圖函式（可選 SD/SEM；誤差棒顏色可 group/black；畫平均值點，不畫原始點）----
make_lineplot <- function(sum_df, raw_df, err = c("SD","SEM"),
                          y_lab = "Absorbance (a.u.)",
                          palette = pal_groups,
                          tag_suffix = "D0_raw",
                          errorbar_color = errorbar_color_mode,   # "group" 或 "black"
                          errbar_lwd = errorbar_linewidth,
                          mean_pt_size = mean_point_size,
                          mean_pt_stroke = mean_point_stroke) {
  
  err <- match.arg(err)
  err_col <- if (err == "SD") "sd" else "sem"
  dfp <- sum_df %>% dplyr::mutate(err = .data[[err_col]])
  
  # 平均值點 shape 對應各組（依目前 factor levels）
  grp_lv <- if (is.factor(dfp$Group)) levels(dfp$Group) else unique(dfp$Group)
  base_shapes <- c(15, 18, 17, 16, 8, 7, 3, 4, 23, 24)  # 方/菱/三角/圓/星..高辨識序列
  shape_map <- setNames(rep(base_shapes, length.out = length(grp_lv)), grp_lv)
  
  # y 軸上界
  ymax <- max(dfp$mean + dfp$err, na.rm = TRUE)
  base_top <- ifelse(is.finite(ymax) && ymax > 0, ymax, 1)
  ylim_top <- base_top * 1.12
  
  p <- ggplot(dfp, aes(x = Day_label, y = mean, group = Group, color = Group)) +
    geom_line(linewidth = 0.6)
  
  # 誤差棒：依設定決定顏色
  if (identical(errorbar_color, "black")) {
    p <- p + geom_errorbar(aes(ymin = mean - err, ymax = mean + err),
                           width = 0.15, linewidth = errbar_lwd, color = "black")
  } else {
    p <- p + geom_errorbar(aes(ymin = mean - err, ymax = mean + err),
                           width = 0.15, linewidth = errbar_lwd)  # 承接 color=Group
  }
  
  # 各組「平均值」點（不畫 raw data）
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
  return(p)
}

save_both_errplots <- function(long_cur, mode_tag, y_lab) {
  sum_cur <- summary_by_group_day(long_cur)
  p_sd  <- make_lineplot(sum_cur, long_cur, err = "SD",  y_lab = y_lab,
                         tag_suffix = paste0(mode_tag, "_SD"))
  p_sem <- make_lineplot(sum_cur, long_cur, err = "SEM", y_lab = y_lab,
                         tag_suffix = paste0(mode_tag, "_SEM"))
  
  save_all_formats(p_sd,  paste0(figure_tag, "_", attr(p_sd,  "tag_suffix")))
  save_all_formats(p_sem, paste0(figure_tag, "_", attr(p_sem, "tag_suffix")))
  
  pairwise_tests_by_day(
    long_cur,
    mode_tag = mode_tag,
    out_dir_path = out_dir,
    figure_tag_prefix = figure_tag
  )
  
  print(p_sd); print(p_sem)
}

# ---- 7) 四種繪圖方式 ----
# A) 從 D0 開始，不做 normalization
long_A <- long_df %>% dplyr::filter(!is.na(Day_num), Day_num >= 0)
ylab_A <- "Absorbance (a.u.)"
save_both_errplots(long_A, mode_tag = "D0_raw", y_lab = ylab_A)

# B) 從 D0 開始，baseline = D0（1 或 100%）
long_B <- long_df %>% dplyr::filter(!is.na(Day_num), Day_num >= 0)
long_B <- normalize_to_baseline(long_B, baseline_day_num = 0, to_percent = norm_as_percent)
ylab_B <- if (norm_as_percent) "Relative to D0 (%)" else "Relative to D0 (×)"
save_both_errplots(long_B, mode_tag = "D0_norm", y_lab = ylab_B)

# C) 從 D1 開始，不做 normalization
long_C <- long_df %>% dplyr::filter(!is.na(Day_num), Day_num >= 1)
ylab_C <- "Absorbance (a.u.)"
save_both_errplots(long_C, mode_tag = "D1_raw", y_lab = ylab_C)

# D) 從 D1 開始，baseline = D1（1 或 100%）
long_D <- long_df %>% dplyr::filter(!is.na(Day_num), Day_num >= 1)
long_D <- normalize_to_baseline(long_D, baseline_day_num = 1, to_percent = norm_as_percent)
ylab_D <- if (norm_as_percent) "Relative to D1 (%)" else "Relative to D1 (×)"
save_both_errplots(long_D, mode_tag = "D1_norm", y_lab = ylab_D)

message("Done. Files saved to: ", normalizePath(out_dir))
