# 可選：設定工作資料夾
setwd("C:/Users/danny/Documents/R_project/CSN6_PTEN_plotting/fig_5g_tumor_growth_Jin_figure_data")

# =========================================================
# Nature journal style line plots for Tumor Growth assay
# (Jin_figure_data：表內已給各組各日 mean 與 SD)
# - 不設定/不繪製 Day 0；其他天照 input
# - 無個別 raw 點；僅畫各組每日平均與 SD 誤差棒
# - 輸出：.tiff (LZW, 600 dpi), .pdf (向量), .jpg (600 dpi)
# =========================================================

# ---- 0) 基本設定 ----
input_path <- "fig_5g_tumor_growth_Jin_figure_data_cleaned.xlsx"
sheet_name <- 1
figure_tag <- "Fig_5G_TumorGrowth_JinFig_line"
out_dir    <- "figure_outputs"

# 尺寸與造型（與先前一致）
fig_width_in        <- 3.35
fig_height_in       <- 3.0
errorbar_color_mode <- "group"   # "group" 或 "black"
errorbar_linewidth  <- 0.10
mean_point_size     <- 2.8
mean_point_stroke   <- 0.30
y_label             <- "Tumor volume (a.u.)"

# 固定圖例順序與配色
desired_order <- c("CMV5","PTEN","PTEN+CSN6","CSN6")
pal_groups_named <- c(
  "CMV5"      = "#E41A1C",
  "PTEN"      = "#FFB300",
  "PTEN+CSN6" = "#A6D854",
  "CSN6"      = "#4DBBD5"
)

# ---- 1) 套件 ----
pkgs <- c("readxl","tidyverse","ggpubr")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 2) 讀取 & 轉為「每組-每日-mean/sd」的長表 ----
stopifnot(file.exists(input_path))
df_raw <- readxl::read_excel(path = input_path, sheet = sheet_name)

nm <- names(df_raw); nm_lc <- tolower(nm)
group_candidates <- nm[nm_lc %in% c("group","condition","treatment","genotype","line","cellline")]
group_col <- if (length(group_candidates) >= 1) group_candidates[1] else NA_character_

# 檢測「雙列表頭」：第 1 列為 day labels，欄名含 "tumor" 與 "sd"
special_two_row <- (tolower(as.character(df_raw[[1]][1])) %in% c("condition","group")) &&
  any(grepl("tumor", nm_lc)) &&
  any(nm_lc == "sd")

if (special_two_row) {
  idx_mean_start <- which(grepl("tumor", nm_lc))[1]
  idx_sd_start   <- which(nm_lc == "sd")[1]
  mean_idx <- idx_mean_start:(idx_sd_start - 1)
  sd_idx   <- idx_sd_start:ncol(df_raw)
  
  mean_day_labels <- as.character(dplyr::slice(df_raw, 1)[, mean_idx][1, ])
  sd_day_labels   <- as.character(dplyr::slice(df_raw, 1)[, sd_idx][1, ])
  mean_day_num <- readr::parse_number(tolower(mean_day_labels))
  sd_day_num   <- readr::parse_number(tolower(sd_day_labels))
  
  df_me <- df_raw[-1, c(1, mean_idx)]
  names(df_me) <- c("Group", paste0("day_", mean_day_num))
  df_sd <- df_raw[-1, c(1, sd_idx)]
  names(df_sd) <- c("Group", paste0("day_", sd_day_num))
  
  long_me <- df_me %>%
    tidyr::pivot_longer(-Group, names_to = "Day_col", values_to = "mean") %>%
    dplyr::mutate(Day_num = readr::parse_number(Day_col)) %>%
    dplyr::select(Group, Day_num, mean)
  
  long_sd <- df_sd %>%
    tidyr::pivot_longer(-Group, names_to = "Day_col", values_to = "sd") %>%
    dplyr::mutate(Day_num = readr::parse_number(Day_col)) %>%
    dplyr::select(Group, Day_num, sd)
  
  long_sum <- long_me %>%
    dplyr::left_join(long_sd, by = c("Group","Day_num")) %>%
    dplyr::mutate(
      Group = as.character(Group),
      mean  = suppressWarnings(as.numeric(mean)),
      sd    = suppressWarnings(as.numeric(sd))
    )
} else {
  # 一般情況：找「日」的平均欄 (不含 "sd") 與對應的 SD 欄
  day_mean_cols <- nm[grepl("^\\s*day\\s*\\d+\\s*$", nm_lc) |
                        (grepl("^\\s*day\\s*\\d+", nm_lc) & !grepl("sd", nm_lc))]
  if (length(day_mean_cols) == 0)
    stop("偵測不到日次的平均欄（如 'day 1'）。請檢查檔案。")
  
  day_sd_cols <- nm[grepl("^\\s*day\\s*\\d+.*sd", nm_lc) | grepl("sd.*day\\s*\\d+", nm_lc)]
  if (length(day_sd_cols) == 0)
    stop("偵測不到對應的 SD 欄位（如 'day 1 sd' 或 'sd day 1'）。")
  
  # 取 Group 欄
  if (!is.na(group_col)) {
    df_me <- df_raw %>% dplyr::select(all_of(c(group_col, day_mean_cols)))
    df_sd <- df_raw %>% dplyr::select(all_of(c(group_col, day_sd_cols)))
  } else {
    df_me <- df_raw %>% dplyr::mutate(Group = "Group") %>% dplyr::select(Group, all_of(day_mean_cols))
    df_sd <- df_raw %>% dplyr::mutate(Group = "Group") %>% dplyr::select(Group, all_of(day_sd_cols))
  }
  
  long_me <- df_me %>%
    tidyr::pivot_longer(-Group, names_to = "Day_raw", values_to = "mean") %>%
    dplyr::mutate(Day_num = readr::parse_number(tolower(Day_raw))) %>%
    dplyr::select(Group, Day_num, mean)
  
  long_sd <- df_sd %>%
    tidyr::pivot_longer(-Group, names_to = "Day_raw_sd", values_to = "sd") %>%
    dplyr::mutate(Day_num = readr::parse_number(tolower(Day_raw_sd))) %>%
    dplyr::select(Group, Day_num, sd)
  
  long_sum <- long_me %>%
    dplyr::left_join(long_sd, by = c("Group","Day_num")) %>%
    dplyr::mutate(
      Group = as.character(Group),
      mean  = suppressWarnings(as.numeric(mean)),
      sd    = suppressWarnings(as.numeric(sd))
    )
}

# ---- 2.1) 不繪製 Day 0；建立 Day_label 因子順序 ----
long_sum <- long_sum %>%
  dplyr::filter(!is.na(Day_num) & Day_num > 0) %>%   # ← 直接排除 0
  dplyr::mutate(Day_label = paste0("D", Day_num))

day_levels <- long_sum %>%
  dplyr::distinct(Day_label, Day_num) %>%
  dplyr::arrange(Day_num) %>% dplyr::pull(Day_label)
long_sum <- long_sum %>%
  dplyr::mutate(Day_label = factor(Day_label, levels = day_levels))

# ---- 3) 固定 Group 顯示順序與配色 ----
long_sum <- long_sum %>%
  dplyr::mutate(
    Group = dplyr::recode(
      Group,
      "CMV"         = "CMV5",
      "PTEN + CSN6" = "PTEN+CSN6",
      "CSN6 + PTEN" = "PTEN+CSN6",
      .default = Group
    )
  )
levs <- desired_order[desired_order %in% unique(long_sum$Group)]
long_sum <- long_sum %>% dplyr::mutate(Group = factor(Group, levels = levs))
pal_groups <- pal_groups_named[levels(long_sum$Group)]

# ---- 4) 存檔工具 ----
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

# ---- 5) 作圖（僅 SD；不畫 raw 點）----
make_lineplot_sum <- function(sum_df,
                              y_lab = y_label,
                              palette = pal_groups,
                              tag_suffix = "SD",
                              errorbar_color = errorbar_color_mode,
                              errbar_lwd = errorbar_linewidth,
                              mean_pt_size = mean_point_size,
                              mean_pt_stroke = mean_point_stroke) {
  
  dfp <- sum_df %>% dplyr::mutate(err = sd)
  
  grp_lv <- if (is.factor(dfp$Group)) levels(dfp$Group) else unique(dfp$Group)
  base_shapes <- c(15,18,17,16,8,7,3,4,23,24)   # 平均值點樣式（不畫 raw）
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
                           width = 0.15, linewidth = errbar_lwd)  # 承接 color=Group
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

# ---- 6) 產生圖並輸出（僅 SD）----
p_SD <- make_lineplot_sum(long_sum, tag_suffix = "NoD0_SD")
save_all_formats(p_SD, paste0(figure_tag, "_", attr(p_SD, "tag_suffix")))
print(p_SD)

message("Done. Files saved to: ", normalizePath(out_dir))
