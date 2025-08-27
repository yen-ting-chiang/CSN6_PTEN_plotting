setwd("C:/Users/danny/Documents/R_project/CSN6_PTEN_plotting/fig_5e_wound_healing")

# =========================================================
# Nature journal style barplot for wound-healing assay
# - 讀取 /mnt/data/fig_5e_wound_healing_cleaned.xlsx（欄：CMV, PTEN, PTEN+CSN6）
# - 產生兩種圖：SD 誤差棒 & SEM 誤差棒
# - 疊加原始點、Welch's t-test（Holm 校正）顯著性括號與標籤
# - 輸出：.tiff (LZW, 600 dpi), .pdf (向量), .jpg (600 dpi)
# - 另輸出 pairwise t 檢定之 raw p-value 與 Holm adjusted p-value 的 CSV
# =========================================================

# ---- 0) 基本設定 ----
input_path  <- "fig_5e_wound_healing_cleaned.xlsx"
sheet_name  <- 1
y_label     <- "Wound-healing metric (a.u.)"
figure_tag  <- "Fig_5E_WoundHealing_barplot"   # 會自動加上 _SD / _SEM
out_dir     <- "figure_outputs"

# Nature 單欄寬 ≈ 85 mm（3.35 in）
fig_width_in  <- 3.35
fig_height_in <- 2.4

# ---- 1) 套件安裝 / 載入 ----
pkgs <- c("readxl","tidyverse","rstatix","ggpubr")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 2) 讀取資料並整理長格式 ----
stopifnot(file.exists(input_path))
raw_df <- readxl::read_excel(path = input_path, sheet = sheet_name)

expected_cols <- c("CMV","PTEN","PTEN+CSN6")
if (!all(expected_cols %in% names(raw_df))) {
  stop("資料欄位需包含：", paste(expected_cols, collapse = ", "))
}

long_df <- raw_df %>%
  dplyr::select(all_of(expected_cols)) %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = "Group", values_to = "Value") %>%
  dplyr::mutate(Group = factor(Group, levels = expected_cols)) %>%
  dplyr::filter(!is.na(Value))

# ---- 3) 摘要統計（Mean, SD, SEM）----
summary_df <- long_df %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    n    = dplyr::n(),
    mean = mean(Value),
    sd   = sd(Value),
    sem  = sd / sqrt(n),
    .groups = "drop"
  )


# ---- 4) 成對比較（Welch’s t-test, Holm 校正）----
comparisons <- list(c("CMV","PTEN"),
                    c("PTEN","PTEN+CSN6"),
                    c("CMV","PTEN+CSN6"))

ttest_raw <- rstatix::pairwise_t_test(
  data = long_df,
  formula = Value ~ Group,
  comparisons = comparisons,
  p.adjust.method = "holm",
  pool.sd = FALSE,       # <-- 使用 Welch（不假設等變異）
  paired = FALSE,
  detailed = TRUE
)

# 匯出 CSV（raw p 與 Holm adjusted p）
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ttest_export <- ttest_raw %>%
  dplyr::transmute(
    comparison = paste(group1, "vs", group2),
    group1, group2, n1, n2, statistic, df,
    p_value = p,
    p_adj_holm = p.adj,
    signif_adj = p.adj.signif,
    conf_low = conf.low, conf_high = conf.high
  )
readr::write_csv(ttest_export,
                 file.path(out_dir, paste0(figure_tag, "_pairwise_ttest_pvalues.csv")))

# 顯示在圖上的標籤（使用 adjusted p）
ttest_for_plot <- ttest_raw %>%
  dplyr::mutate(
    label = rstatix::p_format(p.adj, digits = 2, accuracy = 0.0001, add.p = TRUE)
  )




# ---- 5) 配色（色盲友善）----
cols <- c("CMV" = "#4C78A8", "PTEN" = "#F58518", "PTEN+CSN6" = "#54A24B")

# ---- 6) 輔助函式：存檔 ----
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


# ---- 7) 作圖函式（可選 SD 或 SEM；可選是否顯示 p 值）----
make_barplot <- function(summary_df, long_df, ttest_df,
                         err_type = c("SD","SEM"),
                         show_p   = FALSE) {
  err_type <- match.arg(err_type)
  err_col  <- if (err_type == "SD") "sd" else "sem"
  df_plot  <- summary_df %>% dplyr::mutate(err = .data[[err_col]])
  
  # y 上界：若不畫 p 值，留一點空白；若畫，依註記行數自動墊高
  y_max    <- max(df_plot$mean + dplyr::coalesce(df_plot$err, 0), na.rm = TRUE)
  base_top <- ifelse(is.finite(y_max) && y_max > 0, y_max, 1)
  step     <- 0.08 * base_top
  
  ttest_df2 <- ttest_df %>%
    dplyr::mutate(y.position = base_top + (dplyr::row_number() - 1) * step)
  
  ylim_top <- if (show_p && nrow(ttest_df2) > 0) {
    base_top + step * (nrow(ttest_df2) + 0.8)
  } else {
    base_top * 1.10
  }
  
  p <- ggplot(df_plot, aes(x = Group, y = mean, fill = Group)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.4) +
    geom_errorbar(aes(ymin = mean - err, ymax = mean + err),
                  width = 0.16, linewidth = 0.4) +
    geom_jitter(data = long_df, aes(x = Group, y = Value),
                width = 0.10, height = 0, size = 1.6, stroke = 0.2,
                shape = 21, color = "black", fill = "white", alpha = 0.9) +
    scale_fill_manual(values = cols, guide = "none") +
    labs(x = NULL, y = y_label) +
    coord_cartesian(ylim = c(0, ylim_top)) +
    theme_classic(base_size = 8, base_family = "sans") +
    theme(
      axis.text.x  = element_text(size = 8, color = "black"),
      axis.text.y  = element_text(size = 8, color = "black"),
      axis.title.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
      axis.line    = element_line(linewidth = 0.4, color = "black"),
      axis.ticks   = element_line(linewidth = 0.4, color = "black"),
      plot.margin  = margin(t = 2, r = 4, b = 2, l = 4)
    )
  
  if (show_p && nrow(ttest_df2) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      data = ttest_df2,
      label = "label", y.position = "y.position",
      xmin = "group1", xmax = "group2",
      tip.length = 0.01, size = 2.6, bracket.size = 0.4
    )
  }
  return(p)
}


# ---- 8) 產生兩張圖：SD 與 SEM（不顯示 p 值）----
plot_SD  <- make_barplot(summary_df, long_df, ttest_for_plot, err_type = "SD",  show_p = FALSE)
plot_SEM <- make_barplot(summary_df, long_df, ttest_for_plot, err_type = "SEM", show_p = FALSE)

# ---- 9) 存檔 ----
save_all_formats(plot_SD,  paste0(figure_tag, "_SD"))
save_all_formats(plot_SEM, paste0(figure_tag, "_SEM"))

message("Done. Files saved to: ", normalizePath(out_dir))
print(plot_SD)
print(plot_SEM)
