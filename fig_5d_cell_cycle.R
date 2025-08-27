# 可選：設定工作資料夾
setwd("C:/Users/danny/Documents/R_project/CSN6_PTEN_plotting/fig_5d_cell_cycle")

# =========================================================
# Nature journal style grouped barplot for Cell Cycle
# - 讀取 fig_5d_cell_cycle_cleaned.xlsx
# - x 軸 = Condition；填色 = Stage（sub G1 / G1 / S / G2M）
# - 產生 SD 與 SEM 兩種誤差棒版本（圖上不顯示 p 值）
# - 輸出：.tiff (LZW, 600 dpi), .pdf (向量), .jpg (600 dpi)
# =========================================================

# ---- 0) 基本設定 ----
input_path   <- "fig_5d_cell_cycle_cleaned.xlsx"
sheet_name   <- 1
y_label      <- "Percentage of cell cycle stages (%)"
figure_tag   <- "Fig_5D_CellCycle_grouped"
out_dir      <- "figure_outputs"

# 若檔案為「長表」請在此指定欄名；若為「寬表」留空即可由腳本自動偵測
condition_col <- NULL
stage_col     <- NULL
value_col     <- NULL

# 寬表時預期的階段欄位（與 Excel 欄名一致即可）
stage_cols_wide <- c("sub G1","G1","S","G2M")

# Nature 單欄寬 ≈ 85 mm（3.35 in）
fig_width_in  <- 3.35
fig_height_in <- 2.4

# ---- 1) 套件安裝 / 載入 ----
pkgs <- c("readxl","tidyverse","ggpubr")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 2) 讀取資料 ----
stopifnot(file.exists(input_path))
df <- readxl::read_excel(path = input_path, sheet = sheet_name)

# ---- 3) 轉長表：支援長/寬兩種輸入格式 ----
is_long <- !is.null(condition_col) && !is.null(stage_col) && !is.null(value_col) &&
  all(c(condition_col, stage_col, value_col) %in% names(df))

if (is_long) {
  long_df <- df %>%
    dplyr::select(all_of(c(condition_col, stage_col, value_col))) %>%
    dplyr::rename(Condition = !!condition_col, Stage = !!stage_col, Value = !!value_col)
} else {
  # 自動偵測寬表：找一欄類別（Condition），其餘數值（Stage）
  if (is.null(condition_col)) {
    guess <- names(df)[tolower(names(df)) %in% c("condition","group","treatment","genotype","line","cellline")]
    condition_col <- if (length(guess) >= 1) guess[1] else names(df)[!purrr::map_lgl(df, is.numeric)][1]
  }
  stopifnot(condition_col %in% names(df))
  num_cols <- names(df)[purrr::map_lgl(df, is.numeric)]
  stage_cols <- if (length(intersect(stage_cols_wide, names(df))) > 0) {
    intersect(stage_cols_wide, names(df))
  } else num_cols
  if (length(stage_cols) < 1) stop("找不到任何數值欄位可作為 Stage。")
  
  long_df <- df %>%
    dplyr::select(all_of(c(condition_col, stage_cols))) %>%
    tidyr::pivot_longer(cols = all_of(stage_cols), names_to = "Stage", values_to = "Value") %>%
    dplyr::rename(Condition = !!condition_col)
}

# ---- 3.1) 正規化 Condition 名稱並固定 x 軸順序 ----
# 固定順序：CMV → PTEN → CSN6+PTEN → CSN6
long_df <- long_df %>%
  dplyr::mutate(
    Condition = dplyr::recode(
      as.character(Condition),
      "PTEN+CSN6" = "CSN6+PTEN",
      "PTEN + CSN6" = "CSN6+PTEN",
      "CSN6 + PTEN" = "CSN6+PTEN",
      "CMV5" = "CMV",
      .default = as.character(Condition)
    )
  )
cond_order <- c("CMV","PTEN","CSN6+PTEN","CSN6")
cond_present <- intersect(cond_order, unique(long_df$Condition))
long_df <- long_df %>% dplyr::mutate(Condition = factor(Condition, levels = cond_present))

# ---- 3.2) Stage 命名與順序（把 G2/M 或 G2-M 視為 G2M）----
long_df <- long_df %>%
  dplyr::mutate(
    Stage = dplyr::recode(
      as.character(Stage),
      "G2/M" = "G2M",
      "G2-M" = "G2M",
      "Sub-G1" = "sub G1",
      "Sub G1" = "sub G1",
      .default = as.character(Stage)
    )
  )
stage_order_pref <- c("sub G1","G1","S","G2M")
stage_levels <- unique(c(intersect(stage_order_pref, unique(long_df$Stage)),
                         setdiff(unique(long_df$Stage), stage_order_pref)))
long_df <- long_df %>%
  dplyr::filter(!is.na(Value)) %>%
  dplyr::mutate(Stage = factor(Stage, levels = stage_levels))

# ---- 4) 摘要統計（Mean, SD, SEM）----
summary_df <- long_df %>%
  dplyr::group_by(Condition, Stage) %>%
  dplyr::summarise(
    n    = dplyr::n(),
    mean = mean(Value),
    sd   = sd(Value),
    sem  = sd / sqrt(n),
    .groups = "drop"
  )

# ---- 5) 調色：NPG（Nature Publishing Group）palette 對應到各 Stage ----
stage_palette_npg <- c(
  "sub G1" = "#E64B35",  # Vermilion
  "G1"     = "#4DBBD5",  # Sky blue
  "S"      = "#00A087",  # Bluish green
  "G2M"    = "#3C5488"   # Indigo
)
# 若還有其他 Stage，從備用色補齊
fallback_colors <- c("#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85")
pal_vec <- setNames(rep(NA_character_, length(levels(summary_df$Stage))), levels(summary_df$Stage))
for (s in names(pal_vec)) {
  pal_vec[s] <- if (s %in% names(stage_palette_npg)) stage_palette_npg[s] else fallback_colors[1]
  if (!(s %in% names(stage_palette_npg))) fallback_colors <- fallback_colors[-1]
}

# ---- 6) 輔助：輸出多格式 ----
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

# ---- 7) 作圖函式（群組並排；SD/SEM 任選；白色點、靠近柱體且避開 error bar）----
make_grouped_bar <- function(summary_df, long_df, err_type = c("SD","SEM")) {
  err_type <- match.arg(err_type)
  err_col  <- if (err_type == "SD") "sd" else "sem"
  df_plot  <- summary_df %>% dplyr::mutate(err = .data[[err_col]])
  
  y_max    <- max(df_plot$mean + dplyr::coalesce(df_plot$err, 0), na.rm = TRUE)
  base_top <- ifelse(is.finite(y_max) && y_max > 0, y_max, 1)
  ylim_top <- base_top * 1.12
  
  dodge <- position_dodge(width = 0.70)
  
  p <- ggplot(df_plot, aes(x = Condition, y = mean, fill = Stage)) +
    geom_col(position = dodge, width = 0.60, color = "black", linewidth = 0.4) +
    geom_errorbar(aes(ymin = mean - err, ymax = mean + err),
                  position = dodge, width = 0.16, linewidth = 0.4) +
    # 資料點：用 shape=Stage 提供 dodge 分組，並以小幅左右抖動避開中央誤差棒
    geom_point(
      data = long_df,
      aes(x = Condition, y = Value, shape = Stage),
      position = position_jitterdodge(
        jitter.width = 0.6,   # 0.04 ~ 0.06 皆可；越小越靠近柱中央
        jitter.height = 0,
        dodge.width  = 0.70    # 必須與上方 dodge 一致
      ),
      size = 0.8, stroke = 0.1,
      fill = "white", color = "black", alpha = 0.9
    ) +
    scale_shape_manual(values = rep(21, length(levels(summary_df$Stage)))) +
    guides(shape = "none") +
    scale_fill_manual(values = pal_vec, name = NULL) +
    labs(x = NULL, y = y_label) +
    coord_cartesian(ylim = c(0, ylim_top)) +
    theme_classic(base_size = 8, base_family = "sans") +
    theme(
      axis.text.x  = element_text(size = 8, color = "black"),
      axis.text.y  = element_text(size = 8, color = "black"),
      axis.title.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
      axis.line    = element_line(linewidth = 0.4, color = "black"),
      axis.ticks   = element_line(linewidth = 0.4, color = "black"),
      legend.position = "right",
      plot.margin  = margin(t = 2, r = 4, b = 2, l = 4)
    )
  return(p)
}

# ---- 8) 產生與輸出圖檔（SD 與 SEM）----
plot_SD  <- make_grouped_bar(summary_df, long_df, err_type = "SD")
plot_SEM <- make_grouped_bar(summary_df, long_df, err_type = "SEM")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save_all_formats(plot_SD,  paste0(figure_tag, "_SD"))
save_all_formats(plot_SEM, paste0(figure_tag, "_SEM"))

message("Done. Files saved to: ", normalizePath(out_dir))
print(plot_SD)
print(plot_SEM)

