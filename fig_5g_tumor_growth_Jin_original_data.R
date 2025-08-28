# 可選：設定工作資料夾
setwd("C:/Users/danny/Documents/R_project/CSN6_PTEN_plotting/fig_5g_tumor_growth_Jin_original_data")

# =========================================================
# Nature journal style line plots for Tumor Growth assay
# (Jin original data：表內已給各組各日 mean 與 SD)
# - Day 0 一律設為 0 並繪製；其餘天照 input
# - 輸出 SD（一定）與 SEM（若有 n 才能計算）兩種誤差棒
# - 逐日 Welch t-test（若有 n 才能計算）並以 Holm 校正，輸出 CSV
# - 輸出：.tiff (LZW, 600 dpi), .pdf (向量), .jpg (600 dpi)
# =========================================================

# ---- 0) 基本設定 ----
input_path <- "fig_5g_tumor_growth_Jin_original_data_cleaned.xlsx"
sheet_name <- 1
figure_tag <- "Fig_5G_TumorGrowth_Jin_line"
out_dir    <- "figure_outputs"

# 尺寸與造型（沿用你先前設定）
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

# ---- 2) 讀取 & 轉為「每組-每日-mean/sd(/n)」的長表 ----
stopifnot(file.exists(input_path))
df_raw <- readxl::read_excel(path = input_path, sheet = sheet_name)

nm <- names(df_raw); nm_lc <- tolower(nm)
group_candidates <- nm[nm_lc %in% c("group","condition","treatment","genotype","line","cellline")]
group_col <- if (length(group_candidates) >= 1) group_candidates[1] else NA_character_

# 特別處理「兩層表頭」：第一列是 day 標籤，欄名含 "tumor volume" 與 "sd"
special_jin_layout <- (tolower(as.character(df_raw[[1]][1])) %in% c("condition","group")) &&
  any(grepl("tumor", nm_lc)) &&
  any(nm_lc == "sd")

if (special_jin_layout) {
  idx_mean_start <- which(grepl("tumor", nm_lc))[1]
  idx_sd_start   <- which(nm_lc == "sd")[1]
  mean_idx <- idx_mean_start:(idx_sd_start - 1)
  sd_idx   <- idx_sd_start:ncol(df_raw)
  
  # 第 1 列是 day 標籤
  mean_day_labels <- as.character(dplyr::slice(df_raw, 1)[, mean_idx][1, ])
  sd_day_labels   <- as.character(dplyr::slice(df_raw, 1)[, sd_idx][1, ])
  mean_day_num <- readr::parse_number(tolower(mean_day_labels))
  sd_day_num   <- readr::parse_number(tolower(sd_day_labels))
  
  # 資料列從第 2 列開始
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
      sd    = suppressWarnings(as.numeric(sd)),
      n     = NA_real_               # 若檔案未提供 n
    )
} else {
  # 一般長/寬表結構（保留原自動偵測）
  is_long <- any(nm_lc %in% c("day","time","timepoint")) &&
    any(grepl("\\b(mean|avg|value|vol|volume)\\b", nm_lc)) &&
    any(grepl("\\b(sd|stdev|std)\\b", nm_lc))
  
  if (is_long) {
    day_col  <- nm[which(nm_lc %in% c("day","time","timepoint"))][1]
    mean_col <- nm[grep("\\b(mean|avg|value|vol|volume)\\b", nm_lc)][1]
    sd_col   <- nm[grep("\\b(sd|stdev|std)\\b", nm_lc)][1]
    n_col    <- if (any(nm_lc %in% c("n","replicates","count","size","sample","samples"))) {
      nm[which(nm_lc %in% c("n","replicates","count","size","sample","samples"))][1]
    } else NA_character_
    
    long_sum <- df_raw %>%
      { if (!is.na(group_col)) dplyr::select(., all_of(c(group_col, day_col, mean_col, sd_col, n_col))) else dplyr::select(., all_of(c(day_col, mean_col, sd_col, n_col))) } %>%
      dplyr::rename(Group = !!group_col, Day_raw = !!day_col, mean = !!mean_col, sd = !!sd_col) %>%
      { if (!is.na(n_col)) dplyr::rename(., n = !!n_col) else dplyr::mutate(., n = NA_real_) } %>%
      dplyr::mutate(Group = if (is.null(.$Group)) "Group" else as.character(Group),
                    Day_num = readr::parse_number(tolower(Day_raw))) %>%
      dplyr::select(Group, Day_num, mean, sd, n)
  } else {
    day_mean_cols <- nm[grepl("^\\s*day\\s*\\d+\\s*$", nm_lc) |
                          (grepl("^\\s*day\\s*\\d+", nm_lc) & !grepl("sd", nm_lc))]
    if (length(day_mean_cols) == 0) stop("偵測不到日次的平均欄（如 'day 1'）。")
    day_sd_cols <- nm[grepl("^\\s*day\\s*\\d+.*sd", nm_lc) | grepl("sd.*day\\s*\\d+", nm_lc)]
    day_n_cols  <- nm[grepl("^\\s*day\\s*\\d+.*(n|rep)", nm_lc) | grepl("(n|rep).*day\\s*\\d+", nm_lc)]
    
    if (!is.na(group_col)) df_me <- df_raw %>% dplyr::select(all_of(c(group_col, day_mean_cols))) else
      df_me <- df_raw %>% dplyr::mutate(Group = "Group") %>% dplyr::select(Group, all_of(day_mean_cols))
    long_me <- df_me %>% tidyr::pivot_longer(-Group, names_to = "Day_raw", values_to = "mean")
    
    if (length(day_sd_cols) > 0) {
      if (!is.na(group_col)) df_sd <- df_raw %>% dplyr::select(all_of(c(group_col, day_sd_cols))) else
        df_sd <- df_raw %>% dplyr::mutate(Group = "Group") %>% dplyr::select(Group, all_of(day_sd_cols))
      long_sd <- df_sd %>%
        tidyr::pivot_longer(-Group, names_to = "Day_raw_sd", values_to = "sd") %>%
        dplyr::mutate(Day_num = readr::parse_number(tolower(Day_raw_sd))) %>%
        dplyr::select(Group, Day_num, sd)
    } else stop("偵測不到 SD 欄位（如 'day 1 sd' 或 'day1_sd'）。")
    
    if (length(day_n_cols) > 0) {
      if (!is.na(group_col)) df_n <- df_raw %>% dplyr::select(all_of(c(group_col, day_n_cols))) else
        df_n <- df_raw %>% dplyr::mutate(Group = "Group") %>% dplyr::select(Group, all_of(day_n_cols))
      long_n <- df_n %>%
        tidyr::pivot_longer(-Group, names_to = "Day_raw_n", values_to = "n") %>%
        dplyr::mutate(Day_num = readr::parse_number(tolower(Day_raw_n))) %>%
        dplyr::select(Group, Day_num, n)
    } else long_n <- NULL
    
    long_sum <- long_me %>%
      dplyr::mutate(Day_num = readr::parse_number(tolower(Day_raw))) %>%
      dplyr::select(Group, Day_num, mean) %>%
      dplyr::left_join(long_sd, by = c("Group","Day_num")) %>%
      { if (!is.null(long_n)) dplyr::left_join(., long_n, by = c("Group","Day_num")) else dplyr::mutate(., n = NA_real_) }
  }
}

# ---- 2.1) 一致化 Day 標籤、補 D0=0（修正 n 造成的錯誤）----
long_sum <- long_sum %>%
  dplyr::mutate(Day_label = paste0("D", Day_num),
                Group = as.character(Group))

if (!any(long_sum$Day_num == 0)) {
  d0 <- long_sum %>%
    dplyr::distinct(Group) %>%
    dplyr::mutate(
      Day_num   = 0,
      Day_label = "D0",
      mean      = 0,
      sd        = 0,
      n         = NA_real_  # 這裡固定設 NA，避免 mutate(n = n) 的錯誤
    )
  long_sum <- dplyr::bind_rows(d0, long_sum)
} else {
  long_sum <- long_sum %>%
    dplyr::mutate(
      mean = ifelse(Day_num == 0, 0, mean),
      sd   = ifelse(Day_num == 0, 0, sd)
    )
}

# 重新設定 Day 因子順序（含 D0）
day_levels <- long_sum %>%
  dplyr::distinct(Day_label, Day_num) %>%
  dplyr::arrange(Day_num) %>% 
  dplyr::pull(Day_label)
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

# ---- 4) 預先計算 SEM（僅在 n 有值時可用）----
long_sum <- long_sum %>%
  dplyr::mutate(
    sem = dplyr::if_else(is.na(n) | n <= 1, NA_real_, sd / sqrt(n))
  )

# ---- 5) 作圖函式（使用摘要 mean + 指定誤差欄位）----
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

make_lineplot_sum <- function(sum_df, err_col = c("sd","sem"),
                              y_lab = y_label,
                              palette = pal_groups,
                              tag_suffix = "SD",
                              errorbar_color = errorbar_color_mode,
                              errbar_lwd = errorbar_linewidth,
                              mean_pt_size = mean_point_size,
                              mean_pt_stroke = mean_point_stroke) {
  err_col <- match.arg(err_col)
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

# ---- 6) 逐日 Welch t-test（僅在有 n 時可計算；用摘要統計公式）----
welch_from_summary <- function(m1, s1, n1, m2, s2, n2) {
  se2 <- (s1^2)/n1 + (s2^2)/n2
  t   <- (m1 - m2) / sqrt(se2)
  df  <- se2^2 / ( ((s1^2)/n1)^2/(n1-1) + ((s2^2)/n2)^2/(n2-1) )
  p   <- 2 * stats::pt(-abs(t), df)
  list(t = t, df = df, p = p)
}

pairwise_tests_from_summary <- function(sum_df, mode_tag, out_dir_path, figure_tag_prefix) {
  if (all(is.na(sum_df$n))) {
    message("未提供 n，無法由摘要統計計算 t 檢定；將略過 CSV 輸出。")
    return(invisible(NULL))
  }
  days <- sum_df %>% dplyr::distinct(Day_label) %>% dplyr::pull(Day_label)
  out <- list()
  for (dlv in as.character(days)) {
    sub <- sum_df %>% dplyr::filter(Day_label == dlv) %>% dplyr::select(Group, mean, sd, n)
    if (nrow(sub) < 2 || any(is.na(sub$n))) next
    grps <- as.character(sub$Group)
    combs <- utils::combn(grps, 2, simplify = FALSE)
    for (cc in combs) {
      g1 <- cc[1]; g2 <- cc[2]
      r1 <- sub[sub$Group == g1, ][1, ]; r2 <- sub[sub$Group == g2, ][1, ]
      w  <- welch_from_summary(r1$mean, r1$sd, r1$n, r2$mean, r2$sd, r2$n)
      out[[length(out)+1]] <- data.frame(
        Day_label = dlv,
        group1 = g1, group2 = g2,
        n1 = r1$n, n2 = r2$n,
        statistic = w$t, df = w$df, p = w$p,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out)) {
    message("無可用的 n 來計算 t 檢定。")
    return(invisible(NULL))
  }
  res <- dplyr::bind_rows(out) %>%
    dplyr::group_by(Day_label) %>%
    dplyr::mutate(p.adj = stats::p.adjust(p, method = "holm")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p.adj.signif = ggpubr::pvalue_to_asterisk(p.adj))
  if (!dir.exists(out_dir_path)) dir.create(out_dir_path, recursive = TRUE)
  readr::write_csv(
    res %>% dplyr::transmute(
      Day_label,
      comparison = paste(group1, "vs", group2),
      group1, group2, n1, n2, statistic, df, p_value = p, p_adj_holm = p.adj,
      signif_adj = p.adj.signif
    ),
    file.path(out_dir_path, paste0(figure_tag_prefix, "_D0is0_pairwise_ttest_pvalues.csv"))
  )
  invisible(res)
}

# ---- 7) 畫圖與輸出 ----
# SD（一定輸出）
p_SD <- make_lineplot_sum(long_sum, err_col = "sd",  tag_suffix = "D0is0_SD")
save_all_formats(p_SD, paste0(figure_tag, "_", attr(p_SD, "tag_suffix")))
print(p_SD)

# SEM（僅在 n 有值時輸出）
if (any(!is.na(long_sum$sem))) {
  p_SEM <- make_lineplot_sum(long_sum, err_col = "sem", tag_suffix = "D0is0_SEM")
  save_all_formats(p_SEM, paste0(figure_tag, "_", attr(p_SEM, "tag_suffix")))
  print(p_SEM)
} else {
  message("未提供 n，無法計算 SEM；僅輸出 SD 圖。")
}

# 逐日 Welch t-test（僅在 n 有值時）
pairwise_tests_from_summary(long_sum, mode_tag = "D0is0",
                            out_dir_path = out_dir, figure_tag_prefix = figure_tag)

message("Done. Files saved to: ", normalizePath(out_dir))
