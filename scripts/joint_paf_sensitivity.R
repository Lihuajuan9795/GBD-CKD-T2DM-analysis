# =====================================================================
#  joint_paf_sensitivity.R
#  高 BMI 与低体力活动「联合人群归因分数（joint PAF）」估计方法敏感性分析
#  ---------------------------------------------------------------------
#  产出（对应 Supplementary Figures S1–S6 与 Supplementary Table S3）：
#    S1 / S4  DALYs / Deaths：1990–2023 年三种方法的联合 PAF 趋势
#    S2 / S5  DALYs / Deaths：2023 年另两种方法相对主分析的差值（百分点）
#    S3 / S6  DALYs / Deaths：2023 年三种方法的联合 PAF 估计值
#    Table S3 2023 年全球及 5 个 SDI 五分位的三种方法比较
#
#  方法对应稿件 2.4（联合归因负担）与局限性说明：
#    本次使用的 GBD 2023 数据提取中**不含**两风险联合的原生（combined-risk）
#    估计，因此联合 PAF 只能由两个单风险 PAF 推导。联合估计依赖连接两风险的
#    因果结构假设，故比较三种组合方法：
#      (1) 乘性补集法  PAF = 1 − (1 − PAF_bmi)(1 − PAF_pa)     ← 主分析
#      (2) 加法求和法  PAF = PAF_bmi + PAF_pa
#      (3) 保守最大值法 PAF = max(PAF_bmi, PAF_pa)
#    三法之间的差异按**方法学敏感性**报告，不作为实质性结果。
# =====================================================================

rm(list = ls())

# ---------------------------------------------------------------------
# 0. 配置
# ---------------------------------------------------------------------
GBD_DIR <- "D:/2023同等学力相关/2025年开题相关/需上交资料/提交给鞠老师资料（原始数据、论文终版）/下载数据"
DATA_DIR <- "D:/2023同等学力相关/2025年开题相关/需上交资料/提交给鞠老师资料（原始数据、论文终版）/数据分析版/糖尿病肾病数据/糖尿病肾病数据"
OUT_DIR <- "D:/2023同等学力相关/2025年开题相关/论文写作/投稿论文/output/joint_paf"

FILE_BMI_REGION  <- file.path(GBD_DIR, "第一批数据/第一批数据/全球＋5SDI+7超级地区+21地区_type 2_High body-mass index.csv")
FILE_PA_REGION   <- file.path(GBD_DIR, "第三批数据/第三批数据/全球＋5SDI+7超级地区+21地区_type 2_Low physical activity.csv")
FILE_TOTAL_SDI   <- file.path(DATA_DIR, "5sdi+type2  both-年龄分组-反事实联合分析.csv")
FILE_AGE_STANDARD <- file.path(DATA_DIR, "2023_age_stand.csv")

START_YEAR <- 1990
END_YEAR   <- 2023
CAUSE      <- "Chronic kidney disease due to diabetes mellitus type 2"
RISK_BMI   <- "High body-mass index"
RISK_PA    <- "Low physical activity"

# 分析范围：全球 + 5 个 SDI 五分位
ANALYSIS_LOCS <- c("Global", "High SDI", "High-middle SDI", "Middle SDI",
                   "Low-middle SDI", "Low SDI")

# ---------------------------------------------------------------------
# 1. 依赖包
# ---------------------------------------------------------------------
pkgs <- c("data.table", "dplyr", "tidyr", "readr", "ggplot2", "scales", "patchwork")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---------------------------------------------------------------------
# 2. 标准人口与工具函数
# ---------------------------------------------------------------------
clean_age_name <- function(x) {
  x <- as.character(x); x <- trimws(x)
  x <- gsub(" years", "", x, fixed = TRUE)
  x <- gsub(" year",  "", x, fixed = TRUE)
  x <- gsub(" to ", "-", x, fixed = TRUE)
  x <- gsub(" plus", "+", x, fixed = TRUE)
  x
}

AGE_GROUPS <- c("25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95+")

age_stand <- readr::read_csv(FILE_AGE_STANDARD, show_col_types = FALSE)
names(age_stand) <- trimws(names(age_stand))
if (!"age_name" %in% names(age_stand)) names(age_stand)[names(age_stand) == "age"] <- "age_name"
wcol <- grep("^w_stand", names(age_stand), value = TRUE)[1]

weight_df <- age_stand %>%
  dplyr::mutate(age_name = clean_age_name(age_name)) %>%
  dplyr::filter(age_name %in% AGE_GROUPS) %>%
  dplyr::mutate(age_name = factor(age_name, levels = AGE_GROUPS)) %>%
  dplyr::group_by(age_name) %>%
  dplyr::summarise(w_raw = mean(.data[[wcol]], na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(w = w_raw / sum(w_raw, na.rm = TRUE)) %>%
  dplyr::arrange(age_name)
stopifnot(nrow(weight_df) == length(AGE_GROUPS))

read_gbd <- function(path) {
  data.table::fread(path, showProgress = FALSE) %>%
    data.table::as.data.table() %>%
    dplyr::mutate(age_name = clean_age_name(age_name), year = as.integer(year),
                  val = as.numeric(val))
}

# ---------------------------------------------------------------------
# 3. 读入数据
# ---------------------------------------------------------------------
cat("读入年龄别归因率与总负担 ...\n")
attr_df <- dplyr::bind_rows(
    read_gbd(FILE_BMI_REGION),
    read_gbd(FILE_PA_REGION)
  ) %>%
  dplyr::filter(sex_name == "Both", age_name %in% AGE_GROUPS,
                metric_name == "Rate", cause_name == CAUSE,
                rei_name %in% c(RISK_BMI, RISK_PA),
                location_name %in% ANALYSIS_LOCS,
                year >= START_YEAR, year <= END_YEAR)

total_df <- read_gbd(FILE_TOTAL_SDI) %>%
  dplyr::filter(sex_name == "Both", age_name %in% AGE_GROUPS, metric_name == "Rate",
                location_name %in% ANALYSIS_LOCS, year >= START_YEAR, year <= END_YEAR)

get_rate <- function(risk_name) {
  attr_df %>%
    dplyr::filter(rei_name == risk_name) %>%
    dplyr::group_by(location_name, year, age_name, measure_name) %>%
    dplyr::summarise(rate = sum(val, na.rm = TRUE), .groups = "drop")
}
get_total <- function() {
  total_df %>%
    dplyr::group_by(location_name, year, age_name, measure_name) %>%
    dplyr::summarise(total_rate = sum(val, na.rm = TRUE), .groups = "drop")
}

# ---------------------------------------------------------------------
# 4. 计算三种联合 PAF（年龄别 → ASR 层面）
# ---------------------------------------------------------------------
build_joint <- function(measure_target) {
  d <- get_total() %>%
    dplyr::filter(measure_name == !!measure_target) %>%
    dplyr::left_join(get_rate(RISK_BMI) %>% dplyr::filter(measure_name == !!measure_target) %>%
                       dplyr::select(-measure_name) %>% dplyr::rename(rate_bmi = rate),
                     by = c("location_name", "year", "age_name")) %>%
    dplyr::left_join(get_rate(RISK_PA) %>% dplyr::filter(measure_name == !!measure_target) %>%
                       dplyr::select(-measure_name) %>% dplyr::rename(rate_pa = rate),
                     by = c("location_name", "year", "age_name")) %>%
    dplyr::mutate(
      rate_bmi = ifelse(is.na(rate_bmi), 0, rate_bmi),
      rate_pa  = ifelse(is.na(rate_pa),  0, rate_pa),

      # 单风险 PAF（截断到 0–1）
      paf_bmi = dplyr::if_else(total_rate > 0, pmin(pmax(rate_bmi / total_rate, 0), 1), NA_real_),
      paf_pa  = dplyr::if_else(total_rate > 0, pmin(pmax(rate_pa  / total_rate, 0), 1), NA_real_),

      # (1) 乘性补集法（主分析）
      paf_mult = 1 - (1 - paf_bmi) * (1 - paf_pa),
      # (2) 加法求和法
      paf_add  = pmin(paf_bmi + paf_pa, 1),
      # (3) 保守最大值法
      paf_max  = pmax(paf_bmi, paf_pa, na.rm = TRUE)
    )
  d
}

# 年龄别率 → ASR 层面的联合 PAF
asr_paf <- function(d) {
  w_df <- data.frame(age_name = factor(AGE_GROUPS, levels = AGE_GROUPS), w = weight_df$w)
  d2 <- d %>% dplyr::mutate(age_name = factor(age_name, levels = AGE_GROUPS)) %>%
    dplyr::left_join(w_df, by = "age_name")

  num <- d2 %>%
    dplyr::group_by(location_name, year) %>%
    dplyr::summarise(
      total_asr = sum(total_rate * w, na.rm = TRUE),
      joint_mult = sum(total_rate * paf_mult * w, na.rm = TRUE),
      joint_add  = sum(total_rate * paf_add  * w, na.rm = TRUE),
      joint_max  = sum(total_rate * paf_max  * w, na.rm = TRUE),
      w_bmi = sum(total_rate * paf_bmi * w, na.rm = TRUE),
      w_pa  = sum(total_rate * paf_pa  * w, na.rm = TRUE),
      .groups = "drop"
    )

  num %>%
    dplyr::mutate(
      paf_mult = dplyr::if_else(total_asr > 0, joint_mult / total_asr, NA_real_),
      paf_add  = dplyr::if_else(total_asr > 0, joint_add  / total_asr, NA_real_),
      paf_max  = dplyr::if_else(total_asr > 0, joint_max  / total_asr, NA_real_),
      paf_bmi  = dplyr::if_else(total_asr > 0, w_bmi      / total_asr, NA_real_),
      paf_pa   = dplyr::if_else(total_asr > 0, w_pa       / total_asr, NA_real_),
      additive_minus_mult_pp = (paf_add - paf_mult) * 100,
      max_minus_mult_pp      = (paf_max - paf_mult) * 100
    ) %>%
    dplyr::select(location_name, year, paf_mult, paf_add, paf_max, paf_bmi, paf_pa,
                  additive_minus_mult_pp, max_minus_mult_pp)
}

# ---------------------------------------------------------------------
# 5. 输出 S1–S6 与 Table S3
# ---------------------------------------------------------------------
export_measure <- function(measure_target, tag) {
  d  <- build_joint(measure_target)
  ap <- asr_paf(d)

  # S1 / S4：趋势
  readr::write_csv(
    ap %>% dplyr::select(location_name, year, paf_mult, paf_add, paf_max, paf_bmi, paf_pa),
    file.path(OUT_DIR, sprintf("Figure_S%s_joint_PAF_trend_%s_sourcedata.csv", tag$fig_trend, tag$lab))
  )

  # S2 / S5：2023 年差值（百分点）
  readr::write_csv(
    ap %>% dplyr::filter(year == END_YEAR) %>%
      dplyr::select(location_name, paf_mult, paf_add, paf_max,
                    additive_minus_mult_pp, max_minus_mult_pp),
    file.path(OUT_DIR, sprintf("Figure_S%s_joint_PAF_2023_diff_%s_sourcedata.csv", tag$fig_diff, tag$lab))
  )

  # S3 / S6：2023 年三种方法估计值
  readr::write_csv(
    ap %>% dplyr::filter(year == END_YEAR) %>%
      dplyr::select(location_name, paf_mult, paf_add, paf_max),
    file.path(OUT_DIR, sprintf("Figure_S%s_joint_PAF_2023_values_%s_sourcedata.csv", tag$fig_val, tag$lab))
  )

  # Table S3 用
  ap %>% dplyr::filter(year == END_YEAR) %>% dplyr::mutate(measure = tag$mlab)
}

tab_s3 <- dplyr::bind_rows(
  export_measure("DALYs (Disability-Adjusted Life Years)",
                 list(fig_trend = 1, fig_diff = 2, fig_val = 3, lab = "DALY", mlab = "DALYs")),
  export_measure("Deaths",
                 list(fig_trend = 4, fig_diff = 5, fig_val = 6, lab = "Death", mlab = "Deaths"))
)

table_s3 <- tab_s3 %>%
  dplyr::transmute(
    Measure = measure,
    `Location/SDI group` = location_name,
    `High BMI PAF (%)`       = round(paf_bmi * 100, 2),
    `Low PA PAF (%)`         = round(paf_pa  * 100, 2),
    `Joint PAF, multiplicative (%)`     = round(paf_mult * 100, 2),
    `Joint PAF, additive (%)`           = round(paf_add  * 100, 2),
    `Joint PAF, conservative max (%)`   = round(paf_max  * 100, 2)
  )
readr::write_csv(table_s3, file.path(OUT_DIR, "Supplementary_Table_S3_joint_PAF_methods.csv"))

# ---------------------------------------------------------------------
# 6. 作图（S1/S4 趋势、S3/S6 柱状图）
# ---------------------------------------------------------------------
for (mlab in c("DALYs", "Deaths")) {
  lab <- if (mlab == "DALYs") "DALY" else "Death"
  tr <- readr::read_csv(
    file.path(OUT_DIR, sprintf("Figure_S%s_joint_PAF_trend_%s_sourcedata.csv",
                               if (mlab == "DALYs") 1 else 4, lab)),
    show_col_types = FALSE) %>%
    tidyr::pivot_longer(cols = c(paf_mult, paf_add, paf_max),
                        names_to = "method", values_to = "paf") %>%
    dplyr::mutate(method = dplyr::recode(method,
                                         paf_mult = "Multiplicative complement (main)",
                                         paf_add  = "Additive sum",
                                         paf_max  = "Conservative maximum"))

  p <- ggplot2::ggplot(tr, ggplot2::aes(year, paf, colour = method, linetype = method)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(~ location_name) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Year", y = "Joint PAF (age-standardised)", colour = "", linetype = "",
                  title = sprintf("Joint PAF for %s by estimation method, 1990–2023", mlab)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  ggplot2::ggsave(file.path(OUT_DIR, sprintf("Figure_S%s_joint_PAF_trend_%s.png",
                                             if (mlab == "DALYs") 1 else 4, lab)),
                  p, width = 10, height = 6, dpi = 300)

  # ---- S3 / S6：2023 年三种方法估计值（分组柱状图） ----
  fig_val <- if (mlab == "DALYs") 3 else 6
  vals <- readr::read_csv(
    file.path(OUT_DIR, sprintf("Figure_S%s_joint_PAF_2023_values_%s_sourcedata.csv", fig_val, lab)),
    show_col_types = FALSE) %>%
    tidyr::pivot_longer(cols = c(paf_mult, paf_add, paf_max),
                        names_to = "method", values_to = "paf") %>%
    dplyr::mutate(method = factor(dplyr::recode(method,
                                                paf_mult = "Multiplicative complement (main)",
                                                paf_add  = "Additive sum",
                                                paf_max  = "Conservative maximum"),
                                  levels = c("Multiplicative complement (main)",
                                             "Additive sum", "Conservative maximum")))

  pv <- ggplot2::ggplot(vals, ggplot2::aes(location_name, paf, fill = method)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(x = "", y = sprintf("Joint PAF, 2023 (%s, age-standardised)", mlab),
                  fill = "", title = sprintf("Joint PAF in 2023 by estimation method (%s)", mlab)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  ggplot2::ggsave(file.path(OUT_DIR, sprintf("Figure_S%s_joint_PAF_2023_%s.png", fig_val, lab)),
                  pv, width = 9, height = 5.5, dpi = 300)
}

cat("\n输出目录：", OUT_DIR, "\n"); print(list.files(OUT_DIR))
