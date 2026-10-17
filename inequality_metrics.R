# =====================================================================
#  inequality_metrics.R
#  跨国社会经济不平等分析：Slope Index of Inequality (SII) 与
#  Concentration Index (CI)，附集中曲线（concentration curve）
#  ---------------------------------------------------------------------
#  产出：
#    - Figure 6        1990 / 2023 年集中曲线源数据
#    - Figure 6 SII/CI 汇总表（1990 / 2023 × 四项指标 × 高 BMI / 低体力活动）
#    - Supplementary Table S2（SII / CI 全套，含按 SDI 五分位的 2023 年结果）
#
#  方法对应稿件 2.4.4：
#    * SII：以国家按 SDI 排序后的人口累积相对秩为自变量、
#          人口数为权重的加权线性回归斜率，表示 SDI 分布两端 ASR 的绝对差。
#    * CI ：2 × Σ(pop × ASR × (rank − 0.5)) / Σ(pop × ASR)，相对不平等测度。
#    * 正值 = 负担集中于高 SDI 国家；负值 = 集中于低 SDI 国家。
#    * CI 的 95% CI 用 Bootstrap（默认 1000 次）获得。
# =====================================================================

rm(list = ls())

# ---------------------------------------------------------------------
# 0. 配置
# ---------------------------------------------------------------------
GBD_DIR <- "D:/2023同等学力相关/2025年开题相关/需上交资料/提交给鞠老师资料（原始数据、论文终版）/下载数据"
DATA_DIR <- "D:/2023同等学力相关/2025年开题相关/需上交资料/提交给鞠老师资料（原始数据、论文终版）/数据分析版/糖尿病肾病数据/糖尿病肾病数据"
OUT_DIR <- "D:/2023同等学力相关/2025年开题相关/论文写作/投稿论文/output/inequality"

FILE_BMI_CTRY <- file.path(GBD_DIR, "第一批数据/第一批数据/204国家_type 2_High body-mass index.csv")
FILE_PA_CTRY  <- file.path(GBD_DIR, "第三批数据/第三批数据/204国家_type 2_Low physical activity.csv")
FILE_POP_PATTERN <- file.path(DATA_DIR, "GBD_population_2023", "*.csv")
FILE_SDI <- file.path(DATA_DIR, "SDI2023_匹配地区名称后.xlsx")

YEARS    <- c(1990, 2023)
CAUSE    <- "Chronic kidney disease due to diabetes mellitus type 2"
RISK_BMI <- "High body-mass index"
RISK_PA  <- "Low physical activity"

MEASURES <- c(
  DALYs  = "DALYs (Disability-Adjusted Life Years)",
  Deaths = "Deaths",
  YLDs   = "YLDs (Years Lived with Disability)",
  YLLs   = "YLLs (Years of Life Lost)"
)

NON_COUNTRY <- c(
  "Global", "High SDI", "High-middle SDI", "Middle SDI", "Low-middle SDI", "Low SDI",
  "Sub-Saharan Africa", "High-income", "Latin America and Caribbean",
  "Central Europe, Eastern Europe, and Central Asia",
  "Southeast Asia, East Asia, and Oceania",
  "High-income Asia Pacific", "Central Asia", "East Asia", "South Asia", "Southeast Asia",
  "Australasia", "Caribbean", "Central Europe", "Eastern Europe", "Western Europe",
  "Andean Latin America", "Central Latin America", "Southern Latin America",
  "Tropical Latin America", "North Africa and Middle East", "High-income North America",
  "Oceania", "Central Sub-Saharan Africa", "Eastern Sub-Saharan Africa",
  "Southern Sub-Saharan Africa", "Western Sub-Saharan Africa"
)

N_BOOT <- 1000
set.seed(20231017)

# ---------------------------------------------------------------------
# 1. 依赖包
# ---------------------------------------------------------------------
pkgs <- c("data.table", "dplyr", "tidyr", "readr", "readxl", "stringr", "ggplot2")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

read_gbd <- function(path) {
  dat <- data.table::fread(path, showProgress = FALSE)
  data.table::as.data.table(dat) %>%
    dplyr::mutate(year = as.integer(year), val = as.numeric(val),
                  lower = as.numeric(lower), upper = as.numeric(upper))
}

# ---------------------------------------------------------------------
# 2. 读入数据
# ---------------------------------------------------------------------
cat("读入 204 国家 ASR ...\n")
asr <- dplyr::bind_rows(read_gbd(FILE_BMI_CTRY), read_gbd(FILE_PA_CTRY)) %>%
  dplyr::filter(sex_name == "Both", age_name == "Age-standardized",
                metric_name == "Rate", cause_name == CAUSE, year %in% YEARS) %>%
  dplyr::transmute(
    measure = stringr::str_extract(measure_name, "^[A-Za-z]+"),
    risk = dplyr::if_else(rei_name == RISK_BMI, "High BMI", "Low physical activity"),
    location_name, year, val
  )

# 国家人口（全年龄，Both）
cat("读入人口 ...\n")
pop_files <- sort(Sys.glob(FILE_POP_PATTERN))
if (length(pop_files) == 0)  # 兼容 CSV 位于子目录的情况
  pop_files <- sort(list.files(file.path(DATA_DIR, "GBD_population_2023"),
                               pattern = "\.csv$", recursive = TRUE, full.names = TRUE))
pop <- dplyr::bind_rows(lapply(pop_files, read_gbd)) %>%
  dplyr::filter(measure_name == "Population", metric_name == "Number",
                sex_name == "Both", year %in% YEARS) %>%
  dplyr::group_by(location_name, year) %>%
  dplyr::summarise(population = sum(val, na.rm = TRUE), .groups = "drop")

# SDI（宽表 → 长表）
cat("读入 SDI ...\n")
sdi_raw <- readxl::read_excel(FILE_SDI)
names(sdi_raw) <- trimws(names(sdi_raw))
if (!"location_name" %in% names(sdi_raw)) names(sdi_raw)[names(sdi_raw) == "Location"] <- "location_name"
sdi_long <- sdi_raw %>%
  tidyr::pivot_longer(cols = -location_name, names_to = "year", values_to = "sdi") %>%
  dplyr::mutate(year = suppressWarnings(as.integer(year))) %>%
  dplyr::filter(!is.na(year), year %in% YEARS)

# 合并
dat_all <- asr %>%
  dplyr::left_join(pop,      by = c("location_name", "year")) %>%
  dplyr::left_join(sdi_long, by = c("location_name", "year")) %>%
  dplyr::filter(!is.na(val), !is.na(sdi), !is.na(population), population > 0, val >= 0)

missing_sdi <- setdiff(unique(asr$location_name), unique(sdi_long$location_name))
if (length(missing_sdi) > 0) {
  warning("以下地区在 SDI 文件中缺名，已从不平等分析中剔除：",
          paste(head(missing_sdi, 10), collapse = "; "))
}

# ---------------------------------------------------------------------
# 3. 核心函数
# ---------------------------------------------------------------------
add_relative_rank <- function(dat) {
  dat %>%
    dplyr::arrange(sdi) %>%
    dplyr::mutate(
      pop_total = sum(population, na.rm = TRUE),
      pop_cum   = cumsum(population),
      relative_rank = (pop_cum - population / 2) / pop_total
    )
}

calc_sii <- function(dat) {
  dat2 <- dat %>% dplyr::filter(!is.na(val), !is.na(relative_rank),
                                !is.na(population), population > 0)
  if (nrow(dat2) < 5) {
    return(data.frame(SII = NA_real_, SII_LCI = NA_real_, SII_UCI = NA_real_,
                      P_value = NA_real_, n_country = nrow(dat2)))
  }
  fit <- stats::lm(val ~ relative_rank, data = dat2, weights = population)
  ci  <- suppressMessages(stats::confint(fit))
  data.frame(
    SII      = as.numeric(stats::coef(fit)[2]),
    SII_LCI  = as.numeric(ci["relative_rank", 1]),
    SII_UCI  = as.numeric(ci["relative_rank", 2]),
    P_value  = as.numeric(summary(fit)$coefficients["relative_rank", "Pr(>|t|)"]),
    n_country = nrow(dat2)
  )
}

calc_ci_index <- function(dat) {
  dat2 <- dat %>% dplyr::filter(!is.na(val), !is.na(relative_rank),
                                !is.na(population), population > 0, val >= 0)
  if (nrow(dat2) < 5) return(NA_real_)
  numerator   <- 2 * sum(dat2$population * dat2$val * (dat2$relative_rank - 0.5), na.rm = TRUE)
  denominator <- sum(dat2$population * dat2$val, na.rm = TRUE)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  numerator / denominator
}

boot_ci_index <- function(dat, n_boot = N_BOOT) {
  dat2 <- dat %>% dplyr::filter(!is.na(val), !is.na(sdi),
                                !is.na(population), population > 0, val >= 0)
  if (nrow(dat2) < 8) return(c(low = NA_real_, high = NA_real_))
  vals <- replicate(n_boot, {
    idx <- sample(seq_len(nrow(dat2)), replace = TRUE)
    d <- add_relative_rank(dat2[idx, , drop = FALSE])
    calc_ci_index(d)
  })
  out <- stats::quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE)
  c(low = as.numeric(out[1]), high = as.numeric(out[2]))
}

make_concentration_curve <- function(dat) {
  dat2 <- dat %>%
    dplyr::arrange(sdi) %>%
    dplyr::filter(!is.na(val), !is.na(population), population > 0, val >= 0) %>%
    dplyr::mutate(
      pop_total     = sum(population, na.rm = TRUE),
      burden_total  = sum(population * val, na.rm = TRUE),
      pop_cum       = cumsum(population),
      burden_cum    = cumsum(population * val),
      cum_pop_pct    = pop_cum / pop_total,
      cum_burden_pct = burden_cum / burden_total
    )
  dat2 %>% dplyr::select(location_name, sdi, val, population, cum_pop_pct, cum_burden_pct)
}

# ---------------------------------------------------------------------
# 4. 计算 SII / CI（1990 与 2023，四项指标，两个风险）
# ---------------------------------------------------------------------
cat("计算 SII / CI ...\n")
summary_rows <- list()
curve_rows   <- list()

for (mname in names(MEASURES)) {
  m <- MEASURES[[mname]]
  for (rk in c("High BMI", "Low physical activity")) {
    for (yr in YEARS) {
      sub <- dat_all %>% dplyr::filter(measure == mname, risk == rk, year == yr)
      if (nrow(sub) < 10) next
      sub <- add_relative_rank(sub)

      s <- calc_sii(sub)
      ci_est <- calc_ci_index(sub)
      ci_boot <- boot_ci_index(sub)

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        measure = mname, risk = rk, year = yr,
        SII = s$SII, SII_LCI = s$SII_LCI, SII_UCI = s$SII_UCI,
        SII_P = s$P_value,
        CI = ci_est, CI_LCI = ci_boot[["low"]], CI_UCI = ci_boot[["high"]],
        n_country = s$n_country,
        stringsAsFactors = FALSE
      )

      cv <- make_concentration_curve(sub) %>%
        dplyr::mutate(measure = mname, risk = rk, year = yr)
      curve_rows[[length(curve_rows) + 1]] <- cv
    }
  }
}

sii_ci_summary <- dplyr::bind_rows(summary_rows)
curve_data     <- dplyr::bind_rows(curve_rows)

# ---------------------------------------------------------------------
# 5. 按 SDI 五分位分组的 2023 年结果（Supplementary Table S2 下段）
# ---------------------------------------------------------------------
sdi_quint <- dat_all %>%
  dplyr::filter(year == 2023) %>%
  dplyr::group_by(measure, risk) %>%
  dplyr::mutate(sdi_quintile = dplyr::ntile(sdi, 5)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    sdi_group = factor(sdi_quintile, levels = 1:5,
                       labels = c("Low SDI", "Low-middle SDI", "Middle SDI",
                                  "High-middle SDI", "High SDI"))
  )

sdi_group_summary <- sdi_quint %>%
  dplyr::group_by(measure, risk, sdi_group) %>%
  dplyr::summarise(
    n_country = dplyr::n(),
    median_asr = stats::median(val, na.rm = TRUE),
    mean_asr   = mean(val, na.rm = TRUE),
    total_pop  = sum(population, na.rm = TRUE),
    .groups = "drop"
  )

# ---------------------------------------------------------------------
# 6. 输出
# ---------------------------------------------------------------------
fig6 <- curve_data %>%
  dplyr::transmute(measure, risk, year, cum_pop_pct, cum_burden_pct) %>%
  dplyr::arrange(measure, risk, year, cum_pop_pct)
readr::write_csv(fig6, file.path(OUT_DIR,
  "Figure6_Concentration_curves_inequality_1990_2023_sourcedata.csv"))

readr::write_csv(sii_ci_summary, file.path(OUT_DIR,
  "Figure6_Concentration_curves_inequality_1990_2023_SII_CI_summary.csv"))

readr::write_csv(sdi_group_summary, file.path(OUT_DIR,
  "Supplementary_Table_S2_SII_CI_by_SDI_quintile_2023.csv"))

# ---------------------------------------------------------------------
# 7. 快速绘图（可选）
# ---------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- ggplot2::ggplot(fig6, ggplot2::aes(x = cum_pop_pct, y = cum_burden_pct,
                                          colour = factor(year))) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey60") +
    ggplot2::facet_grid(risk ~ measure, scales = "free") +
    ggplot2::labs(x = "Cumulative share of population ranked by SDI",
                  y = "Cumulative share of attributable burden",
                  colour = "Year") +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(OUT_DIR, "Figure6_concentration_curves.png"),
                  p, width = 11, height = 7, dpi = 300)
}

cat("\n=== SII / CI 汇总 ===\n"); print(sii_ci_summary)
cat("\n输出目录：", OUT_DIR, "\n"); print(list.files(OUT_DIR))
