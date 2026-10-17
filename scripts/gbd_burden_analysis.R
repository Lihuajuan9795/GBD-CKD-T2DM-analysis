# =====================================================================
#  gbd_burden_analysis.R
#  Global burden of CKD-T2DM attributable to high BMI and low physical activity
#  ---------------------------------------------------------------------
#  产出：
#    - Figure 1  全球 / 5 个 SDI 分组年龄标化率（ASR）时间趋势 1990–2023
#    - Figure 2  1990→2023 绝对归因数三因素分解（人口增长 / 老龄化 / 流行病学改变）
#    - Figure 3  2023 年年龄-性别别归因率
#    - Figure 4  21 个 GBD 区域年均绝对变化（annualised change）
#    - Table 1  + Supplementary Table S1（ASR 与 EAPC 全套）
#    - Supplementary Figure（204 国家 EAPC）
#
#  方法对应稿件 2.4.1（EAPC）、2.4.2（分解）、2.4.5（制图数据准备）
#  单风险 ASR 直接取 GBD 提供的 "Age-standardized" 行；
#  联合效应无 GBD 原生估计，按稿件口径由两个单风险 PAF 以乘性补集法推导
#  （见 joint_paf_sensitivity.R 的三方法敏感性比较）。
# =====================================================================

rm(list = ls())

# ---------------------------------------------------------------------
# 0. 配置：把下面的路径改成你本机的实际路径
# ---------------------------------------------------------------------
GBD_DIR <- "D:/2023同等学力相关/2025年开题相关/需上交资料/提交给鞠老师资料（原始数据、论文终版）/下载数据"
DATA_DIR <- "D:/2023同等学力相关/2025年开题相关/需上交资料/提交给鞠老师资料（原始数据、论文终版）/数据分析版/糖尿病肾病数据/糖尿病肾病数据"
OUT_DIR <- "D:/2023同等学力相关/2025年开题相关/论文写作/投稿论文/output/gbd_burden"

FILE_BMI_REGION <- file.path(GBD_DIR, "第一批数据/第一批数据/全球＋5SDI+7超级地区+21地区_type 2_High body-mass index.csv")
FILE_PA_REGION  <- file.path(GBD_DIR, "第三批数据/第三批数据/全球＋5SDI+7超级地区+21地区_type 2_Low physical activity.csv")
FILE_BMI_CTRY   <- file.path(GBD_DIR, "第一批数据/第一批数据/204国家_type 2_High body-mass index.csv")
FILE_PA_CTRY    <- file.path(GBD_DIR, "第三批数据/第三批数据/204国家_type 2_Low physical activity.csv")

# 总负担（反事实：去除该风险后的总率）——用于计算 PAF 与联合效应
FILE_TOTAL_REGION <- file.path(DATA_DIR, "22地区-type2 both-年龄分组-反事实联合分析.csv")
FILE_TOTAL_SDI    <- file.path(DATA_DIR, "5sdi+type2  both-年龄分组-反事实联合分析.csv")

# 人口与标准人口
FILE_POP_PATTERN  <- file.path(DATA_DIR, "GBD_population_2023", "*.csv")
FILE_AGE_STANDARD <- file.path(DATA_DIR, "2023_age_stand.csv")

# ---------------------------------------------------------------------
# 1. 分析参数
# ---------------------------------------------------------------------
START_YEAR <- 1990
END_YEAR   <- 2023
CAUSE      <- "Chronic kidney disease due to diabetes mellitus type 2"
RISK_BMI   <- "High body-mass index"
RISK_PA    <- "Low physical activity"

MEASURES <- c(
  Deaths = "Deaths",
  DALYs  = "DALYs (Disability-Adjusted Life Years)",
  YLDs   = "YLDs (Years Lived with Disability)",
  YLLs   = "YLLs (Years of Life Lost)"
)

GBD_21_REGIONS <- c(
  "High-income Asia Pacific", "Central Asia", "East Asia", "South Asia", "Southeast Asia",
  "Australasia", "Caribbean", "Central Europe", "Eastern Europe", "Western Europe",
  "Andean Latin America", "Central Latin America", "Southern Latin America",
  "Tropical Latin America", "North Africa and Middle East", "High-income North America",
  "Oceania", "Central Sub-Saharan Africa", "Eastern Sub-Saharan Africa",
  "Southern Sub-Saharan Africa", "Western Sub-Saharan Africa"
)

SDI_5 <- c("High SDI", "High-middle SDI", "Middle SDI", "Low-middle SDI", "Low SDI")

NON_COUNTRY <- c(
  "Global", SDI_5,
  "Sub-Saharan Africa", "High-income", "Latin America and Caribbean",
  "Central Europe, Eastern Europe, and Central Asia",
  "Southeast Asia, East Asia, and Oceania",
  GBD_21_REGIONS
)

# ---------------------------------------------------------------------
# 2. 依赖包
# ---------------------------------------------------------------------
pkgs <- c("data.table", "dplyr", "tidyr", "readr", "stringr", "ggplot2", "scales")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---------------------------------------------------------------------
# 3. 工具函数
# ---------------------------------------------------------------------
clean_age_name <- function(x) {
  x <- as.character(x); x <- trimws(x)
  x <- gsub(" years", "", x, fixed = TRUE)
  x <- gsub(" year",  "", x, fixed = TRUE)
  x <- gsub(" to ", "-", x, fixed = TRUE)
  x <- gsub(" plus", "+", x, fixed = TRUE)
  x
}

clean_sex_name <- function(x) {
  x <- as.character(x); x <- trimws(x)
  dplyr::case_when(
    tolower(x) == "male"   ~ "Male",
    tolower(x) == "female" ~ "Female",
    tolower(x) == "both"   ~ "Both",
    TRUE ~ x
  )
}

standardize_cols <- function(dat) {
  dat <- data.table::as.data.table(dat)
  ren <- c(location = "location_name", sex = "sex_name", age = "age_name",
           measure = "measure_name", metric = "metric_name", year_id = "year")
  for (old in names(ren)) {
    if (old %in% names(dat) && !(ren[[old]] %in% names(dat))) {
      data.table::setnames(dat, old, ren[[old]])
    }
  }
  dat %>%
    dplyr::mutate(
      location_name = ifelse(location_name == "Russia", "Russian Federation", location_name),
      sex_name  = clean_sex_name(sex_name),
      age_name  = clean_age_name(age_name),
      year      = as.integer(year),
      val       = as.numeric(val),
      lower     = as.numeric(lower),
      upper     = as.numeric(upper)
    )
}

read_gbd <- function(path) {
  if (!file.exists(path)) stop("找不到文件：", path)
  dat <- data.table::fread(path, showProgress = FALSE)
  standardize_cols(dat)
}

# ---------------- 标准人口（GBD 2023 世界标准人口，仅保留 ≥25 岁） -----------
age_stand <- readr::read_csv(FILE_AGE_STANDARD, show_col_types = FALSE)
names(age_stand) <- trimws(names(age_stand))
if (!"age_name" %in% names(age_stand)) names(age_stand)[names(age_stand) == "age"] <- "age_name"
age_stand <- age_stand %>% dplyr::mutate(age_name = clean_age_name(age_name))

AGE_GROUPS <- c("25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95+")

weight_df <- age_stand %>%
  dplyr::filter(age_name %in% AGE_GROUPS) %>%
  dplyr::mutate(age_name = factor(age_name, levels = AGE_GROUPS)) %>%
  dplyr::group_by(age_name) %>%
  dplyr::summarise(w_raw = mean(.data[[grep("^w_stand", names(age_stand), value = TRUE)[1]]]),
                   .groups = "drop") %>%
  dplyr::mutate(w = w_raw / sum(w_raw, na.rm = TRUE)) %>%
  dplyr::arrange(age_name)

stopifnot(nrow(weight_df) == length(AGE_GROUPS))
W_STAND <- weight_df$w                       # BAPC 与 ASR 加权都用它
cat("标准人口权重（≥25 岁，15 组）之和 =", sum(W_STAND), "\n")

# ---------------- EAPC：ln(ASR) = a + b*year -----------------------------
calc_eapc <- function(dat) {
  dat <- dat %>%
    dplyr::filter(year >= START_YEAR, year <= END_YEAR, is.finite(val), val > 0) %>%
    dplyr::arrange(year)
  if (nrow(dat) < 5) {
    return(data.frame(EAPC = NA_real_, LCI = NA_real_, UCI = NA_real_,
                      P_value = NA_real_, n_years = nrow(dat)))
  }
  dat$y <- log(dat$val)
  mod <- stats::lm(y ~ year, data = dat)
  co  <- summary(mod)[["coefficients"]]
  beta <- co["year", "Estimate"]; se <- co["year", "Std. Error"]; p <- co["year", "Pr(>|t|)"]
  data.frame(
    EAPC    = (exp(beta) - 1) * 100,
    LCI     = (exp(beta - 1.96 * se) - 1) * 100,
    UCI     = (exp(beta + 1.96 * se) - 1) * 100,
    P_value = p,
    n_years = nrow(dat)
  )
}

# ---------------- 由年龄别率计算 ASR -------------------------------------
asr_from_agespec <- function(rate_df, value_col = "rate") {
  rate_df <- rate_df %>% dplyr::mutate(age_name = factor(clean_age_name(age_name), levels = AGE_GROUPS))
  w_df <- data.frame(age_name = factor(AGE_GROUPS, levels = AGE_GROUPS), w = W_STAND)
  rate_df %>%
    dplyr::left_join(w_df, by = "age_name") %>%
    dplyr::summarise(asr = sum(.data[[value_col]] * w, na.rm = TRUE), .groups = "drop")
}

# ---------------- 联合归因率（乘性补集法 / 主分析）-------------------------
# joint_rate = r_bmi + r_pa - r_bmi * r_pa / r_total
# 等价于：PAF_joint = 1 - (1 - PAF_bmi)(1 - PAF_pa)，再乘总率
add_joint_rate <- function(df) {
  df %>%
    dplyr::mutate(
      rate_bmi   = ifelse(is.na(rate_bmi), 0, rate_bmi),
      rate_pa    = ifelse(is.na(rate_pa),  0, rate_pa),
      paf_bmi    = dplyr::if_else(total_rate > 0, pmin(pmax(rate_bmi / total_rate, 0), 1), NA_real_),
      paf_pa     = dplyr::if_else(total_rate > 0, pmin(pmax(rate_pa  / total_rate, 0), 1), NA_real_),
      joint_paf  = 1 - (1 - paf_bmi) * (1 - paf_pa),
      joint_rate = total_rate * joint_paf
    )
}

# ---------------- 三因素分解（GBD 标准算法，Shapley 型）-------------------
# D = Σ_a pop_a * rate_a
# 分解为 Aging（年龄结构）、Population（人口规模）、Epidemiology（年龄别率）
decompose_3factor <- function(pop1, pop2, rate1, rate2) {
  a <- as.numeric(pop1); b <- as.numeric(pop2)
  r <- as.numeric(rate1) / 1e5; s <- as.numeric(rate2) / 1e5
  p1 <- sum(a); p2 <- sum(b)
  if (p1 <= 0 || p2 <= 0) return(NULL)
  a1 <- a / p1; a2 <- b / p2

  aging <- (sum(a2 * p1 * r) + sum(a2 * p2 * s)) / 3 + (sum(a2 * p1 * s) + sum(a2 * p2 * r)) / 6 -
           (sum(a1 * p1 * r) + sum(a1 * p2 * s)) / 3 - (sum(a1 * p1 * s) + sum(a1 * p2 * r)) / 6
  pop   <- (sum(a1 * p2 * r) + sum(a2 * p2 * s)) / 3 + (sum(a1 * p2 * s) + sum(a2 * p2 * r)) / 6 -
           (sum(a1 * p1 * r) + sum(a2 * p1 * s)) / 3 - (sum(a1 * p1 * s) + sum(a2 * p1 * r)) / 6
  epi   <- (sum(a1 * p1 * s) + sum(a2 * p2 * s)) / 3 + (sum(a1 * p2 * s) + sum(a2 * p1 * s)) / 6 -
           (sum(a1 * p1 * r) + sum(a2 * p2 * r)) / 3 - (sum(a1 * p2 * r) + sum(a2 * p1 * r)) / 6

  D1 <- sum(a * r); D2 <- sum(b * s)
  data.frame(
    D1990 = D1, D2023 = D2, total_change = D2 - D1,
    Aging_abs = aging, Population_abs = pop, Epidemiology_abs = epi,
    Aging_pct = aging / (D2 - D1) * 100,
    Population_pct = pop / (D2 - D1) * 100,
    Epidemiology_pct = epi / (D2 - D1) * 100,
    residual_check = (aging + pop + epi) - (D2 - D1)
  )
}

# ---------------------------------------------------------------------
# 4. 读入数据
# ---------------------------------------------------------------------
cat("读入 GBD 归因数据 ...\n")
attr_region <- dplyr::bind_rows(
  read_gbd(FILE_BMI_REGION),
  read_gbd(FILE_PA_REGION)
)
attr_ctry <- dplyr::bind_rows(
  read_gbd(FILE_BMI_CTRY),
  read_gbd(FILE_PA_CTRY)
)

# 总负担（反事实）文件：无 rei_name 列
total_region <- dplyr::bind_rows(read_gbd(FILE_TOTAL_REGION), read_gbd(FILE_TOTAL_SDI))
if (!"rei_name" %in% names(total_region)) total_region$rei_name <- NA_character_

cat("读入人口数据 ...\n")
pop_files <- sort(Sys.glob(FILE_POP_PATTERN))
if (length(pop_files) == 0)  # 兼容 CSV 位于子目录的情况
  pop_files <- sort(list.files(file.path(DATA_DIR, "GBD_population_2023"),
                               pattern = "\.csv$", recursive = TRUE, full.names = TRUE))
if (length(pop_files) == 0) stop("未找到人口文件：", FILE_POP_PATTERN)
population <- dplyr::bind_rows(lapply(pop_files, read_gbd)) %>%
  dplyr::filter(measure_name == "Population", metric_name == "Number")

# 人口文件为单岁年龄，聚合到 5 岁组（25-29 ... 95+）
age_map <- data.frame(
  age_num = c(0:99, 100),
  age_name = c(rep("<5 years", 5),
               rep(c("5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
                     "45-49","50-54","55-59","60-64","65-69","70-74","75-79",
                     "80-84","85-89","90-94"), each = 5),
               rep("95+ years", 1))
)
population <- population %>%
  dplyr::mutate(age_num = suppressWarnings(as.integer(age_name))) %>%
  dplyr::filter(!is.na(age_num)) %>%
  dplyr::left_join(age_map, by = "age_num") %>%
  dplyr::mutate(age_name = clean_age_name(age_name.y)) %>%
  dplyr::group_by(location_name, sex_name, year, age_name) %>%
  dplyr::summarise(pop = sum(val, na.rm = TRUE), .groups = "drop")

# ---------------------------------------------------------------------
# 5. Figure 1：全球与 5 个 SDI 分组的 ASR 趋势
# ---------------------------------------------------------------------
asr_series <- attr_region %>%
  dplyr::filter(
    sex_name == "Both",
    age_name == "Age-standardized",
    metric_name == "Rate",
    cause_name == CAUSE,
    location_name %in% c("Global", SDI_5),
    year >= START_YEAR & year <= END_YEAR
  ) %>%
  dplyr::transmute(
    measure = stringr::str_extract(measure_name, "^[A-Za-z]+"),
    risk = dplyr::if_else(rei_name == RISK_BMI, "High BMI", "Low physical activity"),
    location_name, year, val, lower, upper
  )

readr::write_csv(asr_series, file.path(OUT_DIR, "Figure1_Global_ASR_trend_HighBMI_LowPA_1990_2023_sourcedata.csv"))

# ---------------------------------------------------------------------
# 6. EAPC（Table 1 / Supplementary Table S1 / 204 国家补充图）
# ---------------------------------------------------------------------
cat("计算 EAPC ...\n")
asr_all <- dplyr::bind_rows(
  attr_ctry   %>% dplyr::filter(sex_name == "Both", age_name == "Age-standardized", metric_name == "Rate"),
  attr_region %>% dplyr::filter(sex_name == "Both", age_name == "Age-standardized", metric_name == "Rate")
) %>%
  dplyr::transmute(
    measure = stringr::str_extract(measure_name, "^[A-Za-z]+"),
    risk = dplyr::if_else(rei_name == RISK_BMI, "High BMI", "Low physical activity"),
    location_name, year, val, lower, upper
  ) %>%
  dplyr::distinct(measure, risk, location_name, year, .keep_all = TRUE)

eapc_table <- asr_all %>%
  dplyr::group_by(measure, risk, location_name) %>%
  dplyr::group_modify(~ calc_eapc(.x)) %>%
  dplyr::ungroup()

endpoint <- asr_all %>%
  dplyr::filter(year %in% c(START_YEAR, END_YEAR)) %>%
  dplyr::transmute(measure, risk, location_name,
                   year = paste0("val_", year),
                   l = sprintf("%.2f (%.2f–%.2f)", val, lower, upper)) %>%
  tidyr::pivot_wider(names_from = year, values_from = l)

eapc_out <- endpoint %>%
  dplyr::left_join(eapc_table, by = c("measure", "risk", "location_name")) %>%
  dplyr::mutate(
    EAPC_fmt = sprintf("%.2f (%.2f–%.2f)", EAPC, LCI, UCI),
    location_type = dplyr::case_when(
      location_name == "Global" ~ "Global",
      location_name %in% SDI_5  ~ "SDI quintile",
      location_name %in% GBD_21_REGIONS ~ "GBD region",
      TRUE ~ "Country/territory"
    )
  )

readr::write_csv(eapc_out, file.path(OUT_DIR, "Table1_and_Supplementary_Table_S1_ASR_EAPC_full.csv"))

# Table 1（正文）：仅 Global + 5 个 SDI，DALYs 与 Deaths
table1 <- eapc_out %>%
  dplyr::filter(location_name %in% c("Global", SDI_5), measure %in% c("DALYs", "Deaths"))
readr::write_csv(table1, file.path(OUT_DIR, "Table1_global_SDI_ASR_EAPC.csv"))

# 204 国家 EAPC 补充图数据
countries <- setdiff(unique(asr_all$location_name), NON_COUNTRY)
cat("识别到国家/地区数 =", length(countries), "（应为 204）\n")
readr::write_csv(
  eapc_out %>% dplyr::filter(location_type == "Country/territory", measure == "DALYs"),
  file.path(OUT_DIR, "Supplementary_Figure_EAPC_204countries_DALYs_sourcedata.csv")
)

# ---------------------------------------------------------------------
# 7. Figure 4：21 个 GBD 区域年均绝对变化（含联合效应）
# ---------------------------------------------------------------------
cat("计算 21 区域年均变化 ...\n")
build_joint_asr_one <- function(loc, measure_target) {
  yr <- c(START_YEAR, END_YEAR)
  get_rate <- function(risk) {
    attr_region %>%
      dplyr::filter(location_name == loc, sex_name == "Both", year %in% yr,
                    age_name %in% AGE_GROUPS, metric_name == "Rate",
                    measure_name == !!measure_target, cause_name == CAUSE, rei_name == risk) %>%
      dplyr::group_by(year, age_name) %>% dplyr::summarise(rate = sum(val, na.rm = TRUE), .groups = "drop")
  }
  tot <- total_region %>%
    dplyr::filter(location_name == loc, sex_name == "Both", year %in% yr,
                  age_name %in% AGE_GROUPS, metric_name == "Rate",
                  measure_name == !!measure_target, cause_name == CAUSE) %>%
    dplyr::group_by(year, age_name) %>% dplyr::summarise(total_rate = sum(val, na.rm = TRUE), .groups = "drop")

  d <- tot %>%
    dplyr::left_join(get_rate(RISK_BMI) %>% dplyr::rename(rate_bmi = rate),
                     by = c("year", "age_name")) %>%
    dplyr::left_join(get_rate(RISK_PA)  %>% dplyr::rename(rate_pa  = rate),
                     by = c("year", "age_name")) %>%
    add_joint_rate()

  d %>%
    dplyr::select(year, age_name, rate_bmi, rate_pa, joint_rate) %>%
    tidyr::pivot_longer(cols = c(rate_bmi, rate_pa, joint_rate),
                        names_to = "risk", values_to = "rate") %>%
    dplyr::mutate(risk = dplyr::recode(risk,
                                       rate_bmi = "High BMI",
                                       rate_pa = "Low physical activity",
                                       joint_rate = "Joint")) %>%
    dplyr::group_by(risk, year) %>%
    dplyr::group_modify(~ asr_from_agespec(.x, "rate")) %>%
    dplyr::ungroup()
}

fig4_list <- list()
for (m in MEASURES) {
  for (loc in GBD_21_REGIONS) {
    x <- try(build_joint_asr_one(loc, m), silent = TRUE)
    if (inherits(x, "try-error") || is.null(x) || nrow(x) == 0) next
    x <- x %>%
      dplyr::mutate(measure = names(MEASURES)[which(MEASURES == m)], region = loc) %>%
      dplyr::select(measure, region, risk, year, asr)
    fig4_list[[length(fig4_list) + 1]] <- x
  }
}
fig4 <- dplyr::bind_rows(fig4_list) %>%
  tidyr::pivot_wider(names_from = year, values_from = asr) %>%
  dplyr::mutate(mean_annual_change_per_100000_per_year =
                  (as.numeric(.data[["2023"]]) - as.numeric(.data[["1990"]])) / (END_YEAR - START_YEAR)) %>%
  dplyr::select(measure, region, risk, mean_annual_change_per_100000_per_year)

readr::write_csv(fig4, file.path(OUT_DIR, "Figure4_Mean_annual_change_21GBD_regions_1990_2023_sourcedata.csv"))

# ---------------------------------------------------------------------
# 8. Figure 3：2023 年年龄-性别别归因率
# ---------------------------------------------------------------------
cat("整理年龄-性别别率 ...\n")
fig3 <- attr_region %>%
  dplyr::filter(location_name == "Global", sex_name %in% c("Male", "Female"),
                year == END_YEAR, age_name %in% AGE_GROUPS,
                metric_name == "Rate", cause_name == CAUSE) %>%
  dplyr::transmute(
    measure = stringr::str_extract(measure_name, "^[A-Za-z]+"),
    risk = dplyr::if_else(rei_name == RISK_BMI, "High BMI", "Low physical activity"),
    sex = sex_name, age_group = age_name,
    rate_per_100000 = val, lower = lower, upper = upper
  )

readr::write_csv(fig3, file.path(OUT_DIR, "Figure3_Age_sex_specific_rates_2023_HighBMI_LowPA_sourcedata.csv"))

# ---------------------------------------------------------------------
# 9. Figure 2：三因素分解（Global + 21 区域）
# ---------------------------------------------------------------------
cat("进行三因素分解 ...\n")
get_pop_vec <- function(loc, yr) {
  population %>%
    dplyr::filter(location_name == loc, sex_name == "Both", year == yr,
                  age_name %in% AGE_GROUPS) %>%
    dplyr::mutate(age_name = factor(age_name, levels = AGE_GROUPS)) %>%
    dplyr::arrange(age_name) %>% dplyr::pull(pop)
}
get_rate_vec <- function(loc, yr, measure_target, risk = NULL) {
  d <- attr_region %>%
    dplyr::filter(location_name == loc, sex_name == "Both", year == yr,
                  age_name %in% AGE_GROUPS, metric_name == "Rate",
                  measure_name == !!measure_target, cause_name == CAUSE)
  if (!is.null(risk)) d <- d %>% dplyr::filter(rei_name == risk)
  d <- d %>% dplyr::group_by(age_name) %>% dplyr::summarise(rate = sum(val, na.rm = TRUE), .groups = "drop")
  d %>% dplyr::mutate(age_name = factor(age_name, levels = AGE_GROUPS)) %>%
    dplyr::arrange(age_name) %>% dplyr::pull(rate)
}
get_joint_vec <- function(loc, yr, measure_target) {
  tot <- total_region %>%
    dplyr::filter(location_name == loc, sex_name == "Both", year == yr,
                  age_name %in% AGE_GROUPS, metric_name == "Rate",
                  measure_name == !!measure_target, cause_name == CAUSE) %>%
    dplyr::group_by(age_name) %>% dplyr::summarise(total_rate = sum(val, na.rm = TRUE), .groups = "drop")
  b <- attr_region %>%
    dplyr::filter(location_name == loc, sex_name == "Both", year == yr,
                  age_name %in% AGE_GROUPS, metric_name == "Rate",
                  measure_name == !!measure_target, cause_name == CAUSE, rei_name == RISK_BMI) %>%
    dplyr::group_by(age_name) %>% dplyr::summarise(rate_bmi = sum(val, na.rm = TRUE), .groups = "drop")
  p <- attr_region %>%
    dplyr::filter(location_name == loc, sex_name == "Both", year == yr,
                  age_name %in% AGE_GROUPS, metric_name == "Rate",
                  measure_name == !!measure_target, cause_name == CAUSE, rei_name == RISK_PA) %>%
    dplyr::group_by(age_name) %>% dplyr::summarise(rate_pa = sum(val, na.rm = TRUE), .groups = "drop")
  d <- tot %>% dplyr::left_join(b, by = "age_name") %>% dplyr::left_join(p, by = "age_name") %>% add_joint_rate()
  d %>% dplyr::mutate(age_name = factor(age_name, levels = AGE_GROUPS)) %>%
    dplyr::arrange(age_name) %>% dplyr::pull(joint_rate)
}

decomp_rows <- list()
for (m in MEASURES) {
  for (loc in c("Global", GBD_21_REGIONS)) {
    pop1 <- get_pop_vec(loc, START_YEAR); pop2 <- get_pop_vec(loc, END_YEAR)
    if (length(pop1) == 0 || length(pop2) == 0) next
    for (rt in c("High BMI", "Low physical activity", "Joint")) {
      r1 <- switch(rt,
                   "High BMI"              = get_rate_vec(loc, START_YEAR, m, RISK_BMI),
                   "Low physical activity" = get_rate_vec(loc, START_YEAR, m, RISK_PA),
                   "Joint"                 = get_joint_vec(loc, START_YEAR, m))
      r2 <- switch(rt,
                   "High BMI"              = get_rate_vec(loc, END_YEAR, m, RISK_BMI),
                   "Low physical activity" = get_rate_vec(loc, END_YEAR, m, RISK_PA),
                   "Joint"                 = get_joint_vec(loc, END_YEAR, m))
      if (length(r1) != length(pop1) || length(r2) != length(pop2)) next
      out <- decompose_3factor(pop1, pop2, r1, r2)
      if (is.null(out)) next
      out$measure <- names(MEASURES)[which(MEASURES == m)]
      out$location <- loc
      out$risk <- rt
      decomp_rows[[length(decomp_rows) + 1]] <- out
    }
  }
}
decomp <- dplyr::bind_rows(decomp_rows)
readr::write_csv(decomp, file.path(OUT_DIR, "Figure2_Decomposition_HighBMI_LowPA_1990_2023_sourcedata.csv"))
cat("分解残差最大值（应 < 1e-6）：", max(abs(decomp$residual_check), na.rm = TRUE), "\n")

# ---------------------------------------------------------------------
# 10. 完成
# ---------------------------------------------------------------------
cat("\n全部输出已写入：", OUT_DIR, "\n")
print(list.files(OUT_DIR))
