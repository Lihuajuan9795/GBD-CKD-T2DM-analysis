# =====================================================================
#  case_control_paf.R
#  本地县级医院病例对照研究：高 BMI / 低体力活动与 DKD 的
#  比值比（OR）与人群归因分数（PAF，Miettinen 公式）
#  ---------------------------------------------------------------------
#  产出：
#    - Table 2   基线特征（318 例 DKD vs 1,460 例非 DKD 的 T2DM 对照）
#    - Table 3   本地 PAF 与 GBD 2023 中国估计的比较（Figure 8 源数据）
#    - Figure S7 按性别、年龄分层的 PAF 森林图源数据
#    - Figure 8  本地 PAF vs GBD 2023 中国 PAF
#
#  方法对应稿件 2.5：
#    * 病例：主要出院诊断为 ICD-10 E11.2（2 型糖尿病伴肾的并发症）的住院患者；
#      对照：同期主要诊断为 E11 但无 DKD 的住院患者。
#    * 分析单位为**住院人次**（1,583 名患者贡献 1,778 人次，162 人有 ≥2 次住院）；
#      未做剔除重复住院的敏感性分析，已在局限性中说明。
#    * 高 BMI：BMI ≥ 25 kg/m²；低体力活动：标准化问卷判定。
#    * 完整病例分析（complete-case）：高 BMI 229 / 1054；低体力活动 231 / 1053。
#    * 年龄/性别缺失时由身份证号补全（第 7–14 位出生日期，第 17 位奇偶判性别）。
#    * 协变量由 15 个出院诊断字段按 ICD-10 前缀 + 名称关键词判定，无缺失。
#    * 非条件 logistic 回归，校正年龄、性别、高血压、冠心病、心衰、
#      脑卒中后遗症、脂肪肝、贫血与住院天数。
#    * PAF = Pe × (OR − 1) / [1 + Pe × (OR − 1)]，Pe 为对照暴露率；
#      95% CI 用 Bootstrap（1,000 次，分层重抽样）。
#    * **主分析 PAF 用未校正 OR**（与已存档的 PAF分析结果.xlsx 一致）；
#      校正 OR 及其 PAF 作为敏感性分析，写入
#      Table3_sensitivity_adjusted_OR_PAF.csv，并在稿件 2.5 / 3.8 中报告。
#
#  ⚠ 本分析使用去标识化的回顾性常规诊疗数据（伦理批件号 2025gxry11，
#    已豁免知情同意）。原始个体数据不得公开存储。
# =====================================================================

rm(list = ls())

# ---------------------------------------------------------------------
# 0. 配置
# ---------------------------------------------------------------------
LOCAL_DIR <- "D:/2023同等学力相关/2025年开题相关/结果分析/验证数据"
OUT_DIR   <- "D:/2023同等学力相关/2025年开题相关/论文写作/投稿论文/output/case_control"

FILE_RAW     <- file.path(LOCAL_DIR, "糖尿病患者2023-2024原始数据.xlsx")
FILE_CASES   <- file.path(LOCAL_DIR, "2型-DKD糖尿病肾病 3筛改（分析版）.xlsx")
FILE_CONTROLS <- file.path(LOCAL_DIR, "非DKD的糖尿病患者2023-2024  5筛（分析版）.xlsx")

# GBD 2023 中国（2023 年）PAF，用于 Table 3 / Figure 8 的外部比较
GBD_CHINA_PAF <- data.frame(
  exposure = c("High BMI", "Low physical activity"),
  gbd_deaths_paf = c(32.5, 8.1), gbd_deaths_l = c(14.1, 3.1), gbd_deaths_u = c(56.8, 13.2),
  gbd_dalys_paf  = c(37.5, 7.3), gbd_dalys_l  = c(16.2, 2.8), gbd_dalys_u  = c(63.2, 12.2),
  stringsAsFactors = FALSE
)

N_BOOT <- 1000
set.seed(20250101)

# ---------------------------------------------------------------------
# 1. 依赖包
# ---------------------------------------------------------------------
pkgs <- c("dplyr", "readxl", "tidyr", "readr", "stringr", "broom", "ggplot2", "scales")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---------------------------------------------------------------------
# 2. 读入数据
# ---------------------------------------------------------------------
raw <- readxl::read_excel(FILE_RAW, sheet = 1)
cases    <- readxl::read_excel(FILE_CASES,    sheet = 1)
controls <- readxl::read_excel(FILE_CONTROLS, sheet = 1)

# 原表部分列名带尾随空格（例如 "BMI "），统一清理
for (nm in c("raw", "cases", "controls")) {
  x <- get(nm); names(x) <- trimws(names(x)); assign(nm, x)
}

cat("原始住院记录：", nrow(raw), " 条\n")
cat("病例表：", nrow(cases), " | 对照表：", nrow(controls), "\n")

norm_time <- function(x) as.character(x)

raw <- raw %>%
  dplyr::mutate(admit_key = norm_time(!!dplyr::sym("入院时间")),
                disch_key = norm_time(!!dplyr::sym("出院时间")))
cases <- cases %>%
  dplyr::mutate(admit_key = norm_time(!!dplyr::sym("入院时间")),
                disch_key = norm_time(!!dplyr::sym("出院时间")))
controls <- controls %>%
  dplyr::mutate(admit_key = norm_time(!!dplyr::sym("入院时间")),
                disch_key = norm_time(!!dplyr::sym("出院时间")))

# ---------------------------------------------------------------------
# 3. 病例 / 对照标识
#   匹配键 = (入院时间, 出院时间)。对照表的「住院号」是 1..n 的序号，
#   不是真实住院号，不能用于匹配。若一条记录同时出现在两表，按 DKD 优先。
# ---------------------------------------------------------------------
key <- function(d) paste0(d$admit_key, "||", d$disch_key)

case_keys <- key(cases)
ctrl_keys <- key(controls)
dup_keys  <- intersect(case_keys, ctrl_keys)
if (length(dup_keys) > 0) {
  cat("同时出现在病例表与对照表的记录：", length(dup_keys), " 条，按 DKD 优先从对照剔除。\n")
}

ctrl_keep <- controls[!ctrl_keys %in% dup_keys, ]

case_flag <- data.frame(k = case_keys,    group = "DKD",     stringsAsFactors = FALSE)
ctrl_flag <- data.frame(k = key(ctrl_keep), group = "Non-DKD", stringsAsFactors = FALSE)
flag <- dplyr::bind_rows(case_flag, ctrl_flag)

# 暴露与人口学信息（来自筛查后的病例/对照表）
expo <- dplyr::bind_rows(
  cases     %>% dplyr::mutate(group = "DKD")     %>% dplyr::select(admit_key, disch_key, group,
                                                                    性别, 年龄, BMI, 是否低体力活动),
  ctrl_keep %>% dplyr::mutate(group = "Non-DKD") %>% dplyr::select(admit_key, disch_key, group,
                                                                    性别, 年龄, BMI, 是否低体力活动)
)

dat <- raw %>%
  dplyr::mutate(k = key(raw)) %>%
  dplyr::left_join(flag, by = "k") %>%
  dplyr::filter(!is.na(group)) %>%
  dplyr::left_join(expo %>% dplyr::select(admit_key, disch_key, 性别, 年龄, BMI, 是否低体力活动),
                   by = c("admit_key", "disch_key"))

cat("匹配成功：", nrow(dat), " 条（目标 1,778 = 318 + 1,460）\n")
print(table(dat$group))

# 唯一患者数 / 重复住院情况
# 注意：「住院号」每次住院都不同（1,778 条记录 = 1,778 个住院号），不能用于计数患者；
#       唯一患者应以「身份证号」为准（1,583 名患者，162 人有 ≥2 次住院）。
pid <- if ("身份证号" %in% names(dat)) as.character(dat$身份证号) else as.character(dat$住院号)
stopifnot(sum(is.na(pid)) == 0)
n_patients <- length(unique(pid))
n_repeat   <- sum(table(pid) >= 2)
cat("唯一患者数：", n_patients, " | 有 ≥2 次住院的患者数：", n_repeat, "\n")

# ---------------------------------------------------------------------
# 4. 由身份证号补全年龄与性别
# ---------------------------------------------------------------------
parse_id <- function(id) {
  id <- as.character(id)
  ok <- nchar(id) >= 17
  birth <- ifelse(ok, substr(id, 7, 14), NA_character_)
  sexdig <- ifelse(ok, suppressWarnings(as.integer(substr(id, 17, 17))), NA_integer_)
  data.frame(birth = birth, sex_from_id = ifelse(is.na(sexdig), NA_character_,
                                                 ifelse(sexdig %% 2 == 1, "男", "女")),
             stringsAsFactors = FALSE)
}

idinfo <- parse_id(dat$身份证号)
ref_year <- suppressWarnings(
  as.numeric(stringr::str_sub(dat$disch_key, 1, 4))
) %>% max(na.rm = TRUE)

dat <- dat %>%
  dplyr::mutate(
    sex_raw = as.character(性别),
    age_raw = suppressWarnings(as.numeric(年龄)),
    birth   = as.character(idinfo$birth),
    sex_id  = idinfo$sex_from_id,
    age_from_id = ref_year - suppressWarnings(as.numeric(stringr::str_sub(birth, 1, 4))),
    sex = dplyr::coalesce(ifelse(is.na(sex_raw) | sex_raw == "", NA_character_, sex_raw), sex_id),
    age = dplyr::coalesce(age_raw, age_from_id)
  )
cat("年龄缺失：", sum(is.na(dat$age)), " | 性别缺失：", sum(is.na(dat$sex)), "\n")

# ---------------------------------------------------------------------
# 5. 暴露与协变量
# ---------------------------------------------------------------------
dat <- dat %>%
  dplyr::mutate(
    bmi_val = suppressWarnings(as.numeric(as.character(BMI))),
    high_bmi = dplyr::if_else(!is.na(bmi_val), bmi_val >= 25, NA),
    low_pa   = dplyr::case_when(
      是否低体力活动 %in% c("是", "yes", "Yes", "Y", "1") ~ TRUE,
      是否低体力活动 %in% c("否", "no",  "No",  "N", "0") ~ FALSE,
      TRUE ~ NA
    ),
    los = suppressWarnings(as.numeric(as.character(住院天数))),
    case = ifelse(group == "DKD", 1, 0)
  )

# 15 个出院诊断字段
dx_code_cols <- sprintf("出院第%d诊断编码", 1:15)
dx_name_cols <- sprintf("出院第%d诊断名称", 1:15)
dx_code_cols <- intersect(dx_code_cols, names(dat))
dx_name_cols <- intersect(dx_name_cols, names(dat))

has_any <- function(df, pattern_code = NULL, pattern_name = NULL) {
  hit <- rep(FALSE, nrow(df))
  if (!is.null(pattern_code)) {
    for (cc in dx_code_cols) {
      v <- as.character(df[[cc]])
      hit <- hit | (!is.na(v) & grepl(pattern_code, v, ignore.case = TRUE))
    }
  }
  if (!is.null(pattern_name)) {
    for (cc in dx_name_cols) {
      v <- as.character(df[[cc]])
      hit <- hit | (!is.na(v) & grepl(pattern_name, v, ignore.case = TRUE))
    }
  }
  as.integer(hit)
}

dat <- dat %>%
  dplyr::mutate(
    hypertension = has_any(., "^(I10|I11|I13|I15)", "高血压"),
    chd          = has_any(., "^I2[0-5]",           "冠心病|缺血性心脏病|心绞痛|心肌梗死"),
    heart_fail   = has_any(., "^I50",               "心衰|心力衰竭"),
    stroke_seq   = has_any(., "^I69",               "脑梗死後遗症|脑梗死后遗症|脑出血后遗症|脑血管病后遗症"),
    fatty_liver  = has_any(., "^(K76\\.0|K74\\.6)", "脂肪肝"),
    anaemia      = has_any(., "^D(5[0-9]|6[0-4])",  "贫血")
  )

# ---------------------------------------------------------------------
# 6. Table 2：基线特征
# ---------------------------------------------------------------------
fmt_p <- function(p) ifelse(is.na(p), "–", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

baseline_row <- function(label, x_case, x_ctrl, test) {
  data.frame(Characteristic = label,
             `DKD cases (n=318)` = x_case,
             `Non-DKD controls (n=1460)` = x_ctrl,
             `P value` = test, stringsAsFactors = FALSE, check.names = FALSE)
}

tab2 <- list()
tab2[[1]] <- baseline_row("Hospital admissions, n",
                          sum(dat$case == 1), sum(dat$case == 0), "–")

tt <- try(stats::t.test(age ~ case, data = dat), silent = TRUE)
tab2[[2]] <- baseline_row("Age at admission, years, mean ± SD",
                          sprintf("%.1f ± %.1f", mean(dat$age[dat$case == 1], na.rm = TRUE),
                                  stats::sd(dat$age[dat$case == 1], na.rm = TRUE)),
                          sprintf("%.1f ± %.1f", mean(dat$age[dat$case == 0], na.rm = TRUE),
                                  stats::sd(dat$age[dat$case == 0], na.rm = TRUE)),
                          fmt_p(if (inherits(tt, "try-error")) NA else tt$p.value))

sex_tab <- table(dat$sex, dat$case)
ct <- try(stats::chisq.test(sex_tab), silent = TRUE)
tab2[[3]] <- baseline_row("Male, n (%)",
                          sprintf("%d (%.1f)", sum(dat$sex == "男" & dat$case == 1, na.rm = TRUE),
                                  mean(dat$sex[dat$case == 1] == "男", na.rm = TRUE) * 100),
                          sprintf("%d (%.1f)", sum(dat$sex == "男" & dat$case == 0, na.rm = TRUE),
                                  mean(dat$sex[dat$case == 0] == "男", na.rm = TRUE) * 100),
                          fmt_p(if (inherits(ct, "try-error")) NA else ct$p.value))

add_expo_row <- function(label, var) {
  v <- dat[[var]]
  a <- mean(v[dat$case == 1], na.rm = TRUE) * 100
  b <- mean(v[dat$case == 0], na.rm = TRUE) * 100
  na <- sum(!is.na(v) & dat$case == 1); nb <- sum(!is.na(v) & dat$case == 0)
  tb <- table(factor(v, levels = c(TRUE, FALSE)), dat$case)
  ct2 <- try(stats::chisq.test(tb), silent = TRUE)
  baseline_row(label,
               sprintf("%d/%d (%.1f)", sum(v[dat$case == 1], na.rm = TRUE), na, a),
               sprintf("%d/%d (%.1f)", sum(v[dat$case == 0], na.rm = TRUE), nb, b),
               fmt_p(if (inherits(ct2, "try-error")) NA else ct2$p.value))
}
tab2[[4]] <- add_expo_row("High BMI (BMI ≥25 kg/m²), n/N (%)", "high_bmi")
tab2[[5]] <- add_expo_row("Low physical activity, n/N (%)", "low_pa")

add_comorb_row <- function(label, var) {
  v <- dat[[var]]
  ct3 <- try(stats::chisq.test(table(v, dat$case)), silent = TRUE)
  baseline_row(label,
               sprintf("%d (%.1f)", sum(v[dat$case == 1], na.rm = TRUE), mean(v[dat$case == 1]) * 100),
               sprintf("%d (%.1f)", sum(v[dat$case == 0], na.rm = TRUE), mean(v[dat$case == 0]) * 100),
               fmt_p(if (inherits(ct3, "try-error")) NA else ct3$p.value))
}
for (nm in c("hypertension", "chd", "heart_fail", "stroke_seq", "fatty_liver", "anaemia")) {
  lab <- c(hypertension = "Hypertension, n (%)", chd = "Coronary heart disease, n (%)",
           heart_fail = "Heart failure, n (%)", stroke_seq = "Stroke sequelae, n (%)",
           fatty_liver = "Fatty liver, n (%)", anaemia = "Anaemia, n (%)")[nm]
  tab2[[length(tab2) + 1]] <- add_comorb_row(lab, nm)
}

tt2 <- try(stats::t.test(los ~ case, data = dat), silent = TRUE)
tab2[[length(tab2) + 1]] <- baseline_row(
  "Length of stay, days, mean ± SD",
  sprintf("%.1f ± %.1f", mean(dat$los[dat$case == 1], na.rm = TRUE), stats::sd(dat$los[dat$case == 1], na.rm = TRUE)),
  sprintf("%.1f ± %.1f", mean(dat$los[dat$case == 0], na.rm = TRUE), stats::sd(dat$los[dat$case == 0], na.rm = TRUE)),
  fmt_p(if (inherits(tt2, "try-error")) NA else tt2$p.value))

table2 <- dplyr::bind_rows(tab2)
readr::write_csv(table2, file.path(OUT_DIR, "Table2_local_baseline_characteristics.csv"))
cat("\n=== Table 2 ===\n"); print(table2, row.names = FALSE)

# ---------------------------------------------------------------------
# 7. OR 与 PAF
# ---------------------------------------------------------------------
ADJ_VARS <- c("age", "sex", "hypertension", "chd", "heart_fail",
              "stroke_seq", "fatty_liver", "anaemia", "los")

fit_or <- function(d, exposure, adjusted = TRUE) {
  d2 <- d %>% dplyr::filter(!is.na(.data[[exposure]]))
  if (adjusted) {
    rhs <- paste(c(exposure, ADJ_VARS), collapse = " + ")
    d2 <- d2 %>% dplyr::filter(!is.na(age), !is.na(sex), !is.na(los))
    if (nrow(d2) < 30) return(NULL)
    m <- stats::glm(stats::as.formula(paste("case ~", rhs)), data = d2, family = stats::binomial())
  } else {
    if (nrow(d2) < 30) return(NULL)
    m <- stats::glm(stats::as.formula(paste("case ~", exposure)), data = d2, family = stats::binomial())
  }
  s <- summary(m)$coefficients
  or <- exp(stats::coef(m)[exposure])
  ci <- suppressMessages(stats::confint(m))[exposure, ]
  data.frame(OR = or, OR_LCI = exp(ci[1]), OR_UCI = exp(ci[2]),
             P_value = s[exposure, "Pr(>|z|)"], n_case = sum(d2$case == 1),
             n_ctrl = sum(d2$case == 0))
}

# Miettinen PAF：Pe = 对照暴露率；稿件 Table 3 报告的是**未校正 OR** 的 PAF
miettinen_paf <- function(pe, or) {
  x <- pe * (or - 1)
  x / (1 + x)
}

bootstrap_paf <- function(d, exposure, n_boot = N_BOOT) {
  d2 <- d %>% dplyr::filter(!is.na(.data[[exposure]]))
  idx_case <- which(d2$case == 1); idx_ctrl <- which(d2$case == 0)
  vals <- replicate(n_boot, {
    s_case <- sample(idx_case, length(idx_case), replace = TRUE)
    s_ctrl <- sample(idx_ctrl, length(idx_ctrl), replace = TRUE)
    sb <- d2[c(s_case, s_ctrl), ]
    if (nrow(sb) == 0 || length(unique(sb$case)) < 2) return(NA_real_)
    tab <- table(sb[[exposure]], sb$case)
    if (any(dim(tab) < 2)) return(NA_real_)
    a <- tab["TRUE", "1"]; b <- tab["TRUE", "0"]      # 暴露：病例 / 对照
    cc <- tab["FALSE", "1"]; dd <- tab["FALSE", "0"]
    if ((a + cc) == 0 || (b + dd) == 0) return(NA_real_)
    pe <- b / (b + dd)
    or <- (a * dd) / (cc * b)
    if (!is.finite(or) || or <= 0) return(NA_real_)
    miettinen_paf(pe, or) * 100
  })
  stats::quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE)
}

# 校正模型的 PAF 及其 Bootstrap CI：Pe 取（重抽样后的）对照暴露率，
# OR 取该次重样本上重新拟合的多因素 logistic 回归的暴露系数。
# 稿件口径：主分析用**未校正 OR** 的 PAF，校正 OR 的 PAF 作为敏感性分析。
bootstrap_paf_glm <- function(d, exposure, adjusted = TRUE, n_boot = N_BOOT) {
  use <- c(exposure, if (adjusted) ADJ_VARS)
  d2  <- d %>% dplyr::filter(stats::complete.cases(dplyr::select(., dplyr::all_of(use))))
  idx_case <- which(d2$case == 1); idx_ctrl <- which(d2$case == 0)
  vals <- replicate(n_boot, {
    sb <- d2[c(sample(idx_case, length(idx_case), replace = TRUE),
               sample(idx_ctrl, length(idx_ctrl), replace = TRUE)), ]
    if (length(unique(sb$case)) < 2) return(NA_real_)
    pe <- mean(sb[[exposure]][sb$case == 0], na.rm = TRUE)
    if (!is.finite(pe) || pe <= 0 || pe >= 1) return(NA_real_)
    m  <- try(stats::glm(stats::as.formula(paste("case ~",
                   paste(c(exposure, if (adjusted) ADJ_VARS), collapse = " + "))),
                 data = sb, family = stats::binomial()), silent = TRUE)
    if (inherits(m, "try-error")) return(NA_real_)
    or <- exp(stats::coef(m)[[exposure]])
    if (!is.finite(or) || or <= 0) return(NA_real_)
    miettinen_paf(pe, or) * 100
  })
  stats::quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE)
}

exposures <- c(high_bmi = "High BMI", low_pa = "Low physical activity")
paf_rows <- list(); or_rows <- list()

for (ex in names(exposures)) {
  d2 <- dat %>% dplyr::filter(!is.na(.data[[ex]]))
  a <- sum(d2[[ex]] & d2$case == 1); b <- sum(d2[[ex]] & d2$case == 0)
  cc <- sum(!d2[[ex]] & d2$case == 1); dd <- sum(!d2[[ex]] & d2$case == 0)
  pe <- b / (b + dd)
  or_crude <- (a * dd) / (cc * b)
  paf_pt <- miettinen_paf(pe, or_crude) * 100

  bt <- bootstrap_paf(dat, ex)                       # 未校正 OR 的 PAF CI（主分析）
  or_un <- fit_or(dat, ex, adjusted = FALSE)
  or_ad <- fit_or(dat, ex, adjusted = TRUE)

  # 敏感性分析：用校正 OR 重新计算 PAF 及其 Bootstrap CI
  pe_adj <- mean(dat[[ex]][dat$case == 0], na.rm = TRUE)
  paf_adj_pt <- miettinen_paf(pe_adj, or_ad$OR) * 100
  bt_adj <- bootstrap_paf_glm(dat, ex, adjusted = TRUE)

  paf_rows[[ex]] <- data.frame(
    Exposure = exposures[[ex]],
    `Local PAF (%; 95% CI)` = sprintf("%.1f (%.1f to %.1f)", paf_pt, bt[[1]], bt[[2]]),
    `PAF point estimate (%)` = paf_pt, `PAF LCI (%)` = bt[[1]], `PAF UCI (%)` = bt[[2]],
    `Unadjusted OR (95% CI)` = sprintf("%.3f (%.3f–%.3f)", or_un$OR, or_un$OR_LCI, or_un$OR_UCI),
    `Adjusted OR (95% CI)`   = sprintf("%.3f (%.3f–%.3f)", or_ad$OR, or_ad$OR_LCI, or_ad$OR_UCI),
    `P value (unadjusted)` = or_un$P_value,
    `P value (adjusted)` = or_ad$P_value,
    `PAF from adjusted OR (%)`   = paf_adj_pt,
    `PAF adj LCI (%)`            = bt_adj[[1]],
    `PAF adj UCI (%)`            = bt_adj[[2]],
    `Pe (control exposure %)` = pe * 100,
    `n cases / n controls analysed` = sprintf("%d / %d", or_un$n_case, or_un$n_ctrl),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  or_rows[[ex]] <- data.frame(Exposure = exposures[[ex]],
                              OR_unadj = or_un$OR, LCI = or_un$OR_LCI, UCI = or_un$OR_UCI,
                              OR_adj = or_ad$OR, LCI_adj = or_ad$OR_LCI, UCI_adj = or_ad$OR_UCI)
}

paf_tab <- dplyr::bind_rows(paf_rows)
or_tab  <- dplyr::bind_rows(or_rows)
readr::write_csv(or_tab, file.path(OUT_DIR, "Local_OR_unadjusted_and_adjusted.csv"))
cat("\n=== OR（未校正 / 校正） ===\n"); print(or_tab, row.names = FALSE)

# Table 3：本地 PAF vs GBD 2023 中国
table3 <- paf_tab %>%
  dplyr::select(Exposure, `Local PAF (%; 95% CI)`) %>%
  dplyr::left_join(
    GBD_CHINA_PAF %>% dplyr::rename(Exposure = exposure) %>%
      dplyr::transmute(
        Exposure,
        `GBD 2023 China deaths PAF (%; 95% UI)` = sprintf("%.1f (%.1f to %.1f)",
                                                          gbd_deaths_paf, gbd_deaths_l, gbd_deaths_u),
        `GBD 2023 China DALYs PAF (%; 95% UI)`  = sprintf("%.1f (%.1f to %.1f)",
                                                          gbd_dalys_paf, gbd_dalys_l, gbd_dalys_u)
      ),
    by = "Exposure"
  )
readr::write_csv(table3, file.path(OUT_DIR, "Table3_local_PAF_vs_GBD_China.csv"))

# 敏感性分析：校正 OR 对应的 PAF（稿件 2.5 / 3.8 报告）
sens <- paf_tab %>%
  dplyr::transmute(
    Exposure,
    `Unadjusted OR (95% CI)`,
    `PAF from unadjusted OR, main analysis (%)` = sprintf("%.1f (%.1f to %.1f)",
                                                          `PAF point estimate (%)`,
                                                          `PAF LCI (%)`, `PAF UCI (%)`),
    `Adjusted OR (95% CI)`,
    `PAF from adjusted OR, sensitivity (%)` = sprintf("%.1f (%.1f to %.1f)",
                                                      `PAF from adjusted OR (%)`,
                                                      `PAF adj LCI (%)`, `PAF adj UCI (%)`)
  )
readr::write_csv(sens, file.path(OUT_DIR, "Table3_sensitivity_adjusted_OR_PAF.csv"))
cat("\n=== 敏感性分析：校正 OR 的 PAF ===\n"); print(sens, row.names = FALSE)
cat("\n=== Table 3 ===\n"); print(table3, row.names = FALSE)

# Figure 8 源数据
fig8 <- dplyr::bind_rows(
  paf_tab %>% dplyr::transmute(
    figure = "Figure 8", risk = Exposure, outcome = "Local (hospital case-control)",
    source = "Local", paf = `PAF point estimate (%)`, lci = `PAF LCI (%)`, uci = `PAF UCI (%)`),
  GBD_CHINA_PAF %>% dplyr::transmute(
    figure = "Figure 8", risk = exposure, outcome = "Deaths",
    source = "GBD 2023 China", paf = gbd_deaths_paf, lci = gbd_deaths_l, uci = gbd_deaths_u),
  GBD_CHINA_PAF %>% dplyr::transmute(
    figure = "Figure 8", risk = exposure, outcome = "DALYs",
    source = "GBD 2023 China", paf = gbd_dalys_paf, lci = gbd_dalys_l, uci = gbd_dalys_u)
)
readr::write_csv(fig8, file.path(OUT_DIR, "Figure8_Local_PAF_vs_GBD_China_2023_sourcedata.csv"))

# ---------------------------------------------------------------------
# 8. Figure S7：按性别与年龄分层的 PAF
# ---------------------------------------------------------------------
dat <- dat %>% dplyr::mutate(
  age_group = dplyr::case_when(
    is.na(age) ~ NA_character_,
    age < 60  ~ "<60",
    age < 75  ~ "60–74",
    TRUE      ~ "≥75"
  )
)

strat_one <- function(d2, ex, label) {
  if (is.null(d2) || nrow(d2) == 0) return(NULL)
  v <- d2[[ex]]
  a <- sum(v & d2$case == 1, na.rm = TRUE); b <- sum(v & d2$case == 0, na.rm = TRUE)
  cc <- sum(!v & d2$case == 1, na.rm = TRUE); dd <- sum(!v & d2$case == 0, na.rm = TRUE)
  if ((b + dd) == 0 || (a + cc) == 0 || (cc * b) == 0) return(NULL)
  pe <- b / (b + dd); or <- (a * dd) / (cc * b)
  data.frame(stratum = label, exposure = exposures[[ex]],
             paf = miettinen_paf(pe, or) * 100, or = or,
             n_case = a + cc, n_ctrl = b + dd, stringsAsFactors = FALSE)
}

strat_rows <- list()
for (ex in names(exposures)) {
  for (lv in c("男", "女")) {
    lbl <- if (lv == "男") "Sex: Male" else "Sex: Female"
    r <- strat_one(dat %>% dplyr::filter(sex == lv, !is.na(.data[[ex]])), ex, lbl)
    if (!is.null(r)) strat_rows[[length(strat_rows) + 1]] <- r
  }
  for (ag in c("<60", "60–74", "≥75")) {
    r <- strat_one(dat %>% dplyr::filter(age_group == ag, !is.na(.data[[ex]])), ex,
                   paste("Age:", ag))
    if (!is.null(r)) strat_rows[[length(strat_rows) + 1]] <- r
  }
}
strat <- dplyr::bind_rows(strat_rows)
readr::write_csv(strat, file.path(OUT_DIR, "Figure_S7_stratified_PAF_forest_sourcedata.csv"))
cat("\n=== 分层 PAF ===\n"); print(strat, row.names = FALSE)

# 森林图
if (nrow(strat) > 0) {
  p <- ggplot2::ggplot(strat, ggplot2::aes(x = paf, y = stratum, colour = exposure)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(0.5)) +
    ggplot2::facet_wrap(~ exposure, scales = "free_y") +
    ggplot2::labs(x = "Local PAF (%)", y = "", colour = "") +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(OUT_DIR, "Figure_S7_stratified_PAF_forest.png"),
                  p, width = 7, height = 5, dpi = 300)
}

cat("\n输出目录：", OUT_DIR, "\n"); print(list.files(OUT_DIR))
cat("\n稿件口径提醒：Table 3 / Figure 8 的主分析 PAF 由**未校正 OR** 算得",
    "（2×2 表 a/b/c/d 与 PAF分析结果.xlsx 一致）；",
    "校正 OR 及其 PAF 见 Table3_sensitivity_adjusted_OR_PAF.csv，",
    "在 2.5 与 3.8 中作为敏感性分析报告。\n")
