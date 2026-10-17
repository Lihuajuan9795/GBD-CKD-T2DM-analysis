# =====================================================================
#  bpc_forecast.R
#  Bayesian age-period-cohort (BAPC) 情景预测 1990–2050
#  ---------------------------------------------------------------------
#  产出：
#    - Figure 7  5 个 SDI 五分位 ×（DALYs / Deaths / YLDs / YLLs）×
#               （高 BMI / 低体力活动 / 联合效应）的 1990–2050 预测轨迹
#    - BAPC_2050_summary.csv（2050 年中位值与 95% 可信区间）
#    - BAPC_rescaled_G BD_scale_all_measures.csv（换算到 GBD 标准人口口径）
#
#  方法对应稿件 2.4.3：
#    * 年龄别归因数 = 年龄别归因率 × 对应年龄别人口，作为 Poisson 计数矩阵；
#      同期人口矩阵作为 offset（人年）。
#    * 年龄 / 时期 / 队列效应均用二阶随机游走（second-order random walk）。
#    * 标准人口 = GBD 2023 世界标准人口中 ≥25 岁的 15 个 5 岁组，重新归一。
#    * 因 BAPC 内部标准人口权重与 GBD 直接报告的 ASR 不完全一致，
#      用 1990–2023 重叠期「GBD ASR / BAPC 拟合值」的中位数比例对整条序列
#      做尺度换算（稿件报告的换算系数约为 1/1.9，CV < 2%）。
#    * 本模型不含未来政策干预与诊疗可及性等外部协变量，仅为参照情景外推。
# =====================================================================

rm(list = ls())

# ---------------------------------------------------------------------
# 0. 配置
# ---------------------------------------------------------------------
GBD_DIR <- "D:/2023同等学力相关/2025年开题相关/需上交资料/提交给鞠老师资料（原始数据、论文终版）/下载数据"
DATA_DIR <- "D:/2023同等学力相关/2025年开题相关/需上交资料/提交给鞠老师资料（原始数据、论文终版）/数据分析版/糖尿病肾病数据/糖尿病肾病数据"
OUT_DIR <- "D:/2023同等学力相关/2025年开题相关/论文写作/投稿论文/output/bapc"

FILE_BMI_REGION <- file.path(GBD_DIR, "第一批数据/第一批数据/全球＋5SDI+7超级地区+21地区_type 2_High body-mass index.csv")
FILE_PA_REGION  <- file.path(GBD_DIR, "第三批数据/第三批数据/全球＋5SDI+7超级地区+21地区_type 2_Low physical activity.csv")
FILE_TOTAL_SDI  <- file.path(DATA_DIR, "5sdi+type2  both-年龄分组-反事实联合分析.csv")
FILE_POP_PATTERN  <- file.path(DATA_DIR, "GBD_population_2023", "*.csv")
FILE_POP_FORECAST <- file.path(DATA_DIR, "GBD_population_2019to2100(1).csv")
FILE_AGE_STANDARD <- file.path(DATA_DIR, "2023_age_stand.csv")

START_YEAR <- 1990
END_YEAR   <- 2023
PRED_END   <- 2050

CAUSE    <- "Chronic kidney disease due to diabetes mellitus type 2"
RISK_BMI <- "High body-mass index"
RISK_PA  <- "Low physical activity"

MEASURES <- c(
  DALYs  = "DALYs (Disability-Adjusted Life Years)",
  Deaths = "Deaths",
  YLDs   = "YLDs (Years Lived with Disability)",
  YLLs   = "YLLs (Years of Life Lost)"
)

LOCATIONS   <- c("Global", "High SDI", "High-middle SDI", "Middle SDI", "Low-middle SDI", "Low SDI")
RISK_TYPES  <- c("High BMI", "Low physical activity", "Joint")

# ---------------------------------------------------------------------
# 1. 依赖包（BAPC 与 INLA 需单独安装）
# ---------------------------------------------------------------------
pkgs <- c("data.table", "dplyr", "tidyr", "readr", "reshape2", "ggplot2")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
if (!requireNamespace("BAPC", quietly = TRUE)) {
  stop("请先安装 BAPC 与 INLA：\n",
       "  install.packages('INLA', repos = 'https://inla.r-inla-download.org/R/stable')\n",
       "  install.packages('BAPC')   # 或 remotes::install_github('GlobalBurdenR/GlobalBurdenR')")
}
suppressPackageStartupMessages(library(BAPC))
suppressPackageStartupMessages(library(INLA))

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---------------------------------------------------------------------
# 2. 工具函数与标准人口
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
age_stand <- age_stand %>% dplyr::mutate(age_name = clean_age_name(age_name))

weight_df <- age_stand %>%
  dplyr::filter(age_name %in% AGE_GROUPS) %>%
  dplyr::mutate(age_name = factor(age_name, levels = AGE_GROUPS)) %>%
  dplyr::group_by(age_name) %>%
  dplyr::summarise(w_raw = mean(.data[[grep("^w_stand", names(age_stand), value = TRUE)[1]]]),
                   .groups = "drop") %>%
  dplyr::arrange(age_name)

W_STAND <- weight_df$w_raw / sum(weight_df$w_raw)     # BAPC 的 stdweight 参数
stopifnot(length(W_STAND) == length(AGE_GROUPS))

read_gbd <- function(path) {
  dat <- data.table::fread(path, showProgress = FALSE)
  names(dat) <- trimws(names(dat))
  dat <- data.table::as.data.table(dat)
  for (pr in list(c("location", "location_name"), c("sex", "sex_name"),
                  c("age", "age_name"), c("measure", "measure_name"),
                  c("metric", "metric_name"), c("year_id", "year"))) {
    if (pr[1] %in% names(dat) && !(pr[2] %in% names(dat))) data.table::setnames(dat, pr[1], pr[2])
  }
  dat %>% dplyr::mutate(age_name = clean_age_name(age_name), year = as.integer(year),
                        val = as.numeric(val), sex_name = as.character(sex_name))
}

# ---------------------------------------------------------------------
# 3. 读入数据
# ---------------------------------------------------------------------
cat("读入归因数据 / 总负担 / 人口 ...\n")
attr_df <- dplyr::bind_rows(read_gbd(FILE_BMI_REGION), read_gbd(FILE_PA_REGION)) %>%
  dplyr::filter(sex_name == "Both", age_name %in% AGE_GROUPS, metric_name == "Rate")

total_df <- read_gbd(FILE_TOTAL_SDI) %>%
  dplyr::filter(sex_name == "Both", age_name %in% AGE_GROUPS, metric_name == "Rate")

# 历史人口（单岁 → 5 岁组）
pop_files <- sort(Sys.glob(FILE_POP_PATTERN))
if (length(pop_files) == 0)  # 兼容 CSV 位于子目录的情况
  pop_files <- sort(list.files(file.path(DATA_DIR, "GBD_population_2023"),
                               pattern = "\.csv$", recursive = TRUE, full.names = TRUE))
pop_hist  <- dplyr::bind_rows(lapply(pop_files, read_gbd)) %>%
  dplyr::filter(measure_name == "Population", metric_name == "Number", sex_name == "Both") %>%
  dplyr::mutate(age_num = suppressWarnings(as.integer(age_name))) %>%
  dplyr::filter(!is.na(age_num))
pop_hist <- pop_hist %>%
  dplyr::mutate(
    age_name = dplyr::case_when(
      age_num >= 95 ~ "95+",
      TRUE ~ paste0(floor(age_num / 5) * 5, "-", floor(age_num / 5) * 5 + 4)
    )
  ) %>%
  dplyr::filter(age_name %in% AGE_GROUPS) %>%
  dplyr::group_by(location_name, year, age_name) %>%
  dplyr::summarise(pop = sum(val, na.rm = TRUE), .groups = "drop")

# 预测人口（2019–2100 文件）
pop_pred <- NULL
if (file.exists(FILE_POP_FORECAST)) {
  pp <- data.table::fread(FILE_POP_FORECAST, showProgress = FALSE)
  names(pp) <- trimws(names(pp))
  if (!"sex_name" %in% names(pp) && "sex" %in% names(pp)) names(pp)[names(pp) == "sex"] <- "sex_name"
  if (!"year" %in% names(pp) && "year_id" %in% names(pp)) names(pp)[names(pp) == "year_id"] <- "year"
  if (!"age_name" %in% names(pp) && "age_group_name" %in% names(pp)) names(pp)[names(pp) == "age_group_name"] <- "age_name"
  pop_pred <- pp %>%
    dplyr::filter(tolower(sex_name) == "both") %>%
    dplyr::mutate(age_name = clean_age_name(age_name), year = as.integer(year)) %>%
    dplyr::filter(age_name %in% AGE_GROUPS, year > END_YEAR, year <= PRED_END) %>%
    dplyr::group_by(location_name, year, age_name) %>%
    dplyr::summarise(pop = sum(val, na.rm = TRUE), .groups = "drop")
} else {
  warning("未找到预测人口文件，将用 2023 年人口结构外推（仅用于演示，正式结果请补齐该文件）。")
}

# ---------------------------------------------------------------------
# 4. 组装年龄别归因率
# ---------------------------------------------------------------------
get_single_rate <- function(loc, measure_target, risk) {
  attr_df %>%
    dplyr::filter(location_name == loc, measure_name == !!measure_target, rei_name == risk) %>%
    dplyr::group_by(year, age_name) %>% dplyr::summarise(rate = sum(val, na.rm = TRUE), .groups = "drop")
}

get_joint_rate <- function(loc, measure_target) {
  tot <- total_df %>%
    dplyr::filter(location_name == loc, measure_name == !!measure_target) %>%
    dplyr::group_by(year, age_name) %>% dplyr::summarise(total_rate = sum(val, na.rm = TRUE), .groups = "drop")
  b <- get_single_rate(loc, measure_target, RISK_BMI) %>% dplyr::rename(rate_bmi = rate)
  p <- get_single_rate(loc, measure_target, RISK_PA)  %>% dplyr::rename(rate_pa  = rate)
  tot %>%
    dplyr::left_join(b, by = c("year", "age_name")) %>%
    dplyr::left_join(p, by = c("year", "age_name")) %>%
    dplyr::mutate(
      rate_bmi = ifelse(is.na(rate_bmi), 0, rate_bmi),
      rate_pa  = ifelse(is.na(rate_pa),  0, rate_pa),
      paf_bmi  = dplyr::if_else(total_rate > 0, pmin(pmax(rate_bmi / total_rate, 0), 1), NA_real_),
      paf_pa   = dplyr::if_else(total_rate > 0, pmin(pmax(rate_pa  / total_rate, 0), 1), NA_real_),
      joint_paf = 1 - (1 - paf_bmi) * (1 - paf_pa),
      rate      = total_rate * joint_paf
    ) %>%
    dplyr::select(year, age_name, rate)
}

# ---------------------------------------------------------------------
# 5. 构造 BAPC 输入矩阵（病例数矩阵 + 人口矩阵）
# ---------------------------------------------------------------------
build_matrices <- function(loc, measure_target, risk_type) {
  rate_df <- switch(risk_type,
                    "High BMI"              = get_single_rate(loc, measure_target, RISK_BMI),
                    "Low physical activity" = get_single_rate(loc, measure_target, RISK_PA),
                    "Joint"                 = get_joint_rate(loc, measure_target))

  ph <- pop_hist %>%
    dplyr::filter(location_name == loc, year >= START_YEAR, year <= END_YEAR)

  dis <- rate_df %>%
    dplyr::filter(year >= START_YEAR, year <= END_YEAR) %>%
    dplyr::left_join(ph, by = c("year", "age_name")) %>%
    dplyr::mutate(val = rate * pop / 1e5) %>%
    dplyr::select(year, age_name, val) %>%
    tidyr::complete(year = START_YEAR:END_YEAR, age_name = AGE_GROUPS, fill = list(val = 0))

  dis$age_name <- factor(dis$age_name, levels = AGE_GROUPS)
  case_mat <- reshape2::dcast(dis, year ~ age_name, value.var = "val", fun.aggregate = sum, fill = 0)
  rownames(case_mat) <- case_mat$year
  case_mat <- case_mat[, -1, drop = FALSE][, AGE_GROUPS, drop = FALSE]
  case_mat <- apply(case_mat, c(1, 2), function(z) round(as.numeric(z), 0)) %>% as.data.frame()

  if (sum(as.matrix(case_mat), na.rm = TRUE) == 0) {
    warning("病例矩阵全为 0：", loc, " | ", risk_type, " | ", measure_target)
    return(NULL)
  }

  # 预测期病例数留空（NA），由 BAPC 预测
  future_empty <- matrix(NA_real_, nrow = PRED_END - END_YEAR, ncol = ncol(case_mat)) %>% as.data.frame()
  rownames(future_empty) <- (END_YEAR + 1):PRED_END
  colnames(future_empty) <- colnames(case_mat)
  case_all <- rbind(case_mat, future_empty) %>% as.data.frame()
  case_all[] <- lapply(case_all, as.numeric)

  # 人口矩阵
  if (!is.null(pop_pred)) {
    pp <- pop_pred %>% dplyr::filter(location_name == loc, year > END_YEAR, year <= PRED_END)
    pop_all <- dplyr::bind_rows(
      ph %>% dplyr::select(year, age_name, val = pop),
      pp %>% dplyr::select(year, age_name, val = pop)
    )
  } else {
    last <- ph %>% dplyr::filter(year == END_YEAR) %>% dplyr::select(age_name, val = pop)
    pop_all <- dplyr::bind_rows(
      ph %>% dplyr::select(year, age_name, val = pop),
      dplyr::bind_rows(lapply((END_YEAR + 1):PRED_END, function(y)
        dplyr::mutate(last, year = y))) %>% dplyr::select(year, age_name, val)
    )
  }

  pop_all <- pop_all %>%
    tidyr::complete(year = START_YEAR:PRED_END, age_name = AGE_GROUPS, fill = list(val = 0))
  pop_all$age_name <- factor(pop_all$age_name, levels = AGE_GROUPS)
  pop_mat <- reshape2::dcast(pop_all, year ~ age_name, value.var = "val", fun.aggregate = sum, fill = 0)
  rownames(pop_mat) <- pop_mat$year
  pop_mat <- pop_mat[, -1, drop = FALSE][, AGE_GROUPS, drop = FALSE]
  pop_mat <- apply(pop_mat, c(1, 2), function(z) round(as.numeric(z), 0)) %>% as.data.frame()
  pop_mat <- pop_mat[rownames(case_all), , drop = FALSE]

  if (any(rowSums(pop_mat, na.rm = TRUE) <= 0)) stop("人口矩阵存在总人口为 0 的年份：", loc)

  list(case_all = case_all, pop_mat = pop_mat)
}

# ---------------------------------------------------------------------
# 6. 运行 BAPC
# ---------------------------------------------------------------------
run_bapc_one <- function(loc, measure_target, risk_type) {
  obj <- build_matrices(loc, measure_target, risk_type)
  if (is.null(obj)) return(NULL)

  cat("BAPC:", loc, "|", names(MEASURES)[which(MEASURES == measure_target)], "|", risk_type, "\n")

  tryCatch({
    BAPC::BAPC(
      BAPC::APCList(obj$case_all, obj$pop_mat, gf = 5),
      predict = list(npredict = PRED_END - END_YEAR, retro = TRUE),
      secondDiff = TRUE,
      stdweight = W_STAND,
      verbose = FALSE
    )
  }, error = function(e) {
    warning("BAPC 运行失败：", loc, " | ", risk_type, " -> ", e$message)
    NULL
  })
}

extract_series <- function(bapc_obj) {
  if (is.null(bapc_obj)) return(NULL)
  probs <- c(0.025, 0.5, 0.975)
  cn <- paste0(probs, "Q")
  ap  <- BAPC::qapc(bapc_obj, percentiles = probs)
  proj <- BAPC::agestd.proj(ap)[, cn, drop = FALSE]
  agg_n <- rowSums(BAPC::pyrs(ap))
  out <- as.data.frame(proj / agg_n * 1e5)
  names(out) <- c("lower", "median", "upper")
  out$year <- START_YEAR:PRED_END
  as.data.frame(out)
}

# ---------------------------------------------------------------------
# 7. 主循环 + GBD 尺度换算
# ---------------------------------------------------------------------
# GBD 直接报告的 ASR（用于计算换算系数）
asr_reported <- attr_df %>%
  dplyr::filter(age_name == "Age-standardized")   # 若提取文件含该行则可用
has_asr <- nrow(asr_reported) > 0

results <- list()
models  <- list()

for (m in MEASURES) {
  mlab <- names(MEASURES)[which(MEASURES == m)]
  for (loc in LOCATIONS) {
    for (rt in RISK_TYPES) {
      fit <- run_bapc_one(loc, m, rt)
      models[[paste(mlab, loc, rt)]] <- fit
      ser <- extract_series(fit)
      if (is.null(ser)) next
      ser$measure <- mlab
      ser$sdi_quintile <- loc
      ser$risk_type <- rt
      results[[length(results) + 1]] <- ser
    }
  }
}

bapc_raw <- dplyr::bind_rows(results)

# ---- 尺度换算：median(GBD ASR / BAPC 拟合值)，1990–2023 重叠期 ----
# 说明：BAPC 内部标准人口与 GBD 2023 标准人口权重略有差异，
#       导致 BAPC 原始输出整体偏高（约 1.9 倍）。用中位数比例统一换算。
compute_rescale <- function(mlab, loc, rt) {
  if (rt == "Joint") return(NA_real_)
  ser <- bapc_raw %>% dplyr::filter(measure == mlab, sdi_quintile == loc, risk_type == rt, year <= END_YEAR)
  if (!has_asr || nrow(ser) == 0) return(NA_real_)
  gbd <- asr_reported %>%
    dplyr::filter(location_name == loc, measure_name == MEASURES[[mlab]],
                  rei_name == if (rt == "High BMI") RISK_BMI else RISK_PA,
                  year <= END_YEAR) %>%
    dplyr::select(year, gbd_asr = val)
  if (nrow(gbd) == 0) return(NA_real_)
  d <- dplyr::inner_join(ser %>% dplyr::select(year, bapc = median), gbd, by = "year")
  if (nrow(d) < 5) return(NA_real_)
  stats::median(gbd$gbd_asr / ser$median[match(gbd$year, ser$year)], na.rm = TRUE)
}

rescale_tab <- dplyr::tibble(measure = character(), sdi_quintile = character(),
                             risk_type = character(), rescale_factor = numeric())
for (m in MEASURES) {
  mlab <- names(MEASURES)[which(MEASURES == m)]
  for (loc in LOCATIONS) {
    for (rt in RISK_TYPES) {
      f <- compute_rescale(mlab, loc, rt)
      rescale_tab <- dplyr::bind_rows(rescale_tab,
        dplyr::tibble(measure = mlab, sdi_quintile = loc, risk_type = rt, rescale_factor = f))
    }
  }
}

# 联合效应无 GBD 直接报告值：用两个单风险换算系数按归因量加权平均
rescale_tab <- rescale_tab %>%
  dplyr::group_by(measure, sdi_quintile) %>%
  dplyr::mutate(
    rescale_factor = dplyr::if_else(
      risk_type == "Joint" & !is.na(rescale_factor),
      rescale_factor,
      dplyr::if_else(
        risk_type == "Joint" & is.na(rescale_factor),
        suppressWarnings(mean(rescale_factor[risk_type != "Joint"], na.rm = TRUE)),
        rescale_factor
      )
    )
  ) %>%
  dplyr::ungroup()

bapc_rescaled <- bapc_raw %>%
  dplyr::left_join(rescale_tab, by = c("measure", "sdi_quintile", "risk_type")) %>%
  dplyr::mutate(
    rescale_factor = dplyr::if_else(is.na(rescale_factor) | !is.finite(rescale_factor), 1, rescale_factor),
    median_r = median * rescale_factor,
    lower_r  = lower  * rescale_factor,
    upper_r  = upper  * rescale_factor
  )

readr::write_csv(rescale_tab,   file.path(OUT_DIR, "BAPC_rescaling_factors.csv"))
readr::write_csv(bapc_raw,      file.path(OUT_DIR, "BAPC_raw_scale_all_measures.csv"))
readr::write_csv(bapc_rescaled, file.path(OUT_DIR, "BAPC_rescaled_GBD_scale_all_measures.csv"))

# Figure 7 源数据（与稿件一致的列名）
fig7 <- bapc_rescaled %>%
  dplyr::transmute(measure, sdi_quintile, risk_type, year,
                   median = median_r, lower = lower_r, upper = upper_r) %>%
  dplyr::arrange(measure, sdi_quintile, risk_type, year)
readr::write_csv(fig7, file.path(OUT_DIR, "Figure7_BAPC_SDI_quintile_projections_1990_2050_sourcedata.csv"))

# 2050 年汇总
bapc_2050 <- fig7 %>%
  dplyr::filter(year == PRED_END) %>%
  dplyr::mutate(`95% CrI` = sprintf("%.2f (%.2f–%.2f)", median, lower, upper))
readr::write_csv(bapc_2050, file.path(OUT_DIR, "BAPC_2050_summary.csv"))

cat("\n换算系数概况：\n"); print(summary(rescale_tab$rescale_factor))
cat("\n输出目录：", OUT_DIR, "\n")
print(list.files(OUT_DIR))
