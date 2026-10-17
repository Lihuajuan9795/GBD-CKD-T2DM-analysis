# =====================================================================
#  mapping_global.R
#  2023 年全球国家层面年龄标化归因率（ASR）choropleth 地图
#  ---------------------------------------------------------------------
#  产出：
#    - Figure 5  四项指标 ×（高 BMI / 低体力活动）的 2023 年世界地图
#      （png / pdf，以及对应的 _sourcedata.csv）
#
#  方法对应稿件 2.4.5：
#    * 使用 GBD 2023 的国家边界绘制；
#    * 配色采用等间距或分位数分级；
#    * 为小岛屿国家与人口稠密地区加插图以提高可读性。
#
#  ⚠ 地图合规要求（务必遵守）：
#    1. 底图必须使用符合国家标准的边界数据，建议从
#       「自然资源部标准地图服务系统」（https://bzdt.ch.mnr.gov.cn/）
#       下载 1:100 万世界国家地图（含南海诸岛插图版）后指定 SHAPEFILE 路径；
#       或使用 GBD 官方发布的、已按我国标准处理的边界文件。
#    2. 中国台湾、中国香港、中国澳门必须作为中国的一部分绘制，
#       不得单独着色或单独标注为国家。本脚本在读入后强制执行合并
#       （见 enforce_map_compliance()）。
#    3. 必须保留南海诸岛插图（脚本默认输出右下角南海插图）。
#    4. 不得使用未经审核的境外在线底图服务（Google / OSM / Mapbox 等）。
# =====================================================================

rm(list = ls())

# ---------------------------------------------------------------------
# 0. 配置
# ---------------------------------------------------------------------
GBD_DIR <- "D:/2023同等学力相关/2025年开题相关/需上交资料/提交给鞠老师资料（原始数据、论文终版）/下载数据"
DATA_DIR <- "D:/2023同等学力相关/2025年开题相关/需上交资料/提交给鞠老师资料（原始数据、论文终版）/数据分析版/糖尿病肾病数据/糖尿病肾病数据"
OUT_DIR <- "D:/2023同等学力相关/2025年开题相关/论文写作/投稿论文/output/maps"

FILE_BMI_CTRY <- file.path(GBD_DIR, "第一批数据/第一批数据/204国家_type 2_High body-mass index.csv")
FILE_PA_CTRY  <- file.path(GBD_DIR, "第三批数据/第三批数据/204国家_type 2_Low physical activity.csv")

# 优先使用合规底图；没有则退回 GBD_maps.RData（仅内部核对用，不建议投稿）
SHAPEFILE <- ""                                    # <- 填入合规世界国家边界 shp/geojson 路径
GBD_MAPS_RDATA <- file.path(DATA_DIR, "GBD_maps.RData")

MAP_YEAR   <- 2023
CAUSE      <- "Chronic kidney disease due to diabetes mellitus type 2"
RISK_BMI   <- "High body-mass index"
RISK_PA    <- "Low physical activity"

MEASURES <- c(
  DALYs  = "DALYs (Disability-Adjusted Life Years)",
  Deaths = "Deaths",
  YLDs   = "YLDs (Years Lived with Disability)",
  YLLs   = "YLLs (Years of Life Lost)"
)

N_BREAKS      <- 7
BREAK_METHOD  <- "quantile"    # "quantile" 或 "equal"

# ---------------------------------------------------------------------
# 1. 依赖包
# ---------------------------------------------------------------------
pkgs <- c("data.table", "dplyr", "tidyr", "readr", "stringr", "ggplot2", "sf", "scales", "patchwork")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---------------------------------------------------------------------
# 2. 读入底图并执行合规处理
# ---------------------------------------------------------------------
load_world <- function() {
  if (nzchar(SHAPEFILE) && file.exists(SHAPEFILE)) {
    world <- sf::st_read(SHAPEFILE, quiet = TRUE)
  } else if (file.exists(GBD_MAPS_RDATA)) {
    warning("未指定合规底图 SHAPEFILE，退回 GBD_maps.RData；投稿前请替换为国家标准地图数据。")
    env <- new.env(); load(GBD_MAPS_RDATA, envir = env)
    if (!exists("world_GBD", envir = env)) stop("GBD_maps.RData 中没有 world_GBD 对象。")
    world <- sf::st_as_sf(get("world_GBD", envir = env))
  } else {
    stop("请指定 SHAPEFILE（自然资源部标准地图）或提供 GBD_maps.RData。")
  }
  world
}

enforce_map_compliance <- function(world) {
  nm <- names(world)

  # 名称字段（不同数据源字段名不同）
  name_col <- nm[grepl("name|Name|NAME|admin|ADMIN", nm)][1]
  iso_col  <- nm[grepl("iso|ISO|id|ID|GBD|LOC", nm)][1]

  # --- 合规要点 1：中国台湾 / 中国香港 / 中国澳门归入中国 ---
  tw_keys <- c("TWN", "TAIWAN", "Taiwan", "台湾", "TW")
  hk_keys <- c("HKG", "HONG KONG", "Hong Kong", "香港")
  mo_keys <- c("MAC", "MACAO", "Macao", "Macau", "澳门")

  if (!is.na(name_col)) {
    vals <- as.character(world[[name_col]])
    vals[vals %in% c(tw_keys, hk_keys, mo_keys)] <- "China"
    world[[name_col]] <- vals
  }
  if (!is.na(iso_col)) {
    v2 <- as.character(world[[iso_col]])
    v2[v2 %in% c("TWN", "TW", "HKG", "HK", "MAC", "MO")] <- "CHN"
    world[[iso_col]] <- v2
  }

  # 若存在独立要素，按 China 做一次几何 dissolve（合并为一个要素）
  if (!is.na(name_col)) {
    idx <- which(as.character(world[[name_col]]) == "China")
    if (length(idx) > 1) {
      geo <- sf::st_union(sf::st_geometry(world[idx, ]))
      world <- world[-idx, ]
      cn <- world[idx[1], ]
      sf::st_geometry(cn) <- geo
      world <- rbind(world, cn)
    }
  }

  world
}

world <- enforce_map_compliance(load_world())

# --- 合规自检：不应再有独立的中国台湾要素 ---
name_col_chk <- names(world)[grepl("name|Name|NAME|admin|ADMIN", names(world))][1]
if (!is.na(name_col_chk)) {
  bad <- as.character(world[[name_col_chk]]) %in% c("Taiwan", "TAIWAN", "台湾", "Hong Kong", "Macao")
  if (any(bad)) warning("底图仍存在单独的中国台湾/香港/澳门要素，请更换为合规底图后再出图！")
}

# ---------------------------------------------------------------------
# 3. 读入 2023 年国家 ASR
# ---------------------------------------------------------------------
read_gbd <- function(path) {
  dat <- data.table::fread(path, showProgress = FALSE)
  data.table::as.data.table(dat) %>%
    dplyr::mutate(year = as.integer(year), val = as.numeric(val))
}

asr_map <- dplyr::bind_rows(read_gbd(FILE_BMI_CTRY), read_gbd(FILE_PA_CTRY)) %>%
  dplyr::filter(sex_name == "Both", age_name == "Age-standardized",
                metric_name == "Rate", cause_name == CAUSE, year == MAP_YEAR) %>%
  dplyr::transmute(
    measure = stringr::str_extract(measure_name, "^[A-Za-z]+"),
    risk = dplyr::if_else(rei_name == RISK_BMI, "High BMI", "Low physical activity"),
    location_name, location_id, ASR = val
  )

# ---------------------------------------------------------------------
# 4. 分级配色
# ---------------------------------------------------------------------
make_bins <- function(x, n = N_BREAKS, method = BREAK_METHOD) {
  x <- x[is.finite(x)]
  br <- if (identical(method, "quantile")) {
    unique(stats::quantile(x, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE))
  } else {
    seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 1)
  }
  br <- unique(as.numeric(br))
  labels <- paste(signif(br[-length(br)], 3), signif(br[-1], 3), sep = " to ")
  list(breaks = br, labels = labels)
}

PAL <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
         "#4292C6", "#2171B5", "#084594")

# ---------------------------------------------------------------------
# 5. 绘图（主图 + 南海诸岛插图）
# ---------------------------------------------------------------------
plot_one_map <- function(dat_sf, title, legend_title) {
  base <- ggplot2::ggplot(dat_sf) +
    ggplot2::geom_sf(ggplot2::aes(fill = bin_label), colour = "grey65", linewidth = 0.05) +
    ggplot2::scale_fill_manual(values = PAL, na.value = "grey90", drop = FALSE,
                               name = legend_title) +
    ggplot2::coord_sf(xlim = c(-180, 208), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "bottom",
      legend.key.size = grid::unit(0.35, "cm"),
      legend.text = ggplot2::element_text(size = 6),
      legend.title = ggplot2::element_text(size = 7, face = "bold")
    )

  # --- 合规要点 2：南海诸岛插图（右下角叠加） ---
  inset <- ggplot2::ggplot(dat_sf) +
    ggplot2::geom_sf(ggplot2::aes(fill = bin_label), colour = "grey65", linewidth = 0.05) +
    ggplot2::scale_fill_manual(values = PAL, na.value = "grey90", drop = FALSE,
                               guide = "none") +
    ggplot2::coord_sf(xlim = c(105, 125), ylim = c(3, 25), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "grey40", fill = NA, linewidth = 0.4))

  # 用 patchwork 把南海插图叠加到主图右下角
  base + patchwork::inset_element(inset, left = 0.74, right = 0.99,
                                  bottom = 0.02, top = 0.40)
}

# ---------------------------------------------------------------------
# 6. 主循环
# ---------------------------------------------------------------------
source_out <- list()
plots <- list()

for (mname in names(MEASURES)) {
  for (rk in c("High BMI", "Low physical activity")) {
    sub <- asr_map %>% dplyr::filter(measure == mname, risk == rk)
    if (nrow(sub) == 0) next

    bb <- make_bins(sub$ASR)
    sub <- sub %>%
      dplyr::mutate(bin_label = cut(ASR, breaks = bb$breaks, labels = bb$labels,
                                    include.lowest = TRUE)) %>%
      dplyr::mutate(bin_label = as.character(bin_label))

    # 与底图连接：优先用 location_id（GBD 边界文件），否则用名称
    join_col <- if ("Location.ID" %in% names(world)) "Location.ID" else
                names(world)[grepl("name|Name|NAME", names(world))][1]

    if (join_col == "Location.ID") {
      map_sf <- world %>% dplyr::left_join(sub %>% dplyr::select(location_id, ASR, bin_label),
                                           by = c("Location.ID" = "location_id"))
    } else {
      map_sf <- world %>% dplyr::left_join(sub %>% dplyr::select(location_name, ASR, bin_label),
                                           by = setNames("location_name", join_col))
    }

    tag <- paste0(mname, "_", gsub(" ", "", rk))
    title <- sprintf("%s attributable to %s, %d (ASR per 100,000)", mname, rk, MAP_YEAR)

    source_out[[tag]] <- sub %>%
      dplyr::mutate(figure = "Figure 5", map_name = location_name) %>%
      dplyr::select(figure, measure, risk, map_name, ASR, bin_label)

    p <- plot_one_map(map_sf, title, paste0("ASR per 100,000\n", mname))
    plots[[tag]] <- p

    ggplot2::ggsave(file.path(OUT_DIR, paste0("Figure5_World_map_", tag, "_", MAP_YEAR, ".png")),
                    p, width = 10, height = 5.6, dpi = 300)
    ggplot2::ggsave(file.path(OUT_DIR, paste0("Figure5_World_map_", tag, "_", MAP_YEAR, ".pdf")),
                    p, width = 10, height = 5.6)
  }
}

readr::write_csv(dplyr::bind_rows(source_out),
                 file.path(OUT_DIR, "Figure5_World_maps_2023_HighBMI_LowPA_sourcedata.csv"))

cat("\n输出目录：", OUT_DIR, "\n"); print(list.files(OUT_DIR))
cat("\n提示：投稿前请确认南海诸岛插图完整、中国台湾已并入中国。\n")
