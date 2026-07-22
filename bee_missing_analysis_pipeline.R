
# ============================================================================
# BeeDisease 누락 분석용
# ----------------------------------------------------------------------------
# 기존 bee_spatial_analysis.R / bee_spatial_analysis2.R 의 로직을 재정리하여
# 아래 누락 분석을 한 번에 수행할 수 있도록 만든 독립 스크립트입니다.
#
# 포함 분석
#   1) 기술통계 + 상관분석 표 (Table sbv4 후보)
#   2) 공간 진단/모형적합도 비교표 (OLS vs CAR)
#   3) 시군구 집계 SBV 발생 지도 (전체 + 종별 facet)
#   4) 최종 CAR 모형의 공간효과 지도
#
# 이번 파일에서 의도적으로 제외한 것
#   - sum_sun / sum_rain / water / forest 상위 구간의 종별 비교
#     (해당 부분은 연구일지 결과와 함께 나중에 일괄 반영)
#
# 실행 전 확인
#   - shapefile(BND_SIGUNGU_PG.shp 등) 경로를 cfg$shp_path 에 지정하세요.
#   - 업로드된 CSV 파일들은 /mnt/data 기준 fallback 경로가 이미 들어 있습니다.
#
# 실행 예시
#   source("bee_missing_analysis_pipeline.R")
#   cfg$shp_path <- "data/raw/BND_SIGUNGU_PG.shp"
#   run_missing_analysis(cfg)
# ============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

# ---------------------------
# 0) 패키지 로드
# ---------------------------
required_pkgs <- c(
  "dplyr", "tidyr", "readr", "stringr", "purrr", "tibble",
  "lubridate", "sf", "spdep", "spatialreg", "ggplot2", "scales", "geosphere"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "다음 패키지 설치가 필요합니다: ",
    paste(missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(lubridate)
  library(sf)
  library(spdep)
  library(spatialreg)
  library(ggplot2)
  library(scales)
  library(geosphere)
})

sf::sf_use_s2(FALSE)

# ---------------------------
# 1) 공용 설정
# ---------------------------
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

pick_existing <- function(paths, required = TRUE) {
  paths <- unique(paths)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hits <- paths[file.exists(paths)]
  if (length(hits) > 0) return(normalizePath(hits[1], winslash = "/", mustWork = FALSE))
  if (required) {
    stop("파일을 찾지 못했습니다. 후보 경로: ", paste(paths, collapse = " | "), call. = FALSE)
  }
  NA_character_
}

cfg <- list(
  # 입력 파일
  cerana_path = pick_existing(c(
    file.path(getwd(), "6_Apis_cerana_with_nearby.csv"),
    file.path(getwd(), "data", "final", "6_Apis_cerana_with_nearby.csv"),
    "/mnt/data/6_Apis_cerana_with_nearby.csv"
  )),
  mellifera_path = pick_existing(c(
    file.path(getwd(), "6_Apis_mellifera_with_nearby.csv"),
    file.path(getwd(), "data", "final", "6_Apis_mellifera_with_nearby.csv"),
    "/mnt/data/6_Apis_mellifera_with_nearby.csv"
  )),
  weather_path = pick_existing(c(
    file.path(getwd(), "weather.csv"),
    file.path(getwd(), "data", "raw", "weather.csv"),
    "/mnt/data/weather.csv"
  )),
  station_path = pick_existing(c(
    file.path(getwd(), "weather_station_info.csv"),
    file.path(getwd(), "data", "raw", "weather_station_info.csv"),
    "/mnt/data/weather_station_info.csv"
  )),
  landcover_path = pick_existing(c(
    file.path(getwd(), "환경부_환경공간정보_년도별 대분류토지피복통계 현황_20230901.csv"),
    file.path(getwd(), "data", "raw", "환경부_환경공간정보_년도별 대분류토지피복통계 현황_20230901.csv"),
    "/mnt/data/환경부_환경공간정보_년도별 대분류토지피복통계 현황_20230901.csv"
  )),
  shp_path = pick_existing(c(
    Sys.getenv("SIGUNGU_SHP", unset = ""),
    file.path(getwd(), "BND_SIGUNGU_PG.shp"),
    file.path(getwd(), "data", "raw", "BND_SIGUNGU_PG.shp"),
    file.path(getwd(), "shape", "BND_SIGUNGU_PG.shp"),
    file.path(getwd(), "shp", "BND_SIGUNGU_PG.shp"),
    "/mnt/data/BND_SIGUNGU_PG.shp"
  ), required = FALSE),

  # 분석 옵션
  weather_year = 2018,
  response_transform = "raw",   # "raw" or "log1p"
  include_grass_in_forest = TRUE,
  include_wetland_in_water = TRUE,

  # 출력 폴더
  out_dir = file.path(getwd(), "out_missing_analysis")
)

# ---------------------------
# 2) 입출력 헬퍼
# ---------------------------
safe_read_csv <- function(path, encodings = c("CP949", "UTF-8", "EUC-KR")) {
  last_err <- NULL
  for (enc in encodings) {
    out <- try(
      readr::read_csv(path, locale = locale(encoding = enc), show_col_types = FALSE),
      silent = TRUE
    )
    if (!inherits(out, "try-error")) return(as_tibble(out))
    last_err <- out
  }
  stop("CSV 로드 실패: ", path, "\n", as.character(last_err), call. = FALSE)
}

norm_text <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "\u3000", "")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x
}

canon_sido <- function(x) {
  x <- norm_text(x)
  dplyr::case_when(
    str_detect(x, "서울") ~ "서울",
    str_detect(x, "부산") ~ "부산",
    str_detect(x, "대구") ~ "대구",
    str_detect(x, "인천") ~ "인천",
    str_detect(x, "광주") ~ "광주",
    str_detect(x, "대전") ~ "대전",
    str_detect(x, "울산") ~ "울산",
    str_detect(x, "세종") ~ "세종",
    str_detect(x, "경기") ~ "경기",
    str_detect(x, "강원") ~ "강원",
    str_detect(x, "충청북|충북") ~ "충북",
    str_detect(x, "충청남|충남") ~ "충남",
    str_detect(x, "전라북|전북") ~ "전북",
    str_detect(x, "전라남|전남") ~ "전남",
    str_detect(x, "경상북|경북") ~ "경북",
    str_detect(x, "경상남|경남") ~ "경남",
    str_detect(x, "제주") ~ "제주",
    TRUE ~ x
  )
}

canon_sigungu <- function(x) {
  x <- norm_text(x)
  x <- stringr::str_replace_all(x, "\\(.*?\\)", "")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x
}

sigungu_root <- function(x) {
  x <- canon_sigungu(x)
  out <- ifelse(
    str_detect(x, "^.+시\\s+.+구$"),
    str_replace(x, "^(.*?시)\\s+.*$", "\\1"),
    x
  )
  out <- ifelse(
    str_detect(out, "^.+시.+구$") & !str_detect(out, "\\s"),
    str_replace(out, "^(.*?시).*$", "\\1"),
    out
  )
  out
}

code_to_sido <- function(code) {
  code <- substr(as.character(code), 1, 2)
  dplyr::case_when(
    code == "11" ~ "서울",
    code == "26" ~ "부산",
    code == "27" ~ "대구",
    code == "28" ~ "인천",
    code == "29" ~ "광주",
    code == "30" ~ "대전",
    code == "31" ~ "울산",
    code == "36" ~ "세종",
    code == "41" ~ "경기",
    code == "42" ~ "강원",
    code == "43" ~ "충북",
    code == "44" ~ "충남",
    code == "45" ~ "전북",
    code == "46" ~ "전남",
    code == "47" ~ "경북",
    code == "48" ~ "경남",
    code == "50" ~ "제주",
    TRUE ~ NA_character_
  )
}

# ---------------------------
# 3) shapefile 로드
# ---------------------------
detect_sigungu_cols <- function(sig_sf) {
  nms <- names(sig_sf)
  code_col <- nms[str_detect(nms, regex("^(SIG_CD|SIGUNGU_CD|SGG_CD|ADM_CD|LAWD_CD)$", ignore_case = TRUE))]
  name_col <- nms[str_detect(nms, regex("^(SIG_KOR_NM|SIG_NM|SIGUNGU_NM|SGG_NM)$", ignore_case = TRUE))]
  sido_col <- nms[str_detect(nms, regex("^(CTP_KOR_NM|SIDO_NM|PROV_NM|SIDO)$", ignore_case = TRUE))]

  list(
    code_col = code_col[1] %||% NA_character_,
    name_col = name_col[1] %||% NA_character_,
    sido_col = sido_col[1] %||% NA_character_
  )
}

load_sigungu_sf <- function(shp_path) {
  if (is.na(shp_path) || !file.exists(shp_path)) {
    stop(
      "시군구 shapefile 을 찾지 못했습니다. cfg$shp_path 또는 SIGUNGU_SHP 환경변수를 지정하세요.",
      call. = FALSE
    )
  }

  sig_sf <- sf::st_read(shp_path, options = "ENCODING=CP949", quiet = TRUE)
  if (is.na(st_crs(sig_sf))) stop("시군구 shapefile 에 CRS 정보가 없습니다.", call. = FALSE)
  if (!all(st_is_valid(sig_sf))) sig_sf <- st_make_valid(sig_sf)

  cols <- detect_sigungu_cols(sig_sf)
  if (is.na(cols$code_col) || is.na(cols$name_col)) {
    stop(
      "시군구 shapefile 에서 코드/이름 컬럼을 찾지 못했습니다.\n",
      "names(sig_sf): ", paste(names(sig_sf), collapse = ", "),
      call. = FALSE
    )
  }

  sig_sf <- sig_sf %>%
    mutate(
      SIG_CD = as.character(.data[[cols$code_col]]),
      SIGUNGU_NM = canon_sigungu(.data[[cols$name_col]]),
      SIDO_NM = if (!is.na(cols$sido_col) && cols$sido_col %in% names(sig_sf)) {
        canon_sido(.data[[cols$sido_col]])
      } else {
        code_to_sido(.data[[cols$code_col]])
      },
      SIGUNGU_ROOT = sigungu_root(SIGUNGU_NM),
      key_exact = paste(SIDO_NM, SIGUNGU_NM, sep = "::"),
      key_root = paste(SIDO_NM, SIGUNGU_ROOT, sep = "::")
    ) %>%
    select(SIG_CD, SIDO_NM, SIGUNGU_NM, SIGUNGU_ROOT, key_exact, key_root, geometry)

  cent <- st_point_on_surface(sig_sf)
  cent_ll <- st_transform(cent, 4326)
  xy <- st_coordinates(cent_ll)
  sig_sf$centroid_lon <- xy[, 1]
  sig_sf$centroid_lat <- xy[, 2]

  sig_sf
}

# ---------------------------
# 4) 날씨 요약 (기존 코드 정의 유지)
#   - 2018년 6~8월
#   - mean_temp : 평균기온 평균
#   - sum_rain  : 일최다강수량 평균
#   - sum_sun   : 합계일조시간 평균
#   - ratio_sun : 일조율 평균
# ---------------------------
load_summer_weather <- function(weather_path, station_path, year_keep = 2018) {
  wx <- safe_read_csv(weather_path)
  stn <- safe_read_csv(station_path)

  required_wx <- c("지점", "일시", "평균기온(°C)", "일최다강수량(mm)", "합계일조시간(hr)", "일조율")
  required_stn <- c("지점", "위도", "경도")
  miss_wx <- setdiff(required_wx, names(wx))
  miss_stn <- setdiff(required_stn, names(stn))
  if (length(miss_wx) > 0) stop("weather.csv 필수 컬럼 누락: ", paste(miss_wx, collapse = ", "), call. = FALSE)
  if (length(miss_stn) > 0) stop("weather_station_info.csv 필수 컬럼 누락: ", paste(miss_stn, collapse = ", "), call. = FALSE)

  wx2 <- wx %>%
    mutate(
      일시 = suppressWarnings(lubridate::ymd(`일시`)),
      year = year(일시),
      month = month(일시)
    ) %>%
    filter(!is.na(일시), year == year_keep, month %in% c(6, 7, 8)) %>%
    group_by(`지점`) %>%
    summarise(
      mean_temp = mean(`평균기온(°C)`, na.rm = TRUE),
      sum_rain = mean(`일최다강수량(mm)`, na.rm = TRUE),
      sum_sun = mean(`합계일조시간(hr)`, na.rm = TRUE),
      ratio_sun = mean(`일조율`, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(station_id = `지점`)

  stn2 <- stn %>%
    transmute(
      station_id = `지점`,
      station_lat = as.numeric(`위도`),
      station_lon = as.numeric(`경도`)
    )

  list(weather_summer = wx2, station_info = stn2)
}

nearest_station_id <- function(lat, lon, stations_df) {
  if (nrow(stations_df) == 0) return(NA_integer_)
  d <- geosphere::distHaversine(
    matrix(c(lon, lat), ncol = 2),
    matrix(c(stations_df$station_lon, stations_df$station_lat), ncol = 2)
  )
  stations_df$station_id[which.min(d)]
}

attach_climate_to_sigungu <- function(sig_sf, weather_summer, station_info) {
  sig_env <- sig_sf %>%
    st_drop_geometry() %>%
    transmute(SIG_CD, centroid_lat, centroid_lon) %>%
    mutate(
      nearest_station = purrr::map2_dbl(
        centroid_lat, centroid_lon,
        ~ nearest_station_id(.x, .y, station_info)
      )
    ) %>%
    left_join(weather_summer, by = c("nearest_station" = "station_id")) %>%
    select(SIG_CD, nearest_station, mean_temp, sum_rain, sum_sun, ratio_sun)

  left_join(sig_sf, sig_env, by = "SIG_CD")
}

# ---------------------------
# 5) 토지피복 로드
#   - forest = 산림 + 초지(옵션)
#   - water  = 수역 + 습지(옵션)
#   - exact / root key 둘 다 만들어서 행정구역 변경에 대응
# ---------------------------
load_landcover <- function(landcover_path,
                           include_grass_in_forest = TRUE,
                           include_wetland_in_water = TRUE) {
  lc <- safe_read_csv(landcover_path)

  req <- c("자료년도", "국가구분", "시도", "시군구",
           "시가화건조지역", "농업지역", "산림지역",
           "초지", "습지", "수역", "합계")
  miss <- setdiff(req, names(lc))
  if (length(miss) > 0) stop("landcover 필수 컬럼 누락: ", paste(miss, collapse = ", "), call. = FALSE)

  lc2 <- lc %>%
    filter(`국가구분` %in% c("대한민국", "남한")) %>%
    filter(`자료년도` == max(`자료년도`, na.rm = TRUE)) %>%
    mutate(
      SIDO_NM = canon_sido(`시도`),
      SIGUNGU_NM = canon_sigungu(`시군구`),
      SIGUNGU_ROOT = sigungu_root(SIGUNGU_NM),
      total_area = as.numeric(`합계`),
      forest_raw = as.numeric(`산림지역`) + if (include_grass_in_forest) as.numeric(`초지`) else 0,
      agri_raw = as.numeric(`농업지역`),
      urban_raw = as.numeric(`시가화건조지역`),
      water_raw = as.numeric(`수역`) + if (include_wetland_in_water) as.numeric(`습지`) else 0,
      forest = forest_raw / total_area,
      agri = agri_raw / total_area,
      urban = urban_raw / total_area,
      water = water_raw / total_area,
      key_exact = paste(SIDO_NM, SIGUNGU_NM, sep = "::"),
      key_root = paste(SIDO_NM, SIGUNGU_ROOT, sep = "::")
    ) %>%
    select(SIDO_NM, SIGUNGU_NM, SIGUNGU_ROOT, key_exact, key_root, forest, agri, urban, water)

  exact_tbl <- lc2 %>%
    group_by(key_exact) %>%
    summarise(
      forest = mean(forest, na.rm = TRUE),
      agri = mean(agri, na.rm = TRUE),
      urban = mean(urban, na.rm = TRUE),
      water = mean(water, na.rm = TRUE),
      .groups = "drop"
    )

  root_tbl <- lc2 %>%
    group_by(key_root) %>%
    summarise(
      forest_root = mean(forest, na.rm = TRUE),
      agri_root = mean(agri, na.rm = TRUE),
      urban_root = mean(urban, na.rm = TRUE),
      water_root = mean(water, na.rm = TRUE),
      .groups = "drop"
    )

  list(exact = exact_tbl, root = root_tbl)
}

attach_landcover_to_sigungu <- function(sig_sf, lc_list) {
  sig2 <- sig_sf %>%
    left_join(lc_list$exact, by = "key_exact") %>%
    left_join(lc_list$root, by = "key_root") %>%
    mutate(
      forest = coalesce(forest, forest_root),
      agri = coalesce(agri, agri_root),
      urban = coalesce(urban, urban_root),
      water = coalesce(water, water_root)
    ) %>%
    select(-forest_root, -agri_root, -urban_root, -water_root)

  na_rows <- sum(
    !stats::complete.cases(st_drop_geometry(sig2)[, c("forest", "agri", "urban", "water")])
  )
  message(sprintf("[INFO] landcover NA 행 수: %d / %d", na_rows, nrow(sig2)))
  sig2
}

# ---------------------------
# 6) 벌 데이터 로드 및 시군구 집계
# ---------------------------
load_bee_points <- function(path, species_label) {
  df <- safe_read_csv(path)
  req <- c("Latitude", "Longitude", "발생두수")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) stop("벌 데이터 필수 컬럼 누락: ", paste(miss, collapse = ", "), call. = FALSE)

  df %>%
    transmute(
      Latitude = as.numeric(Latitude),
      Longitude = as.numeric(Longitude),
      발생두수 = as.numeric(`발생두수`),
      species = species_label
    ) %>%
    filter(!is.na(Latitude), !is.na(Longitude), !is.na(발생두수))
}

join_points_to_sigungu <- function(points_df, sig_sf) {
  pts <- st_as_sf(points_df, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
  pts <- st_transform(pts, st_crs(sig_sf))

  joined <- suppressWarnings(
    st_join(pts, sig_sf[, c("SIG_CD", "SIDO_NM", "SIGUNGU_NM")], left = TRUE, join = st_intersects)
  )

  # 경계선/해안선 근처 포인트는 nearest feature 로 보정
  miss_idx <- which(is.na(joined$SIG_CD))
  if (length(miss_idx) > 0) {
    nn <- st_nearest_feature(joined[miss_idx, ], sig_sf)
    joined$SIG_CD[miss_idx] <- sig_sf$SIG_CD[nn]
    joined$SIDO_NM[miss_idx] <- sig_sf$SIDO_NM[nn]
    joined$SIGUNGU_NM[miss_idx] <- sig_sf$SIGUNGU_NM[nn]
  }

  st_drop_geometry(joined)
}

aggregate_point_counts <- function(cerana_df, mellifera_df, sig_sf) {
  all_pts <- bind_rows(cerana_df, mellifera_df)
  pts_sig <- join_points_to_sigungu(all_pts, sig_sf)

  by_species <- pts_sig %>%
    group_by(SIG_CD, species) %>%
    summarise(
      sbv_count = sum(발생두수, na.rm = TRUE),
      apiary_n = n(),
      .groups = "drop"
    )

  by_total <- pts_sig %>%
    group_by(SIG_CD) %>%
    summarise(
      sbv_count = sum(발생두수, na.rm = TRUE),
      apiary_n = n(),
      .groups = "drop"
    )

  list(by_species = by_species, by_total = by_total)
}

make_analysis_dataset <- function(sig_sf, counts_total, counts_species, response_transform = "raw") {
  total_df <- sig_sf %>%
    left_join(counts_total, by = "SIG_CD") %>%
    mutate(
      sbv_count = replace_na(sbv_count, 0),
      apiary_n = replace_na(apiary_n, 0),
      response = if (response_transform == "log1p") log1p(sbv_count) else sbv_count
    ) %>%
    st_drop_geometry()

  species_df <- tidyr::expand_grid(
    SIG_CD = sig_sf$SIG_CD,
    species = c("cerana", "mellifera")
  ) %>%
    left_join(sig_sf %>% st_drop_geometry(), by = "SIG_CD") %>%
    left_join(counts_species, by = c("SIG_CD", "species")) %>%
    mutate(
      sbv_count = replace_na(sbv_count, 0),
      apiary_n = replace_na(apiary_n, 0)
    )

  list(total = total_df, by_species = species_df)
}

# ---------------------------
# 7) 기술통계 + 상관분석 표
# ---------------------------
cor_one_var <- function(df, response, var_name) {
  dat <- df %>%
    select(all_of(c(response, var_name))) %>%
    filter(if_all(everything(), ~ !is.na(.x)))

  if (nrow(dat) < 5) {
    return(tibble(
      variable = var_name,
      n = nrow(dat),
      mean = NA_real_,
      sd = NA_real_,
      min = NA_real_,
      median = NA_real_,
      max = NA_real_,
      pearson_r = NA_real_,
      pearson_p = NA_real_,
      spearman_rho = NA_real_,
      spearman_p = NA_real_
    ))
  }

  x <- dat[[var_name]]
  y <- dat[[response]]
  p1 <- suppressWarnings(cor.test(y, x, method = "pearson"))
  p2 <- suppressWarnings(cor.test(y, x, method = "spearman", exact = FALSE))

  tibble(
    variable = var_name,
    n = nrow(dat),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    pearson_r = unname(p1$estimate),
    pearson_p = p1$p.value,
    spearman_rho = unname(p2$estimate),
    spearman_p = p2$p.value
  )
}

make_desc_corr_table <- function(df, response = "sbv_count") {
  vars <- c("mean_temp", "sum_sun", "ratio_sun", "sum_rain", "forest", "agri", "urban", "water")
  vars <- vars[vars %in% names(df)]
  bind_rows(lapply(vars, function(v) cor_one_var(df, response = response, var_name = v)))
}

# ---------------------------
# 8) LaTeX tabular 저장
# ---------------------------
escape_latex <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\\\", "\\\\textbackslash{}")
  x <- str_replace_all(x, "([#$%&_{}])", "\\\\\\1")
  x <- str_replace_all(x, "~", "\\\\textasciitilde{}")
  x <- str_replace_all(x, "\\^", "\\\\textasciicircum{}")
  x
}

write_latex_table <- function(df, file, caption = NULL, label = NULL, digits = 3) {
  df2 <- df
  for (j in seq_along(df2)) {
    if (is.numeric(df2[[j]])) {
      df2[[j]] <- ifelse(is.na(df2[[j]]), "", formatC(df2[[j]], digits = digits, format = "f"))
    } else {
      df2[[j]] <- ifelse(is.na(df2[[j]]), "", as.character(df2[[j]]))
    }
  }

  align <- paste(c("l", rep("r", ncol(df2) - 1)), collapse = "")
  body_lines <- apply(df2, 1, function(row) paste(escape_latex(row), collapse = " & "))
  body_lines <- paste0(body_lines, " \\\\")

  lines <- c(
    "\\begin{table}[t]",
    if (!is.null(caption)) paste0("\\caption{", escape_latex(caption), "}") else NULL,
    if (!is.null(label)) paste0("\\label{", label, "}") else NULL,
    "\\centering",
    paste0("\\begin{tabular}{", align, "}"),
    "\\hline",
    paste(paste0("\\textbf{", escape_latex(names(df2)), "}"), collapse = " & "),
    "\\\\",
    "\\hline",
    body_lines,
    "\\hline",
    "\\end{tabular}",
    "\\end{table}"
  )

  writeLines(lines, con = file, useBytes = TRUE)
}

write_outputs <- function(df, out_csv, out_tex = NULL, caption = NULL, label = NULL, digits = 3) {
  readr::write_csv(df, out_csv)
  if (!is.null(out_tex)) {
    write_latex_table(df, out_tex, caption = caption, label = label, digits = digits)
  }
}

# ---------------------------
# 9) 공간 모형 적합 및 진단
# ---------------------------
make_nb_objects <- function(sig_sf_subset) {
  nb <- spdep::poly2nb(sig_sf_subset, queen = TRUE)
  list(
    nb = nb,
    listw_W = spdep::nb2listw(nb, style = "W", zero.policy = TRUE),
    listw_B = spdep::nb2listw(nb, style = "B", zero.policy = TRUE)
  )
}

safe_moran <- function(x, listw) {
  if (length(x) < 5) return(tibble(I = NA_real_, p_value = NA_real_))
  mt <- try(spdep::moran.test(x, listw, zero.policy = TRUE), silent = TRUE)
  if (inherits(mt, "try-error")) return(tibble(I = NA_real_, p_value = NA_real_))
  tibble(
    I = unname(mt$estimate[[1]]),
    p_value = mt$p.value
  )
}

make_X_from_model <- function(fit, newdata) {
  b <- coef(fit)
  ff <- formula(fit)
  X <- model.matrix(ff, data = newdata)

  if ("(Intercept)" %in% names(b) && !("(Intercept)" %in% colnames(X))) {
    X <- cbind("(Intercept)" = rep(1, nrow(X)), X)
  }

  miss <- setdiff(names(b), colnames(X))
  if (length(miss) > 0) {
    Z <- matrix(0, nrow = nrow(X), ncol = length(miss))
    colnames(Z) <- miss
    X <- cbind(X, Z)
  }

  X <- X[, names(b), drop = FALSE]
  list(X = X, b = as.numeric(b))
}

extract_spatial_component <- function(fit, data) {
  xb <- make_X_from_model(fit, data)
  trend_fit <- drop(xb$X %*% xb$b)

  response_fit <- try(as.numeric(stats::fitted(fit)), silent = TRUE)
  if (inherits(response_fit, "try-error") || length(response_fit) != nrow(data)) {
    response_fit <- fit$fitted.values %||% try(as.numeric(stats::predict(fit)), silent = TRUE)
  }
  if (inherits(response_fit, "try-error") || length(response_fit) != nrow(data)) {
    stop("CAR fitted values 를 추출하지 못했습니다.", call. = FALSE)
  }

  tibble(
    trend_fit = as.numeric(trend_fit),
    fitted_response = as.numeric(response_fit),
    spatial_effect = as.numeric(response_fit) - as.numeric(trend_fit),
    residual = data[[all.vars(formula(fit))[1]]] - as.numeric(response_fit)
  )
}

fit_ols_car_pair <- function(sig_sf_full, data_full, response, predictors, model_name) {
  needed <- unique(c("SIG_CD", response, predictors))
  dat <- data_full %>%
    select(any_of(needed)) %>%
    filter(if_all(all_of(c(response, predictors)), ~ !is.na(.x)))

  if (nrow(dat) < 20) {
    stop(sprintf("[%s] 유효 관측치가 너무 적습니다: n=%d", model_name, nrow(dat)), call. = FALSE)
  }

  sf_dat <- sig_sf_full %>% semi_join(dat, by = "SIG_CD")
  dat <- dat %>% arrange(match(SIG_CD, sf_dat$SIG_CD))
  sf_dat <- sf_dat %>% arrange(match(SIG_CD, dat$SIG_CD))

  nb_obj <- make_nb_objects(sf_dat)
  fml <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))

  ols_fit <- stats::lm(fml, data = dat)
  car_fit <- spatialreg::spautolm(
    formula = fml,
    data = dat,
    listw = nb_obj$listw_B,
    family = "CAR",
    zero.policy = TRUE
  )
  car_fit$call$formula <- fml

  moran_y <- safe_moran(dat[[response]], nb_obj$listw_W)
  moran_ols <- safe_moran(residuals(ols_fit), nb_obj$listw_W)
  moran_car <- safe_moran(residuals(car_fit), nb_obj$listw_W)

  ll_ols <- as.numeric(logLik(ols_fit))
  ll_car <- as.numeric(logLik(car_fit))
  lr_stat <- 2 * (ll_car - ll_ols)
  lr_p <- pchisq(lr_stat, df = 1, lower.tail = FALSE)

  diag_row <- tibble(
    model = model_name,
    response = response,
    predictors = paste(predictors, collapse = ", "),
    n = nrow(dat),
    moran_y_I = moran_y$I,
    moran_y_p = moran_y$p_value,
    moran_ols_resid_I = moran_ols$I,
    moran_ols_resid_p = moran_ols$p_value,
    moran_car_resid_I = moran_car$I,
    moran_car_resid_p = moran_car$p_value,
    lambda = car_fit$lambda %||% NA_real_,
    lr_stat = lr_stat,
    lr_p = lr_p,
    logLik_ols = ll_ols,
    logLik_car = ll_car,
    AIC_ols = AIC(ols_fit),
    AIC_car = AIC(car_fit)
  )

  list(
    ols = ols_fit,
    car = car_fit,
    data = dat,
    sf = sf_dat,
    nb = nb_obj$nb,
    listw_W = nb_obj$listw_W,
    listw_B = nb_obj$listw_B,
    diagnostics = diag_row
  )
}

# ---------------------------
# 10) 파스텔 스타일
# ---------------------------
pal_species <- c(
  cerana = "#B9D6CC",
  mellifera = "#E6C3B5"
)

pal_seq <- c("#F7F4EE", "#ECE6F2", "#D9E4F2", "#C5D8EA", "#A8C4DC", "#83A8C8")
pal_div_low <- "#BFD5E6"
pal_div_mid <- "#F8F5F0"
pal_div_high <- "#E5B8A8"

base_theme_pastel <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "#E9E7E2", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      plot.title = element_text(face = "bold", size = 14, color = "#4E4A46"),
      plot.subtitle = element_text(size = 10, color = "#6A655F"),
      legend.title = element_text(size = 10, color = "#4E4A46"),
      legend.text = element_text(size = 9, color = "#5A554F"),
      strip.text = element_text(face = "bold", color = "#4E4A46"),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

save_count_map_total <- function(sig_sf_counts, out_file) {
  p <- ggplot(sig_sf_counts) +
    geom_sf(aes(fill = sbv_count), color = "#F7F6F2", linewidth = 0.15) +
    scale_fill_gradientn(
      colours = pal_seq,
      trans = "sqrt",
      labels = scales::comma,
      na.value = "#F3F0EA",
      name = "SBV 발생두수"
    ) +
    labs(
      title = "시군구 집계 SBV 발생 분포",
      subtitle = "전체 발생두수 기준"
    ) +
    coord_sf(datum = NA) +
    base_theme_pastel()

  ggsave(out_file, p, width = 8.5, height = 10, dpi = 320, bg = "white")
  invisible(p)
}

save_count_map_by_species <- function(sig_sf, species_df, out_file) {
  long_sf <- sig_sf %>%
    select(SIG_CD, geometry) %>%
    left_join(species_df, by = "SIG_CD")

  p <- ggplot(long_sf) +
    geom_sf(aes(fill = sbv_count), color = "#F7F6F2", linewidth = 0.12) +
    facet_wrap(~ species, ncol = 2) +
    scale_fill_gradientn(
      colours = pal_seq,
      trans = "sqrt",
      labels = scales::comma,
      na.value = "#F3F0EA",
      name = "SBV 발생두수"
    ) +
    labs(
      title = "시군구 집계 SBV 발생 분포",
      subtitle = "종별 발생두수"
    ) +
    coord_sf(datum = NA) +
    base_theme_pastel()

  ggsave(out_file, p, width = 12, height = 6.8, dpi = 320, bg = "white")
  invisible(p)
}

save_spatial_effect_map <- function(sig_sf_effect, out_file) {
  lim <- max(abs(sig_sf_effect$spatial_effect), na.rm = TRUE)
  if (!is.finite(lim) || is.na(lim)) lim <- 1

  p <- ggplot(sig_sf_effect) +
    geom_sf(aes(fill = spatial_effect), color = "#F7F6F2", linewidth = 0.15) +
    scale_fill_gradient2(
      low = pal_div_low,
      mid = pal_div_mid,
      high = pal_div_high,
      midpoint = 0,
      limits = c(-lim, lim),
      labels = scales::number_format(accuracy = 0.1),
      na.value = "#F3F0EA",
      name = "공간효과"
    ) +
    labs(
      title = "최종 CAR 모형의 공간효과",
      subtitle = "fitted - Xβ"
    ) +
    coord_sf(datum = NA) +
    base_theme_pastel()

  ggsave(out_file, p, width = 8.5, height = 10, dpi = 320, bg = "white")
  invisible(p)
}

# ---------------------------
# 11) QC 요약
# ---------------------------
write_qc_summary <- function(sig_sf, total_df, out_file) {
  qc <- c(
    sprintf("시군구 수: %d", nrow(sig_sf)),
    sprintf("총 발생두수 합계: %s", format(sum(total_df$sbv_count, na.rm = TRUE), big.mark = ",")),
    sprintf("기상 complete rows: %d", sum(complete.cases(total_df[, c("mean_temp", "sum_sun", "sum_rain", "ratio_sun")]))),
    sprintf("토지피복 complete rows: %d", sum(complete.cases(total_df[, c("forest", "agri", "urban", "water")]))),
    sprintf("기후+토지피복 complete rows: %d", sum(complete.cases(total_df[, c("sum_sun", "sum_rain", "water")]))),
    sprintf("response transform: %s", cfg$response_transform)
  )
  writeLines(qc, con = out_file, useBytes = TRUE)
}

# ---------------------------
# 12) 메인 파이프라인
# ---------------------------
run_missing_analysis <- function(cfg = cfg) {
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

  message("[1/7] 시군구 shp 로드")
  sig_sf <- load_sigungu_sf(cfg$shp_path)

  message("[2/7] 기상 자료 로드 및 centroid 기반 매핑")
  wx <- load_summer_weather(cfg$weather_path, cfg$station_path, year_keep = cfg$weather_year)
  sig_sf <- attach_climate_to_sigungu(sig_sf, wx$weather_summer, wx$station_info)

  message("[3/7] 토지피복 로드 및 시군구 조인")
  lc <- load_landcover(
    cfg$landcover_path,
    include_grass_in_forest = cfg$include_grass_in_forest,
    include_wetland_in_water = cfg$include_wetland_in_water
  )
  sig_sf <- attach_landcover_to_sigungu(sig_sf, lc)

  message("[4/7] 벌 포인트 로드 및 시군구 집계")
  cerana_df <- load_bee_points(cfg$cerana_path, "cerana")
  mellifera_df <- load_bee_points(cfg$mellifera_path, "mellifera")
  counts <- aggregate_point_counts(cerana_df, mellifera_df, sig_sf)

  message("[5/7] 분석 데이터셋 구성")
  analysis_dat <- make_analysis_dataset(
    sig_sf,
    counts_total = counts$by_total,
    counts_species = counts$by_species,
    response_transform = cfg$response_transform
  )

  total_df <- analysis_dat$total
  species_df <- analysis_dat$by_species

  # CSV 저장용 분석 데이터셋
  readr::write_csv(total_df, file.path(cfg$out_dir, "municipal_analysis_dataset.csv"))
  write_qc_summary(sig_sf, total_df, file.path(cfg$out_dir, "qc_summary.txt"))

  message("[6/7] 기술통계/상관분석 표 생성")
  sbv4_tbl <- make_desc_corr_table(total_df, response = "sbv_count")
  write_outputs(
    sbv4_tbl,
    out_csv = file.path(cfg$out_dir, "table_sbv4_desc_corr.csv"),
    out_tex = file.path(cfg$out_dir, "table_sbv4_desc_corr.tex"),
    caption = "Descriptive statistics and correlations with municipal SBV counts.",
    label = "tb:sbv4",
    digits = 3
  )

  message("[7/7] OLS/CAR 적합 및 공간 진단")
  model_results <- list()
  diag_list <- list()

  try({
    model_results$climate <- fit_ols_car_pair(
      sig_sf_full = sig_sf,
      data_full = total_df,
      response = "response",
      predictors = c("sum_sun", "sum_rain", "mean_temp", "ratio_sun"),
      model_name = "climate_only"
    )
    diag_list[[length(diag_list) + 1]] <- model_results$climate$diagnostics
  }, silent = TRUE)

  try({
    model_results$landcover <- fit_ols_car_pair(
      sig_sf_full = sig_sf,
      data_full = total_df,
      response = "response",
      predictors = c("forest", "agri", "urban", "water"),
      model_name = "landcover_only"
    )
    diag_list[[length(diag_list) + 1]] <- model_results$landcover$diagnostics
  }, silent = TRUE)

  # 연구일지 10/01 기준 최종 변수 후보: sum_sun, sum_rain, water
  try({
    model_results$final <- fit_ols_car_pair(
      sig_sf_full = sig_sf,
      data_full = total_df,
      response = "response",
      predictors = c("sum_sun", "sum_rain", "water"),
      model_name = "final_combined"
    )
    diag_list[[length(diag_list) + 1]] <- model_results$final$diagnostics
  }, silent = TRUE)

  if (length(diag_list) == 0) {
    stop("모든 공간 모형 적합이 실패했습니다. shapefile/landcover/weather 조인을 먼저 확인하세요.", call. = FALSE)
  }

  diag_tbl <- bind_rows(diag_list)
  write_outputs(
    diag_tbl,
    out_csv = file.path(cfg$out_dir, "table_spatial_diagnostics.csv"),
    out_tex = file.path(cfg$out_dir, "table_spatial_diagnostics.tex"),
    caption = "Spatial diagnostics and goodness-of-fit comparison between OLS and CAR models.",
    label = "tb:spatialdiag",
    digits = 4
  )

  # 지도 출력
  sig_count_total <- sig_sf %>%
    left_join(total_df %>% select(SIG_CD, sbv_count), by = "SIG_CD")
  save_count_map_total(
    sig_count_total,
    file.path(cfg$out_dir, "fig_municipal_sbv_count_total.png")
  )

  save_count_map_by_species(
    sig_sf,
    species_df %>% select(SIG_CD, species, sbv_count),
    file.path(cfg$out_dir, "fig_municipal_sbv_count_by_species.png")
  )

  if (!is.null(model_results$final)) {
    spatial_comp <- extract_spatial_component(model_results$final$car, model_results$final$data)

    final_effect_sf <- sig_sf %>%
      left_join(
        model_results$final$data %>%
          select(SIG_CD) %>%
          bind_cols(spatial_comp %>% select(spatial_effect)),
        by = "SIG_CD"
      )

    save_spatial_effect_map(
      final_effect_sf,
      file.path(cfg$out_dir, "fig_spatial_effect_final_car.png")
    )

    saveRDS(model_results$final$car, file.path(cfg$out_dir, "final_car_model.rds"))
    sf::st_write(final_effect_sf, file.path(cfg$out_dir, "municipal_sf_with_outputs.gpkg"), delete_dsn = TRUE, quiet = TRUE)
  }

  message("[DONE] 출력 폴더: ", normalizePath(cfg$out_dir, winslash = "/", mustWork = FALSE))
  invisible(list(
    sig_sf = sig_sf,
    total_df = total_df,
    species_df = species_df,
    diagnostics = diag_tbl,
    sbv4 = sbv4_tbl,
    models = model_results
  ))
}

# ---------------------------
# 13) 자동 실행은 하지 않음
# ---------------------------
message("함수 로드 완료. shapefile 경로를 확인한 뒤 run_missing_analysis(cfg) 를 실행하세요.")

# ---------------------------
# 14) 삭제 후 콘솔에서 실행
# ---------------------------
source("bee_missing_analysis_pipeline.R")
cfg$shp_path <- "data/raw/BND_SIGUNGU_PG.shp"
run_missing_analysis(cfg)