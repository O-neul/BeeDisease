# 분석 자료가 들어 있는 프로젝트 폴더를 자동으로 선택합니다.
# 사용하는 PC의 경로가 다르면 아래 후보에 해당 경로를 추가하세요.
PROJECT_DIR_CANDIDATES <- c(
  "C:/.Soyeon/DS/BeeDisease",
  "C:/Users/haa/Downloads/BeeDisease/BeeDisease"
)
PROJECT_DIR <- PROJECT_DIR_CANDIDATES[dir.exists(PROJECT_DIR_CANDIDATES)][1]
if (is.na(PROJECT_DIR)) {
  stop(
    "프로젝트 폴더를 찾지 못했습니다. PROJECT_DIR_CANDIDATES에 실제 경로를 입력하세요."
  )
}
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(janitor)
  library(sf)
  library(spdep)
  library(spatialreg)
  library(ggplot2)
  library(car)
  library(geosphere)
  library(lubridate)
})

# ---------------------------
# 0) 공통 유틸
# ---------------------------
require_pkg <- function(pkgs) {
  miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(miss)) stop("패키지 설치 필요: ", paste(miss, collapse=", "))
}
require_pkg(c("dplyr","readr","stringr","janitor","sf","spdep","spatialreg",
              "ggplot2","tidyr","geosphere","lubridate"))

guess_delim_from_lines <- function(x, n = 20) {
  lines <- head(x, n)
  cand  <- c("," = ",", "\t" = "\t", ";" = ";", "|" = "|")
  counts <- sapply(cand, function(d) sum(str_count(lines, fixed(d))))
  names(which.max(counts))[1]
}

# CP949/UTF-8/구분자/엑셀까지 탄력 로더
load_table_resilient <- function(path) {
  if (!file.exists(path)) stop(sprintf("파일 없음: %s", path))
  message(sprintf("[INFO] 파일 크기: %s bytes", file.size(path)))
  raw <- read_file_raw(path)
  enc_guess <- try(guess_encoding(raw[1:min(length(raw), 200000)]), silent = TRUE)
  enc <- "CP949"
  if (!inherits(enc_guess, "try-error") && nrow(enc_guess) > 0) enc <- enc_guess$encoding[1]
  message(sprintf("[INFO] 인코딩 추정: %s", enc))
  
  lines <- read_lines(path, locale = locale(encoding = enc), progress = FALSE)
  if (length(lines) > 0 && !str_detect(lines[length(lines)], "\r$|\n$")) {
    lines[length(lines)] <- paste0(lines[length(lines)], "\n")
  }
  
  delim <- guess_delim_from_lines(lines, n = 20)
  message(sprintf("[INFO] 구분자 추정: '%s'", ifelse(delim == "\t", "\\t", delim)))
  
  tf <- tempfile(fileext = ".csv"); write_lines(lines, tf)
  
  df <- try(read_delim(tf, delim = delim, locale = locale(encoding = enc),
                       na = c("", "NA", "-", "NaN"), guess_max = 100000,
                       show_col_types = FALSE),
            silent = TRUE)
  if (!inherits(df, "try-error") && nrow(df) > 0 && ncol(df) > 1) return(df)
  
  if (requireNamespace("data.table", quietly = TRUE)) {
    message("[INFO] readr 실패 → data.table::fread 시도")
    df2 <- try(data.table::fread(tf, encoding = enc, sep = delim,
                                 na.strings = c("", "NA", "-", "NaN")),
               silent = TRUE)
    if (!inherits(df2, "try-error") && nrow(df2) > 0 && ncol(df2) > 1) return(as.data.frame(df2))
  }
  
  if (requireNamespace("readxl", quietly = TRUE)) {
    message("[INFO] CSV 실패 → 엑셀 시도(readxl)")
    fmt <- try(readxl::excel_format(path), silent = TRUE)
    if (!inherits(fmt, "try-error") && !is.na(fmt) && fmt != "unknown") {
      df3 <- try(readxl::read_excel(path), silent = TRUE)
      if (!inherits(df3, "try-error") && nrow(df3) > 0) return(as.data.frame(df3))
    }
  }
  stop("테이블 로드 실패: 인코딩/구분자/파일형식 문제")
}

# 이름 정규화(공백/괄호/전각 공백 제거)
norm_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\u3000", "", x)      # 전각 공백 제거
  x <- gsub("\\s+", "", x)        # 모든 공백 제거
  x <- gsub("\\(.*\\)$", "", x)   # 말미 괄호 제거
  x
}

# ---------------------------
# 1) 데이터 경로
# ---------------------------
PATH_LANDCOVER <- "data/raw/환경부_환경공간정보_년도별 대분류토지피복통계 현황_20230901.csv"
PATH_SIG_SHP   <- "data/raw/BND_SIGUNGU_PG.shp"
PATH_BEE_C     <- "data/final/6_Apis_cerana_with_nearby.csv"
PATH_BEE_M     <- "data/final/6_Apis_mellifera_with_nearby.csv"
PATH_WEATHER   <- "data/raw/weather.csv"
PATH_STATIONS  <- "data/raw/weather_station_info.csv"

# 기존 기후 분석 코드와 동일하게 2018년 여름 자료를 사용합니다.
# 논문 분석 설정이 달랐다면 이 값만 변경하세요.
WEATHER_YEAR <- 2018

# ---------------------------
# 2) 기후자료 준비 및 벌 관측치에 결합
# ---------------------------
prepare_summer_climate <- function(weather_path = PATH_WEATHER,
                                   stations_path = PATH_STATIONS,
                                   year_keep = WEATHER_YEAR) {
  if (!file.exists(weather_path)) stop("기상 파일이 없습니다: ", weather_path)
  if (!file.exists(stations_path)) stop("기상관측소 파일이 없습니다: ", stations_path)

  weather <- read.csv(weather_path, fileEncoding = "CP949", check.names = TRUE)
  stations <- read.csv(stations_path, fileEncoding = "CP949", check.names = TRUE)

  req_weather <- c(
    "지점", "일시", "평균기온..C.", "일최다강수량.mm.",
    "합계일조시간.hr.", "일조율"
  )
  req_stations <- c("지점", "위도", "경도")

  miss_weather <- setdiff(req_weather, names(weather))
  miss_stations <- setdiff(req_stations, names(stations))

  if (length(miss_weather)) {
    stop(
      "weather.csv에 필요한 열이 없습니다: ",
      paste(miss_weather, collapse = ", "),
      "\n현재 열 이름: ", paste(names(weather), collapse = ", ")
    )
  }
  if (length(miss_stations)) {
    stop(
      "weather_station_info.csv에 필요한 열이 없습니다: ",
      paste(miss_stations, collapse = ", ")
    )
  }

  weather$일시 <- lubridate::ymd(weather$일시)

  summer_climate <- weather |>
    dplyr::filter(!is.na(일시)) |>
    dplyr::mutate(
      year = lubridate::year(일시),
      month = lubridate::month(일시)
    ) |>
    dplyr::filter(year == year_keep, month %in% c(6, 7, 8)) |>
    dplyr::group_by(지점) |>
    dplyr::summarise(
      mean_temp = mean(평균기온..C., na.rm = TRUE),
      sum_rain = mean(일최다강수량.mm., na.rm = TRUE),
      sum_sun = mean(합계일조시간.hr., na.rm = TRUE),
      ratio_sun = mean(일조율, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::transmute(
      station_id = as.character(지점),
      mean_temp, sum_sun, ratio_sun, sum_rain
    )

  stations_clean <- stations |>
    dplyr::transmute(
      station_id = as.character(지점),
      lat = suppressWarnings(as.numeric(위도)),
      lon = suppressWarnings(as.numeric(경도))
    ) |>
    tidyr::drop_na(station_id, lat, lon)

  if (nrow(summer_climate) == 0) {
    stop("선택한 연도의 여름 기후자료가 없습니다. WEATHER_YEAR를 확인하세요: ", year_keep)
  }

  message(
    "[OK] 기후자료 준비: year=", year_keep,
    ", 기상 관측소 수=", nrow(summer_climate)
  )

  list(climate = summer_climate, stations = stations_clean)
}

nearest_station_id <- function(lat, lon, stations_df) {
  lat <- suppressWarnings(as.numeric(lat))
  lon <- suppressWarnings(as.numeric(lon))
  if (is.na(lat) || is.na(lon)) return(NA_character_)

  d <- geosphere::distHaversine(
    matrix(c(lon, lat), ncol = 2),
    as.matrix(stations_df[, c("lon", "lat")])
  )
  as.character(stations_df$station_id[which.min(d)])
}

attach_climate_to_points <- function(bee_path, climate_tbl, stations_tbl,
                                     lon_col = "Longitude",
                                     lat_col = "Latitude") {
  if (!file.exists(bee_path)) stop("벌 자료가 없습니다: ", bee_path)
  bee <- read.csv(bee_path, check.names = TRUE)

  if (!all(c(lon_col, lat_col) %in% names(bee))) {
    stop("벌 자료에 좌표 열이 없습니다: ", lon_col, ", ", lat_col)
  }

  climate_vars <- c("mean_temp", "sum_sun", "ratio_sun", "sum_rain")
  if (all(climate_vars %in% names(bee))) {
    message("[INFO] 벌 자료에 기후 변수가 이미 있어 기존 값을 사용합니다: ", bee_path)
    return(bee)
  }

  # 불완전하게 남아 있는 기존 기후 열이 있으면 제거한 후 다시 결합합니다.
  bee <- bee |> dplyr::select(-dplyr::any_of(c(climate_vars, "nearest_station")))

  bee$nearest_station <- vapply(
    seq_len(nrow(bee)),
    function(i) nearest_station_id(
      bee[[lat_col]][i], bee[[lon_col]][i], stations_tbl
    ),
    character(1)
  )

  out <- bee |>
    dplyr::left_join(
      climate_tbl,
      by = c("nearest_station" = "station_id")
    )

  na_by_var <- vapply(
    climate_vars,
    function(v) sum(is.na(out[[v]])),
    integer(1)
  )
  message(
    "[INFO] 기후 결합 후 NA: ",
    paste(paste0(names(na_by_var), "=", na_by_var), collapse = ", ")
  )
  out
}

# ---------------------------
# 3) 환경부 토지피복 로더(한글 컬럼 기반, 최신 연도, 비율 변환)
#    forest=산림(+초지 옵션), water=수역(+습지 옵션)
# ---------------------------
load_landcover_kmoe <- function(path_landcover,
                                use_grass_in_forest = TRUE,
                                use_wetland_in_water = TRUE,
                                to_ratio = TRUE) {
  lc_raw <- load_table_resilient(path_landcover)
  req <- c("자료년도","국가구분","시도","시군구",
           "시가화건조지역","농업지역","산림지역",
           "초지","습지","수역","합계")
  miss <- setdiff(req, names(lc_raw))
  if (length(miss)) stop("landcover 필수 컬럼 누락: ", paste(miss, collapse=", "))
  
  yr_max <- suppressWarnings(max(lc_raw$자료년도, na.rm = TRUE))
  lc <- lc_raw %>%
    filter(국가구분 %in% c("대한민국","남한")) %>%
    filter(자료년도 == yr_max)
  
  forest_base <- lc[["산림지역"]]
  if (use_grass_in_forest && "초지" %in% names(lc)) forest_base <- forest_base + lc[["초지"]]
  
  water_base <- lc[["수역"]]
  if (use_wetland_in_water && "습지" %in% names(lc)) water_base <- water_base + lc[["습지"]]
  
  if (to_ratio) {
    total <- as.numeric(lc[["합계"]])
    forest_v <- as.numeric(forest_base) / total
    agri_v   <- as.numeric(lc[["농업지역"]]) / total
    urban_v  <- as.numeric(lc[["시가화건조지역"]]) / total
    water_v  <- as.numeric(water_base) / total
  } else {
    forest_v <- as.numeric(forest_base)
    agri_v   <- as.numeric(lc[["농업지역"]])
    urban_v  <- as.numeric(lc[["시가화건조지역"]])
    water_v  <- as.numeric(water_base)
  }
  
  lc_out <- lc %>%
    transmute(
      SIGUNGU_NM = norm_key(시군구),
      forest = forest_v, agri = agri_v, urban = urban_v, water = water_v
    ) %>%
    group_by(SIGUNGU_NM) %>%
    summarise(
      forest = mean(forest, na.rm=TRUE),
      agri   = mean(agri,   na.rm=TRUE),
      urban  = mean(urban,  na.rm=TRUE),
      water  = mean(water,  na.rm=TRUE),
      .groups = "drop"
    )
  message("[OK] Landcover 준비 완료: 연도=", yr_max, ", 스케일=", ifelse(to_ratio,"비율(0-1)","원단위"))
  lc_out
}

# ---------------------------
# 3) 시군구 경계 로드
# ---------------------------
load_sigungu_sf <- function(path_shp = PATH_SIG_SHP) {
  sf <- st_read(path_shp, options = "ENCODING=CP949", quiet = TRUE)
  if (!all(st_is_valid(sf))) sf <- st_make_valid(sf)
  if (is.na(st_crs(sf))) stop("시군구 shp에 CRS 없음. st_set_crs로 EPSG 지정 필요.")
  sf
}

# ---------------------------
# 4) shp 키 자동 탐지 + 포인트→landcover 안전 조인
# ---------------------------
detect_sigungu_keys <- function(sig_sf) {
  nms <- names(sig_sf)
  name_pat <- "(SIG(_KOR)?_?NM|SIGUNGU_?NM|SGG_?NM|EMD(C)?_?NM|SIG_NM|SIGNM|시군구.?명|시군구.?이름)"
  code_pat <- "(SIG(_?CD)?|SIGUNGU_?CD|SGG_?CD|ADM_?CD|LAWD_?CD|CODE|코드)"
  name_col <- nms[grepl(name_pat, nms, ignore.case = TRUE)]
  code_col <- nms[grepl(code_pat, nms, ignore.case = TRUE)]
  if (length(name_col) == 0 && "SIG_KOR_NM" %in% nms) name_col <- "SIG_KOR_NM"
  if (length(name_col) == 0 && "SIG_NM"     %in% nms) name_col <- "SIG_NM"
  if (length(code_col) == 0 && "SIG_CD"     %in% nms) code_col <- "SIG_CD"
  list(
    name_col = if (length(name_col)) name_col[1] else NA_character_,
    code_col = if (length(code_col)) code_col[1] else NA_character_
  )
}

attach_landcover_to_points <- function(bee_input, sigungu_sf, landcover_tbl,
                                       lon_col="Longitude", lat_col="Latitude") {
  bee <- if (is.data.frame(bee_input)) {
    as.data.frame(bee_input)
  } else if (is.character(bee_input) && length(bee_input) == 1) {
    read.csv(bee_input, check.names = TRUE)
  } else {
    stop("bee_input은 data.frame 또는 파일 경로여야 합니다.")
  }
  stopifnot(all(c(lon_col, lat_col) %in% names(bee)))
  
  keys <- detect_sigungu_keys(sigungu_sf)
  nm_col <- keys$name_col
  cd_col <- keys$code_col
  if (is.na(nm_col) && is.na(cd_col)) {
    stop("shp에서 시군구 이름/코드 컬럼을 찾지 못했습니다. names(sigungu_sf) 확인 요망.")
  }
  
  pts <- sf::st_as_sf(bee, coords = c(lon_col, lat_col), crs = 4326, remove = FALSE)
  pts <- sf::st_transform(pts, sf::st_crs(sigungu_sf))
  sgg2 <- suppressWarnings(sf::st_buffer(sigungu_sf, 0))
  
  pick_cols <- c(na.omit(c(nm_col, cd_col)))
  joined <- sf::st_join(pts, sgg2[, pick_cols, drop = FALSE],
                        join = sf::st_intersects, left = TRUE)
  df <- sf::st_drop_geometry(joined)
  
  if (!is.na(nm_col) && nm_col %in% names(df)) df[[nm_col]] <- norm_key(df[[nm_col]])
  landcover_tbl$SIGUNGU_NM <- norm_key(landcover_tbl$SIGUNGU_NM)
  
  # 1순위: 이름으로 조인
  if (!is.na(nm_col) && nm_col %in% names(df)) {
    df1 <- dplyr::left_join(df, landcover_tbl, by = setNames("SIGUNGU_NM", nm_col))
    na_lc1 <- rowSums(is.na(df1[, c("forest","agri","urban","water")]))
    if (sum(na_lc1 > 0) < nrow(df1)) {
      message(sprintf("[INFO] 이름기반 조인: landcover NA %d / %d", sum(na_lc1 > 0), nrow(df1)))
      return(df1)
    } else {
      message("[WARN] 이름기반 조인 결과가 전부 NA — 코드 매핑으로 재시도")
    }
  }
  
  # 2순위: 코드→이름 매핑 후 조인
  if (!is.na(cd_col) && !is.na(nm_col) &&
      cd_col %in% names(sigungu_sf) && nm_col %in% names(sigungu_sf)) {
    xwalk <- sigungu_sf |>
      sf::st_drop_geometry() |>
      dplyr::select(!!nm_col, !!cd_col) |>
      dplyr::mutate(
        !!nm_col := norm_key(.data[[nm_col]]),
        !!cd_col := as.character(.data[[cd_col]])
      ) |>
      dplyr::distinct()
    
    if (cd_col %in% names(df)) {
      df$.__sgg_code__ <- as.character(df[[cd_col]])
      df2 <- dplyr::left_join(df, xwalk, by = setNames(cd_col, ".__sgg_code__"))
      df2$.__sgg_name__ <- if (nm_col %in% names(df2)) df2[[nm_col]] else NA_character_
      if (nm_col %in% names(df2)) {
        df2$.__sgg_name__ <- ifelse(is.na(df2$.__sgg_name__), df2[[nm_col]], df2$.__sgg_name__)
      }
      df2 <- dplyr::left_join(
        df2,
        dplyr::rename(landcover_tbl, .__sgg_name__ = SIGUNGU_NM),
        by = ".__sgg_name__"
      )
      na_lc2 <- rowSums(is.na(df2[, c("forest","agri","urban","water")]))
      message(sprintf("[INFO] 코드↔이름 매핑 후 조인: landcover NA %d / %d",
                      sum(na_lc2 > 0), nrow(df2)))
      return(df2)
    }
  }
  
  message("[FAIL] landcover 조인 실패. shp 컬럼/값을 확인하세요.")
  if (!is.na(nm_col) && nm_col %in% names(df)) utils::str(head(df[[nm_col]]))
  if (!is.na(cd_col) && cd_col %in% names(df)) utils::str(head(df[[cd_col]]))
  df
}

# ---------------------------
# 5) 이웃/가중치 및 통합 CAR 적합
# ---------------------------
build_neighbors <- function(df, lon_col="Longitude", lat_col="Latitude",
                            crs_proj=5181, d2_m=6000, style="B") {
  sf_pts  <- st_as_sf(df, coords = c(lon_col, lat_col), crs = 4326)
  sf_proj <- suppressWarnings(st_transform(sf_pts, crs_proj))
  coords  <- st_coordinates(sf_proj)
  nb      <- dnearneigh(coords, d1 = 0, d2 = d2_m)
  iso     <- which(card(nb) == 0)
  if (length(iso)) {
    message(sprintf("[WARN] 고립점 %d개(이웃 0). zero.policy=TRUE로 처리.", length(iso)))
  }
  listw <- nb2listw(nb, style = style, zero.policy = TRUE)
  list(nb=nb, listw=listw, sf_proj=sf_proj)
}

# 논문의 기존 두 모형과 직접 비교할 수 있도록,
# Table 4의 기후 변수와 Table 5의 토지피복 변수를 한 모형에 동시에 적합합니다.
# urban은 토지피복 조성자료의 기준범주로 제외합니다.
fit_car_both <- function(
    df,
    d2_m = 6000,
    keep_lc = c("forest", "agri", "water"),
    keep_clim = c("mean_temp", "sum_sun", "ratio_sun", "sum_rain")
) {
  required <- c("발생두수", "species", "Longitude", "Latitude")
  missing_required <- setdiff(required, names(df))
  if (length(missing_required)) {
    stop("필수 변수가 없습니다: ", paste(missing_required, collapse=", "))
  }

  lc_all <- c("forest", "agri", "urban", "water")
  clim_all <- c("mean_temp", "sum_sun", "ratio_sun", "sum_rain")

  cand_lc <- intersect(keep_lc, intersect(lc_all, names(df)))
  cand_clim <- intersect(keep_clim, intersect(clim_all, names(df)))

  missing_lc <- setdiff(keep_lc, names(df))
  missing_clim <- setdiff(keep_clim, names(df))
  if (length(missing_lc)) {
    stop("통합 데이터에 토지피복 변수가 없습니다: ", paste(missing_lc, collapse=", "))
  }
  if (length(missing_clim)) {
    stop("통합 데이터에 기후 변수가 없습니다: ", paste(missing_clim, collapse=", "))
  }

  used_vars <- c(cand_lc, cand_clim)
  message("[MODEL] 사용 landcover: ", paste(cand_lc, collapse=", "))
  message("[MODEL] 사용 climate: ", paste(cand_clim, collapse=", "))

  fml_obj <- stats::as.formula(
    paste0("발생두수 ~ (", paste(used_vars, collapse=" + "), ") * species")
  )
  message("[MODEL] 공식: ", deparse(fml_obj))

  need <- c(required, used_vars)
  dat <- df |>
    dplyr::select(dplyr::all_of(need)) |>
    dplyr::mutate(
      발생두수 = suppressWarnings(readr::parse_number(as.character(발생두수))),
      Longitude = suppressWarnings(readr::parse_number(as.character(Longitude))),
      Latitude = suppressWarnings(readr::parse_number(as.character(Latitude))),
      dplyr::across(
        dplyr::all_of(used_vars),
        ~ suppressWarnings(readr::parse_number(as.character(.x)))
      ),
      species = factor(as.character(species), levels=c("cerana", "mellifera"))
    ) |>
    tidyr::drop_na(dplyr::all_of(need))

  if (nrow(dat) < 20) {
    stop(sprintf(
      "완전사례 선택 후 유효 행이 너무 적습니다 (n=%d). species 값과 변수 결측치를 확인하세요.",
      nrow(dat)
    ))
  }

  species_n <- table(dat$species)
  if (any(species_n == 0)) {
    stop("두 종 중 하나의 유효 관측치가 0개입니다. species 값을 확인하세요.")
  }
  message("[MODEL] 분석 표본: cerana=", species_n[["cerana"]],
          ", mellifera=", species_n[["mellifera"]],
          ", total=", nrow(dat))

  zero_var <- used_vars[vapply(dat[used_vars], function(x) {
    is.na(stats::sd(x)) || stats::sd(x) == 0
  }, logical(1))]
  if (length(zero_var)) {
    stop("변동이 없는 설명변수가 있습니다: ", paste(zero_var, collapse=", "))
  }

  g <- build_neighbors(dat, d2_m=d2_m, style="B")

  fit <- spatialreg::spautolm(
    formula=fml_obj,
    data=dat,
    listw=g$listw,
    family="CAR",
    zero.policy=TRUE
  )
  fit$call$formula <- fml_obj

  b <- stats::coef(fit)
  me_fun <- function(base) {
    base_eff <- if (base %in% names(b)) unname(b[base]) else NA_real_
    int_names <- c(
      paste0(base, ":speciesmellifera"),
      paste0("speciesmellifera:", base)
    )
    int_match <- intersect(int_names, names(b))
    int_eff <- if (length(int_match)) unname(b[int_match[1]]) else 0
    list(
      cerana=base_eff,
      mellifera=base_eff + int_eff,
      species_difference=int_eff
    )
  }
  me_list <- lapply(used_vars, me_fun)
  names(me_list) <- used_vars

  list(
    fit=fit,
    data=dat,
    listw=g$listw,
    variables=used_vars,
    marginal_effects=me_list
  )
}

# ---------------------------
# 6) 진단/시각화 헬퍼
# ---------------------------
diagnose_both <- function(df) {
  nm <- names(df)
  lon <- nm[grepl("^(Longitude|lon|long|x|경도)$", nm, ignore.case=TRUE)][1]
  lat <- nm[grepl("^(Latitude|lat|y|위도)$",     nm, ignore.case=TRUE)][1]
  cat("# rows:", nrow(df), "\n")
  cat("# NAs by var:\n")
  na_cnt <- sapply(
    c("발생두수", "species", "forest", "agri", "urban", "water",
      "mean_temp", "sum_sun", "ratio_sun", "sum_rain", lon, lat),
    function(x) if (x %in% nm) sum(is.na(df[[x]])) else NA_integer_
  )
  print(na_cnt)
  x <- suppressWarnings(as.numeric(df[[lon]])); y <- suppressWarnings(as.numeric(df[[lat]]))
  cat("# coord summary:\n"); print(summary(x)); print(summary(y))
  out_kr <- which(x<124 | x>132 | y<33 | y>39)
  cat("# out of Korea bbox:", length(out_kr), "\n")
}

make_X_from_model <- function(fit, newdata) {
  b  <- coef(fit)
  ff <- formula(fit)
  X  <- model.matrix(ff, data = newdata)
  if ("(Intercept)" %in% names(b) && !("(Intercept)" %in% colnames(X))) {
    X <- cbind("(Intercept)"=1, X)
  }
  miss <- setdiff(names(b), colnames(X))
  if (length(miss)) {
    Z <- matrix(0, nrow=nrow(X), ncol=length(miss)); colnames(Z) <- miss
    X <- cbind(X, Z)
  }
  X <- X[, names(b), drop=FALSE]
  list(X=X, b=as.numeric(b))
}

# ================================================================
# 7) 실행 파이프라인
# ================================================================
# 7-1) 기후, Landcover, Shapefile 로드
weather_obj <- prepare_summer_climate(
  PATH_WEATHER, PATH_STATIONS, year_keep = WEATHER_YEAR
)
landcover_tbl <- load_landcover_kmoe(
  PATH_LANDCOVER,
  use_grass_in_forest = TRUE,
  use_wetland_in_water = TRUE,
  to_ratio = TRUE
)
sig_sf <- load_sigungu_sf(PATH_SIG_SHP)

# 7-2) 벌 포인트에 기후자료를 먼저 붙이고, 같은 행에 토지피복을 추가
cerana_climate <- attach_climate_to_points(
  PATH_BEE_C, weather_obj$climate, weather_obj$stations
)
mellifera_climate <- attach_climate_to_points(
  PATH_BEE_M, weather_obj$climate, weather_obj$stations
)

cerana_df <- attach_landcover_to_points(
  cerana_climate, sig_sf, landcover_tbl
) |>
  dplyr::mutate(species = "cerana")

mellifera_df <- attach_landcover_to_points(
  mellifera_climate, sig_sf, landcover_tbl
) |>
  dplyr::mutate(species = "mellifera")

# 7-3) 토지피복과 기후변수를 모두 가진 통합 데이터 생성
both <- dplyr::bind_rows(cerana_df, mellifera_df) |>
  dplyr::mutate(species = factor(species, levels = c("cerana", "mellifera")))

message("[CHECK] 통합 데이터 열: ", paste(names(both), collapse = ", "))
diagnose_both(both)


# 7-4) 논문 Table 4 + Table 5 변수의 통합 CAR 모형
# 객체 이름을 res_combined로 통일했습니다.
res_combined <- fit_car_both(
  both,
  d2_m=6000,
  keep_lc=c("forest", "agri", "water"),
  keep_clim=c("mean_temp", "sum_sun", "ratio_sun", "sum_rain")
)

# 핵심 결과
print(summary(res_combined$fit))
message("[Marginal Effects: landcover + climate combined CAR]")
print(res_combined$marginal_effects)

# 비교에 필요한 진단값
message("[MODEL DIAGNOSTICS]")
print(c(
  n=nrow(res_combined$data),
  lambda=res_combined$fit$lambda,
  AIC=stats::AIC(res_combined$fit),
  logLik=as.numeric(stats::logLik(res_combined$fit))
))

message("[RESIDUAL MORAN TEST]")
print(spdep::moran.test(
  residuals(res_combined$fit),
  res_combined$listw,
  zero.policy=TRUE
))

# 콘솔 결과를 텍스트 파일로 저장: 논문 결과와 비교할 때 이 파일만 보내도 됩니다.
capture.output(
  {
    cat("=== Combined CAR model summary ===\n")
    print(summary(res_combined$fit))
    cat("\n=== Variables ===\n")
    print(res_combined$variables)
    cat("\n=== Marginal effects ===\n")
    print(res_combined$marginal_effects)
    cat("\n=== Model diagnostics ===\n")
    print(c(
      n=nrow(res_combined$data),
      lambda=res_combined$fit$lambda,
      AIC=stats::AIC(res_combined$fit),
      logLik=as.numeric(stats::logLik(res_combined$fit))
    ))
    cat("\n=== Residual Moran test ===\n")
    print(spdep::moran.test(
      residuals(res_combined$fit),
      res_combined$listw,
      zero.policy=TRUE
    ))
  },
  file="combined_CAR_results_climate_landcover.txt"
)
message("[SAVE] combined_CAR_results_climate_landcover.txt 저장 완료: ", getwd())

# 예측값(Xβ) 시각화
tmp <- make_X_from_model(res_combined$fit, res_combined$data)
res_combined$data$yhat <- drop(tmp$X %*% tmp$b)

print(
  ggplot(res_combined$data, aes(x=species, y=yhat)) +
    geom_boxplot(outlier.shape=NA, width=0.5) +
    labs(
      title="CAR(Xβ) 예측 분포 — 토지피복·기후 통합 모형",
      x="종", y="예측 발생두수"
    ) +
    theme_minimal()
)
