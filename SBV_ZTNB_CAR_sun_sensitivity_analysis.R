# ======================================================================
# SBV 발생두수 민감도 분석: zero-truncated negative-binomial CAR
#
# 목적
# 1) 기존 논문과 같은 기후/토지피복 분리 모형을 count 자료에 맞게 재적합
# 2) sum_sun과 ratio_sun의 동시 포함에 따른 계수 불안정성 확인
# 3) sum_sun 단독 모형과 ratio_sun 단독 모형으로 결과의 강건성 평가
# 4) MCMC 진단, IRR, 종별 결합효과, 공간효과 및 조건부 예측표 자동 저장
#
# 반응변수: 발생두수 원자료(count), 발생 기록이 있는 지점만 사용(Y > 0)
# 우도: zero-truncated negative binomial
# 공간효과: 6 km 거리 기반 exact sparse CAR(escar)
# ======================================================================

PROJECT_DIR <- "C:/.Soyeon/DS/BeeDisease"
if (!dir.exists(PROJECT_DIR)) {
  stop("PROJECT_DIR가 존재하지 않습니다. 파일 상단의 PROJECT_DIR를 수정하세요: ", PROJECT_DIR)
}
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(sf)
  library(spdep)
  library(geosphere)
  library(lubridate)
  library(ggplot2)
  library(tibble)
  library(brms)
  library(posterior)
  library(loo)
  library(Matrix)
})

# ---------------------------
# 0) 실제 프로젝트 경로
# ---------------------------
PATH_FLORAL    <- "data/processed/3_Floral_Source_ID.csv"
PATH_WEATHER   <- "data/raw/weather.csv"
PATH_STATIONS  <- "data/raw/weather_station_info.csv"
PATH_LANDCOVER <- "data/raw/환경부_환경공간정보_년도별 대분류토지피복통계 현황_20230901.csv"
PATH_SIG_SHP   <- "data/raw/BND_SIGUNGU_PG.shp"
PATH_BEE_C     <- "data/final/6_Apis_cerana_with_nearby.csv"
PATH_BEE_M     <- "data/final/6_Apis_mellifera_with_nearby.csv"

# ---------------------------
# 1) 기후·밀원 데이터 함수
# ---------------------------
load_floral <- function(path = PATH_FLORAL) {
  floral <- read.csv(path, stringsAsFactors = FALSE) %>%
    rename_with(tolower) %>%
    rename(
      species   = 종,
      latitude  = 위도,
      longitude = 경도,
      abundance = 카운트
    )

  stopifnot(all(c("id", "species", "latitude", "longitude", "abundance") %in% names(floral)))
  floral
}

prepare_summer_temp <- function(weather_path = PATH_WEATHER,
                                stations_path = PATH_STATIONS,
                                year_keep = 2018) {
  weather  <- read.csv(weather_path, fileEncoding = "CP949")
  stations <- read.csv(stations_path, fileEncoding = "CP949")

  weather$일시 <- ymd(weather$일시)

  weather <- weather %>%
    filter(!is.na(일시)) %>%
    mutate(
      year = year(일시),
      month = month(일시),
      season = case_when(
        month %in% c(12, 1, 2) ~ "winter",
        month %in% c(3, 4, 5)  ~ "spring",
        month %in% c(6, 7, 8)  ~ "summer",
        TRUE                   ~ "fall"
      )
    ) %>%
    filter(year == year_keep)

  summer_temp <- weather %>%
    group_by(지점, season) %>%
    summarise(
      mean_temp = mean(평균기온..C., na.rm = TRUE),
      sum_rain  = mean(일최다강수량.mm., na.rm = TRUE),
      sum_sun   = mean(합계일조시간.hr., na.rm = TRUE),
      ratio_sun = mean(일조율, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(season == "summer") %>%
    rename(station_id = 지점)

  stations_clean <- stations %>%
    select(지점, 위도, 경도) %>%
    rename(station_id = 지점, lat = 위도, lon = 경도)

  list(summer_temp = summer_temp, stations = stations_clean)
}

nearest_station_id <- function(lat, lon, stations_df) {
  d <- distHaversine(
    matrix(c(lon, lat), ncol = 2),
    matrix(c(stations_df$lon, stations_df$lat), ncol = 2)
  )
  stations_df$station_id[which.min(d)]
}

nearby_floral_ids <- function(lat, lon, floral_df, threshold_km = 6) {
  d <- distHaversine(
    matrix(c(lon, lat), ncol = 2),
    as.matrix(floral_df[, c("longitude", "latitude")])
  ) / 1000

  floral_df$id[d <= threshold_km]
}

flower_stats <- function(ids, floral_df, prefer_list) {
  ids <- ids[!is.na(ids) & ids != ""]
  matched <- floral_df %>% filter(id %in% ids)

  preferred    <- matched %>% filter(species %in% prefer_list)
  nonpreferred <- matched %>% filter(!species %in% prefer_list)

  tibble(
    prefer_count    = nrow(preferred),
    prefer_abund    = sum(preferred$abundance, na.rm = TRUE),
    nonprefer_count = nrow(nonpreferred),
    nonprefer_abund = sum(nonpreferred$abundance, na.rm = TRUE)
  )
}

cerana_preferred <- c(
  '밤나무','산밤나무','약밤나무','배추','유채','머귀나무','왕초피','초피나무',
  '산초나무','벼','피마자','다래','쥐다래','개다래','섬다래','감나무','고욤나무',
  '개옺나무','검양옺나무','산검양옺나무','옺나무','붉나무','광대싸리'
)

mellifera_preferred <- c(
  '족제비싸리','아까지나무','갈참나무','굴참나무','떡갈나무','종가시나무',
  '물참나무','상수리나무','신갈나무','졸참나무','가시나무','붉가시나무',
  '갈졸참나무','떡갈참나무','떡신갈나무','주름잎','누운주름잎','덩굴장미',
  '목향장미','생열귀나무','용가시나무','인가목','해당화','고추','애기똥풀',
  '다래','쥐다래','개다래','섬다래','감나무','고욤나무','개옺나무','검양옺나무',
  '산검양옺나무','옺나무','붉나무','광대싸리'
)

prepare_bee_dataset <- function(bee_path, floral_df, summer_temp, stations_df,
                                prefer_list, threshold_km = 6) {
  bee <- read.csv(bee_path)
  stopifnot(all(c("Latitude", "Longitude", "발생두수") %in% names(bee)))

  if (!"nearby_floral_ids" %in% names(bee)) {
    bee$nearby_floral_ids <- lapply(seq_len(nrow(bee)), function(i) {
      nearby_floral_ids(
        lat = bee$Latitude[i],
        lon = bee$Longitude[i],
        floral_df = floral_df,
        threshold_km = threshold_km
      )
    })
  }

  bee$nearest_station <- apply(
    bee[, c("Latitude", "Longitude")],
    1,
    function(row) {
      nearest_station_id(
        as.numeric(row["Latitude"]),
        as.numeric(row["Longitude"]),
        stations_df
      )
    }
  )

  bee2 <- bee %>%
    left_join(summer_temp, by = c("nearest_station" = "station_id"))

  if (!"season" %in% names(bee2)) bee2$season <- "summer"

  nearby_cols <- grep("^nearby_", names(bee2), value = TRUE)

  flower_vars <- bind_rows(lapply(seq_len(nrow(bee2)), function(i) {
    if (length(nearby_cols) == 1 && nearby_cols[1] == "nearby_floral_ids") {
      ids <- bee2$nearby_floral_ids[[i]]
    } else {
      vals <- as.vector(unlist(bee2[i, nearby_cols], use.names = FALSE))
      ids <- vals[!is.na(vals) & vals != ""]
    }
    flower_stats(ids, floral_df, prefer_list)
  }))

  bind_cols(bee2, flower_vars)
}

# ---------------------------
# 2) 토지피복 데이터 함수
# ---------------------------
guess_delim_from_lines <- function(x, n = 20) {
  lines <- head(x, n)
  cand <- c("," = ",", "\t" = "\t", ";" = ";", "|" = "|")
  counts <- sapply(cand, function(d) sum(str_count(lines, fixed(d))))
  names(which.max(counts))[1]
}

load_table_resilient <- function(path) {
  if (!file.exists(path)) stop("파일 없음: ", path)

  raw <- read_file_raw(path)
  enc_guess <- try(guess_encoding(raw[1:min(length(raw), 200000)]), silent = TRUE)
  enc <- "CP949"
  if (!inherits(enc_guess, "try-error") && nrow(enc_guess) > 0) {
    enc <- enc_guess$encoding[1]
  }

  lines <- read_lines(path, locale = locale(encoding = enc), progress = FALSE)
  delim <- guess_delim_from_lines(lines)

  out <- read_delim(
    path,
    delim = delim,
    locale = locale(encoding = enc),
    na = c("", "NA", "-", "NaN"),
    guess_max = 100000,
    show_col_types = FALSE
  )

  as.data.frame(out)
}

norm_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\u3000", "", x)
  x <- gsub("\\s+", "", x)
  x <- gsub("\\(.*\\)$", "", x)
  x
}

load_landcover_kmoe <- function(path_landcover = PATH_LANDCOVER,
                                use_grass_in_forest = TRUE,
                                use_wetland_in_water = TRUE,
                                to_ratio = TRUE) {
  lc_raw <- load_table_resilient(path_landcover)

  req <- c(
    "자료년도", "국가구분", "시도", "시군구", "시가화건조지역", "농업지역",
    "산림지역", "초지", "습지", "수역", "합계"
  )

  miss <- setdiff(req, names(lc_raw))
  if (length(miss)) stop("landcover 필수 컬럼 누락: ", paste(miss, collapse = ", "))

  lc_raw$자료년도 <- suppressWarnings(as.numeric(lc_raw$자료년도))
  yr_max <- max(lc_raw$자료년도, na.rm = TRUE)

  lc <- lc_raw %>%
    filter(grepl("대한민국|남한", 국가구분)) %>%
    filter(자료년도 == yr_max)

  forest_base <- as.numeric(lc[["산림지역"]])
  if (use_grass_in_forest) forest_base <- forest_base + as.numeric(lc[["초지"]])

  water_base <- as.numeric(lc[["수역"]])
  if (use_wetland_in_water) water_base <- water_base + as.numeric(lc[["습지"]])

  total <- as.numeric(lc[["합계"]])

  if (to_ratio) {
    forest_v <- forest_base / total
    agri_v   <- as.numeric(lc[["농업지역"]]) / total
    urban_v  <- as.numeric(lc[["시가화건조지역"]]) / total
    water_v  <- water_base / total
  } else {
    forest_v <- forest_base
    agri_v   <- as.numeric(lc[["농업지역"]])
    urban_v  <- as.numeric(lc[["시가화건조지역"]])
    water_v  <- water_base
  }

  lc %>%
    transmute(
      SIGUNGU_NM = norm_key(시군구),
      forest = forest_v,
      agri = agri_v,
      urban = urban_v,
      water = water_v
    ) %>%
    group_by(SIGUNGU_NM) %>%
    summarise(
      across(c(forest, agri, urban, water), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}

load_sigungu_sf <- function(path_shp = PATH_SIG_SHP) {
  sig <- st_read(path_shp, options = "ENCODING=CP949", quiet = TRUE)
  if (!all(st_is_valid(sig))) sig <- st_make_valid(sig)
  if (is.na(st_crs(sig))) stop("시군구 shp에 CRS가 없습니다.")
  sig
}

detect_sigungu_keys <- function(sig_sf) {
  nms <- names(sig_sf)
  name_pat <- "(SIG(_KOR)?_?NM|SIGUNGU_?NM|SGG_?NM|SIG_NM|SIGNM|시군구.?명|시군구.?이름)"
  code_pat <- "(SIG(_?CD)?|SIGUNGU_?CD|SGG_?CD|ADM_?CD|LAWD_?CD|CODE|코드)"

  name_col <- nms[grepl(name_pat, nms, ignore.case = TRUE)]
  code_col <- nms[grepl(code_pat, nms, ignore.case = TRUE)]

  list(
    name_col = if (length(name_col)) name_col[1] else NA_character_,
    code_col = if (length(code_col)) code_col[1] else NA_character_
  )
}

# 기존 attach_landcover_to_points()는 파일을 다시 읽어서 기후 파생변수를 잃는다.
# 이 함수는 이미 기후변수가 붙은 df_cerana/df_mellifera에 토지피복을 결합한다.
attach_landcover_to_df <- function(bee_df, sigungu_sf, landcover_tbl,
                                   lon_col = "Longitude", lat_col = "Latitude") {
  stopifnot(all(c(lon_col, lat_col) %in% names(bee_df)))

  keys <- detect_sigungu_keys(sigungu_sf)
  nm_col <- keys$name_col

  if (is.na(nm_col)) {
    stop("시군구 shp에서 시군구 명칭 컬럼을 찾지 못했습니다: ",
         paste(names(sigungu_sf), collapse = ", "))
  }

  bee_df <- bee_df %>% mutate(.row_id_join = row_number())

  pts <- st_as_sf(
    bee_df,
    coords = c(lon_col, lat_col),
    crs = 4326,
    remove = FALSE
  ) %>%
    st_transform(st_crs(sigungu_sf))

  sig_use <- sigungu_sf[, nm_col, drop = FALSE]

  joined <- st_join(
    pts,
    suppressWarnings(st_buffer(sig_use, 0)),
    join = st_intersects,
    left = TRUE
  ) %>%
    st_drop_geometry()

  # 경계 위 점이 두 행정구역과 겹쳐 행이 중복되는 경우 첫 매칭만 유지
  duplicated_n <- sum(duplicated(joined$.row_id_join))
  if (duplicated_n > 0) {
    warning("행정경계 중복 매칭 행: ", duplicated_n, "개. 각 점의 첫 매칭을 유지합니다.")
    joined <- joined %>%
      group_by(.row_id_join) %>%
      slice(1) %>%
      ungroup()
  }

  joined[[nm_col]] <- norm_key(joined[[nm_col]])
  landcover_tbl <- landcover_tbl %>%
    mutate(SIGUNGU_NM = norm_key(SIGUNGU_NM))

  out <- left_join(
    joined,
    landcover_tbl,
    by = setNames("SIGUNGU_NM", nm_col)
  ) %>%
    arrange(.row_id_join) %>%
    select(-.row_id_join)

  message(
    "[Landcover join] NA 행: ",
    sum(!complete.cases(out[, c("forest", "agri", "urban", "water")])),
    " / ", nrow(out)
  )

  out
}

# ---------------------------
# 3) 기후 + 토지피복 통합 데이터 생성
# ---------------------------
floral <- load_floral()
wx <- prepare_summer_temp(year_keep = 2018)

# 실제 기존 코드에서 쓰는 객체와 동일한 단계
# df_cerana / df_mellifera에는 발생두수, 좌표, 기후, prefer_count 등이 들어간다.
df_cerana <- prepare_bee_dataset(
  PATH_BEE_C,
  floral_df = floral,
  summer_temp = wx$summer_temp,
  stations_df = wx$stations,
  prefer_list = cerana_preferred,
  threshold_km = 6
) %>%
  mutate(species = "cerana")

df_mellifera <- prepare_bee_dataset(
  PATH_BEE_M,
  floral_df = floral,
  summer_temp = wx$summer_temp,
  stations_df = wx$stations,
  prefer_list = mellifera_preferred,
  threshold_km = 6
) %>%
  mutate(species = "mellifera")

landcover_tbl <- load_landcover_kmoe(
  PATH_LANDCOVER,
  use_grass_in_forest = TRUE,
  use_wetland_in_water = TRUE,
  to_ratio = TRUE
)

sig_sf <- load_sigungu_sf(PATH_SIG_SHP)

cerana_full <- attach_landcover_to_df(df_cerana, sig_sf, landcover_tbl)
mellifera_full <- attach_landcover_to_df(df_mellifera, sig_sf, landcover_tbl)

both_full <- bind_rows(cerana_full, mellifera_full) %>%
  filter(season == "summer") %>%
  mutate(species = factor(species, levels = c("cerana", "mellifera")))

# 실제 모델에 사용할 변수 확인
MODEL_VARS_RAW <- c(
  "mean_temp", "sum_sun", "ratio_sun", "sum_rain",
  "forest", "agri", "water"
)

REQUIRED_VARS <- c(
  "발생두수", "species", "Longitude", "Latitude", MODEL_VARS_RAW
)

missing_vars <- setdiff(REQUIRED_VARS, names(both_full))
if (length(missing_vars)) {
  stop("통합 데이터에 없는 변수: ", paste(missing_vars, collapse = ", "))
}

cat("\n[통합 데이터 변수 확인]\n")
print(intersect(c(
  "발생두수", "species", "Longitude", "Latitude",
  "prefer_count", "mean_temp", "sum_sun", "ratio_sun", "sum_rain",
  "forest", "agri", "urban", "water"
), names(both_full)))


# ---------------------------
# 4) 발생 기록이 있는 지점만을 위한 zero-truncated CAR 모형
# ---------------------------
# INLA 설치본마다 zero-truncated likelihood 이름이 달라지는 문제를 피하기 위해
# brms의 공식 trunc(lb = 0) 문법을 사용한다.
#
# 분석모형:
#   발생두수 | 발생두수 > 0 ~ Poisson 또는 Negative Binomial
#   log(mu_i) = 환경변수 + species + 상호작용 + CAR 공간효과
#
# 주의:
#   이 분석은 SBV가 발생할 확률을 설명하지 않는다.
#   발생 기록이 있는 지점 사이에서 발생 규모가 얼마나 다른지를 설명한다.

# ---------------------------
# 4-1) MCMC 실행 설정
# ---------------------------
N_CHAINS  <- 4
N_ITER    <- 2000
N_WARMUP  <- 1000
N_CORES   <- min(N_CHAINS, max(1, parallel::detectCores(logical = FALSE)))
RANDOM_SEED <- 20260721

# cmdstanr와 CmdStan이 설치되어 있으면 빠른 cmdstanr를 사용하고,
# 아니면 CRAN brms 기본 구성에서 사용하기 쉬운 rstan으로 실행한다.
choose_brms_backend <- function() {
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    cmdstan_ok <- try(cmdstanr::cmdstan_path(), silent = TRUE)
    if (!inherits(cmdstan_ok, "try-error") && nzchar(cmdstan_ok)) {
      return("cmdstanr")
    }
  }
  "rstan"
}

BRMS_BACKEND <- choose_brms_backend()
message("[brms backend] ", BRMS_BACKEND)

options(mc.cores = N_CORES)
if (BRMS_BACKEND == "rstan") {
  rstan::rstan_options(auto_write = TRUE)
}

# ---------------------------
# 4-2) 분석자료와 공간 인접행렬 준비
# ---------------------------
prepare_truncated_car_data <- function(df, d2_m = 6000) {
  x_raw <- c(
    "mean_temp", "sum_sun", "ratio_sun", "sum_rain",
    "forest", "agri", "water"
  )

  need <- c("발생두수", "species", "Longitude", "Latitude", x_raw)

  dat <- df %>%
    select(all_of(need)) %>%
    drop_na() %>%
    mutate(
      species = factor(species, levels = c("cerana", "mellifera")),
      # Stan 코드에는 ASCII 변수명만 전달한다.
      # 한글 반응변수명을 prior/formula 내부에서 직접 사용하면
      # 일부 Stan 컴파일러에서 invalid character 오류가 발생할 수 있다.
      sbv_count = as.numeric(.data[["발생두수"]])
    )

  if (nrow(dat) < 20) {
    stop("완전사례가 너무 적습니다: n = ", nrow(dat))
  }
  if (any(is.na(dat$species))) {
    stop("species에 cerana/mellifera 외 값이 있습니다.")
  }
  if (any(dat$sbv_count <= 0)) {
    stop(
      "zero-truncated 모형에는 발생두수 > 0인 행만 들어가야 합니다. ",
      "0 이하인 행 수: ", sum(dat$sbv_count <= 0)
    )
  }
  if (any(abs(dat$sbv_count - round(dat$sbv_count)) > 1e-8)) {
    stop("발생두수에 정수가 아닌 값이 있습니다. 원래 count를 사용해야 합니다.")
  }
  dat$sbv_count <- as.integer(round(dat$sbv_count))

  message("[양의 발생자료 확인] n = ", nrow(dat))
  message("[양의 발생자료 확인] 최소 발생두수 = ", min(dat$sbv_count))
  message("[양의 발생자료 확인] 0인 관측치 수 = ", sum(dat$sbv_count == 0))

  # 원변수는 유지하고 모형용 z-score 변수를 추가한다.
  scale_info <- vector("list", length(x_raw))
  names(scale_info) <- x_raw

  for (v in x_raw) {
    m <- mean(dat[[v]])
    s <- sd(dat[[v]])
    if (!is.finite(s) || s == 0) {
      stop("분산이 0인 변수: ", v)
    }

    dat[[paste0(v, "_z")]] <- (dat[[v]] - m) / s
    scale_info[[v]] <- c(center = m, scale = s)
  }

  scale_table <- bind_rows(lapply(names(scale_info), function(v) {
    tibble(
      variable = v,
      center = unname(scale_info[[v]]["center"]),
      scale = unname(scale_info[[v]]["scale"])
    )
  }))

  # 기존 코드와 동일하게 EPSG:5181에서 6km 거리 이웃을 만든다.
  pts <- st_as_sf(
    dat,
    coords = c("Longitude", "Latitude"),
    crs = 4326,
    remove = FALSE
  ) %>%
    st_transform(5181)

  coords <- st_coordinates(pts)
  nb <- dnearneigh(coords, d1 = 0, d2 = d2_m)

  n <- nrow(dat)
  W <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    if (length(nb[[i]]) > 0) {
      W[i, nb[[i]]] <- 1
    }
  }
  W <- 1 * ((W + t(W)) > 0)
  diag(W) <- 0

  original_neighbor_count <- rowSums(W)
  isolated <- which(original_neighbor_count == 0)

  # brms의 proper CAR는 모든 위치가 최소 한 개의 이웃을 가져야 안정적으로 적합된다.
  # 원래 6km 이웃이 없는 지점에 한해서만 가장 가까운 양의거리 지점 하나를 연결한다.
  # 6km 이내 이웃이 있던 지점의 관계는 변경하지 않는다.
  fallback_links <- vector("list", length(isolated))

  if (length(isolated) > 0) {
    for (k in seq_along(isolated)) {
      i <- isolated[k]
      d2 <- rowSums((coords - matrix(coords[i, ], nrow = n, ncol = 2, byrow = TRUE))^2)
      d2[i] <- Inf

      # 같은 좌표의 중복 관측치는 거리 0이므로 먼저 제외한다.
      positive <- which(is.finite(d2) & d2 > 0)
      if (length(positive) > 0) {
        j <- positive[which.min(d2[positive])]
      } else {
        # 모든 좌표가 동일한 극단적 경우에만 임의의 다른 행을 연결한다.
        j <- setdiff(seq_len(n), i)[1]
      }

      W[i, j] <- 1
      W[j, i] <- 1

      fallback_links[[k]] <- tibble(
        isolated_id = i,
        linked_id = j,
        distance_m = sqrt(d2[j])
      )
    }
  }

  fallback_table <- if (length(fallback_links)) {
    bind_rows(fallback_links)
  } else {
    tibble(isolated_id = integer(), linked_id = integer(), distance_m = numeric())
  }

  final_neighbor_count <- rowSums(W)
  if (any(final_neighbor_count == 0)) {
    stop("최근접 이웃 보완 후에도 이웃이 없는 지점이 있습니다.")
  }

  # CAR grouping factor의 level과 W의 행·열 이름을 정확히 맞춘다.
  spatial_levels <- as.character(seq_len(n))
  dat$spatial_id <- factor(spatial_levels, levels = spatial_levels)
  dimnames(W) <- list(spatial_levels, spatial_levels)
  storage.mode(W) <- "double"

  neighbor_summary <- tibble(
    n_observations = n,
    distance_band_m = d2_m,
    isolated_under_distance_band = length(isolated),
    nearest_neighbor_fallback_links = nrow(fallback_table),
    min_neighbors_before_fallback = min(original_neighbor_count),
    median_neighbors_before_fallback = median(original_neighbor_count),
    max_neighbors_before_fallback = max(original_neighbor_count),
    min_neighbors_after_fallback = min(final_neighbor_count),
    median_neighbors_after_fallback = median(final_neighbor_count),
    max_neighbors_after_fallback = max(final_neighbor_count)
  )

  message("[6km 이웃 0개였던 관측치] ", length(isolated))
  message("[최근접 이웃 보완 후 이웃 0개] ", sum(final_neighbor_count == 0))

  list(
    data = dat,
    W = W,
    nb = nb,
    scale_table = scale_table,
    scale_info = scale_info,
    isolated = isolated,
    fallback_table = fallback_table,
    neighbor_summary = neighbor_summary
  )
}

spatial_data <- prepare_truncated_car_data(both_full, d2_m = 6000)
dat_zt <- spatial_data$data
W_car <- spatial_data$W

# 결과 파일에서는 원래 한글 변수와 ASCII 분석 변수를 모두 유지한다.
stopifnot(all(dat_zt$sbv_count == dat_zt[["발생두수"]]))


# ======================================================================
# 5) 분리·대체 모형 민감도 분석
# ======================================================================

OUT_DIR <- "out_ztnb_car_sensitivity"
FIT_DIR <- file.path(OUT_DIR, "fits")
TAB_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
TXT_DIR <- file.path(OUT_DIR, "summaries")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TXT_DIR, showWarnings = FALSE, recursive = TRUE)

# 같은 파일명이 존재하면 brms가 저장된 모형을 재사용한다.
# 모형식/자료/사전분포가 바뀌었을 때 자동 재적합하려면 on_change를 사용한다.
FILE_REFIT_MODE <- "on_change"
ADAPT_DELTA <- 0.99
MAX_TREEDEPTH <- 15
N_PRED_DRAWS <- 1000L

# ---------------------------
# 5-1) 일조 지표 상관성 확인
# ---------------------------
sun_cor_pearson <- suppressWarnings(
  cor.test(dat_zt$sum_sun, dat_zt$ratio_sun, method = "pearson")
)
sun_cor_spearman <- suppressWarnings(
  cor.test(dat_zt$sum_sun, dat_zt$ratio_sun, method = "spearman", exact = FALSE)
)

sun_correlation <- tibble(
  comparison = "sum_sun vs ratio_sun",
  n = sum(complete.cases(dat_zt$sum_sun, dat_zt$ratio_sun)),
  pearson_r = unname(sun_cor_pearson$estimate),
  pearson_p = sun_cor_pearson$p.value,
  spearman_rho = unname(sun_cor_spearman$estimate),
  spearman_p = sun_cor_spearman$p.value,
  pairwise_vif_approx = 1 / (1 - unname(sun_cor_pearson$estimate)^2)
)

write.csv(
  sun_correlation,
  file.path(TAB_DIR, "sun_indicator_correlation.csv"),
  row.names = FALSE
)

# ---------------------------
# 5-2) 분석 모형 정의
# ---------------------------
# climate_both: 기존 기후 모형과 동일하게 두 일조 지표를 동시 포함
# climate_sum_only: 누적 일조시간만 포함
# climate_ratio_only: 일조율만 포함
# landcover: 기존 토지피복 모형과 동일한 변수 구성
model_specs <- list(
  climate_both = c(
    "mean_temp_z", "sum_sun_z", "ratio_sun_z", "sum_rain_z"
  ),
  climate_sum_only = c(
    "mean_temp_z", "sum_sun_z", "sum_rain_z"
  ),
  climate_ratio_only = c(
    "mean_temp_z", "ratio_sun_z", "sum_rain_z"
  ),
  landcover = c(
    "forest_z", "agri_z", "water_z"
  )
)

make_car_formula <- function(predictors) {
  rhs <- paste(predictors, collapse = " + ")
  f <- stats::as.formula(
    paste0(
      "sbv_count | trunc(lb = 0) ~ species * (",
      rhs,
      ") + car(W_car, gr = spatial_id, type = 'escar')"
    )
  )
  brms::bf(f, decomp = "QR")
}

intercept_prior_mean <- log(median(dat_zt$sbv_count))

nb_priors <- c(
  brms::set_prior("normal(0, 1)", class = "b"),
  brms::set_prior(
    sprintf("normal(%.12f, 1.5)", intercept_prior_mean),
    class = "Intercept"
  ),
  brms::set_prior("exponential(1)", class = "shape")
)

fit_one_ztnb_car <- function(model_name, predictors, seed) {
  message("\n[적합/불러오기] ", model_name)

  formula_i <- make_car_formula(predictors)
  fit_file <- file.path(FIT_DIR, paste0(model_name, "_ztnb_car"))

  fit <- brms::brm(
    formula = formula_i,
    data = dat_zt,
    data2 = list(W_car = W_car),
    family = brms::negbinomial(link = "log"),
    prior = nb_priors,
    backend = BRMS_BACKEND,
    chains = N_CHAINS,
    iter = N_ITER,
    warmup = N_WARMUP,
    cores = N_CORES,
    seed = seed,
    control = list(
      adapt_delta = ADAPT_DELTA,
      max_treedepth = MAX_TREEDEPTH
    ),
    save_pars = brms::save_pars(all = TRUE),
    file = fit_file,
    file_refit = FILE_REFIT_MODE,
    refresh = 100
  )

  capture.output(
    summary(fit),
    file = file.path(TXT_DIR, paste0(model_name, "_summary.txt"))
  )

  fit
}

fits <- vector("list", length(model_specs))
names(fits) <- names(model_specs)

for (i in seq_along(model_specs)) {
  fits[[i]] <- fit_one_ztnb_car(
    model_name = names(model_specs)[i],
    predictors = model_specs[[i]],
    seed = RANDOM_SEED + 100L + i
  )
}

# ---------------------------
# 5-3) MCMC 수렴 및 sampler 진단
# ---------------------------
model_diagnostics <- function(fit, model_name) {
  d <- posterior::summarise_draws(
    posterior::as_draws_array(fit)
  ) |>
    as.data.frame()

  np <- brms::nuts_params(fit)

  max_rhat <- max(d$rhat, na.rm = TRUE)
  min_bulk <- min(d$ess_bulk, na.rm = TRUE)
  min_tail <- min(d$ess_tail, na.rm = TRUE)
  divergences <- sum(
    np$Parameter == "divergent__" & np$Value == 1,
    na.rm = TRUE
  )
  treedepth_hits <- sum(
    np$Parameter == "treedepth__" & np$Value >= MAX_TREEDEPTH,
    na.rm = TRUE
  )

  tibble(
    model = model_name,
    n = nrow(dat_zt),
    draws = posterior::ndraws(fit),
    max_rhat = max_rhat,
    min_bulk_ess = min_bulk,
    min_tail_ess = min_tail,
    divergences = divergences,
    treedepth_hits = treedepth_hits,
    converged = (
      max_rhat < 1.01 &&
        min_bulk >= 400 &&
        min_tail >= 400 &&
        divergences == 0 &&
        treedepth_hits == 0
    )
  )
}

diagnostics_table <- bind_rows(lapply(names(fits), function(nm) {
  model_diagnostics(fits[[nm]], nm)
}))

write.csv(
  diagnostics_table,
  file.path(TAB_DIR, "model_diagnostics.csv"),
  row.names = FALSE
)

# ---------------------------
# 5-4) 고정효과·IRR 표
# ---------------------------
extract_fixed_effects <- function(fit, model_name) {
  fx <- as.data.frame(
    brms::fixef(fit, probs = c(0.025, 0.975))
  ) |>
    tibble::rownames_to_column("term")

  fx |>
    transmute(
      model = model_name,
      term = term,
      estimate = Estimate,
      posterior_sd = Est.Error,
      lower_95 = Q2.5,
      upper_95 = Q97.5,
      IRR = exp(Estimate),
      IRR_lower_95 = exp(Q2.5),
      IRR_upper_95 = exp(Q97.5),
      credible_95 = Q2.5 > 0 | Q97.5 < 0
    )
}

fixed_effects_table <- bind_rows(lapply(names(fits), function(nm) {
  extract_fixed_effects(fits[[nm]], nm)
}))

write.csv(
  fixed_effects_table,
  file.path(TAB_DIR, "fixed_effects_and_IRR.csv"),
  row.names = FALSE
)

# ---------------------------
# 5-5) 종별 환경효과 결합
# ---------------------------
# 기준종 cerana의 효과 = 환경변수 주효과
# mellifera의 효과 = 주효과 + speciesmellifera 상호작용
# 종 간 기울기 차이 = 상호작용항
find_draw_column <- function(draw_names, candidates, label) {
  hit <- candidates[candidates %in% draw_names]
  if (length(hit) == 0) {
    stop(
      "사후표본에서 계수를 찾을 수 없습니다 [", label, "]: ",
      paste(candidates, collapse = ", ")
    )
  }
  hit[1]
}

summarise_effect_draws <- function(x) {
  q <- unname(stats::quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  tibble(
    estimate = mean(x, na.rm = TRUE),
    median = q[2],
    lower_95 = q[1],
    upper_95 = q[3],
    IRR = exp(mean(x, na.rm = TRUE)),
    IRR_lower_95 = exp(q[1]),
    IRR_upper_95 = exp(q[3]),
    prob_positive = mean(x > 0, na.rm = TRUE),
    prob_negative = mean(x < 0, na.rm = TRUE),
    direction_probability = max(
      mean(x > 0, na.rm = TRUE),
      mean(x < 0, na.rm = TRUE)
    ),
    credible_95 = q[1] > 0 | q[3] < 0
  )
}

extract_species_slopes <- function(fit, model_name, predictors) {
  draws <- posterior::as_draws_df(fit)
  draw_names <- names(draws)

  bind_rows(lapply(predictors, function(v) {
    main_col <- find_draw_column(
      draw_names,
      paste0("b_", v),
      paste0(model_name, ":", v, ":main")
    )

    interaction_col <- find_draw_column(
      draw_names,
      c(
        paste0("b_speciesmellifera:", v),
        paste0("b_", v, ":speciesmellifera")
      ),
      paste0(model_name, ":", v, ":interaction")
    )

    beta_cerana <- draws[[main_col]]
    beta_difference <- draws[[interaction_col]]
    beta_mellifera <- beta_cerana + beta_difference

    bind_rows(
      summarise_effect_draws(beta_cerana) |>
        mutate(species_effect = "cerana"),
      summarise_effect_draws(beta_mellifera) |>
        mutate(species_effect = "mellifera"),
      summarise_effect_draws(beta_difference) |>
        mutate(species_effect = "mellifera_minus_cerana")
    ) |>
      mutate(
        model = model_name,
        predictor = v,
        .before = 1
      )
  }))
}

species_slopes_table <- bind_rows(lapply(names(fits), function(nm) {
  extract_species_slopes(
    fits[[nm]],
    model_name = nm,
    predictors = model_specs[[nm]]
  )
}))

write.csv(
  species_slopes_table,
  file.path(TAB_DIR, "species_specific_environmental_effects.csv"),
  row.names = FALSE
)

# ---------------------------
# 5-6) CAR·과산포 모수
# ---------------------------
extract_distribution_spatial <- function(fit, model_name) {
  draws <- posterior::as_draws_df(fit)
  pars <- intersect(c("car", "sdcar", "shape"), names(draws))

  bind_rows(lapply(pars, function(p) {
    summarise_effect_draws(draws[[p]]) |>
      transmute(
        model = model_name,
        parameter = p,
        estimate,
        median,
        lower_95,
        upper_95,
        prob_positive,
        prob_negative,
        direction_probability,
        credible_95
      )
  }))
}

spatial_distribution_table <- bind_rows(lapply(names(fits), function(nm) {
  extract_distribution_spatial(fits[[nm]], nm)
}))

write.csv(
  spatial_distribution_table,
  file.path(TAB_DIR, "spatial_and_dispersion_parameters.csv"),
  row.names = FALSE
)

# ---------------------------
# 5-7) 일조효과 강건성 표와 자동 판단
# ---------------------------
sun_robustness_table <- species_slopes_table |>
  filter(
    species_effect == "mellifera_minus_cerana",
    predictor %in% c("sum_sun_z", "ratio_sun_z"),
    model %in% c(
      "climate_both",
      "climate_sum_only",
      "climate_ratio_only"
    )
  ) |>
  select(
    model, predictor, estimate, lower_95, upper_95,
    IRR, IRR_lower_95, IRR_upper_95,
    direction_probability, credible_95
  )

write.csv(
  sun_robustness_table,
  file.path(TAB_DIR, "sun_interaction_robustness.csv"),
  row.names = FALSE
)

get_one_sun_result <- function(model_name_value, predictor_value) {
  
  x <- sun_robustness_table |>
    dplyr::filter(
      .data$model == .env$model_name_value,
      .data$predictor == .env$predictor_value
    )
  
  if (nrow(x) != 1L) {
    
    available <- sun_robustness_table |>
      dplyr::select(.data$model, .data$predictor) |>
      dplyr::distinct() |>
      dplyr::arrange(.data$model, .data$predictor)
    
    stop(
      "일조 강건성 행은 정확히 1개여야 합니다: ",
      model_name_value, " / ", predictor_value,
      " (검색된 행 수 = ", nrow(x), ")\n",
      "현재 사용 가능한 조합:\n",
      paste(
        paste0(available$model, " / ", available$predictor),
        collapse = "\n"
      )
    )
  }
  
  x
}

sum_both <- get_one_sun_result("climate_both", "sum_sun_z")
ratio_both <- get_one_sun_result("climate_both", "ratio_sun_z")
sum_single <- get_one_sun_result("climate_sum_only", "sum_sun_z")
ratio_single <- get_one_sun_result("climate_ratio_only", "ratio_sun_z")

sun_decision <- if (
  isTRUE(sum_single$credible_95) &&
    isTRUE(ratio_single$credible_95) &&
    sign(sum_single$estimate) == sign(ratio_single$estimate)
) {
  paste0(
    "두 일조 지표의 단독 모형에서 종 상호작용이 같은 방향으로 95% 신용구간 기준을 충족하였다. ",
    "일조 관련 종별 차이는 지표 선택에 비교적 강건한 것으로 판단할 수 있다."
  )
} else if (
  xor(isTRUE(sum_single$credible_95), isTRUE(ratio_single$credible_95))
) {
  paste0(
    "두 일조 지표 중 하나의 단독 모형에서만 종 상호작용이 95% 신용구간 기준을 충족하였다. ",
    "일조 관련 결론은 지표 선택에 민감하므로, 유의한 결과만 선택하지 말고 사전에 정한 주 지표를 사용해야 한다."
  )
} else if (
  !isTRUE(sum_single$credible_95) && !isTRUE(ratio_single$credible_95)
) {
  paste0(
    "sum_sun 및 ratio_sun 단독 모형 모두에서 종 상호작용의 95% 신용구간이 0을 포함하였다. ",
    "따라서 고일조 조건의 종별 차이를 핵심 결론으로 유지할 근거가 강건하지 않다."
  )
} else {
  paste0(
    "두 단독 모형의 종 상호작용이 모두 95% 신용구간 기준을 충족했지만 방향이 서로 다르다. ",
    "일조효과의 생물학적 해석이 지표에 따라 달라지므로 핵심 결론으로 사용하기 어렵다."
  )
}

simultaneous_model_note <- if (
  sign(sum_both$estimate) != sign(ratio_both$estimate) &&
    !isTRUE(sum_both$credible_95) &&
    !isTRUE(ratio_both$credible_95)
) {
  paste0(
    "두 지표를 동시에 포함한 모형에서 sum_sun과 ratio_sun의 종 상호작용 계수가 반대 방향이고 ",
    "두 신용구간 모두 0을 포함하였다. 공유 정보로 인해 개별 계수의 불확실성이 증가한 패턴과 일치한다."
  )
} else {
  paste0(
    "두 지표 동시 모형의 계수 방향과 신용구간은 sun_interaction_robustness.csv에서 확인한다."
  )
}

sun_decision_text <- c(
  "[일조 지표 상관성]",
  sprintf(
    "Pearson r = %.4f, Spearman rho = %.4f, 근사 VIF = %.3f",
    sun_correlation$pearson_r,
    sun_correlation$spearman_rho,
    sun_correlation$pairwise_vif_approx
  ),
  "",
  "[동시 포함 모형]",
  simultaneous_model_note,
  "",
  "[단독 모형 강건성 판단]",
  sun_decision,
  "",
  "주의: 모형 선택은 유의한 항의 개수가 아니라 자료 구조, 수렴, 사전 지정된 연구 질문을 기준으로 한다."
)

writeLines(
  sun_decision_text,
  file.path(TXT_DIR, "sun_robustness_interpretation.txt"),
  useBytes = TRUE
)

# ---------------------------
# 5-8) posterior_epred를 사용하지 않는 빠른 절단 음이항 예측
# ---------------------------
make_draw_ids <- function(fit, n_draws = N_PRED_DRAWS) {
  total <- posterior::ndraws(fit)
  n_use <- min(as.integer(n_draws), total)
  unique(as.integer(round(seq.int(1, total, length.out = n_use))))
}

fast_ztnb_epred <- function(fit, newdata, draw_ids) {
  mu <- brms::posterior_linpred(
    fit,
    newdata = newdata,
    transform = TRUE,
    re_formula = NA,
    draw_ids = draw_ids
  )

  if (is.null(dim(mu))) {
    mu <- matrix(mu, ncol = 1)
  }

  shape_all <- posterior::as_draws_matrix(fit, variable = "shape")
  if (ncol(shape_all) != 1) {
    stop("shape 사후표본 열을 하나로 특정할 수 없습니다.")
  }
  shape <- as.numeric(shape_all[draw_ids, 1])
  shape_mat <- matrix(
    shape,
    nrow = length(shape),
    ncol = ncol(mu),
    byrow = FALSE
  )

  log_p0 <- shape_mat * (
    log(shape_mat) - log(shape_mat + mu)
  )
  p0 <- exp(log_p0)

  mu / pmax(1 - p0, .Machine$double.eps)
}

make_newdata_template <- function(n) {
  out <- dat_zt[rep(1, n), , drop = FALSE]

  z_vars <- grep("_z$", names(out), value = TRUE)
  for (v in z_vars) out[[v]] <- 0

  out$species <- factor(
    rep("cerana", n),
    levels = levels(dat_zt$species)
  )
  out$spatial_id <- factor(
    rep(levels(dat_zt$spatial_id)[1], n),
    levels = levels(dat_zt$spatial_id)
  )
  out
}

summarise_prediction_columns <- function(x) {
  tibble(
    estimate = colMeans(x, na.rm = TRUE),
    lower_95 = apply(
      x, 2, stats::quantile,
      probs = 0.025, na.rm = TRUE, names = FALSE
    ),
    upper_95 = apply(
      x, 2, stats::quantile,
      probs = 0.975, na.rm = TRUE, names = FALSE
    )
  )
}

make_conditional_curve <- function(fit, focal_var, model_name) {
  focal_values <- seq(
    stats::quantile(dat_zt[[focal_var]], 0.02, na.rm = TRUE),
    stats::quantile(dat_zt[[focal_var]], 0.98, na.rm = TRUE),
    length.out = 80
  )

  key <- tidyr::expand_grid(
    species = factor(
      c("cerana", "mellifera"),
      levels = levels(dat_zt$species)
    ),
    focal_value = focal_values
  )

  nd <- make_newdata_template(nrow(key))
  nd$species <- key$species
  nd[[focal_var]] <- key$focal_value

  draw_ids <- make_draw_ids(fit)
  epred <- fast_ztnb_epred(fit, nd, draw_ids)
  sm <- summarise_prediction_columns(epred)

  out <- bind_cols(
    tibble(
      model = model_name,
      focal_variable = focal_var,
      species = key$species,
      focal_value = key$focal_value
    ),
    sm
  )

  csv_path <- file.path(
    TAB_DIR,
    paste0(model_name, "_", focal_var, "_conditional_curve.csv")
  )
  write.csv(out, csv_path, row.names = FALSE)

  p <- ggplot(
    out,
    aes(
      x = focal_value,
      y = estimate,
      group = species,
      linetype = species
    )
  ) +
    geom_ribbon(
      aes(ymin = lower_95, ymax = upper_95, fill = species),
      alpha = 0.20,
      colour = NA
    ) +
    geom_line(linewidth = 0.9) +
    scale_y_log10() +
    labs(
      x = paste0(focal_var, " (standardized)"),
      y = "Expected SBV count conditional on Y > 0 (log scale)",
      linetype = "species",
      fill = "species",
      title = paste0(model_name, ": ", focal_var)
    ) +
    theme_classic()

  ggsave(
    filename = file.path(
      FIG_DIR,
      paste0(model_name, "_", focal_var, "_conditional_curve_log.png")
    ),
    plot = p,
    width = 8,
    height = 5.5,
    dpi = 300
  )

  out
}

curve_sum <- make_conditional_curve(
  fits$climate_sum_only,
  "sum_sun_z",
  "climate_sum_only"
)

curve_ratio <- make_conditional_curve(
  fits$climate_ratio_only,
  "ratio_sun_z",
  "climate_ratio_only"
)

# ---------------------------
# 5-9) 낮음·평균·높음·상위 10%에서 종별 예측 대비
# ---------------------------
summarise_vector <- function(x) {
  q <- unname(stats::quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  tibble(
    estimate = mean(x, na.rm = TRUE),
    median = q[2],
    lower_95 = q[1],
    upper_95 = q[3]
  )
}

species_contrast_at_values <- function(
    fit,
    focal_var,
    focal_values,
    value_labels,
    model_name
) {
  if (length(focal_values) != length(value_labels)) {
    stop("focal_values와 value_labels의 길이가 다릅니다.")
  }

  draw_ids <- make_draw_ids(fit)

  bind_rows(lapply(seq_along(focal_values), function(i) {
    nd <- make_newdata_template(2)
    nd$species <- factor(
      c("cerana", "mellifera"),
      levels = levels(dat_zt$species)
    )
    nd[[focal_var]] <- focal_values[i]

    ep <- fast_ztnb_epred(fit, nd, draw_ids)
    cerana <- ep[, 1]
    mellifera <- ep[, 2]
    difference <- mellifera - cerana
    ratio <- mellifera / cerana

    bind_cols(
      tibble(
        model = model_name,
        focal_variable = focal_var,
        value_label = value_labels[i],
        focal_value_z = focal_values[i]
      ),
      summarise_vector(cerana) |>
        rename_with(~ paste0("cerana_", .x)),
      summarise_vector(mellifera) |>
        rename_with(~ paste0("mellifera_", .x)),
      summarise_vector(difference) |>
        rename_with(~ paste0("difference_", .x)),
      summarise_vector(ratio) |>
        rename_with(~ paste0("ratio_", .x)),
      tibble(
        prob_mellifera_gt_cerana = mean(mellifera > cerana),
        ratio_credible_above_1 =
          stats::quantile(ratio, 0.025, na.rm = TRUE) > 1,
        difference_credible_above_0 =
          stats::quantile(difference, 0.025, na.rm = TRUE) > 0
      )
    )
  }))
}

sum_q90_z <- unname(stats::quantile(dat_zt$sum_sun_z, 0.90, na.rm = TRUE))
ratio_q90_z <- unname(stats::quantile(dat_zt$ratio_sun_z, 0.90, na.rm = TRUE))

sum_contrast <- species_contrast_at_values(
  fits$climate_sum_only,
  focal_var = "sum_sun_z",
  focal_values = c(-1, 0, 1, sum_q90_z),
  value_labels = c("-1 SD", "mean", "+1 SD", "upper 10% cutoff"),
  model_name = "climate_sum_only"
)

ratio_contrast <- species_contrast_at_values(
  fits$climate_ratio_only,
  focal_var = "ratio_sun_z",
  focal_values = c(-1, 0, 1, ratio_q90_z),
  value_labels = c("-1 SD", "mean", "+1 SD", "upper 10% cutoff"),
  model_name = "climate_ratio_only"
)

sun_species_contrasts <- bind_rows(sum_contrast, ratio_contrast)

write.csv(
  sun_species_contrasts,
  file.path(TAB_DIR, "sun_species_predicted_contrasts.csv"),
  row.names = FALSE
)

sun_cutoffs <- tibble(
  variable = c("sum_sun", "ratio_sun"),
  raw_upper_10_cutoff = c(
    unname(stats::quantile(dat_zt$sum_sun, 0.90, na.rm = TRUE)),
    unname(stats::quantile(dat_zt$ratio_sun, 0.90, na.rm = TRUE))
  ),
  standardized_upper_10_cutoff = c(sum_q90_z, ratio_q90_z)
)

write.csv(
  sun_cutoffs,
  file.path(TAB_DIR, "sun_upper_10_percent_cutoffs.csv"),
  row.names = FALSE
)

# ---------------------------
# 5-10) 논문용 핵심 결과표
# ---------------------------
# 각 분리 모형에서 고정효과만 정리한다.
paper_fixed_effects <- fixed_effects_table |>
  filter(model %in% c("climate_sum_only", "climate_ratio_only", "landcover")) |>
  arrange(model, term)

write.csv(
  paper_fixed_effects,
  file.path(TAB_DIR, "paper_candidate_fixed_effects.csv"),
  row.names = FALSE
)

paper_species_effects <- species_slopes_table |>
  filter(model %in% c("climate_sum_only", "climate_ratio_only", "landcover")) |>
  arrange(model, predictor, species_effect)

write.csv(
  paper_species_effects,
  file.path(TAB_DIR, "paper_candidate_species_effects.csv"),
  row.names = FALSE
)

# ---------------------------
# 5-11) 분석 실행 요약 저장
# ---------------------------
analysis_summary <- c(
  "SBV zero-truncated negative-binomial CAR sensitivity analysis",
  "",
  paste0("Analysis n = ", nrow(dat_zt)),
  paste0("Minimum SBV count = ", min(dat_zt$sbv_count)),
  paste0("Maximum SBV count = ", max(dat_zt$sbv_count)),
  paste0("Backend = ", BRMS_BACKEND),
  paste0("Chains = ", N_CHAINS),
  paste0("Iterations = ", N_ITER),
  paste0("Warmup = ", N_WARMUP),
  paste0("adapt_delta = ", ADAPT_DELTA),
  paste0("max_treedepth = ", MAX_TREEDEPTH),
  "",
  "Models:",
  paste0(
    "- ", names(model_specs), ": ",
    vapply(model_specs, paste, collapse = " + ", FUN.VALUE = character(1))
  ),
  "",
  "Interpretation population:",
  "This analysis concerns SBV count magnitude among occurrence-recorded observations (Y > 0), not the probability of occurrence.",
  "",
  "Sun robustness conclusion:",
  sun_decision
)

writeLines(
  analysis_summary,
  file.path(TXT_DIR, "analysis_summary.txt"),
  useBytes = TRUE
)

# 준비자료도 재사용할 수 있도록 저장한다.
saveRDS(
  list(
    dat_zt = dat_zt,
    W_car = W_car,
    scale_info = spatial_data$scale_info,
    scale_table = spatial_data$scale_table,
    neighbor_summary = spatial_data$neighbor_summary,
    fallback_table = spatial_data$fallback_table
  ),
  file.path(OUT_DIR, "prepared_analysis_data_and_W.rds")
)

cat("\n============================================================\n")
cat("분석 완료\n")
cat("============================================================\n")
cat("\n[수렴 진단]\n")
print(diagnostics_table)
cat("\n[일조 지표 상관]\n")
print(sun_correlation)
cat("\n[일조 상호작용 강건성]\n")
print(sun_robustness_table)
cat("\n[자동 판단]\n")
cat(sun_decision, "\n")
cat("\n결과 폴더: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("주요 파일:\n")
cat("- tables/model_diagnostics.csv\n")
cat("- tables/fixed_effects_and_IRR.csv\n")
cat("- tables/species_specific_environmental_effects.csv\n")
cat("- tables/sun_interaction_robustness.csv\n")
cat("- tables/sun_species_predicted_contrasts.csv\n")
cat("- summaries/sun_robustness_interpretation.txt\n")
cat("- figures/climate_sum_only_sum_sun_z_conditional_curve_log.png\n")
cat("- figures/climate_ratio_only_ratio_sun_z_conditional_curve_log.png\n")
