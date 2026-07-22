# ======================================================================
# SBV 발생두수: 발생 지점(Y>0) 대상 zero-truncated Poisson/NB CAR 모형
# - 기존 프로젝트에 실제로 존재하는 파일 경로와 변수명을 사용함
# - 반응변수: 발생두수 (원자료 count; 로그변환하지 않음)
# - 기후변수: mean_temp, sum_sun, ratio_sun, sum_rain
# - 토지피복: forest, agri, water (urban은 기준범주로 제외)
# - 종: species (cerana 기준, mellifera 비교)
# - 공간구조: EPSG:5181, 0 < 거리 <= 6,000m, brms exact sparse CAR
# - 발생 기록이 있는 지점만 포함하므로 0에서 절단된 count likelihood를 사용
# - Stan 컴파일 안정성을 위해 반응변수는 모형 내부에서 sbv_count로 이름을 바꿈
# - 제공된 코드에 조사봉군수/노출량 변수가 없으므로 offset은 사용하지 않음
# ======================================================================

setwd("C:/.Soyeon/DS/BeeDisease")
setwd("C:/Users/haa/Downloads/BeeDisease/BeeDisease")

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

# ---------------------------
# 4-3) zero-truncated Poisson/NB CAR 적합
# ---------------------------
# urban은 토지피복 구성비의 기준범주로 제외한다.
# decomp = "QR"은 상관이 높은 설명변수로 인한 샘플링 불안정을 줄인다.
formula_zt <- bf(
  sbv_count | trunc(lb = 0) ~
    species * (
      mean_temp_z + sum_sun_z + ratio_sun_z + sum_rain_z +
      forest_z + agri_z + water_z
    ) +
    car(W_car, gr = spatial_id, type = "escar"),
  decomp = "QR"
)

# prior() 안에 median(dat_zt$...) 같은 R 식을 직접 넣으면,
# brms가 그 식을 Stan 코드 문자열로 그대로 넘길 수 있다.
# 특히 한글 변수명이 포함되면 Stan lexer에서 invalid character 오류가 난다.
# 따라서 R에서 먼저 숫자를 계산한 뒤 ASCII 문자열 prior로 전달한다.
intercept_prior_mean <- log(median(dat_zt$sbv_count))

common_priors <- c(
  set_prior("normal(0, 1)", class = "b"),
  set_prior(
    sprintf("normal(%.12f, 1.5)", intercept_prior_mean),
    class = "Intercept"
  )
)

nb_priors <- c(
  common_priors,
  set_prior("exponential(1)", class = "shape")
)

message(
  "[Intercept prior] normal(",
  format(intercept_prior_mean, digits = 6),
  ", 1.5)"
)

# 실제 사용 가능한 prior class를 실행 전에 확인하고 파일로 저장한다.
dir.create("out_zt_car_brms", showWarnings = FALSE)

prior_check_poisson <- get_prior(
  formula_zt,
  data = dat_zt,
  data2 = list(W_car = W_car),
  family = poisson(link = "log")
)

prior_check_nb <- get_prior(
  formula_zt,
  data = dat_zt,
  data2 = list(W_car = W_car),
  family = negbinomial(link = "log")
)

write.csv(
  prior_check_poisson,
  "out_zt_car_brms/available_priors_zero_truncated_poisson.csv",
  row.names = FALSE
)
write.csv(
  prior_check_nb,
  "out_zt_car_brms/available_priors_zero_truncated_nb.csv",
  row.names = FALSE
)

# 컴파일 전에 생성될 Stan 코드를 확인·저장한다.
# 이 단계가 통과하면 prior나 formula에 한글 식별자가 Stan 코드로 유입되지 않은 것이다.
stan_code_poisson <- stancode(
  formula_zt,
  data = dat_zt,
  data2 = list(W_car = W_car),
  family = poisson(link = "log"),
  prior = common_priors
)
stan_code_nb <- stancode(
  formula_zt,
  data = dat_zt,
  data2 = list(W_car = W_car),
  family = negbinomial(link = "log"),
  prior = nb_priors
)

writeLines(
  stan_code_poisson,
  "out_zt_car_brms/generated_zero_truncated_poisson.stan",
  useBytes = TRUE
)
writeLines(
  stan_code_nb,
  "out_zt_car_brms/generated_zero_truncated_nb.stan",
  useBytes = TRUE
)

if (grepl("발생두수", stan_code_poisson, fixed = TRUE) ||
    grepl("발생두수", stan_code_nb, fixed = TRUE)) {
  stop("생성된 Stan 코드에 한글 반응변수명이 남아 있습니다.")
}

message("[Stan 코드 사전검사 완료] ASCII 반응변수 sbv_count 사용")

message("[적합 시작] Zero-truncated Poisson CAR")
fit_ztp <- brm(
  formula = formula_zt,
  data = dat_zt,
  data2 = list(W_car = W_car),
  family = poisson(link = "log"),
  prior = common_priors,
  backend = BRMS_BACKEND,
  chains = N_CHAINS,
  iter = N_ITER,
  warmup = N_WARMUP,
  cores = N_CORES,
  seed = RANDOM_SEED,
  control = list(adapt_delta = 0.97, max_treedepth = 13),
  save_pars = save_pars(all = TRUE),
  file = "out_zt_car_brms/zero_truncated_poisson_car_fit"
)

message("[적합 시작] Zero-truncated Negative-binomial CAR")
fit_ztnb <- brm(
  formula = formula_zt,
  data = dat_zt,
  data2 = list(W_car = W_car),
  family = negbinomial(link = "log"),
  prior = nb_priors,
  backend = BRMS_BACKEND,
  chains = N_CHAINS,
  iter = N_ITER,
  warmup = N_WARMUP,
  cores = N_CORES,
  seed = RANDOM_SEED + 1,
  control = list(adapt_delta = 0.97, max_treedepth = 13),
  save_pars = save_pars(all = TRUE),
  file = "out_zt_car_brms/zero_truncated_negative_binomial_car_fit"
)

# ---------------------------
# 4-4) LOO/WAIC와 계수표
# ---------------------------
fit_ztp <- add_criterion(fit_ztp, criterion = c("loo", "waic"))
fit_ztnb <- add_criterion(fit_ztnb, criterion = c("loo", "waic"))

extract_model_criteria <- function(fit, model_name) {
  loo_obj <- fit$criteria$loo
  waic_obj <- fit$criteria$waic

  tibble(
    model = model_name,
    elpd_loo = loo_obj$estimates["elpd_loo", "Estimate"],
    elpd_loo_se = loo_obj$estimates["elpd_loo", "SE"],
    LOOIC = -2 * loo_obj$estimates["elpd_loo", "Estimate"],
    LOOIC_se = 2 * loo_obj$estimates["elpd_loo", "SE"],
    WAIC = waic_obj$estimates["waic", "Estimate"],
    WAIC_se = waic_obj$estimates["waic", "SE"],
    pareto_k_over_0_7 = sum(loo_obj$diagnostics$pareto_k > 0.7, na.rm = TRUE),
    pareto_k_over_1 = sum(loo_obj$diagnostics$pareto_k > 1, na.rm = TRUE)
  )
}

model_compare <- bind_rows(
  extract_model_criteria(fit_ztp, "Zero-truncated Poisson CAR"),
  extract_model_criteria(fit_ztnb, "Zero-truncated Negative-binomial CAR")
) %>%
  arrange(LOOIC)

loo_pairwise <- as.data.frame(
  loo::loo_compare(
    list(
      Zero_truncated_Poisson_CAR = fit_ztp$criteria$loo,
      Zero_truncated_NB_CAR = fit_ztnb$criteria$loo
    )
  )
) %>%
  tibble::rownames_to_column("model")

print(loo_pairwise)

make_fixed_effect_table <- function(fit, model_name) {
  fx <- as.data.frame(fixef(fit, probs = c(0.025, 0.975))) %>%
    rownames_to_column("term")

  fx %>%
    transmute(
      model = model_name,
      term = term,
      mean = Estimate,
      sd = Est.Error,
      `0.025quant` = Q2.5,
      `0.975quant` = Q97.5,
      IRR = exp(Estimate),
      IRR_lower = exp(Q2.5),
      IRR_upper = exp(Q97.5)
    )
}

fixed_effects <- bind_rows(
  make_fixed_effect_table(fit_ztp, "Zero-truncated Poisson CAR"),
  make_fixed_effect_table(fit_ztnb, "Zero-truncated Negative-binomial CAR")
)

# ---------------------------
# 4-5) 발생 기록 조건부 예측평균
# ---------------------------
make_prediction_table <- function(fit, prefix) {
  pred <- as.data.frame(fitted(
    fit,
    scale = "response",
    summary = TRUE,
    probs = c(0.025, 0.5, 0.975)
  ))

  expected_names <- c("Estimate", "Est.Error", "Q2.5", "Q50", "Q97.5")
  if (ncol(pred) != length(expected_names)) {
    stop("fitted() 결과 열 수가 예상과 다릅니다: ", paste(names(pred), collapse = ", "))
  }
  names(pred) <- paste0(prefix, c(
    "_pred_mean", "_pred_sd", "_pred_lower", "_pred_median", "_pred_upper"
  ))
  pred
}

pred_ztp <- make_prediction_table(fit_ztp, "ztp")
pred_ztnb <- make_prediction_table(fit_ztnb, "ztnb")

prediction_data <- bind_cols(dat_zt, pred_ztp, pred_ztnb)

best_model <- model_compare$model[1]
if (best_model == "Zero-truncated Poisson CAR") {
  prediction_data$selected_pred_mean <- prediction_data$ztp_pred_mean
} else {
  prediction_data$selected_pred_mean <- prediction_data$ztnb_pred_mean
}

p_pred <- ggplot(prediction_data, aes(x = species, y = selected_pred_mean)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  labs(
    title = paste0(best_model, ": 종별 예측 SBV 발생두수"),
    subtitle = "SBV 발생 기록이 있는 지점(Y>0)에 조건부인 예측평균",
    x = "꿀벌 종",
    y = "예측 발생두수 E(Y | Y>0)"
  ) +
  theme_minimal()

print(p_pred)

ggsave(
  "out_zt_car_brms/selected_model_species_prediction.png",
  p_pred,
  width = 8,
  height = 5,
  dpi = 300
)

# ---------------------------
# 5) 결과 저장
# ---------------------------
write.csv(
  fixed_effects,
  "out_zt_car_brms/zero_truncated_car_fixed_effects_IRR.csv",
  row.names = FALSE
)

write.csv(
  model_compare,
  "out_zt_car_brms/zero_truncated_poisson_vs_nb_car.csv",
  row.names = FALSE
)

write.csv(
  loo_pairwise,
  "out_zt_car_brms/zero_truncated_poisson_vs_nb_loo_compare.csv",
  row.names = FALSE
)

write.csv(
  spatial_data$scale_table,
  "out_zt_car_brms/standardization_information.csv",
  row.names = FALSE
)

write.csv(
  spatial_data$neighbor_summary,
  "out_zt_car_brms/spatial_neighbor_summary.csv",
  row.names = FALSE
)

write.csv(
  spatial_data$fallback_table,
  "out_zt_car_brms/isolated_site_nearest_neighbor_links.csv",
  row.names = FALSE
)

write.csv(
  prediction_data %>%
    select(
      발생두수, species, Longitude, Latitude,
      mean_temp, sum_sun, ratio_sun, sum_rain,
      forest, agri, water,
      starts_with("ztp_pred_"),
      starts_with("ztnb_pred_"),
      selected_pred_mean
    ),
  "out_zt_car_brms/zero_truncated_car_predictions.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    analysis_population = "SBV occurrence-recorded sites only (발생두수 > 0)",
    response = "발생두수",
    likelihood_poisson = "Poisson truncated below or equal to 0 via brms trunc(lb=0)",
    likelihood_negative_binomial = "Negative binomial truncated below or equal to 0 via brms trunc(lb=0)",
    fixed_effects = paste(
      "species * (mean_temp_z + sum_sun_z + ratio_sun_z + sum_rain_z +",
      "forest_z + agri_z + water_z)"
    ),
    neighbor_definition = "0 < distance <= 6000 m in EPSG:5181",
    isolated_site_handling = paste0(
      "Only sites with zero neighbors under 6 km were linked to their nearest positive-distance site; n=",
      length(spatial_data$isolated)
    ),
    spatial_model = "brms exact sparse CAR (escar)",
    backend = BRMS_BACKEND,
    chains = N_CHAINS,
    iter = N_ITER,
    warmup = N_WARMUP,
    seed_poisson = RANDOM_SEED,
    seed_negative_binomial = RANDOM_SEED + 1
  ),
  "out_zt_car_brms/model_specification.csv",
  row.names = FALSE
)

capture.output(
  summary(fit_ztp),
  file = "out_zt_car_brms/zero_truncated_poisson_car_summary.txt"
)

capture.output(
  summary(fit_ztnb),
  file = "out_zt_car_brms/zero_truncated_negative_binomial_car_summary.txt"
)

saveRDS(
  fit_ztp,
  "out_zt_car_brms/zero_truncated_poisson_car_fit.rds"
)

saveRDS(
  fit_ztnb,
  "out_zt_car_brms/zero_truncated_negative_binomial_car_fit.rds"
)

cat("\n[모형 비교: LOOIC/WAIC가 낮을수록 우수]\n")
print(model_compare)

cat("\n[LOO 직접 비교]\n")
print(loo_pairwise)

cat("\n[고정효과와 IRR]\n")
print(fixed_effects)

cat("\n[해석 범위]\n")
cat("이 모형은 SBV가 발생할 확률을 분석하지 않습니다.\n")
cat("SBV 발생이 기록된 지점들 중에서 발생두수가 얼마나 큰지를 분석합니다.\n")
cat("IRR은 다른 조건과 CAR 공간효과를 통제했을 때 기저 count 강도의 배수효과입니다.\n")
cat("최종 모형은 out_zt_car_brms/zero_truncated_poisson_vs_nb_car.csv에서 LOOIC가 낮은 모형입니다.\n")
