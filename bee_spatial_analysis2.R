setwd("C:/Users/haa/Downloads/BeeDisease/BeeDisease")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(sf)
})

# -----------------------------
# 0) 경로 설정 — 환경에 맞게 수정
# -----------------------------
PATH_LANDCOVER_CSV <- "data/raw/환경부_환경공간정보_년도별 대분류토지피복통계 현황_20230901.csv"
PATH_SIG_BOUNDARY  <- "data/raw/BND_SIGUNGU_PG.shp"

# -----------------------------
# 1) 토지피복 집계 로더 (연도 맞춤)
# -----------------------------
load_landcover_stats <- function(path_csv = PATH_LANDCOVER_CSV, year_keep = 2018) {
  readr::read_csv(path_csv, locale = locale(encoding = "UTF-8")) %>%
    select(자료년도, 국가구분, 시도, 시군구,
           시가화건조지역, 농업지역, 산림지역, 초지, 습지, 나지, 수역, 합계) %>%
    filter(자료년도 == year_keep) %>%
    mutate(across(c(시도, 시군구), ~ str_squish(as.character(.x)))) %>%
    mutate(
      urban_share   = 100 * 시가화건조지역/합계,
      agri_share    = 100 * 농업지역/합계,
      forest_share  = 100 * 산림지역/합계,
      grass_share   = 100 * 초지/합계,
      wetland_share = 100 * 습지/합계,
      barren_share  = 100 * 나지/합계,
      water_share   = 100 * 수역/합계
    ) %>%
    select(시도, 시군구, ends_with("_share"))
}

lc_stats <- load_landcover_stats(PATH_LANDCOVER_CSV, year_keep = 2018)

# 시군구만 남겨 중복 해소(여러 시도에 같은 시군구명이 있을 때 평균 or 합계를 선택)
# 보통 비율변수이니 평균이 합리적입니다.
lc_stats_by_sig <- lc_stats %>%
  dplyr::group_by(시군구) %>%
  dplyr::summarise(dplyr::across(dplyr::ends_with("_share"), mean, na.rm = TRUE), .groups = "drop")


# -----------------------------
# 2) 시군구 경계 로더 (컬럼명 표준화)
#   - 당신 경계 파일의 실제 컬럼명에 맞춰 rename만 조정하면 됨
# -----------------------------
load_sig_boundary <- function(path_sig = PATH_SIG_BOUNDARY) {
  sig <- st_read(path_sig, quiet = TRUE) %>% st_make_valid()
  sig %>%
    rename(
      SIG_NM = SIGUNGU_NM  # 시군구 이름
    ) %>%
    mutate(
      SIG_NM = stringr::str_squish(as.character(SIG_NM))
    )
}

# -----------------------------
# 3) 포인트 데이터에 시군구 라벨 & landcover 조인
# -----------------------------
attach_landcover_to_points <- function(df_points,  # Longitude, Latitude 필요
                                       lon_col = "Longitude", lat_col = "Latitude",
                                       sig_sf, lc_stats_by_sig) {
  stopifnot(all(c(lon_col, lat_col) %in% names(df_points)))
  
  pts <- sf::st_as_sf(df_points, coords = c(lon_col, lat_col), crs = 4326)
  if (sf::st_crs(sig_sf) != sf::st_crs(pts)) {
    pts <- sf::st_transform(pts, sf::st_crs(sig_sf))
  }
  
  # sig_sf에서 존재하는 열만 선택(여기서는 SIG_NM만 있을 가능성 높음)
  keep_cols <- intersect(c("SIDO_NM","SIG_NM"), names(sig_sf))
  sig_min   <- sig_sf[, c(keep_cols, "geometry")]
  
  # 포인트-인-폴리곤
  pts2 <- sf::st_join(pts, sig_min, left = TRUE)
  
  out <- pts2 |>
    sf::st_drop_geometry()
  
  # 시군구명 만들기
  if ("SIG_NM" %in% names(out)) {
    out <- out |> dplyr::mutate(시군구 = SIG_NM) |> dplyr::select(-SIG_NM)
  } else {
    stop("경계 데이터에 SIG_NM(=시군구명) 컬럼이 없습니다.")
  }
  
  # landcover(시군구 단위) 조인
  out <- out |>
    dplyr::left_join(lc_stats_by_sig, by = "시군구")
  
  out
}


# -----------------------------
# 4) 실행: 기존 데이터에 landcover 열만 '추가'
#   - cerana_summer / mellifera_summer / both 가 이미 있음
# -----------------------------
# (1) 입력 데이터 확인
required_objs <- c("cerana_summer","mellifera_summer","both")
missing <- required_objs[!vapply(required_objs, exists, logical(1))]
if (length(missing)) {
  stop(sprintf("다음 객체가 메모리에 없습니다: %s\n기존 스크립트를 먼저 실행해 %s를 만들어 주세요.",
               paste(missing, collapse=", "), paste(required_objs, collapse=", ")))
}

# (2) 참조 데이터 로드
sig_sf   <- load_sig_boundary(PATH_SIG_BOUNDARY)
lc_stats <- load_landcover_stats(PATH_LANDCOVER_CSV, year_keep = 2018)

# (3) 조인 실행 — 여기서부터 결과물만 새로 생김
cerana_summer_lc    <- attach_landcover_to_points(cerana_summer,    sig_sf = sig_sf, lc_stats = lc_stats)
mellifera_summer_lc <- attach_landcover_to_points(mellifera_summer, sig_sf = sig_sf, lc_stats = lc_stats)
both_lc             <- attach_landcover_to_points(both,              sig_sf = sig_sf, lc_stats = lc_stats)

# -----------------------------
# 5) 체크: landcover 열 존재 여부/결측률 간단 점검
# -----------------------------
message("[CHECK] landcover 변수 열 머리:")
print(intersect(names(both_lc), c("urban_share","agri_share","forest_share","grass_share",
                                  "wetland_share","barren_share","water_share")))

message("[CHECK] 결측률(%) — both_lc:")
na_rate <- colMeans(is.na(both_lc[, c("시도","시군구",
                                      "urban_share","agri_share","forest_share","grass_share",
                                      "wetland_share","barren_share","water_share")]))*100
print(round(na_rate, 2))










library(dplyr)
library(sf)
library(geosphere)
library(lubridate)
library(spdep)
library(spatialreg)

# ======================================================================
# [공용함수]
# ======================================================================

# 1-1. Floral 로딩/정리
load_floral <- function(path = "data/processed/3_Floral_Source_ID.csv") {
  floral <- read.csv(path, stringsAsFactors = FALSE) %>%
    rename_with(tolower) %>%
    rename(
      species   = 종,
      latitude  = 위도,
      longitude = 경도,
      abundance = 카운트
    )
  stopifnot(all(c("id","species","latitude","longitude","abundance") %in% names(floral)))
  floral
}

# 1-2. 날씨 요약 (여름 평균기온, 관측소 정보와 조인)
prepare_summer_temp <- function(weather_path = "data/raw/weather.csv",
                                stations_path = "data/raw/weather_station_info.csv",
                                year_keep = 2018) {
  weather  <- read.csv(weather_path, fileEncoding = "CP949")
  stations <- read.csv(stations_path, fileEncoding = "CP949")
  weather$일시 <- ymd(weather$일시)
  weather <- weather %>% filter(!is.na(일시)) %>%
    mutate(year = year(일시), month = month(일시)) %>%
    filter(year == {{year_keep}}) %>%
    mutate(season = case_when(
      month %in% c(12,1,2) ~ "winter",
      month %in% c(3,4,5)  ~ "spring",
      month %in% c(6,7,8)  ~ "summer",
      TRUE                 ~ "fall"
    ))
  
  seasonal_avg_temp <- weather %>%
    group_by(지점, season) %>%
    summarise(
      mean_temp = mean(평균기온..C., na.rm = TRUE),
      sum_rain  = mean(일최다강수량.mm., na.rm = TRUE),
      sum_sun  = mean(합계일조시간.hr., na.rm = TRUE),
      ratio_sun  = mean(일조율, na.rm = TRUE),
      .groups = "drop") %>%
    filter(season == "summer") %>%
    rename(station_id = 지점)
  
  stations_clean <- stations %>%
    select(지점, 위도, 경도) %>%
    rename(station_id = 지점, lat = 위도, lon = 경도)
  
  list(summer_temp = seasonal_avg_temp, stations = stations_clean)
}

# 1-3. 한 지점에서 가장 가까운 관측소 ID
nearest_station_id <- function(lat, lon, stations_df) {
  d <- distHaversine(matrix(c(lon, lat), ncol = 2),
                     matrix(c(stations_df$lon, stations_df$lat), ncol = 2))
  stations_df$station_id[which.min(d)]
}

# 1-4. 벌 군집 주변 6km 밀원 ID (필요 시)
nearby_floral_ids <- function(lat, lon, floral_df, threshold_km = 6) {
  d <- distHaversine(matrix(c(lon, lat), ncol = 2),
                     as.matrix(floral_df[, c("longitude","latitude")])) / 1000
  floral_df$id[d <= threshold_km]
}

# 1-5. 선호/비선호 요약 변수 생성
flower_stats <- function(ids, floral_df, prefer_list) {
  ids <- ids[!is.na(ids) & ids != ""]
  matched <- floral_df %>% filter(id %in% ids)
  
  preferred    <- matched %>% filter(species %in% prefer_list)
  nonpreferred <- matched %>% filter(!species %in% prefer_list)
  
  tibble(
    prefer_count     = nrow(preferred),
    prefer_abund     = sum(preferred$abundance, na.rm = TRUE),
    nonprefer_count  = nrow(nonpreferred),
    nonprefer_abund  = sum(nonpreferred$abundance, na.rm = TRUE)
  )
}

# 1-6. 이웃/가중치 구성 (0~6km)
build_neighbors <- function(df, lon_col = "Longitude", lat_col = "Latitude",
                            crs_proj = 5181, d2_m = 6000, style = "W") {
  sf_pts  <- st_as_sf(df, coords = c(lon_col, lat_col), crs = 4326)
  sf_proj <- st_transform(sf_pts, crs = crs_proj)
  coords  <- st_coordinates(sf_proj)
  nb      <- dnearneigh(coords, d1 = 0, d2 = d2_m)
  listw   <- nb2listw(nb, style = style, zero.policy = TRUE)
  list(sf_pts = sf_pts, sf_proj = sf_proj, nb = nb, listw = listw)
}

# 1-7. SAR / CAR 적합 도우미
fit_sar <- function(y, X, listw) {
  lagsarlm(formula = y ~ ., data = X, listw = listw, zero.policy = TRUE)
}
fit_car <- function(formula, data, listw_sym) {
  spautolm(formula, data = data, listw = listw_sym, family = "CAR", zero.policy = TRUE)
}

# ======================================================================
# [종별 데이터]
# ======================================================================

# cerana 기준 선호종
cerana_preferred <- c('밤나무','산밤나무','약밤나무','배추','유채','머귀나무','왕초피','초피나무',
                      '산초나무','벼','피마자','다래','쥐다래','개다래','섬다래','감나무','고욤나무',
                      '개옺나무','검양옺나무','산검양옺나무','옺나무','붉나무','광대싸리')

mellifera_preferred <- c( '족제비싸리', '아까지나무', '갈참나무', '굴참나무', '떡갈나무', '종가시나무',
                          '물참나무', '상수리나무', '신갈나무', '졸참나무', '가시나무', '붉가시나무',
                          '갈졸참나무', '떡갈참나무', '떡신갈나무', '주름잎', '누운주름잎', '덩굴장미',
                          '목향장미', '생열귀나무', '용가시나무', '인가목', '해당화', '고추', '애기똥풀',
                          '다래', '쥐다래', '개다래', '섬다래', '감나무', '고욤나무', '개옺나무', '검양옺나무',
                          '산검양옺나무', '옺나무', '붉나무', '광대싸리' )

# 종별 데이터 전처리(재사용 가능)
prepare_bee_dataset <- function(bee_path, floral_df, summer_temp, stations_df,
                                prefer_list = cerana_preferred,
                                threshold_km = 6) {
  
  bee <- read.csv(bee_path)
  
  # 주변 밀원 ID가 파일에 이미 있으면 그걸 쓰고, 없으면 계산
  if (!"nearby_floral_ids" %in% names(bee)) {
    bee$nearby_floral_ids <- lapply(1:nrow(bee), function(i) {
      nearby_floral_ids(
        lat = bee$Latitude[i],
        lon = bee$Longitude[i],
        floral_df = floral_df,
        threshold_km = threshold_km
      )
    })
  }
  
  # 가장 가까운 관측소
  bee$nearest_station <- apply(
    bee[, c("Latitude","Longitude")], 1, \(row) {
      nearest_station_id(as.numeric(row["Latitude"]), as.numeric(row["Longitude"]), stations_df)
    }
  )
  
  # 여름 평균기온 결합
  bee2 <- bee %>%
    left_join(summer_temp, by = c("nearest_station" = "station_id"))
  
  # 시즌 정보가 별도 없으면 여름만으로 태깅(필요 시 수정)
  if (!"season" %in% names(bee2)) bee2$season <- "summer"
  
  # 선호/비선호 파생변수
  nearby_cols <- grep("^nearby_", names(bee2), value = TRUE)
  
  flower_vars <- dplyr::bind_rows(
    lapply(seq_len(nrow(bee2)), function(i) {
      # nearby_*가 하나의 list 열(nearby_floral_ids)인 경우
      if (length(nearby_cols) == 1 && nearby_cols[1] == "nearby_floral_ids") {
        ids <- bee2$nearby_floral_ids[[i]]
      } else {
        # nearby_1, nearby_2 ... 형태를 허용
        vals <- as.vector(unlist(bee2[i, nearby_cols], use.names = FALSE))
        ids  <- vals[!is.na(vals) & vals != ""]
      }
      flower_stats(ids, floral_df, prefer_list)
    })
  )
  
  bind_cols(bee2, flower_vars)
}

# ======================================================================
# [종별/전체 모델링]
# ======================================================================

# 3-1. 공용 리소스
floral <- load_floral()
wx <- prepare_summer_temp()
summer_temp  <- wx$summer_temp
stations_tbl <- wx$stations

# 3-2. 종별 데이터 준비
df_cerana    <- prepare_bee_dataset("data/final/6_Apis_cerana_with_nearby.csv",
                                    floral, summer_temp, stations_tbl,
                                    prefer_list = cerana_preferred, threshold_km = 6) %>%
  mutate(species = "cerana")

df_mellifera <- prepare_bee_dataset("data/final/6_Apis_mellifera_with_nearby.csv",
                                    floral, summer_temp, stations_tbl,
                                    # mellifera 선호종이 다르면 전용 리스트를 넣으세요.
                                    prefer_list = mellifera_preferred, threshold_km = 6) %>%
  mutate(species = "mellifera")

# 3-3. 여름만 필터 및 결측치 제거
cerana_summer    <- df_cerana    %>% filter(season == "summer") %>%
  tidyr::drop_na(발생두수, prefer_count, mean_temp, sum_sun, ratio_sun, sum_rain)

mellifera_summer <- df_mellifera %>% filter(season == "summer") %>%
  tidyr::drop_na(발생두수, prefer_count, mean_temp, sum_sun, ratio_sun, sum_rain)


# 3-4. 이웃/가중치 생성(종별)
g_c  <- build_neighbors(cerana_summer)
g_m  <- build_neighbors(mellifera_summer)

# 3-5. SAR (종별)
y_c <- cerana_summer$발생두수
X_c <- cerana_summer[, c("prefer_count", "sum_sun")]  # 변수셋 교체 가능
sar_c <- fit_sar(y_c, X_c, g_c$listw)
summary(sar_c)

y_m <- mellifera_summer$발생두수
X_m <- mellifera_summer[, c("prefer_count", "ratio_sun","sum_rain")]
sar_m <- fit_sar(y_m, X_m, g_m$listw)
summary(sar_m)

# 3-6. CAR (종별)
listw_sym_c <- nb2listw(g_c$nb, style = "B", zero.policy = TRUE)
listw_sym_m <- nb2listw(g_m$nb, style = "B", zero.policy = TRUE)

car_c <- fit_car(
  발생두수 ~ prefer_count +
    mean_temp + sum_sun + ratio_sun + sum_rain,
  data = cerana_summer, listw_sym = listw_sym_c
)
summary(car_c)

car_m <- fit_car(
  발생두수 ~ prefer_count +
    sum_sun  + sum_rain,
  data = mellifera_summer, listw_sym = listw_sym_m
)
summary(car_m)

# 3-7. 두 종 결합 + 합동 모델
both <- bind_rows(cerana_summer, mellifera_summer) %>%
  mutate(species = factor(species)) %>%
  tidyr::drop_na(발생두수, prefer_count, mean_temp, sum_sun, ratio_sun, sum_rain)

# 결합 데이터에서 이웃은 "혼합 네트워크"로 구성(두 종 위치가 한 평면에 같이 존재)
g_both <- build_neighbors(both)
listw_both    <- g_both$listw
listw_sym_both<- nb2listw(g_both$nb, style = "B", zero.policy = TRUE)

# SAR (합동) — 종 효과 + 상호작용까지
X_both <- both %>%
  transmute(prefer_count, mean_temp, species)

sar_both <- lagsarlm(
  발생두수 ~ (prefer_count + mean_temp + sum_sun + ratio_sun + sum_rain)*species,
  data = both, listw = listw_both, zero.policy = TRUE
)
summary(sar_both)

# CAR (합동)
car_both <- spautolm(
  발생두수 ~ (prefer_count +
            mean_temp + sum_sun + ratio_sun + sum_rain)*species,
  data = both, listw = listw_sym_both, family = "CAR", zero.policy = TRUE
)
summary(car_both)











library(dplyr)
library(sf)
library(spdep)
library(spatialreg)

# 1) 모형에 쓸 변수 목록 (water_share는 다중공선성 회피 위해 제외)
vars_needed <- c(
  "발생두수","prefer_count","mean_temp","sum_sun","ratio_sun","sum_rain",
  "urban_share","agri_share","forest_share","grass_share","wetland_share","barren_share",
  "species", "Longitude","Latitude"
)

# 2) 완전자료만 필터
both_lc_clean <- both_lc %>%
  select(any_of(vars_needed)) %>%
  tidyr::drop_na()

# species 팩터 레벨 원본과 동일하게 고정(상호작용 설계행렬 안정)
if (is.factor(both$species)) {
  both_lc_clean$species <- factor(both_lc_clean$species, levels = levels(both$species))
} else {
  both_lc_clean$species <- factor(both_lc_clean$species)
}

# 3) 이 "깨끗한" 행들로 인접/가중치 다시 생성
g_both_lc <- build_neighbors(
  df      = both_lc_clean,
  lon_col = "Longitude",
  lat_col = "Latitude",
  crs_proj = 5181,   # 기존과 동일
  d2_m     = 6000,
  style    = "W"
)
listw_sym_both_lc <- nb2listw(g_both_lc$nb, style = "B", zero.policy = TRUE)

# 4) CAR 적합 (상호작용 포함)
car_both_lc <- spautolm(
  발생두수 ~ (prefer_count +
            mean_temp + sum_sun + ratio_sun + sum_rain +
            urban_share + agri_share + forest_share +
            grass_share + wetland_share + barren_share) * species,
  data = both_lc_clean,
  listw = listw_sym_both_lc,
  family = "CAR",
  zero.policy = TRUE,
  na.action = na.exclude   # 안전장치 (잔차/적합치 정렬 유지)
)
summary(car_both_lc)
