setwd("C:/Users/haa/Downloads/BeeDisease/BeeDisease")

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

# ======================================================================
# [일조량 상위 집단 분석]
# ======================================================================

suppressPackageStartupMessages({
  if (!"ggplot2" %in% (.packages())) library(ggplot2)
  if (!"dplyr"   %in% (.packages())) library(dplyr)
})

stopifnot(exists("both"), exists("car_both"))

# -----------------------------
# 0) 헬퍼: 적합 terms로 디자인행렬 생성 + 열 정렬/보강
# -----------------------------
make_X_from_model <- function(fit, newdata) {
  b <- coef(fit)              # 회귀계수(람다 제외), named numeric
  ff <- formula(fit)
  
  # 1) 디자인행렬 생성 (수치형 matrix 보장)
  X <- model.matrix(ff, data = newdata)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  # 2) (Intercept) 보강: b에는 있는데 X에 없으면 1-열을 'matrix'로 추가
  if ("(Intercept)" %in% names(b) && !("(Intercept)" %in% colnames(X))) {
    X <- cbind("(Intercept)" = rep(1, nrow(X)), X)
    X <- as.matrix(X)
    storage.mode(X) <- "double"
  }
  
  # 3) 빠진 열(주로 상호작용/더미) 0-열을 'matrix'로 추가
  missing_cols <- setdiff(names(b), colnames(X))
  if (length(missing_cols) > 0) {
    Z <- matrix(0, nrow = nrow(X), ncol = length(missing_cols))
    colnames(Z) <- missing_cols
    X <- cbind(X, Z)
    X <- as.matrix(X)
    storage.mode(X) <- "double"
  }
  
  # 4) 계수 순서대로 정렬 + 최종 수치형 matrix 보장
  X <- X[, names(b), drop = FALSE]
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  # 5) 계수도 수치형 벡터로
  b_vec <- as.numeric(b)      # 이름 순서는 위에서 맞췄으므로 OK
  
  list(X = X, b = b_vec)
}

# -----------------------------
# 1) 일조량 상위 집단 정의
# -----------------------------
q_top  <- 0.90  # 상위 10% 컷 (필요시 조정)
thr_sun <- quantile(both$sum_sun, q_top, na.rm = TRUE)
message(sprintf("[INFO] sum_sun %.0f%% 컷오프 = %.3f", (1 - q_top) * 100, thr_sun))

high_sun <- both %>% filter(sum_sun >= thr_sun)

# 팩터 레벨 고정(모델과 동일 레벨 유지)
high_sun$species <- factor(high_sun$species, levels = levels(both$species))

message(sprintf("[INFO] 상위 집단 크기: %d / %d (%.1f%%)",
                nrow(high_sun), nrow(both), 100 * nrow(high_sun) / nrow(both)))
print(table(high_sun$species))

# -----------------------------
# 2) 관측값 기준 종별 요약 + 탐색적 검정
# -----------------------------
obs_stats <- high_sun %>%
  group_by(species) %>%
  summarise(
    n      = n(),
    mean_y = mean(발생두수, na.rm = TRUE),
    sd_y   = sd(발생두수, na.rm = TRUE),
    med_y  = median(발생두수, na.rm = TRUE),
    q25_y  = quantile(발생두수, 0.25, na.rm = TRUE),
    q75_y  = quantile(발생두수, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
print(obs_stats)

if (length(unique(high_sun$species)) == 2) {
  message("[탐색적] t-test:")
  print(try(t.test(발생두수 ~ species, data = high_sun), silent = TRUE))
  message("[탐색적] Wilcoxon:")
  print(try(wilcox.test(발생두수 ~ species, data = high_sun, exact = FALSE), silent = TRUE))
}

# 관측값 분포
ggplot(high_sun, aes(x = species, y = 발생두수)) +
  geom_boxplot() +
  labs(title = sprintf("관측값 분포 (sum_sun 상위 %.0f%%, 임계=%.2f)", (1 - q_top) * 100, thr_sun),
       x = "종", y = "발생두수")

# -----------------------------
# 3) CAR 모형 기반 예측값 (E[y|X] = Xβ)
# -----------------------------
tmp <- make_X_from_model(car_both, high_sun)
high_sun$yhat <- drop(tmp$X %*% tmp$b)

pred_stats <- high_sun %>%
  group_by(species) %>%
  summarise(
    n        = n(),
    mean_hat = mean(yhat, na.rm = TRUE),
    sd_hat   = sd(yhat,   na.rm = TRUE),
    med_hat  = median(yhat, na.rm = TRUE),
    q25_hat  = quantile(yhat, 0.25, na.rm = TRUE),
    q75_hat  = quantile(yhat, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
print(pred_stats)

# 예측값 분포
ggplot(high_sun, aes(x = species, y = yhat)) +
  geom_boxplot() +
  labs(title = sprintf("모형 예측값 분포 (sum_sun 상위 %.0f%%, 임계=%.2f)", (1 - q_top) * 100, thr_sun),
       x = "종", y = "예측 발생두수 (CAR, Xβ)")

# -----------------------------
# 4) Counterfactual: 같은 관측에서 species만 바꿔보기
# -----------------------------
lvl <- levels(both$species)
high_sun_cer <- high_sun; high_sun_cer$species <- factor("cerana",    levels = lvl)
high_sun_mel <- high_sun; high_sun_mel$species <- factor("mellifera", levels = lvl)

tmp_cer <- make_X_from_model(car_both, high_sun_cer)
tmp_mel <- make_X_from_model(car_both, high_sun_mel)

high_sun$yhat_cer_cf <- drop(tmp_cer$X %*% tmp_cer$b)
high_sun$yhat_mel_cf <- drop(tmp_mel$X %*% tmp_mel$b)
high_sun$delta_species <- high_sun$yhat_mel_cf - high_sun$yhat_cer_cf

delta_summary <- high_sun %>%
  summarise(
    n         = n(),
    mean_diff = mean(delta_species, na.rm = TRUE),
    med_diff  = median(delta_species, na.rm = TRUE),
    q25_diff  = quantile(delta_species, 0.25, na.rm = TRUE),
    q75_diff  = quantile(delta_species, 0.75, na.rm = TRUE)
  )
print(delta_summary)

ggplot(high_sun, aes(x = "", y = delta_species)) +
  geom_boxplot() +
  coord_flip() +
  labs(title = "Counterfactual: 같은 셋에서 species만 변경 (mellifera - cerana)",
       x = NULL, y = "예측 차이")

# -----------------------------
# 5) 일조량 한계효과(계수) 확인
# -----------------------------
b <- coef(car_both)
me_cerana     <- unname(b["sum_sun"])
me_mellifera  <- unname(b["sum_sun"] + b["sum_sun:speciesmellifera"])
message(sprintf("[일조량 한계효과] cerana=%.3f, mellifera=%.3f", me_cerana, me_mellifera))

# -----------------------------
# 6) 시각화
# -----------------------------

# delta_species 히스토그램
ggplot(high_sun, aes(x = delta_species)) +
  geom_histogram(binwidth = 5, color = "black", fill = "skyblue") +
  labs(title = "Counterfactual 차이 분포 (mellifera - cerana)",
       x = "예측 차이", y = "빈도")

# delta_species 밀도곡선
ggplot(high_sun, aes(x = delta_species)) +
  geom_density(fill = "orange", alpha = 0.4) +
  labs(title = "Counterfactual 차이 분포 (mellifera - cerana)",
       x = "예측 차이", y = "밀도")

# -----------------------------
# 7) (선택) 결과 저장
# -----------------------------
# dir.create("out", showWarnings = FALSE)
# write.csv(obs_stats,  "out/highsun_obs_stats.csv",  row.names = FALSE)
# write.csv(pred_stats, "out/highsun_pred_stats.csv", row.names = FALSE)
# write.csv(high_sun[, c("species","발생두수","yhat","yhat_cer_cf","yhat_mel_cf","delta_species")],
#           "out/highsun_rowwise_preds.csv", row.names = FALSE)