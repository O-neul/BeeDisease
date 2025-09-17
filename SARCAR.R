setwd("C:/.Soyeon/BeeDisease")

library(dplyr)

########## floral data ##########

# 데이터 불러오기
bee_cluster <- read.csv("data/final/6_Apis_cerana_with_nearby.csv")
floral <- read.csv("data/processed/3_Floral_Source_ID.csv", stringsAsFactors = FALSE)

# 컬럼명 통일 및 재정렬
colnames(floral) <- tolower(colnames(floral))
floral <- floral %>%
  rename(
    species = 종,
    latitude = 위도,
    longitude = 경도,
    abundance = 카운트
  )

# 거리 계산 함수: 각 꿀벌 군집 위치에 대해 6km 이내 밀원 찾기
find_nearby_floral_ids <- function(lat, lon, floral_df, threshold_km = 6) {
  bee_coord <- matrix(c(lon, lat), ncol = 2)
  floral_coords <- floral_df %>% select(longitude, latitude) %>% as.matrix()
  
  dists <- distHaversine(bee_coord, floral_coords) / 1000  # m → km
  nearby_ids <- floral_df$id[dists <= threshold_km]
  return(nearby_ids)
}

# 꿀벌 군집별 nearby 밀원 ID 리스트 생성
bee_cluster$nearby_floral_ids <- lapply(1:nrow(bee_cluster), function(i) {
  find_nearby_floral_ids(
    lat = bee_cluster$Latitude[i],
    lon = bee_cluster$Longitude[i],
    floral_df = floral,
    threshold_km = 6
  )
})

# 시각화
library(geosphere)
library(sf)
library(ggplot2)
library(purrr)

# 1. sf 객체로 변환
floral_subset <- floral %>% filter(id %in% unlist(bee_cluster$nearby_floral_ids))
bee_subset <- bee_cluster %>% filter(lengths(nearby_floral_ids) > 0)

# sf 객체로 변환
floral_sf <- st_as_sf(floral_subset, coords = c("longitude", "latitude"), crs = 4326)
bee_sf <- st_as_sf(bee_subset, coords = c("Longitude", "Latitude"), crs = 4326)

connection_lines <- map_dfr(1:nrow(bee_subset), function(i) {
  bee_lat <- bee_subset$Latitude[i]
  bee_lon <- bee_subset$Longitude[i]
  ids <- bee_subset$nearby_floral_ids[[i]]
  
  # 해당 군집에 연결된 밀원만 필터링
  linked_florals <- floral %>% filter(id %in% ids)
  
  if (nrow(linked_florals) == 0) return(NULL)
  
  # 각 연결을 LINESTRING으로 생성
  lines <- lapply(1:nrow(linked_florals), function(j) {
    coords <- matrix(c(
      bee_lon, bee_lat,
      linked_florals$longitude[j], linked_florals$latitude[j]
    ), ncol = 2, byrow = TRUE)
    
    st_linestring(coords)
  })
  
  # sf 객체로 변환
  sf::st_sf(
    geometry = sf::st_sfc(lines),
    crs = 4326
  )
})

# sf 포인트 객체 생성
bee_sf <- st_as_sf(bee_subset, coords = c("Longitude", "Latitude"), crs = 4326)
floral_sf <- st_as_sf(floral_subset, coords = c("longitude", "latitude"), crs = 4326)

# 시각화
ggplot() +
  geom_sf(data = connection_lines, color = "gray70", size = 0.3, alpha = 0.7) +
  geom_sf(data = floral_sf, color = "darkgreen", size = 1.5, alpha = 0.8) +
  geom_sf(data = bee_sf, color = "blue", size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "6km 이내 밀원 연결 시각화",
       subtitle = "파랑: 꿀벌 군집 / 초록: 밀원 / 회색선: 연결선")

########## weather data cleaning ##########

library(dplyr)
library(lubridate)
library(geosphere)  # 거리 계산용

# 데이터 불러오기
weather <- read.csv("data/raw/weather.csv", fileEncoding = "CP949")
stations <- read.csv("data/raw/weather_station_info.csv", fileEncoding = "CP949")
# 날짜 파싱 및 연/월 추출
weather$일시 <- ymd(weather$일시)
sum(is.na(weather$일시))  # 실패 개수 확인
weather <- weather[!is.na(weather$일시), ]
weather$year <- year(weather$일시)
weather$month <- month(weather$일시)
# 2018년만 필터링
weather_2018 <- weather %>% filter(year == 2018)
# 계절 구분
weather_2018$season <- case_when(
  weather_2018$month %in% c(12, 1, 2) ~ "winter",
  weather_2018$month %in% c(3, 4, 5)  ~ "spring",
  weather_2018$month %in% c(6, 7, 8)  ~ "summer",
  weather_2018$month %in% c(9, 10, 11) ~ "fall"
)
# 계절별 평균기온 계산 (지점별)
seasonal_avg_temp <- weather_2018 %>%
  group_by(지점, season) %>%
  summarise(mean_temp = mean(평균기온..C., na.rm=TRUE)) %>%
  ungroup()

# 여름만 필터링
summer_temp <- seasonal_avg_temp %>% filter(season == "summer")

# 관측소 위도/경도 정리
stations_clean <- stations %>%
  select(지점, 위도, 경도) %>%
  rename(station_id = 지점, lat = 위도, lon = 경도)

# --------- 군집 정보 불러오기 ---------
bee_cluster <- read.csv("data/final/6_Apis_cerana_with_nearby.csv")

# 각 군집에서 가장 가까운 관측소 매칭
get_nearest_station <- function(lat, lon, stations) {
  dists <- distHaversine(matrix(c(lon, lat), ncol = 2),
                         matrix(c(stations$lon, stations$lat), ncol = 2))
  stations$station_id[which.min(dists)]
}

bee_cluster$nearest_station <- apply(
  bee_cluster[, c("Latitude", "Longitude")],
  1,
  function(row) {
    get_nearest_station(as.numeric(row["Latitude"]), as.numeric(row["Longitude"]), stations_clean)
  }
)

# 여름 기온 데이터만 병합
bee_cluster_temp <- merge(bee_cluster, summer_temp,
                          by.x = "nearest_station", by.y = "지점")
# 결과 예시 확인
head(bee_cluster_temp)

########## floral data cleaning ##########

# 선호 종 목록 (cerana 기준)
cerana_preferred <- c('밤나무', '산밤나무', '약밤나무', '배추', '유채', '머귀나무',
                      '왕초피', '초피나무', '산초나무', '벼', '피마자',
                      '다래', '쥐다래', '개다래', '섬다래',
                      '감나무', '고욤나무', '개옺나무', '검양옺나무', '산검양옺나무',
                      '옺나무', '붉나무', '광대싸리')
# 밀원 군집 데이터 (3_Floral_Source_ID.csv)
floral <- read.csv("data/processed/3_Floral_Source_ID.csv", stringsAsFactors = FALSE)
names(floral) <- tolower(names(floral))  # 컬럼명 정규화 (예: id, species, abundance)
floral <- floral %>%
  rename(
    species = 종,
    latitude = 위도,
    longitude = 경도,
    abundance = 카운트
  )
# 데이터 점검
stopifnot(all(c("id", "species", "abundance") %in% names(floral)))
# nearby_* 컬럼 추출
nearby_cols <- grep("^nearby_", names(bee_cluster_temp), value = TRUE)
# 파생 변수 계산 함수
get_flower_stats <- function(row, floral_df, prefer_list) {
  ids <- as.character(row[nearby_cols])
  ids <- ids[!is.na(ids) & ids != ""]
  
  matched <- floral_df %>% filter(id %in% ids)
  
  preferred <- matched %>% filter(species %in% prefer_list)
  nonpreferred <- matched %>% filter(!species %in% prefer_list)
  
  return(data.frame(
    prefer_count = nrow(preferred),
    prefer_abund = sum(preferred$abundance, na.rm = TRUE),
    nonprefer_count = nrow(nonpreferred),
    nonprefer_abund = sum(nonpreferred$abundance, na.rm = TRUE)
  ))
}
# 모든 행에 대해 반복 적용
flower_vars <- bind_rows(
  apply(bee_cluster_temp, 1, get_flower_stats, floral_df = floral, prefer_list = cerana_preferred)
)
# 최종 병합
cluster_final <- bind_cols(bee_cluster_temp, flower_vars)

########## neighbor ##########

library(sf)
library(spdep)

# 여름 군집만 필터링
summer_data <- cluster_final[cluster_final$season == "summer", ]

# 좌표 추출
coords <- summer_data[, c("Longitude", "Latitude")]

sf_points <- st_as_sf(summer_data, coords = c("Longitude", "Latitude"), crs = 4326)
sf_points_proj <- st_transform(sf_points, crs = 5181)
proj_coords <- st_coordinates(sf_points_proj)

# 4. dnearneigh로 0~6000m 이웃 계산
nb <- dnearneigh(proj_coords, d1 = 0, d2 = 6000)

# 5. 검토
summary(nb)

########## 인접행렬 시각화 ##########

# 좌표 매트릭스 추출
coord_matrix <- st_coordinates(sf_points_proj)

# 시각화
plot(st_geometry(sf_points_proj), pch = 20, cex = 0.6, col = "darkgray", main = "꿀벌 군집 인접 네트워크 (6km)")
plot(nb, coord_matrix, add = TRUE, col = "blue", lwd = 1)

########## SAR modelling ##########

library(spatialreg)

# 1. 인접 리스트 → 가중치 행렬 생성
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# 2. 종속변수와 설명변수 지정
y <- cluster_final$발생두수
X <- cluster_final[, c("prefer_count", "mean_temp")]

# 3. SAR 모델 적합
sar_model <- lagsarlm(y ~ ., data = X, listw = listw, zero.policy = TRUE)

# 4. 결과 확인
summary(sar_model)

########## CAR modelling ##########

listw_sym <- nb2listw(nb, style = "B", zero.policy = TRUE)

car_model <- spautolm(y ~ prefer_count + prefer_abund +
                        nonprefer_count + nonprefer_abund + mean_temp,
                      data = cluster_final,
                      listw = listw_sym,
                      family = "CAR",  # <-- 핵심
                      zero.policy = TRUE)

summary(car_model)
