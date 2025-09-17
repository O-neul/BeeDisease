setwd("C:/.Soyeon/BeeDisease")
setwd("C:/Users/haa/Downloads/BeeDisease/BeeDisease")

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(geosphere)
  library(lubridate)
  library(spdep)
  library(spatialreg)
  library(ggplot2)
  library(readr)
})

# ---------------------------------------------------------------
# 1. 데이터 로드
# ---------------------------------------------------------------

# 1-1. 밀원
floral <- read.csv("data/processed/3_Floral_Source_ID.csv", stringsAsFactors = FALSE) %>%
  rename_with(tolower) %>%
  rename(
    species   = 종,
    latitude  = 위도,
    longitude = 경도,
    abundance = 카운트
  )

# 1-2. 날씨
weather  <- read.csv("data/raw/weather.csv", fileEncoding = "CP949")
stations <- read.csv("data/raw/weather_station_info.csv", fileEncoding = "CP949")

# 1-3. 행정구역 (시군구 경계)
sig_sf <- st_read("data/raw/BND_SIGUNGU_PG.shp", options = "ENCODING=CP949")

# 1-4. 토지피복 (환경부 제공 CSV)
landcover <- read.delim(
  "data/raw/환경부_환경공간정보_년도별 대분류토지피복통계 현황_20230901.csv",
  fileEncoding = "CP949"
)

print(names(landcover))
str(landcover)
head(landcover)

# ===============================================================
# Landcover CSV "복구 로더" + 매핑 + both CAR 모델 (새 파일)
# ===============================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(janitor)
  library(sf); library(spdep); library(spatialreg); library(ggplot2)
})

# ---------- 0) 유틸: 구분자 추정 ----------
guess_delim_from_lines <- function(x, n = 10) {
  lines <- head(x, n)
  cand  <- c("," = ",", "\t" = "\t", ";" = ";", "|" = "|")
  counts <- sapply(cand, function(d) sum(str_count(lines, fixed(d))))
  names(which.max(counts))[1]
}

# ---------- 1) 유틸: 탄력적 CSV/텍스트 로더 ----------
load_table_resilient <- function(path) {
  if (!file.exists(path)) stop(sprintf("파일 없음: %s", path))
  message(sprintf("[INFO] 파일 크기: %s bytes", file.size(path)))
  
  # 1) 원시 바이트 + 인코딩 추정
  raw <- read_file_raw(path)
  enc_guess <- try(guess_encoding(raw[1:min(length(raw), 200000)]), silent = TRUE)
  enc <- "CP949"
  if (!inherits(enc_guess, "try-error") && nrow(enc_guess) > 0) {
    enc <- enc_guess$encoding[1]
  }
  message(sprintf("[INFO] 인코딩 추정: %s", enc))
  
  # 2) 라인 단위로 읽고 마지막 줄 보정(미완성 줄 → 개행 추가)
  lines <- read_lines(path, locale = locale(encoding = enc), progress = FALSE)
  if (length(lines) == 0) {
    message("[WARN] 라인이 0개 → 다른 구분자/인코딩 재시도 필요")
  }
  if (length(lines) > 0 && !str_detect(lines[length(lines)], "\r$|\n$")) {
    lines[length(lines)] <- paste0(lines[length(lines)], "\n")
  }
  
  # 3) 구분자 추정
  delim <- guess_delim_from_lines(lines, n = 20)
  message(sprintf("[INFO] 구분자 추정: '%s'", ifelse(delim == "\t", "\\t", delim)))
  
  # 4) 임시 파일로 정리해 읽기
  tf <- tempfile(fileext = ".csv")
  write_lines(lines, tf)
  
  # 5) readr로 시도
  df <- try(read_delim(tf, delim = delim, locale = locale(encoding = enc),
                       na = c("", "NA", "-", "NaN"), guess_max = 100000), silent = TRUE)
  if (!inherits(df, "try-error") && nrow(df) > 0 && ncol(df) > 1) {
    message("[OK] read_delim 성공")
    return(df)
  }
  
  # 6) data.table::fread 백업
  if (requireNamespace("data.table", quietly = TRUE)) {
    message("[INFO] readr 실패 → data.table::fread 시도")
    df2 <- try(data.table::fread(tf, encoding = enc, sep = delim, na.strings = c("", "NA", "-", "NaN")), silent = TRUE)
    if (!inherits(df2, "try-error") && nrow(df2) > 0 && ncol(df2) > 1) {
      return(as.data.frame(df2))
    }
  }
  
  # 7) 엑셀로 저장된 파일 가능성 체크
  if (requireNamespace("readxl", quietly = TRUE)) {
    message("[INFO] CSV 실패 → 엑셀 시도(readxl)")
    fmt <- try(readxl::excel_format(path), silent = TRUE)
    if (!inherits(fmt, "try-error") && !is.na(fmt) && fmt != "unknown") {
      df3 <- try(readxl::read_excel(path), silent = TRUE)
      if (!inherits(df3, "try-error") && nrow(df3) > 0) {
        return(as.data.frame(df3))
      }
    }
  }
  
  stop("테이블 로드 실패: 인코딩/구분자/파일형식 문제")
}

# ---------- 2) 토지피복 로드 + 컬럼 자동 매핑 ----------
load_landcover_and_map <- function(path_landcover) {
  lc <- load_table_resilient(path_landcover) %>% as.data.frame()
  message("[INFO] 원본 열이름:")
  print(names(lc))
  
  # 열이름 정리(한/영/공백/특수문자 → 스네이크케이스)
  lc2 <- lc %>% janitor::clean_names()
  
  # 시군구 코드 후보: sigungu_cd, sig_cd, adm_cd, lawd_cd, ... (한글도 대응)
  code_col <- names(lc2)[str_detect(names(lc2),
                                    "(sigungu|si_gun_gu|sgg|adm|lawd|code|cd|코드)")]
  if (length(code_col) == 0) {
    # 이름 컬럼으로 붙일 수도 있으니, 시군구명 후보도 따로 보관
    name_col <- names(lc2)[str_detect(names(lc2), "(sigungu|si_gun_gu|sgg|시군구|시_군_구|행정구역).*?(nm|name|명)?")]
  } else {
    name_col <- NULL
  }
  message(sprintf("[INFO] 코드 후보: %s", paste(code_col, collapse = ", ")))
  if (!is.null(name_col)) message(sprintf("[INFO] 명칭 후보: %s", paste(name_col, collapse = ", ")))
  
  # 대분류 토지피복 후보(한/영 키워드)
  forest_cols <- names(lc2)[str_detect(names(lc2), "(forest|산림|임야|초지)")]
  agri_cols   <- names(lc2)[str_detect(names(lc2), "(agri|농업|경작|논|밭|과수)")]
  urban_cols  <- names(lc2)[str_detect(names(lc2), "(urban|시가|도시|공업|상업|주거|교통)")]
  water_cols  <- names(lc2)[str_detect(names(lc2), "(water|수계|하천|호소|습지|저수지)")]
  message("[INFO] 추정된 토지피복 컬럼:")
  print(list(forest = forest_cols, agri = agri_cols, urban = urban_cols, water = water_cols))
  
  # 비율/면적이 여러 열이면 평균(또는 합계)로 축약
  summariser <- function(.df, cols) {
    if (length(cols) == 0) return(NA_real_)
    rowMeans(as.matrix(.df[, cols, drop = FALSE]), na.rm = TRUE)
  }
  
  # 키 식별: 코드 우선, 없으면 명칭으로.
  if (length(code_col) >= 1) {
    key <- code_col[1]
    lc_out <- lc2 %>%
      mutate(
        forest = summariser(., forest_cols),
        agri   = summariser(., agri_cols),
        urban  = summariser(., urban_cols),
        water  = summariser(., water_cols)
      ) %>%
      group_by(.data[[key]]) %>%
      summarise(
        forest = mean(forest, na.rm = TRUE),
        agri   = mean(agri,   na.rm = TRUE),
        urban  = mean(urban,  na.rm = TRUE),
        water  = mean(water,  na.rm = TRUE),
        .groups = "drop"
      ) %>%
      rename(SIGUNGU_CD = all_of(key))
  } else if (!is.null(name_col) && length(name_col) >= 1) {
    key <- name_col[1]
    lc_out <- lc2 %>%
      mutate(
        forest = summariser(., forest_cols),
        agri   = summariser(., agri_cols),
        urban  = summariser(., urban_cols),
        water  = summariser(., water_cols)
      ) %>%
      group_by(.data[[key]]) %>%
      summarise(
        forest = mean(forest, na.rm = TRUE),
        agri   = mean(agri,   na.rm = TRUE),
        urban  = mean(urban,  na.rm = TRUE),
        water  = mean(water,  na.rm = TRUE),
        .groups = "drop"
      ) %>%
      rename(SIGUNGU_NM = all_of(key))
  } else {
    stop("시군구 코드/명칭 컬럼을 찾지 못했습니다. names(lc2)를 확인해 매핑 규칙을 추가하세요.")
  }
  
  lc_out
}

# ---------- 3) 시군구 경계 로드 ----------
load_sigungu_sf <- function(path_shp = "data/raw/BND_SIGUNGU_PG.shp") {
  sf <- st_read(path_shp, options = "ENCODING=CP949", quiet = TRUE)
  sf
}

# ---------- 4) bee 두 종 결합 + 토지피복 join ----------
# ---- 교체용 함수: CRS 정렬 + 유효성 보정 포함 ----
build_both_with_landcover <- function(bee_c_path, bee_m_path,
                                      landcover_tbl, sigungu_sf) {
  # 1) bee 읽기
  bee_c <- read.csv(bee_c_path)
  bee_m <- read.csv(bee_m_path)
  
  # 2) sf 변환 (WGS84로 만듦)
  c_sf <- st_as_sf(bee_c, coords = c("Longitude","Latitude"), crs = 4326, remove = FALSE)
  m_sf <- st_as_sf(bee_m, coords = c("Longitude","Latitude"), crs = 4326, remove = FALSE)
  
  # 3) 시군구 유효성 보정 + CRS 정렬
  sgg <- sigungu_sf
  if (!all(st_is_valid(sgg))) sgg <- st_make_valid(sgg)
  if (is.na(st_crs(sgg))) stop("시군구 shp에 CRS가 없습니다. st_set_crs로 올바른 EPSG를 지정하세요.")
  
  c_sf <- st_transform(c_sf, st_crs(sgg))
  m_sf <- st_transform(m_sf, st_crs(sgg))
  
  # 4) 조인 컬럼 선택
  sig_cols <- names(sgg)
  has_cd <- any(grepl("(SIGUNGU|SGG).*?(CD|code|cd|코드)", sig_cols, ignore.case = TRUE))
  has_nm <- any(grepl("(SIGUNGU|SGG|시군구).*?(NM|NAME|명)", sig_cols, ignore.case = TRUE))
  pick_cols <- character(0)
  if (has_cd) pick_cols <- c(pick_cols, sig_cols[grepl("(SIGUNGU|SGG).*?(CD|code|cd|코드)", sig_cols, ignore.case = TRUE)][1])
  if (has_nm) pick_cols <- c(pick_cols, sig_cols[grepl("(SIGUNGU|SGG|시군구).*?(NM|NAME|명)", sig_cols, ignore.case = TRUE)][1])
  if (length(pick_cols) == 0) pick_cols <- sig_cols[1]
  
  # 5) 공간조인
  sgg2 <- suppressWarnings(st_buffer(sgg, 0))
  c_join <- st_join(c_sf, sgg2[, pick_cols, drop = FALSE], join = st_intersects, left = TRUE)
  m_join <- st_join(m_sf, sgg2[, pick_cols, drop = FALSE], join = st_intersects, left = TRUE)
  
  c_df <- st_drop_geometry(c_join)
  m_df <- st_drop_geometry(m_join)
  
  # 6) landcover 키 결정 (코드 우선, 없으면 명칭)
  if ("SIGUNGU_CD" %in% names(landcover_tbl) &&
      any(grepl("CD|code|cd|코드", names(c_df), ignore.case = TRUE))) {
    key_lc <- "SIGUNGU_CD"
    key_sf <- names(c_df)[grepl("CD|code|cd|코드", names(c_df), ignore.case = TRUE)][1]
  } else if ("SIGUNGU_NM" %in% names(landcover_tbl) &&
             any(grepl("NM|NAME|명", names(c_df), ignore.case = TRUE))) {
    key_lc <- "SIGUNGU_NM"
    key_sf <- names(c_df)[grepl("NM|NAME|명", names(c_df), ignore.case = TRUE)][1]
  } else {
    stop("시군구 키 매칭 실패: landcover_tbl과 시군구 조인 컬럼을 확인하세요.")
  }
  
  c_df2 <- dplyr::left_join(c_df, landcover_tbl, by = setNames(key_lc, key_sf))
  m_df2 <- dplyr::left_join(m_df, landcover_tbl, by = setNames(key_lc, key_sf))
  
  both <- dplyr::bind_rows(
    c_df2 %>% dplyr::mutate(species = "cerana"),
    m_df2 %>% dplyr::mutate(species = "mellifera")
  )
  
  # 7) 진단 메시지
  if (any(is.na(both[[key_sf]]))) {
    message("[WARN] 시군구 조인이 일부 NA입니다. 좌표/경계 확인 요망.")
  }
  if (any(is.na(both$forest) | is.na(both$agri) | is.na(both$urban) | is.na(both$water))) {
    message("[WARN] landcover NA 존재. 키/연도/코드 불일치 확인 요망.")
  }
  both
}

# ---------- 5) 공간가중행렬 + both CAR ----------
fit_both_car_with_lc <- function(both_df) {
  stopifnot(is.data.frame(both_df))
  df <- both_df
  
  # 1) 좌표 열 자동 탐색 (영문/한글/대소문자)
  guess_lon <- names(df)[grepl("^(longitude|lon|long|x|경도)$", names(df), ignore.case = TRUE)][1]
  guess_lat <- names(df)[grepl("^(latitude|lat|y|위도)$", names(df), ignore.case = TRUE)][1]
  if (is.na(guess_lon) || is.na(guess_lat)) {
    stop("경도/위도 열을 찾지 못했습니다. (예: Longitude/Latitude). names(df)로 실제 열 이름을 확인하세요.")
  }
  
  # 2) 좌표 숫자형 강제 + 결측/비정상 제거
  df[[guess_lon]] <- suppressWarnings(as.numeric(df[[guess_lon]]))
  df[[guess_lat]] <- suppressWarnings(as.numeric(df[[guess_lat]]))
  before_n <- nrow(df)
  df <- df %>%
    filter(is.finite(.data[[guess_lon]]), is.finite(.data[[guess_lat]]))
  if (nrow(df) < before_n) {
    message(sprintf("[INFO] 좌표 비정상/결측 %d행 제거", before_n - nrow(df)))
  }
  
  # 3) 필수 기본 변수 확인
  if (!("발생두수" %in% names(df))) {
    stop("종속변수 '발생두수'가 데이터에 없습니다.")
  }
  if (!("species" %in% names(df))) {
    stop("'species' 변수가 없습니다. build_both_with_landcover() 과정 확인.")
  }
  
  # 4) landcover 후보와 기후/밀원 후보 중 실제 존재하는 것만 사용
  lc_candidates   <- c("forest","agri","urban","water")
  clim_candidates <- c("prefer_count","mean_temp","sum_sun","ratio_sun","sum_rain")
  
  lc_used   <- intersect(lc_candidates,   names(df))
  clim_used <- intersect(clim_candidates, names(df))
  
  if (length(lc_used) == 0) {
    stop("토지피복 변수(forest/agri/urban/water)가 하나도 없습니다. landcover 조인 단계를 확인하세요.")
  }
  
  # 5) 동적 포뮬러 구성: 발생두수 ~ (clim + lc)*species
  rhs_terms <- c(lc_used, clim_used, "species")
  rhs       <- paste(rhs_terms, collapse = " + ")
  fml_txt   <- paste0("발생두수 ~ (", paste(c(lc_used, clim_used), collapse = " + "), ")*species")
  fml       <- as.formula(fml_txt)
  
  message("[MODEL] 사용 변수: ", paste(rhs_terms, collapse = ", "))
  message("[MODEL] 공식: ", fml_txt)
  
  # 6) 분석용 데이터 구성: 필요한 열만 추출 + 결측 제거
  need_for_drop <- unique(c("발생두수", "species", lc_used, clim_used, guess_lon, guess_lat))
  dat <- df[, need_for_drop, drop = FALSE] %>%
    tidyr::drop_na("발생두수", "species", all_of(lc_used), all_of(guess_lon), all_of(guess_lat)) %>%
    mutate(species = factor(species))
  
  if (nrow(dat) < 10) {
    stop(sprintf("유효 행이 너무 적습니다(n=%d). 좌표/조인/결측을 확인하세요.", nrow(dat)))
  }
  
  # 7) sf 변환 → 투영(5181) → 이웃행렬
  sf_pts  <- st_as_sf(dat, coords = c(guess_lon, guess_lat), crs = 4326)
  sf_proj <- suppressWarnings(st_transform(sf_pts, 5181))
  coords  <- st_coordinates(sf_proj)
  if (!is.numeric(coords[,1]) || !is.numeric(coords[,2])) {
    stop("좌표가 숫자가 아닙니다. 원본 좌표 열을 확인하세요.")
  }
  
  nb  <- dnearneigh(coords, d1 = 0, d2 = 6000)
  lwB <- nb2listw(nb, style = "B", zero.policy = TRUE)
  
  # 8) CAR 적합
  fit <- spautolm(fml, data = dat, listw = lwB, family = "CAR", zero.policy = TRUE)
  
  # 9) 결과 + 한계효과(landcover만)
  b <- coef(fit)
  me_fun <- function(base) {
    base_eff <- unname(b[base])
    int_nm   <- paste0(base, ":speciesmellifera")
    int_eff  <- ifelse(int_nm %in% names(b), unname(b[int_nm]), 0)
    list(cerana = base_eff, mellifera = base_eff + int_eff)
  }
  me_list <- lapply(lc_used, me_fun)
  names(me_list) <- lc_used
  
  list(
    fit = fit,
    data = dat,
    coords_epsg = 5181,
    lon_col = guess_lon,
    lat_col = guess_lat,
    vars_used = list(landcover = lc_used, climate = clim_used),
    marginal_effects = me_list
  )
}

# ===================== 실행 파이프라인 =========================
# 1) 토지피복 읽기(복구 포함)
landcover_tbl <- load_landcover_and_map("data/raw/환경부_환경공간정보_년도별 대분류토지피복통계 현황_20230901.csv")
message("[CHECK] landcover 테이블 미리보기"); print(head(landcover_tbl))

# 2) 시군구 경계
sig_sf <- load_sigungu_sf("data/raw/BND_SIGUNGU_PG.shp")
message("[DEBUG] sig_sf CRS: "); print(st_crs(sig_sf))

# 3) bee + landcover 결합
both <- build_both_with_landcover(
  bee_c_path = "data/final/6_Apis_cerana_with_nearby.csv",
  bee_m_path = "data/final/6_Apis_mellifera_with_nearby.csv",
  landcover_tbl = landcover_tbl,
  sigungu_sf = sig_sf
)

message("[CHECK] both 미리보기"); print(head(both))
# both가 만들어진 뒤에 바로 실행
diagnose_both <- function(both) {
  nm <- names(both)
  lon <- nm[grepl("^(Longitude|lon|long|x|경도)$", nm, ignore.case=TRUE)][1]
  lat <- nm[grepl("^(Latitude|lat|y|위도)$",     nm, ignore.case=TRUE)][1]
  if (is.na(lon) || is.na(lat)) stop("좌표 열을 찾지 못했습니다. both의 열 이름을 확인하세요.")
  
  cat("# rows:", nrow(both), "\n")
  cat("# NAs by var:\n")
  na_cnt <- sapply(c("발생두수","species","forest","agri","urban","water", lon, lat),
                   function(x) if (x %in% nm) sum(is.na(both[[x]])) else NA_integer_)
  print(na_cnt)
  
  # 좌표 수치형/범위 점검
  x <- suppressWarnings(as.numeric(both[[lon]]))
  y <- suppressWarnings(as.numeric(both[[lat]]))
  cat("# coord summary:\n"); print(summary(x)); print(summary(y))
  out_kr <- which(x<124 | x>132 | y<33 | y>39) # 한반도 대략 범위
  cat("# out of Korea bbox:", length(out_kr), "\n")
}
diagnose_both(both)
# 4) 모델
res <- fit_both_car_with_lc(both)
summary(res$fit)

# 5) 토지피복 한계효과(기본종 cerana 기준, mellifera는 상호작용 더함)
b <- coef(res$fit)
me <- function(base, inter_name) {
  c_base <- unname(b[base]); c_int <- ifelse(paste0(base,":speciesmellifera") %in% names(b),
                                             unname(b[paste0(base,":speciesmellifera")]), 0)
  list(cerana = c_base, mellifera = c_base + c_int)
}
effects <- list(
  forest = me("forest"),
  agri   = me("agri"),
  urban  = me("urban"),
  water  = me("water")
)
message("[Marginal Effects]"); print(effects)