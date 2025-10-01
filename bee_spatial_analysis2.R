setwd("C:/.Soyeon/BeeDisease")
setwd("C:/Users/haa/Downloads/BeeDisease/BeeDisease")

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
})

# ---------------------------
# 0) 공통 유틸
# ---------------------------
require_pkg <- function(pkgs) {
  miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(miss)) stop("패키지 설치 필요: ", paste(miss, collapse=", "))
}
require_pkg(c("dplyr","readr","stringr","janitor","sf","spdep","spatialreg","ggplot2"))

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

# ---------------------------
# 2) 환경부 토지피복 로더(한글 컬럼 기반, 최신 연도, 비율 변환)
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

attach_landcover_to_points <- function(bee_path, sigungu_sf, landcover_tbl,
                                       lon_col="Longitude", lat_col="Latitude") {
  bee <- read.csv(bee_path)
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
# 5) 이웃/가중치 및 CAR 적합
# ---------------------------
build_neighbors <- function(df, lon_col="Longitude", lat_col="Latitude",
                            crs_proj=5181, d2_m=6000, style="B") {
  sf_pts  <- st_as_sf(df, coords = c(lon_col, lat_col), crs = 4326)
  sf_proj <- suppressWarnings(st_transform(sf_pts, crs_proj))
  coords  <- st_coordinates(sf_proj)
  nb      <- dnearneigh(coords, d1 = 0, d2 = d2_m)
  iso     <- which(card(nb) == 0)
  if (length(iso)) message(sprintf("[WARN] 고립점 %d개(이웃 0). zero.policy=TRUE로 처리.", length(iso)))
  listw   <- nb2listw(nb, style = style, zero.policy = TRUE)
  list(nb=nb, listw=listw, sf_proj=sf_proj)
}

# 🔧 수정: drop_na 후 데이터(dat)로 내부에서 이웃행렬 재생성
fit_car_both <- function(df, d2_m = 6000) {
  cand_lc   <- intersect(c("forest","agri","urban","water"), names(df))
  cand_clim <- intersect(c("prefer_count","mean_temp","sum_sun","ratio_sun","sum_rain"), names(df))
  stopifnot("발생두수" %in% names(df), "species" %in% names(df))
  if (length(cand_lc) == 0) stop("토지피복 변수가 없습니다.")
  
  fml_txt <- paste0("발생두수 ~ (", paste(c(cand_lc, cand_clim), collapse = " + "), ")*species")
  message("[MODEL] 공식: ", fml_txt)
  fml_obj <- stats::as.formula(fml_txt)
  
  need <- unique(c("발생두수","species","Longitude","Latitude", cand_lc, cand_clim))
  dat <- df |>
    dplyr::select(dplyr::any_of(need)) |>
    tidyr::drop_na(dplyr::any_of(c("발생두수","species","Longitude","Latitude", cand_lc))) |>
    dplyr::mutate(species = factor(species, levels = c("cerana","mellifera")))
  if (nrow(dat) < 20) stop(sprintf("유효 행이 너무 적습니다 (n=%d).", nrow(dat)))
  
  # dat로 이웃행렬 생성 (차원/순서 일치)
  g <- build_neighbors(dat, d2_m = d2_m, style = "B")
  listw_sym <- g$listw
  
  # formula 객체를 직접 전달
  fit <- spatialreg::spautolm(
    formula = fml_obj,
    data    = dat,
    listw   = listw_sym,
    family  = "CAR",
    zero.policy = TRUE
  )
  # 요 줄로 summary() 재평가 문제도 차단
  fit$call$formula <- fml_obj
  
  # 한계효과(landcover만)
  b <- coef(fit)
  me_fun <- function(base) {
    base_eff <- unname(b[base])
    int_nm   <- paste0(base, ":speciesmellifera")
    int_eff  <- ifelse(int_nm %in% names(b), unname(b[int_nm]), 0)
    list(cerana = base_eff, mellifera = base_eff + int_eff)
  }
  me_list <- lapply(cand_lc, me_fun); names(me_list) <- cand_lc
  
  list(fit = fit, data = dat, marginal_effects = me_list)
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
  na_cnt <- sapply(c("발생두수","species","forest","agri","urban","water", lon, lat),
                   function(x) if (x %in% nm) sum(is.na(df[[x]])) else NA_integer_)
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
# 7-1) Landcover, Shapefile 로드
landcover_tbl <- load_landcover_kmoe(PATH_LANDCOVER, use_grass_in_forest=TRUE,
                                     use_wetland_in_water=TRUE, to_ratio=TRUE)
sig_sf <- load_sigungu_sf(PATH_SIG_SHP)

# 7-2) 벌 포인트 + 시군구 + landcover 조인 (이름 우선, 실패시 코드 매핑)
cerana_df    <- attach_landcover_to_points(PATH_BEE_C, sig_sf, landcover_tbl) %>%
  mutate(species = if ("species" %in% names(.)) species else "cerana")
mellifera_df <- attach_landcover_to_points(PATH_BEE_M, sig_sf, landcover_tbl) %>%
  mutate(species = if ("species" %in% names(.)) species else "mellifera")

# 7-3) 결합 + 진단
both <- bind_rows(cerana_df, mellifera_df) %>%
  mutate(species = factor(species, levels = c("cerana","mellifera")))
diagnose_both(both)

# 7-4) CAR 적합 + 요약 + 한계효과 (이웃행렬은 함수 내부에서 생성)
res <- fit_car_both(both, d2_m = 6000)
print(summary(res$fit))

message("[Marginal Effects]")
print(res$marginal_effects)

# 7-5) 예측값(Xβ) 시각화
tmp <- make_X_from_model(res$fit, res$data)
res$data$yhat <- drop(tmp$X %*% tmp$b)

ggplot(res$data, aes(x = species, y = yhat)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  labs(title = "CAR(Xβ) 예측 분포 (종별)", x = "종", y = "예측 발생두수") +
  theme_minimal()

message("== 해석 체크 ==")
message("1) summary(res$fit)에서 p-value/CI 확인")
message("2) landcover는 조성자료: 기준범주 하나 드롭하거나 ILR 변환 고려")
message("3) 스케일: 비율(0-1). 10%p 변화 효과 = 계수×0.1")
message("4) count 자료: log1p(발생두수) 대안 또는 GLMM/INLA도 검토")


#######################
# 변수선택
fit_car_both <- function(df, d2_m = 6000, keep_lc = NULL, keep_clim = NULL) {
  # 후보군
  cand_lc_all   <- intersect(c("forest","agri","urban","water"), names(df))
  cand_clim_all <- intersect(c("prefer_count","mean_temp","sum_sun","ratio_sun","sum_rain"), names(df))
  
  # 사용 변수(선택적 필터링)
  cand_lc   <- cand_lc_all
  cand_clim <- cand_clim_all
  if (!is.null(keep_lc))   cand_lc   <- intersect(cand_lc_all, keep_lc)
  if (!is.null(keep_clim)) cand_clim <- intersect(cand_clim_all, keep_clim)
  if (length(cand_lc) == 0) stop("토지피복 선택 결과가 비었습니다. keep_lc 인자를 확인하세요.")
  
  stopifnot("발생두수" %in% names(df), "species" %in% names(df))
  
  used_vars <- c(cand_lc, cand_clim)
  message("[MODEL] 사용 landcover: ", paste(cand_lc, collapse=", "))
  if (length(cand_clim)) message("[MODEL] 사용 climate: ", paste(cand_clim, collapse=", "))
  message("[MODEL] 공식: 발생두수 ~ (", paste(used_vars, collapse = " + "), ")*species")
  
  fml_obj <- stats::as.formula(paste0("발생두수 ~ (", paste(used_vars, collapse = " + "), ")*species"))
  
  need <- unique(c("발생두수","species","Longitude","Latitude", used_vars))
  dat <- df |>
    dplyr::select(dplyr::any_of(need)) |>
    tidyr::drop_na(dplyr::any_of(c("발생두수","species","Longitude","Latitude", cand_lc))) |>
    dplyr::mutate(species = factor(species, levels = c("cerana","mellifera")))
  if (nrow(dat) < 20) stop(sprintf("유효 행이 너무 적습니다 (n=%d).", nrow(dat)))
  
  # dat로 이웃행렬 생성
  g <- build_neighbors(dat, d2_m = d2_m, style = "B")
  listw_sym <- g$listw
  
  fit <- spatialreg::spautolm(
    formula = fml_obj,
    data    = dat,
    listw   = listw_sym,
    family  = "CAR",
    zero.policy = TRUE
  )
  fit$call$formula <- fml_obj
  
  # 한계효과(선택된 landcover만)
  b <- coef(fit)
  me_fun <- function(base) {
    base_eff <- unname(b[base])
    int_nm   <- paste0(base, ":speciesmellifera")
    int_eff  <- ifelse(int_nm %in% names(b), unname(b[int_nm]), 0)
    list(cerana = base_eff, mellifera = base_eff + int_eff)
  }
  me_list <- lapply(cand_lc, me_fun); names(me_list) <- cand_lc
  
  list(fit = fit, data = dat, marginal_effects = me_list)
}

# 7-4) CAR 적합: forest + water 만 (기후변수 전부 제외)
res <- fit_car_both(
  both,
  d2_m = 6000,
  keep_lc   = c("forest","water", "agri"),
  keep_clim = character(0)   # ← 기후변수 강제 제외
)

print(summary(res$fit))

message("[Marginal Effects: forest & water only]")
print(res$marginal_effects)

# 예측값(Xβ) 박스플롯
tmp <- make_X_from_model(res$fit, res$data)
res$data$yhat <- drop(tmp$X %*% tmp$b)

ggplot(res$data, aes(x = species, y = yhat)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  labs(title = "CAR(Xβ) 예측 분포 (종별) — forest & water only",
       x = "종", y = "예측 발생두수") +
  theme_minimal()