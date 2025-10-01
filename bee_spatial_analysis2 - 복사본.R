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
  library(tidyr)
  if (!requireNamespace("car", quietly = TRUE)) install.packages("car")
  library(car)
})

# ---------------------------
# 0) 설정값(편집 가능)
# ---------------------------
PATH_LANDCOVER <- "data/raw/환경부_환경공간정보_년도별 대분류토지피복통계 현황_20230901.csv"
PATH_SIG_SHP   <- "data/raw/BND_SIGUNGU_PG.shp"
PATH_BEE_C     <- "data/final/6_Apis_cerana_with_nearby.csv"
PATH_BEE_M     <- "data/final/6_Apis_mellifera_with_nearby.csv"

# 변수 선택 옵션
VIF_THRESHOLD  <- 5
DROP_BASELINE  <- TRUE
BASELINE_VAR   <- "urban" # 드롭할 기준범주
STANDARDIZE_X  <- TRUE    # landcover 연속변수 z-score
NEIGHBOR_R_M   <- 8000    # 이웃 반경

# ---------------------------
# 1) 공통 유틸/로더
# ---------------------------
require_pkg <- function(pkgs) {
  miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(miss)) stop("패키지 설치 필요: ", paste(miss, collapse=", "))
}
require_pkg(c("dplyr","readr","stringr","janitor","sf","spdep","spatialreg","ggplot2","car","tidyr"))

guess_delim_from_lines <- function(x, n = 20) {
  lines <- head(x, n)
  cand  <- c("," = ",", "\t" = "\t", ";" = ";", "|" = "|")
  counts <- sapply(cand, function(d) sum(str_count(lines, fixed(d))))
  names(which.max(counts))[1]
}
load_table_resilient <- function(path) {
  if (!file.exists(path)) stop(sprintf("파일 없음: %s", path))
  message(sprintf("[INFO] 파일 크기: %s bytes", file.size(path)))
  raw <- read_file_raw(path)
  enc_guess <- try(guess_encoding(raw[1:min(length(raw), 200000)]), silent = TRUE)
  enc <- "CP949"; if (!inherits(enc_guess, "try-error") && nrow(enc_guess) > 0) enc <- enc_guess$encoding[1]
  message(sprintf("[INFO] 인코딩 추정: %s", enc))
  lines <- read_lines(path, locale = locale(encoding = enc), progress = FALSE)
  if (length(lines) > 0 && !str_detect(lines[length(lines)], "\r$|\n$")) lines[length(lines)] <- paste0(lines[length(lines)], "\n")
  delim <- guess_delim_from_lines(lines, n = 20); message(sprintf("[INFO] 구분자 추정: '%s'", ifelse(delim == "\t", "\\t", delim)))
  tf <- tempfile(fileext = ".csv"); write_lines(lines, tf)
  df <- try(read_delim(tf, delim = delim, locale = locale(encoding = enc), na = c("", "NA", "-", "NaN"),
                       guess_max = 100000, show_col_types = FALSE), silent = TRUE)
  if (!inherits(df, "try-error") && nrow(df) > 0 && ncol(df) > 1) return(df)
  if (requireNamespace("data.table", quietly = TRUE)) {
    message("[INFO] readr 실패 → data.table::fread 시도")
    df2 <- try(data.table::fread(tf, encoding = enc, sep = delim, na.strings = c("", "NA", "-", "NaN")), silent = TRUE)
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
norm_key <- function(x) {
  x <- trimws(as.character(x)); x <- gsub("\u3000", "", x); x <- gsub("\\s+", "", x); x <- gsub("\\(.*\\)$", "", x); x
}

load_landcover_kmoe <- function(path_landcover, use_grass_in_forest=TRUE, use_wetland_in_water=TRUE, to_ratio=TRUE) {
  lc_raw <- load_table_resilient(path_landcover)
  req <- c("자료년도","국가구분","시도","시군구","시가화건조지역","농업지역","산림지역","초지","습지","수역","합계")
  miss <- setdiff(req, names(lc_raw)); if (length(miss)) stop("landcover 필수 컬럼 누락: ", paste(miss, collapse=", "))
  lc <- lc_raw %>% filter(grepl("대한민국|남한", 국가구분)) %>% filter(자료년도 == max(자료년도, na.rm=TRUE))
  forest_base <- lc[["산림지역"]]; if (use_grass_in_forest && "초지" %in% names(lc)) forest_base <- forest_base + lc[["초지"]]
  water_base  <- lc[["수역"]];     if (use_wetland_in_water && "습지" %in% names(lc)) water_base  <- water_base  + lc[["습지"]]
  if (to_ratio) {
    total <- as.numeric(lc[["합계"]])
    forest_v <- as.numeric(forest_base)/total; agri_v <- as.numeric(lc[["농업지역"]])/total
    urban_v  <- as.numeric(lc[["시가화건조지역"]])/total; water_v <- as.numeric(water_base)/total
  } else {
    forest_v <- as.numeric(forest_base); agri_v <- as.numeric(lc[["농업지역"]]); urban_v <- as.numeric(lc[["시가화건조지역"]]); water_v <- as.numeric(water_base)
  }
  out <- lc %>% transmute(SIGUNGU_NM = norm_key(시군구), forest=forest_v, agri=agri_v, urban=urban_v, water=water_v) %>%
    group_by(SIGUNGU_NM) %>% summarise(across(c(forest,agri,urban,water), ~mean(., na.rm=TRUE)), .groups="drop")
  message("[OK] Landcover 준비 완료: 연도=", max(lc$자료년도, na.rm=TRUE), ", 스케일=", ifelse(to_ratio,"비율(0-1)","원단위"))
  out
}
load_sigungu_sf <- function(path_shp = PATH_SIG_SHP) {
  sf <- st_read(path_shp, options = "ENCODING=CP949", quiet = TRUE)
  if (!all(st_is_valid(sf))) sf <- st_make_valid(sf)
  if (is.na(st_crs(sf))) stop("시군구 shp에 CRS 없음. st_set_crs로 EPSG 지정 필요.")
  sf
}
detect_sigungu_keys <- function(sig_sf) {
  nms <- names(sig_sf)
  name_pat <- "(SIG(_KOR)?_?NM|SIGUNGU_?NM|SGG_?NM|EMD(C)?_?NM|SIG_NM|SIGNM|시군구.?명|시군구.?이름)"
  code_pat <- "(SIG(_?CD)?|SIGUNGU_?CD|SGG_?CD|ADM_?CD|LAWD_?CD|CODE|코드)"
  name_col <- nms[grepl(name_pat, nms, ignore.case = TRUE)]
  code_col <- nms[grepl(code_pat, nms, ignore.case = TRUE)]
  if (length(name_col) == 0 && "SIG_KOR_NM" %in% nms) name_col <- "SIG_KOR_NM"
  if (length(name_col) == 0 && "SIG_NM"     %in% nms) name_col <- "SIG_NM"
  if (length(code_col) == 0 && "SIG_CD"     %in% nms) code_col <- "SIG_CD"
  list(name_col = if (length(name_col)) name_col[1] else NA_character_,
       code_col = if (length(code_col)) code_col[1] else NA_character_)
}
attach_landcover_to_points <- function(bee_path, sigungu_sf, landcover_tbl, lon_col="Longitude", lat_col="Latitude") {
  bee <- read.csv(bee_path); stopifnot(all(c(lon_col, lat_col) %in% names(bee)))
  keys <- detect_sigungu_keys(sigungu_sf); nm_col <- keys$name_col; cd_col <- keys$code_col
  if (is.na(nm_col) && is.na(cd_col)) stop("shp에서 시군구 이름/코드 컬럼을 찾지 못했습니다.")
  pts <- st_as_sf(bee, coords = c(lon_col, lat_col), crs = 4326, remove = FALSE) |> st_transform(st_crs(sigungu_sf))
  joined <- st_join(pts, suppressWarnings(st_buffer(sigungu_sf, 0))[, c(na.omit(c(nm_col, cd_col))), drop = FALSE],
                    join = st_intersects, left = TRUE) |> st_drop_geometry()
  if (!is.na(nm_col) && nm_col %in% names(joined)) joined[[nm_col]] <- norm_key(joined[[nm_col]])
  landcover_tbl$SIGUNGU_NM <- norm_key(landcover_tbl$SIGUNGU_NM)
  # 1) 이름 조인
  if (!is.na(nm_col) && nm_col %in% names(joined)) {
    df1 <- left_join(joined, landcover_tbl, by = setNames("SIGUNGU_NM", nm_col))
    na1 <- rowSums(is.na(df1[, c("forest","agri","urban","water")]))
    if (sum(na1 > 0) < nrow(df1)) { message(sprintf("[INFO] 이름기반 조인: landcover NA %d / %d", sum(na1 > 0), nrow(df1))); return(df1) }
    message("[WARN] 이름기반 조인 전부 NA — 코드매핑 시도")
  }
  # 2) 코드→이름 매핑
  if (!is.na(cd_col) && !is.na(nm_col) && cd_col %in% names(sigungu_sf) && nm_col %in% names(sigungu_sf)) {
    xwalk <- sigungu_sf |> st_drop_geometry() |> select(!!nm_col, !!cd_col) |>
      mutate(!!nm_col := norm_key(.data[[nm_col]]), !!cd_col := as.character(.data[[cd_col]])) |> distinct()
    if (cd_col %in% names(joined)) {
      joined$.__sgg_code__ <- as.character(joined[[cd_col]])
      df2 <- left_join(joined, xwalk, by = setNames(cd_col, ".__sgg_code__"))
      df2$.__sgg_name__ <- if (nm_col %in% names(df2)) df2[[nm_col]] else NA_character_
      if (nm_col %in% names(df2)) df2$.__sgg_name__ <- ifelse(is.na(df2$.__sgg_name__), df2[[nm_col]], df2$.__sgg_name__)
      df2 <- left_join(df2, rename(landcover_tbl, .__sgg_name__ = SIGUNGU_NM), by = ".__sgg_name__")
      na2 <- rowSums(is.na(df2[, c("forest","agri","urban","water")]))
      message(sprintf("[INFO] 코드↔이름 매핑 후 조인: landcover NA %d / %d", sum(na2 > 0), nrow(df2)))
      return(df2)
    }
  }
  message("[FAIL] landcover 조인 실패. shp 컬럼/값 확인 요망."); joined
}

# ---------------------------
# 2) 진단/이웃/모형 헬퍼
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
build_neighbors <- function(df, lon_col="Longitude", lat_col="Latitude", crs_proj=5181, d2_m=6000, style="B") {
  sf_pts  <- st_as_sf(df, coords = c(lon_col, lat_col), crs = 4326)
  sf_proj <- suppressWarnings(st_transform(sf_pts, crs_proj))
  coords  <- st_coordinates(sf_proj)
  nb      <- dnearneigh(coords, d1 = 0, d2 = d2_m)
  iso     <- which(card(nb) == 0)
  if (length(iso)) message(sprintf("[WARN] 고립점 %d개(이웃 0). zero.policy=TRUE로 처리.", length(iso)))
  listw   <- nb2listw(nb, style = style, zero.policy = TRUE)
  list(nb=nb, listw=listw, sf_proj=sf_proj)
}

# ---------------------------
# 3) 변수 선택(비공간 lm 기반) + CAR 적합
# ---------------------------
drop_baseline_if_needed <- function(df, drop = if (DROP_BASELINE) BASELINE_VAR else NULL) {
  vars <- c("forest","agri","urban","water")
  if (is.null(drop) || !(drop %in% vars)) return(df)
  keep <- setdiff(vars, drop)
  message(sprintf("[DROP] 기준범주로 '%s' 드롭 → 사용: %s", drop, paste(keep, collapse=", ")))
  df[, unique(c("발생두수","species","Longitude","Latitude", keep)), drop=FALSE]
}
compute_vif_lm <- function(df, rhs_vars) {
  ftxt <- paste0("발생두수 ~ ", paste(rhs_vars, collapse = " + "))
  fit  <- lm(as.formula(ftxt), data = df)
  as.numeric(car::vif(fit)) |> setNames(names(car::vif(fit)))
}
vif_step_select <- function(df, base_vars = c("species"), cand_vars, vif_thr = VIF_THRESHOLD) {
  tmp <- df[, unique(c("발생두수","Longitude","Latitude", base_vars, cand_vars)), drop=FALSE]
  tmp <- tidyr::drop_na(tmp, any_of(c("발생두수", base_vars, cand_vars)))
  keep <- cand_vars
  log_rows <- list()
  repeat {
    rhs <- c(base_vars, keep)
    v   <- try(compute_vif_lm(tmp, rhs), silent = TRUE)
    if (inherits(v, "try-error")) break
    v_lc <- v[names(v) %in% keep]
    log_rows[[length(log_rows)+1]] <- data.frame(var=names(v_lc), VIF=as.numeric(v_lc))
    if (length(v_lc) == 0) break
    v_max <- max(v_lc, na.rm = TRUE)
    if (is.finite(v_max) && v_max > vif_thr) {
      worst <- names(which.max(v_lc))
      message(sprintf("[VIF] %s=%.2f > %.1f → 제거", worst, v_max, vif_thr))
      keep <- setdiff(keep, worst)
      if (length(keep) == 0) break
    } else break
  }
  log_df <- if (length(log_rows)) bind_rows(log_rows, .id = "iter") else NULL
  list(selected = keep, log = log_df)
}
fit_car_with_vars <- function(df, vars, d2_m = NEIGHBOR_R_M, standardize = STANDARDIZE_X) {
  stopifnot(length(vars) >= 1)
  # === (중요) 행 식별자 보존: 이후 상위10% 서브셋 매칭용 ===
  dat <- df %>%
    mutate(id_row = row_number()) %>%                              # <-- 추가
    select(any_of(c("id_row","발생두수","species","Longitude","Latitude", vars))) %>%
    drop_na(any_of(c("발생두수","species","Longitude","Latitude", vars))) %>%
    mutate(species = factor(species, levels = c("cerana","mellifera")))
  if (nrow(dat) < 20) stop(sprintf("유효 행이 너무 적습니다 (n=%d).", nrow(dat)))
  if (standardize) {
    num_vars <- vars[vars %in% c("forest","agri","urban","water")]
    if (length(num_vars)) dat[, num_vars] <- lapply(dat[, num_vars, drop=FALSE], scale)
    message("[STD] 연속 landcover 표준화 적용: ", paste(num_vars, collapse=", "))
  }
  g <- build_neighbors(dat, d2_m = d2_m, style = "B")
  listw_sym <- g$listw
  fml_txt <- paste0("발생두수 ~ (", paste(vars, collapse = " + "), ")*species")
  fml_obj <- stats::as.formula(fml_txt)
  fit <- spatialreg::spautolm(formula = fml_obj, data = dat, listw = listw_sym, family = "CAR", zero.policy = TRUE)
  fit$call$formula <- fml_obj
  b <- coef(fit)
  lc_used <- intersect(vars, c("forest","agri","urban","water"))
  me_fun <- function(base) {
    base_eff <- unname(b[base]); int_nm <- paste0(base, ":speciesmellifera")
    int_eff  <- ifelse(int_nm %in% names(b), unname(b[int_nm]), 0)
    list(cerana = base_eff, mellifera = base_eff + int_eff)
  }
  me_list <- if (length(lc_used)) { out <- lapply(lc_used, me_fun); names(out) <- lc_used; out } else list()
  list(fit=fit, data=dat, vars=vars, marginal_effects=me_list)
}

# ---------------------------
# 4) 실행 파이프라인
# ---------------------------
# 4-1) 데이터 로드
landcover_tbl <- load_landcover_kmoe(PATH_LANDCOVER, use_grass_in_forest=TRUE, use_wetland_in_water=TRUE, to_ratio=TRUE)
sig_sf <- load_sigungu_sf(PATH_SIG_SHP)
cerana_df    <- attach_landcover_to_points(PATH_BEE_C, sig_sf, landcover_tbl) %>% mutate(species = if ("species" %in% names(.)) species else "cerana")
mellifera_df <- attach_landcover_to_points(PATH_BEE_M, sig_sf, landcover_tbl) %>% mutate(species = if ("species" %in% names(.)) species else "mellifera")
both <- bind_rows(cerana_df, mellifera_df) %>% mutate(species = factor(species, levels = c("cerana","mellifera")))
diagnose_both(both)

# 4-2) 변수 준비 (기준범주 드롭 → VIF 선택)
base_df <- drop_baseline_if_needed(both, drop = if (DROP_BASELINE) BASELINE_VAR else NULL)
# (중요) 이후 상위 10% 컷은 base_df의 원스케일 agri로 계산한다.
cand_vars <- intersect(c("forest","agri","urban","water"), names(base_df))
message("[INFO] 후보 landcover: ", paste(cand_vars, collapse=", "))

sel <- vif_step_select(base_df, base_vars = "species", cand_vars = cand_vars, vif_thr = VIF_THRESHOLD)
message("[INFO] VIF 선택 결과: ", paste(sel$selected, collapse=", "))
if (!is.null(sel$log)) {
  message("[VIF log] 마지막 반복의 VIF:")
  print(sel$log %>% group_by(var) %>% summarize(lastVIF = dplyr::last(VIF)))
}

# 4-3) 선택 변수로 CAR 적합
res <- fit_car_with_vars(base_df, vars = sel$selected, d2_m = NEIGHBOR_R_M, standardize = STANDARDIZE_X)
print(summary(res$fit))

message("[Marginal Effects]")
print(res$marginal_effects)
##########
# 4-2) 변수 준비 (기준범주만 드롭하고 VIF 선택은 생략)
base_df <- drop_baseline_if_needed(both, drop = if (DROP_BASELINE) BASELINE_VAR else NULL)

# 전체 모형: landcover 전부 사용(단, BASELINE_VAR은 합 1 제약 때문에 제외)
vars_all <- intersect(c("forest","agri","urban","water"), names(base_df))
if (DROP_BASELINE && BASELINE_VAR %in% vars_all) {
  vars_all <- setdiff(vars_all, BASELINE_VAR)
  message(sprintf("[FULL] 기준범주 '%s'만 제외하고 전부 사용 → %s",
                  BASELINE_VAR, paste(vars_all, collapse=", ")))
} else {
  message(sprintf("[FULL] 모든 landcover 사용 → %s", paste(vars_all, collapse=", ")))
}

# 4-3) 전체 변수로 CAR 적합 (VIF 미사용)
res <- fit_car_with_vars(base_df, vars = vars_all, d2_m = NEIGHBOR_R_M, standardize = STANDARDIZE_X)
print(summary(res$fit))

message("[Marginal Effects]")
print(res$marginal_effects)

# 4-4) 예측값(Xβ) 시각화(전체)
make_X_from_model <- function(fit, newdata) {
  b  <- coef(fit); ff <- formula(fit); X <- model.matrix(ff, data = newdata)
  if ("(Intercept)" %in% names(b) && !("(Intercept)" %in% colnames(X))) X <- cbind("(Intercept)"=1, X)
  miss <- setdiff(names(b), colnames(X))
  if (length(miss)) { Z <- matrix(0, nrow=nrow(X), ncol=length(miss)); colnames(Z) <- miss; X <- cbind(X, Z) }
  X <- X[, names(b), drop=FALSE]; list(X=X, b=as.numeric(b))
}
tmp_all <- make_X_from_model(res$fit, res$data)
res$data$yhat <- drop(tmp_all$X %*% tmp_all$b)

ggplot(res$data, aes(x = species, y = yhat)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  labs(title = sprintf("CAR(Xβ) 예측 분포 (선택 변수: %s)", paste(res$vars, collapse=", ")),
       x = "종", y = "예측 발생두수") +
  theme_minimal()

message("== 해석 가이드 ==")
message("* p-value 개선: 기준범주 드롭, VIF 축소, 표준화, 반경 튜닝(", NEIGHBOR_R_M, "m)")
message("* 추가: 기상/밀원 포함, kNN 이웃, Poisson/NegBin 공간모형(INLA/BRMS) 고려")


# ======================================================================
# [수역(water) 상위 집단 분석]
# ======================================================================

stopifnot(exists("base_df"), exists("res"))

# 1) 상위 10% 컷오프 계산 (원자료 water 기준)
q_top   <- 0.90
thr_water <- quantile(base_df$water, q_top, na.rm = TRUE)
message(sprintf("[INFO] water 상위 %.0f%% 컷오프 = %.5f", (1 - q_top) * 100, thr_water))

# 2) 컷오프 충족 행의 전역 id 식별
base_df <- base_df %>% mutate(id_row = row_number())
high_ids_w <- base_df %>%
  filter(!is.na(water) & water >= thr_water) %>%
  pull(id_row)

high_water <- res$data %>%
  filter(id_row %in% high_ids_w) %>%
  mutate(species = factor(species, levels = c("cerana","mellifera")))

message(sprintf("[INFO] 상위 집단 크기(water): %d / %d (%.1f%%)",
                nrow(high_water), nrow(res$data), 100 * nrow(high_water) / nrow(res$data)))
print(table(high_water$species))

# 3) 관측값 기준 종별 요약 + 탐색적 검정
obs_stats_w <- high_water %>%
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
print(obs_stats_w)

if (length(unique(high_water$species)) == 2) {
  message("[탐색적] t-test (water 상위):")
  print(try(t.test(발생두수 ~ species, data = high_water), silent = TRUE))
  message("[탐색적] Wilcoxon (water 상위):")
  print(try(wilcox.test(발생두수 ~ species, data = high_water, exact = FALSE), silent = TRUE))
}

ggplot(high_water, aes(x = species, y = 발생두수)) +
  geom_boxplot() +
  labs(title = sprintf("관측값 분포 (water 상위 %.0f%%, 임계=%.5f)", (1 - q_top) * 100, thr_water),
       x = "종", y = "발생두수")

# 4) CAR 모형 기반 예측값 (E[y|X] = Xβ)
tmp_sub_w <- make_X_from_model(res$fit, high_water)
high_water$yhat <- drop(tmp_sub_w$X %*% tmp_sub_w$b)

pred_stats_w <- high_water %>%
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
print(pred_stats_w)

ggplot(high_water, aes(x = species, y = yhat)) +
  geom_boxplot() +
  labs(title = sprintf("모형 예측값 분포 (water 상위 %.0f%%, 임계=%.5f)", (1 - q_top) * 100, thr_water),
       x = "종", y = "예측 발생두수 (CAR, Xβ)")

# 5) Counterfactual: 같은 관측에서 species만 바꾸기
lvl <- levels(res$data$species)
high_water_cer <- high_water; high_water_cer$species <- factor("cerana",    levels = lvl)
high_water_mel <- high_water; high_water_mel$species <- factor("mellifera", levels = lvl)

tmp_cer_w <- make_X_from_model(res$fit, high_water_cer)
tmp_mel_w <- make_X_from_model(res$fit, high_water_mel)

high_water$yhat_cer_cf <- drop(tmp_cer_w$X %*% tmp_cer_w$b)
high_water$yhat_mel_cf <- drop(tmp_mel_w$X %*% tmp_mel_w$b)
high_water$delta_species <- high_water$yhat_mel_cf - high_water$yhat_cer_cf

delta_summary_w <- high_water %>%
  summarise(
    n         = n(),
    mean_diff = mean(delta_species, na.rm = TRUE),
    med_diff  = median(delta_species, na.rm = TRUE),
    q25_diff  = quantile(delta_species, 0.25, na.rm = TRUE),
    q75_diff  = quantile(delta_species, 0.75, na.rm = TRUE)
  )
print(delta_summary_w)

ggplot(high_water, aes(x = "", y = delta_species)) +
  geom_boxplot() +
  coord_flip() +
  labs(title = "Counterfactual: 같은 셋에서 species만 변경 (mellifera - cerana) — water 상위",
       x = NULL, y = "예측 차이")

# 6) water 한계효과(계수) 확인
b <- coef(res$fit)
if ("water" %in% names(b)) {
  me_cerana_w    <- unname(b["water"])
  me_mellifera_w <- unname(b["water"] + ifelse("water:speciesmellifera" %in% names(b), b["water:speciesmellifera"], 0))
  message(sprintf("[water 한계효과] cerana=%.3f, mellifera=%.3f", me_cerana_w, me_mellifera_w))
} else {
  message("[주의] 최종 선택 변수에 'water'가 포함되지 않아 한계효과를 계산할 수 없습니다.")
}

# 7) (선택) delta 분포 시각화
ggplot(high_water, aes(x = delta_species)) +
  geom_histogram(binwidth = 5, color = "black", fill = "skyblue") +
  labs(title = "Counterfactual 차이 분포 (mellifera - cerana) — water 상위",
       x = "예측 차이", y = "빈도")

ggplot(high_water, aes(x = delta_species)) +
  geom_density(fill = "orange", alpha = 0.4) +
  labs(title = "Counterfactual 차이 분포 (mellifera - cerana) — water 상위",
       x = "예측 차이", y = "밀도")

# 8) (선택) 결과 저장 예시
# dir.create("out_water", showWarnings = FALSE)
# write.csv(obs_stats_w,  "out_water/highwater_obs_stats.csv",  row.names = FALSE)
# write.csv(pred_stats_w, "out_water/highwater_pred_stats.csv", row.names = FALSE)
# write.csv(high_water[, c("species","발생두수","yhat","yhat_cer_cf","yhat_mel_cf","delta_species")],
#           "out_water/highwater_rowwise_preds.csv", row.names = FALSE)


# ======================================================================
# [산림(forest) 상위 집단 분석]
# ======================================================================

stopifnot(exists("base_df"), exists("res"))

# 1) 상위 10% 컷오프 계산 (원자료 forest 기준)
q_top   <- 0.90
thr_forest <- quantile(base_df$forest, q_top, na.rm = TRUE)
message(sprintf("[INFO] forest 상위 %.0f%% 컷오프 = %.5f", (1 - q_top) * 100, thr_forest))

# 2) 컷오프 충족 행의 전역 id 식별
base_df <- base_df %>% mutate(id_row = row_number())
high_ids_f <- base_df %>%
  filter(!is.na(forest) & forest >= thr_forest) %>%
  pull(id_row)

high_forest <- res$data %>%
  filter(id_row %in% high_ids_f) %>%
  mutate(species = factor(species, levels = c("cerana","mellifera")))

message(sprintf("[INFO] 상위 집단 크기(forest): %d / %d (%.1f%%)",
                nrow(high_forest), nrow(res$data), 100 * nrow(high_forest) / nrow(res$data)))
print(table(high_forest$species))

# 3) 관측값 기준 종별 요약 + 탐색적 검정
obs_stats_f <- high_forest %>%
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
print(obs_stats_f)

if (length(unique(high_forest$species)) == 2) {
  message("[탐색적] t-test (forest 상위):")
  print(try(t.test(발생두수 ~ species, data = high_forest), silent = TRUE))
  message("[탐색적] Wilcoxon (forest 상위):")
  print(try(wilcox.test(발생두수 ~ species, data = high_forest, exact = FALSE), silent = TRUE))
}

ggplot(high_forest, aes(x = species, y = 발생두수)) +
  geom_boxplot() +
  labs(title = sprintf("관측값 분포 (forest 상위 %.0f%%, 임계=%.5f)", (1 - q_top) * 100, thr_forest),
       x = "종", y = "발생두수")

# 4) CAR 모형 기반 예측값 (E[y|X] = Xβ)
tmp_sub_f <- make_X_from_model(res$fit, high_forest)
high_forest$yhat <- drop(tmp_sub_f$X %*% tmp_sub_f$b)

pred_stats_f <- high_forest %>%
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
print(pred_stats_f)

ggplot(high_forest, aes(x = species, y = yhat)) +
  geom_boxplot() +
  labs(title = sprintf("모형 예측값 분포 (forest 상위 %.0f%%, 임계=%.5f)", (1 - q_top) * 100, thr_forest),
       x = "종", y = "예측 발생두수 (CAR, Xβ)")

# 5) Counterfactual: 같은 관측에서 species만 바꾸기
lvl <- levels(res$data$species)
high_forest_cer <- high_forest; high_forest_cer$species <- factor("cerana",    levels = lvl)
high_forest_mel <- high_forest; high_forest_mel$species <- factor("mellifera", levels = lvl)

tmp_cer_f <- make_X_from_model(res$fit, high_forest_cer)
tmp_mel_f <- make_X_from_model(res$fit, high_forest_mel)

high_forest$yhat_cer_cf <- drop(tmp_cer_f$X %*% tmp_cer_f$b)
high_forest$yhat_mel_cf <- drop(tmp_mel_f$X %*% tmp_mel_f$b)
high_forest$delta_species <- high_forest$yhat_mel_cf - high_forest$yhat_cer_cf

delta_summary_f <- high_forest %>%
  summarise(
    n         = n(),
    mean_diff = mean(delta_species, na.rm = TRUE),
    med_diff  = median(delta_species, na.rm = TRUE),
    q25_diff  = quantile(delta_species, 0.25, na.rm = TRUE),
    q75_diff  = quantile(delta_species, 0.75, na.rm = TRUE)
  )
print(delta_summary_f)

ggplot(high_forest, aes(x = "", y = delta_species)) +
  geom_boxplot() +
  coord_flip() +
  labs(title = "Counterfactual: 같은 셋에서 species만 변경 (mellifera - cerana) — forest 상위",
       x = NULL, y = "예측 차이")

# 6) forest 한계효과(계수) 확인
b <- coef(res$fit)
if ("forest" %in% names(b)) {
  me_cerana_f    <- unname(b["forest"])
  me_mellifera_f <- unname(b["forest"] + ifelse("forest:speciesmellifera" %in% names(b), b["forest:speciesmellifera"], 0))
  message(sprintf("[forest 한계효과] cerana=%.3f, mellifera=%.3f", me_cerana_f, me_mellifera_f))
} else {
  message("[주의] 최종 선택 변수에 'forest'가 포함되지 않아 한계효과를 계산할 수 없습니다.")
}

# 7) (선택) delta 분포 시각화
ggplot(high_forest, aes(x = delta_species)) +
  geom_histogram(binwidth = 5, color = "black", fill = "skyblue") +
  labs(title = "Counterfactual 차이 분포 (mellifera - cerana) — forest 상위",
       x = "예측 차이", y = "빈도")

ggplot(high_forest, aes(x = delta_species)) +
  geom_density(fill = "orange", alpha = 0.4) +
  labs(title = "Counterfactual 차이 분포 (mellifera - cerana) — forest 상위",
       x = "예측 차이", y = "밀도")

# 8) (선택) 결과 저장 예시
# dir.create("out_forest", showWarnings = FALSE)
# write.csv(obs_stats_f,  "out_forest/highforest_obs_stats.csv",  row.names = FALSE)
# write.csv(pred_stats_f, "out_forest/highforest_pred_stats.csv", row.names = FALSE)
# write.csv(high_forest[, c("species","발생두수","yhat","yhat_cer_cf","yhat_mel_cf","delta_species")],
#           "out_forest/highforest_rowwise_preds.csv", row.names = FALSE)