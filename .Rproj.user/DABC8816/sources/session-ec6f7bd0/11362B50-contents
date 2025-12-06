#' @useDynLib DRsurvCRT, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- expose Surv() / cluster() to formulas ---------------------------------
#' @keywords internal
#' @noRd
.surv_env <- function(parent = parent.frame()) {
  env <- new.env(parent = parent)
  env$Surv    <- survival::Surv
  env$cluster <- survival::cluster
  env
}

# ---- parse Surv(time, status) from LHS --------------------------------------
#' @keywords internal
#' @noRd
.surv_lhs <- function(formula) {
  lhs <- as.character(formula)[2]
  # allow survival::Surv(...)
  lhs <- sub("^[[:alnum:]_.:]*::(Surv\\()", "\\1", lhs)
  inside <- sub("^\\s*Surv\\s*\\((.*)\\)\\s*$", "\\1", lhs)
  parts <- trimws(strsplit(inside, ",", fixed = TRUE)[[1]])
  parts <- gsub("^`|`$", "", parts)
  if (length(parts) < 2L) stop("formula must be Surv(time, status) ~ ...", call. = FALSE)
  list(time = parts[1], status = parts[2])
}

# ---- find cluster() and strip from RHS --------------------------------------
#' @keywords internal
#' @noRd
.extract_cluster <- function(formula) {
  tl  <- attr(terms(formula), "term.labels")
  if (is.null(tl)) tl <- character(0)
  is_cl <- grepl("^cluster\\s*\\(", tl)
  cl_var <- NULL
  if (any(is_cl)) {
    cl_inside <- sub("^cluster\\s*\\((.*)\\)\\s*$", "\\1", tl[is_cl][1])
    cl_var <- trimws(cl_inside)
  }
  rhs_wo_cluster <- paste(tl[!is_cl], collapse = if (length(tl[!is_cl])) " + " else "")
  list(cluster = cl_var, rhs_wo_cluster = rhs_wo_cluster)
}

# ---- default censoring formula builder --------------------------------------
#' @keywords internal
#' @noRd
.build_cens_formula_from <- function(nm, rhs_wo_cluster, clvar = NULL) {
  base <- sprintf("Surv(%s,%s==0) ~ %s", nm$time, nm$status, rhs_wo_cluster)
  if (!is.null(clvar) && nzchar(clvar)) paste0(base, " + cluster(", clvar, ")") else base
}

# ---- model matrices WITHOUT cluster() ---------------------------------------
#' @keywords internal
#' @noRd
.model_mats_safe <- function(formula, cens_formula, data) {
  info <- .extract_cluster(formula)
  fX   <- if (nzchar(info$rhs_wo_cluster)) as.formula(paste("~", info$rhs_wo_cluster)) else ~ 1
  Xs   <- model.matrix(fX, data = data)
  if (ncol(Xs) && colnames(Xs)[1] == "(Intercept)") Xs <- Xs[, -1, drop = FALSE]

  infoC <- .extract_cluster(cens_formula)
  fC    <- if (nzchar(infoC$rhs_wo_cluster)) as.formula(paste("~", infoC$rhs_wo_cluster)) else ~ 1
  Xc    <- model.matrix(fC, data = data)
  if (ncol(Xc) && colnames(Xc)[1] == "(Intercept)") Xc <- Xc[, -1, drop = FALSE]

  list(Xs = Xs, Xc = Xc)
}

# ---- treatment probabilities (cluster-level) --------------------------------
#' @keywords internal
#' @noRd
.compute_trt_probs <- function(data, trt, strata_col) {
  cl_ids <- unique(data[[strata_col]])
  A_first <- vapply(cl_ids, function(g) data[[trt]][which(data[[strata_col]] == g)[1]], numeric(1))
  tab <- prop.table(table(A_first))
  if (length(tab) != 2L) stop("Treatment must be binary at the cluster level.", call. = FALSE)
  c(p0 = as.numeric(tab[1]), p1 = as.numeric(tab[2]))
}

# ---- prob validity -----------------------------------------------------------
#' @keywords internal
#' @noRd
.validate_probs <- function(trt_prob1, trt_prob0, tol = 1e-8) {
  if (!is.numeric(trt_prob1) || !is.numeric(trt_prob0) ||
      length(trt_prob1) != 1L || length(trt_prob0) != 1L) {
    stop("trt_prob1 and trt_prob0 must be numeric scalars.")
  }
  if (trt_prob1 <= 0 || trt_prob1 >= 1 || trt_prob0 <= 0 || trt_prob0 >= 1) {
    stop("trt_prob1 and trt_prob0 must be in (0,1).")
  }
  if (abs(trt_prob1 + trt_prob0 - 1) > tol) {
    stop("trt_prob1 + trt_prob0 must equal 1.")
  }
  invisible(TRUE)
}

# ---- survival-grid sanity warning -------------------------------------------
#' @keywords internal
#' @noRd
.check_S_warn <- function(S, name = "S", cols = c(1, 2)) {
  if (!is.matrix(S)) {
    warning(sprintf("%s should be a numeric matrix [S1, S0, diff].", name), call. = FALSE)
    return(invisible(FALSE))
  }
  cols <- cols[cols %in% seq_len(ncol(S))]
  if (!length(cols)) return(invisible(TRUE))
  zero_prop <- vapply(cols, function(j) mean(S[, j] == 0, na.rm = TRUE), numeric(1))
  if (any(zero_prop >= 0.5)) {
    bad <- cols[zero_prop >= 0.5]
    warning(sprintf("%s: column(s) %s have ≥ 50%% zeros; RMST/SPCE may be unstable.",
                    name, paste(bad, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}
