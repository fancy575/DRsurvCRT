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

# ---- survival-array checks (multi-state) ------------------------------------
#' @keywords internal
#' @noRd
.check_S_array <- function(S, name = "S", tol = 1e-9) {
  if (!is.array(S) || length(dim(S)) != 3L) {
    stop(sprintf("%s must be a 3-D array with dims [time, treatment, state].", name))
  }
  d <- dim(S)
  n_time   <- d[1]
  n_trt    <- d[2]
  n_state  <- d[3]

  dn <- dimnames(S)
  trt_names    <- if (!is.null(dn)) dn[[2]] else NULL
  state_names  <- if (!is.null(dn)) dn[[3]] else NULL

  for (k in seq_len(n_state)) {
    mat_k <- S[, , k, drop = FALSE][, , 1L]
    bad   <- colSums(mat_k < tol, na.rm = TRUE) > (n_time / 2)

    if (any(bad)) {
      bad_trt_idx  <- which(bad)
      trt_label    <- if (!is.null(trt_names)) paste(trt_names[bad_trt_idx], collapse = ", ")
      else paste(bad_trt_idx, collapse = ", ")
      state_label  <- if (!is.null(state_names)) state_names[k] else as.character(k)
      stop(sprintf(
        "%s validity failed: state %s, treatment arm(s) %s have >50%% of time points < %.1e.",
        name, state_label, trt_label, tol
      ))
    }
  }
  invisible(TRUE)
}

# ---- preprocess S: clamp [0,1] and enforce nonincreasing in t ---------------
#' @keywords internal
#' @noRd
.preprocess_S_array <- function(S) {
  S <- pmin(pmax(S, 0), 1)  # clamp to [0,1]
  for (q in seq_len(dim(S)[3])) {         # states
    for (c in seq_len(dim(S)[2])) {       # treatment arms
      for (t in 2:nrow(S)) {              # time (start at second)
        if (S[t, c, q] > S[t - 1, c, q]) {
          S[t, c, q] <- S[t - 1, c, q]
        }
      }
    }
  }
  S
}

# ---- interpolate S on [0, tau] for one tau ----------------------------------
#' @keywords internal
#' @noRd
.interpolate_S_array <- function(S, tau, time) {
  # S: [time × 2 × max_s], time: vector including 0
  if (!tau %in% time) {
    interp_time <- c(time[time <= tau], tau)
  } else {
    interp_time <- time[time <= tau]
  }
  arr <- array(
    apply(S, c(2, 3), function(x) stats::approx(time, x, xout = interp_time, rule = 2)$y),
    dim = c(length(interp_time), dim(S)[2], dim(S)[3])
  )
  list(S = arr, time = interp_time)
}

# ---- RMT-IF integrand integrals for one S (cluster or individual) ----------
#' @keywords internal
#' @noRd
.compute_for_S_array <- function(S, time) {
  max_s <- dim(S)[3]
  # Binary special case: max_s == 1 => NO minus overlap terms (RMST-like)
  if (max_s == 1L) {
    integral_1 <- pracma::trapz(time, S[, 1, 1] * 1)  # S1 * 1
    integral_2 <- pracma::trapz(time, S[, 2, 1] * 1)  # S0 * 1
    return(matrix(c(integral_1, integral_2), nrow = 2L,
                  dimnames = list(c("s1qs0q1","s0qs1q1"), "1")))
  }

  # General multi-state case (max_s >= 2)
  out <- sapply(seq_len(max_s), function(q) {
    if (q == max_s) {
      # last slice, compare to "absorbing" 1
      i1 <- pracma::trapz(time, S[, 1, q] * 1) -
        pracma::trapz(time, S[, 1, q] * S[, 2, q])
      i0 <- pracma::trapz(time, S[, 2, q] * 1) -
        pracma::trapz(time, S[, 2, q] * S[, 1, q])
    } else {
      i1 <- pracma::trapz(time, S[, 1, q] * S[, 2, q + 1]) -
        pracma::trapz(time, S[, 1, q] * S[, 2, q])
      i0 <- pracma::trapz(time, S[, 2, q] * S[, 1, q + 1]) -
        pracma::trapz(time, S[, 2, q] * S[, 1, q])
    }
    c(i1, i0)
  })
  rownames(out) <- c("s1qs0q1","s0qs1q1")
  colnames(out) <- as.character(seq_len(max_s))
  out
}

# ---- compute stage-wise integrals at one tau --------------------------------
#' @keywords internal
#' @noRd
.compute_rmtif_integrals_one_tau <- function(S_c, S_i, time, tau) {
  # S_c, S_i: [time × 2 × max_s] including t=0 row
  Sc_pre  <- .preprocess_S_array(S_c)
  Si_pre  <- .preprocess_S_array(S_i)

  # interpolate to [0, tau]
  tmp_c   <- .interpolate_S_array(Sc_pre, tau = tau, time = time)
  tmp_i   <- .interpolate_S_array(Si_pre, tau = tau, time = time)
  Sc_interp <- tmp_c$S;  time_interp <- tmp_c$time
  Si_interp <- tmp_i$S   # time_interp is same

  Sc_int <- as.data.frame(.compute_for_S_array(Sc_interp, time_interp))
  Si_int <- as.data.frame(.compute_for_S_array(Si_interp, time_interp))

  list(Sc_int = Sc_int, Si_int = Si_int)
}
