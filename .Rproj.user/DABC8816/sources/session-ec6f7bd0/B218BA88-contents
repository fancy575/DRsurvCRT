#' @title Point-estimation core (internal)
#' @description
#' Fits nuisance models on treated subset and the flipped dataset, then calls
#' C++ backends to obtain **full** grids (S and RMST) on the event-time grid.
#' If `tau` is supplied, also computes RMST-at-τ by the linear patch.
#' No printing here; downstream functions slice as needed.
#' @keywords internal
#' @noRd
.DR_est_core <- function(data, formula, cens_formula, intv,
                         method = c("marginal","frailty"),
                         trt_prob = NULL,
                         fit_controls = NULL,
                         e_time = NULL,
                         strata = NULL,
                         tau = NULL) {

  method <- match.arg(method)

  # ensure Surv()/cluster() are visible
  environment(formula) <- .surv_env(parent = environment(formula) %||% parent.frame())
  if (!is.null(cens_formula))
    environment(cens_formula) <- .surv_env(parent = environment(cens_formula) %||% parent.frame())

  nm    <- .surv_lhs(formula)
  info  <- .extract_cluster(formula)
  clvar <- info$cluster %||% strata
  rhs   <- info$rhs_wo_cluster

  if (method == "frailty" && is.null(clvar))
    stop("frailty method requires cluster() in formula or a 'strata' column name.", call. = FALSE)

  # build modeling formulas (attach cluster() when available)
  rhs_surv <- rhs
  if (!is.null(clvar) && nzchar(clvar))
    rhs_surv <- paste(rhs_surv, sprintf("cluster(%s)", clvar), sep = if (nzchar(rhs_surv)) " + " else "")
  if (is.null(cens_formula)) {
    cf <- .build_cens_formula_from(nm, rhs, clvar)
    cens_formula <- as.formula(cf)
  }
  f_event <- as.formula(sprintf("Surv(%s,%s) ~ %s", nm$time, nm$status, rhs_surv))

  # fixed event-time grid
  data <- data[order(data[[nm$time]]), , drop = FALSE]
  et   <- if (is.null(e_time)) sort(unique(data[[nm$time]])) else as.numeric(e_time)

  # cluster id column used by backends
  cl_col <- clvar %||% strata
  if (is.null(cl_col))
    stop("Need a cluster identifier (cluster() in formula or 'strata=').", call. = FALSE)

  # treatment probabilities
  if (is.null(trt_prob)) {
    pr <- .compute_trt_probs(data, trt = intv, strata_col = cl_col)
    p0 <- unname(pr["p0"]); p1 <- unname(pr["p1"])
  } else {
    stopifnot(is.numeric(trt_prob), length(trt_prob) == 2L, all(is.finite(trt_prob)))
    p0 <- trt_prob[1]; p1 <- trt_prob[2]
  }
  .validate_probs(p1, p0)

  # fits: treated subset, then flip and refit
  sub1 <- data[data[[intv]] == 1, , drop = FALSE]
  if (method == "marginal") {
    fit_e1 <- survival::coxph(f_event, data = sub1)
    fit_c1 <- survival::coxph(cens_formula, data = sub1)
  } else {
    if (is.null(fit_controls)) {
      fit_controls <- frailtyEM::emfrail_control(
        se = FALSE, se_adj = FALSE, ca_test = FALSE, lik_ci = FALSE, zph = FALSE,
        em_control  = list(eps = 1e-6, maxit = 200, fast_fit = TRUE),
        nlm_control = list(iterlim = 50)
      )
    }
    fit_e1 <- frailtyEM::emfrail(f_event, data = sub1,
                                 distribution = frailtyEM::emfrail_dist("gamma", theta = 2),
                                 control = fit_controls)
    fit_c1 <- frailtyEM::emfrail(cens_formula, data = sub1,
                                 distribution = frailtyEM::emfrail_dist("gamma", theta = 2),
                                 control = fit_controls)
  }
  data0 <- data; data0[[intv]] <- 1 - data0[[intv]]
  sub0  <- data0[data0[[intv]] == 1, , drop = FALSE]
  if (method == "marginal") {
    fit_e0 <- survival::coxph(f_event, data = sub0)
    fit_c0 <- survival::coxph(cens_formula, data = sub0)
  } else {
    fit_e0 <- frailtyEM::emfrail(f_event, data = sub0,
                                 distribution = frailtyEM::emfrail_dist("gamma", theta = 2),
                                 control = fit_controls)
    fit_c0 <- frailtyEM::emfrail(cens_formula, data = sub0,
                                 distribution = frailtyEM::emfrail_dist("gamma", theta = 2),
                                 control = fit_controls)
  }

  # design matrices WITHOUT cluster()
  XsXc <- .model_mats_safe(f_event, cens_formula, data)

  # call C++ backends
  if (method == "marginal") {
    raw <- marginal_est(
      ftime       = data[[nm$time]],
      delta       = data[[nm$status]],
      trt         = data[[intv]],
      strata      = data[[cl_col]],
      trt_prob1   = p1,
      trt_prob0   = p0,
      censor_cov  = XsXc$Xc,
      surv_cov    = XsXc$Xs,
      censor_fit1 = stats::coef(fit_c1),  censor_fit0 = stats::coef(fit_c0),
      surv_fit1   = stats::coef(fit_e1),  surv_fit0   = stats::coef(fit_e0),
      e_time      = et,
      RMST_cal    = TRUE
    )
  } else {
    raw <- frailty_est(
      ftime       = data[[nm$time]],
      delta       = data[[nm$status]],
      trt         = data[[intv]],
      strata      = data[[cl_col]],
      trt_prob1   = p1,
      trt_prob0   = p0,
      censor_cov  = XsXc$Xc,
      surv_cov    = XsXc$Xs,
      censor_fit1 = stats::coef(fit_c1),  censor_fit0  = stats::coef(fit_c0),
      beta_c1     = exp(fit_c1$logtheta), beta_c0      = exp(fit_c0$logtheta),
      surv_fit1   = stats::coef(fit_e1),  surv_fit0    = stats::coef(fit_e0),
      beta_s1     = exp(fit_e1$logtheta), beta_s0      = exp(fit_e0$logtheta),
      e_time      = et,
      RMST_cal    = TRUE
    )
  }

  et_out <- as.numeric(raw$event_time)
  S_cl   <- as.matrix(raw$S_cluster);     rownames(S_cl)  <- et_out
  S_ind  <- as.matrix(raw$S_individual);  rownames(S_ind) <- et_out
  .check_S_warn(S_cl,  "S_cluster")
  .check_S_warn(S_ind, "S_individual")

  # RMST at tau (optional, computed here—NOT in variance)
  RMST_tau_cluster <- RMST_tau_ind <- NULL
  if (!is.null(tau)) {
    idx   <- pmax(1L, findInterval(tau, et_out))
    tdiff <- tau - et_out[idx]
    Scl_sub  <- S_cl[idx,  , drop = FALSE]
    Sind_sub <- S_ind[idx, , drop = FALSE]
    RMST_tau_cluster <- as.matrix(raw$RMST_cluster_out)[idx, , drop = FALSE] + sweep(Scl_sub,  1L, tdiff, `*`)
    RMST_tau_ind     <- as.matrix(raw$RMST_ind_out)[idx,     , drop = FALSE] + sweep(Sind_sub, 1L, tdiff, `*`)
    rownames(RMST_tau_cluster) <- rownames(RMST_tau_ind) <- tau
  }

  list(
    trt_prob            = c(p0 = p0, p1 = p1),
    event_time          = et_out,
    S_cluster_full      = S_cl,
    S_ind_full          = S_ind,
    RMST_cluster_full   = as.matrix(raw$RMST_cluster_out),
    RMST_ind_full       = as.matrix(raw$RMST_ind_out),
    RMST_tau_cluster    = RMST_tau_cluster,   # NULL if tau is NULL
    RMST_tau_ind        = RMST_tau_ind,       # NULL if tau is NULL
    cluster_col         = cl_col
  )
}

#' @title Jackknife variance core (internal)
#' @description
#' Leave-one-cluster-out jackknife on a **shared** event-time grid. For SPCE,
#' aggregation is over the full grid; for RMST, we re-use RMST-at-τ from the core
#' (no re-patching here).
#' @keywords internal
#' @noRd
.DR_var_jackknife <- function(data, formula, cens_formula, intv,
                              method = c("marginal","frailty"),
                              trt_prob = NULL,
                              e_time_full = NULL,
                              tau = NULL,
                              strata = NULL,
                              estimand = c("SPCE","RMST")) {
  method   <- match.arg(method)
  estimand <- match.arg(estimand)

  environment(formula) <- .surv_env(parent = environment(formula) %||% parent.frame())
  if (!is.null(cens_formula))
    environment(cens_formula) <- .surv_env(parent = environment(cens_formula) %||% parent.frame())

  nm <- .surv_lhs(formula)
  data <- data[order(data[[nm$time]]), , drop = FALSE]

  if (is.null(e_time_full)) e_time_full <- sort(unique(data[[nm$time]]))
  e_time_full <- as.numeric(e_time_full)

  info  <- .extract_cluster(formula)
  clvar <- info$cluster %||% strata
  if (is.null(clvar)) stop("Need cluster() in formula or 'strata' for jackknife.", call. = FALSE)

  clusters <- unique(data[[clvar]])
  K <- length(clusters)
  JKfac <- (K - 1) / K

  if (estimand == "SPCE") {
    nT <- length(e_time_full)
    agg_cl  <- array(NA_real_, dim = c(nT, 3, K))
    agg_ind <- array(NA_real_, dim = c(nT, 3, K))
  } else {
    if (is.null(tau) || !length(tau))
      stop("Provide non-empty 'tau' for RMST jackknife.", call. = FALSE)
    nL <- length(tau)
    agg_cl  <- array(NA_real_, dim = c(nL, 3, K))
    agg_ind <- array(NA_real_, dim = c(nL, 3, K))
  }

  for (m in seq_along(clusters)) {
    keep <- data[data[[clvar]] != clusters[m], , drop = FALSE]
    est_m <- .DR_est_core(keep, formula, cens_formula, intv, method,
                          trt_prob, fit_controls = NULL,
                          e_time = e_time_full, strata = clvar, tau = tau)
    if (estimand == "SPCE") {
      agg_cl[ , , m]  <- est_m$S_cluster_full
      agg_ind[ , , m] <- est_m$S_ind_full
    } else {
      agg_cl[ , , m]  <- est_m$RMST_tau_cluster
      agg_ind[ , , m] <- est_m$RMST_tau_ind
    }
  }

  mean_cl  <- apply(agg_cl,  c(1,2), mean, na.rm = TRUE)
  mean_ind <- apply(agg_ind, c(1,2), mean, na.rm = TRUE)
  r_cl  <- sweep(agg_cl,  c(1,2), mean_cl,  "-")
  r_ind <- sweep(agg_ind, c(1,2), mean_ind, "-")
  var_cl  <- JKfac * apply(r_cl^2,  c(1,2), sum, na.rm = TRUE)
  var_ind <- JKfac * apply(r_ind^2, c(1,2), sum, na.rm = TRUE)
  cov_cl  <- JKfac * apply(r_cl[,1,] * r_cl[,2,], 1, sum, na.rm = TRUE)
  cov_ind <- JKfac * apply(r_ind[,1,] * r_ind[,2,], 1, sum, na.rm = TRUE)

  list(
    Cluster    = cbind(var_cl, cov = cov_cl),
    Individual = cbind(var_ind, cov = cov_ind),
    event_time = e_time_full
  )
}
