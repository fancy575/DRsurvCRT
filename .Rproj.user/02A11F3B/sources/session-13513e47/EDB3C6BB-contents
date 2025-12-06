#' @title Multi-state stage-specific survival core (internal)
#' @description
#' status ∈ {0,1,...,max_s}, 0 = censoring. For each state s ≥ 1, we:
#'   - keep rows with status == 0 or status ≥ s
#'   - define event = 1{status != 0}
#'   - fit nuisance models for A=1 and the flipped dataset
#'   - call C++ backends (marginal_est or frailty_est) on a common event-time grid
#'
#' Returns S_stage_cluster/S_stage_ind: [time_full × 2 × max_s] with a t=0 row of 1s.
#'
#' @keywords internal
#' @noRd
.DR_stage_surv_core <- function(data, formula, cens_formula, intv,
                                method = c("marginal","frailty"),
                                trt_prob = NULL,
                                fit_controls = NULL,
                                e_time = NULL) {

  method <- match.arg(method)

  # ensure Surv()/cluster() are visible
  environment(formula) <- .surv_env(parent = environment(formula) %||% parent.frame())
  if (!is.null(cens_formula))
    environment(cens_formula) <- .surv_env(parent = environment(cens_formula) %||% parent.frame())

  nm    <- .surv_lhs(formula)
  info  <- .extract_cluster(formula)
  clvar <- info$cluster
  rhs   <- info$rhs_wo_cluster

  if (is.null(clvar) || !nzchar(clvar))
    stop("Outcome formula must include cluster(<id>) for CRT.", call. = FALSE)

  time_var   <- nm$time
  status_var <- nm$status
  cl_col     <- clvar

  # sort by time
  data <- data[order(data[[time_var]]), , drop = FALSE]

  # treatment probabilities at cluster level
  if (is.null(trt_prob)) {
    pr <- .compute_trt_probs(data, trt = intv, strata_col = cl_col)
    p0 <- unname(pr["p0"])
    p1 <- unname(pr["p1"])
  } else {
    stopifnot(is.numeric(trt_prob), length(trt_prob) == 2L, all(is.finite(trt_prob)))
    p0 <- trt_prob[1]
    p1 <- trt_prob[2]
  }
  .validate_probs(trt_prob1 = p1, trt_prob0 = p0)

  # event-time grid (positive times only)
  all_etime <- if (is.null(e_time)) sort(unique(data[[time_var]])) else sort(unique(as.numeric(e_time)))
  all_etime <- all_etime[all_etime > 0]
  if (!length(all_etime))
    stop("No positive event times found.", call. = FALSE)

  max_s <- max(data[[status_var]], na.rm = TRUE)
  if (max_s < 1L)
    stop("Need at least one non-censoring state (status > 0).", call. = FALSE)

  # RHS including cluster() for nuisance fits
  rhs_surv <- rhs
  if (nzchar(rhs_surv))
    rhs_surv <- paste(rhs_surv, sprintf("+ cluster(%s)", cl_col))
  else
    rhs_surv <- sprintf("cluster(%s)", cl_col)

  nT  <- length(all_etime)
  S_c <- array(NA_real_, dim = c(nT, 2L, max_s))
  S_i <- array(NA_real_, dim = c(nT, 2L, max_s))

  for (s in seq_len(max_s)) {
    # stage s: keep censored or status >= s, define event = 1{status != 0}
    sub_dat <- data[
      data[[status_var]] == 0L | data[[status_var]] >= s,
      , drop = FALSE
    ]
    sub_dat$event <- as.integer(sub_dat[[status_var]] != 0L)

    # event and censoring formulas for this stage
    f_event_s <- as.formula(
      sprintf("Surv(%s, event) ~ %s", time_var, rhs_surv)
    )
    f_cens_s <- if (is.null(cens_formula)) {
      # stage-specific censoring model: Surv(time, event==0) ~ rhs_surv
      as.formula(sprintf("Surv(%s, event==0) ~ %s", time_var, rhs_surv))
    } else {
      cens_formula  # if user insists; will still use event as delta in C++ calls
    }

    environment(f_event_s) <- .surv_env(parent = environment(formula))
    environment(f_cens_s)  <- .surv_env(parent = environment(formula))

    # A=1 subset
    sub1 <- sub_dat[sub_dat[[intv]] == 1L, , drop = FALSE]

    if (method == "marginal") {
      fit_e1 <- survival::coxph(f_event_s, data = sub1)
      fit_c1 <- survival::coxph(f_cens_s,  data = sub1)
    } else {
      # frailty
      if (is.null(fit_controls)) {
        fit_controls <- frailtyEM::emfrail_control(
          se       = FALSE, se_adj = FALSE, ca_test = FALSE,
          lik_ci   = FALSE, zph    = FALSE,
          em_control  = list(eps = 1e-6, maxit = 200, fast_fit = TRUE),
          nlm_control = list(iterlim = 50)
        )
      }
      fit_e1 <- frailtyEM::emfrail(
        f_event_s, data = sub1,
        distribution = frailtyEM::emfrail_dist("gamma", theta = 2),
        control      = fit_controls
      )
      fit_c1 <- frailtyEM::emfrail(
        f_cens_s, data = sub1,
        distribution = frailtyEM::emfrail_dist("gamma", theta = 2),
        control      = fit_controls
      )
    }

    # flip treatment
    sub_dat0 <- sub_dat
    sub_dat0[[intv]] <- 1L - sub_dat0[[intv]]

    sub0 <- sub_dat0[sub_dat0[[intv]] == 1L, , drop = FALSE]

    if (method == "marginal") {
      fit_e0 <- survival::coxph(f_event_s, data = sub0)
      fit_c0 <- survival::coxph(f_cens_s,  data = sub0)
    } else {
      fit_e0 <- frailtyEM::emfrail(
        f_event_s, data = sub0,
        distribution = frailtyEM::emfrail_dist("gamma", theta = 2),
        control      = fit_controls
      )
      fit_c0 <- frailtyEM::emfrail(
        f_cens_s, data = sub0,
        distribution = frailtyEM::emfrail_dist("gamma", theta = 2),
        control      = fit_controls
      )
    }

    # design matrices WITHOUT cluster()
    XsXc <- .model_mats_safe(f_event_s, f_cens_s, sub_dat)

    # call C++ backends
    if (method == "marginal") {
      raw_s <- marginal_est(
        ftime       = sub_dat[[time_var]],
        delta       = sub_dat[["event"]],
        trt         = sub_dat[[intv]],
        strata      = sub_dat[[cl_col]],
        trt_prob1   = p1,
        trt_prob0   = p0,
        censor_cov  = XsXc$Xc,
        surv_cov    = XsXc$Xs,
        censor_fit1 = stats::coef(fit_c1),  censor_fit0 = stats::coef(fit_c0),
        surv_fit1   = stats::coef(fit_e1),  surv_fit0   = stats::coef(fit_e0),
        e_time      = all_etime,
        RMST_cal    = FALSE
      )
    } else {
      raw_s <- frailty_est(
        ftime       = sub_dat[[time_var]],
        delta       = sub_dat[["event"]],
        trt         = sub_dat[[intv]],
        strata      = sub_dat[[cl_col]],
        trt_prob1   = p1,
        trt_prob0   = p0,
        censor_cov  = XsXc$Xc,
        surv_cov    = XsXc$Xs,
        censor_fit1 = stats::coef(fit_c1),  censor_fit0  = stats::coef(fit_c0),
        beta_c1     = exp(fit_c1$logtheta), beta_c0      = exp(fit_c0$logtheta),
        surv_fit1   = stats::coef(fit_e1),  surv_fit0    = stats::coef(fit_e0),
        beta_s1     = exp(fit_e1$logtheta), beta_s0      = exp(fit_e0$logtheta),
        e_time      = all_etime,
        RMST_cal    = FALSE
      )
    }

    S_c[ , , s] <- as.matrix(raw_s$S_cluster)[, 1:2, drop = FALSE]
    S_i[ , , s] <- as.matrix(raw_s$S_individual)[, 1:2, drop = FALSE]
  }

  # validity checks and prepend t = 0 row of 1s
  .check_S_array(S_c, name = "S_stage_cluster")
  .check_S_array(S_i, name = "S_stage_ind")

  time_full <- c(0, all_etime)
  max_s     <- dim(S_c)[3]

  S_c0 <- abind::abind(array(1, dim = c(1, 2, max_s)), S_c, along = 1)
  S_i0 <- abind::abind(array(1, dim = c(1, 2, max_s)), S_i, along = 1)

  list(
    trt_prob        = c(p0 = p0, p1 = p1),
    event_time_core = all_etime,   # positive times only
    event_time_full = time_full,   # includes 0
    S_stage_cluster = S_c0,        # [time_full × 2 × max_s]
    S_stage_ind     = S_i0,
    cluster_col     = cl_col,
    max_state       = max_s
  )
}


#' @title Point-estimation core (internal)
#' @description
#' Wrapper around .DR_stage_surv_core() that returns either:
#'
#' * SPCE (stage-specific survival): S_stage_cluster/S_stage_ind
#' * RMTIF (generalized RMST-type win estimand): RMTIF_cluster/RMTIF_ind
#'   evaluated at every positive event time τ on the grid.
#'
#' Binary case (max_state == 1) is a special case where the RMT-IF
#' reduces to RMST-type integrals (no minus overlap terms).
#'
#' @keywords internal
#' @noRd
.DR_est_core <- function(data, formula, cens_formula, intv,
                         method   = c("marginal","frailty"),
                         estimand = c("SPCE","RMTIF"),
                         trt_prob = NULL,
                         fit_controls = NULL,
                         e_time   = NULL) {

  method   <- match.arg(method)
  estimand <- match.arg(estimand)

  stage <- .DR_stage_surv_core(
    data         = data,
    formula      = formula,
    cens_formula = cens_formula,
    intv         = intv,
    method       = method,
    trt_prob     = trt_prob,
    fit_controls = fit_controls,
    e_time       = e_time
  )

  p0        <- stage$trt_prob["p0"]
  p1        <- stage$trt_prob["p1"]
  time_core <- stage$event_time_core   # positive event times
  time_full <- stage$event_time_full   # 0 + time_core
  S_c0      <- stage$S_stage_cluster
  S_i0      <- stage$S_stage_ind
  cl_col    <- stage$cluster_col
  max_s     <- stage$max_state

  if (estimand == "SPCE") {
    # Just return stage-specific survival
    return(list(
      trt_prob          = c(p0 = p0, p1 = p1),
      event_time        = time_full,       # includes t=0
      S_stage_cluster   = S_c0,            # [time_full × 2 × max_s]
      S_stage_ind       = S_i0,
      RMTIF_cluster     = NULL,
      RMTIF_ind         = NULL,
      stagewise_cluster = NULL,
      stagewise_ind     = NULL,
      cluster_col       = cl_col,
      max_state         = max_s
    ))
  }

  # ---- RMT-IF branch: integrals at each τ in time_core ----------------------
  tau_vec <- time_core
  nL      <- length(tau_vec)

  stagewise_cl <- vector("list", nL)
  stagewise_in <- vector("list", nL)

  R_cl <- matrix(NA_real_, nrow = nL, ncol = 3L)
  R_in <- matrix(NA_real_, nrow = nL, ncol = 3L)
  colnames(R_cl) <- colnames(R_in) <- c("R1","R0","R1-R0")
  rownames(R_cl) <- rownames(R_in) <- tau_vec

  for (k in seq_len(nL)) {
    tau <- tau_vec[k]
    tmp <- .compute_rmtif_integrals_one_tau(S_c0, S_i0, time_full, tau)

    Sc_int <- as.matrix(tmp$Sc_int)  # 2 × max_s (or 2 × 1 in binary case)
    Si_int <- as.matrix(tmp$Si_int)

    # add "sum" col and diff row, like your DR_win_est
    Sc_int <- rbind(
      Sc_int,
      diff = Sc_int[1, ] - Sc_int[2, ]
    )
    Sc_int <- cbind(Sc_int,
                    sum = c(sum(Sc_int[1, ], na.rm = TRUE),
                            sum(Sc_int[2, ], na.rm = TRUE),
                            sum(Sc_int[3, ], na.rm = TRUE)))
    rownames(Sc_int)[1:2] <- c("s1qs0qp1","s0qs1qp1")
    colnames(Sc_int) <- c(paste0("stage_", seq_len(max_s)), "sum")

    Si_int <- rbind(
      Si_int,
      diff = Si_int[1, ] - Si_int[2, ]
    )
    Si_int <- cbind(Si_int,
                    sum = c(sum(Si_int[1, ], na.rm = TRUE),
                            sum(Si_int[2, ], na.rm = TRUE),
                            sum(Si_int[3, ], na.rm = TRUE)))
    rownames(Si_int)[1:2] <- c("s1qs0qp1","s0qs1qp1")
    colnames(Si_int) <- c(paste0("stage_", seq_len(max_s)), "sum")

    stagewise_cl[[k]] <- Sc_int
    stagewise_in[[k]] <- Si_int

    R_cl[k, ] <- c(Sc_int["s1qs0qp1","sum"],
                   Sc_int["s0qs1qp1","sum"],
                   Sc_int["diff","sum"])
    R_in[k, ] <- c(Si_int["s1qs0qp1","sum"],
                   Si_int["s0qs1qp1","sum"],
                   Si_int["diff","sum"])
  }

  list(
    trt_prob          = c(p0 = p0, p1 = p1),
    event_time        = tau_vec,          # τ grid (positive event times)
    S_stage_cluster   = S_c0,             # still available for plotting SPCE-type curves
    S_stage_ind       = S_i0,
    RMTIF_cluster     = R_cl,
    RMTIF_ind         = R_in,
    stagewise_cluster = stagewise_cl,
    stagewise_ind     = stagewise_in,
    cluster_col       = cl_col,
    max_state         = max_s
  )
}

#' @title Jackknife variance core (internal)
#' @description
#' Leave-one-cluster-out jackknife for:
#'   * SPCE: stage-specific survival S_s(t) by arm and diff
#'   * RMTIF: generalized RMST-type win curves R1(τ), R0(τ), R1-R0(τ)
#'
#' Works for method = "marginal" or "frailty".
#'
#' @keywords internal
#' @noRd
.DR_var_jackknife <- function(data, formula, cens_formula, intv,
                              method    = c("marginal","frailty"),
                              estimand  = c("SPCE","RMTIF"),
                              trt_prob  = NULL,
                              fit_controls = NULL,
                              e_time_full  = NULL) {

  method   <- match.arg(method)
  estimand <- match.arg(estimand)

  environment(formula) <- .surv_env(parent = environment(formula) %||% parent.frame())
  if (!is.null(cens_formula))
    environment(cens_formula) <- .surv_env(parent = environment(cens_formula) %||% parent.frame())

  nm   <- .surv_lhs(formula)
  data <- data[order(data[[nm$time]]), , drop = FALSE]

  if (is.null(e_time_full))
    e_time_full <- sort(unique(data[[nm$time]]))
  e_time_full <- as.numeric(e_time_full)
  e_time_full <- e_time_full[e_time_full > 0]

  info  <- .extract_cluster(formula)
  clvar <- info$cluster
  if (is.null(clvar))
    stop("Need cluster() in formula for jackknife.", call. = FALSE)

  clusters <- unique(data[[clvar]])
  K        <- length(clusters)
  JKfac    <- (K - 1) / K

  # full-data fit for shapes
  est0 <- .DR_est_core(
    data         = data,
    formula      = formula,
    cens_formula = cens_formula,
    intv         = intv,
    method       = method,
    estimand     = estimand,
    trt_prob     = trt_prob,
    fit_controls = fit_controls,
    e_time       = e_time_full
  )

  max_s <- est0$max_state

  if (estimand == "SPCE") {
    # est0$S_stage_cluster is [T_full+1 × 2 × max_s], drop t=0 row
    S0_cl  <- est0$S_stage_cluster[-1, , , drop = FALSE]
    S0_ind <- est0$S_stage_ind[-1,  , , drop = FALSE]
    T_full <- dim(S0_cl)[1]

    # collect replicate S1,S0,S1-S0 for each cluster, time, state
    agg_cl  <- array(NA_real_, dim = c(T_full, 3L, max_s, K))
    agg_ind <- array(NA_real_, dim = c(T_full, 3L, max_s, K))

    for (m in seq_along(clusters)) {
      keep <- data[data[[clvar]] != clusters[m], , drop = FALSE]

      est_m <- .DR_est_core(
        data         = keep,
        formula      = formula,
        cens_formula = cens_formula,
        intv         = intv,
        method       = method,
        estimand     = "SPCE",
        trt_prob     = trt_prob,
        fit_controls = fit_controls,
        e_time       = e_time_full
      )

      Scl_m  <- est_m$S_stage_cluster[-1, , , drop = FALSE]
      Sind_m <- est_m$S_stage_ind[-1,  , , drop = FALSE]

      for (s in seq_len(max_s)) {
        S1_cl <- Scl_m[, 1, s]
        S0_cl <- Scl_m[, 2, s]
        S1_in <- Sind_m[, 1, s]
        S0_in <- Sind_m[, 2, s]

        agg_cl[ , , s, m]  <- cbind(S1_cl, S0_cl, S1_cl - S0_cl)
        agg_ind[ , , s, m] <- cbind(S1_in, S0_in, S1_in - S0_in)
      }
    }

    var_stage_cl  <- array(NA_real_, dim = c(T_full, 3L, max_s))
    var_stage_ind <- array(NA_real_, dim = c(T_full, 3L, max_s))

    for (t in seq_len(T_full)) {
      for (j in 1:3) {
        for (s in seq_len(max_s)) {
          v_cl <- agg_cl[t, j, s, ]
          v_in <- agg_ind[t, j, s, ]

          mu_cl <- mean(v_cl, na.rm = TRUE)
          mu_in <- mean(v_in, na.rm = TRUE)

          var_stage_cl[t, j, s]  <- JKfac * sum((v_cl - mu_cl)^2, na.rm = TRUE)
          var_stage_ind[t, j, s] <- JKfac * sum((v_in - mu_in)^2, na.rm = TRUE)
        }
      }
    }

    dimnames(var_stage_cl) <- list(
      time   = e_time_full,
      comp   = c("Var(S1)", "Var(S0)", "Var(S1-S0)"),
      state  = paste0("state_", seq_len(max_s))
    )
    dimnames(var_stage_ind) <- dimnames(var_stage_cl)

    return(list(
      Cluster = list(
        var_stage  = var_stage_cl,
        event_time = e_time_full
      ),
      Individual = list(
        var_stage  = var_stage_ind,
        event_time = e_time_full
      )
    ))
  }

  # ---- RMT-IF jackknife -----------------------------------------------------
  R0_cl  <- est0$RMTIF_cluster   # [T_full × 3]
  T_full <- nrow(R0_cl)

  agg_R_cl  <- array(NA_real_, dim = c(T_full, 3L, K))
  agg_R_ind <- array(NA_real_, dim = c(T_full, 3L, K))

  for (m in seq_along(clusters)) {
    keep <- data[data[[clvar]] != clusters[m], , drop = FALSE]

    est_m <- .DR_est_core(
      data         = keep,
      formula      = formula,
      cens_formula = cens_formula,
      intv         = intv,
      method       = method,
      estimand     = "RMTIF",
      trt_prob     = trt_prob,
      fit_controls = fit_controls,
      e_time       = e_time_full
    )

    Rcl_m  <- est_m$RMTIF_cluster  # [time × 3]
    Rind_m <- est_m$RMTIF_ind

    agg_R_cl[ , , m]  <- Rcl_m
    agg_R_ind[ , , m] <- Rind_m
  }

  mean_R_cl  <- apply(agg_R_cl,  c(1, 2), mean, na.rm = TRUE)
  mean_R_ind <- apply(agg_R_ind, c(1, 2), mean, na.rm = TRUE)

  r_cl  <- sweep(agg_R_cl,  c(1, 2), mean_R_cl,  "-")
  r_ind <- sweep(agg_R_ind, c(1, 2), mean_R_ind, "-")

  var_R_cl  <- JKfac * apply(r_cl^2,  c(1, 2), sum, na.rm = TRUE)
  var_R_ind <- JKfac * apply(r_ind^2, c(1, 2), sum, na.rm = TRUE)

  cov_R_cl  <- JKfac * apply(r_cl[, 1, ] * r_cl[, 2, ], 1, sum, na.rm = TRUE)
  cov_R_ind <- JKfac * apply(r_ind[, 1, ] * r_ind[, 2, ], 1, sum, na.rm = TRUE)

  out_cl  <- cbind(var_R_cl, cov = cov_R_cl)
  out_ind <- cbind(var_R_ind, cov = cov_R_ind)

  rownames(out_cl)  <- e_time_full
  rownames(out_ind) <- e_time_full
  colnames(out_cl)  <- c("Var(R1)","Var(R0)","Var(R1-R0)","Cov(R1,R0)")
  colnames(out_ind) <- colnames(out_cl)

  list(
    Cluster = list(
      var_R      = out_cl,
      event_time = e_time_full
    ),
    Individual = list(
      var_R      = out_ind,
      event_time = e_time_full
    )
  )
}


