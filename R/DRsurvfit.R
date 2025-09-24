#' Doubly-robust survival estimation for CRTs (SPCE or RMST)
#'
#' @description
#' Fits doubly-robust estimators for cluster-randomized trials with right-censored
#' outcomes. SPCE is returned on the **full event-time grid**. RMST is returned at
#' user-specified `tau` (vector). Jackknife variance uses a **shared** event-time grid
#' across leave-one-cluster-out fits and **reuses** RMST-at-τ computed in the core.
#'
#' The returned object includes metadata: final fitted outcome/censoring formulas,
#' the cluster id column, number of clusters, df for jackknife t-intervals (= K-1),
#' sample sizes, and cluster-level treatment split.
#'
#' @param data A data.frame.
#' @param formula Outcome model: e.g., \code{Surv(time, event) ~ W1 + W2 + Z1 + Z2 + cluster(M)}.
#'   - For \code{method = "frailty"}, a \code{cluster(<id>)} term is **required** (or supply \code{strata=}).
#'   - For \code{method = "marginal"}, the cluster term is optional; if you also pass \code{strata=},
#'     \code{cluster(strata)} is appended to both outcome and censoring formulas.
#' @param cens_formula Optional censoring model. If \code{NULL}, reuse the RHS (without \code{cluster()})
#'   and set LHS to \code{Surv(time, event == 0)}. If a cluster id is known, it is appended
#'   as \code{+ cluster(id)}.
#' @param intv Character: cluster-level treatment column (0/1), constant within cluster.
#' @param method \code{"marginal"} or \code{"frailty"}.
#' @param estimand \code{"SPCE"} or \code{"RMST"}.
#' @param tau Numeric vector of RMST landmark times (required for RMST).
#' @param trt_prob Optional length-2 numeric \code{(p0, p1)}. If \code{NULL}, computed from first row per cluster.
#' @param variance \code{"none"} or \code{"jackknife"}.
#' @param fit_controls Optional \code{frailtyEM::emfrail_control()} list (frailty method).
#' @param strata Optional cluster id column name to use if the outcome formula does not include \code{cluster()}.
#' @param spce_summary_times Optional numeric vector of times to show in \code{summary()} for SPCE.
#' @param verbose Logical.
#'
#' @return An object of class \code{DRsurvfit}.
#'         Key fields: \code{event_time}, \code{S_full_cluster}, \code{S_full_ind},
#'         \code{S_cluster}/\code{S_individual} (SPCE), \code{RMST_cluster}/\code{RMST_ind} (RMST),
#'         \code{trt_prob}, \code{var_cluster}/\code{var_ind} (if jackknife),
#'         \code{formula_outcome}, \code{cens_formula}, \code{cluster_col}, \code{n_clusters}, \code{df_jackknife},
#'         \code{n_obs}, \code{n_events}, \code{cluster_trt_counts}.
#'
#' @examples
#' \dontrun{
#' data(dat)
#' fit_spce <- DRsurvfit(
#'   dat,
#'   survival::Surv(time, event) ~ W1 + W2 + Z1 + Z2 + cluster(M),
#'   intv = "A",
#'   method = "marginal",
#'   estimand = "SPCE",
#'   variance = "jackknife",
#'   spce_summary_times = c(0.5, 1, 1.5, 2)
#' )
#' summary(fit_spce, times = c(0.5, 1, 1.5, 2))
#'
#' fit_rmst <- DRsurvfit(
#'   dat,
#'   survival::Surv(time, event) ~ W1 + W2 + Z1 + Z2 + cluster(M),
#'   intv = "A",
#'   method = "frailty",
#'   estimand = "RMST",
#'   tau = c(0.5, 1, 2),
#'   variance = "jackknife"
#' )
#' summary(fit_rmst)
#' }
#' @export
DRsurvfit <- function(data,
                      formula,
                      cens_formula = NULL,
                      intv,
                      method   = c("marginal", "frailty"),
                      estimand = c("SPCE", "RMST"),
                      tau = NULL,
                      trt_prob = NULL,
                      variance = c("none","jackknife"),
                      fit_controls = NULL,
                      strata = NULL,
                      spce_summary_times = NULL,
                      verbose = FALSE) {

  method   <- match.arg(method)
  estimand <- match.arg(estimand)
  variance <- match.arg(variance)
  if (!is.data.frame(data)) stop("'data' must be a data.frame.", call. = FALSE)

  # fit core on full grid; pass tau so core also returns RMST-at-τ (if any)
  est <- .DR_est_core(
    data = data, formula = formula, cens_formula = cens_formula, intv = intv,
    method = method, trt_prob = trt_prob, fit_controls = fit_controls,
    e_time = NULL, strata = strata, tau = if (estimand == "RMST") tau else NULL
  )

  et   <- est$event_time
  S_cl <- est$S_cluster_full
  S_in <- est$S_ind_full

  # --- NEW: gather metadata for the returned object --------------------------
  # We reconstruct the final fitted formulas and cluster info in the same way as the core.
  environment(formula) <- .surv_env(parent = environment(formula) %||% parent.frame())
  if (!is.null(cens_formula))
    environment(cens_formula) <- .surv_env(parent = environment(cens_formula) %||% parent.frame())

  nm    <- .surv_lhs(formula)
  info  <- .extract_cluster(formula)
  clvar <- info$cluster %||% strata
  rhs   <- info$rhs_wo_cluster
  rhs_surv <- rhs
  if (!is.null(clvar) && nzchar(clvar))
    rhs_surv <- paste(rhs_surv, sprintf("cluster(%s)", clvar), sep = if (nzchar(rhs_surv)) " + " else "")
  if (is.null(cens_formula)) {
    cf <- .build_cens_formula_from(nm, rhs, clvar)
    cens_formula_final <- as.formula(cf)
  } else {
    cens_formula_final <- cens_formula
  }
  formula_outcome_final <- as.formula(sprintf("Surv(%s,%s) ~ %s", nm$time, nm$status, rhs_surv))

  cluster_col <- est$cluster_col
  n_obs       <- nrow(data)
  n_events    <- if (nm$status %in% names(data)) sum(as.integer(data[[nm$status]] == 1L), na.rm = TRUE) else NA_integer_
  K           <- length(unique(data[[cluster_col]]))
  df_jack     <- max(1L, K - 1L)

  # cluster-level treatment split by first row per cluster (same rule used to compute trt_prob)
  cl_ids <- unique(data[[cluster_col]])
  A_first <- vapply(cl_ids, function(g) data[[intv]][which(data[[cluster_col]] == g)[1]], numeric(1))
  cluster_trt_counts <- c(n_trt0 = sum(A_first == 0, na.rm = TRUE),
                          n_trt1 = sum(A_first == 1, na.rm = TRUE))
  # ---------------------------------------------------------------------------

  if (estimand == "SPCE") {
    var_out <- if (variance == "jackknife") {
      .DR_var_jackknife(
        data, formula, cens_formula, intv, method,
        trt_prob, e_time_full = et, tau = NULL,
        strata = est$cluster_col, estimand = "SPCE"
      )
    } else NULL

    ans <- list(
      method   = method,
      estimand = "SPCE",
      event_time      = et,
      S_full_cluster  = S_cl,
      S_full_ind      = S_in,
      S_cluster       = S_cl,
      S_individual    = S_in,
      RMST_cluster    = NULL,
      RMST_ind        = NULL,
      trt_prob        = unname(est$trt_prob),
      var_cluster     = if (!is.null(var_out)) var_out$Cluster else NULL,
      var_ind         = if (!is.null(var_out)) var_out$Individual else NULL,
      spce_summary_times = spce_summary_times,

      # NEW: metadata for rich summaries/repro
      formula_outcome = formula_outcome_final,
      cens_formula    = cens_formula_final,
      cluster_col     = cluster_col,
      n_clusters      = K,
      df_jackknife    = df_jack,
      n_obs           = n_obs,
      n_events        = n_events,
      cluster_trt_counts = cluster_trt_counts,
      call            = match.call(),
      jackknife       = identical(variance, "jackknife")
    )

  } else { # RMST
    if (is.null(tau) || !length(tau) || any(!is.finite(tau)))
      stop("'tau' must be a non-empty numeric vector for RMST.", call. = FALSE)

    R_cl <- est$RMST_tau_cluster
    R_in <- est$RMST_tau_ind
    rownames(R_cl) <- rownames(R_in) <- tau

    var_out <- if (variance == "jackknife") {
      .DR_var_jackknife(
        data, formula, cens_formula, intv, method,
        trt_prob, e_time_full = et, tau = tau,
        strata = est$cluster_col, estimand = "RMST"
      )
    } else NULL

    ans <- list(
      method   = method,
      estimand = "RMST",
      event_time      = et,
      S_full_cluster  = S_cl,
      S_full_ind      = S_in,
      S_cluster       = NULL,
      S_individual    = NULL,
      RMST_cluster    = R_cl,
      RMST_ind        = R_in,
      trt_prob        = unname(est$trt_prob),
      var_cluster     = if (!is.null(var_out)) var_out$Cluster else NULL,
      var_ind         = if (!is.null(var_out)) var_out$Individual else NULL,
      tau             = tau,

      # NEW: same metadata for RMST fits
      formula_outcome = formula_outcome_final,
      cens_formula    = cens_formula_final,
      cluster_col     = cluster_col,
      n_clusters      = K,
      df_jackknife    = df_jack,
      n_obs           = n_obs,
      n_events        = n_events,
      cluster_trt_counts = cluster_trt_counts,
      call            = match.call(),
      jackknife       = identical(variance, "jackknife")
    )
  }

  class(ans) <- "DRsurvfit"
  ans
}

# ---- S3: summary -------------------------------------------------------------

#' Summary method for DRsurvfit objects
#'
#' For SPCE: prints S1, S0, and (S1-S0) with t-based confidence intervals at
#' user-specified times (or defaults). df for t-intervals is (#clusters - 1).
#' For RMST: prints R1, R0, and (R1-R0) with t-based CIs at requested tau.
#'
#' @param object A DRsurvfit object.
#' @param digits Number of digits to print.
#' @param times Optional numeric times for SPCE summary (ignored for RMST).
#' @param alpha Confidence level = 1 - alpha for the t-intervals. Default 0.05.
#' @param ... Unused.
#' @export
summary.DRsurvfit <- function(object, digits = 4, times = NULL, alpha = 0.05, ...) {
  # Header
  cat(sprintf("DRsurvfit: method = %s, estimand = %s\n", object$method, object$estimand))
  cat("Treatment probs (p0, p1): ",
      paste(signif(object$trt_prob, digits), collapse = ", "), "\n", sep = "")
  if (!is.null(object$formula_outcome)) {
    cat("Outcome model:   ", deparse(object$formula_outcome), "\n", sep = "")
  } else if (!is.null(object$formula)) {
    cat("Outcome model:   ", deparse(object$formula), "\n", sep = "")
  }
  if (!is.null(object$cens_formula)) {
    cat("Censoring model: ", deparse(object$cens_formula), "\n", sep = "")
  }
  if (!is.null(object$cluster_col)) {
    cat("Cluster id col:  ", object$cluster_col, "\n", sep = "")
  }
  if (!is.null(object$n_clusters)) {
    cat("Clusters (M):    ", object$n_clusters, "\n", sep = "")
  }
  if (!is.null(object$n_obs)) {
    cat("Obs (N):         ", object$n_obs, "\n", sep = "")
  }
  cat("\n")

  # formatter: x (lcl, ucl)
  fmt_ci <- function(est, se, tcrit) {
    if (anyNA(se) || is.na(tcrit)) {
      return(sprintf("%s", formatC(est, digits = digits, format = "fg")))
    }
    lcl <- est - tcrit * se
    ucl <- est + tcrit * se
    sprintf("%s (%s, %s)",
            formatC(est, digits = digits, format = "fg"),
            formatC(lcl, digits = digits, format = "fg"),
            formatC(ucl, digits = digits, format = "fg"))
  }

  # prints a block for SPCE/RMST with CIs if Vmat available
  print_block <- function(Mat, Vmat, rown, label, colnames_default, alpha, tcrit) {
    if (is.null(Mat)) return(invisible(NULL))
    cat(label, "\n")
    if (is.null(colnames(Mat)) || any(colnames(Mat) == "")) {
      colnames(Mat) <- colnames_default[seq_len(ncol(Mat))]
    }
    if (!is.null(Vmat)) {
      # Vmat columns: Var(col1), Var(col2), Var(col3), cov(S1,S0)  (we use 1:3)
      se1   <- sqrt(pmax(0, Vmat[, 1L, drop = TRUE]))
      se0   <- sqrt(pmax(0, Vmat[, 2L, drop = TRUE]))
      sedif <- sqrt(pmax(0, Vmat[, 3L, drop = TRUE]))
      c1 <- fmt_ci(Mat[, 1L, drop = TRUE], se1,   tcrit)
      c2 <- fmt_ci(Mat[, 2L, drop = TRUE], se0,   tcrit)
      c3 <- fmt_ci(Mat[, 3L, drop = TRUE], sedif, tcrit)
    } else {
      c1 <- formatC(Mat[, 1L, drop = TRUE], digits = digits, format = "fg")
      c2 <- formatC(Mat[, 2L, drop = TRUE], digits = digits, format = "fg")
      c3 <- formatC(Mat[, 3L, drop = TRUE], digits = digits, format = "fg")
    }
    out <- cbind(
      `Col1 (LCL, UCL)` = c1,
      `Col0 (LCL, UCL)` = c2,
      `ColDiff (LCL, UCL)` = c3
    )
    # Relabel column headers to the desired names
    colnames(out) <- c(
      paste0(colnames_default[1], " (LCL, UCL)"),
      paste0(colnames_default[2], " (LCL, UCL)"),
      paste0(colnames_default[3], " (LCL, UCL)")
    )
    rownames(out) <- rown
    print(noquote(out))
    if (!is.null(Vmat) && !is.na(tcrit)) {
      df <- attr(tcrit, "df")
      if (!is.null(df)) cat(sprintf("  t-intervals with df = %d, alpha = %.3f\n\n", df, alpha))
      else               cat(sprintf("  t-intervals (alpha = %.3f)\n\n", alpha))
    } else if (is.null(Vmat)) {
      cat("  (jackknife variances not available; showing point estimates only)\n\n")
    } else {
      cat("\n")
    }
    invisible(NULL)
  }

  # t critical and df
  K  <- object$n_clusters %||% NA_integer_
  df <- object$df_jackknife %||% if (is.finite(K)) max(1L, K - 1L) else NA_integer_
  tcrit <- if (is.finite(df)) stats::qt(1 - alpha/2, df = df) else NA_real_
  if (is.finite(df)) attr(tcrit, "df") <- df

  if (object$estimand == "SPCE") {
    et <- object$event_time
    # which times to show
    show_times <- times
    if (is.null(show_times) && !is.null(object$spce_summary_times))
      show_times <- object$spce_summary_times
    if (is.null(show_times)) {
      qs <- stats::quantile(et, probs = c(0.10, 0.25, 0.50, 0.75), type = 1)
      show_times <- as.numeric(qs)
    }
    # left-nearest indices on the shared grid
    idx  <- pmax(1L, findInterval(show_times, et))
    rown <- paste0("t=", formatC(et[idx], digits = digits, format = "fg"))

    # slice to selected times
    Scl  <- if (!is.null(object$S_cluster))    object$S_cluster[idx, , drop = FALSE]    else NULL
    Sind <- if (!is.null(object$S_individual)) object$S_individual[idx, , drop = FALSE] else NULL
    Vcl  <- if (!is.null(object$var_cluster))  object$var_cluster[idx, , drop = FALSE]  else NULL
    Vind <- if (!is.null(object$var_ind))      object$var_ind[idx, , drop = FALSE]      else NULL

    # print SPCE blocks
    print_block(Scl,  Vcl,  rown, "Cluster-level SPCE:",
                c("S1","S0","S1-S0"), alpha, tcrit)
    print_block(Sind, Vind, rown, "Individual-level SPCE:",
                c("S1","S0","S1-S0"), alpha, tcrit)

  } else { # RMST
    # row labels
    rown <- NULL
    if (!is.null(object$landmarks)) {
      rown <- paste0("tau=", formatC(as.numeric(object$landmarks), digits = digits, format = "fg"))
    } else if (!is.null(object$tau)) {
      rown <- paste0("tau=", formatC(as.numeric(object$tau), digits = digits, format = "fg"))
    }

    # pull matrices
    Rcl  <- object$RMST_cluster
    Rind <- object$RMST_ind
    Vcl  <- if (!is.null(object$var_cluster)) object$var_cluster else NULL
    Vind <- if (!is.null(object$var_ind))     object$var_ind     else NULL

    # Ensure expected column names
    if (!is.null(Rcl)  && (is.null(colnames(Rcl))  || any(colnames(Rcl)  == ""))) colnames(Rcl)  <- c("R1","R0","R1-R0")
    if (!is.null(Rind) && (is.null(colnames(Rind)) || any(colnames(Rind) == ""))) colnames(Rind) <- c("R1","R0","R1-R0")

    # print RMST blocks with t-based CIs
    print_block(Rcl,  Vcl,  rown, "Cluster-level RMST:",
                c("R1","R0","R1-R0"), alpha, tcrit)
    print_block(Rind, Vind, rown, "Individual-level RMST:",
                c("R1","R0","R1-R0"), alpha, tcrit)
  }

  invisible(object)
}


# ---- S3: plot ---------------------------------------------------------------

#' Plot method for DRsurvfit objects
#'
#' @description
#' Plots survival curves (treatment \code{S1} solid, control \code{S0} dashed)
#' for either the cluster- or individual-level estimand. If \code{tau} is present
#' (RMST fit), vertical lines are drawn at \code{tau}.
#'
#' @param x A \code{DRsurvfit} object.
#' @param level \code{"cluster"} or \code{"individual"}.
#' @param ... Unused.
#'
#' @return The input object \code{x} (invisibly).
#'
#' @examples
#' \dontrun{
#' data(dat)
#' fit_spce <- DRsurvfit(
#'   dat, survival::Surv(time, event) ~ W1 + W2 + Z1 + Z2 + cluster(M),
#'   intv = "A", method = "marginal", estimand = "SPCE"
#' )
#' plot(fit_spce, level = "cluster")
#'
#' fit_rmst <- DRsurvfit(
#'   dat, survival::Surv(time, event) ~ W1 + W2 + Z1 + Z2 + cluster(M),
#'   intv = "A", method = "frailty", estimand = "RMST", tau = c(0.5,1,2)
#' )
#' plot(fit_rmst, level = "individual")
#' }
#' @export
plot.DRsurvfit <- function(x, level = c("cluster", "individual"), ...) {
  level <- match.arg(level)
  mat <- if (level == "cluster") x$S_full_cluster else x$S_full_ind
  if (is.null(mat)) {
    warning("No survival curves available to plot.", call. = FALSE)
    return(invisible(x))
  }
  yy1 <- mat[, 1L]; yy0 <- mat[, 2L]
  plot(x$event_time, yy1, type = "l", xlab = "time", ylab = "Survival",
       main = sprintf("%s-level survival curves (%s)", level, x$method))
  lines(x$event_time, yy0, lty = 2)
  legend("topright", legend = c("Treatment (S1)","Control (S0)"), lty = c(1,2), bty = "n")

  # draw tau if present (RMST)
  tau <- x$tau %||% NULL
  if (!is.null(tau)) abline(v = tau, lty = 3)

  invisible(x)
}
