#' Doubly-robust multi-state survival estimation for CRTs (SPCE or RMT-IF)
#'
#' @description
#' Fits doubly-robust estimators for cluster-randomized trials with right-censored,
#' possibly multi-state outcomes.
#'
#' The outcome is specified as \code{Surv(time, status)}, where
#' \code{status ∈ {0,1,2,...,S_max}} and \code{status = 0} denotes censoring.
#' Values \code{1,2,...,S_max} are ordered states, with the largest state typically
#' representing an absorbing "worst" outcome (e.g., death).
#'
#' The function supports two estimands:
#' \itemize{
#'   \item \strong{SPCE}: stage-specific survival probabilities
#'     \eqn{S_s(t)} for each state \eqn{s=1,\dots,S_{\max}} at all event times.
#'   \item \strong{RMTIF}: a generalized win-based restricted mean time estimand
#'     constructed from the multi-state survival surfaces. When \code{status} is
#'     binary (\code{0/1}), this reduces to an RMST-type quantity (no explicit
#'     RMST target times are required; the estimand is evaluated on the full
#'     event-time grid).
#' }
#'
#' Jackknife variance is computed via leave-one-cluster-out re-fitting on a shared
#' event-time grid:
#' \itemize{
#'   \item For \code{estimand = "SPCE"}: variances of \eqn{S_{1}(t)}, \eqn{S_{0}(t)},
#'         and \eqn{S_{1}(t) - S_{0}(t)} at each time and state.
#'   \item For \code{estimand = "RMTIF"}: variances and covariance of
#'         \eqn{R_{1}(\tau)}, \eqn{R_{0}(\tau)}, and \eqn{R_{1}(\tau) - R_{0}(\tau)}
#'         at each event time \eqn{\tau}.
#' }
#'
#' The returned object includes metadata needed for summaries and plotting:
#' final fitted outcome/censoring formulas, the cluster id column, number of clusters,
#' degrees of freedom for jackknife t-intervals (= M - 1), sample sizes, and the
#' cluster-level treatment split.
#'
#' @param data A \code{data.frame}.
#' @param formula Outcome model: e.g.,
#'   \code{Surv(time, status) ~ W1 + W2 + Z1 + Z2 + cluster(M)}.
#'   The left-hand side must be \code{Surv(time, status)} with
#'   \code{status ∈ {0,1,2,...}} and \code{0} indicating censoring.
#'
#'   The right-hand side \emph{must} include a \code{cluster(<id>)} term specifying
#'   the cluster id for CRTs. All other covariates may be individual- or cluster-level.
#'
#' @param cens_formula Optional censoring model. If \code{NULL}, the censoring model
#'   is built automatically from the outcome formula by:
#'   \itemize{
#'     \item reusing the RHS (excluding \code{cluster()});
#'     \item using LHS \code{Surv(time, event == 0)};
#'     \item appending \code{+ cluster(id)} when a cluster id is present.
#'   }
#'   If supplied, \code{cens_formula} is used as-is for all stage-specific fits,
#'   but the DR estimating equations still use the stage-specific \code{event}
#'   indicator as described above.
#'
#' @param intv Character: name of the cluster-level treatment column (0/1),
#'   constant within cluster.
#'
#' @param method \code{"marginal"} or \code{"frailty"}.
#'   \itemize{
#'     \item \code{"marginal"}: fits \code{survival::coxph} models with
#'           \code{cluster(<id>)} robust variance.
#'     \item \code{"frailty"}: fits \code{frailtyEM::emfrail} gamma-frailty models
#'           for outcome and censoring.
#'   }
#'
#' @param estimand \code{"SPCE"} or \code{"RMTIF"}.
#'   \itemize{
#'     \item \code{"SPCE"}: returns stage-specific survival arrays
#'           \code{S_stage_cluster} and \code{S_stage_ind} with dimensions
#'           \code{[time × 2 × S_max]}.
#'     \item \code{"RMTIF"}: returns the generalized win-based RMST-type estimand
#'           at each event time, along with stage-wise contributions. For a binary
#'           status, this reduces to an RMST-type quantity.
#'   }
#'
#' @param trt_prob Optional length-2 numeric vector \code{(p0, p1)} giving the
#'   cluster-level treatment probabilities for arms 0 and 1. If \code{NULL}, they
#'   are computed as the empirical proportion of first-row treatment assignments
#'   per cluster.
#'
#' @param variance \code{"none"} or \code{"jackknife"} for variance estimation.
#'
#' @param fit_controls Optional \code{frailtyEM::emfrail_control()} list, used only
#'   when \code{method = "frailty"}. If \code{NULL}, default fast-fitting controls
#'   are used (no standard errors from the frailtyEM fits are required here).
#'
#' @param verbose Logical; currently unused but kept for future verbosity options.
#'
#' @return An object of class \code{"DRsurvfit"} with fields depending on
#'   \code{estimand}:
#'
#'   \describe{
#'     \item{Common:}{
#'       \itemize{
#'         \item \code{method}: fitted method (\code{"marginal"} or \code{"frailty"}).
#'         \item \code{estimand}: requested estimand (\code{"SPCE"} or \code{"RMTIF"}).
#'         \item \code{trt_prob}: numeric vector \code{c(p0, p1)}.
#'         \item \code{event_time}: time grid:
#'           \itemize{
#'             \item SPCE: all event times including \code{0}.
#'             \item RMTIF: positive event times \eqn{\tau} at which the RMT-IF is evaluated.
#'           }
#'         \item \code{max_state}: maximum observed non-zero status.
#'         \item \code{cluster_col}: name of the cluster id column.
#'         \item \code{n_clusters}: number of clusters (\eqn{M}).
#'         \item \code{df_jackknife}: jackknife degrees of freedom (\eqn{M - 1}).
#'         \item \code{n_obs}: total number of observations.
#'         \item \code{n_events}: total number of non-censoring observations
#'               (\code{status != 0}).
#'         \item \code{cluster_trt_counts}: counts of treated and control clusters
#'               \code{c(n_trt0, n_trt1)} based on first row per cluster.
#'         \item \code{formula_outcome}: fully reconstructed outcome formula.
#'         \item \code{cens_formula}: final censoring formula used.
#'         \item \code{call}: the matched call.
#'         \item \code{jackknife}: logical indicating whether jackknife variance
#'               was computed.
#'       }
#'     }
#'
#'     \item{If \code{estimand = "SPCE"}:}{
#'       \itemize{
#'         \item \code{S_stage_cluster}: 3D array \code{[time × 2 × S_max]} with
#'               stage-specific cluster-level survival:
#'               \code{S_stage_cluster[ , 1, s]} = \eqn{S_1^{(s)}(t)},
#'               \code{S_stage_cluster[ , 2, s]} = \eqn{S_0^{(s)}(t)}.
#'         \item \code{S_stage_ind}: analogous individual-level survival array.
#'         \item \code{var_stage_cluster}: jackknife variances for
#'               \eqn{S_1^{(s)}(t)}, \eqn{S_0^{(s)}(t)}, and
#'               \eqn{S_1^{(s)}(t) - S_0^{(s)}(t)} as a 3D array
#'               \code{[time × 3 × S_max]} with dimension names
#'               \code{comp = c("Var(S1)","Var(S0)","Var(S1-S0)")}, when
#'               \code{variance = "jackknife"}; otherwise \code{NULL}.
#'         \item \code{var_stage_ind}: analogous individual-level variance array.
#'       }
#'     }
#'
#'     \item{If \code{estimand = "RMTIF"}:}{
#'       \itemize{
#'         \item \code{RMTIF_cluster}: matrix \code{[time × 3]} with columns
#'               \code{c("R1","R0","R1-R0")} giving the cluster-level RMT-IF curves
#'               at each event time \eqn{\tau}.
#'         \item \code{RMTIF_ind}: analogous individual-level RMT-IF matrix.
#'         \item \code{stagewise_cluster}: list of length \code{length(event_time)};
#'               each element is a \code{3 × (S_max+1)} matrix of stage-wise
#'               contributions with rows
#'               \code{c("s1qs0qp1","s0qs1qp1","diff")} and columns
#'               \code{c("stage_1",...,"stage_Smax","sum")}.
#'         \item \code{stagewise_ind}: analogous individual-level list.
#'         \item \code{var_rmtif_cluster}: jackknife variance/covariance matrix
#'               \code{[time × 4]} with columns
#'               \code{c("Var(R1)","Var(R0)","Var(R1-R0)","Cov(R1,R0)")},
#'               when \code{variance = "jackknife"}; otherwise \code{NULL}.
#'         \item \code{var_rmtif_ind}: analogous individual-level matrix.
#'         \item \code{S_stage_cluster}, \code{S_stage_ind}: the underlying
#'               stage-specific survival arrays are also returned for convenience.
#'       }
#'     }
#'   }
#'
#' @examples
#' \dontrun{
#' data(dat)
#'
#' ## Multi-state SPCE (stage-specific survival)
#' fit_spce <- DRsurvfit(
#'   datm,
#'   survival::Surv(time, status) ~ W1 + W2 + Z1 + Z2 + cluster(cluster),
#'   intv    = "trt",
#'   method  = "marginal",
#'   estimand = "SPCE",
#'   variance = "jackknife"
#' )
#'
#' ## Multi-state RMT-IF (binary reduces to RMST-type)
#' fit_rmtif <- DRsurvfit(
#'   datm,
#'   survival::Surv(time, status) ~ W1 + W2 + Z1 + Z2 + cluster(cluster),
#'   intv    = "trt",
#'   method  = "frailty",
#'   estimand = "RMTIF",
#'   variance = "jackknife"
#' )
#' }
#' @export
DRsurvfit <- function(data,
                      formula,
                      cens_formula = NULL,
                      intv,
                      method   = c("marginal", "frailty"),
                      estimand = c("SPCE", "RMTIF"),
                      trt_prob = NULL,
                      variance = c("none", "jackknife"),
                      fit_controls = NULL,
                      verbose = FALSE) {

  method   <- match.arg(method)
  estimand <- match.arg(estimand)
  variance <- match.arg(variance)

  if (!is.data.frame(data))
    stop("'data' must be a data.frame.", call. = FALSE)

  ## --- fit core on full grid -----------------------------------------------
  est <- .DR_est_core(
    data         = data,
    formula      = formula,
    cens_formula = cens_formula,
    intv         = intv,
    method       = method,
    estimand     = estimand,
    trt_prob     = trt_prob,
    fit_controls = fit_controls,
    e_time       = NULL
  )

  et        <- est$event_time
  S_stage_c <- est$S_stage_cluster   # [time × 2 × max_state] or NULL
  S_stage_i <- est$S_stage_ind
  max_s     <- est$max_state
  cluster_col <- est$cluster_col

  ## --- reconstruct final formulas and cluster info -------------------------
  environment(formula) <- .surv_env(parent = environment(formula) %||% parent.frame())
  if (!is.null(cens_formula))
    environment(cens_formula) <- .surv_env(parent = environment(cens_formula) %||% parent.frame())

  nm   <- .surv_lhs(formula)
  info <- .extract_cluster(formula)
  clvar <- info$cluster
  rhs   <- info$rhs_wo_cluster

  if (is.null(clvar) || !nzchar(clvar))
    stop("Outcome formula must include cluster(<id>) for CRT.", call. = FALSE)

  rhs_surv <- rhs
  if (nzchar(rhs_surv))
    rhs_surv <- paste(rhs_surv, sprintf("+ cluster(%s)", clvar))
  else
    rhs_surv <- sprintf("cluster(%s)", clvar)

  if (is.null(cens_formula)) {
    cf <- .build_cens_formula_from(nm, rhs, clvar)
    cens_formula_final <- as.formula(cf)
  } else {
    cens_formula_final <- cens_formula
  }

  formula_outcome_final <- as.formula(
    sprintf("Surv(%s,%s) ~ %s", nm$time, nm$status, rhs_surv)
  )

  n_obs <- nrow(data)
  n_events <- if (nm$status %in% names(data)) {
    sum(as.integer(data[[nm$status]] != 0L), na.rm = TRUE)
  } else {
    NA_integer_
  }

  K       <- length(unique(data[[cluster_col]]))
  df_jack <- max(1L, K - 1L)

  ## cluster-level treatment split (first row per cluster)
  cl_ids <- unique(data[[cluster_col]])
  A_first <- vapply(cl_ids, function(g)
    data[[intv]][which(data[[cluster_col]] == g)[1]],
    numeric(1L))
  cluster_trt_counts <- c(
    n_trt0 = sum(A_first == 0, na.rm = TRUE),
    n_trt1 = sum(A_first == 1, na.rm = TRUE)
  )

  ## --- jackknife variance --------------------------------------------------
  var_out <- if (identical(variance, "jackknife")) {
    .DR_var_jackknife(
      data         = data,
      formula      = formula,
      cens_formula = cens_formula,
      intv         = intv,
      method       = method,
      estimand     = estimand,
      trt_prob     = trt_prob,
      fit_controls = fit_controls,
      e_time_full  = et[et > 0]
    )
  } else {
    NULL
  }

  ## --- assemble return object ---------------------------------------------
  if (estimand == "SPCE") {
    var_stage_cluster <- if (!is.null(var_out)) var_out$Cluster$var_stage else NULL
    var_stage_ind     <- if (!is.null(var_out)) var_out$Individual$var_stage else NULL

    ans <- list(
      method   = method,
      estimand = "SPCE",
      event_time        = et,            # includes 0
      max_state         = max_s,
      S_stage_cluster   = S_stage_c,
      S_stage_ind       = S_stage_i,
      RMTIF_cluster     = NULL,
      RMTIF_ind         = NULL,
      stagewise_cluster = NULL,
      stagewise_ind     = NULL,
      trt_prob          = unname(est$trt_prob),
      var_stage_cluster = var_stage_cluster,
      var_stage_ind     = var_stage_ind,

      ## metadata
      formula_outcome   = formula_outcome_final,
      cens_formula      = cens_formula_final,
      cluster_col       = cluster_col,
      n_clusters        = K,
      df_jackknife      = df_jack,
      n_obs             = n_obs,
      n_events          = n_events,
      cluster_trt_counts = cluster_trt_counts,
      call              = match.call(),
      jackknife         = identical(variance, "jackknife")
    )

  } else {  # RMTIF
    var_rmtif_cluster <- if (!is.null(var_out)) var_out$Cluster$var_R else NULL
    var_rmtif_ind     <- if (!is.null(var_out)) var_out$Individual$var_R else NULL

    ans <- list(
      method   = method,
      estimand = "RMTIF",
      event_time        = et,                # τ grid > 0
      max_state         = max_s,
      S_stage_cluster   = S_stage_c,
      S_stage_ind       = S_stage_i,
      RMTIF_cluster     = est$RMTIF_cluster,
      RMTIF_ind         = est$RMTIF_ind,
      stagewise_cluster = est$stagewise_cluster,
      stagewise_ind     = est$stagewise_ind,
      trt_prob          = unname(est$trt_prob),
      var_rmtif_cluster = var_rmtif_cluster,
      var_rmtif_ind     = var_rmtif_ind,

      ## metadata
      formula_outcome   = formula_outcome_final,
      cens_formula      = cens_formula_final,
      cluster_col       = cluster_col,
      n_clusters        = K,
      df_jackknife      = df_jack,
      n_obs             = n_obs,
      n_events          = n_events,
      cluster_trt_counts = cluster_trt_counts,
      call              = match.call(),
      jackknife         = identical(variance, "jackknife")
    )
  }

   class(ans) <- "DRsurvfit"
  ans
}


#' Summary method for DRsurvfit objects (multi-state SPCE / RMT-IF)
#'
#' @description
#' Produces tabular summaries for multi-state doubly-robust estimators:
#' \itemize{
#'   \item For \code{estimand = "SPCE"}: stage-specific survival probabilities
#'         \eqn{S_s(t)} at selected times \eqn{\tau}, with optional jackknife
#'         t-based confidence intervals.
#'   \item For \code{estimand = "RMTIF"}: RMT-IF curves \eqn{R_1(\tau)},
#'         \eqn{R_0(\tau)}, and \eqn{R_1(\tau) - R_0(\tau)} at the same set
#'         of \eqn{\tau}, again with optional jackknife t-based intervals.
#' }
#'
#' The same argument \code{tau} is used for both estimands. If \code{tau} is
#' \code{NULL}, the function uses the 25\%, 50\%, and 75\% quantiles of the
#' event-time grid (excluding time 0 if present).
#'
#' @param object A \code{DRsurvfit} object.
#' @param level Character: \code{"cluster"} or \code{"individual"} level summary.
#' @param tau Optional numeric vector of times at which to summarize both
#'   \code{SPCE} and \code{RMTIF}. If \code{NULL}, the 25\%, 50\%, and 75\%
#'   quantiles of the event-time grid are used.
#' @param states Optional integer vector of states to summarize for
#'   \code{estimand = "SPCE"}. Defaults to all states \code{1:object$max_state}.
#' @param digits Number of digits to print for estimates and confidence limits.
#' @param alpha Nominal type I error for the intervals; coverage is
#'   \code{1 - alpha}. Default is \code{0.05}, giving 95\% confidence intervals.
#' @param ... Unused.
#'
#' @return The input object \code{object}, invisibly.
#' @export
summary.DRsurvfit <- function(object,
                              level  = c("cluster", "individual"),
                              tau    = NULL,
                              states = NULL,
                              digits = 4,
                              alpha  = 0.05,
                              ...) {

  level <- match.arg(level)

  ## --- header -------------------------------------------------------------- ##
  cat(sprintf("DRsurvfit: method = %s, estimand = %s\n",
              object$method, object$estimand))
  if (!is.null(object$trt_prob)) {
    cat("Treatment probs (p0, p1): ",
        paste(signif(object$trt_prob, digits), collapse = ", "),
        "\n", sep = "")
  }
  if (!is.null(object$formula_outcome)) {
    cat("Outcome model:   ", deparse(object$formula_outcome), "\n", sep = "")
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
  if (!is.null(object$n_events)) {
    cat("Events (status != 0): ", object$n_events, "\n", sep = "")
  }
  cat("\n")

  ## --- helper: format estimate + CI ---------------------------------------- ##
  fmt_ci <- function(est, se, tcrit) {
    if (is.null(se) || anyNA(se) || is.na(tcrit)) {
      return(formatC(est, digits = digits, format = "fg"))
    }
    lcl <- est - tcrit * se
    ucl <- est + tcrit * se
    sprintf("%s (%s, %s)",
            formatC(est, digits = digits, format = "fg"),
            formatC(lcl, digits = digits, format = "fg"),
            formatC(ucl, digits = digits, format = "fg"))
  }

  ## --- t critical and df --------------------------------------------------- ##
  K  <- object$n_clusters %||% NA_integer_
  df <- object$df_jackknife %||%
    if (is.finite(K)) max(1L, K - 1L) else NA_integer_
  tcrit <- if (is.finite(df)) stats::qt(1 - alpha / 2, df = df) else NA_real_
  if (is.finite(df)) attr(tcrit, "df") <- df

  et <- object$event_time

  ## --- choose tau if NULL -------------------------------------------------- ##
  if (is.null(tau)) {
    et_pos <- et[et > 0]
    if (length(et_pos)) {
      qs <- stats::quantile(et_pos, probs = c(0.25, 0.50, 0.75), type = 1)
      tau <- as.numeric(qs)
    } else {
      tau <- unique(et)
    }
  }
  tau <- sort(unique(as.numeric(tau)))
  idx <- pmax(1L, findInterval(tau, et))
  row_lab <- paste0("t=", formatC(et[idx], digits = digits, format = "fg"))

  ## ======================================================================== ##
  ## SPCE: stage-specific survival S_s(t)
  ## ======================================================================== ##
  if (identical(object$estimand, "SPCE")) {
    S_arr <- if (level == "cluster") object$S_stage_cluster else object$S_stage_ind
    V_arr <- if (level == "cluster") object$var_stage_cluster else object$var_stage_ind

    if (is.null(S_arr) || length(dim(S_arr)) != 3L) {
      warning("SPCE summaries unavailable: missing stage-specific survival arrays.", call. = FALSE)
      return(invisible(object))
    }

    max_s <- object$max_state %||% dim(S_arr)[3]
    if (is.null(max_s) || is.na(max_s)) max_s <- dim(S_arr)[3]

    if (is.null(states)) {
      states <- seq_len(max_s)
    } else {
      states <- intersect(states, seq_len(max_s))
      if (!length(states)) {
        stop("No valid states to summarize.", call. = FALSE)
      }
    }

    for (s in states) {
      cat(sprintf("Stage-specific SPCE: state %d (%s-level)\n", s, level))

      ## slice survival: [length(tau) × 2]
      S_slice <- S_arr[idx, , s, drop = FALSE]
      S1 <- S_slice[, 1L, drop = TRUE]
      S0 <- S_slice[, 2L, drop = TRUE]
      Sd <- S1 - S0

      if (!is.null(V_arr)) {
        # V_arr: [time × 3 × state]
        V_slice <- V_arr[idx, , s, drop = FALSE]
        se1   <- sqrt(pmax(0, V_slice[, 1L, drop = TRUE]))
        se0   <- sqrt(pmax(0, V_slice[, 2L, drop = TRUE]))
        sedif <- sqrt(pmax(0, V_slice[, 3L, drop = TRUE]))
        c1 <- fmt_ci(S1, se1, tcrit)
        c2 <- fmt_ci(S0, se0, tcrit)
        c3 <- fmt_ci(Sd, sedif, tcrit)
      } else {
        c1 <- formatC(S1, digits = digits, format = "fg")
        c2 <- formatC(S0, digits = digits, format = "fg")
        c3 <- formatC(Sd, digits = digits, format = "fg")
      }

      out <- cbind(
        `S1 (LCL, UCL)`     = c1,
        `S0 (LCL, UCL)`     = c2,
        `S1-S0 (LCL, UCL)`  = c3
      )
      rownames(out) <- row_lab
      print(noquote(out))
      if (!is.null(V_arr) && !is.na(tcrit)) {
        df_here <- attr(tcrit, "df")
        if (!is.null(df_here)) {
          cat(sprintf("  t-intervals with df = %d, alpha = %.3f\n\n", df_here, alpha))
        } else {
          cat(sprintf("  t-intervals (alpha = %.3f)\n\n", alpha))
        }
      } else if (is.null(V_arr)) {
        cat("  (jackknife variances not available; showing point estimates only)\n\n")
      } else {
        cat("\n")
      }
    }

    return(invisible(object))
  }

  ## ======================================================================== ##
  ## RMTIF: win-based RMST-type curves over event time
  ## ======================================================================== ##
  if (identical(object$estimand, "RMTIF")) {
    Rmat <- if (level == "cluster") object$RMTIF_cluster else object$RMTIF_ind
    Vmat <- if (level == "cluster") object$var_rmtif_cluster else object$var_rmtif_ind

    if (is.null(Rmat) || !is.matrix(Rmat) || ncol(Rmat) < 2L) {
      warning("RMTIF summaries unavailable: missing RMTIF matrices.", call. = FALSE)
      return(invisible(object))
    }

    R1 <- Rmat[idx, 1L, drop = TRUE]
    R0 <- Rmat[idx, 2L, drop = TRUE]
    Rd <- Rmat[idx, 3L, drop = TRUE]

    if (!is.null(Vmat)) {
      se1   <- sqrt(pmax(0, Vmat[idx, 1L, drop = TRUE]))
      se0   <- sqrt(pmax(0, Vmat[idx, 2L, drop = TRUE]))
      sedif <- sqrt(pmax(0, Vmat[idx, 3L, drop = TRUE]))
      c1 <- fmt_ci(R1, se1, tcrit)
      c2 <- fmt_ci(R0, se0, tcrit)
      c3 <- fmt_ci(Rd, sedif, tcrit)
    } else {
      c1 <- formatC(R1, digits = digits, format = "fg")
      c2 <- formatC(R0, digits = digits, format = "fg")
      c3 <- formatC(Rd, digits = digits, format = "fg")
    }

    cat(sprintf("RMT-IF summary (%s-level)\n", level))
    out <- cbind(
      `R1 (LCL, UCL)`       = c1,
      `R0 (LCL, UCL)`       = c2,
      `R1-R0 (LCL, UCL)`    = c3
    )
    rownames(out) <- row_lab
    print(noquote(out))
    if (!is.null(Vmat) && !is.na(tcrit)) {
      df_here <- attr(tcrit, "df")
      if (!is.null(df_here)) {
        cat(sprintf("  t-intervals with df = %d, alpha = %.3f\n\n", df_here, alpha))
      } else {
        cat(sprintf("  t-intervals (alpha = %.3f)\n\n", alpha))
      }
    } else if (is.null(Vmat)) {
      cat("  (jackknife variances not available; showing point estimates only)\n\n")
    } else {
      cat("\n")
    }

    cat("Stage-wise RMT-IF decompositions are available in:\n")
    cat(sprintf("  object$stagewise_%s[[k]] for each time index k.\n",
                if (level == "cluster") "cluster" else "ind"))
    return(invisible(object))
  }

  warning("Unknown estimand in DRsurvfit object.", call. = FALSE)
  invisible(object)
}



#' Plot method for DRsurvfit objects (SPCE / RMT-IF)
#'
#' @description
#' Produces plots for multi-state doubly-robust estimators:
#' \itemize{
#'   \item For \code{estimand = "SPCE"}: for each state \eqn{s}, plots the
#'         difference curve \eqn{S_{1,s}(t) - S_{0,s}(t)} with jackknife
#'         t-based confidence bands over time.
#'   \item For \code{estimand = "RMTIF"}: plots the overall RMT-IF difference
#'         curve \eqn{R_1(t) - R_0(t)} (sum of stage-wise contributions) with
#'         jackknife t-based confidence bands over time.
#' }
#'
#' The argument \code{tau} is a truncation time: if supplied, the plot is
#' restricted to \code{event_time <= tau}. If \code{tau} is \code{NULL}, the
#' full event-time grid is used.
#'
#' @param x A \code{DRsurvfit} object.
#' @param level Character: \code{"cluster"} or \code{"individual"}.
#' @param states Optional integer vector of states to plot when
#'   \code{estimand = "SPCE"}. Defaults to all states
#'   \code{1:object$max_state}.
#' @param tau Optional numeric truncation time. If non-\code{NULL}, only
#'   event times \code{<= max(tau)} are plotted. If \code{NULL}, all event
#'   times are plotted.
#' @param alpha Nominal type I error for the intervals; coverage is
#'   \code{1 - alpha}. Default is \code{0.05}.
#' @param ... Unused; included for S3 consistency.
#'
#' @return The input object \code{x}, invisibly.
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs theme_minimal
plot.DRsurvfit <- function(x,
                           level  = c("cluster", "individual"),
                           states = NULL,
                           tau    = NULL,
                           alpha  = 0.05,
                           ...) {

  level <- match.arg(level)
  et <- x$event_time

  if (is.null(et) || !length(et)) {
    warning("No event_time grid stored in DRsurvfit object; nothing to plot.",
            call. = FALSE)
    return(invisible(x))
  }

  ## ----- truncation by tau ------------------------------------------------- ##
  if (!is.null(tau)) {
    tau_star <- max(as.numeric(tau), na.rm = TRUE)
    idx_time <- which(et <= tau_star)
    if (!length(idx_time)) {
      stop("No event_time <= tau; cannot plot.", call. = FALSE)
    }
  } else {
    idx_time <- seq_along(et)
    tau_star <- max(et, na.rm = TRUE)
  }
  et_plot <- et[idx_time]

  ## ----- t critical using df from object ----------------------------------- ##
  K  <- x$n_clusters %||% NA_integer_
  df <- x$df_jackknife %||%
    if (is.finite(K)) max(1L, K - 1L) else NA_integer_
  tcrit <- if (is.finite(df)) stats::qt(1 - alpha / 2, df = df) else NA_real_

  ## ======================================================================== ##
  ## SPCE: plot only S1 - S0 per state with CI
  ## ======================================================================== ##
  if (identical(x$estimand, "SPCE")) {
    S_arr <- if (level == "cluster") x$S_stage_cluster else x$S_stage_ind
    V_arr <- if (level == "cluster") x$var_stage_cluster else x$var_stage_ind

    if (is.null(S_arr) || length(dim(S_arr)) != 3L) {
      warning("SPCE plot unavailable: missing stage-specific survival arrays.",
              call. = FALSE)
      return(invisible(x))
    }

    ## subset time dimension by idx_time
    S_arr <- S_arr[idx_time, , , drop = FALSE]
    if (!is.null(V_arr)) {
      V_arr <- V_arr[idx_time, , , drop = FALSE]
    }

    max_s <- x$max_state %||% dim(S_arr)[3]
    if (is.null(max_s) || is.na(max_s)) max_s <- dim(S_arr)[3]

    if (is.null(states)) {
      states <- seq_len(max_s)
    } else {
      states <- intersect(states, seq_len(max_s))
      if (!length(states)) {
        stop("No valid states to plot.", call. = FALSE)
      }
    }

    ## build plotting data.frame: one row per time × state
    df_list <- vector("list", length(states))
    for (i in seq_along(states)) {
      s <- states[i]

      S1 <- S_arr[, 1L, s, drop = TRUE]
      S0 <- S_arr[, 2L, s, drop = TRUE]
      diff <- S1 - S0

      if (!is.null(V_arr)) {
        var_diff <- V_arr[, 3L, s, drop = TRUE]  # third col: var(diff)
        se_diff  <- sqrt(pmax(0, var_diff))
        if (is.finite(tcrit)) {
          lcl <- diff - tcrit * se_diff
          ucl <- diff + tcrit * se_diff
        } else {
          lcl <- ucl <- NA_real_
        }
      } else {
        lcl <- ucl <- NA_real_
      }

      df_list[[i]] <- data.frame(
        time  = et_plot,
        diff  = diff,
        lcl   = lcl,
        ucl   = ucl,
        state = factor(s, levels = states,
                       labels = paste0("state ", states))
      )
    }

    df_plot <- do.call(rbind, df_list)

    p <- ggplot2::ggplot(df_plot,
                         ggplot2::aes(x = time, y = diff,
                                      color = state, fill = state)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lcl, ymax = ucl),
        alpha = if (!all(is.na(df_plot$lcl))) 0.20 else 0,
        color = NA
      ) +
      ggplot2::geom_line() +
      ggplot2::labs(
        x     = "time",
        y     = expression(S[1](t) - S[0](t)),
        color = "State",
        fill  = "State",
        title = sprintf("Stage-specific SPCE differences (%s-level)", level),
        subtitle = if (is.finite(df))
          sprintf("Truncated at t <= %.3f; t-based %.1f%% CIs, df = %d",
                  tau_star, 100 * (1 - alpha), df)
        else
          sprintf("Truncated at t <= %.3f; CIs not available", tau_star)
      ) +
      ggplot2::theme_minimal()

    print(p)
    return(invisible(x))
  }

  ## ======================================================================== ##
  ## RMT-IF: plot only R1 - R0 over time with CI
  ## ======================================================================== ##
  if (identical(x$estimand, "RMTIF")) {
    Rmat <- if (level == "cluster") x$RMTIF_cluster else x$RMTIF_ind
    Vmat <- if (level == "cluster") x$var_rmtif_cluster else x$var_rmtif_ind

    if (is.null(Rmat) || !is.matrix(Rmat) || ncol(Rmat) < 2L) {
      warning("RMT-IF plot unavailable: missing RMT-IF matrices.",
              call. = FALSE)
      return(invisible(x))
    }

    ## subset time dimension
    Rmat <- Rmat[idx_time, , drop = FALSE]
    if (!is.null(Vmat)) {
      Vmat <- Vmat[idx_time, , drop = FALSE]
    }

    ## difference curve
    if (ncol(Rmat) >= 3L) {
      diff <- Rmat[, 3L, drop = TRUE]
    } else {
      diff <- Rmat[, 1L, drop = TRUE] - Rmat[, 2L, drop = TRUE]
    }

    if (!is.null(Vmat) && ncol(Vmat) >= 3L) {
      var_diff <- Vmat[, 3L, drop = TRUE]
      se_diff  <- sqrt(pmax(0, var_diff))
      if (is.finite(tcrit)) {
        lcl <- diff - tcrit * se_diff
        ucl <- diff + tcrit * se_diff
      } else {
        lcl <- ucl <- NA_real_
      }
    } else {
      lcl <- ucl <- NA_real_
    }

    df_plot <- data.frame(
      time = et_plot,
      diff = diff,
      lcl  = lcl,
      ucl  = ucl
    )

    p <- ggplot2::ggplot(df_plot,
                         ggplot2::aes(x = time, y = diff)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lcl, ymax = ucl),
        alpha = if (!all(is.na(df_plot$lcl))) 0.20 else 0,
        color = NA
      ) +
      ggplot2::geom_line() +
      ggplot2::labs(
        x = "time",
        y = expression(R[1](t) - R[0](t)),
        title = sprintf("RMT-IF difference (%s-level)", level),
        subtitle = if (is.finite(df))
          sprintf("Truncated at t <= %.3f; t-based %.1f%% CIs, df = %d",
                  tau_star, 100 * (1 - alpha), df)
        else
          sprintf("Truncated at t <= %.3f; CIs not available", tau_star)
      ) +
      ggplot2::theme_minimal()

    print(p)
    return(invisible(x))
  }

  warning("Unknown estimand in DRsurvfit object; nothing plotted.", call. = FALSE)
  invisible(x)
}
