#' Simulated single-state CRT survival data
#'
#' @description
#' A simulated cluster-randomized trial (CRT) dataset with a \emph{single}
#' terminal event type and right censoring. Each row corresponds to one
#' individual in a cluster. The data are intended for illustrating the
#' doubly-robust survival estimators in \pkg{DRsurvCRT} for the simple
#' two-state setting (alive vs. terminal event).
#'
#' @format A data frame with \code{n} rows and 8 variables:
#' \describe{
#'   \item{cluster}{Integer cluster identifier.}
#'   \item{trt}{Cluster-level treatment indicator; \code{0} = control,
#'     \code{1} = intervention. Constant within cluster.}
#'   \item{W1}{Cluster-level baseline binary covariate (e.g., center-level
#'     characteristic).}
#'   \item{W2}{Cluster-level baseline continuous covariate.}
#'   \item{Z1}{Individual-level continuous baseline covariate.}
#'   \item{Z2}{Individual-level binary baseline covariate.}
#'   \item{time}{Observed follow-up time (event or censoring time).}
#'   \item{event}{Event indicator; \code{1} = terminal event,
#'     \code{0} = right censored.}
#' }
#'
#' @details
#' The dataset was generated from a cluster-randomized design with
#' covariate-dependent hazards and administrative censoring. It is primarily
#' used to demonstrate calls such as \code{DRsurvfit(..., estimand = "SPCE")}
#' and the corresponding variance and plotting methods in the single-state
#' setting.
#'
#' @examples
#' data(dats)
#'
#' ## quick look
#' head(dats)
#' table(dats$cluster, dats$trt)
#'
#' ## Example use with DRsurvfit (marginal Cox working model)
#' \dontrun{
#' fit_spce <- DRsurvfit(
#'   data    = dats,
#'   formula = survival::Surv(time, event) ~ W1 + W2 + Z1 + Z2 + cluster(cluster),
#'   intv    = "trt",
#'   method  = "marginal",
#'   estimand = "SPCE",
#'   variance = "jackknife"
#' )
#'
#' summary(fit_spce)
#' plot(fit_spce, level = "cluster")
#' }
#'
#' @keywords datasets
"dats"

#' Simulated multi-state CRT survival data (long format)
#'
#' @description
#' A simulated cluster-randomized trial (CRT) dataset with a \emph{multi-state}
#' prioritized endpoint and right censoring, stored in long format. Each
#' individual can experience up to several ordered event types, and contributes
#' multiple rows with increasing \code{time}. The data are designed to
#' illustrate the multi-state SPCE and RMT-IF estimators in \pkg{DRsurvCRT}.
#'
#' @format A data frame with \code{N} rows and 9 variables:
#' \describe{
#'   \item{id}{Individual identifier. Each subject may appear on multiple rows.}
#'   \item{time}{Observed time of the record (e.g., transition time or censoring
#'     time). Times are non-decreasing within an \code{id}.}
#'   \item{event}{Multi-state event indicator. The coding is:
#'     \itemize{
#'       \item \code{0} = still at-risk / no new transition (often the last
#'         row for a censored subject);
#'       \item \code{1}, \code{2}, \code{3} = entry into progressively
#'         more severe event states (e.g., non-fatal event, progression,
#'         death), in increasing clinical priority.
#'     }
#'     Not every subject experiences all states; some may be censored earlier.}
#'   \item{cluster}{Integer cluster identifier.}
#'   \item{trt}{Cluster-level treatment indicator; \code{0} = control,
#'     \code{1} = intervention. Constant within cluster.}
#'   \item{W1}{Cluster-level baseline binary covariate.}
#'   \item{W2}{Cluster-level baseline continuous covariate.}
#'   \item{Z1}{Individual-level continuous baseline covariate.}
#'   \item{Z2}{Individual-level binary baseline covariate.}
#' }
#'
#' @details
#' For each subject \code{id}, rows are ordered by \code{time}. If an individual
#' transitions to a new state (e.g., from state 0 to 1, or from 1 to 2, etc.),
#' the corresponding transition time is recorded with \code{event > 0}. Later
#' rows for the same subject may reflect further transitions or censoring.
#'
#' This structure is suitable for constructing state-specific transformed
#' survival functions \eqn{S_{a,s}(t)} (for treatment arm \eqn{a} and
#' state \eqn{s}) and for computing stage-wise RMT-IF contributions for
#' prioritized composite endpoints. It is the canonical example dataset used
#' by the multi-state core and jackknife routines in \pkg{DRsurvCRT}.
#'
#' @examples
#' data(datm)
#'
#' ## quick look
#' head(datm)
#' table(datm$cluster, datm$trt)
#' table(datm$event)
#'
#' ## Example: fit multi-state SPCE / RMT-IF (schematic)
#' \dontrun{
#' fit_ms <- DRsurvfit(
#'   data    = datm,
#'   formula = survival::Surv(time, event) ~ W1 + W2 + Z1 + Z2 + cluster(cluster),
#'   intv    = "trt",
#'   method  = "frailty",
#'   estimand = "SPCE",           # or "RMTIF"
#'   variance = "jackknife"
#' )
#'
#' summary(fit_ms, level = "cluster")
#' plot(fit_ms, level = "cluster")
#' }
#'
#' @keywords datasets
"datm"
