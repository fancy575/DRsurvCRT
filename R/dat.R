#' Example CRT survival dataset
#'
#' `dat` is a toy dataset used to demonstrate \code{DRsurvfit}. It mimics a
#' cluster-randomized trial (CRT) with right-censored survival outcomes.
#'
#' @details
#' - Cluster identifier: \strong{M}
#' - Cluster-level treatment: \strong{A} (binary, constant within each \code{M})
#' - Cluster-level covariates: \strong{W1}, \strong{W2} (constant within each \code{M})
#' - Individual-level covariates: \strong{Z1}, \strong{Z2}
#' - Survival outcome: \strong{time} (follow-up time), \strong{event} (1 = event, 0 = censored)
#'
#' In a valid CRT, the cluster-level covariates \code{W1}, \code{W2} and the cluster treatment
#' \code{A} do not vary within clusters; i.e., each cluster \code{M} has a single value of \code{A}, \code{W1}, \code{W2}.
#'
#' @format A data frame with rows corresponding to individuals and the following columns:
#' \describe{
#'   \item{M}{Cluster ID (integer or factor).}
#'   \item{A}{Cluster-level treatment assignment (0/1), \emph{constant within each} \code{M}.}
#'   \item{W1}{Cluster-level covariate, \emph{constant within each} \code{M}.}
#'   \item{W2}{Cluster-level covariate, \emph{constant within each} \code{M}.}
#'   \item{Z1}{Individual-level covariate.}
#'   \item{Z2}{Individual-level covariate.}
#'   \item{time}{Follow-up time.}
#'   \item{event}{Event indicator (1 = event, 0 = censored).}
#' }
#'
#' @examples
#' # Load the dataset
#' data(dat)
#'
#' # Check cluster-constancy of A, W1, W2 within M (should all be TRUE)
#' all(tapply(dat$A, dat$M, function(x) length(unique(x)) == 1))
#' all(tapply(dat$W1, dat$M, function(x) length(unique(x)) == 1))
#' all(tapply(dat$W2, dat$M, function(x) length(unique(x)) == 1))
#'
#' @usage data(dat)
#' @keywords datasets
"dat"
