#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// ---- helpers ---------------------------------------------------------------

static inline void clamp01_monotone(mat& S){ // cols: 0=S1,1=S0,2=diff
  const uword t = S.n_rows;
  // clamp to [0,1]
  S.cols(0,1).for_each([](double& v){ if(v>1.0) v=1.0; else if(v<0.0) v=0.0; });
  // enforce non-increasing
  for(uword k=1;k<t;++k){
    if(S(k,0) > S(k-1,0)) S(k,0) = S(k-1,0);
    if(S(k,1) > S(k-1,1)) S(k,1) = S(k-1,1);
  }
  // recompute diff
  S.col(2) = S.col(0) - S.col(1);
}

static inline uword lower_idx(const vec& tgrid, double tt){
  // first index j with tgrid[j] >= tt
  uword lo=0, hi=tgrid.n_elem;
  while(lo<hi){ uword md=(lo+hi)>>1; if(tgrid[md] < tt) lo=md+1; else hi=md; }
  if(lo==tgrid.n_elem) lo = tgrid.n_elem-1;
  return lo;
}

// arm-specific Breslow baselines with arm indicator (mimics R: multiply by z)
// builds dl_s, dl_c, and their cumulative Lambda_s, Lambda_c in O(N+T)
static void baselines_arm(
    const vec& e_time,
    const vec& ftime, const vec& delta, const vec& z_arm,
    const vec& w_e,   const vec& w_c,
    vec& dl_s, vec& dl_c, vec& Lam_s, vec& Lam_c)
{
  const int t = (int)e_time.n_elem;
  const int N = (int)ftime.n_elem;

  // map subjects to grid times & collect indices per time
  std::vector< std::vector<uword> > ids_at_time(t);
  ids_at_time.reserve(t);
  for(int i=0;i<N;++i){
    uword j = lower_idx(e_time, ftime[i]);
    ids_at_time[j].push_back((uword)i);
  }

  vec add_we(t, fill::zeros), add_wc(t, fill::zeros);
  vec dN(t, fill::zeros), dNc(t, fill::zeros);

  for(int j=0;j<t;++j){
    const auto &v = ids_at_time[j];
    for(uword k=0;k<v.size();++k){
      uword i = v[k];
      if(z_arm[i] > 0.5){                 // in this arm
        add_we[j] += w_e[i];
        add_wc[j] += w_c[i];
        if(delta[i]==1.0) dN[j]++; else dNc[j]++;
      }
    }
  }

  // reverse-cumsum denominators (risk sets in this arm)
  vec den_e(t, fill::zeros), den_c(t, fill::zeros);
  if(t>0){
    den_e[t-1]=add_we[t-1]; den_c[t-1]=add_wc[t-1];
    for(int j=t-2;j>=0;--j){ den_e[j]=den_e[j+1]+add_we[j]; den_c[j]=den_c[j+1]+add_wc[j]; }
  }

  // Breslow increments (guard zeros)
  auto sdiv = [](double a,double b){ return (b<=0.0 || !arma::is_finite(b)) ? 0.0 : (a/b); };
  dl_s.set_size(t); dl_c.set_size(t);
  for(int j=0;j<t;++j){ dl_s[j] = sdiv(dN[j], den_e[j]); dl_c[j] = sdiv(dNc[j], den_c[j]); }

  Lam_s = cumsum(dl_s);
  Lam_c = cumsum(dl_c);
}

// stream H, K, and the martingale integral using LEFT-LIMIT K(t-)
// returns (S_ind(t), S_cl(t)) for an arm
static std::pair< vec, vec > accumulate_arm(
    const vec& e_time,
    const vec& ftime, const vec& delta, const vec& z_arm, double pi_arm,
    const vec& w_e, const vec& w_c,
    const vec& dl_s, const vec& dl_c,
    const uvec& cl_index, const uvec& cl_size)
{
  const int t = (int)e_time.n_elem;
  const int N = (int)ftime.n_elem;
  const uword M = cl_size.n_elem;
  const double eps = 1e-12;

  // map subjects to grid times once
  std::vector< std::vector<uword> > ids_at_time(t);
  for(int i=0;i<N;++i) ids_at_time[ lower_idx(e_time, ftime[i]) ].push_back((uword)i);

  // state vectors
  vec H(N, fill::ones), K(N, fill::ones), q_mart(N, fill::zeros);
  uvec atRisk(N, fill::ones);

  vec S_ind(t, fill::zeros);
  mat S_cl_sum(M, t, fill::zeros);

  for(int j=0;j<t;++j){
    // update H to t_j (right-continuous): H *= exp(- w_e * dl_s[j])
    if(dl_s[j]!=0.0) H %= exp( - dl_s[j] * w_e );

    // capture K(t^-) BEFORE updating K at t_j
    vec K_before = K;

    // update K to t_j: K *= exp(- w_c * dl_c[j])
    if(dl_c[j]!=0.0) K %= exp( - dl_c[j] * w_c );

    // q increment using LEFT-LIMIT denom: (dN_c - dl_c * w_c * Y)/(K^- * H)
    vec denomKH = K_before % H;
    for(int i=0;i<N;++i){
      if(!arma::is_finite(denomKH[i]) || denomKH[i] < eps) denomKH[i] = eps;
    }

    // drift: -dl_c * w_c * Y /(K^- * H)
    if(dl_c[j]!=0.0){
      for(int i=0;i<N;++i) if(atRisk[i]) q_mart[i] += (- dl_c[j] * w_c[i]) / denomKH[i];
    }
    // jump for censored at j: +1/(K^- * H)
    for(uword k=0;k<ids_at_time[j].size();++k){
      uword i = ids_at_time[j][k];
      if(delta[i]==0.0) q_mart[i] += 1.0 / denomKH[i];
    }

    // column at time j: z/(π K^-)*Y  - (z-π)/π * H  +  z/π * H * q
    vec denomKminus = K_before;
    for(int i=0;i<N;++i){
      if(!arma::is_finite(denomKminus[i]) || denomKminus[i] < eps) denomKminus[i] = eps;
    }
    vec term1  = (z_arm/pi_arm) % (conv_to<vec>::from(atRisk) % (1.0 / denomKminus)); // IPCW with K^-
    vec term2  = ((z_arm - pi_arm)/pi_arm) % H;                                       // augmentation
    vec term3  = (z_arm/pi_arm) % (H % q_mart);                                       // martingale integral
    vec col    = term1 - term2 + term3;

    // individual mean
    S_ind[j] = mean(col);

    // cluster sums (sum within cluster; later turn into mean)
    for(int i=0;i<N;++i) S_cl_sum(cl_index[i], j) += col[i];

    // subjects at time j leave risk set now
    for(uword k=0;k<ids_at_time[j].size();++k) atRisk[ ids_at_time[j][k] ] = 0u;
  }

  // average of cluster means: (1/M) * Σ_m ( (1/n_m) Σ_{i∈m} col_i )
  vec S_cl(t, fill::zeros);
  for(uword m=0;m<M;++m){
    if(cl_size[m]>0) S_cl += S_cl_sum.row(m).t() / (double)cl_size[m];
  }
  S_cl /= (double)M;

  return std::make_pair(S_ind, S_cl);
}

// ---- main ------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List marginal_est(const arma::vec& ftime, const arma::vec& delta, const arma::vec trt,
                       const arma::vec& strata, const double& trt_prob1, const double& trt_prob0,
                       const arma::mat& censor_cov, const arma::mat& surv_cov,
                       const arma::vec& censor_fit1, const arma::vec& censor_fit0, 
                       const arma::vec& surv_fit1,   const arma::vec& surv_fit0,
                       const arma::vec& e_time,
                       const bool RMST_cal = true)  // <-- new flag (default: TRUE)
{
  const int N = (int)ftime.n_elem;
  const int t = (int)e_time.n_elem;

  // cluster indexing 0..M-1
  vec strata_M = arma::unique(strata);
  uword M = strata_M.n_elem;
  uvec cl_index(N);
  for (int i = 0; i < N; ++i) {
    cl_index[i] = (uword) index_min(abs(strata_M - strata[i]));
  }
  uvec cl_size(M, fill::zeros);
  for (int i = 0; i < N; ++i) cl_size[cl_index[i]]++;

  // arm weights exp(η)
  vec w_e1 = exp(surv_cov * surv_fit1);
  vec w_c1 = exp(censor_cov * censor_fit1);
  vec w_e0 = exp(surv_cov * surv_fit0);
  vec w_c0 = exp(censor_cov * censor_fit0);

  // baselines per arm
  vec dl_s1, dl_c1, Lam_s1, Lam_c1;
  vec dl_s0, dl_c0, Lam_s0, Lam_c0;
  baselines_arm(e_time, ftime, delta, trt,       w_e1, w_c1, dl_s1, dl_c1, Lam_s1, Lam_c1);
  baselines_arm(e_time, ftime, delta, 1.0 - trt, w_e0, w_c0, dl_s0, dl_c0, Lam_s0, Lam_c0);

  // accumulate S(t) for each arm
  auto S1 = accumulate_arm(e_time, ftime, delta, trt,       trt_prob1, w_e1, w_c1, dl_s1, dl_c1, cl_index, cl_size);
  auto S0 = accumulate_arm(e_time, ftime, delta, 1.0 - trt, trt_prob0, w_e0, w_c0, dl_s0, dl_c0, cl_index, cl_size);

  // pack outputs (t × 3)
  mat out_S_ind(t, 3), out_S_cluster(t, 3);
  out_S_ind.col(0)      = S1.first;    out_S_ind.col(1)      = S0.first;    out_S_ind.col(2)      = S1.first - S0.first;
  out_S_cluster.col(0)  = S1.second;   out_S_cluster.col(1)  = S0.second;   out_S_cluster.col(2)  = S1.second - S0.second;

  // monotone clamp in [0,1]
  clamp01_monotone(out_S_ind);
  clamp01_monotone(out_S_cluster);

  // If RMST is not requested, skip computations and return NULL placeholders
  if (!RMST_cal) {
    Rcpp::List out = Rcpp::List::create(
      _["S_individual"]      = out_S_ind,
      _["S_cluster"]         = out_S_cluster,
      _["Cluster_size"]      = conv_to<vec>::from(cl_size),
      _["event_time"]        = e_time
    );
    out["RMST_cluster_out"] = R_NilValue;
    out["RMST_ind_out"]     = R_NilValue;
    return out;
  }

  // ---- RMST via trapezoids (unchanged) ----
  vec prev = zeros<vec>(t);
  if (t > 1) prev(span(1, t - 1)) = e_time(span(0, t - 2));
  vec dt = e_time - prev;
  mat dtX = repmat(dt, 1, 3);

  rowvec init; init << 1.0 << 1.0 << 0.0;
  mat Scl_prev, Sind_prev;
  if (t > 1) {
    Scl_prev = out_S_cluster.rows(0, t - 2); Scl_prev.insert_rows(0, init);
    Sind_prev = out_S_ind.rows(0, t - 2);    Sind_prev.insert_rows(0, init);
  } else {
    Scl_prev.set_size(1, 3); Scl_prev.row(0) = init;
    Sind_prev.set_size(1, 3); Sind_prev.row(0) = init;
  }

  mat RMST_cluster = (Scl_prev + out_S_cluster) % dtX / 2.0;
  mat RMST_ind     = (Sind_prev + out_S_ind)     % dtX / 2.0;

  mat RMST_cluster_out = cumsum(RMST_cluster, 0);
  mat RMST_ind_out     = cumsum(RMST_ind,     0);

  return Rcpp::List::create(
    _["S_individual"]      = out_S_ind,
    _["S_cluster"]         = out_S_cluster,
    _["Cluster_size"]      = conv_to<vec>::from(cl_size),
    _["event_time"]        = e_time,
    _["RMST_cluster_out"]  = RMST_cluster_out,
    _["RMST_ind_out"]      = RMST_ind_out
  );
}