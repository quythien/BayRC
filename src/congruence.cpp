#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>
using namespace Rcpp;

// ============================================================================
// ITERATION-LEVEL JACCARD INDEX FUNCTIONS
// ============================================================================

// Helper: Compute Jaccard index for a single iteration
inline double compute_jaccard_single(const IntegerVector& A, const IntegerVector& B) {
  int n = A.size();
  int a = 0;  // both rhythmic
  int b = 0;  // A only (loss)
  int c = 0;  // B only (gain)
  
  for (int i = 0; i < n; i++) {
    if (A[i] == 1 && B[i] == 1) a++;
    else if (A[i] == 1 && B[i] == 0) b++;
    else if (A[i] == 0 && B[i] == 1) c++;
  }
  
  int union_size = a + b + c;
  if (union_size == 0) return NA_REAL;
  return static_cast<double>(a) / union_size;
}

// Helper: Compute expected Jaccard under independence
inline double compute_null_jaccard(int n_A, int n_B, int G) {
  if (G == 0) return 0.0;
  double exp_inter = (static_cast<double>(n_A) * n_B) / G;
  double exp_union = n_A + n_B - exp_inter;
  if (exp_union <= 0.0) return 0.0;
  return exp_inter / exp_union;
}

// Main function: Compute iteration-level Jaccard (averaged adjusted Jaccard)
// [[Rcpp::export]]
List compute_iteration_jaccard_cpp(IntegerMatrix rho_A, IntegerMatrix rho_B) {
  int G = rho_A.nrow();  // number of genes
  int K = rho_A.ncol();  // number of iterations
  
  if (rho_B.nrow() != G || rho_B.ncol() != K) {
    stop("Matrices must have same dimensions");
  }
  
  NumericVector jaccard_obs(K);
  NumericVector jaccard_exp(K);
  NumericVector jaccard_adj(K);
  
  for (int k = 0; k < K; k++) {
    IntegerVector A = rho_A(_, k);
    IntegerVector B = rho_B(_, k);
    
    // Count confusion matrix elements
    int a = 0, b = 0, c = 0;
    for (int i = 0; i < G; i++) {
      if (A[i] == 1 && B[i] == 1) a++;
      else if (A[i] == 1 && B[i] == 0) b++;
      else if (A[i] == 0 && B[i] == 1) c++;
    }
    
    // Observed Jaccard
    int union_size = a + b + c;
    if (union_size == 0) {
      jaccard_obs[k] = NA_REAL;
      jaccard_exp[k] = NA_REAL;
      jaccard_adj[k] = NA_REAL;
      continue;
    }
    jaccard_obs[k] = static_cast<double>(a) / union_size;
    
    // Expected Jaccard under independence
    int n_A = a + b;  // total rhythmic in A
    int n_B = a + c;  // total rhythmic in B
    jaccard_exp[k] = compute_null_jaccard(n_A, n_B, G);
    
    // Null-adjusted Jaccard
    if (std::abs(jaccard_exp[k]) < 1.0) {
      jaccard_adj[k] = (jaccard_obs[k] - jaccard_exp[k]) / (1.0 - jaccard_exp[k]);
    } else {
      jaccard_adj[k] = NA_REAL;
    }
  }
  
  // Compute means (removing NAs)
  double sum_obs = 0.0, sum_adj = 0.0;
  int count_obs = 0, count_adj = 0;
  
  for (int k = 0; k < K; k++) {
    if (!NumericVector::is_na(jaccard_obs[k])) {
      sum_obs += jaccard_obs[k];
      count_obs++;
    }
    if (!NumericVector::is_na(jaccard_adj[k])) {
      sum_adj += jaccard_adj[k];
      count_adj++;
    }
  }
  
  double mean_obs = (count_obs > 0) ? sum_obs / count_obs : NA_REAL;
  double mean_adj = (count_adj > 0) ? sum_adj / count_adj : NA_REAL;
  
  return List::create(
    _["jaccard_obs"] = jaccard_obs,
    _["jaccard_exp"] = jaccard_exp,
    _["jaccard_adj"] = jaccard_adj,
    _["mean_jaccard_obs"] = mean_obs,
    _["mean_jaccard_adj"] = mean_adj
  );
}

// Permutation test for iteration-level Jaccard
// [[Rcpp::export]]
NumericVector permutation_jaccard_cpp(IntegerMatrix rho_A,
                                      IntegerMatrix rho_B,
                                      int n_perm,
                                      int seed = 12345) {
  int G = rho_A.nrow();
  int K = rho_A.ncol();
  
  NumericVector null_dist(n_perm);
  std::mt19937 rng(seed);
  std::vector<int> idx(G);
  std::iota(idx.begin(), idx.end(), 0);
  
  for (int perm = 0; perm < n_perm; perm++) {
    // Permute gene labels (consistent across all iterations)
    std::shuffle(idx.begin(), idx.end(), rng);
    
    // Create permuted matrix
    IntegerMatrix rho_B_perm(G, K);
    for (int k = 0; k < K; k++) {
      for (int i = 0; i < G; i++) {
        rho_B_perm(i, k) = rho_B(idx[i], k);
      }
    }
    
    // Compute mean adjusted Jaccard for this permutation
    double sum_adj = 0.0;
    int count = 0;
    
    for (int k = 0; k < K; k++) {
      IntegerVector A = rho_A(_, k);
      IntegerVector B_perm = rho_B_perm(_, k);
      
      int a = 0, b = 0, c = 0;
      for (int i = 0; i < G; i++) {
        if (A[i] == 1 && B_perm[i] == 1) a++;
        else if (A[i] == 1 && B_perm[i] == 0) b++;
        else if (A[i] == 0 && B_perm[i] == 1) c++;
      }
      
      int union_size = a + b + c;
      if (union_size == 0) continue;
      
      double J_obs = static_cast<double>(a) / union_size;
      int n_A = a + b;
      int n_B = a + c;
      double J_exp = compute_null_jaccard(n_A, n_B, G);
      
      if (std::abs(J_exp) < 1.0) {
        double J_adj = (J_obs - J_exp) / (1.0 - J_exp);
        sum_adj += J_adj;
        count++;
      }
    }
    
    null_dist[perm] = (count > 0) ? sum_adj / count : NA_REAL;
  }
  
  return null_dist;
}

// Bootstrap for iteration-level Jaccard
// [[Rcpp::export]]
List bootstrap_jaccard_cpp(IntegerMatrix rho_A,
                          IntegerMatrix rho_B,
                          int n_boot,
                          int seed = 12345) {
  int G = rho_A.nrow();
  int K = rho_A.ncol();
  
  NumericVector boot_obs(n_boot);
  NumericVector boot_adj(n_boot);
  
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> dist(0, G - 1);
  
  for (int b = 0; b < n_boot; b++) {
    // Resample genes (with replacement, across all iterations)
    IntegerMatrix rho_A_boot(G, K);
    IntegerMatrix rho_B_boot(G, K);
    
    for (int i = 0; i < G; i++) {
      int idx = dist(rng);
      for (int k = 0; k < K; k++) {
        rho_A_boot(i, k) = rho_A(idx, k);
        rho_B_boot(i, k) = rho_B(idx, k);
      }
    }
    
    // Compute mean Jaccard for bootstrap sample
    double sum_obs = 0.0, sum_adj = 0.0;
    int count_obs = 0, count_adj = 0;
    
    for (int k = 0; k < K; k++) {
      IntegerVector A = rho_A_boot(_, k);
      IntegerVector B = rho_B_boot(_, k);
      
      int a = 0, b = 0, c = 0;
      for (int i = 0; i < G; i++) {
        if (A[i] == 1 && B[i] == 1) a++;
        else if (A[i] == 1 && B[i] == 0) b++;
        else if (A[i] == 0 && B[i] == 1) c++;
      }
      
      int union_size = a + b + c;
      if (union_size == 0) continue;
      
      double J_obs = static_cast<double>(a) / union_size;
      sum_obs += J_obs;
      count_obs++;
      
      int n_A = a + b;
      int n_B = a + c;
      double J_exp = compute_null_jaccard(n_A, n_B, G);
      
      if (std::abs(J_exp) < 1.0) {
        double J_adj = (J_obs - J_exp) / (1.0 - J_exp);
        sum_adj += J_adj;
        count_adj++;
      }
    }
    
    boot_obs[b] = (count_obs > 0) ? sum_obs / count_obs : NA_REAL;
    boot_adj[b] = (count_adj > 0) ? sum_adj / count_adj : NA_REAL;
  }
  
  return List::create(
    _["boot_obs"] = boot_obs,
    _["boot_adj"] = boot_adj
  );
}

// ============================================================================
// GAIN/LOSS FUNCTIONS FOR ITERATION-LEVEL ANALYSIS
// ============================================================================

// [[Rcpp::export]]
List compute_gain_loss_iterations_cpp(IntegerMatrix rho_A, IntegerMatrix rho_B) {
  int G = rho_A.nrow();
  int K = rho_A.ncol();
  
  NumericVector exp_conserved(K);
  NumericVector exp_gained(K);
  NumericVector exp_lost(K);
  
  for (int k = 0; k < K; k++) {
    IntegerVector A = rho_A(_, k);
    IntegerVector B = rho_B(_, k);
    
    int conserved = 0, gained = 0, lost = 0;
    for (int i = 0; i < G; i++) {
      if (A[i] == 1 && B[i] == 1) conserved++;
      else if (A[i] == 1 && B[i] == 0) lost++;
      else if (A[i] == 0 && B[i] == 1) gained++;
    }
    
    exp_conserved[k] = conserved;
    exp_gained[k] = gained;
    exp_lost[k] = lost;
  }
  
  // Compute means
  double mean_conserved = mean(exp_conserved);
  double mean_gained = mean(exp_gained);
  double mean_lost = mean(exp_lost);
  
  // Compute indices
  double total_union = mean_conserved + mean_gained + mean_lost;
  
  double gain_index, loss_index, gain_loss_ratio;
  
  if (total_union == 0.0) {
    gain_index = 0.0;
    loss_index = 0.0;
    gain_loss_ratio = NA_REAL;
  } else {
    gain_index = mean_gained / total_union;
    loss_index = mean_lost / total_union;
    
    if (loss_index == 0.0) {
      gain_loss_ratio = (gain_index == 0.0) ? NA_REAL : INFINITY;
    } else {
      gain_loss_ratio = gain_index / loss_index;
    }
  }
  
  return List::create(
    _["expected_conserved"] = mean_conserved,
    _["expected_gained"] = mean_gained,
    _["expected_lost"] = mean_lost,
    _["gain_index"] = gain_index,
    _["loss_index"] = loss_index,
    _["gain_loss_ratio"] = gain_loss_ratio
  );
}

// Permutation test for gain-loss ratio
// [[Rcpp::export]]
double gainloss_pvalue_iterations_cpp(IntegerMatrix rho_A,
                                      IntegerMatrix rho_B,
                                      int n_perm,
                                      int seed = 12345) {
  int G = rho_A.nrow();
  int K = rho_A.ncol();
  
  // Compute observed gain-loss ratio
  List obs_result = compute_gain_loss_iterations_cpp(rho_A, rho_B);
  double ratio_obs = obs_result["gain_loss_ratio"];
  
  if (NumericVector::is_na(ratio_obs) || ratio_obs <= 0.0) return NA_REAL;
  
  double dev_obs = std::fabs(std::log(ratio_obs));
  
  // Permutation test
  std::mt19937 rng(seed);
  std::vector<int> idx(G);
  std::iota(idx.begin(), idx.end(), 0);
  
  int extreme = 0;
  
  for (int perm = 0; perm < n_perm; perm++) {
    std::shuffle(idx.begin(), idx.end(), rng);
    
    // Create permuted matrix
    IntegerMatrix rho_B_perm(G, K);
    for (int k = 0; k < K; k++) {
      for (int i = 0; i < G; i++) {
        rho_B_perm(i, k) = rho_B(idx[i], k);
      }
    }
    
    // Compute gain-loss ratio for permutation
    List perm_result = compute_gain_loss_iterations_cpp(rho_A, rho_B_perm);
    double ratio_perm = perm_result["gain_loss_ratio"];
    
    if (NumericVector::is_na(ratio_perm) || ratio_perm <= 0.0) continue;
    
    double dev_perm = std::fabs(std::log(ratio_perm));
    if (dev_perm >= dev_obs) extreme++;
  }
  
  return (1.0 + extreme) / (1.0 + n_perm);
}