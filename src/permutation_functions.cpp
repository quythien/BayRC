#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_null_congruence_cpp(
    IntegerVector pathway_idx,        // 1-based indices of genes in pathway
    LogicalVector rhythmic_A,         // Full rhythmic_A vector for all genes
    LogicalVector rhythmic_B,         // Full rhythmic_B vector for all genes
    int n_perm = 1000,
    int seed = 123
) {
  // Set seed for reproducibility
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed = base_env["set.seed"];
  set_seed(seed);
  
  int n_genes_total = rhythmic_B.size();
  int pathway_size = pathway_idx.size();
  
  // Convert 1-based R indices to 0-based C++ indices
  IntegerVector idx_cpp = pathway_idx - 1;
  
  // Pre-allocate result vector
  NumericVector null_scores(n_perm);
  
  // Create a working copy of rhythmic_B for shuffling
  LogicalVector rhythmic_B_shuffled = clone(rhythmic_B);
  
  // Perform permutation loop
  for (int perm = 0; perm < n_perm; perm++) {
    
    // Shuffle rhythmic_B_shuffled using Fisher-Yates algorithm
    for (int i = n_genes_total - 1; i > 0; i--) {
      int j = rand() % (i + 1);  // Random index from 0 to i
      // Swap elements
      bool temp = rhythmic_B_shuffled[i];
      rhythmic_B_shuffled[i] = rhythmic_B_shuffled[j];
      rhythmic_B_shuffled[j] = temp;
    }
    
    // Calculate congruence for this permutation
    int both_count = 0;
    int either_count = 0;
    
    for (int i = 0; i < pathway_size; i++) {
      int gene_idx = idx_cpp[i];
      
      bool in_A = rhythmic_A[gene_idx];
      bool in_B = rhythmic_B_shuffled[gene_idx];
      
      if (in_A && in_B) {
        both_count++;
      }
      if (in_A || in_B) {
        either_count++;
      }
    }
    
    // Calculate congruence score (both / either)
    if (either_count == 0) {
      null_scores[perm] = NA_REAL;  // Return NA if no rhythmic genes
    } else {
      null_scores[perm] = (double)both_count / (double)either_count;
    }
  }
  
  return null_scores;
}


// [[Rcpp::export]]
NumericVector compute_null_congruence_cpp_v2(
    IntegerVector pathway_idx,
    LogicalVector rhythmic_A,
    LogicalVector rhythmic_B,
    int n_perm = 1000
) {
  
  int n_genes_total = rhythmic_B.size();
  int pathway_size = pathway_idx.size();
  
  // Convert 1-based R indices to 0-based C++ indices
  IntegerVector idx_cpp = pathway_idx - 1;
  
  // Pre-allocate result vector
  NumericVector null_scores(n_perm);
  
  // Create a working copy of rhythmic_B for shuffling
  LogicalVector rhythmic_B_shuffled = clone(rhythmic_B);
  
  // Perform permutation loop
  for (int perm = 0; perm < n_perm; perm++) {
    
    // Shuffle using R's RNG (better quality randomness)
    IntegerVector shuffle_order = Rcpp::sample(n_genes_total, n_genes_total, false) - 1;
    
    // Apply shuffle
    for (int i = 0; i < n_genes_total; i++) {
      rhythmic_B_shuffled[i] = rhythmic_B[shuffle_order[i]];
    }
    
    // Calculate congruence for this permutation
    int both_count = 0;
    int either_count = 0;
    
    for (int i = 0; i < pathway_size; i++) {
      int gene_idx = idx_cpp[i];
      
      bool in_A = rhythmic_A[gene_idx];
      bool in_B = rhythmic_B_shuffled[gene_idx];
      
      if (in_A && in_B) {
        both_count++;
      }
      if (in_A || in_B) {
        either_count++;
      }
    }
    
    // Calculate congruence score (both / either)
    if (either_count == 0) {
      null_scores[perm] = NA_REAL;
    } else {
      null_scores[perm] = (double)both_count / (double)either_count;
    }
  }
  
  return null_scores;
}