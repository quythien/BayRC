################################################################################
# Shared pathway-summary printer for the case-study scripts.
#
# Five of the case studies carried their own identical copy of this function
# and four more called it without defining it, so those four stopped partway
# through their enrichment section. The definition lives here now; the scripts
# source this file instead of repeating it.
################################################################################

print_pathway_summary <- function(result_obj,
                                  filter_by = c("q", "p"),
                                  cutoff = 0.2) {

  filter_by <- match.arg(filter_by)   # ensure "p" or "q"

  # Extract table
  df <- result_obj$results

  # Compute expected union and ratio
  df$Expected_Union <- df$Expected_N_Gain + df$Expected_N_Loss + df$Expected_N_Conserved
  df$Union_Size_Ratio <- df$Expected_Union / df$size

  # Decide filtering method
  if (filter_by == "q") {
    df_sig <- df[df$padj < cutoff, ]
  } else {
    df_sig <- df[df$pval < cutoff, ]
  }

  # Print formatted output
  apply(df_sig, 1, function(x) {
    cat(sprintf(
      "%s (Expected gain = %.1f; Expected loss = %.1f; Expected conserved = %.1f; Expected union = %.1f; Expected union / size = %.1f / %d = %.3f; Gain/Loss ratio = %.3f; p value = %.4g; q value = %.4g)\n\n",
      x["pathway"],
      as.numeric(x["Expected_N_Gain"]),
      as.numeric(x["Expected_N_Loss"]),
      as.numeric(x["Expected_N_Conserved"]),
      as.numeric(x["Expected_Union"]),
      as.numeric(x["Expected_Union"]),
      as.numeric(x["size"]),
      as.numeric(x["Union_Size_Ratio"]),
      as.numeric(x["Gain_Loss_Ratio_Arithmetic"]),
      as.numeric(x["pval"]),
      as.numeric(x["padj"])
    ))
  })
}
