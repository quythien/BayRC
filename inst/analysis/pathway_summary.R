################################################################################
# Shared pathway-summary printer for the case-study scripts.
################################################################################

print_pathway_summary <- function(result_obj,
                                  filter_by = c("q", "p"),
                                  cutoff = 0.2) {

  filter_by <- match.arg(filter_by)

  df <- result_obj$results

  # expected union of gain, loss and conserved, and its share of pathway size
  df$Expected_Union <- df$Expected_N_Gain + df$Expected_N_Loss + df$Expected_N_Conserved
  df$Union_Size_Ratio <- df$Expected_Union / df$size

  if (filter_by == "q") {
    df_sig <- df[df$padj < cutoff, ]
  } else {
    df_sig <- df[df$pval < cutoff, ]
  }

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
