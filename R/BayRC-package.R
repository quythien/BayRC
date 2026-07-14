#' @description
#' BayRC (Bayesian Rhythmicity Comparison) is a unified statistical framework
#' for comparing and interpreting circadian rhythms across biological
#' conditions: age, disease state, tissue, species, or sex. It jointly
#' infers gene-level rhythmicity and phase, computes posterior probabilities
#' of rhythmic and phase concordance, and classifies rhythmic gain, loss,
#' conservation, and direction-specific phase shifts under Bayesian false
#' discovery rate (BFDR) control. It further supports pathway-level
#' enrichment and genome-wide concordance scoring, providing a unified,
#' uncertainty-aware framework for comparative circadian analysis across
#' tissues, species, and disease contexts.
#'
#' @section Main functions:
#' \describe{
#'   \item{MCMC core (paper Sec. 2.1)}{
#'     \itemize{
#'       \item \code{\link{CB_init_single}}
#'       \item \code{\link{CBt_init_single}}
#'       \item \code{\link{CB_MCMC_single_rj_slice}}
#'       \item \code{\link{CB_getAllEst}}
#'       \item \code{\link{CBt_sim_data}}
#'       \item \code{\link{Cosinor_fit}}
#'       \item \code{\link{circular_HDI}}
#'       \item \code{\link{circular_median}}
#'     }
#'   }
#'   \item{Gene-level biomarker detection with BFDR control (paper Sec. 2.2)}{
#'     \itemize{
#'       \item \code{\link{match_symbols}}
#'       \item \code{\link{bfdr_from_posterior}}
#'       \item \code{\link{detect_rhy}}
#'       \item \code{\link{summarize_bay}}
#'       \item \code{\link{transition_classify}}
#'       \item \code{\link{phase_infer}}
#'     }
#'   }
#'   \item{Pathway-level rhythmic enrichment (paper Sec. 2.3)}{
#'     \itemize{
#'       \item \code{\link{pathSelect}}
#'       \item \code{\link{plot_heatmap}}
#'       \item \code{\link{multi_conservation_pathway}}
#'       \item \code{\link{multi_conservation_pathway_bootstrap}}
#'     }
#'   }
#'   \item{Genome-wide concordance summary (paper Sec. 2.4)}{
#'     \itemize{
#'       \item \code{\link{multi_conservation}}
#'       \item \code{\link{pairwise_concordance}}
#'     }
#'   }
#'   \item{Cross-species alignment}{
#'     \itemize{
#'       \item \code{\link{match_homologs}}
#'       \item \code{\link{merge_mcmc}}
#'     }
#'   }
#' }
#' @importFrom Rcpp evalCpp
#' @import stats
#' @import graphics
#' @import grDevices
#' @import utils
#' @importFrom circular circular median.circular mean.circular sd.circular
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_bar geom_histogram
#'   geom_density geom_segment geom_ribbon geom_errorbar geom_tile
#'   facet_wrap facet_grid theme theme_bw theme_minimal theme_void
#'   scale_fill_manual scale_color_manual scale_fill_gradient
#'   scale_x_continuous scale_y_continuous scale_color_gradient
#'   expansion coord_polar element_rect element_line
#'   labs xlab ylab ggtitle element_text element_blank
#'   unit margin ggsave alpha vars
#' @importFrom dplyr filter mutate select arrange group_by summarise
#'   left_join inner_join bind_rows bind_cols rename pull n
#'   ungroup distinct case_when
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map_dfr map map_dbl
#' @importFrom RColorBrewer brewer.pal
#' @useDynLib BayRC, .registration = TRUE
"_PACKAGE"
