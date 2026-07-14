#' @keywords internal
#' @importFrom Rcpp evalCpp
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
