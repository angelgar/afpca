#' Firing rates for 25 Neurons in the Motor Cortex
#'
#' This dataset comes from a novel neuroscience research to understand the
#' dynamics between motor cortex activation and its involvement in skilled
#' movements. It contains firing rates for 25 neurons coming from the motor cortex
#' of a mouse while reaching for a food pellet. These are neuron-specific averages
#' across 147 experimental trials.
#'
#' @format A data frame with 25 rows and two columns including:
#' \describe{
#'   \item{neuron_no}{Integer, Neuron ID}
#'   \item{activation}{A tidyfun tfd object containing firing rates for each neuron across 174 timepoints, in spikes / ms}
#'   ...
#' }
#'
#' @references
#'
#' Sauerbrei, Britton A., et al. "Cortical pattern generation during dexterous
#' movement is input-driven." \emph{Nature} 577.7790 (2020): 386-391.
#'
#' tidyfun manual: tidyfun.github.io/tidyfun
#'
"firing_rates_data"
