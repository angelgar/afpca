## Code used to transform the raw firing rates from a data.frame to a tidyfun object.

library(tidyverse)
library(tidyfun)

firing_rates_data_raw <- read_csv("~/Desktop/firing_rates_data.csv") %>% as.matrix()

firing_rates_data <- tibble(neuron_no = 1:25)

firing_rates_data$activation <- tfd(firing_rates_data_raw, arg = seq(0,1, l = 174))

usethis::use_data(firing_rates_data, overwrite = TRUE)
