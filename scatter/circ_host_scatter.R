library(ggplot2)
library(reshape2)
library(dplyr)


circ_prop <- read.delim('./data/circ_host_prop2.txt')
circ_prop$log_host_tpm <- log10(circ_prop$value_y)
circ_prop$log_circ_tpm <- log10(circ_prop$value_x)
circ_prop$log2_prop <- log2(circ_prop$host_prop)

select_data <- filter(circ_prop, variable == 'Cerebellum')
f_circ_prop <- filter(circ_prop, value_x >= 0.5)

p <- ggplot(circ_prop, aes(log_host_tpm, log_circ_tpm, color=variable)) +
  geom_point() +
  facet_wrap(~variable, nrow = 5)
p
