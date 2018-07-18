library(ggplot2)
library(reshape2)
library(dplyr)


circ_prop <- read.delim('./data/circ_host_prop.txt')
circ_prop$log2_host_tpm <- log2(circ_prop$value_y)
circ_prop$log2_prop <- log2(circ_prop$host_prop)

select_data <- filter(circ_prop, variable == 'Cerebellum')
f_circ_prop <- filter(circ_prop, value_x > 1)

p <- ggplot(f_circ_prop, aes(log2_prop, log2_host_tpm, color=variable)) +
  geom_point() +
  facet_wrap(~variable, nrow = 5)
p
