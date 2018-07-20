library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)


circ_prop <- read.delim('./data/circ_host_prop.txt')
circ_prop$log_host_tpm <- log10(circ_prop$value_y)
circ_prop$log_circ_tpm <- log10(circ_prop$value_x)
circ_prop$log2_prop <- log2(circ_prop$host_prop)

select_data <- filter(circ_prop, variable == 'Cerebellum')
f_circ_prop <- filter(circ_prop, value_x >= 0.5)

p <- ggplot(circ_prop, aes(log_host_tpm, log_circ_tpm, color=variable)) +
  geom_point() +
  facet_wrap(~variable, nrow = 5)
p


circ_exp <- fread('./data/circRNA.exp.norm.group.txt')
m_circ_exp <- melt(as.data.frame(circ_exp))
f_m_circ_exp <- filter(m_circ_exp, value > 0)
f_m_circ_exp$log_tpm <- log10(f_m_circ_exp$value + 1)
p <- ggplot(f_m_circ_exp, aes(log_tpm, fill=variable, color=variable)) +
  geom_density(alpha=0.5) +
  facet_wrap(~variable)
p
ggsave(filename = 'circRNA_exp_distribution.png', plot = p, width = 16, 
  height = 10, type = "cairo-png")
