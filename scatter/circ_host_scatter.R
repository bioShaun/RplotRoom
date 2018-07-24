library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)


circ_prop <- fread('./data/ciri.circ_host.tpm.all.txt')
circ_prop$log_host_tpm <- log10(circ_prop$host_tpm)
circ_prop$log_circ_tpm <- log10(circ_prop$circRNA_tpm)
circ_prop$log_ratio <- log10(circ_prop$junction_reads_ratio)

p <- ggplot(circ_prop, aes(log_host_tpm, log_circ_tpm, color=group_id)) +
  geom_point() +
  facet_wrap(~group_id, nrow = 5) +
  xlab('log10 Host TPM') + ylab('log10 circRNA TPM')
p

ggsave('CIRI_tpm_compare.png', plot = p, width = 18, height = 12,
  dpi = 300, type = "cairo")
ggsave('CIRI_tpm_compare.pdf', plot = p, width = 18, height = 12,
  device = cairo_pdf)

p <- ggplot(circ_prop, aes(log_host_tpm, log_ratio, color=group_id)) +
  geom_point() +
  facet_wrap(~group_id, nrow = 5) +
  xlab('log10 Host TPM') + ylab('log10 circRNA ratio')
p
ggsave('CIRI_ratio_compare.png', plot = p, width = 18, height = 12,
  dpi = 300, type = "cairo")
ggsave('CIRI_ratio_compare.pdf', plot = p, width = 18, height = 12,
  device = cairo_pdf)
