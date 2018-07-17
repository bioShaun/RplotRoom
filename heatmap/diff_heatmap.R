library(pheatmap)
library(RColorBrewer)

genes <- c('lncRNA', 'TUCP', 'miRNA', 'circRNA', 'protein_coding')

diff_heatmap <- function(name) {
  diff_matrix = file.path('./data', paste(name, 'diff.matrix.txt', sep='.'))
  diff_matrix <- read.delim(diff_matrix, row.names = 1)
  pdf(paste(name, 'tissue.pairwise.diff.genes.pdf', sep='.'), width = 12, height = 12, onefile = F)
  pheatmap(diff_matrix, cluster_rows = F, cluster_cols = F, 
    display_numbers=T, number_format='%d',
    color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
    border_color = F, number_color = "red",
    legend=F)
  dev.off()
  
  png(paste(name, 'tissue.pairwise.diff.genes.png', sep='.'), width = 12, height = 12,
    units = 'in', res = 300)
  pheatmap(diff_matrix, cluster_rows = F, cluster_cols = F, 
    display_numbers=T, number_format='%d',
    color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
    border_color = F, number_color = "red",
    legend=F)
  dev.off()
}

sapply(genes, diff_heatmap)
