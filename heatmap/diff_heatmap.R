library(pheatmap)
library(superheat)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(corrplot)
library(reshape2)
library(ggplot2)

#genes <- c('lncRNA', 'TUCP', 'miRNA', 'circRNA', 'protein_coding')
genes <- c('miRNA')

gene_num_df <- read.table('./data/all.detected_genes.all_samples.txt', header = T)
name <- 'miRNA'
num_stats='detected_number'
gene_type_num <- filter(gene_num_df, gene_biotype == name)
rownames(gene_type_num) <- gene_type_num$tissue

# color
color_map_file <- './data/tissue_group.color.txt'
color_cfg <- read.table(color_map_file, header = T,
  comment.char = "", stringsAsFactors = F)
group_col <- color_cfg[, c('group', 'color')]
group_col <- group_col[!duplicated(group_col),]
rownames(group_col) <- group_col$group

color1_part1 <- colorRampPalette(brewer.pal(n = 7, name = "Blues"))(21)
color1_part2 <- colorRampPalette(brewer.pal(n = 9, name = "Blues")[8:9])(40)
color1 <- c(color1_part1, color1_part2)
color2_part1 <- colorRampPalette(viridis::viridis(10)[1:8])(21)
color2_part2 <- colorRampPalette(viridis::viridis(10)[9:10])(40)
color2 <- c(color2_part1, color2_part2)
out_suf2 <- 'viridis'
out_suf1 <- 'blues'

color_pal=color1
diff_heatmap <- function(name, num_stats='detected_number', 
  color_pal=color1, suffix = out_suf1) {
  diff_matrix_file = file.path('./data', paste(name, 'diff.matrix.txt', sep='.'))
  diff_matrix <- read.delim(diff_matrix_file, row.names = 1, check.names = F)
  total_num <- gene_type_num[rownames(diff_matrix), num_stats]
  # total_num <- filter(gene_num_df, gene_type == name)[, num_stats]
  diff_matrix_por <- diff_matrix / total_num
  color_order <- group_col[colnames(diff_matrix),]$color
  pdf(paste(name, suffix, 'tissue.pairwise.diff.genes.pdf', sep='.'), width = 12, height = 12, onefile = F)
  # pheatmap(diff_matrix, cluster_rows = F, cluster_cols = F, 
  #   display_numbers=T, number_format='%.2f',
  #   color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
  #   border_color = F, number_color = "red",
  #   legend=F)
  # diff_matrix_por <- as.matrix(diff_matrix_por)
  # diff_matrix_por[lower.tri(diff_matrix_por)] <- NA
  # diff_matrix <- as.matrix(diff_matrix)
  # diff_matrix[lower.tri(diff_matrix)] <- NA
  superheat(diff_matrix_por, 
    heat.lim = c(0, 0.6),
    heat.pal.values = c(0, 0.5, 1),
    heat.pal = color_pal,
    X.text = diff_matrix,
    X.text.size = 2,
    bottom.label.text.angle = 90,
    bottom.label.text.size = 3,
    left.label.text.size = 3,
    bottom.label.col = color_order,
    bottom.label.size = 0.2,
    left.label.col = 'white',
    extreme.values.na = FALSE,
    legend.breaks = seq(0, 0.6, 0.1))
 
  dev.off()

  png(paste(name, suffix, 'tissue.pairwise.diff.genes.png', sep='.'), width = 12, height = 12,
    units = 'in', res = 300)
  # pheatmap(diff_matrix, cluster_rows = F, cluster_cols = F, 
  #   display_numbers=T, number_format='%d',
  #   color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
  #   border_color = F, number_color = "red",
  #   legend=F)
  superheat(diff_matrix_por, 
    heat.lim = c(0, 0.6),
    heat.pal.values = c(0, 0.5, 1),
    heat.pal = color_pal,
    X.text = as.matrix(diff_matrix),
    X.text.size = 2,
    bottom.label.text.angle = 90,
    bottom.label.text.size = 3,
    left.label.text.size = 3,
    bottom.label.col = color_order,
    bottom.label.size = 0.2,
    left.label.col = 'white',
    extreme.values.na = FALSE,
    legend.breaks = seq(0, 0.6, 0.1))
  dev.off()
}

#diff_heatmap('TUCP')
sapply(genes, diff_heatmap)

# try corrplot
# colmat <- colorRampPalette(c("red", "white", "blue"))
# 
# corrplot(as.matrix(diff_matrix_por), type="upper",
#   addCoef.col = "grey", method = "color", 
#   is.corr = F,cl.lim = c(0, 1), col = colmat(200)
# )

# try ggplot2
# (https://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/)
diff_matrix <- read.delim(diff_matrix_file, check.names = F)
m_diff_matrix_por <- melt(diff_matrix)
colnames(m_diff_matrix_por) <- c('comp1', 'comp2', 'number')
p <- ggplot(m_diff_matrix_por, aes(comp1, comp2)) +
  geom_tile(aes(fill = number), color='white') +
  scale_fill_gradient(low = "white", high = "steelblue")
p
