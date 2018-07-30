library(pheatmap)
library(superheat)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(corrplot)
library(reshape2)
library(ggplot2)
library(tibble)

#genes <- c('lncRNA', 'TUCP', 'miRNA', 'circRNA', 'protein_coding')
genes <- c('miRNA')

gene_num_df <- read.table('./data/gene_num.txt', header = T)
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

theme_heatmap <- function(base_size = 9) {
  theme_bw(base_size = base_size) + theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(), 
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      size = base_size* 0.8,
      angle = 45,
      hjust = 0,
      face = 'bold'
    ),
    axis.text.y = element_text(
      face = 'bold',
      size = base_size* 0.8
    )
  )
} 

melt_diff_matrix <- function(diff_matrix, plot_order, type='up'){
  if (type == 'up') {
    diff_matrix[lower.tri(diff_matrix)] <- NA
  } else if (type == 'down') {
    diff_matrix[upper.tri(diff_matrix)] <- NA
  } else {
    diff_matrix <- diff_matrix
  }
  diff_matrix <- data.frame(diff_matrix, check.names = F)
  diff_matrix$group_id <- rownames(diff_matrix)
  m_diff_matrix_por <- melt(diff_matrix, id.vars = 'group_id')
  
  if (type == 'down')  {
    colnames(m_diff_matrix_por) <- c('comp2', 'comp1', 'number')
  } else {
    colnames(m_diff_matrix_por) <- c('comp1', 'comp2', 'number')
  }
  
  m_diff_matrix_por$comp1 <- factor(m_diff_matrix_por$comp1,
    levels = rev(plot_order))
  m_diff_matrix_por$comp2 <- factor(m_diff_matrix_por$comp2,
    levels = plot_order)
  m_diff_matrix_por <- m_diff_matrix_por[!(is.na(m_diff_matrix_por$number)),]
  return(m_diff_matrix_por)
}


plot_heatmap <- function(m_diff_matrix_por, m_diff_matrix){
  p <- ggplot(m_diff_matrix_por, aes(comp2, comp1)) +
    geom_tile(aes(fill = number), color='white') +
    # scale_fill_gradient(low = "white", high = "steelblue") +
    scale_fill_gradientn(colors = color1, limits=c(0,1)) +
    geom_text(data=m_diff_matrix, aes(comp2, comp1, label = number), size=1.8 )
  
  p <- p + labs(x = "",y = "") + 
    scale_x_discrete(expand = c(0, 0), position = "top") +
    scale_y_discrete(expand = c(0, 0), position = 'right') 
  
  p <- p + theme_heatmap()
  
  return(p)
}


diff_heatmap_updown <- function(name, num_stats='diff_num', plot='all') {
  diff_matrix_file = file.path('./data', paste(name, 'diff.matrix.txt', sep='.'))
  diff_matrix <- read.delim(diff_matrix_file, check.names = F, row.names = 1)
  # gene_type_num <- filter(gene_num_df, gene_biotype == name)
  # rownames(gene_type_num) <- gene_type_num$tissue
  # total_num <- gene_type_num[rownames(diff_matrix), num_stats]
  total_num <- filter(gene_num_df, gene_type == name)[, num_stats]
  plot_order <- rownames(diff_matrix_por)
 
  
  if (plot == 'updown') {
    each_type = 'all'
    diff_matrix_por <- data.frame(t(t(diff_matrix) / total_num), check.names = F)
    out_diff_matrix_por <- rownames_to_column(diff_matrix_por, var = 'Group') 
    write.table(out_diff_matrix_por, file = paste(name, 'diff2', num_stats, 'matrix.txt', sep = '.'),
      sep = '\t', quote = F, row.names = F)    
    t_diff_matrix_por <- melt_diff_matrix(diff_matrix_por, plot_order, type=each_type)
    t_diff_matrix <- melt_diff_matrix(diff_matrix, plot_order, type=each_type)
    lp <- plot_heatmap(t_diff_matrix_por, t_diff_matrix)
    ggsave(filename = paste(name, 'diff.genes.png', sep = '.'),
      width = 8, height = 8)
    ggsave(filename = paste(name, 'diff.genes.pdf', sep = '.'),
      width = 8, height = 8)    
  } else {
    total_diff_matrix <- diff_matrix + t(diff_matrix)
    total_diff_matrix_por <- total_diff_matrix / total_num
    out_total_diff_matrix_por <- rownames_to_column(total_diff_matrix_por, var = 'Group')
    write.table(out_total_diff_matrix_por, file = paste(name, 'total.diff2', num_stats, 'matrix.txt', sep = '.'),
      sep = '\t', quote = F, row.names = F) 
    each_type = 'up'
    total_m_diff_matrix_por <- melt_diff_matrix(total_diff_matrix_por, plot_order, type=each_type)
    total_m_diff_matrix <- melt_diff_matrix(total_diff_matrix, plot_order, type=each_type)
    tlp <- plot_heatmap(total_m_diff_matrix_por, total_m_diff_matrix)
    ggsave(filename = paste(name, 'total.diff2', num_stats, 'png', sep = '.'),
      width = 8, height = 8)
    ggsave(filename = paste(name, 'total.diff2', num_stats, 'pdf', sep = '.'),
      width = 8, height = 8)     
  }
  

}


genes <- c('lncRNA', 'TUCP', 'miRNA', 'circRNA', 'protein_coding')
#genes <- c('miRNA')

diff_heatmap_updown('miRNA')
sapply(genes, diff_heatmap_updown)








