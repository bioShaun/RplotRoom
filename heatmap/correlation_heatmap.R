# plot heatmap using both pheatmap & superheat


library(pheatmap)
library(dplyr)
library(gplots)
library(superheat)
library(scales)

options(stringsAsFactors = F)

name <- 'test'

cor_file <- 'correlation_heatmap.data'
cor_df <- read.delim(cor_file, row.names = 1, check.names = F)

tissue_sample <- 'tissue_sample'
tissue_df <- read.delim(tissue_sample, header = F)
tissue_df <- filter(tissue_df, V2 %in% colnames(cor_df))

color_map_file <- 'color_theme_tissue'
color_cfg <- read.table(color_map_file, header = F,
  comment.char = "", stringsAsFactors = F)
group_colors <- color_cfg$V2
names(group_colors) <- color_cfg$V1

ann_color = data.frame(group = tissue_df$V1)
rownames(ann_color) <- tissue_df$V2
ann_colors = list(group = group_colors)


sample_num = length(colnames(cor_df))
heatmap_width <- (sample_num - 5)/5 + 7
heatmap_heigh <- (sample_num - 5)/5 + 6
fontsize = (sample_num - 5)/10 + 7
cellwidth <- (heatmap_width - 1) * 50/sample_num


pdf(paste(name, "correlation.heatmap.pheatmap.pdf", sep = "."), width = (heatmap_width), height = (heatmap_heigh ), onefile = F)
pheatmap(cor_df, annotation_col = ann_color, annotation_colors = ann_colors,
  annotation_names_row = F, annotation_names_col = F, treeheight_row = 0,
  border_color = NA, cellwidth = cellwidth, cellheight = cellwidth, fontsize = fontsize)
dev.off()


png(paste(name, "correlation.heatmap.pheatmap.png", sep = "."), width = (heatmap_width *
    300), height = (heatmap_heigh * 300), res = 300)
pheatmap(cor_df, annotation_col = ann_color, annotation_colors = ann_colors,
  annotation_names_row = F, annotation_names_col = F,treeheight_row = 0,
  border_color = NA, cellwidth = cellwidth, cellheight = cellwidth,
  fontsize = fontsize)
dev.off()

sh_colors <- group_colors[tissue_df$V1]
names(sh_colors) <- tissue_df$V2


superheatmap <- superheat(cor_df, 
  bottom.label.text.angle = 270,
  left.label.col = "white",
  left.label.text.size = 2,
  bottom.label.text.size = 2,
  col.dendrogram = TRUE,
  pretty.order.rows = TRUE,
  print.plot = F)

sample_in_row <- rownames(cor_df)[superheatmap$order.cols]
row_col <- sh_colors[sample_in_row]

plot_cor_df <- cor_df[, sample_in_row]

png(paste(name, "correlation.heatmap.superheat.png", sep = "."), 
  height = 900, width = 900)
superheat(plot_cor_df, 
  bottom.label.text.angle = 270,
  left.label.col = "white",
  left.label.text.size = 2,
  bottom.label.text.size = 2,
  col.dendrogram = TRUE,
  pretty.order.rows = TRUE,
  bottom.label.col = row_col)
dev.off()

pdf(paste(name, "correlation.heatmap.superheat.pdf", sep = "."), 
  height = 12, width = 12, onefile = F)
superheat(plot_cor_df, 
  bottom.label.text.angle = 270,
  left.label.col = "white",
  left.label.text.size = 2,
  bottom.label.text.size = 2,
  col.dendrogram = TRUE,
  pretty.order.rows = TRUE,
  bottom.label.col = row_col)
dev.off()
