suppressMessages(library(argparser))
suppressMessages(library(edgeR))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))

options(stringsAsFactors = F)

p <- arg_parser("perform differential analysis")
p <- add_argument(
  p, '--counts', 
  help = 'counts file')
p <- add_argument(
  p, '--tpm_table', 
  help = 'quantification normalized tpm table')
p <- add_argument(
  p, '--compare', 
  help = 'compare name.', 
  default = 'Case:Control')
p <- add_argument(
  p, '--samples', 
  help = 'samples to perform diff analysis, example: case_1,case_2:control1,control2')
p <- add_argument(
  p, '--out_dir', 
  help = 'output directory')
p <- add_argument(
  p, '--qvalue', 
  help = 'diff gene qvalue cutoff', 
  default = 0.05)
p <- add_argument(
  p, '--logfc', 
  help = 'diff gene logfc cutoff',  
  default = 1)
p <- add_argument(
  p, '--dispersion', 
  help = 'for analysis with no replicates',  
  default = 0.1)
p <- add_argument(
  p, '--lib_size', 
  help = 'user provided library size', 
  default = NULL)
argv <- parse_args(p)


check_file <- function(input_file) {
  if (! file.exists(input_file)) {
    stop(paste(input_file, 'Not EXISTS!'))
  }
}

check_samples <- function(samples, table_df) {
  non_exist_samples <- samples[(! samples %in% colnames(table_df))]
  if (length(non_exist_samples) > 0) {
    stop(paste(non_exist_samples, 'Not in expression table!'))
  }
}

# read aguments
counts_file <- argv$counts
tpm_table <- argv$tpm_table
compare <- argv$compare
samples <- argv$samples
outdir <- argv$out_dir
qvalue <- argv$qvalue
logfc <- argv$logfc
dispersion <- argv$dispersion
lib_size <- argv$lib_size

## read expression tables
dir.create(outdir, showWarnings = FALSE)
check_file(counts_file)
check_file(tpm_table)
gene_tpm_matrix <- fread(tpm_table, check.names = F)
gene_counts <- fread(counts_file, check.names = F)
gene_tpm_matrix <- data.frame(gene_tpm_matrix)
gene_counts <- data.frame(gene_counts)

## load compare and sample information
each_pair <- unlist(strsplit(compare, split = ":"))
con1_sample <- unlist(strsplit(strsplit2(samples, split=':')[1], split=','))
con2_sample <- unlist(strsplit(strsplit2(samples, split=':')[2], split=','))
samples <- c(con1_sample, con2_sample)
check_samples(samples, gene_counts)
check_samples(samples, gene_tpm_matrix)

## extract analysis sample's expression table
each_pair_cts <- gene_counts[, samples]
rownames(each_pair_cts) <- gene_counts[,1]
each_pair_cts <- each_pair_cts[rowSums(each_pair_cts)>=2,]

## load counts into edgeR DEGList and normalization
conditions = factor(c(rep(each_pair[1], length(con1_sample)),
  rep(each_pair[2], length(con2_sample))))
if (is.na(lib_size)) {
    y <- DGEList(each_pair_cts, group=conditions)
} else {
    lib_size_df <- read.delim(lib_size, stringsAsFactors = F,
                              header = F, row.names = 1)
    ep_lib_size <- lib_size_df[c(con1_sample, con2_sample),]
    y <- DGEList(each_pair_cts, group=conditions, lib.sizes <- ep_lib_size)
}
y <- calcNormFactors(y)

## diff analysis
if (length(con1_sample) > 1 && length(con2_sample) > 1) {
  y <- estimateDisp(y)
  et <- exactTest(y, pair=rev(each_pair))
} else {
  et <- exactTest(y, pair=rev(each_pair), dispersion=dispersion)
}

## merge DE table with tpm table and output
tTags <- topTags(et,n=NULL)
new_tTags <- tTags$table
new_tTags <- new_tTags[, !(names(new_tTags) %in% c("logCPM"))]
each_pair_matrix <- gene_tpm_matrix[, samples]
rownames(each_pair_matrix) <- gene_tpm_matrix[, 1]
merged_df <- merge(each_pair_matrix, new_tTags, by.x = 0, by.y = 0, all.y = T)
sorted_merged_df <- arrange(merged_df, FDR)
colnames(sorted_merged_df)[1] <- 'Gene_ID'
out_file_name_prefix <- paste(outdir, '/', each_pair[1], '_vs_', each_pair[2], sep = '')
write.table(sorted_merged_df,
            file=paste(out_file_name_prefix, 'edgeR.DE_results.txt', sep = '.'),
            sep='\t', quote=F, row.names=F)


## up & down regulated genes results ##
# diff_genes <- c(diff_genes, up_regulate_df$Gene_ID, down_regulate_df$Gene_ID)
# up_regulate_df <- filter(sorted_merged_df, logFC >= logfc, FDR <= qvalue)
# down_regulate_df <- filter(sorted_merged_df, logFC <= -(logfc), FDR <= qvalue)
# diff_genes <- c()
# up_regulate_name_prefix <- paste(outdir, '/', compare, '.',each_pair[1], '-UP',  sep = '')
# down_regulate_name_prefix <- paste(outdir, '/', compare, '.',each_pair[2], '-UP', sep = '')
# if (dim(up_regulate_df)[1] > 0) {
#     write.table(up_regulate_df,
#                 file=paste(up_regulate_name_prefix, 'edgeR.DE_results.txt', sep = '.'),
#                 sep='\t', quote=F, row.names=F)
#     write(as.character(up_regulate_df$Gene_ID),
#           file = paste(up_regulate_name_prefix, 'edgeR.DE_results.diffgenes.txt', sep = '.'),
#           sep = '\n')
# }
# if (dim(down_regulate_df)[1] > 0) {
#     write.table(down_regulate_df,
#                 file=paste(down_regulate_name_prefix, 'edgeR.DE_results.txt', sep = '.'),
#                 sep='\t', quote=F, row.names=F)
#     write(as.character(down_regulate_df$Gene_ID),
#           file = paste(down_regulate_name_prefix, 'edgeR.DE_results.diffgenes.txt', sep = '.'),
#           sep = '\n')
# }
# 
# ## write diff gene list
# if (length(diff_genes) > 0) {
#     write(as.character(diff_genes),
#           file = paste(out_file_name_prefix, 'ALL.edgeR.DE_results.diffgenes.txt', sep = '.'),
#           sep = '\n')
# }
