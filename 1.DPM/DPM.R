# merge P-values from different sources using directional and non-directional methods
#The code is sourced from "Directional integration and pathway enrichment analysis for multi-omics data"

library(optparse)
library(ActivePathways)
library(dplyr)
library(openxlsx)

plot_theme = function(...) {
  theme_bw() +
    theme(    
      plot.title = element_text(size = 22),
      plot.caption = element_text(size = 12),
      plot.subtitle = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.text.x = element_text(size = 12,
                                 angle = 90, hjust = 1, vjust=0.5, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      ...
    )
}


opt = list()
opt$combined_pval_file = 'pvals.tsv'
opt$combined_fc_file = 'fcs.tsv'
opt$merged_pvalues_file = 'merged_pvalues.tsv'


# load the dfs
gene_expression = read.xlsx('全部多组学.xlsx', sheet="PPGs_FPKM",colNames = T,rowNames = T)
protein_expression = read.xlsx('全部多组学.xlsx', sheet="PPGs_蛋白丰度",colNames = T,rowNames = T)
phosphorylation = read.xlsx('全部多组学.xlsx', sheet="PPGs_磷酸化",colNames = T,rowNames = T)

# Identify common row names across all data frames
common_rows <- Reduce(intersect, list(rownames(gene_expression), rownames(protein_expression), rownames(phosphorylation)))

# subset to common_rows
gene_expression = gene_expression[common_rows,]
protein_expression = protein_expression[common_rows,]
phosphorylation = phosphorylation[common_rows,]

# separate into p_value and fold_change dfs
pval_df = gene_expression[,'pvalue',drop=FALSE]
colnames(pval_df) = 'rna'
pval_df$protein = protein_expression$pvalue
pval_df$phosphorylation = phosphorylation$pvalue


fc_df = gene_expression[,'log2FC',drop=FALSE]
colnames(fc_df) = 'rna'
fc_df$protein = protein_expression$log2FC
fc_df$phosphorylation = phosphorylation$log2FC

fc_df[is.na(fc_df)] = 0
pval_df[is.na(pval_df)] = 1

# merge pvals using brown's method
browns_df = merge_p_values(as.matrix(pval_df), method="Brown")

# repeat using DPM
dpm_df = merge_p_values(as.matrix(pval_df), method='DPM', scores_direction = as.matrix(fc_df), constraints_vector = c(1,1,1))

# combine and save to file
res_df = data.frame(Brown = browns_df, DPM = dpm_df)
write.table(res_df, file=opt$merged_pvalues_file, sep='\t', quote=FALSE)


# save of dfs to files
write.table(pval_df, file=opt$combined_pval_file, sep='\t', quote=FALSE)
write.table(fc_df, file=opt$combined_fc_file, sep='\t', quote=FALSE)