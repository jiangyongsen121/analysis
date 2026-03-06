#The code is sourced from "Directional integration and pathway enrichment analysis for multi-omics data"
library(ActivePathways)
library(dplyr)
gene_expression_file = 'overlap_rna.txt'
protein_expression_file = 'overlap_protein.txt'
phosphylation_file = 'overlap_phosphylation.txt'
combined_pval_file = 'pvals.tsv'
combined_fc_file = 'fcs.tsv'
merged_pvalues_file = 'merged_pvalues.tsv'

# load the dfs
gene_expression = read.csv(gene_expression_file, sep='\t', row.names=1)
protein_expression = read.csv(protein_expression_file, sep='\t', row.names=1)
phosphylation = read.csv(phosphylation_file, sep='\t', row.names=1)

# Identify common row names across all data frames
common_rows <- Reduce(intersect, list(rownames(gene_expression), rownames(protein_expression), rownames(phosphylation)))

# subset to common_rows
gene_expression = gene_expression[common_rows,]
protein_expression = protein_expression[common_rows,]
phosphylation = phosphylation[common_rows,]

# separate into p_value and fold_change dfs
pval_df = gene_expression[,'pvalue',drop=FALSE]
colnames(pval_df) = 'rna'
pval_df$protein = protein_expression$pValue
pval_df$phosphylation = phosphylation$pValue

fc_df = gene_expression[,'log2FoldChange',drop=FALSE]
colnames(fc_df) = 'rna'
fc_df$protein = protein_expression$log2FC
fc_df$phosphylation = phosphylation$Log2FC
fc_df[is.na(fc_df)] = 0
pval_df[is.na(pval_df)] = 1
# merge pvals using brown's method
browns_df = merge_p_values(as.matrix(pval_df), method="Brown")
# repeat using DPM
dpm_df = merge_p_values(as.matrix(pval_df), 
                        method='DPM', 
                        scores_direction = as.matrix(fc_df), 
                        constraints_vector = c(-1,-1,-1)) 
# combine and save to file
res_df = data.frame(Brown = browns_df, DPM = dpm_df)
write.table(res_df, file=merged_pvalues_file, sep='\t', quote=FALSE)
# save of dfs to files
write.table(pval_df, file=combined_pval_file, sep='\t', quote=FALSE)
write.table(fc_df, file=combined_fc_file, sep='\t', quote=FALSE)

library("ggplot2")
library("data.table")
library('dplyr')
library('ggpubr')
library("ggrepel")
library("ActivePathways")

# load differential expression p-values
df_pvals <- read.table("pvals.tsv",sep="\t",header=TRUE)
scores <- data.frame(row.names = row.names(df_pvals), protein=df_pvals$protein,
                     rna=df_pvals$rna, phosphylation=df_pvals$methylation)
scores <- as.matrix(scores)
scores[is.na(scores)] <- 1

# load differential expression log2fc values
df_fcs <- read.table("fcs.tsv", sep = "\t", header=TRUE)
scores_dir <- data.frame(row.names = row.names(df_fcs), protein=df_fcs$protein,
                         rna = df_fcs$rna, phosphylation = df_fcs$methylation)
scores_dir <- as.matrix(scores_dir)
scores_dir[is.na(scores_dir)] <- 0

gene_list <- read.table('KEGG_BP_gene_list.txt',header = F)
scores <- scores[match(gene_list$V1,row.names(scores)),]
scores_dir <- scores_dir[match(gene_list$V1,row.names(scores_dir)),]
constraints_vector <- c(1,1,1)
ap_background <- rownames(scores)
fname_GMT = "KEGG_BP_phospho.gmt"
res_brown <- ActivePathways(scores=scores,
                            background=ap_background,
                            merge_method="Brown",
                            gmt=fname_GMT,
                            cytoscape_file_tag="Brown_",
                            custom_colors=c("#b2abd2","#fdb863","#66c2a5"),
                            color_integrated_only="#bababa",
                            geneset_filter=c(10,750),
                            cutoff=0.10,
                            significant=0.05
)
res_dpm <- ActivePathways(scores=scores,
                          scores_direction=scores_dir,
                          constraints_vector=constraints_vector,
                          background=ap_background,
                          merge_method="DPM",
                          gmt=fname_GMT,
                          cytoscape_file_tag="DPM_",
                          custom_colors=c("#b2abd2","#fdb863","#66c2a5"),
                          color_integrated_only="#bababa",
                          geneset_filter=c(10,750),
                          cutoff=0.10,
                          significant=0.05
)

# save pathway results
export_as_CSV(res_brown,"res_brown.csv")
export_as_CSV(res_dpm,"res_dpm.csv")
save(res_brown,res_dpm,scores,scores_dir,file ="panelF_Rdata.Rdata")

### aggregate the res_brown pathways with the res_dpm pathways, while 
### keeping track of lost/gained/maintained pathways between the two methods. 

# 1) Write the aggregated subgroups txt file
col_colors <- c("#b2abd2","#fdb863","#66c2a5","#bababa")
tests <- c('protein','rna','phosphylation','combined')
all_pathways <- data.frame(term_id = c(res_brown$term_id,res_dpm$term_id))
evidence <- append(res_brown$evidence, res_dpm$evidence)
all_pathways$evidence <- evidence
head(all_pathways)

# retain the evidence contribution from each method (contributions for shared pathways will be the same)
sub_df = all_pathways %>% group_by(term_id) %>% summarise(evidence = list(evidence))
sub_df$evidence <- lapply(sub_df$evidence,'[[',1)

col_significance <- sub_df
evidence_columns = do.call(rbind, lapply(col_significance$evidence,
                                         function(x) 0+(tests %in% x)))
colnames(evidence_columns) = tests
col_significance = cbind(col_significance[,"term_id"], evidence_columns)
head(col_significance)

#check for lost/gained/maintained pathways between methods
lostp <- res_brown$term_id[!res_brown$term_id %in% res_dpm$term_id]
gainedp <- res_dpm$term_id[!res_dpm$term_id %in% res_brown$term_id]
sharedp <- res_brown$term_id[res_brown$term_id %in% res_dpm$term_id]
col_significance$directional_impact <- 0
col_significance[col_significance$term_id %in% lostp,]$directional_impact <- 1
col_significance[col_significance$term_id %in% gainedp,]$directional_impact <- 2


col_significance <- as.data.table(col_significance)
instruct.str <- paste('piechart:',
                      ' attributelist="', 
                      paste(tests, collapse=','),
                      '" colorlist="', 
                      paste(col_colors, collapse=','), 
                      '" showlabels=FALSE', sep='')
col_significance[, "instruct" := instruct.str]
utils::write.table(col_significance, 
                   file=paste0("COMBINED__", "subgroups.txt"), 
                   row.names=FALSE, 
                   sep="\t", 
                   quote=FALSE)

# 2) Write the aggregated pathways txt file
df_txtpathways <- data.frame(term_id = c(res_brown$term_id,res_dpm$term_id),
                             term_name = c(res_brown$term_name,res_dpm$term_name),
                             adjusted_p_val = c(res_brown$adjusted_p_val,res_dpm$adjusted_p_val))

sub_df2 = df_txtpathways %>% group_by(term_id,term_name) %>% summarise(adjusted_p_val = min(adjusted_p_val))
sub_df2
pathways_txt <- as.data.table(sub_df2)

utils::write.table(pathways_txt, 
                   file=paste0("COMBINED__", "pathways.txt"), 
                   row.names=FALSE, 
                   sep="\t", 
                   quote=FALSE)

# 3) Write the aggregated gmt file
gmt_main <- read.GMT("KEGG_BP_phospho.gmt")
gmt_main <- gmt_main[pathways_txt$term_id]
write.GMT(gmt_main,paste0("COMBINED__","pathways.gmt"))
