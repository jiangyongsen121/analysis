#GSEA plot
library(GseaVis)
Path = "drought_GSEA/"
intergrated <- readGseaFile(filePath = Path)
setid <- c("PHENYLPROPANOID BIOSYNTHESIS",
           "MONOLIGNOL BIOSYNTHESIS",
           "FLAVONOID BIOSYNTHESIS")
# multiple terms
pdf('gsea_plot.pdf',width=8,height=6)
gseaNb(filePath = Path,
       geneSetID = setid,
       curveCol = ggsci::pal_npg()(4),
       addPval = T,
       pvalX = 1,
       pvalY = 1,
       rankCol=c('red'))
dev.off()
