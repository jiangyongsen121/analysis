library(data.table)
library(tictoc) ##计算R语言命令运行的时间
library(lfmm)
library(openxlsx)

#读取环境数据
env <- read.xlsx('env.xlsx',sheet ="env",rowNames = F)
env$Longitude <- as.numeric(env$Longitude)
env$Latitude <- as.numeric(env$Latitude)

#scale
library(dplyr)
env_site <- env %>%
  distinct(Longitude, Latitude, Altitude, .keep_all = TRUE)
env_scaled <- env_site %>%
  select(-ID, -Longitude, -Latitude, -Altitude) %>%
  scale() %>%
  as.data.frame()
env_site_scaled <- cbind(env_site[ , c("Longitude","Latitude","Altitude")], env_scaled)
env_final <- env %>%
  select(ID, Longitude, Latitude, Altitude) %>%
  left_join(env_site_scaled, by = c("Longitude","Latitude","Altitude"))
row.names(env_final) <- env_final[,1] 
env_final <- env_final[,-1]


genotype <- fread("genotype.txt",h = T, stringsAsFactors = F)
genotype <- as.data.frame(genotype)
row.names(genotype) <- genotype[,1] 
genotype <- genotype[,-1]
genotype <- genotype[match(row.names(env_final), row.names(genotype)),]
identical(rownames(genotype), row.names(env_final))

fwrite(genotype, "snp.imupted.maf0.05.geno0.2_het0.9_LD0.8.lfmm", col.names = F, row.names = F, sep = "\t", na = "9")
LEA::lfmm2geno("snp.imupted.maf0.05.geno0.2_het0.9_LD0.8.lfmm")

#PCA
pc <- LEA::pca('snp.imupted.maf0.05.geno0.2_het0.9_LD0.8.lfmm', scale = TRUE)
# Tracy-Widom tests on all eigenvalues to identify significant components
tw <- LEA::tracy.widom(pc)
print(tw$pvalues[1:10])
# Scree plot
windows()
plot(tw$eigenvalues, main = "Scree plot",
     ylab = "Eigenvalues", xlab = "PCs", t = "b")

# Running LFMM
tic()
mod.lfmm <- lfmm::lfmm_ridge(Y = as.matrix(genotype),
                             X = env_final[,4:28],
                             K = 3)
# P-values
pv <- lfmm::lfmm_test(Y = as.matrix(genotype), X = env_final[,4:28], lfmm = mod.lfmm,
                      calibrate = "gif")
toc()
pvalues <- pv$calibrated.pvalue
dim(pvalues)
head(pvalues)

# False discovery rate control with q-values  -----------------------------
head(pvalues)
my_qvalue <- function(x) {
  q <- qvalue::qvalue(x)
  q <- q$qvalues
  return(q)
}
qvalues <- apply(pvalues, 2, my_qvalue)
head(qvalues)
# Let's add SNP information (name, chromosome, physical position)  -- use the map file --
map <- read.table('289snp.imupted.maf0.05.geno0.2_het0.9_LD0.8.map')
qvalues <- as.data.frame(qvalues)
qvalues$SNP <- as.character(map$V2)
qvalues$CHR <- map$V1
qvalues$BP <- map$V4
head(qvalues)

# Let's decide a q-value cut-off of 1%
qcut <- 0.05
# LFMM results
lfmm_qvalcut <- function (qvalues = NULL, cutoff = NULL) {
  
  qvalues <- qvalues
  cutoff <- cutoff
  res <- list()
  e <- qvalues[, 1:(ncol(qvalues)-3)]
  s <- qvalues[, (ncol(qvalues)-2):ncol(qvalues)]
  
  for (i in 1:ncol(e)) {
    j <- which(e[, i] < cutoff)
    
    if (length(j) == 0) {
      t <- data.frame("ENV" = colnames(e)[i], 
                      "qval" = NA, 
                      "SNP" = NA, 
                      "CHR" = NA, 
                      "BP" = NA)
      res[[i]] <- t
    } else {
      t <- data.frame("ENV" = colnames(e)[i], 
                      "qval" = e[j, i], 
                      "SNP" = s[j, 1], 
                      "CHR" = s[j, 2], 
                      "BP" = s[j, 3])
      res[[i]] <- t
    }
  }
  res <- do.call(rbind.data.frame, res)
  res <- na.omit(res)
  return(res)
}
lfmm.res <- lfmm_qvalcut(qvalues = qvalues, cutoff = qcut)
print(lfmm.res)
