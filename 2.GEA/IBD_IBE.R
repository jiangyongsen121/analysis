#The code is sourced from "Genomic insights into local adaptation and future climate-induced vulnerability of a keystone forest tree in East Asia"

library(vegan)
library(geosphere)
env=read.xlsx('env.xlsx',sheet = "Sheet 1",rowNames = T)
fst_neutral=read.table('neutral_p_dis.mat',header=F)
fst_dro=read.table('Adrought_SNP_p_dis.mat',header=F)
row.names(fst_neutral) <- fst_neutral$V1
fst_neutral <- fst_neutral[,-1]
colnames(fst_neutral) <- row.names(fst_neutral)
row.names(fst_dro) <- fst_dro$V1
fst_dro <- fst_dro[,-1]
colnames(fst_dro) <- row.names(fst_dro)
env=env[row.names(fst_dro),]

library(dplyr)
env$ID=row.names(env)
env_site <- env %>%
  distinct(Longitude, Latitude, Altitude, .keep_all = TRUE)
env_scaled <- env_site %>%
  dplyr::select(-ID, -Longitude, -Latitude, -Altitude) %>%
  scale() %>%
  as.data.frame()
env_site_scaled <- cbind(env_site[ , c("Longitude","Latitude","Altitude")], env_scaled)
env_final <- env %>%
  dplyr::select(ID, Longitude, Latitude, Altitude) %>%
  left_join(env_site_scaled, by = c("Longitude","Latitude","Altitude"))
row.names(env_final) <- env_final[,1] 
env_final <- env_final[,-1]

### env distance
ENV<- env_final[,4:24]
env_dist=vegdist(ENV,method="euclidean",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)
write.table(as.matrix(env_dist),file="env_dist.txt",quote=F,sep="\t",row.names = T,col.names = T)

#### geographical distance
dist=as.data.frame(env_final[,c(1:2)])
muer.dists = distm(dist, fun=distVincentyEllipsoid)
rownames(muer.dists)=colnames(muer.dists) =row.names(env_final)
write.table(as.matrix(muer.dists),file="geo_dist.txt",quote=F,sep="\t",row.names = T,col.names = T)

#### geo-env dist pearson
env_dist=as.matrix(env_dist)
mantel(muer.dists,env_dist,method="pearson",permutations=999)

##### IBD
mantel(fst_dro,muer.dists,method="pearson",permutations=999)

##### IBE
mantel(fst_dro,env_dist,method="pearson",permutations=999)

##### IBE partial
mantel.partial(fst_dro,env_dist,muer.dists,method="pearson",permutations=999)

##### IBD partial
mantel.partial(fst_dro,muer.dists,env_dist,method="pearson",permutations=999)

###############  plot ######################
library(ggplot2)
library(cowplot)
trans <- function(raw_data){
  out_data=data.frame(raw_data[,1])
  colnames(out_data)="value"
  for (i in 2:108){
    temp=data.frame(raw_data[i:108,i])
    colnames(temp)="value"
    out_data=rbind(out_data,temp)
  }
  out_data=na.omit(out_data)
  return(out_data)
}
plot_data=as.data.frame(matrix(nrow=5886,ncol=0))
plot_data$geo_dist=trans(muer.dists)$value
plot_data$env_dist=trans(env_dist)$value
plot_data$fst_neu=trans(fst_neutral)$value
plot_data$fst_dro=trans(fst_dro)$value
colnames(plot_data)=c("geo_dist","env_dist","fstneu","fstdro")
write.csv(plot_data,file="NW_plot_data.csv",quote=F,row.names = F)

p1=ggplot(plot_data)+
  
  geom_smooth(aes(x=geo_dist, y=fstneu),alpha=0.7,formula = y ~ x, method = lm,se=T,level=0.95,color="#464647", fill="#DCE2E6",size = 1.5,fullrange = F) +
  geom_smooth(aes(x=geo_dist, y=fstdro),alpha=0.7,formula = y ~ x, method = lm,se=T,level=0.95,color="#CB2552",fill="#F5E0CA",fullrange = F, size = 1.5) +
  labs(x = "Geographical Distance",y = expression(italic(F)[italic(ST)]/(1-italic(F)[italic(ST)])),size = 5.5)+
  panel_border(color = "black", size = 0.6, linetype = 1, remove = FALSE)+
  theme_bw()+
  theme(text=element_text(family="sans"),
        axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
        axis.text.x=element_text(size=12,colour = "black"),
        axis.text.y=element_text(size=12,colour = "black"), 
        plot.title = element_text(
          size = 15L,
          hjust = 0
        ),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        panel.background=element_rect(fill="white"),
        plot.background = element_rect(fill = "white"),
        #axis.line.x=element_line(colour="black"),
        #axis.line.y=element_line(colour="black"),
        #panel.border=element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm"))

p2=ggplot(plot_data)+
 
  geom_smooth(aes(x=env_dist, y=fstneu),alpha=0.7,formula = y ~ x, method = lm,se=T,level=0.95,color="#464647", fill="#DCE2E6",size = 1.5,fullrange = F) +
  geom_smooth(aes(x=env_dist, y=fstdro),alpha=0.7,formula = y ~ x, method = lm,se=T,level=0.95,color="#CB2552",fill="#F5E0CA",fullrange = F, size = 1.5) +
  labs(x = "Environment Distance",y = expression(italic(F)[italic(ST)]/(1-italic(F)[italic(ST)])),size = 5.5)+
  panel_border(color = "black", size = 0.6, linetype = 1, remove = FALSE)+
  theme_bw()+
  theme(text=element_text(family="sans"),
        axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
        axis.text.x=element_text(size=12,colour = "black"),
        axis.text.y=element_text(size=12,colour = "black"), 
        plot.title = element_text(
          size = 15L,
          hjust = 0
        ),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        panel.background=element_rect(fill="white"),
        plot.background = element_rect(fill = "white"),
        #axis.line.x=element_line(colour="black"),
        #axis.line.y=element_line(colour="black"),
        #panel.border=element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm"))

all=plot_grid(p1,p2,align ="v",labels=c("c","d"),label_size = 20,label_fontfamily = "sans",label_fontface = 1,ncol=1)
ggsave(all,file="NW_IBD_IBE_nodot.pdf",width=5.5,height=8)

