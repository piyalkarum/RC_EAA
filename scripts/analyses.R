
# COMPLETE ANALYES : KARUNARATHNE ET AL. 2023 ####
#...............................................................................
# HYBRIDIZATION MEDIATED RANGE EXPANSION AND CLIMATE RESILIENCE IN TWO KEYSTONE
# SPECIES OF BOREAL FORESTS
# ..............................................................................

# This script contains all the major analyses presented in the publication
# "Hybridization mediated niche expansion and climate change resiliance in two
# keystone stree species of boreal forests."

## 1. RDA ####

library(pegas)
library(ggplot2)
library(raster)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(vegan)
library(qvalue)
library(ggVennDiagram)
library(corrplot)
library(rgeos)
library(adespatial)
library(geodist)
library(robust)

out_pd<-"./output"
dt<-gsub("-","_",substr(Sys.time(),1,10))
source("./src/rdadapt.R")
source("./src/adaptive_index.R")

AF<-readRDS("./data/allele.freq_pop.rds")
AllF<-t(AF)
freq_mean <- colMeans(AllF)
AllF <- AllF[,-which(freq_mean>=0.95 | freq_mean<=0.05)] # remove MAF < 5%
#environmental data (WC)
Env.t<-read.table("./data/env_wc_noGbif.txt",h=T)
# get geo distance svd
gd<-svd(as.dist(geodist(Env.t[,1:2],measure = "geodesic")))$d
## Standardization of the variables
Env<-data.frame(Env.t[,-c(1:5)],geod=gd)
Env<-Env[match(colnames(AF),Env.t$Source),]
Env <- scale(Env, center=TRUE, scale=TRUE)
row.names(Env) <- c(Env.t$Source)
Variables <- data.frame(Env[,-ncol(Env)])
#Variable selection
ev1<-Variables
ccr<-cor(ev1)
corrplot::corrplot(ccr)
diag(ccr)<-0
thr<-range(abs(ccr))[2]
while(thr>0.6){
  ccr<-cor(ev1)
  diag(ccr)<-0
  thr<-range(abs(ccr))[2]
  ev1<-ev1[,-(which.max(colSums(abs(ccr))))]
}
lyrs<-colnames(ev1)
#Variables<-Variables[,lyrs]
#bio8, bio9, bio18

evars<-lyrs
eq<-paste0(evars,collapse = "+")
RDAfull<-rda(as.formula(paste0("AllF ~",eq)),Variables)
(vf<-vif.cca(RDAfull))
#keep only the vars with VIF < 5
while(max(vf)>1.1){
  evars<-names(vf)
  evars<-evars[-which.max(vf)]
  eq<-paste0(evars,collapse = "+")
  RDAfull<-rda(as.formula(paste0("AllF ~",eq)),Variables)
  (vf<-vif.cca(RDAfull))
}

### 1.1 Partial RDA ######################
pRDAfull <- rda(as.formula(paste0("AllF ~ geod +",eq)),  data.frame(Env))
RsquareAdj(pRDAfull)
#anova(pRDAfull)
RDA_env <- rda(as.formula(paste0("AllF ~",eq,"+ Condition(geod)")),  data.frame(Env))
RsquareAdj(RDA_env)
#genome scan
rdadapt_env<-rdadapt(RDA_env, 2)
thres_env <- 0.05/length(rdadapt_env$p.values)## P-values threshold after Bonferroni correction
## Identifying the loci that are below the p-value threshold
outliers <- data.frame(Loci = colnames(AllF)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(AllF)[which(rdadapt_env$p.values<thres_env)], split = "\\."), function(x) x[1])))
outliers <- outliers[order(outliers$contig, outliers$p.value),]## Top hit outlier per contig
outliers_rdadapt_env <- as.character(outliers$Loci[!duplicated(outliers$contig)])
locus_scores <- scores(RDA_env, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci

TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "All outliers"
TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp")) # pull the biplot scores

ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", 4, 2)) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

#### p-value plot ################
snp.sort<-read.table("./data/snp_list_3_ordered.txt",h=F)
cl2<-hcl.colors(n=length(unique(snp.sort$V4)),palette = "Emrld")
allsnp<-data.frame(snp.name=colnames(AllF),rdadapt_env)
allsnp.sorted<-allsnp[match(paste0(snp.sort[,1],".",snp.sort[,2]),allsnp$snp.name),]
plot(-log10(allsnp.sorted$p.values)~rownames(snp.sort),pch=19,cex=.2,col=cl2[snp.sort$V4],xlab="SNP",ylab="log10(RDA - p.values)")
points(-log10(allsnp.sorted$p.values)~rownames(snp.sort),pch=19,cex=.2,col=ifelse(allsnp.sorted$p.values<thres_env,2,NA))
#####################################

#### QQ plot of p-values ####
Outliers <- rep("Neutral", length(colnames(AllF)))
Outliers[colnames(AllF)%in%outliers$Loci] <- "All outliers"
Outliers[colnames(AllF)%in%outliers_rdadapt_env] <- "Top outliers"
Outliers <- factor(Outliers, levels = c("Neutral", "All outliers", "Top outliers"))
write.table(outliers,paste0(out_pd,"PA_RDA_outliers.txt"),row.names=F,quote=F,sep="\t")

TAB_manhatan <- data.frame(pos = 1:length(colnames(AllF)), pvalues = rdadapt_env$p.values,Outliers = Outliers)
TAB_manhatan<-TAB_manhatan[order(TAB_manhatan$pvalues,decreasing = T),2:3]
TAB_manhatan$pos<-seq(0,1,length.out=length(TAB_manhatan$pvalues))

ggplot(data = TAB_manhatan) +
  geom_point(aes(x=pos, y=-log10(pvalues), col = Outliers), size=1.4) +
  scale_color_manual(values = c("gray90", 4, 2)) +
  xlab("Loci") + ylab("-log10(p.values)") +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  facet_wrap(~"QQ plot of p-values", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(), legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

#### Running a simple RDA model ####
RDA_env_unconstrained <- rda(as.formula(paste0("AllF ~ ",eq)),  Variables)
rdadapt_env_unconstrained <- rdadapt(RDA_env_unconstrained, 2) ## Running the rdadapt function
thres_env <- 0.05/length(rdadapt_env_unconstrained$p.values)## Setting the p-value threshold

#### Identifying the outliers for the simple RDA ####
outliers_unconstrained <- data.frame(Loci = colnames(AllF)[which(rdadapt_env_unconstrained$p.values<thres_env)], p.value = rdadapt_env_unconstrained$p.values[which(rdadapt_env_unconstrained$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(AllF)[which(rdadapt_env_unconstrained$p.values<thres_env)], split = "\\."), function(x) x[1])))
outliers_unconstrained <- outliers_unconstrained[order(outliers_unconstrained$contig, outliers_unconstrained$p.value),]
outliers_rdadapt_env_unconstrained <- as.character(outliers_unconstrained$Loci[!duplicated(outliers_unconstrained$contig)])
## For all the outliers
list_outliers_RDA_all <- list(RDA_constrained = as.character(outliers$Loci), RDA_unconstrained = as.character(outliers_unconstrained$Loci))
## Only for the top hit locus per contig
list_outliers_RDA_top <- list(RDA_constrained = outliers_rdadapt_env, RDA_unconstrained = outliers_rdadapt_env_unconstrained)
## common outliers
common_outliers_RDA_top <- Reduce(intersect, list_outliers_RDA_top)
## Adaptively enriched RDA
RDA_outliers <- rda(as.formula(paste0("AllF[,common_outliers_RDA_top] ~", eq)),  Variables)
ss<-summary(RDA_outliers)$concont$importance #take [2] and [5]
#### RDA biplot ####
TAB_loci <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="bp"))
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*3, y=RDA2*3), colour = 2, size = 2, alpha = 0.8) + #"#F9A242FF"
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 3.5, family = "Times", colour=4) +
  xlab(paste0("RDA 1 (",round(ss[2]*100,0),"%)")) + ylab(paste0("RDA 2 (",round(ss[5]*100,0),"%)")) +
  facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

### 1.2 adaptive index predictions ###################
#### Data ####
# download the WorldClim 2.1 bioclimatic dataset at 2.5 AS resolution
# at <https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_2-5m_esri.zip>
# and save it to the data directory
efl<-list.files("./data/wc1.4",pattern = ".tif",full.names = T)
en<-stack(efl)
names(en)<-colnames(Variables)
ras<-dropLayer(en,sort(which(!names(en)%in%evars)))
ran<-extent(-10,140,35,90)
ras<-crop(ras,ran)
rng<-shapefile("./data/PAandPO_combined_map.shp")
#rng<-rng[rng$LEVEL1_COD==1 | rng$LEVEL1_COD==3,]
#rng<-crop(rng,ran)
plot(rng,border="grey90")
points(Env.t[,c(3,2)],pch=19,cex=0.5,col="blue")

admin <- ne_countries(scale = "medium", returnclass = "sf")
scale_env <- attr(Env, 'scaled:scale')
names(scale_env)<-colnames(Env)
center_env <- attr(Env, 'scaled:center')
names(center_env)<-colnames(Env)
res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 2, env_pres = ras,range=rng, method = "loadings", scale_env = scale_env, center_env = center_env)
writeRaster(res_RDA_proj_current$RDA1,paste0(out_pd,"RDA1_prediction.tif"))

RDA_proj <- rasterToPoints(res_RDA_proj_current$RDA1)
RDA_proj[,3] <- (RDA_proj[,3]-min(RDA_proj[,3]))/(max(RDA_proj[,3])-min(RDA_proj[,3]))

#### Adaptive genetic turnover projected across range for RDA1 and RDA2 indexes ####
TAB_RDA <- as.data.frame(RDA_proj)
colnames(TAB_RDA)[3] <- "value"
#TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))
saveRDS(TAB_RDA,paste0(out_pd,"RDA1_table.rds"))

#### Extract RDA values for comparison ####
RDAvals<-data.frame(Env.t[,1:5])
for(i in seq_along(res_RDA_proj_current)){
  tm<-res_RDA_proj_current[[i]]
  tm<-extract(tm,Env.t[,1:2])
  RDAvals<-cbind(RDAvals,tm)
}
colnames(RDAvals)[6:7]<-c("RDA1","RDA2")

#### RDA plot ####
pdf(paste0(out_pd,"/RDA1_.pdf"))
ggplot(data = TAB_RDA) +
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) +
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "C", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
  geom_sf(data = admin, fill=NA, size=0.1) +
  coord_sf(xlim = c(-10, 111), ylim = c(35, 75), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))
dev.off()

### for RDA2
RDA_proj2 <- rasterToPoints(res_RDA_proj_current$RDA2)
RDA_proj2[,3] <- (RDA_proj2[,3]-min(RDA_proj2[,3]))/(max(RDA_proj2[,3])-min(RDA_proj2[,3]))

## Adaptive genetic turnover projected across range for RDA1 and RDA2 indexes
TAB_RDA2 <- as.data.frame(RDA_proj2)
colnames(TAB_RDA2)[3] <- "value"
#TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))
saveRDS(TAB_RDA2,paste0(out_pd,"RDA2_table.rds"))

pdf(paste0(out_pd,"/RDA2_.pdf"))
ggplot(data = TAB_RDA2) +
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) +
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "C", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
  geom_sf(data = admin, fill=NA, size=0.1) +
  coord_sf(xlim = c(-10, 111), ylim = c(35, 75), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))
dev.off()


## 2. Environmental clusters ####
### 2.1 Clustering based on adaptively enriched RDA space ####
library(terra)
library(factoextra)

setwd("/Users/piyalkaru/Desktop/UPP/RUS_CLIN/github/RC_EAA")
env<-read.table("./data/env.2.txt",h=T)
hstat <- data.table::fread("./data/env_wc_noGbif2.txt", h = TRUE)[,3:5]
aes<-rast("./output/RDA1_prediction.tif")
wm<-vect("./data/level3/level3.shp")
rng<-ext(-10,98,35,75)
wm2<-crop(wm,rng)
# remove multicolinearity
ev1<-scale(env[,-c(1:5)])
ccr<-cor(ev1)
diag(ccr)<-0
thr<-range(abs(ccr))[2]
while(thr>0.5){
  ccr<-cor(ev1)
  diag(ccr)<-0
  thr<-range(abs(ccr))[2]
  ev1<-ev1[,-(which.max(colSums(abs(ccr))))]
}
lyrs<-colnames(ev1)
s.env<-ev1[,lyrs]
row.names(s.env)<-env$Source

# extract RDA values
crd<-env[,c("Long","Lat")]
env2<-data.frame(env[,1:5],extract(aes,crd))

hstat<-hstat[match(env2$Source,hstat$Source),]
hybs<-data.frame(hstat[hstat$stat=="pt_hy",])
tt<-data.frame(s.env,eas=env2$RDA_pres_1)
tt<-tt[match(hybs$Source,rownames(tt)),] ## hybrids only

rr<-data.frame(lon=env2$Long,lat=env2$Lat,s.env,esa=env2$RDA_pres_1)
rr2<-rr[-match(hybs$Source,rownames(rr)),]

#### Finding optimal clusters ####
fviz_nbclust(rr[,-c(1,2)],kmeans,method = "wss",k.max = 20)
fviz_nbclust(rr2[,-c(1,2)],kmeans,method="wss")

set.seed(999)
cent<-5
cl<-hcl.colors(cent,"berlin")
cl2<-hcl.colors(n=10,palette = "red-blue")
km2<-kmeans(rr2[,-c(1,2)],centers = cent)
par(mai=c(0,0,0,0))
plot(wm2,col="grey90",border=F,axes=F,buffer=F,bty="n")
points(rr2[,c(1,2)],pch=19,cex=1.8,col=cl[km2$cluster])

#### Hybrid zone clusters ####
km3<-kmeans(tt,centers = 3)
points(vect(hybs,geom=c("Long","Lat"),crs=crs(wm2)),pch=1,cex=1.8,col=rev(cl)[km3$cluster])
text(env[,c("Long","Lat")],labels=env$Source,cex=.6,adj=c(1,0))

saveRDS(list(nonhyb=km2,hyb=km3),"./output/kmeans_clustering2.rds")


## 3. Niche optima and centroid ####
### 3.1. Niche overlap (per Karunarathne et al. 2018) --------------------------
env<-read.table("./data/env.2.txt",h=T)
#### variable selection #####
# This additional step is necessary since the significant variables for the combined species is different from individual species
ev<-scale(env[,-c(1:5)])
ev1<-ev
ccr<-cor(ev1)
diag(ccr)<-0
thr<-range(abs(ccr))[2]
while(thr>0.5){
  ccr<-cor(ev1)
  diag(ccr)<-0
  thr<-range(abs(ccr))[2]
  ev1<-ev1[,-(which.max(colSums(abs(ccr))))]
}
lyrs<-colnames(ev1)

dat<-env
padat<-rowSums(log10(env[env$Species=="PA",lyrs]))
podat<-rowSums(log10(env[env$Species=="PO",lyrs]))
## kernal density
paden<-density(padat,adjust=1)
poden<-density(podat,adjust=1)
drn<-range(c(range(paden$x),range(poden$x)))
## plotting
cld<-rCNV:::makeTransparent(c(1,2),alpha = .6)
plot(0,xlim=drn,ylim=c(0,6),xlab="Environmental gradient",ylab="Probability density of occurrence",xaxt="n",bty="n")
polygon(paden,col=cld[1])
polygon(poden,col=cld[2])
legend("topleft",pch=19,col=cld,legend=c("P.abies","P.obovata"),bty="n")

### 3.2 Niche Centroid ----------------------------------------------
library(factoextra)
library(vegan)
library(ade4)

setwd("/Users/piyalkaru/Desktop/UPP/RUS_CLIN/github/RC_EAA")
env<-read.table("./data/env.2.txt",h=T)
aes<-rast("./output/RDA1_prediction.tif")
papo<-vect("./data/PAandPO_combined_map.shp")
wm<-vect("./data/level3/level3.shp")
hstat <- data.frame(data.table::fread("./data/env_wc_noGbif2.txt", h = TRUE))[,3:5]
hstat<-hstat[match(env$Source,hstat$Source),]
#flag hybrids
h.pop<-hstat[hstat$stat=="pt_hy",3]
crd<-env[,c("Long","Lat")]
env2<-data.frame(env[,1:5],env[,lyrs],extract(aes,crd))

# WO hybrids
env2<-data.frame(env2[!env2$Source %in% h.pop,])
pc.dat<-env2[,c(6:8,10)]
rownames(pc.dat)<-env2$Source

nch.pca<-dudi.pca(pc.dat,scannf = F, nf=6)
cl<-hcl.colors(10,"batlow")
fviz_pca_ind(nch.pca,col.ind = "cos2",gradient.cols=cl,repel=T)
fviz_pca_ind(nch.pca,col.ind = factor(env2$Species),palette = c(2,4),addEllipses = T,repel = T)
s.class(nch.pca$li,fac=factor(all.nch$Species),col=c(2,4))

# With hybrids
env2<-data.frame(env[,1:5],env[,lyrs],extract(aes,crd))
pc.dat<-env2[,c(6:8,10)]
rownames(pc.dat)<-env2$Source

nch.pca<-dudi.pca(pc.dat,scannf = F, nf=6)
fviz_pca_ind(nch.pca,col.ind = "cos2",gradient.cols=cl,repel=T)
stats<-env2$Species
stats[env0$stat=="pt_hy"]<-"pt_hy"
#PA PO niche only
fviz_pca_ind(nch.pca,col.ind = factor(env2$Species),palette = c("#1FA7EA","#EE8900"),addEllipses = T,repel = T)
# add hybrid niche as a third ellipse
fviz_pca_ind(nch.pca,col.ind = factor(stats),palette = c("#1FA7EA","#EE8900","red"),addEllipses = T,repel = T)
s.class(nch.pca$li,fac=factor(env$Species),col=c(2,4))

## 4. SDM with dismo ####
library(raster)
library(dismo)

efl<-list.files("./data/wc1.4/curr/tif",pattern = ".bil",full.names = T)
crd<-read.table("./data/coords_w_hyb_and_gbif_oct28.txt",h=T)
wm<-shapefile("./data/level3/level3.shp")
papo<-shapefile("./data/PAandPO_combined_map.shp")

pres<-crd[crd$stat=="pt_nonhy",]
rng1<-wm[wm$LEVEL2_COD%in%c(10:14,30),]
rng1<-rng1[!rng1$LEVEL3_COD%in%c("FOR","ICE","IRE","GRB","SVA","DEN"),]
## extract bioclim data for populations
bio<-stack(efl)
names(bio)<-paste0("Bio",1:19)
env<-data.frame(pres,extract(bio,pres[,1:2])[,-1])

### 4.1 for P. abies ####
env1<-env[env$Species=="PA",]
env1<-na.omit(env1)
#### Check for autocorrelation ####
ev1<-scale(env1[,-c(1:5)])
ccr<-cor(ev1)
diag(ccr)<-0
thr<-range(abs(ccr))[2]
while(thr>0.6){
  ccr<-cor(ev1)
  diag(ccr)<-0
  thr<-range(abs(ccr))[2]
  ev1<-ev1[,-(which.max(colSums(abs(ccr))))]
}
lyrs<-colnames(ev1)

rng<-extent(-10,130,35,90)
enrast<-bio[[lyrs]]
enrast<-mask(enrast,rng1)
enrast<-crop(enrast,rng)
set.seed(123)

pa<-env1[env1$Species=="PA",1:2]
mx<-maxent(enrast,pa,nbg=nrow(env1)*10)
pt<-predict(enrast,mx)
writeRaster(pt,paste0("./output/",unique(env1$Species),"_dismoMaxEnt_",substr(Sys.time(),1,10),".tif"))


### 4.2 for P. obovata ####
library(raster)
library(dismo)

efl<-list.files("./data/wc1.4/curr/tif",pattern = ".bil",full.names = T)
crd<-read.table("./data/coords_w_hyb_and_gbif_oct28.txt",h=T)
wm<-shapefile("./data/level3/level3.shp")
papo<-shapefile("./data/PAandPO_combined_map.shp")

pres<-crd[crd$stat=="pt_nonhy",]
rng1<-wm[wm$LEVEL2_COD%in%c(10:14,30),]
rng1<-rng1[!rng1$LEVEL3_COD%in%c("FOR","ICE","IRE","GRB","SVA","DEN"),]
## extract bioclim data for populations
bio<-stack(efl)
names(bio)<-paste0("Bio",1:19)
env<-data.frame(pres,extract(bio,pres[,1:2])[,-1])

env1<-env[env$Species=="PO",]
env1<-na.omit(env1)
#### check for autocorrelation ####
ev1<-scale(env1[,-c(1:5)])
ccr<-cor(ev1)
diag(ccr)<-0
thr<-range(abs(ccr))[2]
while(thr>0.6){
  ccr<-cor(ev1)
  diag(ccr)<-0
  thr<-range(abs(ccr))[2]
  ev1<-ev1[,-(which.max(colSums(abs(ccr))))]
}
lyrs<-colnames(ev1)

rng<-extent(-10,130,35,90)
enrast<-bio[[lyrs]]
enrast<-mask(enrast,rng1)
enrast<-crop(enrast,rng)
set.seed(123)

pa<-env1[env1$Species=="PO",1:2]
mx<-maxent(enrast,pa,nbg=nrow(env1)*20)
pt<-predict(enrast,mx)
writeRaster(pt,paste0("./output/",unique(env1$Species),"_dismoMaxEnt_",substr(Sys.time(),1,10),".tif"))

## 5 Paleo-SDM DISMO ####
### 5.1 P. abies ####
library(raster)
library(dismo)

efl<-list.files("./data/wc1.4/curr/tif",pattern = ".bil",full.names = T)
crd<-read.table("./data/coords_w_hyb_and_gbif_oct28.txt",h=T)
wm<-shapefile("./data/level3/level3.shp")
papo<-shapefile("./data/PAandPO_combined_map.shp")

pres<-crd[crd$stat=="pt_nonhy",]
bio<-stack(efl)
names(bio)<-paste0("Bio",1:19)
env<-data.frame(pres,extract(bio,pres[,1:2])[,-1])

env1<-env[env$Species=="PA",]
env1<-na.omit(env1)
#### check for autocorrelation ####
ev1<-scale(env1[,-c(1:5)])
ccr<-cor(ev1)
diag(ccr)<-0
thr<-range(abs(ccr))[2]
while(thr>0.6){
  ccr<-cor(ev1)
  diag(ccr)<-0
  thr<-range(abs(ccr))[2]
  ev1<-ev1[,-(which.max(colSums(abs(ccr))))]
}
lyrs<-colnames(ev1)

rng<-extent(-10,130,35,90)
enrast<-bio[[lyrs]]
enrast<-crop(enrast,rng)
set.seed(123)
pa<-env1[env1$Species=="PA",1:2]
mx<-maxent(enrast,pa,nbg=nrow(env1)*10)

saveRDS(mx,paste0("./output/",unique(env1$Species),"MaxEnt_model",substr(Sys.time(),1,10),".rds"))
#### 5.1.1 for LGM -------------------------------------------------------------
lgm<-list.files("./data/wc1.4/lgm",full.names = T,pattern=".tif")
lgm<-stack(lgm)
names(lgm)<-paste0("Bio",1:19)
lgm<-lgm[[lyrs]]
rng2<-extent(-15,170,20,90)
lgm<-crop(lgm,rng2)
lgm<-mask(lgm,wm)
pa.lgm<-predict(lgm,mx)

writeRaster(pa.lgm,paste0("./output/",unique(env1$Species),"_LGM_dismoMaxEnt_",substr(Sys.time(),1,10),".tif"))

#### 5.1.2 MidHolocene ---------------------------------------------------------
mh<-list.files("./data/wc1.4/midhol",full.names = T,pattern=".tif")
mh<-stack(mh)
names(mh)<-paste0("Bio",1:19)
mh<-mh[[lyrs]]
rng2<-extent(-15,170,20,90)
mh<-crop(mh,rng2)
mh<-mask(mh,wm)
pa.mh<-predict(mh,mx)
writeRaster(pa.mh,paste0("./output/",unique(env1$Species),"_MidHol_dismoMaxEnt_",substr(Sys.time(),1,10),".tif"))

#### 5.1.3 For Inter-glacial period (LIG) --------------------------------------
ig<-list.files("./data/wc1.4/intglac",full.names = T,pattern=".bil")
ig<-stack(ig)
names(ig)<-paste0("Bio",1:19)
ig<-ig[[lyrs]]
rng2<-extent(-15,170,20,90)
ig<-crop(ig,rng2)
ig<-mask(ig,wm)
pa.ig<-predict(ig,mx)
writeRaster(pa.ig,paste0("./output/",unique(env1$Species),"_intGl_dismoMaxEnt_",substr(Sys.time(),1,10),".tif"))

#### 5.1.4 Future --------------------------------------------------------------
fu<-list.files("./data/wc1.4/futur",full.names = T,pattern=".tif")
fu<-stack(fu)
names(fu)<-paste0("Bio",1:19)
fu<-fu[[lyrs]]
rng2<-extent(-15,170,20,90)
fu<-crop(fu,rng2)
fu<-mask(fu,wm)
pa.fu<-predict(fu,mx)
writeRaster(pa.fu,paste0("./output/",unique(env1$Species),"_Futr_dismoMaxEnt_",substr(Sys.time(),1,10),".tif"))

### 5.2 P. obovata ####
library(raster)
library(dismo)

efl<-list.files("./data/wc1.4/curr/tif",pattern = ".bil",full.names = T)
crd<-read.table("./data/coords_w_hyb_and_gbif_oct28.txt",h=T)
wm<-shapefile("./data/level3/level3.shp")
papo<-shapefile("./data/PAandPO_combined_map.shp")

pres<-crd[crd$stat=="pt_nonhy",]
bio<-stack(efl)
names(bio)<-paste0("Bio",1:19)
env<-data.frame(pres,extract(bio,pres[,1:2])[,-1])

env1<-env[env$Species=="PO",]
env1<-na.omit(env1)
## check for autocorrelation
ev1<-scale(env1[,-c(1:5)])
ccr<-cor(ev1)
diag(ccr)<-0
thr<-range(abs(ccr))[2]
while(thr>0.6){
  ccr<-cor(ev1)
  diag(ccr)<-0
  thr<-range(abs(ccr))[2]
  ev1<-ev1[,-(which.max(colSums(abs(ccr))))]
}
lyrs<-colnames(ev1)

rng<-extent(-10,130,35,90)
enrast<-bio[[lyrs]]
enrast<-crop(enrast,rng)
set.seed(123)
pa<-env1[env1$Species=="PO",1:2]
mx<-maxent(enrast,pa,nbg=nrow(env1)*10)
saveRDS(mx,paste0("./output/",unique(env1$Species),"MaxEnt_model",substr(Sys.time(),1,10),".rds"))

#### 5.2.1  for LGM ------------------------------------------------------------
lgm<-list.files("./data/wc1.4/lgm",full.names = T,pattern=".tif")
lgm<-stack(lgm)
names(lgm)<-paste0("Bio",1:19)
lgm<-lgm[[lyrs]]
rng2<-extent(-15,170,20,90)
lgm<-crop(lgm,rng2)
lgm<-mask(lgm,wm)
pa.lgm<-predict(lgm,mx)
writeRaster(pa.lgm,paste0("./output/",unique(env1$Species),"_LGM_dismoMaxEnt_",substr(Sys.time(),1,10),".tif"))

#### 5.2.2 MidHolocene ---------------------------------------------------------
mh<-list.files("./data/wc1.4/midhol",full.names = T,pattern=".tif")
mh<-stack(mh)
names(mh)<-paste0("Bio",1:19)
mh<-mh[[lyrs]]
rng2<-extent(-15,170,20,90)
mh<-crop(mh,rng2)
mh<-mask(mh,wm)
pa.mh<-predict(mh,mx)
writeRaster(pa.mh,paste0("./output/",unique(env1$Species),"_MidHol_dismoMaxEnt_",substr(Sys.time(),1,10),".tif"))

#### 5.2.3 For Inter-glacial period (LIG) --------------------------------------
ig<-list.files("./data/wc1.4/intglac",full.names = T,pattern=".bil")
ig<-stack(ig)
names(ig)<-paste0("Bio",1:19)
ig<-ig[[lyrs]]
rng2<-extent(-15,170,20,90)
ig<-crop(ig,rng2)
ig<-mask(ig,wm)
pa.ig<-predict(ig,mx)
writeRaster(pa.ig,paste0("./output/",unique(env1$Species),"_intGl_dismoMaxEnt_",substr(Sys.time(),1,10),".tif"))

#### 5.2.4 Future projection ---------------------------------------------------
fu<-list.files("./data/wc1.4/futur",full.names = T,pattern=".tif")
fu<-stack(fu)
names(fu)<-paste0("Bio",1:19)
fu<-fu[[lyrs]]
rng2<-extent(-15,170,20,90)
fu<-crop(fu,rng2)
fu<-mask(fu,wm)
pa.fu<-predict(fu,mx)
writeRaster(pa.fu,paste0("./output/",unique(env1$Species),"_Futr_dismoMaxEnt_",substr(Sys.time(),1,10),".tif"))


### 5.3 Plotting locally ####
# current
pa.cr<-rast("./output/PA_dismoMaxEnt.tif")
po.cr<-rast("./output/PO_dismoMaxEnt.tif")
# lgm
pa.lgm<-rast("./output/PA_LGM_dismoMaxEnt.tif")
po.lgm<-rast("./output/PO_LGM_dismoMaxEnt.tif")
# mid-holo
pa.mh<-rast("./output/PA_MidHol_dismoMaxEnt.tif")
po.mh<-rast("./output/PO_MidHol_dismoMaxEnt.tif")
# lig
pa.ig<-rast("./output/PA_intGl_dismoMaxEnt.tif")
po.ig<-rast("./output/PO_intGl_dismoMaxEnt.tif")
# future
pa.fu<-rast("./output/PA_Futr_dismoMaxEnt.tif")
po.fu<-rast("./output/PO_Futr_dismoMaxEnt.tif")
# worldmap
wm<-vect("./data/ne_10m_land/ne_10m_land.shp")

par(mfrow=c(2,5))
plot(pa.fu,legend=F,main="P.abies future")
plot(pa.cr,legend=F,main="P.abies current")
plot(pa.mh,legend=F,main="P.abies mid-Holocene")
plot(pa.lgm,legend=F,main="P.abies LGM")
plot(pa.ig,legend=F,main="P.abies Int.Glac")

plot(po.fu,legend=F,main="P.ovovata future")
plot(po.cr,legend=F,main="P.obovata current")
plot(po.mh,legend=F,main="P.obovata mid-Holocene")
plot(po.lgm,legend=F,main="P.obovata LGM")
plot(po.ig,legend=F,main="P.ovovata Int.Glac")


## 6. RONA calculation ---------------------------------------------------------
gen<-readRDS("./data/allele.freq_pop_Sept6_22.rds")
env<-read.table("./data/env_wc_noGbif2.txt",h=T)
fu<-read.table("./data/future_env.txt",h=T)
env<-env[match(env$Source,colnames(gen)),]
fu<-fu[match(fu$Source,colnames(gen)),]

### remove multicolinearity ####
ev1<-scale(env[,-c(1:5)])
ccr<-cor(ev1)
diag(ccr)<-0
thr<-range(abs(ccr))[2]
while(thr>0.6){
  ccr<-cor(ev1)
  diag(ccr)<-0
  thr<-range(abs(ccr))[2]
  ev1<-ev1[,-(which.max(colSums(abs(ccr))))]
}
lyrs<-colnames(ev1)

lyrs<-c("Bio8" , "Bio9" , "Bio18")
present<-env[,lyrs]
future<-fu[,lyrs]

rona<-matrix(NA,nrow=nrow(gen),ncol=ncol(gen))
stats<-matrix(NA,nrow=nrow(gen),ncol=5)
evlist<-list()
for(i in 1:ncol(present)){
  # progress bar
  pb<-txtProgressBar(min=0,max=ncol(present),style=3,width=50,char="=")
  setTxtProgressBar(pb,i)
  for(k in 1:nrow(gen)){
    x<-present[,i]
    y<-gen[k,]
    mod<-lm(y~x)
    msum<-summary(mod)
    rona[k,]<-abs((msum$coefficients[2] * future[,i] + msum$coefficients[1]) - gen[k,])
    stats[k,]<-c(msum$fstatistic[1],msum$r.squared,msum$adj.r.squared,msum$sigma,msum$coefficients[8])
  }
  colnames(rona)<-colnames(gen);rownames(rona)<-rownames(gen)
  colnames(stats)<-c("stat","r.squared","adj.r.squared","sigma","pval")
  rownames(stats)<-rownames(gen)
  stats<-data.frame(stats)
  ## average RONA
  weighted<-apply(rona,2,weighted.mean,w=stats$r.squared)
  unweighted<-colMeans(rona)
  avg.rona<-data.frame(weighted,unweighted)
  ####
  evlist[[i]]<-list(RONA=rona,stats=stats,meanRONA=avg.rona)
}
close(pb)
names(evlist)<-colnames(present)

saveRDS(evlist,paste0("./output/RONA_out_",substr(Sys.time(),1,10),".rds"))

### Check locally ####
rona<-readRDS("./output/RONA/RONA_out.rds")

for(i in seq_along(rona)){
  tmp1<-rona[[i]]$RONA
  tmp2<-rona[[i]]$stats
  
  top.rona<-tmp1[tmp2$pval<0.05,]
  top.stat<-tmp2[tmp2$pval<0.05,]
  
  weighted<-apply(top.rona,2,weighted.mean,w=top.stat$r.squared)
  unweighted<-colMeans(top.rona)
  se<-apply(top.rona,2,function(x){sd(x)/sqrt(length(x))})
  top.avg.rona<-data.frame(weighted,unweighted,se)
  write.table(top.avg.rona,paste0("output/","top.avg.RONA",substr(Sys.time(),1,10),".",names(rona)[i],".txt"),quote=F,sep="\t")
  
}

#### Partial response plots of top two alleles ####
# need coefficients from all lms and find the top outlier allels and get the partial response.
# the new one with coefficients from sept23
rona<-readRDS("output/RONA_out.rds")
gen<-readRDS("data/allele.freq_pop.rds")
env<-read.table("data/env_wc_noGbif2.txt",h=T)
env<-env[match(env$Source,colnames(gen)),]
lyrs<-c("Bio8" , "Bio9" , "Bio18")
present<-env[,lyrs]
# get the top 10 loci for each variable, do lm on the allele frequency and plot the response curve
set.seed(999)
mls_bio<-list()
for(i in seq_along(rona)){
  tm<-rona[[i]]$stats
  tm<-tm[order(tm$coef,decreasing = T),]
  top10<-rownames(tm)[1:10]
  tgen<-gen[match(top10,rownames(gen)),]
  x<-present[,i]
  y<-seq(min(tgen),max(tgen),length.out=length(x))
  plot(x,y,type="n",main=lyrs[i],xlab="Precipitation (mm)",ylab="Allele frequency")
  cl<-hcl.colors(nrow(tgen),palette = "dark3")
  new_x <- data.frame(x = seq(min(x), max(x), length.out =1000))
  for(j in 1:nrow(tgen)){
    y=tgen[j,]
    mod<-glm(y~x)
    new_y<-predict(mod,newdata=new_x)
    lines(new_x$x,new_y,col=cl[j])
  }
}


### RONA map plots ####
library(terra)
source("./scripts/raster_functions.R")
#### Extract weighted and non-weighted RONA ------------------------------------
wm<-vect("./data/level3/level3.shp")
papo<-vect("./data/PAandPO_combined_map.shp")
fl<-list.files("./output/rona_avg",full.names = T,pattern = ".txt")
crd<-read.table("./data/coords_w_hyb_and_gbif_oct28.txt",h=T)

for(k in seq_along(fl)){
  pi<-read.table(fl[k],h=T)
  pi.n<-pi
  bn<-stringr::str_split(fl[1],"\\.")[[1]][4]
  tm<-crd[crd$Source%in%rownames(pi),]
  pi<-cbind(tm[,1:2],pi)
  rownames(pi)<-rownames(pi.n)
  colnames(pi)[1:2]<-c("lon","lat")
  pr<-rast_index(pi,resolution = .25,extent=20,field = "weighted")
  wo<-dup_mean(pr,pi)
  pr2<-rast_index(wo,resolution = .25,extent=20,field = "weighted")
  
  ### interpolate (fill gaps) with nearest neighbor sliding window ######
  mt<-matrix(values(pr2),nrow = dim(pr2)[1],ncol = dim(pr2)[2],byrow = T)
  wind<-50 # window size
  width<-3 # average cells
  for(i in 1:wind){
    # progress bar
    pb <- txtProgressBar(min = 0, max = wind, style = 3, width = 50, char = "=")
    setTxtProgressBar(pb, i)
    #
    # rows
    ht<-apply(mt,1,function(x){
      mps<-NULL
      mpsf<-NULL
      for(j in 1:length(x)){
        #left
        if(j+width<=length(x)){
          if(is.na(x[j])){
            mps[j]<-max(x[j:j+1],na.rm = T)
          } else {mps[j]<-mean(x[j:(j+width)],na.rm=T)}
        }else{
          if(is.na(x[j])){
            mps[j]<-max(x[j:j-1],na.rm = T)
          } else {mps[j]<-mean(x[j:(j-width)],na.rm=T)}
        }
        #right
        if(j-width>0){
          if(is.na(x[j])){
            mpsf[j]<-max(x[j:j-1],na.rm = T)
          } else {mpsf[j]<-mean(x[j:(j-width)],na.rm=T)}
        }else{
          if(is.na(x[j])){
            mpsf[j]<-max(x[j:j+1],na.rm = T)
          } else {mpsf[j]<-mean(x[j:(j+width)],na.rm=T)}
        }
        
        mps[is.infinite(mps)]<-NA
        mpsf[is.infinite(mpsf)]<-NA
        
        out<-NULL
        for(k in seq_along(mps)){
          out[k]<-mean(c(mps[k],mpsf[k]),na.rm=T)
        }
      }
      return(out)
    })
    
    # columns
    vt<-apply(mt,2,function(x){
      mps<-NULL
      mpsf<-NULL
      for(j in 1:length(x)){
        #left
        if(j+width<length(x)){
          if(is.na(x[j])){
            mps[j]<-max(x[j:j+1],na.rm = T)
          } else {mps[j]<-mean(x[j:(j+width)],na.rm=T)}
        }else{
          if(is.na(x[j])){
            mps[j]<-max(x[j:j-1],na.rm = T)
          } else {mps[j]<-mean(x[j:(j-width)],na.rm=T)}
        }
        #right
        if(j-width>0){
          if(is.na(x[j])){
            mpsf[j]<-max(x[j:j-1],na.rm = T)
          } else {mpsf[j]<-mean(x[j:(j-width)],na.rm=T)}
        }else{
          if(is.na(x[j])){
            mpsf[j]<-max(x[j:j+1],na.rm = T)
          } else {mpsf[j]<-mean(x[j:(j+width)],na.rm=T)}
        }
        
        mps[is.infinite(mps)]<-NA
        mpsf[is.infinite(mpsf)]<-NA
        
        out<-NULL
        for(k in seq_along(mps)){
          out[k]<-mean(c(mps[k],mpsf[k]),na.rm=T)
        }
      }
      return(out)
    })
    
    ## mean matrices
    X<-list(t(ht),vt)
    Y<-do.call(cbind,X)
    Y<-array(Y,dim = c(dim(X[[1]]),length(X)))
    mt<-apply(Y,c(1,2),mean,na.rm=T)
    mt[is.nan(mt)]<-NA
  }
  close(pb)
  
  saveRDS(mt,paste0("./output/PA_PO_wgRONA_extrapl_mat_",bn,"_",substr(Sys.time(),1,10),".rds"))
  ## raster value set
  epr<-pr2
  values(epr)<-mt
  epr<-mask(epr,papo)
  writeRaster(epr,paste0("./output/PA_PO_wgRONA_extrapl_",bn,"_",substr(Sys.time(),1,10),".tif"))
}

#### plotting ####
fls<-list.files("./output/RONA/",pattern = ".tif",full.names=T)
wm<-vect("./data/ne_10m_land/ne_10m_land.shp")
for(files in fls){
  epr<-rast(files)
  cl1<-colorRamps::matlab.like2(length(unique(values(epr),na.rm=T)))
  rng<-ext(-5,100,35,75)
  nn<-gsub("PA_PO_wgRONA_extrapl_","",basename(files))
  plot(wm,border="grey85",col="grey99",ext=rng,mar=c(3.1, 3.1, 2.1, 8.1),main=paste0("Risk of Nonadaptiveness",nn))
  plot(epr,add=T,smooth=T,col=cl1)
  
}

## 7. Pi - NUCLEOTIDE DIVERSITY ####
library(terra)
wm<-vect("./data/level3/level3.shp")
papo<-vect("./data/PAandPO_combined_map.shp")
source("./scripts/raster_functions.R")

### Extrapolation ####
# see the "Extrapolation" section for an explanation of the spatial extrapolation
# used here
pi<-read.table("./data/pi_pop.txt",h=T,col.names = c("pop","lon","lat","pi"))
crd<-read.table("./data/coords_w_hyb_and_gbif_oct28.txt",h=T)
tm<-crd[crd$Source%in%pi$pop,]
pr<-rast_index(pi,resolution = .25,extent=20,field="pi")
wo<-dup_mean(pr,pi)
pr2<-rast_index(wo,resolution = .25,extent=20,field="pi")

#### Interpolate (fill gaps) with nearest neighbor sliding window ######
mt<-matrix(values(pr2),nrow = dim(pr2)[1],ncol = dim(pr2)[2],byrow = T)
wind<-50
width<-3
for(i in 1:wind){
  ## progress bar
  pb <- txtProgressBar(min = 0, max = wind, style = 3, width = 50, char = "=")
  setTxtProgressBar(pb, i)
  ##---
  ## rows
  ht<-apply(mt,1,function(x){
    mps<-NULL
    mpsf<-NULL
    for(j in 1:length(x)){
      #left
      if(j+width<=length(x)){
        if(is.na(x[j])){
          mps[j]<-max(x[j:j+1],na.rm = T)
        } else {mps[j]<-mean(x[j:(j+width)],na.rm=T)}
      }else{
        if(is.na(x[j])){
          mps[j]<-max(x[j:j-1],na.rm = T)
        } else {mps[j]<-mean(x[j:(j-width)],na.rm=T)}
      }
      #right
      if(j-width>0){
        if(is.na(x[j])){
          mpsf[j]<-max(x[j:j-1],na.rm = T)
        } else {mpsf[j]<-mean(x[j:(j-width)],na.rm=T)}
      }else{
        if(is.na(x[j])){
          mpsf[j]<-max(x[j:j+1],na.rm = T)
        } else {mpsf[j]<-mean(x[j:(j+width)],na.rm=T)}
      }
      
      mps[is.infinite(mps)]<-NA
      mpsf[is.infinite(mpsf)]<-NA
      
      out<-NULL
      for(k in seq_along(mps)){
        out[k]<-mean(c(mps[k],mpsf[k]),na.rm=T)
      }
    }
    return(out)
  })
  
  ## columns
  vt<-apply(mt,2,function(x){
    mps<-NULL
    mpsf<-NULL
    for(j in 1:length(x)){
      #left
      if(j+width<length(x)){
        if(is.na(x[j])){
          mps[j]<-max(x[j:j+1],na.rm = T)
        } else {mps[j]<-mean(x[j:(j+width)],na.rm=T)}
      }else{
        if(is.na(x[j])){
          mps[j]<-max(x[j:j-1],na.rm = T)
        } else {mps[j]<-mean(x[j:(j-width)],na.rm=T)}
      }
      #right
      if(j-width>0){
        if(is.na(x[j])){
          mpsf[j]<-max(x[j:j-1],na.rm = T)
        } else {mpsf[j]<-mean(x[j:(j-width)],na.rm=T)}
      }else{
        if(is.na(x[j])){
          mpsf[j]<-max(x[j:j+1],na.rm = T)
        } else {mpsf[j]<-mean(x[j:(j+width)],na.rm=T)}
      }
      
      mps[is.infinite(mps)]<-NA
      mpsf[is.infinite(mpsf)]<-NA
      
      out<-NULL
      for(k in seq_along(mps)){
        out[k]<-mean(c(mps[k],mpsf[k]),na.rm=T)
      }
    }
    return(out)
  })
  
  ## mean matrices
  X<-list(t(ht),vt)
  Y<-do.call(cbind,X)
  Y<-array(Y,dim = c(dim(X[[1]]),length(X)))
  mt<-apply(Y,c(1,2),mean,na.rm=T)
  mt[is.nan(mt)]<-NA 
}
close(pb)

saveRDS(mt,paste0("./output/PA_PO_pi_extrapl_mat_",substr(Sys.time(),1,10),".rds"))
## raster value set
epr<-pr2
values(epr)<-mt
epr<-mask(epr,papo)
writeRaster(epr,paste0("./output/PA_PO_pi_extrapl_",substr(Sys.time(),1,10),".tif"))

### Plotting locally ####
epr<-rast("./output/PA_PO_pi_extrapl.tif")
cl1<-colorRamps::matlab.like2(length(unique(values(epr),na.rm=T)))
rng<-ext(-5,100,35,75)
plot(wm,border="grey85",col="grey99",ext=rng,mar=c(3.1, 3.1, 2.1, 8.1),main=expression(paste("Distribution of ", pi)))
plot(epr,add=T,smooth=T,col=cl1)
plot(papo,border="grey80",add=T)


## 8. PCAdapt ####
library(pcadapt)
### 8.1 Data -------------------------------------------------------------------
bd<-read.pcadapt("/.data/RussClin_filtered8.BED.bed",type = "bed")
### 8.2 PCADAPT and plot ####
pc<-pcadapt(bd,K=20)
plot(pc,option="screeplot",K=20)

info<-read.table("./data/sample_vs_pop_names.txt",h=T)
pop<-info$population
scors<-pc$scores

#### Plot to get optimum K -----------------------------------------------------
pdf("./output/PCs_pop3.pdf",h=8,w=12)
cl<-hcl.colors(length(unique(pop)),palette = "dark2")
par(mfrow=c(2,2))
par(mar=c(4,4,2,3))
for(i in 1:19){
  plot(scors[,i+1]~scors[,i],pch=19,cex=0.7,col=cl[factor(pop)],bty="n",xlab=paste0("PC ",i),ylab=paste0("PC ",i+1),typ="n")
  text(scors[,i+1]~scors[,i],labels=info$ID,cex=0.7,col=cl[factor(pop)])
}
dev.off()

#### Plotting for optimum K=5 --------------------------------------------------
sn<-read.table("./data/snp_names_v8.txt",h=T)
sn.sorted<-read.table("./data/sorted.txt",h=F)
id.sorted<-paste0(sn.sorted$V1,".",sn.sorted$V2)
#### Manhatan plot for peaks ---------------------------------------------------
# K=5 is optimal >>> make manhattan plot with peaks
K=5
pc2<-pcadapt(bd,K=K)
pvals<-na.omit(cbind(pc2$pvalues,pc2$chi2.stat,sn))
#sorted snps
p.snp<-paste0(pvals[,3],".",pvals[,4])
pvals<-pvals[match(id.sorted,p.snp),]
pvals$lg<-sn.sorted[,6]
pvals<-na.omit(pvals)
colnames(pvals)<-c("pval","chi2","chr","pos","lg")
ps<-(-as.numeric(pchisq(pvals$chi2,df = K, lower.tail = FALSE, log.p = TRUE)/log(10)))

## sliding window
wnd=30
mps<-NULL
thr<-quantile(ps[ps<quantile(ps,p=.95)],p=.95)
for(j in 1:length(ps)){
  pb <- txtProgressBar(min = 0, max = length(ps), style = 3, width = 50, char = "=")
  setTxtProgressBar(pb, j)
  if(j+wnd<length(ps)){mps[j]<-mean(ps[j:(j+wnd)])}else{mps[j]<-mean(ps[j:(j-wnd)])}
  #Sys.sleep(time = 1)
}
mps[mps<thr]<-0

cl2<-hcl.colors(n=length(unique(pvals$lg)),palette = "Emrld")
plot(ps,pch=19,cex=0.5,col=cl2[pvals$lg],xlab="Alleles",ylab="-log10(p-values)",bty="n",main=paste0("PCAdapt plot with K=",K," (window: ",wnd,")"))
points(mps,type = "l")

## 9. Hybrid Index (HI) ####
### Extrapolation of hybrid index ---------------------------------------------
# plotting hybrid index on extrapolated map
library(terra)
source("./scripts/raster_functions.R")

wm<-vect("./data/level3/level3.shp")
papo<-vect("./data/PAandPO_combined_map.shp")
fl<-"./data/hybrid_pop.txt"
crd<-read.table("./data/coords_w_hyb_and_gbif_oct28.txt",h=T)

hi<-read.table(fl,h=T)
hi.n<-hi
hi<-hi[,c(3,4,2)]
rownames(hi)<-hi.n$Population
colnames(hi)[1:2]<-c("lon","lat")
# raster indexing
pr<-rast_index(hi,resolution = .25,extent=20,field = "hypop")
wo<-dup_mean(pr,hi)
pr2<-rast_index(wo,resolution = .25,extent=20,field = "hypop")

### 9.1 Interpolate (fill gaps) with nearest neighbor sliding window -----------
mt<-matrix(values(pr2),nrow = dim(pr2)[1],ncol = dim(pr2)[2],byrow = T)
wind<-50
width<-3
for(i in 1:wind){
  ## progress bar
  pb <- txtProgressBar(min = 0, max = wind, style = 3, width = 50, char = "=")
  setTxtProgressBar(pb, i)
  ## rows
  ht<-apply(mt,1,function(x){
    mps<-NULL
    mpsf<-NULL
    for(j in 1:length(x)){
      #left
      if(j+width<=length(x)){
        if(is.na(x[j])){
          mps[j]<-max(x[j:j+1],na.rm = T)
        } else {mps[j]<-mean(x[j:(j+width)],na.rm=T)}
      }else{
        if(is.na(x[j])){
          mps[j]<-max(x[j:j-1],na.rm = T)
        } else {mps[j]<-mean(x[j:(j-width)],na.rm=T)}
      }
      #right
      if(j-width>0){
        if(is.na(x[j])){
          mpsf[j]<-max(x[j:j-1],na.rm = T)
        } else {mpsf[j]<-mean(x[j:(j-width)],na.rm=T)}
      }else{
        if(is.na(x[j])){
          mpsf[j]<-max(x[j:j+1],na.rm = T)
        } else {mpsf[j]<-mean(x[j:(j+width)],na.rm=T)}
      }
      
      mps[is.infinite(mps)]<-NA
      mpsf[is.infinite(mpsf)]<-NA
      
      out<-NULL
      for(k in seq_along(mps)){
        out[k]<-mean(c(mps[k],mpsf[k]),na.rm=T)
      }
    }
    return(out)
  })
  
  ## columns
  vt<-apply(mt,2,function(x){
    mps<-NULL
    mpsf<-NULL
    for(j in 1:length(x)){
      #left
      if(j+width<length(x)){
        if(is.na(x[j])){
          mps[j]<-max(x[j:j+1],na.rm = T)
        } else {mps[j]<-mean(x[j:(j+width)],na.rm=T)}
      }else{
        if(is.na(x[j])){
          mps[j]<-max(x[j:j-1],na.rm = T)
        } else {mps[j]<-mean(x[j:(j-width)],na.rm=T)}
      }
      #right
      if(j-width>0){
        if(is.na(x[j])){
          mpsf[j]<-max(x[j:j-1],na.rm = T)
        } else {mpsf[j]<-mean(x[j:(j-width)],na.rm=T)}
      }else{
        if(is.na(x[j])){
          mpsf[j]<-max(x[j:j+1],na.rm = T)
        } else {mpsf[j]<-mean(x[j:(j+width)],na.rm=T)}
      }
      
      mps[is.infinite(mps)]<-NA
      mpsf[is.infinite(mpsf)]<-NA
      
      out<-NULL
      for(k in seq_along(mps)){
        out[k]<-mean(c(mps[k],mpsf[k]),na.rm=T)
      }
    }
    return(out)
  })
  
  ## mean matrices
  X<-list(t(ht),vt)
  Y<-do.call(cbind,X)
  Y<-array(Y,dim = c(dim(X[[1]]),length(X)))
  mt<-apply(Y,c(1,2),mean,na.rm=T)
  mt[is.nan(mt)]<-NA
}
close(pb)

saveRDS(mt,paste0("./output/PA_PO_hybInd_expl_",substr(Sys.time(),1,10),".rds"))
## raster value set
epr<-pr2
values(epr)<-mt
epr<-mask(epr,papo)
writeRaster(epr,paste0("./output/PA_PO_hybInd_expl_",substr(Sys.time(),1,10),".tif"))

#### plotting locally ----------------------------------------------------------
hi_map<-rast("./output/PA_PO_hybInd_expl.tif")
wm<-vect("./data/level3/level3.shp")
cols<-colorRampPalette(c(4,2))
cl1<-cols(length(unique(values(hi_map),na.rm=T)))
cl1<-colorRamps::matlab.like2(length(unique(values(hi_map),na.rm=T)))
rng<-ext(-5,100,35,75)
cl<-hcl.colors(length(unique(values(hi_map),na.rm=T)),"Blue-Red 2")
plot(wm,border="grey85",col="grey99",ext=rng,mar=c(3.1, 3.1, 2.1, 8.1),main="Hybrid Index")
plot(hi_map,add=T,smooth=T,col=cl1,axes=F,box=F,legend=T)


## 10. Niche centrality ####

### 10.1 Pi vs niche optimum ----
pi<-read.table("./data/pi_pop.txt",h=T)
pi_sorted<-pi[match(env$Source,pi[,1]),]


### 10.2. Linear model of Pi vs abundance --------------------------
library(terra)
epr<-rast("./output/PA_PO_pi_extrapl.tif")
pa_dis<-rast("./douput/PA_dismoMaxEnt.tif")
po_dis<-rast("./output/PO_dismoMaxEnt.tif")

env<-read.table("./data/env_wc_noGbif2.txt",h=T)
env<-env[env$stat!="pt_hy",]
coords<-env[,1:4]
pa_coords<-coords[coords$Species=="PA",]
pa_abun<-extract(pa_dis,pa_coords[,1:2])

po_coords<-coords[coords$Species=="PO",]
po_abun<-extract(po_dis,po_coords[,1:2])

pi_pa<-pi[match(pa_coords$Source,pi$pop),]
pi_po<-pi[match(po_coords$Source,pi$pop),]

pa_pi_ab<-data.frame(pi_pa,abund=pa_abun$`PA_dismoMaxEnt_2022-11-08`,species="PA")
po_pi_ab<-data.frame(pi_po,abund=po_abun$`PO_dismoMaxEnt_2022-11-09`,species="PO")

## linear model
summary(pa_lm<-lm(pi_pa$pi~pa_abun$`PA_dismoMaxEnt_2022-11-08`))
summary(po_lm<-lm(pi_po$pi~po_abun$`PO_dismoMaxEnt_2022-11-09`))
# primitive plot
plot(pi_pa$pi~pa_abun$`PA_dismoMaxEnt_2022-11-08`,pch=19, xlab="Species abundance",ylab="Nucleotide diversity")
abline(pa_lm,lty=2)
points(pi_po$pi~po_abun$`PO_dismoMaxEnt_2022-11-09`,col=4,pch=19)
abline(po_lm,col=4,lty=2)

# better plot with p values and lm
library(ggplot2)
library(viridis) 
library(ggpmisc)
data<-rbind(pa_pi_ab,po_pi_ab)
# Create a scatter plot with ggplot and add colors for each species
(gg <- ggplot(data, aes(x = abund, y = pi, color = species)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, formula = y ~ x, fullrange = TRUE) +
    scale_color_manual(values = c("#1FA7EA", "#EE8900")) +
    labs(x = "Species Abundance", y = "Genetic Diversity (π)", title = "Niche optimum Vs. π by Species") +
    theme_minimal() +
    stat_poly_eq(
      aes(label = paste0("p-value =",format(stat(p.value), digits = 2))),
      formula = y ~ x,
      parse = TRUE,
      size = 3,
      label.x = "right",
      label.y = "top"
    ))


## Additional plotting with ggplot ---------------------------------------------
### Plotting maps with sf and same map region for all maps----------------------
#### Pi ------------------------------------------------------------------------
library(ggplot2)
library(terra)
library(stars)

setwd("/Users/piyalkaru/Desktop/UPP/RUS_CLIN/github/RC_EAA")
epr<-rast("./output/PA_PO_pi_extrapl_2023-06-11.tif")
wb<-vect("./data/level3/level3.shp")
rng<-ext(-10,98,35,75)

# Create a color palette from blue to red
color_palette <- colorRampPalette(c(4, 2))
# extrapolated pi raster
epr.c<-crop(epr,wb)
x_eq<-st_as_stars(epr.c)
x_eq<-st_as_sf(x_eq)
names(x_eq)[1]<-"pi"
# world borders
wb_c<-crop(wb,rng)
wb_sf <- st_as_sf(wb_c)
wb_cod3<-wb_sf["LEVEL3_COD"]

bgmap<-ggplot()+
  geom_sf(data=wb_cod3,color="grey90",fill="grey80",alpha=0.9)

pi_map<-bgmap+
  geom_sf(data=x_eq,aes(fill=pi),color=NA)+
  scale_fill_gradientn(colors=color_palette(100),name="Pi") +
  coord_sf(expand = F) + 
  theme_minimal()

print(pi_map)


#### Addding points to the map ---------------
library(ggplot2)
library(terra)
library(stars)

setwd("/Users/piyalkaru/Desktop/UPP/RUS_CLIN/github/RC_EAA")
epr<-rast("./output/PA_PO_pi_extrapl.tif")
env<-read.table("./data/env_wc_noGbif2.txt",h=T)
pi<-read.table("./data/pi_pop.txt",h=T)
pi_sorted<-pi[match(env$Source,pi[,1]),]
wb<-vect("./data/level3/level3.shp")
# Create a color palette from blue to red
color_palette <- colorRampPalette(c(4, 2))
rng<-ext(-10,98,35,75)

epr.c<-crop(epr,wb)
x_eq<-st_as_stars(epr.c)
x_eq<-st_as_sf(x_eq)
names(x_eq)[1]<-"pi"
# world borders
wb_c<-crop(wb,rng)
wb_sf <- st_as_sf(wb_c)
wb_cod3<-wb_sf["LEVEL3_COD"]

bgmap<-ggplot()+
  geom_sf(data=wb_cod3,color="grey90",fill="grey80",alpha=0.9)

pi.point<-vect(pi_sorted[,1:4],geom=c("lon","lat"),crs=crs(wb))
pi.point_sf<-st_as_sf(pi.point)
## test
pi_map<-bgmap+
  geom_sf(data=x_eq,aes(fill=pi),color=NA)+
  scale_fill_gradientn(colors=color_palette(100),name="PI") +
  geom_sf(data = pi.point_sf, aes(color=pi),stroke=1,fill="black") +
  scale_color_gradientn(colors = color_palette(100), name = "Pi points") +
  geom_sf(data = pi.point_sf, shape=1,size=4,color="black") +
  coord_sf(expand = F) + 
  theme_minimal()


#### Hybrid Index --------------------------------------------------------------
library(ggplot2)
library(terra)
library(stars)

setwd("/Users/piyalkaru/Desktop/UPP/RUS_CLIN/github/RC_EAA")
hi_map<-rast("./output/PA_PO_hybInd_expl.tif")
hi<-read.table("./data/hybrid_pop.txt",h=T)
wb<-vect("./data/level3/level3.shp")
env<-read.table("./data/env_wc_noGbif2.txt",h=T)
# Create a color palette from blue to red
color_palette <- colorRampPalette(c("#EE8900","orange","#1FA7EA","dodgerblue"))
rng<-ext(-10,98,35,75)

hi.c<-crop(hi_map,wb)
x_eq<-st_as_stars(hi.c)
x_eq<-st_as_sf(x_eq)
names(x_eq)[1]<-"hi"

wb_sf<-crop(wb,rng)
wb_sf <- st_as_sf(wb_sf)
wb_cod3<-wb_sf["LEVEL3_COD"]

bgmap<-ggplot()+
  geom_sf(data=wb_cod3,color="grey90",fill="grey80",alpha=0.9)

colno<-length(unique(values(hi.c),na.rm=T))
## adding hi points
hi.point<-vect(hi[,1:4],geom=c("Long","Lat"),crs=crs(wb))
hi.point_sf<-st_as_sf(hi.point)

hi_map<-bgmap+
  geom_sf(data=x_eq,aes(fill=hi),color=NA)+
  scale_fill_gradientn(colors=color_palette(colno),name="HI") +
  geom_sf(data = hi.point_sf, aes(color=hypop),stroke=1,fill="black") +
  scale_color_gradientn(colors = color_palette(colno), name = "Hi points") +
  geom_sf(data = hi.point_sf, shape=1,size=4,color="black") +
  coord_sf(expand = F) + 
  theme_minimal()


#### SDM plots -----------------------------------------------------------------
library(ggplot2)
library(terra)
library(stars)

tif_path="sdm_tifs"
fl<-list.files(tif_path,pattern = ".tif",full.names = T)
# worldmap
wb<-vect("level3.shp")
rng<-ext(-10,98,35,75)

for(i in seq_along(fl)){
  po.cr<-rast(fl[i])
  po.cr.c<-crop(po.cr,wb)
  names(po.cr.c)<-"presence"
  po.cr.c<-crop(po.cr.c,rng)
  po.points<-as.points(po.cr.c)
  po.p.df<-cbind(geom(po.points),as.data.frame(po.points))
  po.cr.rasterpoints<-data.frame(po.p.df[,c(3,4,6)])
  
  wb_sf<-crop(wb,rng)
  wb_sf <- st_as_sf(wb_sf)
  wb_cod3<-wb_sf["LEVEL3_COD"]
  
  bgmap<-ggplot()+
    geom_sf(data=wb_cod3,color="grey90",fill="grey80",alpha=0.9)
  
  hi_map<-bgmap+
    geom_tile(data=po.cr.rasterpoints,aes(x = x, y = y, fill = presence)) +
    scale_fill_gradientn(colors=rev(terrain.colors(1000)),name="Niche suitability") +
    geom_sf(data=wb_cod3,color="grey90",fill=NA,alpha=0.9)+
    xlab("Longitude") + ylab("Latitude") + ggtitle(paste0(gsub(".tif","",basename(fl[i]))))+
    coord_sf(expand = F) + 
    theme_minimal()
  
  pdf(paste0(gsub(".tif", "", basename(fl[i])), ".pdf"),
      width = 8.3 ,  # Width in inches (8.3 inches is A5 width)
      height = 5.8 , # Height in inches (5.8 inches is A5 height)
  )
  print(hi_map)
  dev.off()
}




### Elevation and topography details for the maps ------------------------------
# downloaded from https://www.naturalearthdata.com/downloads/10m-raster-data/10m-natural-earth-2/
library(ggplot2)
library(terra)
library(stars)
elv<-rast("NE2_LR_LC.tif")
elv<-elv["NE2_LR_LC_1"]
rng<-ext(-10,98,35,75)

elv_c<-crop(elv,rng)
elv_sf<-st_as_stars(elv_c)
elv_sf<-st_as_sf(elv_sf)
color_palette2 <- colorRampPalette(c("grey70","black"))

topomap<-ggplot() +
  geom_sf(data=elv_sf,aes(fill=ELEV),color=NA) +
  scale_fill_gradientn(colors = color_palette2(100),name="Elevation (m)")+
  theme_minimal()

#### Climatic clusters ---------------------------------------------------------
# use topmap from the elevation and topo maps
setwd("/Users/piyalkaru/Desktop/UPP/RUS_CLIN/github/RC_EAA")
km<-readRDS("./output/kmeans_clustering.rds")
km2<-km$nonhyb
km3<-km$hyb
x<-km2$cluster
x2<-km3$cluster
clustab<-data.frame(pop=c(names(x),names(x2)),cluster=c(km2$cluster,paste0("h",km3$cluster)))

crd<-data.frame(data.table::fread("./data/env_wc_noGbif2.txt", h = TRUE))
rr2<-data.frame(crd[match(clustab$pop,crd$Source),c(1:2,5)],clustab)

# convert rr2 to sf object
rr2_v<-vect(rr2,geom=c("Long","Lat"),crs=crs(elv))
rr2_sf<-st_as_sf(rr2_v)

# cluster labels for the map
cluster_labels <- c("C. P.abies", "S. P.obovata","N. P.abies","S. P.abies", "N. P.obovata","hybrid 1","hybrid 2","hybrid 3")
stat_labels <- c("Hybrids", "Non-hybrids")

topomap<-ggplot() +
  geom_sf(data=elv_sf,aes(fill=ELEV),color=NA) +
  scale_fill_gradientn(colors = color_palette2(100),name="Elevation (m)")+
  theme_minimal()
# Add points to the existing 'topomap' plot and customize legend labels
topomap +
  geom_sf(data = rr2_sf, aes(color = as.factor(cluster), shape = stat,size=1)) +
  scale_color_discrete(name = "Climatic Clusters", labels = cluster_labels) +  
  scale_shape_discrete(name = "Population type", labels = stat_labels) + 
  #geom_sf_text(data = rr2_sf, aes(label = pop), check_overlap = T,position = "jitter")+
  coord_sf(expand = FALSE)+ xlab("Longitude")+ylab("Latitude")+
  theme_minimal()

#### RONA maps -------------------------------
library(ggplot2)
library(terra)
library(stars)

setwd("/Users/piyalkaru/Desktop/UPP/RUS_CLIN/github/RC_EAA")
fls<-list.files("./output/",pattern = "RONA",full.names=T)
hi_map<-rast(fls[3])
wb<-vect("./data/level3/level3.shp")
rng<-ext(-10,98,35,75)

hi.c<-crop(hi_map,wb)
x_eq<-st_as_stars(hi.c)
x_eq<-st_as_sf(x_eq)
names(x_eq)[1]<-"hi"

wb_sf<-crop(wb,rng)
wb_sf <- st_as_sf(wb_sf)
wb_cod3<-wb_sf["LEVEL3_COD"]

bgmap<-ggplot()+
  geom_sf(data=wb_cod3,color="grey90",fill="grey80",alpha=0.9)

colno<-length(unique(values(hi.c),na.rm=T))
cl1<-colorRamps::matlab.like2(colno)

hi_map<-bgmap+
  geom_sf(data=x_eq,aes(fill=hi),color=NA)+
  scale_fill_gradientn(colors=cl1,name="RONA") +
  coord_sf(expand = F) + 
  xlab("Longitude")+ylab("Latitude")+ggtitle("Bio9")+
  theme_minimal()

print(hi_map)


## EXTRAPOLATION ---------------------------------------------------------------
# see the map_extrapolation Rmd or .html
