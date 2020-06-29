library(Momocs)
library(ggplot2)
library(fpc)




#library(cluster)
#library(NbClust)
#library(factoextra)
#setwd('F:/12Dec2016CopyRajaniDellonlySomeFolders/Idaho/Research/Jan2017Work/AnalysisJan12/GeographicLocation/AllTPSDig/')
#setwd('F:/12Dec2016CopyRajaniDellonlySomeFolders/Idaho/Research/Jan2017Work/AnalysisJan12/TitanEarthLakes/TPSDig40LakesMay30/')
#setwd('F:/12Dec2016CopyRajaniDellonlySomeFolders/Idaho/Research/Jan2017Work/AnalysisJan12/TitanEarthLakes/TPSDig/')
#setwd('F:/12Dec2016CopyRajaniDellonlySomeFolders/Idaho/Research/Jan2017Work/AnalysisJan12/TitanEarthLakes/TitanLakes4thJuly/TPSDigOnlyTitanLakes')
#setwd('F:/12Dec2016CopyRajaniDellonlySomeFolders/Idaho/Research/Jan2017Work/AnalysisJan12/AnalysisAugust1st2018/TPSDig')
setwd('F:/12Dec2016CopyRajaniDellonlySomeFolders/Idaho/Research/Jan2017Work/AnalysisJan12/AnalysisAugust1st2018/TPSDigAll')
#setwd('F:/12Dec2016CopyRajaniDellonlySomeFolders/Idaho/Research/Jan2017Work/AnalysisJan12/AnalysisAugust1st2018/StandardLakes')
filenames <- list.files(path=getwd())

numfiles <- length(filenames)

all_csv <- lapply(filenames,import_tps)

temp.out <- lapply(all_csv, function(x) Out(x$coo))
TitanEarthLakes.out=combine(temp.out)




print(TitanEarthLakes.out)
str(TitanEarthLakes.out)
#LakeType <- c(rep("FeiaLacus",4), rep("OneidaLacus",4), rep("RukwaLacus",4), rep("SparrowLacus",4), rep("VanernLacus",4))#, rep("Titan2",38))#, rep("Tec",10), rep("Vol",4), rep("TitanPS",6), rep("TitanLD", 41), rep("TitanPS",30)) #102 lakes
#LakeType <- c(rep("Titan",70),rep("FeiaLacus",4), rep("OneidaLacus",4), rep("RukwaLacus",4), rep("SparrowLacus",4), rep("VanernLacus",4))
LakeType <- c(rep("Titan",13),rep("FeiaLacus",1),rep("Titan",22),rep("OneidaLacus",1),rep("Titan",6), rep("RukwaLacus",1),rep("Titan",3),rep("SparrowLacus",1), rep("Titan",7),rep("VanernLacus",1), rep("Titan",10))
#LakeType <- c(rep("Titan",7),rep("FeiaLacus",1),rep("Titan",22),rep("OneidaLacus",1),rep("Titan",6), rep("RukwaLacus",1),rep("Titan",3),rep("SparrowLacus",1), rep("Titan",7),rep("VanernLacus",1), rep("Titan",10))
#LakeType <- c(rep("Titan",1),rep("VanernLacus",1), rep("Titan",10))
#LakeType <- c(rep("TitanPS",5), rep("Earth",20), rep("TitanPS",6), rep("TitanLD",39), rep("TitanPS",30)) #100 lakes
#LakeType <- c(rep("TitanPS",4), rep("Earth",20), rep("TitanPS",6), rep("TitanLD",10))#, rep("TitanPS",30)) 40 lakes
#LakeType <- c(rep("Titan",1))
LakeType <- data.frame(LakeType)
TitanEarthLakes.out$fac <- LakeType

plot.new()
panel(TitanEarthLakes.out, fac="LakeType", names=NULL, cex.names = 0.7, points=TRUE, points.pch =25, points.cex = 25)

calibrate_harmonicpower(TitanEarthLakes.out, 'efourier')
calibrate_reconstructions(TitanEarthLakes.out,'efourier',(33), range=1:9)
#TitanEarthLakes.f <- efourier(TitanEarthLakes.out, nb.h =6 , norm=TRUE, smooth.it=6)
TitanEarthLakes.f <- efourier(TitanEarthLakes.out, nb.h =14 , norm=TRUE, smooth.it=6)
hcontrib(TitanEarthLakes.f, harm.range = 1:8)
boxplot(TitanEarthLakes.f)



TitanEarthLakes.new <- TitanEarthLakes.f
#TitanEarthLakes.new                            #the fourier coefficient file
negvals <- TitanEarthLakes.new$coe[,30] < 0    # finding the -ve C2 coefficients, since 4 harmonics, C2=10th
print(negvals[])                                  # logical to show the lakes with -ve C2 or 10th coeff


neg_vals_list = c()                           # creates a vector
for (i in 1:81) {                             # no of lakes
  if (negvals[i] == TRUE) {                   # if it is a negative value
    neg_vals_list <- c(neg_vals_list,i)       # this part gets the row numbers/lake ID for the lakes that have -ve C2
  }
  
}

m<-length(neg_vals_list)                      # m is the integer of the number of values with negative C2
print(m)                                      # which is 3 in number
n<-seq(2,56,2)                                # start from 2nd column, go to 16th (as 4 harmonics and 4 coefficients) and sequence by 2 to get even coefficients
print(n)
nlen <- length(n)  
print(nlen)

print(neg_vals_list)                          # row number 2, 3 and 7 have -ve C2s
#print(TitanEarthLakes.new$coe)
for (i in 1:m) {                            #for the -ve C2/10th coefficients ; 3 in this case
  for (j in 1:nlen) {                       #for the number of even coefficients; 2,4,6,8,10,12,14,16
    TitanEarthLakes.new$coe[neg_vals_list[i],n[j]] <- (-1)*(TitanEarthLakes.new$coe[neg_vals_list[i],n[j]])
  }                                         #multiply by -1
}
print(TitanEarthLakes.new$coe)

TitanEarthLakesNew.p <- PCA(TitanEarthLakes.new)
plot(TitanEarthLakesNew.p,1,morpho=FALSE)

TitanEarthLakesNew.p$rotation[,1]

PC1<-TitanEarthLakesNew.p$rotation[,1]
round(PC1,4)
PC2<-TitanEarthLakesNew.p$rotation[,2]
round(PC2,4)  
PC3<-TitanEarthLakesNew.p$rotation[,3]
round(PC3,4)



#clusGap(TitanEarthLakesNew.p,FUNcluster=kmeans(TitanEarthLakesNew.p,3), K.max, B = 100, verbose = TRUE)
#gap_stat <- clusGap(TitanEarthLakesNew.p, FUN = pam, K.max = 10, B = 500,verbose=TRUE)

#print(TitanEarthLakesNew.p$rotation)
# Print the result
#print(gap_stat, method = "firstmax")
#pam.res <- pam(TitanEarthLakesNew.p$rotation, 2)
#pam.res$cluster

# Visualize pam clusters
#fviz_cluster(pam.res, stand = FALSE, geom = "point",
 #            frame.type = "norm")


#you want to select the loadings from your pca data, maybe the first two
#TitanEarthLakesNew.p$x
#Then import ClusGap
#use kmeans or pam and do clustering on the loadings.

#fviz_pca_ind(TitanEarthLakesNew.p, title = "PCA - Titan data", 
#             habillage = "LakeType",  palette = "jco",
#             geom = "point", ggtheme = theme_classic(),
#             legend = "bottom")

par(oma = c(0.2,3,0.2,3))
plot(TitanEarthLakesNew.p, zoom=1,"LakeType",palette = col_solarized, xax=1, yax=2, pch=20, cex=1, 
     chull.filled = FALSE,pos.shp = "range_axes",ellipses=TRUE,conf.ellipses = 0.95,border.shp = col_alpha("#000000", 0.05), lwd.shp = 1.3,
     col.shp = col_alpha("#0000FF", 0.95),axisvar = FALSE,rug = FALSE, labelspoints=TRUE, abbreviate.labelspoints = FALSE, labelsgroups = FALSE, cex.labelsgroups = 1, rect.labelsgroups = TRUE,title = substitute(""),axisnames = FALSE)#,title = substitute("")#,
     #flipx.shp = FALSE, flipy.shp = FALSE,pos.shp = "range_axes",rotate.shp=0,
     #border.shp = col_alpha("#000000", 0.05), lwd.shp = 1.3,
     #col.shp = col_alpha("#0000FF", 0.95),ellipses=TRUE,conf.ellipses = 0.95,labelspoints = FALSE, #col.labelspoints = par("fg"),
     #cex.labelspoints = 0.9, abbreviate.labelspoints = TRUE,
     #labelsgroups = TRUE, cex.labelsgroups = 0.8, rect.labelsgroups = TRUE,
     #abbreviate.labelsgroups = FALSE, color.legend = FALSE, axisnames = FALSE,
     #axisvar = FALSE, unit = FALSE, eigen = FALSE, rug = FALSE,
     #title = substitute("Titan Lakes"), box = TRUE, old.par = FALSE)
mtext("PC1 (30.4%)", side =3, line = -1.5, cex = 1.2)
mtext("PC2 (15.4%)", side = 2, line = 0.25, cex = 1.2)
mtext("PC3 (11.3%)", side =3, line = -1, cex = 1.2)
#mtext("PC2 (15.6%)", side = 2, line = 0, cex = 1.2)
plot(TitanEarthLakesNew.p$x[1:2])
test <- as.data.frame(TitanEarthLakesNew.p$x[,1:2])
#plot(test.df)
# create data frame with scores
#scores = as.data.frame(pca1$x)

# plot of observations
ggplot(data = test, aes(x = PC1, y = PC2, label = rownames(test))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(colour = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot of Titan Lakes")

plot(TitanEarthLakesNew.p,'LakeType')

plot(TitanEarthLakesNew.p, zoom=1,"LakeType",palette = col_solarized, xax=1, yax=2, pch=20, cex=1, ellipses=TRUE,conf.ellipses = 0.95,
     chull.filled =FALSE,labelspoints = TRUE,
     flipx.shp = FALSE, flipy.shp = FALSE,axisnames = TRUE, axisvar = TRUE, pos.shp = "full_axes",rotate.shp=0)

plot(TitanEarthLakesNew.p, "LakeType", xax = 1, yax = 2, points = TRUE,
#col = "#000000", 
pch = 20, cex = 0.5, palette = col_solarized,
center.origin = TRUE, zoom = 1, xlim = NULL, ylim = NULL,
bg = par("bg"), grid = TRUE, nb.grids = 3, morphospace = TRUE,
pos.shp = c("range", "full", "circle", "xy", "range_axes", "full_axes")[5],
amp.shp = 1, size.shp = 1, nb.shp = 12, nr.shp = 6, nc.shp = 5,
rotate.shp = 0, flipx.shp = FALSE, flipy.shp = FALSE, pts.shp = 60,
border.shp = col_alpha("#0000FF", 0.5), lwd.shp = 2,
col.shp = col_alpha("#0000FF", 0.95), stars = FALSE, ellipses = TRUE,
conf.ellipses = 0.95, ellipsesax = FALSE, conf.ellipsesax = c(0.5, 0.9),
lty.ellipsesax = 1, lwd.ellipsesax = sqrt(2), chull = FALSE,
chull.lty = 1, chull.filled = FALSE, chull.filled.alpha = 0.92,
density = FALSE, lev.density = 20, contour = FALSE, lev.contour = 3,
n.kde2d = 100, delaunay = FALSE, loadings = FALSE,
labelspoints = TRUE, col.labelspoints = par("fg"),
cex.labelspoints = 0.6, abbreviate.labelspoints = TRUE,
labelsgroups = TRUE, cex.labelsgroups = 0.8, rect.labelsgroups = FALSE,
abbreviate.labelsgroups = FALSE, color.legend = FALSE, axisnames = FALSE,
axisvar = TRUE, unit = FALSE, eigen = TRUE, rug = TRUE,
title = substitute(x), box = TRUE, old.par = TRUE)


par(oma = c(0.2,3,0.2,3))
CLUST(TitanEarthLakesNew.p, "LakeType" ,type="phylogram",dist_method="euclidean", hclust_method="average",retain=1:6,cex=0.8,palette = col_solarized)
#CLUST(TitanEarthLakesNew.p, "LakeType", type = "unrooted", dist_method = "euclidean",
#hclust_method = "complete", retain = 0.99,palette = col_qual)
set.seed(321463)
temp<-KMEANS(TitanEarthLakesNew.p, centers = 10)
temp$tot.withinss
wss <- sapply(1:10,  function(k){KMEANS(TitanEarthLakesNew.p, centers =k)$tot.withinss})
wss
plot(1:10,wss)
par(oma = c(0.2,3,0.2,3))
KmeansCalc<-KMEANS(TitanEarthLakesNew.p,1, nax = 1:2, pch=27, cex=1.75)
#plot(KmeansCalc, zoom=1,"LakeType",palette = col_solarized, xax=1, yax=2, pch=20, cex=1, ellipses=TRUE,conf.ellipses = 0.95,
#     chull.filled =FALSE,labelspoints = TRUE,
#     flipx.shp = FALSE, flipy.shp = FALSE,axisnames = TRUE, axisvar = TRUE, pos.shp = "full_axes",rotate.shp=0)

KmeansCalc$
KmeansCalc$cluster

####
KmeansCalc2<-KMEANS(TitanEarthLakesNew.p, center = 4, nax = 1:2, pch=12,cex=1)

#, panel=col("blue","green","purple","red","yellow"))
clvecd = c(1,2,3,4)
col=c("blue","green","purple","red")

KmeansCalc$cluster
typeof(TitanEarthLakesNew.p)
typeof(KmeansCalc$cluster)



plotcluster(TitanEarthLakesNew.p
              , KmeansCalc$cluster
              , pch = length(KmeansCalc$cluster)
              , clvecd=c(1,2,3,4)
              , col=c("blue","green","purple","red") )

#vcol <- c("blue","green","purple","red","yellow")
#plot(TitanEarthLakesNew.p, KmeansCalc$cluster, col=vcol[KmeansCalc$cluster])

#clus <- kmeans(x, centers=5)

#plotcluster(x, clus$cluster, col=vcol[clus$cluster])

plotcluster(x, clus$cluster, col=vcol[clus$cluster])
mtext("PC1", side =3, line = -1, cex = 1.2)
mtext("PC2", side = 2, line = 0, cex = 1.2)
KmeansCalc$cluster

#cluster1 = c()                           # creates a vector
#cluster2 = c()
#cluster3 = c()
#for (i in 1:102) {                           # no of lakes
#  if (KmeansCalc$cluster == 1)   {                 # if it is a negative value
#    cluster1 <- KmeansCalc$cluster[i]       # this part gets the row numbers/lake ID for the lakes that have -ve C2
#  }
  

#}

#print(cluster1)

#set.seed(20)
#print(cluster1)
MANOVA(TitanEarthLakesNew.p, "LakeType",retain=6)
MANOVA_PW(TitanEarthLakesNew.p,"LakeType")
Lakes.l <- LDA(TitanEarthLakesNew.p, "LakeType",retain=6)
plot(Lakes.l)
#TitanEarthLakes.f$fac$plop <- factor(rep(letters[1:2], each=51))
#TitanEarthLakes.l <- LDA(PCA(TitanEarthLakes.f), 'plop')
#TitanEarthLakes.l
#plot(TitanEarthLakes.l)


plot(TitanEarthLakesNew.p, ellipses=TRUE)#, ellipsesax = FALSE, pch=c(1, 2), cex=2)
plot(TitanEarthLakesNew.p, chull=TRUE, pos.shp = "full_axes", abbreTeciate.labelsgroups = TRUE, points=FALSE, labelspoints = TRUE)
plot(TitanEarthLakesNew.p, pos.shp="circle", stars=TRUE, chull.filled=TRUE, palette=col_spring)
scree(TitanEarthLakesNew.p)
scree_plot(TitanEarthLakesNew.p)
boxplot(TitanEarthLakesNew.p)
PCcontrib(TitanEarthLakesNew.p)

#to extract how much variance is explained by each PC
var.base <- (TitanEarthLakesNew.p$sdev)^2/sum(TitanEarthLakesNew.p$sdev^2)
var.base*100
round(var.base*100,3)
sum(var.base[1:6])



