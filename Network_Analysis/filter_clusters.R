# this script is used to take results from k-means clustering and
# perform filtering and merging of the results to produce a final network

library( Seurat )
library( WGCNA )

# minimum size of a gene co-expression cluster
minModSize = 25
# minimum correlation between a gene and the cluster centroid
kAEtoStay = 0.3
# threshold for merging clusters with strongly correlated patterns
# the "height" corresponds to 1 - Pearson's r
cutHeight = 0.3

km = readRDS('kmeans.k50.rds')

obj = readRDS('etoh_final_obj_wsex_wgroup.rds')

spn = subset( obj , idents = c('dSPN','iSPN','eSPN') )

impute1 = readRDS('EtOH.SPN.ALRA.seurat.rds')
 
impute1 = NormalizeData(impute1)

datExpr0 = as.matrix( GetAssayData(impute1) )

# Z score
datExpr=sweep(datExpr0,1,apply(datExpr0,1,mean),"-")
indx.sd=(apply(datExpr0,1,sd))==0 # these will produce NAs
datExpr=sweep(datExpr,1,apply(datExpr,1,sd),"/")
datExpr[indx.sd,]=0
if(sum(is.na(datExpr))!=0){print("NAs in exprsMTX.Z Zscores!")}
#kmeans

kAE = cor( t(km$centers) , t(datExpr) )

colors = km$cluster

for( i in 1:length(colors) ) {
  r.i = kAE[ colors[i] , i ] 
  if( r.i < kAEtoStay ) colors[i] = 0
}

size = table( colors )
too.small = as.numeric(names(which(size<minModSize))) 
colors[ colors %in% too.small ] = 0


centers = sapply( sort(unique(colors)) , function(i) 
  colMeans(datExpr[ colors == i , ]) )

colnames(centers) = paste( 'AE' , sort(unique(colors)) , sep = '' )

r = cor( centers )

d = as.dist( 1 - r )

hc = hclust( d , method = 'average' )

cl = cutree( hclust( d , method = 'average' ) , h = cutHeight )


mergeColors = rep(NA,length(colors))
for( i in 1:max(cl) ) {
  idx = as.numeric( gsub( 'AE','', names(cl)[ cl == i ] ))
  mergeColors[ colors %in% idx ] = i
}
mergeColors = mergeColors - 1
names(mergeColors) = names(colors)

# extend to all samples


source('/home/sament/ALRA/alra.R')

spn = NormalizeData(spn)
norm = t(as.matrix(GetAssayData(spn)))
norm = norm[ , rownames(impute1) ]
res = alra( norm )
alra.out = res[[3]]
rownames(alra.out) = colnames(spn)

datExpr2 = alra.out

MEs = moduleEigengenes( datExpr2 , mergeColors )

saveRDS( mergeColors , file = 'k50.merged.clusters.rds' )

saveRDS( MEs , file = 'k50.merged.MEs.rds' )

saveRDS( datExpr2 , file = 'ErOH.ALRA.allSPNs.rds' )

meta = spn@meta.data
meta = cbind( meta , MEs$eigengenes , MEs$averageExpr )
spn@meta.data = meta

saveRDS( spn , file = 'EtOH.allSPNs.wMEs.rds' )


