library( Seurat )

obj = readRDS('EtOH.SPN.ALRA.seurat.rds')
obj = NormalizeData(obj)
datExpr = t( as.matrix(GetAssayData(obj)) )

colors = readRDS('k50.merged.clusters.rds')
mods = list()
for( i in 1:max(colors) ) {
  mods[[i]] = names(colors)[ colors == i ]
}
names(mods) = paste('M',1:max(colors),sep='')

obj = AddModuleScore( obj , mods )

AE = obj@meta.data[,grep('Cluster',colnames(obj@meta.data))]

MM = cor( AE , datExpr )
MM = t(MM)
colnames(MM) = gsub('Cluster','M',colnames(MM))

m = MM[,'M33']
head( m[ order(m,decreasing=T) ] , 100 )

saveRDS( MM , file = 'module_membership.kAE.rds' )

