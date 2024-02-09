
library( Seurat )


obj = readRDS('/local/projects/idea/sament2/etoh/networks/etoh_final_obj_wsex_wgroup.rds')
 
spn = subset( obj , idents = c('dSPN','iSPN','eSPN') )
# excluding samples 4 and 8 due to substantial differences in sequencing depth
spn = subset( spn , sample_number %in% c(3,7) )

### ALRA (Kluger method)

source('/home/sament/ALRA/alra.R')

norm = t(as.matrix(GetAssayData(spn)))
pct = colSums(norm>0)/nrow(norm)
norm = norm[ , pct > 0.03 ]

res = alra( norm )
alra.out = res[[3]]
rownames(alra.out) = colnames(spn)

m = spn@meta.data
alra = CreateSeuratObject( counts = t(alra.out) , meta.data = m )
alra = NormalizeData( alra )
Idents(alra) = Idents(spn)

alra = FindVariableFeatures(alra)
alra = ScaleData( alra )
alra = RunPCA(alra)
alra = RunUMAP( alra , dims = 1:10 )
alra = FindNeighbors( alra , dims = 1:10 )
alra = FindClusters( alra )

markers = FindAllMarkers( alra , only.pos = T , max.cells.per.ident = 1000 )

pdf('SPN.ALRA.UMAP.pdf')
DimPlot( alra , label = T ) + NoLegend()
FeaturePlot( alra , label = T , combine = F , 
	     features = c('Drd1','Drd2','Drd3','Oprm1') )
dev.off()


saveRDS( alra , file = 'EtOH.SPN.ALRA.seurat.rds' )

write.table( markers , quote=F , sep='\t' , row.names=F ,
	     file = 'EtOH.SPN.ALRA.markers.txt' )




