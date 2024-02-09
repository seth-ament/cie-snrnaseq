library( Seurat )


library( biomaRt )
mart = useMart('ensembl')
mart = useDataset( mart = mart , 
		   'mmusculus_gene_ensembl' )
anno = getBM( mart = mart , 
	      attributes = c('mgi_symbol' ,
			     'gene_biotype' ) )


obj = readRDS('etoh_final_obj_wsex_wgroup.rds')

dge = read.delim('../meta/EtOH.mmu.dstr.degs.combined_meta.min_100_cells_per_batch.2023-12-11.txt')


b3 = read.csv( '/autofs/burnsfs/projects-t3/idea/ewild/etoh/Data_2023/combat8/exprFilter/degs_masterfile_EtOH_B3_pseudobulk_combat8_revised_exprFilter.csv' )
colnames(b3)[2] = 'gene'
b3 = b3[,-1]

g = unique( b3$gene )
n = length(g)
logcpm.max = sapply( 1:n , function(i)
  max( b3$logCPM[ b3$gene == g[i] ] ) )
logcpm.max = data.frame(
        gene = g , logcpm.max )

b3 = merge( b3 , logcpm.max , by = 'gene' )
b3$logcpm.diff = b3$logCPM - b3$logcpm.max
b3 = b3[ b3$logcpm.diff > -3.32 , ]

types = unique(dge$celltype)
n = length(types)
dge.filt = list()
for( i in 1:n ) {
  g = b3$gene[ b3$celltype == types[i] ]
  dge.filt[[i]] = dge[ dge$celltype == types[i] &
		       dge$gene %in% g , ]
}
dge.filt = do.call( rbind , dge.filt )

nsig = table( dge.filt$gene[ dge.filt$metaq < 0.05 ] )

nsig = table( b3$gene[ abs(b3$logFC) > 0.25 & 
	               b3$logcpm.max > 6 ] )
nsig = sort( nsig )
nsig = data.frame( gene = names(nsig) ,
		   nsig = as.vector(nsig) )
nsig = merge( anno , nsig , by = 1 )
nsig = nsig[ order( nsig$nsig , decreasing = T ) , ]
 
features = head( nsig$mgi_symbol[ nsig$gene_biotype == 'protein_coding' ] , 20 )
 
levels(obj) = c('dSPN','iSPN','eSPN','IN-PV','IN-SST','IN-CHAT',
		'Astro','Oligo','Poly','Endo','Ependy','Mural','MG')
 
pdf('shared_degs.vioplot.pdf',width=15,height=3)
VlnPlot( obj , 
	 features = c('Ptger4','Egr1','Ctnna3') , 
	 split.by = 'trt' , 
	 pt.size = 0.5 ,
	 split.plot = T ,
	 combine = T )
dev.off()








