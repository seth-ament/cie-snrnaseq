
b3 = read.csv( '/autofs/burnsfs/projects-t3/idea/ewild/etoh/Data_2023/combat8/exprFilter/degs_masterfile_EtOH_B3_pseudobulk_combat8_revised_exprFilter.csv' )

meta = read.delim('../meta/EtOH.mmu.dstr.degs.combined_meta.min_100_cells_per_batch.2023-12-11.txt')

MM = readRDS('module_membership.kAE.rds')

colors = readRDS('k50.merged.clusters.rds' )

dspn = b3[ b3$celltype == 'dspn' , c(2,3) ]
ispn = b3[ b3$celltype == 'ispn' , c(2,3) ]
espn = b3[ b3$celltype == 'espn' , c(2,3) ]

dspn.meta = meta[ meta$celltype == 'dspn' , c(2,6,7) ]
ispn.meta = meta[ meta$celltype == 'ispn' , c(2,6,7) ]
espn.meta = meta[ meta$celltype == 'espn' , c(2,6,7) ]
dspn.meta$logP = -log10(dspn.meta$metap)
ispn.meta$logP = -log10(ispn.meta$metap)
espn.meta$logP = -log10(espn.meta$metap)


attr = merge( colors , dspn , by.x = 0 , by.y = 1 )
attr = merge( attr , dspn.meta , by.x = 1 , by.y = 1 )
attr = merge( attr , ispn , by.x = 1 , by.y = 1 )
attr = merge( attr , ispn.meta , by.x = 1 , by.y = 1 )
attr = merge( attr , espn , by.x = 1 , by.y = 1 )
attr = merge( attr , espn.meta , by.x = 1 , by.y = 1 )
attr = merge( attr , MM , by.x = 1 , by.y = 0 )

colnames(attr) = c('gene','module',
  'logFC.dspn','metaz.dspn','metap.dspn','logP.dspn',
  'logFC.ispn','metaz.ispn','metap.ispn','logP.ispn',
  'logFC.espn','metaz.espn','metap.espn','logP.espn',
  colnames(MM) )

write.table( attr , quote=F , sep='\t' , row.names=F ,
	     file = 'cytoscape.attributes.txt' )


