# combine pseudobulk p-values from batches 1-3

library( dotgen )

# cell counts per batch

b1.counts = read.csv('/autofs/burnsfs/projects-t3/idea/ewild/etoh/Data_2023/combat8/dec2023/Cell_Dist_B1.csv')

b2.counts = read.csv('/autofs/burnsfs/projects-t3/idea/ewild/etoh/Data_2023/combat8/dec2023/Cell_Dist_B2.csv')

b3.counts = read.csv('/autofs/burnsfs/projects-t3/idea/ewild/etoh/Data_2023/combat8/dec2023/Cell_Dist_B3.csv')

b1.counts[,1] = tolower(b1.counts[,1])
b2.counts[,1] = tolower(b2.counts[,1])
b3.counts[,1] = tolower(b3.counts[,1])
b1.counts[,1] = gsub( 'ependyma','ependy',b1.counts[,1] )
b2.counts[,1] = gsub( 'ependyma','ependy',b2.counts[,1] )
b1.counts[,1] = gsub( 'as','astro',b1.counts[,1] )
b2.counts[,1] = gsub( 'as','astro',b2.counts[,1] )
b1.counts[,1] = gsub( 'ol','oligo',b1.counts[,1] )
b2.counts[,1] = gsub( 'ol','oligo',b2.counts[,1] )
b1.counts[,1] = gsub( 'opc','poly',b1.counts[,1] )
b2.counts[,1] = gsub( 'opc','poly',b2.counts[,1] )

counts = merge( b1.counts , b2.counts , by = 1 , all = T )
counts = merge( counts , b3.counts , by = 1 , all = T )
counts[,1] = gsub( 'in\\-','',counts[,1] )

# DEGs in each batch

B1 = read.csv("/autofs/burnsfs/projects-t3/idea/ewild/etoh/B1_pseudo/B1subset_pseudo_degs/degs_alra_masterfile_EtOH_B1subset_pseudobulk.csv")

B2 = read.csv( "/autofs/burnsfs/projects-t3/idea/ewild/etoh/etoh_obj.new/pseudobulk/subset_B2/degs_alra_masterfile_EtOH_objnew_B2only_pseudobulk.csv")

B3 = read.csv( '/autofs/burnsfs/projects-t3/idea/ewild/etoh/Data_2023/combat8/exprFilter/degs_masterfile_EtOH_B3_pseudobulk_combat8_revised_exprFilter.csv' )

types = unique( B3$celltype )
n = length(types)

rownames(counts) = counts[,1]
counts = counts[,-1]
counts = counts[ types , ]
colnames(counts) = c('b1','b2','b3')

meta = list()
for( i in 1:n ) {
  counts.i = as.vector(t(counts[ types[i] , ]))
  weights = counts.i / sum(counts.i)
  weights = weights * sum(counts.i)/counts.i[3]
  weights[ counts.i < 100 ] = 0

  b1.i = B1[ B1$celltype == types[i] , ]
  b2.i = B2[ B2$celltype == types[i] , ]
  b3.i = B3[ B3$celltype == types[i] , ]
  genes = b3.i[,2]
  b1.i = b1.i[ b1.i[,2] %in% genes , ]
  b2.i = b2.i[ b2.i[,2] %in% genes , ]
  b1.i$z1 = zsc( b1.i$PValue , -1*sign( b1.i$logFC ) )
  b2.i$z2 = zsc( b2.i$PValue , -1*sign( b2.i$logFC ) )
  b3.i$z3 = zsc( b3.i$PValue , -1*sign( b3.i$logFC ) )
  m = merge( b1.i[,c('X1','z1')] , b2.i[,c('X1','z2')] , by = 1 , all=T)
  m = merge( m , b3.i[,c('X1','z3')] , by = 1 , all=T )
  rownames(m) = m$X1
  m = m[,-1]
  weighted = matrix( NA , nrow = nrow(m) , ncol = ncol(m) )
  weighted[,1] = m[,1] * weights[1]
  weighted[,2] = m[,2] * weights[2]
  weighted[,3] = m[,3] * weights[3]
  if(nrow(b2.i)==0) weights = c(0,0,1)
  metaz = rowSums( weighted , na.rm = T ) / sqrt(sum(weights))
  metap = 2*pnorm( metaz )
  metap[ metaz > 0 ] = 2*pnorm( metaz[metaz>0] , lower.tail = F )
  metaq = p.adjust( metap , method = 'fdr' )
  df = data.frame( celltype = types[i] ,
		   gene = rownames(m) ,
		   m ,
		   metaz , metap , metaq )
  meta[[i]] = df
}

meta = do.call( rbind , meta )

meta = meta[ order( meta$metap ) , ]

write.table( meta , row.names=F , quote=F , sep='\t' ,
	     file = 'EtOH.mmu.dstr.degs.combined_meta.min_100_cells_per_batch.2023-12-11.txt' )






