library( RRHO )
library( dotgen )

# first perform a meta-analysis of B1 and B2

# how many cells per type

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

B1 = read.csv("/autofs/burnsfs/projects-t3/idea/ewild/etoh/B1_pseudo/B1subset_pseudo_degs/degs_alra_masterfile_EtOH_B1subset_pseudobulk.csv")

B1$score = -log10(B1$PValue) * sign(B1$logFC)

B2 = read.csv( "/autofs/burnsfs/projects-t3/idea/ewild/etoh/etoh_obj.new/pseudobulk/subset_B2/degs_alra_masterfile_EtOH_objnew_B2only_pseudobulk.csv")

B3 = read.csv( '/autofs/burnsfs/projects-t3/idea/ewild/etoh/Data_2023/combat8/exprFilter/degs_masterfile_EtOH_B3_pseudobulk_combat8_revised_exprFilter.csv' )

B3$score = -log10(B3$PValue) * sign(B3$logFC)

types = unique( B1$celltype )
n = length(types)

rownames(counts) = counts[,1]
counts = counts[,-1]
counts = counts[ types , ]
colnames(counts) = c('b1','b2','b3')

counts = counts[,1:2]
counts = counts[ -which( apply(counts,1,max) < 100 ) , ]

types = rownames(counts)
n = length(types)

meta = list()
for( i in 1:n ) {
  counts.i = as.vector(t(counts[ types[i] , ]))
  weights = counts.i / sum(counts.i)
  weights = weights * sum(counts.i)/counts.i[1]
  weights[ counts.i < 100 ] = 0

  b1.i = B1[ B1$celltype == types[i] , ]
  b2.i = B2[ B2$celltype == types[i] , ]
  # b3.i = B3[ B3$celltype == types[i] , ]
  genes = b1.i[,2]
  b1.i = b1.i[ b1.i[,2] %in% genes , ]
  b2.i = b2.i[ b2.i[,2] %in% genes , ]
  b1.i$z1 = zsc( b1.i$PValue , -1*sign( b1.i$logFC ) )
  b2.i$z2 = zsc( b2.i$PValue , -1*sign( b2.i$logFC ) )
  # b3.i$z3 = zsc( b3.i$PValue , -1*sign( b3.i$logFC ) )
  m = merge( b1.i[,c('X1','z1')] , b2.i[,c('X1','z2')] , by = 1 , all=T)
  # m = merge( m , b3.i[,c('X1','z3')] , by = 1 , all=T )
  rownames(m) = m$X1
  m = m[,-1]
  weighted = matrix( NA , nrow = nrow(m) , ncol = ncol(m) )
  weighted[,1] = m[,1] * weights[1]
  weighted[,2] = m[,2] * weights[2]
  # weighted[,3] = m[,3] * weights[3]
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

write.table( meta , row.names=F , quote=F , sep='\t' ,
	     file = 'EtOH.pseudobulk.B12.meta-analysis.txt' )

types = intersect( meta$celltype , B3$celltype )

n = length(types)
rrho.out = list()
for( i in 1:n ) {
  # system( paste( 'mkdir' , types[i] ) )
  setwd( types[i] )
  b1.i = meta[ meta$celltype == types[i] ,
	   c('gene','metaz') ]
  b3.i = B3[ B3$celltype == types[i] ,
	   c('X1','score') ]
  m = merge( b1.i , b3.i , by = 1 )
  l1 = m[,1:2]
  l2 = m[,c(1,3)]
  out.i = RRHO( list1 = l1 ,
	      list2 = l2 ,
	      labels = c('B1','B3') ,
	      alternative = 'enrichment' ,
 	      stepsize = round(nrow(m)/99) ,
	      outputdir = getwd() ,
	      BY = T ,
	      log10.ind = T ,
	      plots = T )
  rrho.out[[i]] = out.i
  cat( types[i] , max(rrho.out[[i]]$hypermat) , '\n' )
  setwd('..')
}

# overall

m = merge( B3 , meta , by.x = c('X1','celltype') , 
	   by.y = c('gene','celltype') )
t = table( m$FDR < 0.05 , m$metaq < 0.05 )
rep = m[ m$FDR < 0.05 & m$metaq < 0.05 , ]

# 497 DEGs across all cell types in the discovery dataset
 
# 45 of these genes are significant at FDR < 0.05 in both the discovery and replication datasets
# odds ratio = 33.9 , P = 6.2e-50

t = table( m$FDR < 0.05 , m$metap < 0.05 )
 
# 41 of these are concordantly down-regulated in both datasets

# at p < 0.05, 141 of the 497 DEGs are reproducible.

for( i in 1:n ) {
  b1.i = meta[ meta$celltype == types[i] , ]
  b3.i = B3[ B3$celltype == types[i] , ]
  u = intersect( b3.i$X1 , b1.i$gene )
  dn3 = b3.i$X1[ b3.i$PValue < 0.05 & b3.i$logFC > 0 ]
  dn1 = b1.i$gene[ meta$metaz > 0 & meta$metap < 0.05 ]
  t = table( u %in% dn3 , u %in% dn1 )


maxp.up = rep(NA,n)
maxp.dn = rep(NA,n)
maxp = rep(NA,n)
minfdr = rep(NA,n)
for( i in 1:n ) {
  h = rrho.out[[i]]$hypermat
  m = dim(h)[1]
  maxp.up[i] = max(h[1:25,1:25])
  maxp.dn[i] = max(h[76:m,76:m])
  maxp[i] = max(h)
  p = as.matrix(10^(-1*h))
  fdr = p.adjust( as.vector(p) , method = 'fdr' )
  minfdr[i] = min(fdr)
}



data.frame( types , minfdr , maxp , maxp.up , maxp.dn )

  types        minfdr       maxp    maxp.up    maxp.dn
1 oligo  1.000000e+00   1.992602  1.9926020   1.322949
2  ispn  2.555009e-98 101.536151 20.1939228 101.536151
3  poly  3.681521e-08  11.299190 11.1329428   3.880878
4    pv  1.382014e-03   6.443917  2.7938127   6.443917
5  dspn 4.382527e-118 121.349546 39.0866681 121.349546
6    mg  7.010482e-01   4.154252  0.3087645   4.154252
7   sst  1.000000e+00   2.876267  1.0922712   2.876267
8 astro  2.824559e-03   6.291457  6.2914574   4.052033
9  espn  2.883499e-03   5.969852  2.8402678   5.969852


