# compare human AUD DEGs in striatum to our data
library( limma )
library( RRHO )

cn = read.csv('/autofs/burnsfs/projects-t3/idea/ewild/etoh/ZillichDEGs/Zillich_caudateDEGs_munged.csv')

put = read.csv('/autofs/burnsfs/projects-t3/idea/ewild/etoh/ZillichDEGs/Zillich_putamenDEGs_munged.csv')

# mouse

meta = read.delim('/local/projects/idea/sament2/etoh/meta/EtOH.mmu.dstr.degs.combined_meta.min_100_cells_per_batch.2023-12-11.txt')
attr = read.delim('/local/projects/idea/sament2/etoh/networks/cytoscape.attributes.txt')
colors = readRDS('/local/projects/idea/sament2/etoh/networks/k50.merged.clusters.rds')

cn$score = -log10(cn$Zill_C_Pvalue) * sign(cn$Zill_C_log2FC)
put$score = -log10(put$Zill_P_Pvalue) * sign(put$Zill_P_log2FC)

# DEG overlap
# RRHO

types = unique( meta$celltype )
n = length(types)
rrho.out = list()
human = cn[,c('gene','score')]
human = human[ duplicated(human[,1]) == F , ]
human = human[ human[,1] != '' , ]
rownames(human) = human[,1]
for( i in 1:n ) {
  cat( types[i] , '...' )
  mouse = meta[ meta$celltype == types[i] , c('gene','metaz') ]
  rownames(mouse) = mouse[,1]
  m = merge( mouse , human )
  isec = intersect( mouse[,1] , human[,1] )
  l1 = m[,1:2]
  l2 = m[,c(1,3)]
  system( paste( 'mkdir' , types[i] ) )
  setwd( types[i] )
  out.i = RRHO( list1 = l1 ,
	        list2 = l2 ,
		alternative = 'enrichment' , 
		plots = T , 
		BY = F ,
		outputdir = getwd() ,
		labels = c( types[i] , 'human' ) ,
		log10.ind = T )
  rrho.out[[i]] = out.i
  setwd('..')
  cat( max(out.i$hypermat) , '\n' )
}

pdf('hypermat.plots.pdf')
for( i in 1:n ) {
  mtx = rrho.out[[i]]$hypermat
  heatmap( x = as.matrix(mtx) ,
	   col = magma(100) ,
	   scale = 'none' ,
	   breaks = seq(0,6,length.out=101) ,
	   Colv=NA,Rowv=NA )
}
dev.off()

# a simple metric
# across all cell types, is there directional overlap

# down-regulated DEGs
  degs = meta$gene[ meta$metaq < 0.05 & 
		   meta$metaz < 0 ]
  idx = which( cn$gene %in% degs )
  p.dn = geneSetTest( index = idx ,
               statistics = cn$score ,
               alternative = 'down' )
  p.up = geneSetTest( index = idx ,
               statistics = cn$score ,
               alternative = 'up' )
  t = table( cn$gene %in% degs , cn$Zill_C_log2FC < 0 )
  test = fisher.test( t , alternative = 'greater' )

# p.dn = 0.0027
# p.up = 0.997
# fisher: OR = 1.53 , p = 0.004
cn[ cn$gene %in% degs & cn$Zill_C_log2FC < 0 , c(2,3,4) ]

          gene Zill_C_log2FC Zill_C_Pvalue
337     Calhm2  -0.100622693   0.001896977
591       Apoe  -0.144732336   0.004380939
673      Degs2  -0.158104435   0.005380850
685     Atp10a  -0.121850068   0.005586645
1178    Calcrl  -0.116926564   0.013538623
1457      Sgcz  -0.081575413   0.019622275
1820    Mfsd4a  -0.095764375   0.028040599
2283   Slc2a13  -0.074880732   0.038950285
2370    Plppr5  -0.067217919   0.041126505
2510      Dlk2  -0.085595015   0.045307705
2643    Tmeff2  -0.070729030   0.049502352

# up-regulated DEGs
  degs = meta$gene[ meta$metaq < 0.05 &
                   meta$metaz > 0 ]
  idx = which( cn$gene %in% degs )
  p.dn = geneSetTest( index = idx ,
               statistics = cn$score ,
               alternative = 'down' )
  p.up = geneSetTest( index = idx ,
               statistics = cn$score ,
               alternative = 'up' )
  t = table( cn$gene %in% degs , cn$Zill_C_log2FC > 0 )
  test = fisher.test( t , alternative = 'greater' )

# p.dn = 0.35
# p.up = 0.65
# fisher: OR = 0.9 , p = 0.72


n = max(colors)
p.dn = rep(NA,n)
p.up = rep(NA,n)
p.fisher = rep(NA,n)
OR.fisher = rep(NA,n)
nUp = rep(NA,n)
nDown = rep(NA,n)
for( i in 1:max(colors) ) {
  mod = names(colors)[ colors == i ]
  idx = which( cn$gene %in% mod )
  p.dn[i] = geneSetTest( index = idx ,
	       statistics = cn$score ,
	       alternative = 'down' )
  p.up[i] = geneSetTest( index = idx ,
               statistics = cn$score ,
               alternative = 'up' )
  t = table( cn$gene %in% mod , cn$Zill_C_log2FC > 0 )
  test = fisher.test( t )
  p.fisher[i] = test$p.value
  OR.fisher[i] = test$estimate
  nUp[i] = t[2,2]
  nDown[i] = t[2,1]
}

res.cn = data.frame( 
  region = 'caudate' ,
  mod = paste( 'M' , 1:n , sep = '' ) ,
  geneSetTest.p.up = p.up ,
  geneSetTest.p.down = p.dn ,
  nGenesUp = nUp ,
  nGenesDown = nDown ,
  fisher.p = p.fisher ,
  fisher.OR = OR.fisher )
res.cn = res.cn[ order( res.cn$fisher.p ) , ]

# M43
i=43
  mod = names(colors)[ colors == i ]
  t = table( cn$gene %in% mod , cn$Zill_C_log2FC < 0 )
test = fisher.test( t )
 
m43.cn = cn[ cn$gene %in% mod , ]
m43.cn[ order( m43.cn$Zill_C_log2FC ) , ]

write.table( res.cn , sep='\t' , quote=F , row.names=F ,
	     file = 'zillich_vs_modules.overlap_table.txt' )

# putamen

p.dn = rep(NA,n)
p.up = rep(NA,n)
p.fisher = rep(NA,n)
OR.fisher = rep(NA,n)
nUp = rep(NA,n)
nDown = rep(NA,n)
for( i in 1:max(colors) ) {
  mod = names(colors)[ colors == i ]
  idx = which( put$gene %in% mod )
  p.dn[i] = geneSetTest( index = idx ,
               statistics = put$score ,
               alternative = 'down' )
  p.up[i] = geneSetTest( index = idx ,
               statistics = put$score ,
               alternative = 'up' )
  t = table( put$gene %in% mod , put$Zill_P_log2FC > 0 )
  test = fisher.test( t )
  p.fisher[i] = test$p.value
  OR.fisher[i] = test$estimate
  nUp[i] = t[2,1]
  nDown[i] = t[2,2]
}

res.put = data.frame(
  region = 'putamen' ,
  mod = paste( 'M' , 1:n , sep = '' ) ,
  geneSetTest.p.up = p.up ,
  geneSetTest.p.down = p.dn ,
  nGenesUp = nUp ,
  nGenesDown = nDown ,
  fisher.p = p.fisher ,
  fisher.OR = OR.fisher )
res.put = res.put[ order( res.put$fisher.p ) , ]

# actually, the distribution of p-values / fold changes in these putamen data look really off (the original ones from the paper). Very imbalanced. I dont think I'll use these data. Just the Caudate, which looks OK.

# which M43 genes are down-regulated in human AUD?

mod = names(colors)[ colors == 43 ]

cn.mod = cn[ cn$gene %in% mod , ]
cn.mod = cn.mod[ order( cn.mod$Zill_C_Pvalue ) , ]
cn.mod = merge( cn.mod , attr , by = 'gene' )
cn.mod = cn.mod[ cn.mod$X != 7461 , ]
cn.mod$mouse.logfc.norm = ( cn.mod$logFC.dspn - 
			   min(cn.mod$logFC.dspn) ) / 
                            min(cn.mod$logFC.dspn)

cn.mod$human.logfc.norm = ( cn.mod$Zill_C_log2FC - 
			   min(cn.mod$Zill_C_log2FC )) / 
                            min(cn.mod$Zill_C_log2FC)

cn.mod$mean.norm = rowMeans( cn.mod[,c('mouse.logfc.norm',
				       'human.logfc.norm')] )
cn.mod$mouserank = rank( cn.mod$logFC.dspn )
cn.mod$humanrank = rank( cn.mod$Zill_C_Pvalue )
cn.mod$meanrank = apply( cn.mod[,c('mouserank','humanrank')] , 1 , median )
highlight = c('Rgs9','Pde10a','Foxp1','Itpr1','Bcl11b','Rarb')
highlight = cn.mod[ cn.mod$mean.norm > -0.1 , ]

pdf('human_vs_mouse.biplot.pdf')
par( bty = 'l' , cex = 1.5 )
plot( x = cn.mod$logFC.dspn ,
      y = cn.mod$Zill_C_log2FC ,
      pch = 19 , col = 'darkgrey' ,
      xlab = 'logFC in Mouse dSPN' , 
      ylab = 'logFC in Human Caudate' )
points( x = highlight$logFC.dspn ,
        y = highlight$Zill_C_log2FC ,
	pch = 19 , col = 'blue' )
text( x = highlight$logFC.dspn - 0.02 ,
        y = highlight$Zill_C_log2FC ,
        highlight$gene , adj = 1 , cex = 0.5 )
abline( h = 0 , v = 0 )
dev.off()

sig = cn.mod[ cn.mod$metap < 0.05 & cn.mod$Zill_C_Pvalue < 0.05 &
	      cn.mod$metaz < 0 & cn.mod$Zill_C_log2FC < 0 , ]


# make a table of Zillich data for mouse DEGs
 
degs = unique( meta$gene[ meta$metaq < 0.05 ] )
celltype.res = list()
types = unique( meta$celltype )
types = types[c(1,4,3,12,10,7,9,11,8,6,2,5)]
n = length(types)
metaz = matrix( NA , ncol = n , nrow = length(degs) )
rownames(metaz) = degs
colnames(metaz) = types
for( i in 1:n ) {
  meta.i = meta[ meta$celltype == types[i] &
		 meta$gene %in% degs , c('gene','metaz') ]
  metaz[ meta.i$gene  , i ] = meta.i$metaz
}

human = cn[ cn$gene %in% degs , ]
  
m = merge( human[,2:5] , metaz , by.x = 1 , by.y = 0 )
m = m[ order( m$Zill_C_Pvalue ) , ]

write.table( m , quote=F , sep='\t' , row.names=F ,
	     file = 'zillich.cn.for_mouse_degs.txt' )




