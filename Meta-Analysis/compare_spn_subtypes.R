# compare the DEGs across MSN subtypes
# venn diagram
# biplots

library( limma )
library( biomaRt )

httr::set_config(httr::config(ssl_verifypeer = FALSE))
mart = useEnsembl( biomart = 'ensembl' ,
		   dataset = 'mmusculus_gene_ensembl' ,
		   mirror = 'useast' )
anno = getBM( mart = mart , attributes = c(
		'ensembl_gene_id',
		'mgi_symbol',
		'gene_biotype' ) )
 
dge = read.delim('EtOH.mmu.dstr.degs.combined_meta.min_100_cells_per_batch.2023-12-11.txt')
 
b3 = read.csv( '/autofs/burnsfs/projects-t3/idea/ewild/etoh/Data_2023/combat8/exprFilter/degs_masterfile_EtOH_B3_pseudobulk_combat8_revised_exprFilter.csv' )
colnames(b3)[2] = 'gene'
b3 = b3[,-1]

u = unique( dge$gene )
n = length(u)
dspn = rep( 0 , n )
dup = dge$gene[ dge$metaq < 0.05 & dge$celltype == 'dspn' & dge$metaz > 0 ]
ddn = dge$gene[ dge$metaq < 0.05 & dge$celltype == 'dspn' & dge$metaz < 0 ]
dspn[ u %in% dup ] = 1
dspn[ u %in% ddn ] = -1
ispn = rep( 0 , n )
iup = dge$gene[ dge$metaq < 0.05 & dge$celltype == 'ispn' & dge$metaz > 0 ]
idn = dge$gene[ dge$metaq < 0.05 & dge$celltype == 'ispn' & dge$metaz < 0 ]
ispn[ u %in% iup ] = 1
ispn[ u %in% idn ] = -1
espn = rep( 0 , n )
eup = dge$gene[ dge$metaq < 0.05 & dge$celltype == 'espn' & dge$metaz > 0 ]
edn = dge$gene[ dge$metaq < 0.05 & dge$celltype == 'espn' & dge$metaz < 0 ]
espn[ u %in% eup ] = 1
espn[ u %in% edn ] = -1

df = data.frame( dspn , ispn , espn )
rownames(df) = u
 
res = vennCounts( df )

  dspn ispn espn Counts
1    0    0    0  20244
2    0    0    1     81
3    0    1    0     29
4    0    1    1      0
5    1    0    0    189
6    1    0    1      6
7    1    1    0     46
8    1    1    1      7

pdf('msn_deg.vann.pdf')
vennDiagram( res )
dev.off()

# are the effects correlated
 
# degs = rownames(df)[ rowSums(df!=0) > 0 ]

degs = unique( dge$gene[ dge$metaq < 0.05 ] )

g = unique( b3$gene )
n = length(g)
logcpm.max = sapply( 1:n , function(i)
  max( b3$logCPM[ b3$gene == g[i] ] ) )
logcpm.max = data.frame(
	gene = g , logcpm.max )

b3 = merge( b3 , logcpm.max , by = 'gene' )
b3$logcpm.diff = b3$logCPM - b3$logcpm.max
b3 = b3[ b3$logcpm.diff > -3.32 , ]

fclist = list()
types = unique( b3$celltype )
types = types[ c(3,4,1,12,7,8,5,2,11,6,9,10) ]
for( i in 1:length(types) ) {
  fc.i = b3[ b3$celltype == types[i] , c('gene','logFC') ]
  rownames(fc.i) = fc.i$gene
  fc.i = fc.i[ degs , ]
  fclist[[i]] = fc.i$logFC
}

logfc.degs.all = do.call( cbind ,  fclist )
colnames(logfc.degs.all) = types
rownames(logfc.degs.all) = degs
r = cor( logfc.degs.all , use = 'na.or.complete' )

n = length(degs)
o = matrix( 0 , 12 , 12 )
for( i in 1:12 ) {
  for( j in 1:12 ) {
    a = dge$gene[ dge$celltype == types[i] &
		  dge$metaq < 0.05 ]
    b = dge$gene[ dge$celltype == types[j] &
		  dge$metaq < 0.05 ]
    if( length(a) == 0 | length(b) == 0 ) next
    o[i,j] = length( intersect( a , b ) )
  }
}

pdf('logfc.deg.corrplot.pdf')
corrplot( r , type = 'upper' , 
	 method = 'circle' , 
	 diag = F  )
corrplot( r , type = 'lower' ,
	 method = 'number' ,
	 diag = T )
dev.off()


dspn = b3[ b3$celltype == 'dspn' , c('gene','logFC') ]
ispn = b3[ b3$celltype == 'ispn' , c('gene','logFC') ]
espn = b3[ b3$celltype == 'espn' , c('gene','logFC') ]
sst = b3[ b3$celltype == 'sst' , c('gene','logFC') ]
astro = b3[ b3$celltype == 'astro' , c('gene','logFC') ]

logfc = merge( dspn , ispn , by = 1 ) 
logfc = merge( logfc , espn , by = 1 )
colnames(logfc) = c('gene','dspn','ispn','espn')

logfc.sig = logfc[ logfc$gene %in% degs , ]
 
not.pc = grepl( '^Gm[0-9]' , logfc.sig[,1] ) |
     grepl( '[0-9]Rik$' , logfc.sig[,1] ) |
     grepl( '^AI[0-9]' , logfc.sig[,1] )
logfc.pc = logfc.sig[ not.pc == F , ]
tots = rowSums( logfc.pc[,2:4] )
# top.dn = head( logfc.pc[ order( tots ) , 1 ] , 10 )
top.dn = unique( c( 
  head( logfc.pc[ order( logfc.pc[,2] ) , 1 ] , 5 ) ,
  head( logfc.pc[ order( logfc.pc[,3] ) , 1 ] , 5 ) ,
  head( logfc.pc[ order( logfc.pc[,4] ) , 1 ] , 5 ) ) ) 
top.up = unique( c( 
  head( logfc.pc[ order( logfc.pc[,2] , decreasing = T ) , 1 ] , 5 ) ,
  head( logfc.pc[ order( logfc.pc[,3] , decreasing = T ) , 1 ] , 5 ) ,
  head( logfc.pc[ order( logfc.pc[,4] , decreasing = T ) , 1 ] , 5 ) ) ) 
 
# top.up = head( logfc.pc[ order( tots , decreasing = T ) , 1 ] , 10 )
idx.up = which( logfc.sig[,1] %in% top.up ) 
idx.dn = which( logfc.sig[,1] %in% top.dn ) 

cor( logfc.sig[,-1] )

#          dspn      ispn      espn
#dspn 1.0000000 0.8647137 0.5615127
#ispn 0.8647137 1.0000000 0.5361310
#espn 0.5615127 0.5361310 1.0000000

logfc.min = min(logfc.sig[,2:4])
logfc.max = max(logfc.sig[,2:4])
# (-1.81,1.09)

pdf('msn_subtype_logfc_biplots_filt.pdf',
       width=15,height=5)
par( mfrow = c(1,3) ,
     mar = c(5,5,1,1) ,
     bty = 'l' , cex = 1 )
# dSPN vs. iSPN
plot( x = logfc.sig[,2] ,
      y = logfc.sig[,3] ,
      xlab = 'dSPN' ,
      ylab = 'iSPN' ,
      xlim = c(-2,1.25) ,
      ylim = c(-2,1.25) ,
      pch = 19 , col = 'darkgrey' )
abline( h = 0 , v = 0 )
points( x = logfc.sig[c(idx.up,idx.dn),2] ,
        y = logfc.sig[c(idx.up,idx.dn),3] ,
	pch = 19 , col = 'blue' )
text( x = logfc.sig[idx.up,2]+0.05 ,
      y = logfc.sig[idx.up,3] ,
      logfc.sig[idx.up,1] , 
      adj = 0 , cex = 0.7 )
text( x = logfc.sig[idx.dn,2]-0.05 ,
      y = logfc.sig[idx.dn,3] ,
      logfc.sig[idx.dn,1] ,
      adj = 1 , cex = 0.7 )
# dSPN vs. eSPN
plot( x = logfc.sig[,2] ,
      y = logfc.sig[,4] ,
      xlab = 'dSPN' ,
      ylab = 'eSPN' ,
      xlim = c(-2,1.25) ,
      ylim = c(-2,1.25) ,
      pch = 19 , col = 'darkgrey' )
abline( h = 0 , v = 0 )
points( x = logfc.sig[c(idx.up,idx.dn),2] ,
        y = logfc.sig[c(idx.up,idx.dn),4] ,
        pch = 19 , col = 'blue' )
text( x = logfc.sig[idx.up,2]+0.05 ,
      y = logfc.sig[idx.up,4] ,
      logfc.sig[idx.up,1] ,
      adj = 0 , cex = 0.7 )
text( x = logfc.sig[idx.dn,2]-0.05 ,
      y = logfc.sig[idx.dn,4] ,
      logfc.sig[idx.dn,1] ,
      adj = 1 , cex = 0.7 )
# iSPN vs. eSPN
plot( x = logfc.sig[,3] ,
      y = logfc.sig[,4] ,
      xlab = 'iSPN' ,
      ylab = 'eSPN' ,
      xlim = c(-2,1.25) ,
      ylim = c(-2,1.25) ,
      pch = 19 , col = 'darkgrey' )
abline( h = 0 , v = 0 )
points( x = logfc.sig[c(idx.up,idx.dn),3] ,
        y = logfc.sig[c(idx.up,idx.dn),4] ,
        pch = 19 , col = 'blue' )
text( x = logfc.sig[idx.up,3]+0.05 ,
      y = logfc.sig[idx.up,4] ,
      logfc.sig[idx.up,1] ,
      adj = 0 , cex = 0.7 )
text( x = logfc.sig[idx.dn,3]-0.05 ,
      y = logfc.sig[idx.dn,4] ,
      logfc.sig[idx.dn,1] ,
      adj = 1 , cex = 0.7 )
dev.off()








