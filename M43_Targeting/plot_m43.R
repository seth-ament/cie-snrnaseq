# plot "rescue" effects on M43 genes

dge = read.delim('GSE149900.voom-limma.DEGs.txt')

colors = readRDS('/local/projects/idea/sament2/etoh/networks/k50.merged.clusters.rds')

m43 = names(colors)[ colors == 43 ]

geno = unique( dge$geno )

genes = intersect( m43 , dge$GeneName )

logFC = list()

for( i in 1:7 ) {
  dge.i = dge[ dge$geno == geno[i] &
	       dge$GeneName %in% m43 , ]
  rownames(dge.i) = dge.i$GeneName
  dge.i = dge.i[ genes , ]
  logFC[[i]] = dge.i$logFC
}
logFC = do.call( cbind , logFC )
rownames(logFC) = genes
colnames(logFC) = paste( geno , '+/-' , sep = '' )

col_palette = colorRampPalette(c('red','white','blue'))(n=299)
breaks = c( seq( -0.5,0,length.out=150 ) ,
	    seq( 0,0.5,length.out=150) )
x = logFC
x[ x > 0.5 ] = 0.5
x[ x < -0.5 ] = -0.5
pdf('m43_rescue.logfc.heatmap.pdf')
heatmap( x = x ,
	 col = col_palette ,
	 breaks = breaks ,
	 scale = 'none' , 
	 cexRow = 0.3 ,
	 mar = c(5,25) )
dev.off()



