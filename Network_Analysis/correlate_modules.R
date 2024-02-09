library( Seurat )
library( WGCNA )

obj = readRDS('etoh.all_types.w_module_scores.rds' )

mod.stats = read.delim('EtOH.pseudobulk.fisher.modules.txt')
mod.stats$logp = -log10( mod.stats$pval )

up = mod.stats[ mod.stats$dir == 'up' , ]
up = up[ order( up$celltype , up$set ) , ]

down = mod.stats[ mod.stats$dir == 'down' , ]
down = down[ order( down$celltype , down$set ) , ]

dir = down$logp > up$logp

logp = down$logp
logp[ dir == F ] = up$logp[ dir == F ]
logp[ dir == T ] = -1*logp[ dir == T ]

logp.mtx = matrix( logp , ncol = 12 , nrow = 45 )
rownames(logp.mtx) = up$set[1:45]
colnames(logp.mtx) = unique( up$celltype )


# cor SPNs
# this version works better

spn = subset( obj , idents = c('dSPN','iSPN','eSPN') )
mtx = spn@meta.data[,grep('Cluster',colnames(spn@meta.data))]
colnames(mtx) = gsub('Cluster','M',colnames(mtx))
r = cor( mtx )
hc = hclust( as.dist(1-r) , method = 'average' )

pdf('mod.hclust.spns.pdf', height = 3.5 , width = 7 )
plot(hc, hang = -9 )
dev.off()

o = hc$labels[ hc$order ]
logp.mtx = t( logp.mtx[ o , ] )

max = 8
x = logp.mtx
x[ x > max ] = max
x[ x < -1*max ] = -1*max
x = x[ rev(c(3,6,5,11,12,2,1,9,10,4,8,7)) , ]

palette = colorRampPalette(c('blue','white','red'))(299)
breaks = c(seq(-1*max,max,length.out=300))

pdf('mod.logp.heatmap.pdf')
heatmap( x = x ,
	 Rowv = NA , 
	 Colv = NA ,
	 scale = 'none' ,
	 margins = c(30,5) ,
	 col = palette ,
	 breaks = breaks )
dev.off()




