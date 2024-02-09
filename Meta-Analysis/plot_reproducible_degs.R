library( dotgen )
library( viridis )

disc = read.csv( '/autofs/burnsfs/projects-t3/idea/ewild/etoh/Data_2023/combat8/exprFilter/degs_masterfile_EtOH_B3_pseudobulk_combat8_revised_exprFilter.csv' )
disc$z = zsc( disc$PValue , -1*sign( disc$logFC ) )
colnames(disc)[2] = 'gene'
disc$dir = sign( disc$z )

repl = read.delim('EtOH.pseudobulk.B12.meta-analysis.txt' )
repl$dir = sign( repl$metaz )

disc$gene.type = paste( disc$gene , disc$celltype , disc$dir )
repl$gene.type = paste( repl$gene , repl$celltype , repl$dir )

sig.disc = disc$gene.type[ disc$FDR < 0.05 ]
sig.repl = repl$gene.type[ repl$metaq < 0.05 ]
sig.both = intersect( sig.disc , sig.repl )

sig.genes = unique( disc$gene[ disc$gene.type %in% sig.both ] )

types.disc = c('dspn','ispn','espn','pv','sst','chat','astro','oligo',
	       'poly','endo','mural','mg')

n = length(sig.genes)
m = length(types.disc)
z.disc = matrix( 0 , ncol = m , nrow = n )
for( i in 1:n ) {
  for( j in 1:m ) {
    z.ij = disc$z[ disc$gene == sig.genes[i] &
			  disc$celltype == types.disc[j] ]
    if( length(z.ij) == 0 ) next
    z.disc[i,j] = z.ij
  }
}
rownames(z.disc) = sig.genes
colnames(z.disc) = paste( types.disc , 'D' )

                  
types.repl = c('dspn','ispn','espn','pv','sst','astro','oligo',
               'poly','mg')

n = length(sig.genes)
m = length(types.repl)
z.repl = matrix( 0 , ncol = m , nrow = n )
for( i in 1:n ) {
  for( j in 1:m ) {
    z.ij = repl$metaz[ repl$gene == sig.genes[i] &
                          repl$celltype == types.repl[j] ]
    if( length(z.ij) == 0 ) next
    z.repl[i,j] = z.ij
  }
}
rownames(z.repl) = sig.genes
colnames(z.repl) = paste( types.repl , 'R' )

z = cbind( z.disc , NA , z.repl )

x = z
x[ x < -4 ] = -4
x[ x > 4 ] = 4
col_palette = colorRampPalette(c('red','white','blue') )(n=299)
breaks = c( seq( -4 , 0 , length.out = 150 ) ,
	    seq( 0, 4 , length.out = 150 ) )
pdf('sig.gene.z.heatmap.pdf')
heatmap( x = x ,
	 col = col_palette ,
	 breaks = breaks ,
	 Colv = NA ,
	 scale = 'none' ,
         mar = c(5,15) )
dev.off()







