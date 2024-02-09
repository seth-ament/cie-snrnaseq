# create volcano plots highlighting reproducible DEGs

library( dotgen )

disc = read.csv( '/autofs/burnsfs/projects-t3/idea/ewild/etoh/Data_2023/combat8/exprFilter/degs_masterfile_EtOH_B3_pseudobulk_combat8_revised_exprFilter.csv' )
disc$z = zsc( disc$PValue , -1*sign( disc$logFC ) )
colnames(disc)[2] = 'gene'
disc$dir = sign( disc$z )

repl = read.delim('EtOH.pseudobulk.B12.meta-analysis.txt' )
repl$dir = sign( repl$metaz )

B1 = read.csv("/autofs/burnsfs/projects-t3/idea/ewild/etoh/B1_pseudo/B1subset_pseudo_degs/degs_alra_masterfile_EtOH_B1subset_pseudobulk.csv")

disc$gene.type = paste( disc$gene , disc$celltype , disc$dir )
repl$gene.type = paste( repl$gene , repl$celltype , repl$dir )

sig.disc = disc$gene.type[ disc$FDR < 0.05 ]
sig.repl = repl$gene.type[ repl$metaq < 0.05 ]
sig.both = intersect( sig.disc , sig.repl )

sig.genes = unique( disc$gene[ disc$gene.type %in% sig.both ] )

# volcano plots

# dSPN

dspn.disc = disc[ disc$celltype == 'dspn' , ]
logp = -log10( dspn.disc$PValue )
logfc = dspn.disc$logFC
idx = dspn.disc$gene %in% sig.genes
thresh = min( logp[ dspn.disc$FDR < 0.05 ] )

pdf('dSPN.discovery.volcano.highlight_shared_degs.pdf')
par( bty = 'l' , cex = 1.5 )
plot( x = logfc ,
      y = logp ,
      xlab = 'logFC' ,
      ylab = '-log10(P)' ,
      pch = 19 ,
      col = 'darkgrey' )
points( x = logfc[ idx ] ,
        y = logp[ idx ] ,
	col = 'blue' ,
	pch = 19 )
text( x = logfc[ idx ]-0.1 ,
      y = logp[ idx ] ,
      dspn.disc$gene[ idx ] ,
      adj = 1 , cex = 0.5 )
abline( h = thresh , lty = 2 )
dev.off()


dspn.repl = B1[ B1$celltype == 'dspn' , ]
dspn.repl$gene = dspn.repl$X1
logp = -log10( dspn.repl$PValue )
logfc = dspn.repl$logFC
idx = dspn.repl$gene %in% sig.genes
thresh = min( logp[ dspn.repl$FDR < 0.05 ] )

pdf('dSPN.replication.volcano.highlight_shared_degs.pdf' ) 
par( bty = 'l' , cex = 1.5 )
plot( x = logfc ,
      y = logp ,
      xlab = 'logFC' ,
      ylab = '-log10(P)' ,
      pch = 19 ,
      col = 'darkgrey' )
points( x = logfc[ idx ] ,
        y = logp[ idx ] ,
        col = 'blue' ,
        pch = 19 )
text( x = logfc[ idx ]-0.1 ,
      y = logp[ idx ] ,
      dspn.repl$gene[ idx ] ,
      adj = 1 , cex = 0.5 )
abline( h = thresh , lty = 2 )
dev.off()


# iSPN

ispn.disc = disc[ disc$celltype == 'ispn' , ]
logp = -log10( ispn.disc$PValue )
logfc = ispn.disc$logFC
idx = ispn.disc$gene %in% sig.genes
thresh = min( logp[ dspn.disc$FDR < 0.05 ] )

pdf('iSPN.discovery.volcano.highlight_shared_degs.pdf')
par( bty = 'l' , cex = 1.5 )
plot( x = logfc ,
      y = logp ,
      xlab = 'logFC' ,
      ylab = '-log10(P)' ,
      pch = 19 ,
      col = 'darkgrey' )
points( x = logfc[ idx ] ,
        y = logp[ idx ] ,
        col = 'blue' ,
        pch = 19 )
text( x = logfc[ idx ]-0.1 ,
      y = logp[ idx ] ,
      ispn.disc$gene[ idx ] ,
      adj = 1 , cex = 0.5 )
abline( h = thresh , lty = 2 )
dev.off()



ispn.repl = B1[ B1$celltype == 'ispn' , ]
ispn.repl$gene = ispn.repl$X1
logp = -log10( ispn.repl$PValue )
logfc = ispn.repl$logFC
idx = ispn.repl$gene %in% sig.genes
thresh = min( logp[ ispn.repl$FDR < 0.05 ] )

pdf('iSPN.replication.volcano.highlight_shared_degs.pdf')
par( bty = 'l' , cex = 1.5 )
plot( x = logfc ,
      y = logp ,
      xlab = 'logFC' ,
      ylab = '-log10(P)' ,
      pch = 19 ,
      col = 'darkgrey' )
points( x = logfc[ idx ] ,
        y = logp[ idx ] ,
        col = 'blue' ,
        pch = 19 )
text( x = logfc[ idx ]-0.1 ,
      y = logp[ idx ] ,
      ispn.repl$gene[ idx ] ,
      adj = 1 , cex = 0.5 )
abline( h = thresh , lty = 2 )
dev.off()












