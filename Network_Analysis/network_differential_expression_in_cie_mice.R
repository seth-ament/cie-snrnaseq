library( limma )

dge = read.delim('/local/projects/idea/sament2/etoh/meta/EtOH.mmu.dstr.degs.combined_meta.min_100_cells_per_batch.2023-12-11.txt')

colors = readRDS('k50.merged.clusters.rds')

types = unique(dge$celltype)
n = length(types)
m = max(colors)
df = list()
for( i in 1:n ) { 
  dge.i = dge[ dge$celltype == types[i] , ]
  dge.i$score = sign(dge.i$metaz) * -log10(dge.i$metap)
  p.dn = rep(NA,m)
  p.up = rep(NA,m)
  nGenes = rep(NA,m)
  for( k in 1:m ) {
      set = names(colors)[colors==k]
      idx = dge.i$gene %in% set
      nGenes[k] = length(set)
      p.up[k] = geneSetTest(
                statistics = dge.i$score ,
                index = idx ,
                alternative = 'up' )
      p.dn[k] = geneSetTest(
                statistics = dge.i$score ,
                index = idx ,
                alternative = 'down' )
  }
      df.i.up = data.frame(
        celltype = types[i] ,
        set = paste('M',1:m,sep='') ,
        dir = 'up' ,
        nGenes ,
        pval = p.up )
      df.i.dn = data.frame(
        celltype = types[i] ,
        set =  paste('M',1:m,sep=''),
        dir = 'down' ,
        nGenes ,
        pval = p.dn )
      df.i = rbind( df.i.up , df.i.dn )
  df[[i]] = df.i
}
df = do.call( rbind , df )

df$fdr = p.adjust( df$pval , method = 'fdr' )
df$bonf = p.adjust( df$pval , method = 'bonf' )

df = df[ order( df$pval ) , ]

write.table( df , quote=F , row.names=F , sep='\t' ,
	     file = 'EtOH.pseudobulk.geneSetTest.modules.txt' )

### fisher ###

types = unique(dge$celltype)
n = length(types)
m = max(colors)
df = list()
for( i in 1:n ) {
  dge.i = dge[ dge$celltype == types[i] , ]
  up = dge.i$gene[ dge.i$metap < 0.01 & dge.i$metaz > 0 ]
  dn = dge.i$gene[ dge.i$metap < 0.01 & dge.i$metaz < 0 ]
  u = names(colors)
  p.dn = rep(NA,m)
  p.up = rep(NA,m)
  or.dn = rep(NA,m)
  or.up = rep(NA,m)
  o.dn = rep(NA,m)
  o.up = rep(NA,m)
  nSetGenes = rep(NA,m)
  nUp = rep(length(up),m) 
  nDown = rep(length(dn),m)
  for( k in 1:m ) {
      set = names(colors)[colors==k]
      nSetGenes[k] = length(set)
      t.up = table( u %in% up , u %in% set )
      test.up = fisher.test( t.up , alternative = 'greater' )
      p.up[k] = test.up$p.value
      or.up[k] = test.up$estimate
      o.up[k] = t.up[2,2]
      t.dn = table( u %in% dn , u %in% set )
      test.dn = fisher.test( t.dn , alternative = 'greater' )
      p.dn[k] = test.dn$p.value
      or.dn[k] = test.dn$estimate
      o.dn[k] = t.dn[2,2]
  }
      df.i.up = data.frame(
        celltype = types[i] ,
        set = paste('M',1:m,sep='') ,
        dir = 'up' ,
	nDEGs = nUp ,
        nSetGenes ,
	nOverlap = o.up ,
	OR = or.up ,
        pval = p.up )
      df.i.dn = data.frame(
        celltype = types[i] ,
        set = paste('M',1:m,sep='') ,
        dir = 'down' ,
        nDEGs = nDown ,
        nSetGenes ,
        nOverlap = o.dn ,
        OR = or.dn ,
        pval = p.dn )
  df.i = rbind( df.i.up , df.i.dn )
  df[[i]] = df.i
}
df = do.call( rbind , df )
df = df[ order( df$pval ) , ]
df$fdr = p.adjust( df$pval , method = 'fdr' )
df$bonf = p.adjust( df$pval , method = 'bonf' )

write.table( df , quote=F , sep='\t' , row.names=F ,
	     file = 'EtOH.pseudobulk.fisher.modules.txt' )




