library( limma )

dge = read.delim('/local/projects/idea/sament2/etoh/meta/EtOH.mmu.dstr.degs.combined_meta.2023-12-08.txt')

# genesets

f = Sys.glob('genesets/*')
f = setdiff( f , "genesets/README.txt" )
f = f[c(1:3,6,8)]

sets = list()
for( i in 1:length(f) ) {
  sets.i = readLines( f[i] )
  sets.i = strsplit( sets.i , split = '\t' )
  setnames = sapply( 1:length(sets.i) , 
		          function(i) sets.i[[i]][1] )
  sets.i = sapply( 1:length(sets.i) ,
		   function(i) sets.i[[i]][-c(1:2)] )
  names(sets.i) = setnames
  sets[[i]] = 
sets.i
}

names(sets) = gsub('genesets/','',f)
names(sets) = gsub('.txt','',names(sets))

# geneSetTest

types = unique( dge$celltype )
m = length(types)

dbs = names(sets)
n = length(dbs)

df = list()
for( i in 1:m ) {
  cat( types[i] , '\n' )
  dge.i = dge[ dge$celltype == types[i] , ]
  dge.i$score = sign(dge.i$metaz) * -log10(dge.i$metap)
  dge.i$g = toupper( dge.i$gene )
  df.i = list()
  for( j in 1:n ) {
    cat( '...' , dbs[j] , '\n' )
    db = sets[[j]]
    o = length(db)
    p.up = rep(NA,o)
    nGenes = rep(NA,o)
    p.up = rep(NA,o)
    p.dn = rep(NA,o)
    for( k in 1:o ) {
      idx = dge.i$g %in% db[[k]]
      nGenes[k] = length(db[[k]])
      p.up[k] = geneSetTest( 
		statistics = dge.i$score ,
		index = idx ,
		alternative = 'up' )
      p.dn[k] = geneSetTest(
	        statistics = dge.i$score ,
		index = idx ,
		alternative = 'down' )
    }
    df.ij.up = data.frame(
      celltype = types[i] ,
      db = dbs[j] ,
      set = names(db) ,
      dir = 'up' ,
      nGenes ,
      pval = p.up )
    df.ij.dn = data.frame(
      celltype = types[i] ,
      db = dbs[j] ,
      set = names(db) ,
      dir = 'down' ,
      nGenes ,
      pval = p.dn )
    df.ij = rbind( df.ij.up , df.ij.dn )
    df.i[[j]] = df.ij
  }
  names(df.i) = dbs
  dfi2 = do.call( rbind , df.i )
  df[[i]] = dfi2
}
names(df) = types
df2 = do.call( rbind , df )
df2 = df2[ order( df2$pval ) , ]
df2$fdr = p.adjust( df2$pval , method = 'fdr' )

write.table( df2 , row.names=F , sep='\t' , quote=F ,
	     file = 'etoh.dstr.meta.genesettest.txt' )

## top genes in top terms
dge$hsa = toupper( dge$gene )
set = sets[['KEGG_2021_Human']][['Glutamatergic synapse']]
set = sets[['GO_Biological_Process_2023']][['Axon Guidance (GO:0007411)']]
set = sets[['GO_Biological_Process_2023']][["Regulation Of Cell Population Proliferation (GO:0042127)"]]

head( dge[ dge$celltype == 'endo' &
           dge$hsa %in% set , c(1,2,6,7,8) ] , 10 )







