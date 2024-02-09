

colors = readRDS('k50.merged.clusters.rds')

f = Sys.glob('../gsea/genesets/*')
f = setdiff( f , "../gsea/genesets/README.txt" )
f = f[c(1:3,6,7)]

sets = list()
for( i in 1:length(f) ) {
  sets.i = readLines( f[i] )
  sets.i = strsplit( sets.i , split = '\t' )
  setnames = sapply( 1:length(sets.i) ,
                          function(i) sets.i[[i]][1] )
  sets.i = sapply( 1:length(sets.i) ,
                   function(i) sets.i[[i]][-c(1:2)] )
  names(sets.i) = setnames
  sets[[i]] = sets.i
}

names(sets) = gsub('../gsea/genesets/','',f)
names(sets) = gsub('.txt','',names(sets))

u2 = unique(unlist( sets ))

dbs = names(sets)
n = length(dbs)
m = max(colors)
u = intersect( u2 , toupper( names(colors) ))

df = list()
for( i in 1:m ) {
  cat( 'Module' , i , '\n' )
  mod = toupper( names(colors)[ colors == i ] )
  df.i = list()
  for( j in 1:n ) {
    cat( '...' , dbs[j] , '\n' )
    db = sets[[j]]
    o = length(db)
    p = rep(NA,o)
    or = rep(NA,o)
    nGenesSet = rep(NA,o)
    nGenesMod = rep( length(intersect(u,mod)) , o )
    isec = rep(NA,o)
    for( k in 1:o ) {  
      set = intersect( u , db[[k]] )
      nGenesSet[k] = length(set)
      if( length(set) < 5 ) next
      t = table( u %in% mod , u %in% set )
      test = fisher.test( t , alternative = 'greater' )
      p[k] = test$p.value
      or[k] = test$estimate
      isec[k] = t[2,2]
    }
    df.ij = data.frame(
		  mod = paste( 'M' , i , sep = '' ) ,
		  db = names(sets)[j] ,
		  set = names(db) ,
		  nGenesSet ,
		  nGenesMod ,
		  nGenesIsec = isec ,
		  OR = or ,
		  P = p )
    df.i[[j]] = df.ij
  }
  names(df.i) = names(sets)
  df[[i]] = df.i
}
names(df) = paste( 'M' , 1:max(colors) , sep = '' )


df2 = list()
for( i in 1:length(df) ) {
  df2[[i]] = do.call( rbind , df[[i]] )
}
df2 = do.call( rbind , df2 )
df2 = df2[ order( df2$P ) , ]
df2$FDR = p.adjust( df2$P , method = 'fdr' )

saveRDS( df2 , file = 'Mod.GOAnno.rds' )

write.table( df2 , sep='\t' , quote=F , row.names=F ,
	     file = 'Mod.GOAnno.txt' )

