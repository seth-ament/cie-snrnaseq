# find the genes that overlap between modules and
# the GO terms that are enriched within them

# enrichment results
enr = readRDS('Mod.GOAnno.rds' )

# modules
colors = readRDS('k50.merged.clusters.rds')

# gene sets

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

# get the genes and add them to the table

sig = enr[ which( enr$P < 0.01 ) , ]
n = nrow(sig)
glist = rep(NA,n)
for( i in 1:n ) {
  set.i = sets[[ sig$db[i] ]][[ sig$set[i] ]]
  mod.i = names(colors)[ colors == gsub('M','',sig$mod)[i] ]
  idx = which( toupper(mod.i) %in% set.i )
  isec = mod.i[ idx ]
  glist[i] = paste( isec , collapse=',' )
}

sig$isec.genes = glist

write.table( sig , row.names=F , quote=F , sep= '\t' ,
	     file = 'Mod.GOAnno.P0.01.with_isec_genes.txt' )




