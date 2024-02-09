library( viridis )

res = read.delim('etoh.dstr.meta.genesettest.with_fdr.txt')

dge = read.delim('/local/projects/idea/sament2/etoh/meta/EtOH.mmu.dstr.degs.combined_meta.2023-12-08.txt')

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
  sets[[i]] = sets.i
}

names(sets) = gsub('genesets/','',f)
names(sets) = gsub('.txt','',names(sets))

# get the significant sets

sig = res[ res$fdr < 0.05 , ]

sig.sets = list()
for( i in 1:nrow(sig) ) {
  type.i = sig$celltype[i]
  db.i = sig$db[i]
  term.i = sig$set[i]
  dir.i = -1
  if( sig$dir[i] == 'up' ) dir.i = 1
  set.i = sets[[db.i]][[term.i]]
  dge.i = dge$gene[ dge$celltype == sig$celltype[i] &
		    sign(dge$metaz) == dir.i ]
  sig.sets[[i]] = intersect( set.i , toupper(dge.i) )
}

# IRR analysis -- drop terms that overlap more significant ones

u = toupper( unique( dge$gene ) )
thresh = 0.5
drop = vector()
for( i in 1:nrow(sig) ) {
  cat( i )
  set.i = sig.sets[[i]]
  if( i == 1 ) next
  totest = setdiff( 1:(i-1) , drop )
  k2.i = rep(NA,length(totest))
  for( j in 1:length(totest) ) {
    set.j = sig.sets[[ totest[j] ]]
    r1 = as.numeric( u %in% set.i )
    r2 = as.numeric( u %in% set.j )
    ratings = cbind( r1 , r2 )
    k2.i[j] = kappa2( ratings )$value
  }
  if( any( k2.i > thresh ) ) {
	drop = append( drop , i )
	cat( 'dropped\n' )
  } 
  if( all( k2.i <= thresh ) )  cat('\n')
}

keep = 1:217 %in% drop == F & sig$fdr < 0.01

keep.terms = unique( sig[ keep , c('db','set') ] )

n = nrow( keep.terms )
types = c('dspn','ispn','espn','pv','sst','chat',
	  'astro','oligo','poly','endo','mural','mg')
m = length(types)
pmat = matrix( NA , nrow = n , ncol = m )
for( i in 1:n ) {
  for( j in 1:m ) {
    up = res$pval[ res$celltype == types[j] &
		   res$db == keep.terms$db[i] &
		   res$set == keep.terms$set[i] &
		   res$dir == 'up' ]
    dn = res$pval[ res$celltype == types[j] &
                   res$db == keep.terms$db[i] &
                   res$set == keep.terms$set[i] &
                   res$dir == 'down' ]
    if( dn < up ) pmat[i,j] = log10(dn)
    if( up < dn ) pmat[i,j] = -log10(up)
  }
}
rownames(pmat) = keep.terms$set
colnames(pmat) = types

col_palette = colorRampPalette(c('red','white','blue'))(n=299)
breaks = c( seq(-4,0,length.out=150) , seq(0,4,length.out=150) )
x = pmat
x[ x > 4 ] = 4
x[ x < -4 ] = -4
pdf('go_enrich_heatmap.pdf',width=8)
heatmap( x = x ,
	 col = col_palette, 
	 breaks = breaks ,
	 scale = 'none' ,
	 Colv = NA ,
	 mar = c(5,25) )
dev.off()






