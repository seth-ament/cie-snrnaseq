library( biomaRt )

# human-mouse orthologs

mart = useMart( 'ensembl' )
mart = useDataset( mart = mart , dataset = 'hsapiens_gene_ensembl' )
syn = getBM( mart = mart ,
	      attributes = c('ensembl_gene_id','hgnc_symbol') )
orth = getBM( mart = mart ,
	      attributes = c('ensembl_gene_id','mmusculus_homolog_associated_gene_name' ) )
orth = merge( syn , orth )

# load relevant files

MM = readRDS('/local/projects/idea/sament2/etoh/networks/module_membership.kAE.rds' )

mods = readRDS('/local/projects/idea/sament2/etoh/networks/k50.merged.clusters.rds')

druggable = read.delim('druggable_genes_table.txt')

f = '/local/projects/idea/sament2/etoh/gsea/genesets/LINCS_L1000_Chem_Pert_Consensus_Sigs.txt'

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
names(sets) = 'LINCS'


m43 = names(mods)[ mods == 43 ]
m43 = MM[m43,'M43']
m43 = data.frame( gene = names(m43) , kAE = m43 )
m43 = merge( orth , m43 , by.x = 3 , by.y = 1 )
m43 = merge( m43 , druggable , by.x = 2 , by.y = 1 )
m43 = m43[ order( m43$kAE , decreasing = T ) , ]
m43.drug = m43

m29 = names(mods)[ mods == 29 ]
m29 = MM[m29,'M29']
m29 = data.frame( gene = names(m29) , kAE = m29 )
m29 = merge( orth , m29 , by.x = 3 , by.y = 1 )
m29 = merge( m29 , druggable , by.x = 2 , by.y = 1 )
m29 = m29[ order( m29$kAE , decreasing = T ) , ]
m29.drug = m29

write.table( m29.drug , quote=F , row.names=F , sep='\t' ,
	     file = 'm29.druggable_genes.txt' )

write.table( m43.drug , quote=F , row.names=F , sep='\t' ,
             file = 'm43.druggable_genes.txt' )



#

m29 = names(mods)[ mods == 29 ]
m29 = MM[m29,'M29']
m29 = data.frame( gene = names(m29) , kAE = m29 )
m29 = merge( orth , m29 , by.x = 3 , by.y = 1 )
m29 = m29$hgnc_symbol

n = length(sets$LINCS)
or = p = o = nset = rep(NA,n)
u = unique( orth$hgnc_symbol[ orth$mmusculus_homolog_associated_gene_name %in% names(mods) ] )
for( i in 1:n ) {
  set = sets$LINCS[[i]]
  nset[i] = length(set)
  t = table( u %in% set , u %in% m29 )
  test = fisher.test( t , alternative = 'greater' )
  p[i] = test$p.value
  or[i] = test$estimate
  o[i] = t[2,2]
}
m29.lincs = data.frame(
  mod = 'M29' ,
  set = names(sets[[1]]) ,
  nGenesSet = nset ,
  nGenesMod = length(m29) ,
  nOverlap = o ,
  OR = or ,
  P = p )
m29.lincs = m29.lincs[ order( m29.lincs$P ) , ]
m29.lincs$FDR = p.adjust( m29.lincs$P , method = 'fdr' )

#  m43
m43 = names(mods)[ mods == 43 ]
m43 = MM[m43,'M43']
m43 = data.frame( gene = names(m43) , kAE = m43 )
m43 = merge( orth , m43 , by.x = 3 , by.y = 1 )
m43 = m43$hgnc_symbol

n = length(sets$LINCS)
or = p = o = nset = rep(NA,n)
u = unique( orth$hgnc_symbol[ orth$mmusculus_homolog_associated_gene_name %in% names(mods) ] )
for( i in 1:n ) {
  set = sets$LINCS[[i]]
  nset[i] = length(set)
  t = table( u %in% set , u %in% m43 )
  test = fisher.test( t , alternative = 'greater' )
  p[i] = test$p.value
  or[i] = test$estimate
  o[i] = t[2,2]
}
m43.lincs = data.frame(
  mod = 'M43' ,
  set = names(sets[[1]]) ,
  nGenesSet = nset ,
  nGenesMod = length(m43) ,
  nOverlap = o ,
  OR = or ,
  P = p )
m43.lincs = m43.lincs[ order( m43.lincs$P ) , ]
m43.lincs$FDR = p.adjust( m43.lincs$P , method = 'fdr' )

write.table( m29.lincs , row.names=F , quote=F , sep='\t' ,
	    file = 'M29.LINCS.enrich.txt' )
 
write.table( m43.lincs , row.names=F , quote=F , sep='\t' ,
            file = 'M43.LINCS.enrich.txt' )





