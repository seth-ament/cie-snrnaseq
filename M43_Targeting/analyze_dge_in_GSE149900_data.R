library( GEOquery )
library( edgeR )

### get modules

colors = readRDS('/local/projects/idea/sament2/etoh/networks/k50.merged.clusters.rds')
MM = readRDS('/local/projects/idea/sament2/etoh/networks/module_membership.kAE.rds')


getGEOSuppFiles( 'GSE149900' )
gse = getGEO('GSE149900')
eset = gse[[1]]
meta0 = pData( eset )
geno = meta0$'genotype:ch1'
sex = meta0$'Sex:ch1'
gene = gsub( '\\+\\/\\-' , '' , geno )
gene[ gene == 'Wild type' ] = 'wild.type'
meta = data.frame(
  sample = rownames( meta0 ) ,
  geno , gene , sex )

counts = read.csv('/autofs/projects-t2/idea/sament2/etoh/druggable/pde10a/GSE149900/GSE149900_rawCounts_WtHet_Striatum_6m.csv.gz')
rownames(counts) = counts$GeneID

anno = counts[,c('GeneID','GeneName','Description','Chromosome','Strand')]

counts = counts[,-c(835:839)]

fpkm = read.csv('/autofs/projects-t2/idea/sament2/etoh/druggable/pde10a/GSE149900/GSE149900_FPKM_WtHet_Striatum_6m.csv.gz')

# which genes from M43 have data?

i=43
  mod = names(colors)[ colors == i ]
  MM.mod = MM[ mod , ]
  ranked = MM.mod[ order( MM.mod[,paste('M',i,sep='')] ,
                          decreasing = T ) , 
			  paste('M',i,sep='') ]
  m43.genes = intersect( names(ranked) , meta$gene )

# DEGs for each gene in M43

y = DGEList( counts = counts )
keep = filterByExpr( y )
y = y[ keep , ]
y = calcNormFactors( y )
y = estimateDisp( y )
design = model.matrix( ~ 0 + gene + sex , data = meta )
colnames(design) = gsub('gene','',colnames(design))
contrasts = makeContrasts(
  Pde10a = Pde10a - wild.type ,
  Foxp1 = Foxp1 - wild.type ,
  Rgs9 = Rgs9 - wild.type ,
  Itpr1 = Itpr1 - wild.type ,
  Rarb = Rarb - wild.type ,
  Bcl11b = Bcl11b - wild.type ,
  Zswim6 = Zswim6 - wild.type ,
  levels = design )

voom = voom( y , design )
fit = lmFit( voom , design )
fit = contrasts.fit( fit , contrasts )
fit = eBayes( fit )
res = list()
for( i in 1:7 ) {
  outp = topTable( fit , coef = i , n = Inf )
  outp$geno = m43.genes[i]
  outp = merge( anno , outp , by.x = 1 , by.y = 0 )
  res[[i]] = outp
}
res = do.call( rbind , res )
res = merge( anno , Pde10a , by.x = 1 , by.y = 0 )

# test for overlaps

u = names(colors)
mod = names(colors)[ colors == 43 ]
nUpInMod = nDownInMod = nUp = nDown = PValue.UpInMod = 
	PValue.DownInMod = OR.UpInMod = 
	OR.DownInMod = PValue.UpDownInMod = OR.UpDownInMod = 
	  rep(NA,7)

rescue = list()
for( i in 1:7 ) {
  up = res$GeneName[ res$geno == m43.genes[i] & 
		     res$logFC > 0 & res$adj.P.Val < 0.05 ]
  dn = res$GeneName[ res$geno == m43.genes[i] & 
		     res$logFC < 0 & res$adj.P.Val < 0.05 ]
  t.dn = table( u %in% mod , u %in% dn )
  test.dn = fisher.test( t.dn , alternative = 'greater' )
  t.up = table( u %in% mod , u %in% up )
  test.up = fisher.test( t.up , alternative = 'greater' )
  df = data.frame(
    genotype = paste( m43.genes[i] , '+/-' , sep = '' ) ,
    Direction = c('Up','Down') ,
    nDEG = c(length(up),length(dn)) ,
    nDEGinMod = c( t.up[2,2] , t.dn[2,2] ) ,
    P = c( test.up$p.value , test.dn$p.value ) ,
    OR = c( test.up$estimate , test.dn$estimate ) )
  rescue[[i]] = df
}
rescue = do.call( rbind , rescue )
rescue = rescue[ order( rescue$P ) , ]

write.table( res , quote=F , sep='\t' , row.names=F ,
	     file = 'GSE149900.voom-limma.DEGs.txt' )

write.table( rescue , quote=F , sep='\t' , row.names=F ,
	     file = 'GSE149900vsM43.overlap.table.txt' )

saveRDS( meta , file = 'GSE149900.sample_metadata.rds' )
saveRDS( counts , file = 'GSE149900.read_counts.rds' )


