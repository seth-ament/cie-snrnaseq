library( Seurat )
library( lme4 )

obj = readRDS('etoh_final_obj_wsex_wgroup.rds')

colors = readRDS('k50.merged.clusters.rds') 

mods = list()
for( i in 1:max(colors) ) {
  mods[[i]] = names(colors)[ colors == i ]
}
names(mods) = paste('M',1:max(colors),sep='')

obj = AddModuleScore( obj , mods )
batch = rep(NA,ncol(obj))
batch[ obj$sample_number %in% c(1,5) ] = 1
batch[ obj$sample_number %in% c(2,6) ] = 2
batch[ obj$sample_number %in% c(3,7) ] = 3
batch[ obj$sample_number %in% c(4,8) ] = 4
obj$batch = batch

saveRDS( obj , file = 'etoh.all_types.w_module_scores.rds' )

spn = subset( obj , idents = c('dSPN','iSPN','eSPN') )

pdf('EtOH.Module.VlnPlots.pdf',width=15,height=7)
VlnPlot( obj , features = paste('Cluster',1:45,sep='') ,
	 combine = F , split = T , split.by = 'trt' , pt.size = 0 )
dev.off()

types = levels(obj)
df = list()
for( i in 1:12 ) {
  cat( types[i] , '\n' )
  sub = subset( obj , idents = types[i] )
  m = sub@meta.data[,c('trt','batch','sex')]
  m$sample = paste( m$trt , m$sex , m$batch )
  pval = rep(NA,45)
  beta = rep(NA,45)
  for( j in 1:45 ) {
    y = sub@meta.data[,paste('Cluster',j,sep='')]
    fit = lmer( y ~ trt + sex + (1|sample) , data = m )
    res = drop1( fit , ~. , test = 'Chisq' )
    pval[j] = res[2,4]
    beta[j] = coef(fit)$sample[1,2]
  }
  df[[i]] = data.frame(
	   celltype = types[i] ,
	   module = paste('M',1:45,sep='') ,
	   beta = beta ,
	   pval = pval )
}
df = do.call( rbind , df )

df = df[ order( df$pval ) , ]
df$fdr = p.adjust( df$pval , method = 'fdr' )

mod.pseudobulk = list()
for( i in 1:length(types) ) {
  sub = subset( obj , idents = types[i] )
  m = sub@meta.data[,c('trt','batch','sex')]
  m$sample = paste( m$trt , m$sex , m$batch )
  m2 = data.frame(
    trt = c( rep('air',8) , rep('etoh',8) ) ,
    batch = rep(c(1:4),4) ,
    sex = rep(c(rep('F',4),rep('M',4)),2) )
  sample_means = list()
  for( j in 1:45 ) {
    y = sub@meta.data[,paste('Cluster',j,sep='')]
    sample_means[[j]] = as.vector( by( y , m$sample , mean ) )
    # fit = lm( sample_means ~ trt + sex + batch , data = m2 )
    # res = drop1( fit , ~. , test = 'F' )
  }
  sample_means = do.call( rbind , sample_means )
  colnames(sample_means) = paste( m2$trt , m2$sex , m2$batch , sep = '_' )
  rownames(sample_means) = paste('M',1:45,sep='')
  mod.pseudobulk[[i]] = sample_means
}
names(mod.pseudobulk) = types[1:12]

m2 = data.frame(
    trt = c( rep('air',8) , rep('etoh',8) ) ,
    batch = rep(c(1:4),4) ,
    sex = rep(c(rep('F',4),rep('M',4)),2) ,
    xpos = c( rep(0.25,8) , rep(0.5,8) ) )
pdf('module_scatter.pdf',width=10,height=8)
for( j in 1:45 ) {
  par( mfrow = c(3,4) , mar = c(2,2,2,1) , 
       bty = 'l' , cex = 1.2 , oma = c(0,2,2,0) )
  for( i in c(4,2,6,10:12,1,3,5,7:9) ) {
    y = mod.pseudobulk[[i]][j,]
    plot( x = jitter(m2$xpos,1) ,
      y = y , xlim = c(0.1,0.65) ,
      pch = 19 , xaxt = 'n' )
    axis( side = 1 , at = c(0.25,0.5) , labels = c('Air','EtOH') )
    mtext( side = 3 , line = 0.5 , cex = 1.2 , types[i] )
    if( i == 1 ) mtext( outer = T , cex = 2 , paste( 'M' , j , sep = '' ) )
    if( i == 1 ) mtext( outer = T , line = 0.5 , cex = 1.2 , 
		       side = 2 ,  'Average Expression' )
  }
}
dev.off()

pdf('module_boxplot.pdf',width=10,height=8)
for( j in 1:45 ) {
  par( mfrow = c(3,4) , mar = c(2,2,2,1) ,
       bty = 'l' , cex = 1.2 , oma = c(0,2,2,0) )
  for( i in c(4,2,6,10:12,1,3,5,7:9) ) {
    y = mod.pseudobulk[[i]][j,]
    boxplot( y ~ m2$trt , xlab = '' , ylab = '' )
    # axis( side = 1 , at = c(0.25,0.5) , labels = c('Air','EtOH') )
    mtext( side = 3 , line = 0.5 , cex = 1.2 , types[i] )
    if( i == 1 ) mtext( outer = T , cex = 2 , paste( 'M' , j , sep = '' ) )
    if( i == 1 ) mtext( outer = T , line = 0.5 , cex = 1.2 ,
                       side = 2 ,  'Average Expression' )
  }
} 
dev.off()






