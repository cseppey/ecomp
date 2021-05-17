#####
# ecomp boostrap overview
#####

rm(list=ls())

require(abind)
require(venn)
require(xtable)

#---
lar1 <- 3.54
lar1.5 <- 5.51
lar2 <- 7.48

# functions ----
# rnk ---
rnk <- function(x,d){
  apply(x, d, function(y) {
    up <- y[1,] + y[2,]
    down <- y[1,] - y[2,]
    
    rnk <- rank(y[1,], na.last='keep')
    for(i in 1:2){
      for(j in (i+1):3){
        if(is.na(up[i]) | is.na(up[j])){
          next
        }
        cond <- (up[i] > down[j] & down[i] < up[j])
        if(cond){
          rnk[c(i,j)] <- min(rnk[c(i,j)])
        }
      }
    }
    rnk <- rank(rnk, na.last='keep', ties.method=ifelse(max(rnk, na.rm=T) == 2, 'min','max'))
    rnk <- c(3:1)[rnk-diff(c(1, min(rnk, na.rm=T)))]
    
    while(min(rnk, na.rm=T) != 1){
      rnk <- rnk-1
    }
    
    return(rnk)  
  })
}

# nb_of taxa top ---
tax_top <- function(x){
  apply(apply(apply(x, c(2,3), function(x) {
    max <- which(grepl('1', x))
    return(ifelse(length(max) == 1, le[max], 'tie'))
  }), 2, function(y) table(factor(y, levels=c(le,'tie')))), 1,
  function(z) return(c(mean=mean(z),se=sd(z)/sqrt(length(z)))))
}

# table of ranks ---
tb_rnk <- function(x){
  
  tb <- apply(apply(x, 2, function(y){
    yy <- y
    if(all(y == 1)) yy <- rep('same',3)
    if(length(which(y == 1)) == 2) yy <- sub(1,'tie1',yy)
    if(length(which(y == 2)) == 2) yy <- sub(2,'tie2',yy)
    return(yy)
  }), 1, function(z) table(factor(z, levels=c('1','1|2','2','3','tie1','tie2','same'))))
  
  colnames(tb) <- le
  
  return(t(tb))
}

# means of means and se
mms <- function(x){
  apply(x, 2, function(y) apply(y, 1, mean))
}

# download ---
dir_boot <- 'Projets/Ecomp/stat/bootstrap/200821/'
dir_out <- 'Projets/Ecomp/stat/out/200821/'

load('Projets/Ecomp/stat/out/200821/saves/prep.Rdata')

# RETREIVE THE OUTPUT OF THE BOOTSTRAPS ####
lst_out <- list.files(dir_boot, pattern='out', recursive=T, full.names=T)

div <- pvs_nmds <- venns <- lst_omni <- fcts <- perc_ids <- NULL

ind <- 1
for(i in lst_out){
  load(i)
  
  # alpha
  div[['alpha']][['rnk']] <- abind(div[['alpha']][['rnk']], mcls_alpha, along=3)
  
  div[['alpha']][['val']] <- abind(div[['alpha']][['val']], sapply(lst_div$Total, function(x) {
    d <- x$div$div
    return(c(mean=mean(d), se=sd(d)/sqrt(length(d))))
  }), along=3)
  
  # beta
  div[['beta']][['rnk']] <- cbind(div[['beta']][['rnk']], mcl_beta)
  
  div[['beta']][['val']] <- abind(div[['beta']][['val']], sapply(lst_div$Total, function(x) {
    return(c(mean(x$bc),sd(x$bc)/sqrt(length(x$bc))))
  }), along=3)
  
  
  # gamma
  div[['gamma']][['rnk']] <- abind(div[['gamma']][['rnk']], test_gamma, along=3)
  
  div[['gamma']][['val']] <- abind(div[['gamma']][['val']], sapply(lst_div$Total, function(x) x$pool[,c('boot','boot.se')]), along=3)
  
  # nmds
  pvs_nmds <- cbind(pvs_nmds, pv_nmds)
  
  # venn
  venns <- cbind(venns, venn)
  
  # omni OTUs
  lst_omni[[ind]] <- omni_otus
  
  # functions
  fcts <- cbind(fcts, mcls_fct)
  
  # pids
  perc_ids <- abind(perc_ids, as.data.frame(pids), along=3)
  
  #---
  ind <- ind+1
}

# some variables
dn_alpha <- dimnames(div$alpha$rnk)
le <- dn_alpha[[1]]
n_lst <- dn_alpha[[2]]
fct <- levels(ass_sl$fct)

# transfo beta
div$beta$rnk <- aperm(abind(div$beta$rnk, along=3), c(1,3,2))


# GET THE OVERALL ####
lapply(div, function(x) {
  return(list(xtable(tb_rnk(x$rnk[,1,])),
         mms(x$val)))
})

# alpha ----
print('alpha')
# value
mms(div$alpha$val)
# nb taxa
tax_top(div$alpha$rnk[,-1,])

# gamma ----
print('gamma')
# value
mms(div$gamma$val)
# nb taxa
tax_top(div$gamma$rnk[,-1,])

# venn ----
vmse <- round(apply(venns, 1, function(x) return(c(mean=mean(x),se=sd(x)/sqrt(length(x))))))
lvenn <- log(vmse[1,])

# pal <- grey.colors(101, start=0.4, end=1)[102-round((log(vmse[1,])+abs(min(log(vmse[1,]))))/
#                                                       max(log(vmse[1,])+abs(min(log(vmse[1,]))))*100+1)]

alpha <- seq(0.4,1, length.out=100)[round((lvenn-min(lvenn))/(max(lvenn)-min(lvenn)) * 99 +1)]

rgb <- col2rgb(pal_env)

cmbn <- combn(le, 2)
rgb_pair <- apply(cmbn, 2, function(x) rowMeans(rgb[,x]))
colnames(rgb_pair) <- apply(cmbn, 2, function(x) paste(x, collapse='_'))

rgb_ext <- rbind(cbind(rgb,rgb_pair,colMeans(rgb))[,c(3,2,6,1,5,4,7)]/255, alpha)

pal <- apply(rgb_ext, 2, function(x) rgb(x[1],x[2],x[3],x[4]))

#---
cairo_ps(paste(dir_out, 'venn_boot.eps', sep=''), width=lar1.5, height=lar1.5)
par(mar=rep(0, 4))

venn(3, snames=le, cexsn=1)

for(i in seq_along(vmse[1,])){
  polygon(getZones(i,3)[[1]], col=pal[i])
  text(getCentroid(getZones(i,3))[[1]][1], getCentroid(getZones(i,3))[[1]][2], cex=0.75,
       labels=paste(vmse[1,i], '±', round(vmse[2,i],2), '\nn=', c(20,20,40, 20,40,40, 60)[i], sep=''))
}

dev.off()

# list of the omnis ---
for(i in c('order','family','genus')){
  print(i)
  
  tax_tot <- lst_tb <- NULL
  for(j in seq_along(lst_omni)){
    lst_tb[[j]] <- tb <- table(droplevels(taxo_sl[lst_omni[[j]],i]))
    tax_tot <- c(tax_tot, names(tb))
  }
  
  omni <- NULL
  for(j in unique(tax_tot)){
    distro <- sapply(lst_tb, function(x) x[j])
    distro[is.na(distro)] <- 0
    omni <- rbind(omni,c(mean=mean(distro), se=sd(distro)/sqrt(length(distro))))
  }
  row.names(omni) <- unique(tax_tot)
  
  print(xtable(omni[order(omni[,1], decreasing=T),]))
}

# functions ----
grp_fct <- gl(4,3,labels=c(fct))
for(i in fct){
  print(i)
  print(xtable(tb_rnk(fcts[grp_fct == i,])))
}

# pid ----
rnk_pid <- rnk(perc_ids,3)
row.names(rnk_pid) <- le

# table of ranks
print('pid')
print(xtable(tb_rnk(rnk_pid)))

# means of means and se
mms(perc_ids)

#####











