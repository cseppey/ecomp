#####
# ecomp prep
#####

rm(list=ls())

require(PMCMR)
require(venn)
require(multcompView)
require(multcomp)
require(plotrix)
require(RColorBrewer)
require(dendextend)
require(foreach)
require(doSNOW)

setwd('~')

source('bin/src/my_prog/R/pie_taxo.r')
source('bin/src/my_prog/R/pie_taxo_single.r')
source('bin/src/my_prog/R/legend_pie_taxo.r')

argv <- commandArgs(trailingOnly=T)

seed <- as.integer(argv[1])
raref <- as.integer(argv[2])
th <- as.integer(argv[3])
recalc <- as.logical(argv[4])

# ##########
# seed <- 1
# raref <- 10000
# th <- 4
# recalc <- F
# ##########

print(paste('seed',seed))
print(paste('raref',raref))
print(paste('threads',th))
print(paste('recalc',recalc))

if(seed < 10){
  boot <- paste0('00',seed)
} else if (seed < 100){
  boot <- paste0('0',seed)
} else {
  boot <- as.character(seed)
}

print(paste('%%%%%%%%', boot, '%%%%%%%%'))

# param variable ####
lar1 <- 3.54
lar1.5 <- 5.51
lar2 <- 7.48

#---
dom_tax <- c(
  'Stramenopiles*',
    'Oomycota','Hyphochytriomyceta',
    'Bicoecea',
    'Labyrinthulomycetes',
    
    'MAST',

    'Synurophyceae',
    'Chrysophyceae',
    'Xanthophyceae',
    'Eustigmatophyceae',
    'Bacillariophyta','Bolidophyceae',
    'Dictyochophyceae','Pelagophyceae',
  
  'Alveolata*',
    'Dinoflagellata*',
      'Syndiniales','Dinophyceae',
      'Apicomplexa',

    'Ciliophora*',
      'Colpodea','Litostomatea','Oligohymenophorea','Phyllopharyngea','Spirotrichea',
  
  'Rhizaria*',
    'Cercozoa*',
      'Granofilosea','Imbricatea','Sarcomonadea','Thecofilosea','Phaeodarea',
      'Vampyrellida', 'Plasmodiophorida',
    'Foraminifera*',
      'Globothalamea','Monothalamids',
      'Radiolaria',
  
  'Telonemia',
####
  'Centroheliozoa',
  'Haptophyta',
####
  'Cryptophyta*',
    'Goniomonadales','Cryptomonadales',
  'Katablepharidaceae',
####
  'Archaeplastida*',
    'Rhodophyta*',

    'Chlorophyta*',
      'Chlorophyceae','Trebouxiophyceae',
      'Mamiellophyceae',
      'Prasino',
    
    'Streptophyta*',
      'Zygnemophyceae','Klebsormidiophyceae',
      'Embryophyceae*',
####
  'Picozoa',
####
####
  'Apusozoa',
  
  'Opisthokonta*',
    'Ichthyosporea','Choanoflagellida','Metazoa*',
    'Fungi*',
    'Nucleariidea',
  
  'Amoebozoa*',
    'Variosea','Lobosa',
    'Mycetozoa',
  
  'Excavata*',
    'Discoba*',
      'Heterolobosea',
      'Kinetoplastida',
      'Diplonemea','Euglenida'
  
) # est-ce qu'on vérifie les seuil pid?

main_taxa <- grep('*', dom_tax, invert=T, fixed=T, value=T)


# DOWNLOAD ####
dir_base <- 'Projets/Ecomp/stat/'
# dir_in   <- paste0(dir_base, 'in/')
dir_save <- paste0(dir_base, 'out/200821/saves/')

dir_boot <- paste0(dir_base, 'bootstrap/200821/bootstrap_', boot, '/')
dir_bave <- paste0(dir_boot, 'saves/')
dir_bout <- paste0(dir_boot, 'out/')

lapply(list(dir_bave,dir_bout), function(x) dir.create(x, F, T))

# community rarefaction
set.seed(seed)

file <- paste0(dir_bave, 'boot_', boot, '.Rdata')
if(file.exists(file) & recalc == F){
  load(file)
} else {
  
  load(paste0(dir_save, 'prep.Rdata'))
  # #####
  # load(file)
  # #####
  
  # rrarefy
  print('rrarefy')

  mr <- mr_sl[rowSums(mr_sl) > raref,]
  e <- env[row.names(mr),]
  mr <- mr[unlist(tapply(row.names(e), list(e$env), function(x) sample(x, 20))),]

  mr_boot <- as.data.frame(rrarefy(mr, raref))
  mr_boot <- droplevels(mr_boot[,colSums(mr_boot) != 0])

  env_boot <- e[row.names(mr_boot),]
  ass_boot  <- droplevels( ass_sl[names(mr_boot),])
  taxo_boot <- droplevels(taxo_sl[names(mr_boot),])
  
  # save
  print('save')
  save(env_boot, mr_boot, ass_boot, taxo_boot, file=file)
  
  exit()
}

# sort the bootstrap according to taxonomy and environments
ass_boot$taxo <- gsub('Hacrobia','Eukaryota_X',ass_boot$taxo)

env_boot$env <- factor(env_boot$env, levels=levels(env_boot$env)[c(2,1,3)])

ord_otu <- order(ass_boot$taxo)
ord_smp <- order(env_boot$env)

env_sort <- env_boot[ord_smp,]
mr_sort <- mr_boot[ord_smp,ord_otu]
ass_sort <- ass_boot[ord_otu,]
taxo_sort <- taxo_boot[ord_otu,]

le <- levels(env_sort$env)

tab_ech <- table(env_sort$env)

# palette
pal_env <- rgb(c(0,0,1),c(0,196,102)/255, c(1,1,0))
# pal_env <- brewer.pal(3,'Set1')[c(2,3,1)]
names(pal_env) <- le


# DIVERSITY AND ABUNDANCE ####
# Diversity calculation ----
print('##### Diversity #####')

# calculation of metrics
n_lst <- c('Total', main_taxa)

file <- paste0(dir_bave, 'lst_div_', boot, '.Rdata')
if(file.exists(file) & recalc == F){
  load(file)
} else {
  
  # cluster
  th2 <- max(c(th/2,1))
  cl2 <- makeSOCKcluster(th2)
  registerDoSNOW(cl2)
  clusterEvalQ(cl2, library(doSNOW))
  clusterEvalQ(cl2, library(foreach))
  clusterEvalQ(cl2, library(vegan))
  
  lst_div <- foreach(i = n_lst, .verbose=T) %dopar% { # on the taxa
    
    print(paste('###',i,'###'))
    
    ind_otu <- 1:ncol(mr_sort)
    if(i != 'Total'){
      ind_otu <- grep(i, ass_sort$taxo)
    }
    
    # cluster
    th22 <- max(c(1,th/2))
    cl22 <- makeSOCKcluster(th22)
    registerDoSNOW(cl22)
    clusterEvalQ(cl22, library(vegan))
    
    lst <- foreach(j = le, .verbose=T) %dopar% { # on each ecosystem
      
      print(paste('---',j,'---'))
      
      ind_smp <- env_sort$env == j
      
      mr <- mr_sort[ind_smp,ind_otu]
      ind0 <- colSums(mr) != 0
      
      if(nrow(mr) > 0 & ncol(as.data.frame(mr[,ind0])) > 1){
        
        mr <- as.data.frame(mr[,ind0])
        
        #---
        set.seed(0)
        
        # alpha ---
        # abundance
        print('abds')
        abds <- rowSums(mr)
        
        # richness
        print('richness')
        rich <- specnumber(mr)
        
        # shannon
        print('shannon')
        div <- diversity(mr)
        
        # beta ---
        # bray-curtis
        print('BC')
        bc <- vegdist(mr)
        bc[is.na(bc)] <- 0
        
        # gamma ---
        # specaccum
        print('specaccum')
        accum <- specaccum(mr)
        
        # specpool
        print('specpool')
        pool <- specpool(mr)
        
      } else {
        abds <- rich <- div <- rep(0, nrow(mr_sort[ind_smp,]))
        
        bc <- NULL
        
        accum <- data.frame(richness=abds,sd=abds)
        
        pool <- NULL
        
      }
      
      gc()
      
      # fill the list
      l <-list(div=data.frame(accuR=accum$richness, accuSD=accum$sd, abds=abds, rich=rich, div=div),
               pool=pool,
               bc=bc)
      
      return(l)
      
    }
    
    stopCluster(cl22)
    
    gc()
    
    names(lst) <- le
    
    return(lst)
    
  }
  
  stopCluster(cl2)
  
  gc()
  
  names(lst_div) <- n_lst

  save(lst_div, file=file)
}

# diversity Total figure
tot <- lst_div$Total

#---
# pdf(paste0(dir_bout, 'diversity.pdf'), width=lar2, height=lar2)
cairo_ps(paste0(dir_bout, 'diversity.eps'), width=lar2, height=lar2)
par(mfrow=c(2,2), mar=c(5,4,2,1), mgp=c(2.5,1,0))

offset <- 0.07
#

# Gamma ----

plot(NA, xlim=c(0,max(sapply(tot, function(x) length(x$div$accuR)))), 
     ylim=range(sapply(tot, function(x) {return(c(min(apply(x$div[,2:1], 1, diff)), 
                                                  sum(x$pool$boot, x$pool$boot.se)))})),
     xlab='number of samples', ylab='number of OTU')
usr <- par('usr')

text(usr[1]-diff(usr[1:2])*offset, usr[4]+diff(usr[3:4])*offset, 'A', xpd=NA, font=2)

for(i in le){
  
  col <- pal_env[[i]]
  
  # specaccum
  accu <- tot[[i]]$div[,c('accuSD','accuR')]
  
  lines(accu$accuR, col=col)
  
  seg <- apply(as.matrix(accu), 1, function(x) return(c(diff(x), sum(x))))  
  segments(seq_along(accu$accuSD), seg[1,], seq_along(accu$accuSD), seg[2,], col=col)
  
  # specpool
  bse <- unlist(tot[[i]]$pool[,c('boot','boot.se')])
  abline(h=bse[1], col=col)
  abline(h=bse[1]+c(1,-1)*bse[2], col=col, lty=2)
  
}

# Alpha ----
alp <- lapply(tot, function(x) x$div$div)

pvs <- posthoc.kruskal.nemenyi.test(alp)$p.value
n_pvs <- apply(expand.grid(le, le), 1, function(x) paste(x, collapse='-'))
pvs <- c(cbind(rbind(rep(NA,ncol(pvs)),pvs),rep(NA,nrow(pvs)+1)))
names(pvs) <- n_pvs
pvs <- pvs[is.na(pvs) == F]

e <- factor(unlist(lapply(seq_along(alp), function(x) rep(le[x],length(alp[[x]])))), levels=le)
alp <- cbind.data.frame(div=unlist(alp), env=e)

mcl_alpha <- multcompLetters2(div~env, pvs, alp, Letters=1:4)$Letters[le]
if(max(nchar(mcl_alpha)) > 1){
  mcl_alpha[nchar(mcl_alpha) == 2] <- sapply(strsplit(mcl_alpha[nchar(mcl_alpha) == 2], ''), function(x) paste(x, collapse='|'))
}

#---
boxplot(div~env, alp, ylab='H', ylim=range(alp$div)*c(1, 1.2), col=pal_env, xaxt='n')
usr <- par('usr')

axis(1, at=seq_along(le), labels=F)
text(1:length(le), usr[3]-diff(usr[3:4])*0.15, labels=le, xpd=NA, cex=0.75, srt=45)

text(usr[1]-diff(usr[1:2])*offset, usr[4]+diff(usr[3:4])*offset, 'B', xpd=NA, font=2)

text(1:length(le), usr[4]-diff(usr[3:4])*0.1, mcl_alpha, offset=0, pos=3)

# Beta ----
bet <- as.data.frame(sapply(tot, function(x) c(x$bc)))

pvs <- posthoc.kruskal.nemenyi.test(bet)$p.value
n_pvs <- apply(expand.grid(le, le), 1, function(x) paste(x, collapse='-'))
pvs <- c(cbind(rbind(rep(NA,ncol(pvs)),pvs),rep(NA,nrow(pvs)+1)))
names(pvs) <- n_pvs
pvs <- pvs[is.na(pvs) == F]

e <- factor(unlist(lapply(seq_along(bet), function(x) rep(le[x],length(bet[[x]])))), levels=le)
bet <- cbind.data.frame(div=unlist(bet), env=e)

mcl_beta <- multcompLetters2(div~env, pvs, bet, Letters=1:4)$Letters[le]
if(max(nchar(mcl_beta)) > 1){
  mcl_beta[nchar(mcl_beta) == 2] <- sapply(strsplit(mcl_beta[nchar(mcl_beta) == 2], ''), function(x) paste(x, collapse='|'))
}

#---
boxplot(div~env, bet, ylab='Bray-Curtis dist', ylim=range(bet$div)*c(1, 1.2), col=pal_env, xaxt='n')
usr <- par('usr')

axis(1, at=seq_along(le), labels=F)
text(1:length(le), usr[3]-diff(usr[3:4])*0.15, labels=le, xpd=NA, cex=0.75, srt=45)

text(usr[1]-diff(usr[1:2])*offset, usr[4]+diff(usr[3:4])*offset, 'C', xpd=NA, font=2)

text(1:length(le), usr[4]-diff(usr[3:4])*0.1, mcl_beta, offset=0, pos=3)

# NMDS ----
set.seed(0)
bc_tot <- vegdist(mr_sort)
nmds <- metaMDS(bc_tot, trace=F)

#---
plot(nmds$points, col=pal_env[as.numeric(env_sort$env)], pch=19)
usr <- par('usr')

text(usr[1]-diff(usr[1:2])*offset, usr[4]+diff(usr[3:4])*offset, 'D', xpd=NA, font=2)

text(usr[1]+diff(usr[1:2])*0.05, usr[4]-diff(usr[3:4])*0.05, paste('stress:', signif(nmds$stress, 2)), 
     cex=0.75, pos=4)

dev.off()

# test
pv_nmds <- NULL
for(i in le){
  ind_rmv <- labels(bc_tot) != i
  bc <- as.dist(as.matrix(bc_tot)[-ind_rmv,-ind_rmv])
  
  nmds <- metaMDS(bc, trace=F)
  pv_nmds[[paste(le[le != i], collapse='-')]] <- envfit(nmds, data.frame(env_sort$env[-ind_rmv]), permutations=10000)$factors$pvals
}

# NMDS Dinophyceae
mr_dino <- mr_sort[,grep('Dinophyceae', ass_sort$taxo)]
mr_dino <- mr_dino[rowSums(mr_dino) != 0,]

set.seed(0)
nmds_dino <- metaMDS(vegdist(mr_dino), trace=F) 

# remove weird smp
out <- unlist(sapply(apply(nmds_dino$points, 2, function(x) return(boxplot(x, plot=F)$out)), names))

mr_dino2 <- mr_dino[row.names(mr_dino) %in% out == F,]
mr_dino2 <- mr_dino2[,colSums(mr_dino2) != 0]

nmds_dino <- metaMDS(vegdist(mr_dino2), trace=F) 

#---
cairo_ps(paste0(dir_bout, 'nmds_dino.eps'), width=lar1.5, height=lar1.5)

plot(nmds_dino$points, col=pal_env[as.numeric(env_sort$env[row.names(env_sort) %in% row.names(mr_dino)])], pch=19)
usr <- par('usr')

# text(nmds_dino$points, labels=row.names(nmds_dino$points))

text(usr[1]+diff(usr[1:2])*0.95, usr[4]-diff(usr[3:4])*0.05, paste('stress:', signif(nmds_dino$stress, 2)), 
     cex=0.75, pos=2)

dev.off()


# Abundance and alpha diversity vs ecosystem for each taxa ----
print('--- alpha per taxa ---')

ind_eps <- ind_tax <- 1
cairo_ps(paste0(dir_bout, 'diversity_per_taxon_', ind_eps, '.eps'), width=13, height=7)
par(mfrow=c(3,6), mar=c(4,4,4,1), mgp=c(2.5,1,0))

lst_mcl <- NULL

for(i in n_lst){
  for(j in c('abds','rich','div')){
    n <- switch(j, 'abds'='abundance', 'rich'='richness', 'div'='diversity')
    met <- lapply(lst_div[[i]], function(x) x$div[[j]])

    if(j != 'div'){
      met <- lapply(met, function(x) decostand(x, 'log'))
    }
    
    e <- factor(unlist(lapply(seq_along(met), function(x) rep(le[x],length(met[[x]])))), levels=le)
    
    met <- cbind.data.frame(div=unlist(met), env=e)
    
    # test
    pvs <- posthoc.kruskal.nemenyi.test(met$div, met$env)$p.value
    n_pvs <- apply(expand.grid(colnames(pvs), 
                               row.names(pvs)), 1, function(x) paste(x, collapse='-'))[-2]
    pvs <-  na.omit(c(pvs))
    names(pvs) <- n_pvs
    
    mcl <- multcompLetters2(div~env, pvs, met, Letters=1:4)$Letters[le]
    if(max(nchar(mcl)) > 1){
      mcl[nchar(mcl) == 2] <- sapply(strsplit(mcl[nchar(mcl) == 2], ''), function(x) paste(x, collapse='|'))
    }
    
    lst_mcl[[i]][[j]] <- mcl
    
    # graf
    boxplot(div~env, met, ylab=ifelse(j == 'div', n, paste('log', n, '+ 1')), ylim=range(met$div)*c(1, 1.2), col=pal_env, names=F)
    usr <- par('usr')
    
    text(1:3, usr[3]-diff(usr[3:4])*0.23, le, srt=45, xpd=NA)
    
    mtext(paste(i,n), 3, 1.5, cex=0.7, font=2)
    
    text(1:3, usr[4]-diff(usr[3:4])*0.1, mcl, offset=0, pos=3)
    
  }
  
  ind_tax <- ind_tax +1
  if(ind_tax %% 6 == 1){
    dev.off()
    
    ind_eps <- ind_eps + 1
    cairo_ps(paste0(dir_bout, 'diversity_per_taxon_', ind_eps, '.eps'), width=13, height=7)
    par(mfrow=c(3,6), mar=c(4,4,4,1), mgp=c(2.5,1,0))
  }
}

dev.off()

mcls_alpha <- as.data.frame(sapply(lst_mcl, function(x) x$div))

# Gamma diversity vs ecosystem for each taxa ----
print('--- gamma per taxa ---')

mar <- c(2,2,4,1.5)
mgp <- c(2.5,1,0)

ind_eps <- ind_tax <- 1

lay <- function(){
  layout(rbind(cbind(rep(1,3), matrix(1:18+2, nrow=3, byrow=T)), c(0,rep(2,6))), height=c(1,1,1,0.3), width=c(0.3,rep(1,6)))
  par(mar=rep(0,4))
  
  plot.new()
  text(0.5,0.5, 'number of OTUs', srt=90)
  
  plot.new()
  text(0.5,0.5, 'number of samples')
  
  par(mar=mar, mgp=mgp)
}

cairo_ps(paste0(dir_bout, 'specaccum_per_taxon_', ind_eps, '.eps'), width=12, height=7)
lay()

for(i in c('Total', main_taxa)){
  
  spec <- lst_div[[i]]
  
  plot(NA, xlim=c(0, max(sapply(spec, function(x) length(x$div$accuR)))), 
       ylim=range(sapply(spec, function(x) {return(c(min(apply(x$div[,2:1], 1, diff)), 
                                                     sum(x$pool$boot, x$pool$boot.se)))})),
       xlab='', ylab='', main=i)
  usr <- par('usr')
  
  for(j in le){
    
    col <- pal_env[[j]]
    
    # specaccum
    accu <- spec[[j]]$div[,c('accuSD','accuR')]
    
    if(sum(accu) == 0) next
    
    lines(accu$accuR, col=col)
    
    seg <- apply(as.matrix(accu), 1, function(x) return(c(diff(x), sum(x))))  
    segments(seq_along(accu$accuSD), seg[1,], seq_along(accu$accuSD), seg[2,], col=col)
    
    # specpool
    bse <- unlist(spec[[j]]$pool[,c('boot','boot.se')])
    abline(h=bse[1], col=col)
    abline(h=bse[1]+c(1,-1)*bse[2], col=col, lty=2)
    
  }
  
  ind_tax <- ind_tax + 1
  
  #---
  if(ind_tax %% 18 == 1){
    dev.off()
    
    ind_eps <- ind_eps + 1
    
    cairo_ps(paste0(dir_bout, 'specaccum_per_taxon_', ind_eps, '.eps'), width=12, height=7)
    lay()
  }
  
}

dev.off()


# VENN ####
print('##### Venn ####')

pa <- sapply(mr_sort, function(x) tapply(x, env_sort$env, function(y) ifelse(sum(y) > 0, 1, 0)))

pa <- as.data.frame(t(pa))

pa_code <- t(c(2^2,2^1,2^0)*t(pa))
venn <- table(rowSums(pa_code))
lvenn <- log(venn)

n_smp <- c(tab_ech[3],tab_ech[2],sum(tab_ech[2:3]),tab_ech[1],
           sum(tab_ech[c(1,3)]), sum(tab_ech[1:2]), sum(tab_ech))
# pal <- grey.colors(101, start=0.4, end=1)[102-round((log(venn)+abs(min(log(venn))))/
#                                     max(log(venn)+abs(min(log(venn))))*100+1)]
alpha <- seq(0.4,1, length.out=100)[round((lvenn-min(lvenn))/(max(lvenn)-min(lvenn)) * 99 +1)]

rgb <- col2rgb(pal_env)

cmbn <- combn(le, 2)
rgb_pair <- apply(cmbn, 2, function(x) rowMeans(rgb[,x]))
colnames(rgb_pair) <- apply(cmbn, 2, function(x) paste(x, collapse='_'))

rgb_ext <- rbind(cbind(rgb,rgb_pair,colMeans(rgb))[,c(3,2,6,1,5,4,7)]/255, alpha)

pal <- apply(rgb_ext, 2, function(x) rgb(x[1],x[2],x[3],x[4]))

#---
cairo_ps(paste(dir_bout, 'venn.eps', sep=''), width=lar1.5, height=lar1.5)
par(mar=rep(0, 4))

venn(3, snames=names(pa), cexsn=1)

for(i in seq_along(venn)){
  
  polygon(getZones(i,3)[[1]], col=pal[i])
  text(getCentroid(getZones(i,3))[[1]][1], getCentroid(getZones(i,3))[[1]][2], cex=0.75,
       labels=paste(venn[i], '\nn=', n_smp[i], sep=''))
  
}

dev.off()

# search taxo omni
omni_otus <- row.names(pa)[rowSums(pa) == 3]


# FUNCTIONAL GROUPS ####
print('##### Functional group ####')
# pie charts ----
print('--- pie chart ---')

fct <- levels(ass_sort$fct)
mr_pa <- decostand(mr_sort, 'pa')

# rearrange Chrysophyceae
ind_cp <- which(taxo_sort$class == 'Chrysophyceae' & ass_sort$fct == 'phototrophic')
ind_ch <- which(taxo_sort$class == 'Chrysophyceae' & ass_sort$fct == 'consumer')
ind_cu <- which(taxo_sort$class == 'Chrysophyceae' & ass_sort$fct == 'unknown')

# reallocation séquence
rs <- cbind(rowSums(mr_sort[,ind_cp]), rowSums(mr_sort[,ind_ch]), rowSums(mr_sort[,ind_cu]))

seq_realloc_c <- t(apply(rs, 1, function(x) {
  y <- x[1:2]
  if(x[1] == 0) {
    y[2] <- x[2] + x[3]
  } else if (x[2] == 0) {
    y[1] <- x[1] + x[3]
  } else {
    y[1] <- x[1] + x[3] * (x[1]/sum(x[1:2]))
    y[2] <- x[2] + x[3] * (x[2]/sum(x[1:2]))
  }
  return(y)
}))

colnames(seq_realloc_c) <- c('Chrysophyceae_p','Chrysophyceae_c')

# reallocation otu
rs <- cbind(rowSums(mr_pa[,ind_cp]), rowSums(mr_pa[,ind_ch]), rowSums(mr_pa[,ind_cu]))

otu_realloc_c <- t(apply(rs, 1, function(x) {
  y <- x[1:2]
  if(x[1] == 0) {
    y[2] <- x[2] + x[3]
  } else if (x[2] == 0) {
    y[1] <- x[1] + x[3]
  } else {
    y[1] <- x[1] + x[3] * (x[1]/sum(x[1:2]))
    y[2] <- x[2] + x[3] * (x[2]/sum(x[1:2]))
  }
  return(y)
}))

colnames(otu_realloc_c) <- c('Chrysophyceae_p','Chrysophyceae_c')


# Dinophyceae ---
ind_din <- which(taxo_sort$class == 'Dinophyceae' & ass_sort$fct == 'unknown')

# reallocation séquence
rs <- rowSums(mr_sort[,ind_din])

seq_realloc_d <- t(sapply(seq_along(rs), function(x){
  if(env_sort$env[x] == 'marine'){
    y <- rs[x]*c(0.42, 1-0.42)
  } else {
    y <- rs[x]*c(0.88, 1-0.88)
  }
  return(y)
}))

colnames(seq_realloc_d) <- c('Dinphyceae_p','Dinphyceae_c')

# reallocation OTU
rs <- rowSums(mr_pa[,ind_din])

otu_realloc_d <- t(sapply(seq_along(rs), function(x){
  if(env_sort$env[x] == 'marine'){
    y <- rs[x]*c(0.42, 1-0.42)
  } else {
    y <- rs[x]*c(0.88, 1-0.88)
  }
  return(y)
}))

colnames(otu_realloc_d) <- c('Dinphyceae_p','Dinphyceae_c')

# supress the Chryso|Dino and replace by the reallocation
old_ind <- c(ind_cp, ind_ch, ind_cu, ind_din)

# mr ---
mr_pie <- cbind.data.frame(rbind.data.frame(mr_sort[,-old_ind], mr_pa[,-old_ind]),
                           rbind.data.frame(seq_realloc_c, otu_realloc_c), rbind.data.frame(seq_realloc_d, otu_realloc_d))

# taxo
mat_new <- matrix(c(c(as.character(unlist(taxo_sort[taxo_sort$class == 'Chrysophyceae',][1,1:3])),
                      'Chrysophyceae_p','Chrysophyceae_p_X','Chrysophyceae_p_XX','Chrysophyceae_p_XXX','Chrysophyceae_p_XXX.sp'),
                    c(as.character(unlist(taxo_sort[taxo_sort$class == 'Chrysophyceae',][1,1:3])),
                      'Chrysophyceae_c','Chrysophyceae_c_X','Chrysophyceae_c_XX','Chrysophyceae_c_XXX','Chrysophyceae_c_XXX.sp'),
                    c(as.character(unlist(taxo_sort[taxo_sort$class == 'Dinophyceae',][1,1:3])),
                      'Dinophyceae_p','Dinophyceae_p_X','Dinophyceae_p_XX','Dinophyceae_p_XXX','Dinophyceae_p_XXX.sp'),
                    c(as.character(unlist(taxo_sort[taxo_sort$class == 'Dinophyceae',][1,1:3])),
                      'Dinophyceae_c','Dinophyceae_c_X','Dinophyceae_c_XX','Dinophyceae_c_XXX','Dinophyceae_c_XXX.sp')),
                  nrow=4, byrow=T, dimnames=list(c('Chrysophyceae_p','Chrysophyceae_c','Dinophyceae_p','Dinophyceae_c'),
                                                 c(colnames(taxo_sort))))

taxo_pie <- rbind.data.frame(as.matrix(taxo_sort[-old_ind,]), mat_new)

# funct
fct_otu_pie <- factor(c(as.character(ass_sort$fct[-old_ind]), 'phototrophic','consumer', 'phototrophic','consumer'), levels=fct)

# arrangement of the pie arguments
selec_smp <- factor(apply(sapply(expand.grid(env_sort$env, c('relative abundance','richness')), as.character), 
                          1, function(x) paste(x, collapse=' ')))
n_selec_smp <- rep(levels(selec_smp),length(fct)+1)
selec_smp <- lapply(rep(levels(selec_smp),length(fct)+1), function(x) which(selec_smp == x))
names(selec_smp) <- n_selec_smp

selec_otu <- sapply(levels(as.factor(fct_otu_pie)), function(x) names(mr_pie)[which(fct_otu_pie == x)])
selec_otu[['overall']] <- names(mr_pie)

selec_otu2 <- NULL
n_selec_otu2 <- gl(length(fct)+1,length(le)*2, labels=names(selec_otu))
for(i in seq_along(n_selec_otu2)){
  selec_otu2[[i]] <- selec_otu[[n_selec_otu2[i]]]
}
names(selec_otu2) <- n_selec_otu2

# pie
pt <- pie_taxo(mr_pie, taxo_pie, 1:4, selec_smp, selec_otu2, 0.01, show=F,
               pal_ini=colorRampPalette(brewer.pal(8,'Set1')[-c(2,5)]),
               new_order=list(phylum=c("Alveolata","Stramenopiles","Rhizaria","Eukaryota_X","Archaeplastida",
                                       "Opisthokonta","Apusozoa","Amoebozoa","Excavata","Protalveolata" )))

# per grp func
for(i in c(fct, 'overall')){
  ind_sort <- which(names(selec_otu2) == i)+4
  
  cairo_ps(paste0(dir_bout, 'pie_', i,'.eps'), width=lar2*1.5, height=lar2*1.5)
  layout(cbind(matrix(c(4,3,2,1,6,5),ncol=2, byrow=T), rep(7,3)), c(1,1,1.5))
  par(mar=c(1.5,2,2,1.5), oma=rep(2, 4), xaxs='i', yaxs='i', xpd=NA)
  
  for(j in ind_sort){
    plot.new()
    title(names(selec_smp)[j-4])
    pie_taxo_single(pie_taxo=pt, sel_smp=j, x=0.5, y=0.5, ray=0.45, cex=0.75, last_tax_text=F)
  }
  
  plot.new()
  legend_pie_taxo(pt, 0.5,0.5, cex=0.7, last_tax_text=F)
  
  dev.off()
}

# relative abundance ----
print('--- rel abu fct ---')
mr <- mr_pie[1:nrow(env_sort),]

mr_funct <- aggregate(t(mr), list(fct_otu_pie), sum)
row.names(mr_funct) <- mr_funct[,1]
mr_funct <- as.data.frame(decostand(t(mr_funct[,-1]), 'total'))

#---
mcls_fct <- NULL
for(i in mr_funct){
  pvs <- posthoc.kruskal.nemenyi.test(i, env_sort$env)$p.value
  n_pvs <- apply(expand.grid(le, le), 1, function(x) paste(x, collapse='-'))
  pvs <- c(cbind(rbind(rep(NA,ncol(pvs)),pvs),rep(NA,nrow(pvs)+1)))
  names(pvs) <- n_pvs
  pvs <- pvs[is.na(pvs) == F]
  
  fct_perc <- data.frame(perc=i, env=env_sort$env)

  mcl <- multcompLetters2(perc~env, pvs, fct_perc, Letters=1:4)$Letters[le]
  if(max(nchar(mcl)) > 1){
    mcl[nchar(mcl) == 2] <- sapply(strsplit(mcl[nchar(mcl) == 2], ''), function(x) paste(x, collapse='|'))
  }
  
  mcls_fct <- c(mcls_fct, mcl)
}

#---
# pdf(paste0(dir_bout,'nemen_fct_relabu.pdf'), width=lar1, height=lar1)
cairo_ps(paste0(dir_bout,'nemen_fct_relabu.eps'), width=lar1, height=lar1)
par(mar=c(5,4,1,1))

lf <- length(fct)
xs <- 1:lf

plot.new()
plot.window(c(0,lf),range(mr_funct)*c(1,1.2))
usr <- par('usr')

box('plot')
abline(v=1:3, lwd=0.7, lty=3)

# bxp
sht <- 0.2
seq <- seq(0+sht, 1-sht, length.out=length(le))
for(i in seq_along(le)){
  boxplot(mr_funct[env_sort$env == le[i],], boxwex=seq[1]*0.75, col=pal_env[i], at=seq[i]+c(0:(lf-1)),
          add=T, names=F, axes=F, lwd=0.5)
}

# axis
ys <- seq(0,1,length.out=6)
axis(2, ys, sprintf("%.1f", ys), las=2, cex.axis=0.7)
mtext('sequences relative abundance', 2, 2.5, cex=0.7)

text(xs-0.5, usr[3]-diff(usr[3:4])*0.20, fct, xpd=NA, cex=0.7, srt=45)
axis(1, c(1:lf)-0.5, labels=NA)

# test
text(apply(expand.grid(seq, 1:length(fct)-1), 1, sum), usr[4]-diff(usr[3:4])*0.15,
     mcls_fct, cex=0.7, pos=3)

dev.off()


# PID VS ECOSYSTEMS ####
print('##### PID vs ecosystems #####')

lst_pid <- list(marine    =ass_sort$pid[colSums(mr_sort[env_sort$env == le[1],]) != 0],
                freshwater=ass_sort$pid[colSums(mr_sort[env_sort$env == le[2],]) != 0],
                soil      =ass_sort$pid[colSums(mr_sort[env_sort$env == le[3],]) != 0])

pvs <- posthoc.kruskal.nemenyi.test(lst_pid)$p.value
n_pvs <- apply(expand.grid(le[1:2], le[2:3]), 1, function(x) paste(x, collapse='-'))[-2]
pvs <-  as.numeric(na.omit(c(pvs)))
names(pvs) <- n_pvs

fac_pid <- sapply(lst_pid, length)
pid <- data.frame(pid=unlist(lst_pid), env=factor(c(rep(le[1],fac_pid[1]),rep(le[2],fac_pid[2]),rep(le[3],fac_pid[3]))))

mcl_pid <- multcompLetters2(pid~env, pvs, pid, Letters=1:4)$Letters[le]
if(max(nchar(mcl_pid)) > 1){
  mcl_pid[nchar(mcl_pid) == 2] <- sapply(strsplit(mcl_pid[nchar(mcl_pid) == 2], ''), function(x) paste(x, collapse='|'))
}

#---
cairo_ps(paste0(dir_bout, 'pid.eps'), width=lar1.5, height=lar1.5)

boxplot(lst_pid, ylim=range(lst_pid)*c(1,1.1), yaxt='n', ylab='percentage of identity', col=pal_env)
usr <- par('usr')

a <- axis(2, labels=F, tick=F)
axis(2, yaxp=c(a[1], 100, length(a[a <= 100])-1))

text(1:3, usr[4] - diff(usr[3:4])*0.1, mcl_pid)

dev.off()

#---
pids <- lapply(lst_pid, function(x) return(c(mean=mean(x), se=sd(x)/sqrt(length(x)))))


# PIE CHART PAR TAXA ####
print('##### pie chart taxo ####')

# make the schematic phylo tree
file <- paste(dir_save, 'tax_pie.csv', sep='')
if(file.exists(file)){
  tax <- read.table(file, sep=',', row.names=1, h=T)
} else {
  tax <- as.data.frame(matrix(NA, ncol=length(dom_tax), nrow=length(dom_tax)))
  dimnames(tax) <- list(dom_tax,dom_tax)
  for(i in 1:nrow(tax)){
    for(j in 1:ncol(tax)){
      if(i == j){
        tax[i,j] <- 0}
      if(i < j){
        tax[i,j] <- ''
      }
    }
  }

  write.table(tax, file=file, sep='\t')
  
  stop()
  
}

# retreive the bootstrap esti of richness
pred <- sapply(lst_div, function(x) sapply(x, function(y) {
  z <- c(y$pool$boot, y$pool$boot.se)
  if(is.null(z)){
    z <- c(NA, NA)
  }
  return(z)
}), simplify='array') # Diversity calculation

# test of the bootstrap
test_gamma <- apply(pred, 3, function(x) {
  up <- x[1,] + x[2,]
  down <- x[1,] - x[2,]
  
  rnk <- rank(x[1,], na.last='keep')
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

row.names(test_gamma) <- le

# dendro ---
tax2 <- tax
row.names(tax2) <- names(tax2) <- sub('.', '', names(tax), fixed=T)
hc <- hclust(as.dist(tax2))
dend <- as.dendrogram(hc)

n1 <- dend[[2]][[2]][[2]][[1]][[1]]
n2 <- dend[[2]][[2]][[2]][[1]][[2]]
dend[[2]][[2]][[2]][[1]][[1]] <- n2
dend[[2]][[2]][[2]][[1]][[2]] <- n1

# palette ---
pal_tax <- unlist(pt$lst_pal) # pie charts
n <- names(pal_tax)
n <- sapply(strsplit(n, '.', fixed=T), '[[', 2)
n <- sub('Filosa-','', n)
names(pal_tax) <- n
pal_tax <- pal_tax[-grep('other', names(pal_tax))]

#---
t1 <- c('Excavata','Amoebozoa','Opisthokonta','Archaeplastida','Cryptista','Haptista',
        'Rhizaria','Stramenopiles','Alveolata')
p1 <- NULL
for(i in 1:length(t1)){if(t1[i] %in% n){p1 <- c(p1, pal_tax[t1[i]])} else {p1 <- c(p1, NA)}}
# p1[c(4,6)] <- 'white'
names(p1) <- t1
rgb <- rowMeans(col2rgb(pal_tax[grep('Haptophyta|Centroheliozoa', names(pal_tax))]))/255
p1['Haptista'] <- rgb(rgb[1],rgb[2],rgb[3])
rgb <- rowMeans(col2rgb(pal_tax[grep('Katablepharidaceae|Cryptophyceae', names(pal_tax))]))/255
p1['Cryptista'] <- rgb(rgb[1],rgb[2],rgb[3])
t1[c(2,6)] <- c('Amoebo-\nzoa','Hap-\ntista')

p1p <- apply(col2rgb(p1), 2, function(x) {
  y <- x/255
  rgb(y[1],y[2],y[3],0.1)
})

#---
t2 <- c('Holozoa','Streptophyta','Chlorophyta','Cercozoa','Endomyxa','Retaria',
        'Ochrophyta','Ciliophora','Dinoflagellata')
p2 <- NULL
for(i in 1:length(t2)){if(t2[i] %in% n){p2 <- c(p2, pal_tax[t2[i]])} else {p2 <- c(p2, NA)}}
names(p2) <- t2
p2['Retaria'] <- pal_tax['Radiolaria']
rgb <- rowMeans(col2rgb(pal_tax[grep('Endomyxa', names(pal_tax))]))/255
p2['Endomyxa'] <- rgb(rgb[1],rgb[2],rgb[3])
rgb <- rowMeans(col2rgb(pal_tax[grep('Choanoflagellida|Ichthyosporea', names(pal_tax))]))/255
p2['Holozoa'] <- rgb(rgb[1],rgb[2],rgb[3])
t2[c(2,5,9)] <- c('Strepto-\nphyta','Endo-\nmyxa','Dinofla-\ngellata')

p2p <- apply(col2rgb(p2), 2, function(x) {
  y <- x/255
  rgb(y[1],y[2],y[3],0.1)
})

t3 <- get_leaves_attr(dend, 'label')

# plot ---
xb1 <- 25
xsp <- 5
cex <- 0.6
x_bp <- -35
shift <- 17
max_pred <- 500
prop_bar <- 0.8
rad <- 1.5

#---
# pdf(paste(dir_bout,'/pie_tax.pdf', sep=''), width=lar1.5, height=10)
cairo_ps(paste(dir_bout,'/pie_tax.eps', sep=''), width=lar1.5, height=10, fallback_resolution=900)
par(mar=c(0,6,0,27), cex=0.75, xpd=NA)

plot(dend, axes=F, horiz=T, leaflab='none')

rect(x_bp,0.5, xb1,nrow(tax)+0.5, col=rgb(1,1,1,0.1), border=NA)

# phylum
y1 <- c(0,4,8, 14,22,25, 28,38,51)+0.5
y2 <- c(4,7,13, 22,25,27, 38,51,59)+0.5
u1 <- union(y1,y2)
segments(x_bp-3*shift+(shift*(1-prop_bar)), u1, xb1, u1)
# rect(x_bp-3*shift+(shift*(1-prop_bar)),  y1, xb1, y2, col=p1p, border=NA)
rect(xb1-xsp, y1, xb1, y2, col=p1, border=NA)
text(xb1-xsp/2, rowMeans(cbind(y1, y2)), labels=t1, srt=90, cex=cex)

# division
y1 <- c(10,15,18, 28,33,35, 38,51,57)+0.5
y2 <- c(13,18,22, 33,35,38, 46,56,59)+0.5
xb2 <- xb1-xsp*1.5
u2 <- union(y1,y2)
u2 <- u2[u2 %in% u1 == F]
segments(x_bp-3*shift+(shift*(1-prop_bar)), u2, xb2, u2, lty=2)
# rect(x_bp-3*shift+(shift*(1-prop_bar)), y1, xb2, y2, col=p2p, border=NA)
rect(xb2-xsp,  y1, xb2, y2, col=p2, border=NA)
text(xb2-xsp/2, rowMeans(cbind(y1, y2)), labels=t2, srt=90, cex=cex)

# classes
r_wdt <- min(strwidth(dom_tax, 'user', cex=cex))*1.3
ys <- seq(0.5, length(hc$labels)+1)
# rect(0, ys[-length(ys)], r_wdt, ys[-1], col=p3, border=F)
text(0, ys[-length(ys)]+0.45, t3, pos=4, cex=cex)

# barplot
x <- x_bp
x_test <- NULL
for(i in 1:3){
  x_seq <- seq(0,max_pred, length.out=3)
  
  # lines(c(x,x-max(x_seq)*(shift/max_pred*prop_bar)), c(0,0))
  
  xs <- x-x_seq*(shift/max_pred*prop_bar)
  x_test <- c(x_test, xs[1])
  
  # segments(xs, rep(-0.1, length(x_seq)), xs, rep(0, length(x_seq)))
  segments(xs, rep(0.5, length(x_seq)), xs, rep(59.5, length(x_seq)), col=c('grey50','grey50',1), lwd=0.5)
  
  text(xs, rep(0, length(x_seq)), labels=x_seq, cex=cex)
  
  # text(x-0.5*shift, ncol(tax)+1, labels=le[i], cex=cex)
  text(xs[2], length(labels(dend))+1, labels=le[i], cex=cex)
  
  x <- x-shift
}

# ---
ind_tax <- 1
agg_tax <- NULL
for(i in labels(dend)){
  if(i %in% dom_tax){
    
    # pie
    mr_tax <- mr_sort[,grep(i, ass_sort$taxo)]
    agg <- aggregate(rowSums(mr_tax), list(env_sort$env), sum)[,2]
    
    floating.pie(r_wdt+(x_bp-r_wdt)/2, ind_tax, agg/tab_ech, radius = rad,
                 col=pal_env[as.logical(agg)], border=F, startpos=-pi/2)
    # floating.pie(-1, ind_tax, agg/tab_ech, radius = 0.25,
    #              col=pal_env[as.logical(agg)], border=F, startpos=-pi/2)
    
    agg_tax <- rbind(agg_tax, agg)
    
    # barplot
    lrg <- 0.2
    x <- x_bp
    for(j in 1:3){
      if(is.na(pred[1,j,i]) == F){
        correc <- (shift/max_pred*0.8)
        xbar <- x-pred[1,j,i]*correc
        se <- pred[2,j,i]*correc
        if(pred[1,j,i] < max_pred){
          rect(x, ind_tax-lrg, xbar, ind_tax+lrg, col=pal_env[j], border=NA)
          segments(xbar-se, ind_tax, xbar+se, ind_tax)
        } else {
          text(x, ind_tax-lrg, labels=paste(round(pred[1,j,i]), '±', round(pred[2,j,i])), adj=c(-0.05,0), cex=cex)
        }
      }
      x <- x - shift
    }
    
    # test bootstrap
    text(x_bp - shift*c(0,1,2) + shift*0.07, ind_tax, test_gamma[,i], cex=cex*0.66)
    
  }
  
  ind_tax <- ind_tax + 1
}

text(x_bp-1.5*shift, -1, 'predicted number of OTUs', cex=cex)

dev.off()




#####

file <- paste0(dir_bave, 'out_', boot, '.Rdata')
save(lst_div, mcls_alpha, mcl_beta, test_gamma, pv_nmds, venn, omni_otus, mcls_fct, pids, pt, pal_env, file=file)






