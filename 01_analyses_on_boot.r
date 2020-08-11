#####
# ecomp prep
#####

rm(list=ls())

require(PMCMR)
# require(venn)
require(multcompView)
require(multcomp)
require(plotrix)
require(RColorBrewer)
# require(betapart)
# require(dendextend)
require(foreach)
require(doSNOW)

source('bin/src/my_prog/R/pie_taxo.r')
source('bin/src/my_prog/R/pie_taxo_single.r')
source('bin/src/my_prog/R/legend_pie_taxo.r')

argv <- commandArgs(trailingOnly=T)

seed <- as.integer(argv[1])
boot <- argv[2]
raref <- as.integer(argv[3])
th <- as.integer(argv[4])

##########
seed <- 1
boot <- '001'
raref <- 10000
th <- 4
##########

#---
# cluster
cl <- makeSOCKcluster(th)
registerDoSNOW(cl)

#---
# param variable

lar1 <- 3.54
lar1.5 <- 5.51
lar2 <- 7.48

dom_tax <- c('Apusozoa',
             'Amoebozoa*',
             'Mycetozoa',
             'Variosea','Lobosa',
             'Opisthokonta*',
             'Fungi*',
             'Ichthyosporea','Nucleariidea',
             'Choanoflagellida',
             'Metazoa*',
             'Excavata*',
             'Heterolobosea',
             'Kinetoplastida',
             'Diplonemea','Euglenida',
             'Goniomonadales','Cryptomonadales',
             'Katablepharidaceae',
             'Picozoa',
             'Archaeplastida*',
             'Rhodophyta*',
             'Prasino',
             'Chlorophyceae','Mamiellophyceae','Trebouxiophyceae',
             'Streptophyta*',
             'Embryophyceae*',
             'Klebsormidiophyceae','Zygnemophyceae',
             'Centroheliozoa',
             'Haptophyta',
             'Telonemia',
             'Rhizaria*',
             'Foraminifera*',
             'Globothalamea','Monothalamids',
             'Radiolaria',
             'Cercozoa*',
             'Vampyrellida', 'Plasmodiophorida',
             'Granofilosea','Imbricatea','Sarcomonadea','Thecofilosea','Phaeodarea',
             'Stramenopiles*',
             'MAST',
             'Labyrinthulomycetes',
             'Oomycota','Hyphochytriomyceta',
             'Bicoecea',
             'Xanthophyceae','Eustigmatophyceae',
             'Chrysophyceae','Synurophyceae',
             'Bacillariophyta','Bolidophyceae',
             'Dictyochophyceae','Pelagophyceae',
             'Alveolata*',
             'Ciliophora*',
             'Colpodea','Spirotrichea','Litostomatea','Oligohymenophorea','Phyllopharyngea',
             'Apicomplexa',
             'Dinophyta*',
             'Dinophyceae','Syndiniales'
) # est-ce qu'on vérifie les seuil pid?

main_taxa <- grep('*', dom_tax, invert=T, fixed=T, value=T)


# DOWNLOAD ####
dir_base <- 'Projets/Ecomp/stat/'
# dir_in   <- paste0(dir_base, 'in/')
dir_save <- paste0(dir_base, 'out/200807/saves/')

dir_boot <- paste0(dir_base, 'bootstrap/200807/bootstrap_', boot, '/')
dir_bave <- paste0(dir_boot, 'saves/')
dir_bout <- paste0(dir_boot, 'out/')

lapply(list(dir_bave,dir_bout), function(x) dir.create(x, F, T))

# community
set.seed(seed)

load(paste0(dir_save, 'prep.Rdata'))

file <- paste0(dir_boot, 'boot_', boot, '.Rdata')
if(file.exists(file)){
  load(file)
} else {
  # rrarefy
  print('rrarefy')

  mr <- mr_sl[rowSums(mr_sl) > raref,]
  e <- env[row.names(mr),]
  mr <- mr[unlist(tapply(row.names(e), list(e$env), function(x) sample(x, 20))),]
  
  # rrarefy
  print('rrarefy')
  mr_boot <- as.data.frame(rrarefy(mr, raref))
  mr_boot <- droplevels(mr_boot[,colSums(mr_boot) != 0])
  
  env_boot <- e[row.names(mr_boot),]
  ass_boot  <- droplevels( ass_sl[names(mr_boot),])
  taxo_boot <- droplevels(taxo_sl[names(mr_boot),])
  
  # save
  print('save')
  save(env_boot, mr_boot, ass_boot, taxo_boot, file=file)
}

# sort the bootstrap according to taxonomy and environments
ass_boot$taxo <- gsub('Hacrobia','Eukaryota_X',ass_boot$taxo)

env_boot$env <- factor(env_boot$env, levels=levels(env_boot$env)[c(2,1,3)])
le <- levels(env_boot$env)

ord_otu <- order(ass_boot$taxo)
ord_smp <- order(env_boot$env)

env_sort <- env_boot[ord_smp,]
mr_sort <- mr_boot[ord_smp,ord_otu]
ass_sort <- ass_boot[ord_otu,]
taxo_sort <- taxo_boot[ord_otu,]

# palette
pal_env <- rgb(c(0,0,1),c(0,196,102)/255, c(1,1,0))
names(pal_env) <- le

# CHECK SOME STUFF ####
# # marine fraction ----
# pal_mar <- pal_sub[1:5]
# 
# tb <- table(droplevels(env_boot$ech[env_boot$env == 'marine']))
# env_mar <- droplevels(env_boot[env_boot$ech %in% names(tb)[tb == 5],])
# 
# mr_mar <- mr_sort[row.names(env_mar),]
# mr_mar <- mr_mar[,colSums(mr_mar) != 0]
# 
# # fake full frac
# ind_frac <- env_mar$sub != 'marine'
# 
# mr_fff <- as.data.frame(rrarefy(mr_mar[ind_frac,], 2500))
# mr_fff <- sapply(mr_fff, function(x) tapply(x, list(env_mar$ech[ind_frac]), sum))
# row.names(mr_fff) <- paste0(row.names(mr_fff), 'f')
# mr_fff <- rbind.data.frame(mr_fff,mr_mar[ind_frac == F,])
# mr_fff <- mr_fff[,colSums(mr_fff) != 0]
# 
# env_fff <- rbind(env_mar[ind_frac == F,],env_mar[ind_frac == F,])
# row.names(env_fff) <- row.names(mr_fff)
# env_fff$sub <- as.character(env_fff$sub)
# env_fff$sub[grep('f', row.names(env_fff))] <- 'marine_F'
# env_fff$sub <- as.factor(env_fff$sub)
# 
# #---
# lst_mr_mar <- list(frac=list(mr=mr_mar,
#                              env=env_mar),
#                    fff =list(mr=mr_fff,
#                              env=env_fff))
# 
# #---
# for(i in names(lst_mr_mar)){
#   
#   mr <- lst_mr_mar[[i]]$mr
#   env <- lst_mr_mar[[i]]$env
#   
#   nmds_mar <- metaMDS(vegdist(mr))
#   
#   #---
#   cairo_ps(paste0(dir_bout, 'NMDS_marine_', i, '.eps'))
#   
#   plot(nmds_mar$points, type='n')
#   ordispider(nmds_mar, env$ech, col='grey')
#   points(nmds_mar$points, pch=19, col=pal_mar[as.numeric(env$sub)])
#   
#   legend('topright', legend=levels(env$sub), col=pal_mar[1:length(levels(env$sub))], pch=19, bty='n')
#   
#   dev.off()
#   
#   if(i == 'fff'){
#     coord1 <- metaMDS(vegdist(mr[env$sub == 'marine',]))$points
#     coord2 <- metaMDS(vegdist(mr[env$sub == 'marine_F',]))$points
#     row.names(coord1) <- row.names(coord2) <- levels(env$ech)
#     
#     pv <- signif(protest(coord1,coord2)$signif, 2)
#     
#     cairo_ps(paste0(dir_bout, 'procruste_marine.eps'))
#     plot(procrustes(coord1,coord2), main=paste('from full to frac\npv =', pv))
#     dev.off()
#   }
#   
# }
# 
# #
# 
# DIVERSITY AND ABUNDANCE ####
# Diversity calculation ----
print('##### Diversity ####')

# calculation of metrics

# cluster
cl2 <- makeSOCKcluster(th/2)
registerDoSNOW(cl2)
clusterEvalQ(cl2, library(doSNOW))
clusterEvalQ(cl2, library(foreach))
clusterEvalQ(cl2, library(vegan))

# calculation of metrics
n_lst <- c('Total', main_taxa)

file <- paste0(dir_bave, 'lst_div_', boot, '.Rdata')
if(file.exists(file)){
  load(file)
} else {
  lst_div <- foreach(i = n_lst, .verbose=T) %dopar% {
    
    print(paste('###',i,'###'))
    
    ind_otu <- 1:ncol(mr_sort)
    if(i != 'Total'){
      ind_otu <- grep(i, ass_sort$taxo)
    }
    
    # cluster
    cl22 <- makeSOCKcluster(th/4)
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
        abds <- div <- rep(0, nrow(mr_sort[ind_smp,]))
        
        bc <- NULL
        
        accum <- data.frame(richness=abds,sd=abds)
        
        pool <- NULL
        
      }
      
      gc()
      
      # fill the list
      l <-list(div=data.frame(accuR=accum$richness, accuSD=accum$sd, abds=abds, div=div),
               pool=pool,
               bc=bc)
      
      return(l)
      
    }
    
    names(lst) <- le
    
    return(lst)
    
  }
  
  names(lst_div) <- n_lst

  save(lst_div, file=file)
}

# diversity Total figure
tot <- lst_div$Total

#---
cairo_ps(paste0(dir_bout, 'diversity.eps'), width=lar2, height=lar2)
par(mfrow=c(2,2), mar=c(5,4,2,1), mgp=c(2.5,1,0))

# Gamma ----

plot(NA, xlim=c(0,max(sapply(tot, function(x) length(x$div$accuR)))), 
     ylim=range(sapply(tot, function(x) {return(c(min(apply(x$div[,2:1], 1, diff)), 
                                                  sum(x$pool$boot, x$pool$boot.se)))})),
     xlab='number of samples', ylab='number of OTU')
usr <- par('usr')

text(usr[1]-diff(usr[1:2])*0.05, usr[4]+diff(usr[3:4])*0.05, 'A', xpd=NA, font=2)

for(i in le){
  
  col <- pal_env[[i]]
  
  # specaccum
  accu <- tot[[i]]$div[,c('accuSD','accuR')]
  
  lines(accu$accuR, col=col)
  
  seg <- apply(as.matrix(accu), 1, function(x) return(c(diff(x), sum(x))))  
  segments(seq_along(accu$accuSD), seg[1,], seq_along(accu$accuSD), seg[2,], col=col)
  
  # specpool
  boot <- unlist(tot[[i]]$pool[,c('boot','boot.se')])
  abline(h=boot[1], col=col)
  abline(h=boot[1]+c(1,-1)*boot[2], col=col, lty=2)
  
  # segment finale
  m <- max(accu$accuR)
  d <- floor(diff(c(m, boot[1])))
  
  x <- nrow(accu)
  segments(x, m, x, m+d, col=col, lty=3)
  
  text(x, m+0.5*d, d, pos=2)  
  
}

# Alpha ----
alp <- lapply(tot, function(x) x$div$div)

pvs <- posthoc.kruskal.nemenyi.test(alp)$p.value
n_pvs <- apply(expand.grid(le, le), 1, function(x) paste(x, collapse='-'))
pvs <- c(cbind(rbind(rep(NA,ncol(pvs)),pvs),rep(NA,nrow(pvs)+1)))
names(pvs) <- n_pvs
pvs <- pvs[is.na(pvs) == F]

mcl <- multcompLetters(pvs)$Letters
mcl <- mcl[c(length(mcl),seq_along(mcl)[-length(mcl)])]

#---
boxplot(alp, ylab='H', ylim=range(alp)*c(1, 1.2), col=pal_env, xaxt='n')
usr <- par('usr')

axis(1, at=seq_along(le), labels=F)
text(1:length(le), usr[3]-diff(usr[3:4])*0.15, labels=le, xpd=NA, cex=0.75, srt=45)

text(usr[1]-diff(usr[1:2])*0.05, usr[4]+diff(usr[3:4])*0.05, 'B', xpd=NA, font=2)

text(1:length(le), usr[4]-diff(usr[3:4])*0.1, mcl, offset=0, pos=3)

# Beta ----
bet <- as.data.frame(sapply(tot, function(x) c(x$bc)))

pvs <- posthoc.kruskal.nemenyi.test(bet)$p.value
n_pvs <- apply(expand.grid(le, le), 1, function(x) paste(x, collapse='-'))
pvs <- c(cbind(rbind(rep(NA,ncol(pvs)),pvs),rep(NA,nrow(pvs)+1)))
names(pvs) <- n_pvs
pvs <- pvs[is.na(pvs) == F]

mcl <- multcompLetters(pvs)$Letters
mcl <- mcl[c(length(mcl),seq_along(mcl)[-length(mcl)])]

#---
ulbet <- unlist(bet)
boxplot(bet, ylab='Bray-Curtis dist', ylim=c(min(ulbet), max(ulbet+diff(range(ulbet))*0.25)), col=pal_env, xaxt='n')
usr <- par('usr')

axis(1, at=seq_along(le), labels=F)
text(1:length(le), usr[3]-diff(usr[3:4])*0.15, labels=le, xpd=NA, cex=0.75, srt=45)

text(usr[1]-diff(usr[1:2])*0.05, usr[4]+diff(usr[3:4])*0.05, 'C', xpd=NA, font=2)

text(1:length(le), usr[4]-diff(usr[3:4])*0.1, mcl, offset=0, pos=3)

# NMDS ----
set.seed(0)
bc_tot <- vegdist(mr_sort)
nmds <- metaMDS(bc_tot)

ef <- envfit(nmds, data.frame(env_sort$sub), permutations=10000)

#---
plot(nmds$points, col=pal_env[as.numeric(env_sort$env)], pch=19)
usr <- par('usr')

text(usr[1]-diff(usr[1:2])*0.05, usr[4]+diff(usr[3:4])*0.05, 'D', xpd=NA, font=2)

text(usr[1]+diff(usr[1:2])*0.05, usr[4]-diff(usr[3:4])*0.05, paste('stress:', signif(nmds$stress, 2)), 
     cex=0.75, pos=4)

dev.off()

# NMDS Dinophyceae
cairo_ps(paste0(dir_bout, 'nmds_dino.eps'), width=lar1.5, height=lar1.5)

mr_dino <- mr_sort[,grep('Dinophyceae', ass_sort$taxo)]
mr_dino <- mr_dino[-which(row.names(mr_sort) %in% c('3034','3084','3078')),]# the two samples make complete outlyar
mr_dino <- mr_dino[row.names(mr_dino) != 0,]

set.seed(0)
nmds_dino <- metaMDS(vegdist(mr_dino)) 

plot(nmds_dino$points, col=pal_env[as.numeric(env_sort$env[row.names(env_sort) %in% row.names(mr_dino)])], pch=19)
usr <- par('usr')

# text(nmds_dino$points, labels=row.names(nmds_dino$points))

text(usr[1]+diff(usr[1:2])*0.95, usr[4]-diff(usr[3:4])*0.05, paste('stress:', signif(nmds_dino$stress, 2)), 
     cex=0.75, pos=2)

dev.off()


# Abundance and alpha diversity vs ecosystem for each taxa ----
# Gamma diversity vs ecosystem for each taxa ----
# VENN ####
# FUNCTIONAL GROUPS ####
# pie charts ----
fct <- levels(ass_sort$fct)
mr_pa <- decostand(mr_sort, 'pa')

# rearrange Chrysophyceae
ind_cp <- which(taxo_sort$class == 'Chrysophyceae' & ass_sort$fct == 'phototrophic')
ind_ch <- which(taxo_sort$class == 'Chrysophyceae' & ass_sort$fct == 'heterotrophic')
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
fct_otu_pie <- c(as.character(ass_sort$fct[-old_ind]), 'phototrophic','heterotrophic', 'phototrophic','heterotrophic')

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
pt <- pie_taxo(mr_pie, taxo_pie, 1:4, selec_smp, selec_otu2, 0.01, show=F, root='Eukaryota')

# per grp func
for(i in c(fct, 'overall')){
  ind_nf <- which(names(selec_otu2) == i)+4
  
  cairo_ps(paste0(dir_bout, 'pie_', i,'.eps'), width=lar2*1.5, height=lar2*1.5)
  layout(cbind(matrix(c(4,3,2,1,6,5),ncol=2, byrow=T), rep(7,3)), c(1,1,1.5))
  par(mar=c(1.5,2,2,1.5), oma=rep(2, 4), xaxs='i', yaxs='i', xpd=NA)
  
  for(j in ind_nf){
    print(j)
    plot.new()
    title(names(selec_smp)[j-4])
    pie_taxo_single(pie_taxo=pt, sel_smp=j, x=0.5, y=0.5, ray=0.45, cex=0.75, last_tax_text=F)
  }
  
  plot.new()
  legend_pie_taxo(pt, 0.5,0.5, cex=0.7, last_tax_text=F)
  
  dev.off()
}

# relative abundance ----

mr <- mr_pie[1:nrow(env_sort),]

n <- 'sequence nb relative abundance'

mr_funct <- aggregate(t(mr), list(fct_otu_pie), sum)
row.names(mr_funct) <- mr_funct[,1]
mr_funct <- as.data.frame(decostand(t(mr_funct[,-1]), 'total'))

#---
mcls <- NULL
for(i in mr_funct){
  pvs <- posthoc.kruskal.nemenyi.test(i, env_sort$env)$p.value
  n_pvs <- apply(expand.grid(le, le), 1, function(x) paste(x, collapse='-'))
  pvs <- c(cbind(rbind(rep(NA,ncol(pvs)),pvs),rep(NA,nrow(pvs)+1)))
  names(pvs) <- n_pvs
  pvs <- pvs[is.na(pvs) == F]
  
  mcl <- multcompLetters(pvs)$Letters
  mcl <- mcl[c(length(mcl),seq_along(mcl)[-length(mcl)])]
  
  mcls <- c(mcls, mcl)
}

#---
cairo_ps(paste0(dir_bout,'nemen_fct_relabu.eps'), width=lar2, height=lar1.5)
par(mar=c(7,4,2,1))

lf <- length(fct)
xs <- 1:lf

plot.new()
plot.window(c(0,lf),range(mr_funct)*c(1,1.2))
usr <- par('usr')
box('plot')

sht <- 0.2
seq <- seq(0+sht, 1-sht, length.out=length(le))
for(i in seq_along(le)){
  boxplot(mr_funct[env_sort$env == le[i],], boxwex=seq[1]*0.75, col=pal_env[i], at=seq[i]+c(0:(lf-1)),
          add=T, names=F, axes=F, lwd=0.5)
}

ys <- seq(0,1,length.out=6)
axis(2, ys, sprintf("%.1f", ys), las=2)

text(xs-0.5, usr[3]-diff(usr[3:4])*0.12, fct, xpd=NA, srt=45)
axis(1, c(1:lf)-0.5, labels=NA)

text(apply(expand.grid(seq, 1:length(fct)-1), 1, sum), usr[4]-diff(usr[3:4])*0.1,
     mcls, pos=3, offset=0, cex=0.6, srt=90)

dev.off()

# PID VS ECOSYSTEMS ####
# pid vs ecosystem for each taxon ----
# pid vs habitats ----
# PIE CHART PER TAXA ####

#
















