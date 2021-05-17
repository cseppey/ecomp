#####
# ecomp prep
#####

rm(list=ls())

require(foreach)
require(doSNOW)

#---
# cluster
cl <- makeSOCKcluster(4)
registerDoSNOW(cl)


#---
# download ####

dir_in <- 'Projets/Ecomp/stat/in/200821/'
dir_out <- 'Projets/Ecomp/stat/out/200821/'
dir_save <- paste0(dir_out, 'saves/') 

dir.create(dir_save, showWarnings=F, recursive=T, mode='775')

funct <- read.table(paste0(dir_in, 'grp_funct.csv'), h=T, sep=',')
                          
recalc <- F

#---
file1 <- paste0(dir_save, 'mr_ass.Rdata')
if(file.exists(file1) & recalc == F){
  load(file1)
} else {
  # to make the chunks, check chunking.sh
  files <- sapply(strsplit(list.files(paste0(dir_in, 'chunks'), '.mr', full.names=T), '.', fixed=T), '[[', 1)
  
  system.time(mr_tot  <- read.table(paste0(files[1], '.mr'),  h=T))
  system.time(ass_tot <- read.table(paste0(files[1], '.ass'), sep='\t'))
  
  for(i in files[2:length(files)]){
    print(i)
    mr_tot <-  cbind.data.frame(mr_tot,  read.table(paste0(i, '.mr'), h=T))
    ass_tot <- rbind.data.frame(ass_tot, read.table(paste0(i, '.ass'), sep='\t'))
  }
  
  # arrange assignation
  names(ass_tot) <- c('OTU_id','pid','e_value','taxo','GB_id','seq')
  row.names(ass_tot) <- paste0('X', sapply(strsplit(as.character(ass_tot$OTU_id), '_'), '[[', 1))
  
  ass_tot$seq <- as.character(ass_tot$seq)
  ass_tot$taxo <- as.character(ass_tot$taxo)
  
  save(mr_tot, ass_tot, file=file1)    
}

# parse in function of length or if assignation too scrap
ass_tot$seq_lgt <- sl <- nchar(ass_tot$seq)

lim_lgt <- c(100,250)

hist(sl, breaks=100)
abline(v=lim_lgt)

lst_ind <- list(full=1:length(sl), 
                pars_lgt=which(sl >= lim_lgt[1] & sl <= lim_lgt[2]),
                pars_miss_tax=which(ass_tot$taxo != ''),
                pars_60=which(ass_tot$pid >= 60))

#---
if(recalc){
  for(i in names(lst_ind)){
    print(i)
    
    cairo_ps(paste0(dir_out, 'pid_vs_evalue_', i, '.eps'))
    
    ass <- ass_tot[lst_ind[[i]][is.na(ass_tot$pid) == F],]
    ssl <- sl[lst_ind[[i]]]
    
    plot(ass$seq_lgt, ass$pid, pch=19, cex=0.1)
    usr <- par('usr')
    
    dev.off()
  }
  
  cairo_ps(paste0(dir_out, 'seq_lgt.eps'))
  
  plot(log(table(ass_tot$seq_lgt)), type='l', ylab='log nb OTUs', xaxt='n')
  axis(1)
  abline(v=lim_lgt)
  
  dev.off()
}

#---
ind_ok <- intersect(intersect(lst_ind$pars_lgt, lst_ind$pars_miss_tax), lst_ind$pars_60)

mr_sl <- mr_tot[,ind_ok]
ass_sl <- droplevels(ass_tot[ind_ok,])

# taxo
taxo_sl <- as.data.frame(t(sapply(strsplit(ass_sl$taxo, ';'), '[', 1:8)))
dimnames(taxo_sl) <- list(row.names(ass_sl),c('reign','phylum','division','class','order','family','genus','species'))

levels(taxo_sl$phylum)[levels(taxo_sl$phylum) == 'Hacrobia'] <- 'Eukaryota_X'

# functions
file2 <- paste0(dir_save, 'fct.Rdata')
if(file.exists(file2) & recalc == F){
  load(file2)
} else {
  chunks <- split(ass_sl, 0:(nrow(ass_sl)-1) %/% 10000)
  print(length(chunks))
  
  fct <- foreach(i = chunks, .verbose=T, .combine='rbind') %dopar% {
    f <- NULL
    for(j in as.character(i$taxo)){
      cnt <- 0
      for(k in strsplit(j, ';')[[1]]){
        ind <- which(k == funct$taxon)
        if(length(ind)){
          f <- rbind(f, c(as.character(funct$funct[ind]), as.character(k)))
          cnt <- 1
          break
        }
      }
      if(cnt == 0){
        f <- rbind(f, c('unknown',j))
      }
    }
    return(f)
  }

  save(fct, file=file2)    
}

fct2 <- as.factor(fct[,1])
fct2 <- factor(fct2, levels=levels(fct2)[c(1,3,2,4)])
ass_sl <- cbind.data.frame(ass_sl, fct=fct2)

# env
rn <- row.names(mr_sl)
env <- data.frame(env=c('freshwater','marine','soil')[as.numeric(substr(rn, 1, 1))],
                  ech=substr(rn, 1, 4), sub=paste0(substr(rn, 1, 1), substr(rn, 5, nchar(rn))))
row.names(env) <- rn

#---
env$env[grepl('1021|1022|1023|1024', env$ech)] <- 'soil'

#---
# removal of the fraction

ind_nfrac <- env$sub %in% c('1','2e','3')

env <- env[ind_nfrac,]
mr_sl <- mr_sl[ind_nfrac,]
mr_sl <- mr_sl[,colSums(mr_sl) != 0]
ass_sl <- ass_sl[names(mr_sl),]
taxo_sl <- taxo_sl[names(mr_sl),]

#---
file3 <- paste0(dir_save, 'prep.Rdata')
save(env, mr_sl, ass_sl, taxo_sl, file=file3)

# display all taxo ####

tax_select <- c('Amoebozoa','Ciliophora')
taxo_sel <- matrix(unlist(strsplit(sort(unique(apply(taxo_sl[grep(paste(tax_select, collapse='|'),
                                                                  ass_sl$taxo),-8], 1, function(x) paste(x, collapse=';')))), 
                                   split=';')), ncol=7, byrow=T)
names(taxo_sel) <- names(taxo_sl)[-8]

write.table(taxo_sel, paste0(dir_out, 'tax_selec.csv'), quote=F, sep='\t', row.names=F, col.names=F)






