#This file is to rerun the previous data with the same functions that eDNA is being run with to limit any introduced error in analysis

# load libraries
library("ape")
library("pegas")
library("reshape2")
library("seqinr")
library("haplotypes")
library("adegenet")
library("poppr")

# set working directory
setwd("/Users/taylorely/Documents/Grad_Work/Chapter1/existing_data_species/short_previous_haplotypes")

#######################################################################################3
#make functions
############
#this is to take a fasta file and convert it into a df with each line being an individual
reformat<-function(filepath){
  prev_species<-read.fasta(filepath, as.string = TRUE)
  prev_species<-as.data.frame(do.call(rbind, prev_species))
  prev_species$samples<-row.names(prev_species)
  colnames(prev_species)<-c("sequences","samples")
  split<-data.frame(do.call("rbind", strsplit(as.character(prev_species$samples), "_", fixed = TRUE)))
  combined<-cbind(prev_species,split)
  colnames(combined)<-c("sequences","samples","site","haplo_number","frequency")
  seqs<-rep(combined$sequences, times=combined$frequency)
  sites<-rep(combined$site, times=combined$frequency)
  species_finished<<-data.frame(cbind(seqs,sites))
}
#####
#get pairwise phist
previous_haplotypes_pairphi<-function(fasta,species,bp,order){
  folder_name<-"/Users/taylorely/Documents/Grad_Work/Chapter1/existing_data_species/short_previous_haplotypes/HaplotypesPairPhi"
  x <- as.list(fasta$seqs)
  site<-fasta$sites
  site<-as.factor(site)
  x <- as.DNAbin(ape::as.alignment(x))
  x <- as.dna(x)
  #calculate pairwise phist
  pst<-pairPhiST(x, site, nperm=1000, negatives=FALSE, showprogbar=FALSE) #calculates distance in code
  #remove pops that it did not occur in
  i_order<-order
  pops<-unique(fasta$sites)
  nopops<-setdiff(i_order, pops)
  if(length(nopops)!=0){
    i<-1
    for (i in 1:length(nopops)){
      i_order<-i_order[(i_order != nopops[i])]
    }}
  #first start with reordering phist
  pst_test<-as.matrix(pst$PhiST) #convert to matrix
  pst_test[upper.tri(pst_test)] <- t(pst_test)[upper.tri(pst_test)] #copy so symmetrical 
  phist_reordered <- reorder_mat(mat = pst_test, order = i_order) #reorder
  #second reorder p-value
  pst_test2<-as.matrix(pst$p) #convert to matrix
  pst_test2[upper.tri(pst_test2)] <- t(pst_test2)[upper.tri(pst_test2)] #copy so symmetrical 
  pvalue_reordered <- reorder_mat(mat = pst_test2, order = i_order) #reorder
  #now merge the two halves together
  phist_reordered[upper.tri(phist_reordered)] <- pvalue_reordered[upper.tri(pvalue_reordered)]
  write.csv(phist_reordered, file = paste(folder_name, species,"_",bp, "_island_pretty_haplotypes_pst.csv", sep=""))
}
# get a nice table with all pairwise phist information
phist_table_previous<-function(folderinput,folderoutput,runname){
  files<-list.files(folderinput)
  phisttable<-data.frame(matrix(ncol = 10, nrow = 0))
  i<-1
  for(i in 1:length(files)){
    #load a single phist table
    #grab phist values
    phistspecies<-read.csv(file= paste0(folderinput,files[i]), stringsAsFactors = F)
    numbcol<-ncol(phistspecies)
    phistnumb<-(numbcol-1)/2
    phistspecies2<-phistspecies[,1:phistnumb]
    phistspecies2<-melt(phistspecies2, id.vars="X")
    variable<-phistspecies2$variable
    names2 <- as.data.frame(str_split_fixed(variable, '.', 7))
    phistspecies2$variable<-names2$V7
    phistspecies2<-na.omit(phistspecies2)
    #grab p-values
    pnumb<-phistnumb+2
    phistspecies3<-phistspecies[,c(1,pnumb:numbcol)]
    phistspecies3<-melt(phistspecies3, id.vars="X")
    phistspecies3<-na.omit(phistspecies3)
    variable2<-phistspecies3$variable
    names3 <- as.data.frame(str_split_fixed(variable2, '.', 3))
    phistspecies3$variable<-names3$V3
    #correct p-values FDR
    phistspecies3$fdr<-p.adjust(p=phistspecies3$value, method = "fdr")
    #correct p-values bonforroni
    phistspecies3$bon<-p.adjust(p=phistspecies3$value, method = "bonferroni")
    #rename for p-value
    colnames(phistspecies3)<-c("X","variable","p-value","fdr p-value","bon p-value")
    #merge together
    phistall <- merge(phistspecies2,phistspecies3,by=c("X","variable"), all = TRUE)
    Species<-files[i]
    phistall<-cbind(Species,phistall)
    #if significant
    phistall$sigp <- with(phistall, ifelse(phistall$`p-value` < 0.05, "yes", "no"))
    phistall$sigfdr <- with(phistall, ifelse(phistall$`fdr p` < 0.05, "yes", "no"))
    phistall$sigbon <- with(phistall, ifelse(phistall$`bon p`< 0.05, "yes", "no"))
    #bind to table
    phisttable<-rbind(phisttable,phistall)
  }
  colnames(phisttable)<-c("Species","Population 1","Population 2","PhiST value","p-value","fdr p-value","bonferroni p-value","sig p-value","sig fdr","sig bon")
  write.csv(phisttable, file = paste(folderoutput,"/",runname,"_phist_table_manyspecies.csv", sep=""))
}

#####
#calculate genetic diversity statistics
statstable<-function(fasta,species,bp){
    folder_name<-"/Users/taylorely/Documents/Grad_Work/Chapter1/existing_data_species/short_previous_haplotypes/BasicStats/"
    pops<-unique(fasta$sites)
    summary_stats<-data.frame(matrix(ncol = 15, nrow = 0))
    x <- as.list(fasta$seqs)
    x <- as.DNAbin(ape::as.alignment(x))
    h <- pegas::haplotype(x)
    N<-nrow(fasta)
    H<-length(unique(fasta$seqs))
    nd<-nuc.div(h,variance=TRUE) 
    nd_var<-paste(round(nd[1],digits = 6),round(nd[2],digits = 6),sep = " ± ")
    hd<-hap.div(h,variance=TRUE)
    hd_var<-paste(round(hd[1],digits = 4),round(hd[2],digits = 6),sep = " ± ")
    tajima<-tajima.test(x) 
    R2<-R2.test(x,B = 1000,quiet = TRUE,plot = FALSE)
    population<-"overall"
    bp_length<-nchar(fasta$seqs[1])
    basicstats<-c(species,population,bp_length,N,H,hd_var,nd_var,round(tajima$D,digits = 3),round(R2$R2[1],digits = 3),tajima$Pval.normal,R2$P.val[1],tajima$Pval.beta,hd[1],nd[1])
    summary_stats<-rbind(summary_stats,basicstats)
    i<-1
    for (i in 1:length(pops)){
    temp <- fasta[fasta$sites==pops[i],]
    x <- as.list(temp$seqs)
    x <- as.DNAbin(ape::as.alignment(x))
    h <- pegas::haplotype(x)
    N<-nrow(temp)
    H<-length(unique(temp$seqs))
    nd<-nuc.div(h,variance=TRUE) 
    nd_var<-paste(round(nd[1],digits = 6),round(nd[2],digits = 6),sep = " ± ")
    hd<-hap.div(h,variance=TRUE)
    hd_var<-paste(round(hd[1],digits = 4),round(hd[2],digits = 6),sep = " ± ")
    tajima<-tajima.test(x) 
    population<-temp$sites[1]
    bp_length<-nchar(temp$seqs[1])
    R2<-R2.test(x,B = 1000,quiet = TRUE,plot = FALSE)
    basicstats<-c(species,population,bp_length,N,H,hd_var,nd_var,round(tajima$D,digits = 3),round(R2$R2[1],digits = 3),tajima$Pval.normal,R2$P.val[1],tajima$Pval.beta,hd[1],nd[1])
    summary_stats<-rbind(summary_stats,basicstats)
    }
    colnames(summary_stats)<-c("Species","Population","bp length","Sample Size","Number of Haplotypes","hap div var","nuc div var","Tajima's D","R2","Tajima's D p-value","R2 p-value","Tajima's D p-value beta","haplotype diversity","nucleotide diversity")
    write.csv(summary_stats, file = paste(folder_name,species,"_",bp, "_basic_seq_stats.csv", sep=""))
}

#####
#make mismatch graphs
MMD_previous<-function(fasta,species,bp){
    y <- as.list(fasta$seqs)
      y <- as.DNAbin(ape::as.alignment(y))
      d <- dist.dna(y, "N")
      maxdist<-max(d)
      if (maxdist > 0){
        if(maxdist < 4) {maxdist<-4}
        lcol = c("lightblue", "darkred")
        lty = c(1, 1)
        jpeg(file=paste(species,"_",bp,"_all_islands_MMD.jpeg",sep = ""), width=11, height=10, units = "in", res=400)
        par(mar=c(5.1, 5.1, 4.1, 2.1))
        h <- hist(d, xlab = "Number of Pairwise Differences",ylab="Frequency", main = "",xaxt = "n",cex.axis = 2, cex.lab = 2,  ylim = c(0,1), freq = FALSE,breaks = c(0:maxdist))
        axis(1, at = c(0:maxdist),cex.axis = 2)
        dd <- density(d, bw = 2)
        lines(dd, col = lcol[1], lty = lty[1],lwd = 4)
        ## by David Winter:
        theta <- mean(d)
        upper <- ceiling(max(d))
        e <- sapply(0:upper, function(i) theta^i / (theta + 1)^(i + 1))
        lines(e, col = lcol[2], lty = lty[2],lwd = 4)
        rug(d)
        psr <- par("usr")
        xx <- psr[2]/2
        yy <- psr[4] * (0.5 + 0.5/par("plt")[4])
        legend(xx, yy, c("Empirical", "Stable expectation"),
               lty = lty, col = lcol, bg = "white", lwd = 4,bty = "n",
               xjust = 0.5, yjust = 0.5, horiz = TRUE, xpd = TRUE,cex=2)
        #legend("topleft", c("Empirical", "Stable expectation"),
        #       lty = 1, col = lcol, bg = "white", bty = "n")
        invisible(list(histogram = h, empirical.density = dd, expected.curve = e))
        dev.off()
      pops<-unique(fasta$sites)
      for (i in 1:length(pops)){
        temppop <- fasta[fasta$sites==pops[i],]
        y <- as.list(temppop$seqs)
        y <- as.DNAbin(ape::as.alignment(y))
        d <- dist.dna(y, "N")
        maxdist<-max(d)
        if(maxdist>0){
          if(maxdist < 4) {maxdist<-4}
          lcol = c("lightblue", "darkred")
          lty = c(1, 1)
          jpeg(file=paste(species,"_",bp,"_",temppop$sites[1],"_MMD.jpeg",sep = ""), width=11, height=10, units = "in", res=400)
          par(mar=c(5.1, 5.1, 4.1, 2.1))
          h <- hist(d, xlab = "Number of Pairwise Differences",ylab="Frequency", main = "",xaxt = "n",cex.axis = 2, ylim = c(0,1), cex.lab = 2, freq = FALSE,breaks = c(0:maxdist))
          axis(1, at = c(0:maxdist),cex.axis = 2)
          dd <- density(d, bw = 2)
          lines(dd, col = lcol[1], lty = lty[1],lwd = 4)
          ## by David Winter:
          theta <- mean(d)
          upper <- ceiling(max(d))
          e <- sapply(0:upper, function(i) theta^i / (theta + 1)^(i + 1))
          lines(e, col = lcol[2], lty = lty[2],lwd = 4)
          rug(d)
          psr <- par("usr")
          xx <- psr[2]/2
          yy <- psr[4] * (0.5 + 0.5/par("plt")[4])
          legend(xx, yy, c("Empirical", "Stable expectation"),
                 lty = lty, col = lcol, bg = "white", lwd = 4,bty = "n",
                 xjust = 0.5, yjust = 0.5, horiz = TRUE, xpd = TRUE,cex=2)
          #legend("topleft", c("Empirical", "Stable expectation"),
          #       lty = 1, col = lcol, bg = "white", bty = "n")
          invisible(list(histogram = h, empirical.density = dd, expected.curve = e))
          dev.off()
        }
      }
    } 
} 

#####
#Convert fasta files to right format for PGDSpider, this is for arlequin
fastafilesperspecies<-function(datafile,species,length,folder){
  #set unique number so that each name is unique
  datafile$ID<-row.names(datafile)
  seqs_fasta<-as.list(datafile$seqs)
  #set names ex: >OTU1 population:Oahu haplo_1
  names_fasta<- as.list(paste(datafile$ID," population:",datafile$sites," ",sep=""))
  #write as fasta files
  write.fasta(sequences=seqs_fasta,names=names_fasta,file.out=paste(folder,species,"_",length,".fasta",sep = ""))
}
#convert to fasta per pop for each species
fastabypop<-function(datafile, species, length, folder){
  pops<-unique(datafile$sites)
  for (i in 1:length(pops)){
    temp2 <- datafile[datafile$sites==pops[i],]
    #set sequences
    seqs_fasta<-as.list(temp2$seqs)
    #set names ex: >OTU1 population:Oahu haplo_1
    names_fasta<- as.list(paste(species,"_",length," population:",temp2$sites, sep=""))
    #write as fasta files
    write.fasta(sequences=seqs_fasta,names=names_fasta,file.out=paste(folder,species,"_",length,"_",temp2$sites[1],".fasta",sep = ""))
  }
}
#fasta to nexus file for DNASP
nexusfromfasta<-function(folderinput,folderoutput){
  files<-list.files(folderinput)
  i<-1
  for(i in 1:length(files)){
    fastaspecies<-read.fasta(file= paste0(folderinput,"/",files[i]))
    write.nexus.data(fastaspecies, file=paste(folderoutput,files[i],"_nexus.nex",sep = ""))
  }
}

#################################################################################3
longer<-"orginal"
shorter<-"short"
mhiorder<-c("Kauai","Oahu","Maui","Hawaii")
hkorder<-c("Kauai","Oahu","Maui","Kona","Hilo")
################

#lutjanus kasmira
lutkas<-reformat("Lutjanus_kasmira_original.fasta")
lutkasname<-"Lutjanus_kasmira"
fastafilesperspecies(datafile = lutkas,
                     runname = lutkasname,
                     folder = "fastaredo/")
lutkas_pegas_short<-previous_haplotypes_pairphi(lutkas,lutkasname,longer,mhiorder)
lutkas_stats_long<-statstable(lutkas,lutkasname,longer)
MMD_previous(lutkas,lutkasname,longer)
#fastabypop<-function(datafile, species, length, folder)
fastabypop(datafile=lutkas,
           species=lutkasname,
           length=longer,
           folder = "fastaredo/")


abuabd_long<-reformat("Abuabd_long.fasta")
abuabd_short<-reformat("Abuabd_short.fasta")
abuabd<-"Abudefduf_abdominalis"
fastafilesperspecies(datafile = abuabd_short,
                     species = abuabd,
                     length = shorter,
                     folder = "fastaredo/")
fastabypop(datafile=abuabd_short,
           species=abuabd,
           length=shorter,
           folder = "fastaredo/")
abuabd_pegas_short<-previous_haplotypes_pairphi(abuabd_short,abuabd,shorter,mhiorder)
abuabd_pegas_long<-previous_haplotypes_pairphi(abuabd_long,abuabd,longer,mhiorder)
abuabd_stats_short<-statstable(abuabd_short,abuabd,shorter)
abuabd_stats_long<-statstable(abuabd_long,abuabd,longer)
MMD_previous(fasta = abuabd_short,
             species = abuabd,
             bp = shorter)
MMD_previous(fasta = abuabd_long,
             species = abuabd,
             bp = longer)

abuvai_short<-reformat("Abd_vaigiensis_kim_tissue_short_MHI.fasta")
abuvai_long<-reformat("Abd_vaigiensis_kim_tissue_long_MHI.fasta")
abuvai<-"Abudefduf_vaigiensis_MHI_nonegatives"
fastafilesperspecies(datafile = abuvai_short,
                     species = abuvai,
                     length = shorter,
                     folder = "fastaredo/")
fastabypop(datafile=abuvai_long,
           species=abuvai,
           length=longer,
           folder = "fastaredo/")
abuvai_pegas_short<-previous_haplotypes_pairphi(abuvai_short,abuvai,shorter,mhiorder)
abuvai_pegas_long<-previous_haplotypes_pairphi(abuvai_long,abuvai,longer,mhiorder)
abuvai_stats_short<-statstable(abuvai_short,abuvai,shorter)
abuvai_stats_long<-statstable(abuvai_long,abuvai,longer)
MMD_previous(fasta = abuvai_short,
             species = abuvai,
             bp = shorter)
MMD_previous(fasta = abuvai_long,
             species = abuvai,
             bp = longer)

acanig_long<-reformat("Acanig_longer.fasta")
acanig_short<-reformat("Acanig_short.fasta")
acanig<-"Acanthurus_nigrofuscus"
fastafilesperspecies(datafile = acanig_short,
                     species = acanig,
                     length = shorter,
                     folder = "fastaredo/")
fastabypop(datafile=acanig_short,
           species=acanig,
           length=shorter,
           folder = "fastaredo/")
acanig_pegas_short<-previous_haplotypes_pairphi(acanig_short,acanig,shorter,mhiorder)
acanig_pegas_long<-previous_haplotypes_pairphi(acanig_long,acanig,longer,mhiorder)
acanig_stats_short<-statstable(acanig_short,acanig,shorter)
acanig_stats_long<-statstable(acanig_long,acanig,longer)
MMD_previous(fasta = acanig_short,
             species = acanig,
             bp = shorter)
MMD_previous(fasta = acanig_long,
             species = acanig,
             bp = longer)


acanigroris<-"Acanthurus_nigroris"
acanigroris_long<-reformat("Acanthurus_nigroris_longer.fasta")
acanigroris_short<-reformat("Acanthurus_nigroris_short.fasta")
fastafilesperspecies(datafile = acanigroris_short,
                     species = acanigroris,
                     length=shorter,
                     folder = "fasta_PGD/")
fastabypop(datafile=acanigroris_short,
           species=acanigroris,
           length=shorter,
           folder = "fastaredo/")
acanigroris_pegas_short<-previous_haplotypes_pairphi(acanigroris_short,acanigroris,shorter,mhiorder)
acanigroris_pegas_long<-previous_haplotypes_pairphi(acanigroris_long,acanigroris,longer,mhiorder)
acanigroris_stats_short<-statstable(acanigroris_short,acanigroris,shorter)
acanigroris_stats_long<-statstable(acanigroris_long,acanigroris,longer)
MMD_previous(fasta = acanigroris_short,
             species = acanigroris,
             bp = shorter)
MMD_previous(fasta = acanigroris_long,
             species = acanigroris,
             bp = longer)


#doesn't overlap so no short region
acaoli_long<-reformat("Acaoli_long_doesntoverlap.fasta")
acaoli<-"Acanthurus_olivaceus"
acaoli_pegas_long<-previous_haplotypes_pairphi(acaoli_long,acaoli,longer)
acaoli_stats_long<-statstable(acaoli_long,acaoli,longer)
MMD_previous(fasta = acaoli_long,
             species = acaoli,
             bp = longer)


cepharg_long<-reformat("Cepharg_longer.fasta")
cepharg_short<-reformat("Cepharg_short.fasta")
cepharg<-"Cephalopholis_argus"
fastafilesperspecies(datafile = cepharg_long,
                     species = cepharg,
                     length=longer,
                     folder = "fasta_PGD/")
fastabypop(datafile=cepharg_long,
           species=cepharg,
           length=longer,
           folder = "fastaredo/")
cepharg_pegas_short<-previous_haplotypes_pairphi(cepharg_short,cepharg,shorter,mhiorder)
cepharg_pegas_long<-previous_haplotypes_pairphi(cepharg_long,cepharg,longer,mhiorder)
cepharg_stats_short<-statstable(cepharg_short,cepharg,shorter)
cepharg_stats_long<-statstable(cepharg_long,cepharg,longer)
MMD_previous(fasta = cepharg_short,
             species = cepharg,
             bp = shorter)
MMD_previous(fasta = cepharg_long,
             species = cepharg,
             bp = longer)


chafre_long<-reformat("Chafre_longer.fasta")
chafre_short<-reformat("Chafre_short.fasta")
charfre<-"Chaetodon_fre"
fastabypop(datafile=chafre_short,
           species=charfre,
           length=shorter,
           folder = "fastaredo/")
charfre_stats_short<-statstable(chafre_short,charfre,shorter)
charfre_stats_long<-statstable(chafre_long,charfre,longer)
MMD_previous(fasta = chafre_short,
             species = charfre,
             bp = shorter)
MMD_previous(fasta = chafre_long,
             species = charfre,
             bp = longer)


chamul_long<-reformat("Chamul_longer.fasta")
chamul_short<-reformat("Chamul_short.fasta")
charmul<-"Chaetodon_multicinctus"
fastafilesperspecies(datafile = chamul_short,
                     species = charmul,
                     length=shorter,
                     folder = "fastaredo/")
fastabypop(datafile=chamul_short,
           species=charmul,
           length=shorter,
           folder = "fastaredo/")
charmul_pegas_short<-previous_haplotypes_pairphi(chamul_short,charmul,shorter,mhiorder)
charmul_pegas_long<-previous_haplotypes_pairphi(chamul_long,charmul,longer,mhiorder)
charmul_stats_short<-statstable(chamul_short,charmul,shorter)
charmul_stats_long<-statstable(chamul_long,charmul,longer)
MMD_previous(fasta = chamul_short,
             species = charmul,
             bp = shorter)
MMD_previous(fasta = chamul_long,
             species = charmul,
             bp = longer)

chamil_long<-reformat("Chamil_longer.fasta")
chamil_short<-reformat("Chamil_short.fasta")
chamil<-"Chaetodon_miliaris"
fastafilesperspecies(datafile = chamil_short,
                     species = chamil,
                     length=shorter,
                     folder = "fastaredo/")
fastabypop(datafile=chamil_short,
           species=chamil,
           length=shorter,
           folder = "fastaredo/")
charmil_pegas_short<-previous_haplotypes_pairphi(chamil_short,charmil,shorter,mhiorder)
charmil_pegas_long<-previous_haplotypes_pairphi(chamil_long,charmil,longer,mhiorder)
charmil_stats_short<-statstable(chamil_short,charmil,shorter)
charmil_stats_long<-statstable(chamil_long,charmil,longer)
MMD_previous(fasta = chamil_short_limited,
             species = charmil_limited,
             bp = shorter)
MMD_previous(fasta = chamil_long_limited,
             species = charmil_limited,
             bp = longer)
#only Hawaii and Kauai
chamil_long_limited<-chamil_long[chamil_long$sites=="Kauai" | chamil_long$sites=="Hawaii",]
chamil_short_limited<-chamil_short[chamil_short$sites=="Kauai" | chamil_short$sites=="Hawaii",]
charmil_limited<-"Chaetodon_miliaris_limited"
statstable(chamil_short_limited,charmil_limited,shorter)
statstable(chamil_long_limited,charmil_limited,longer)


chrvan<-"Chromis_vanderbilti"
chrvan_long<-reformat("chrvan_longer.fasta")
chrvan_short<-reformat("chrvan_shorter.fasta")
chrvan_long$ID<-row.names(chrvan_long)
seqs_fasta<-as.list(chrvan_long$seqs)
#set names ex: >OTU1 population:Oahu haplo_1
names_fasta<- as.list(paste(chrvan_long$ID," population:",chrvan_long$sites," ",sep=""))
#write as fasta files
write.fasta(sequences=seqs_fasta,names=names_fasta,file.out=paste(folder,"Chromis_vanderbilti_long",".fasta",sep = ""))
fastafilesperspecies(datafile=chrvan_short,
           species=chrvan,
           length=shorter,
           folder = "fastaredo/")
fastabypop(datafile=chrvan_short,
           species=chrvan,
           length=shorter,
           folder = "fastaredo/")
previous_haplotypes_pairphi(chrvan_short,chrvan,shorter,mhiorder)
previous_haplotypes_pairphi(chrvan_long,chrvan,longer,mhiorder)
charmil_stats_short<-statstable(chrvan_short,chrvan,shorter)
charmil_stats_long<-statstable(chrvan_long,chrvan,longer)
MMD_previous(fasta = chrvan_short,
             species = chrvan,
             bp = shorter)
MMD_previous(fasta = chrvan_long,
             species = chrvan,
             bp = longer)


ctestr_long<-reformat("Ctestr_long.fasta")
ctestr_short<-reformat("Ctestr_short.fasta")
ctestr<-"Ctenochaetus flavescens"
fastabypop(datafile=ctestr_short,
           species=ctestr,
           length=shorter,
           folder = "fastaredo/")
ctestr_pegas_short<-previous_haplotypes_pairphi(ctestr_short,ctestr,shorter)
ctestr_pegas_long<-previous_haplotypes_pairphi(ctestr_long,ctestr,longer)
ctestr_stats_short<-statstable(ctestr_short,ctestr,shorter)
ctestr_stats_long<-statstable(ctestr_long,ctestr,longer)
MMD_previous(fasta = ctestr_short,
             species = ctestr,
             bp = shorter)
MMD_previous(fasta = ctestr_long,
             species = ctestr,
             bp = longer)


#doesn't overlap much so no shorter region
lutful_long<-reformat("Lutjanus_fulvus_longer.fasta")
lutful<-"Lutjanus_fulvus"
fastafilesperspecies(datafile = lutful_long,
                     species = lutful,
                     length = longer,
                     folder = "fasta_PGD/")
fastabypop(datafile=lutful_long,
           species=lutful,
           length=longer,
           folder = "fastaredo/")
lutful_pegas_long<-previous_haplotypes_pairphi(lutful_long,lutful,longer,mhiorder)
lutful_stats_long<-statstable(lutful_long,lutful,longer)
MMD_previous(fasta = lutful_long_limited,
             species = lutful_limited,
             bp = longer)
#only Hawaii, Oahu, Kauai
lutful_long_limited<-lutful_long[lutful_long$sites=="Hawaii" | lutful_long$sites=="Kauai" | lutful_long$sites=="Oahu",]
lutful_limited<-"Lutjanus_fulvus_limited"
lutful_stats_long<-statstable(lutful_long_limited,lutful_limited,longer)
fastafilesperspecies(datafile = lutful_long_limited,
                     species = lutful_limited,
                     length = longer,
                     folder = "fasta_PGD/")
fastabypop(datafile=lutful_long_limited,
           species=lutful_limited,
           length=longer,
           folder = "fastaredo/")


mullflav_short<-reformat("Mulflav_short.fasta")
mullflav_longer<-reformat("Mulflav_longer.fasta")
mulloidichthys_flav<-"mulloidichthys_flav"
fastafilesperspecies(datafile = mullflav_short,
                     species = mulloidichthys_flav,
                     length = shorter,
                     folder = "fastaredo/")
fastabypop(datafile=mullflav_short,
           species=mulloidichthys_flav,
           length=shorter,
           folder = "fastaredo/")
mull_flav_phist<-previous_haplotypes_pairphi(mullflav_short,mulloidichthys_flav,shorter,mhiorder)
mull_flav_longer_phist<-previous_haplotypes_pairphi(mullflav_longer,mulloidichthys_flav,longer,mhiorder)
mull_flav_stats_short<-statstable(mullflav_short,mulloidichthys_flav,shorter)
mull_flav_stats_longer<-statstable(mullflav_longer,mulloidichthys_flav,longer)
MMD_previous(fasta = mullflav_short,
             species = mulloidichthys_flav,
             bp = shorter)
MMD_previous(fasta = mullflav_longer,
             species = mulloidichthys_flav,
             bp = longer)


zebflav_long<-reformat("Zebrasoma_flavescens_longer.fasta")
zebflav_short<-reformat("Zebrasoma_flavescens_short.fasta")
zebflav<-"Zebrasoma_flavescens_konahilo_long"
fastafilesperspecies(datafile = zebflav_short,
                     species = zebflav,
                     length=shorter,
                     folder = "fastaredo/")
fastabypop(datafile=zebflav_short,
           species=zebflav,
           length=shorter,
           folder = "fastaredo/")
zebflav_pegas_short<-previous_haplotypes_pairphi(zebflav_short,zebflav,shorter)
zebflav_pegas_long<-previous_haplotypes_pairphi(zebflav_long,zebflav,longer)
zebflav_stats_short<-statstable(zebflav_short,zebflav,shorter)
zebflav_stats_long<-statstable(zebflav_long,zebflav,longer)
MMD_previous(fasta = zebflav_short,
             species = zebflav,
             bp = shorter)
MMD_previous(fasta = zebflav_long,
             species = zebflav,
             bp = longer)
zebflav_long <- data.frame(lapply(zebflav_long, function(x) {gsub("Hilo", "Hawaii", x)}))
zebflav_long <- data.frame(lapply(zebflav_long, function(x) {gsub("Kona", "Hawaii", x)}))
zebflav_short <- data.frame(lapply(zebflav_short, function(x) {gsub("Hilo", "Hawaii", x)}))
zebflav_short <- data.frame(lapply(zebflav_short, function(x) {gsub("Kona", "Hawaii", x)}))
zebflav<-"Zebrasoma_flavescens"
zebflav_pegas_short<-previous_haplotypes_pairphi(zebflav_short,zebflav,shorter,hkorder)
zebflav_pegas_long<-previous_haplotypes_pairphi(zebflav_long,zebflav,longer,hkorder)
zebflav_stats_short<-statstable(zebflav_short,zebflav,shorter)
zebflav_stats_long<-statstable(zebflav_long,zebflav,longer)

#no short because different region
triobe<-"Triaenodon obesus"
triobe_long<-reformat("triodon_tissue_long.fasta")
fastafilesperspecies(datafile = triobe_long,
                     species = triobe,
                     length=longer,
                     folder = "fastaredo/")
fastabypop(datafile=triobe_long,
           species=triobe,
           length=longer,
           folder = "fastaredo/")
triobe_stats_long<-statstable(triobe_long,triobe,longer)
MMD_previous(fasta = triobe_long,
             species = triobe,
             bp = longer)
previous_haplotypes_pairphi(triobe_long,triobe,longer,mhiorder)

stenella<-"Stenella longirostris"
stenella_long<-reformat("stenella_tissue_long.fasta")
fastafilesperspecies(datafile = stenella_long,
                     species = stenella,
                     length=longer,
                     folder = "fastaredo/")
fastabypop(datafile=stenella_long,
           species=stenella,
           length=longer,
           folder = "fastaredo/")
stenella_stats_long<-statstable(stenella_long,stenella,longer)
MMD_previous(fasta = stenella_long,
             species = stenella,
             bp = longer)
previous_haplotypes_pairphi(stenella_long,stenella,longer,mhiorder)

carmel<-"Caranx melampygus"
carmel_long<-reformat("Carmel_tissue_long.fasta")
carmel_stats_long<-statstable(carmel_long,carmel,longer)
fastafilesperspecies(datafile = carmel_long,
                     species = carmel,
                     length=longer,
                     folder = "fastaredo/")
fastabypop(datafile=carmel_long,
           species=carmel,
           length=longer,
           folder = "fastaredo/")
MMD_previous(fasta = carmel_long,
             species = carmel,
             bp = longer)
previous_haplotypes_pairphi(carmel_long,carmel,longer,mhiorder)

carign<-"Caranx ignobilis"
carign_long<-reformat("Carign_tissue_long.fasta")
fastafilesperspecies(datafile = carign_long,
                      species = carign,
                     length=longer,
                     folder = "fastaredo/")
fastabypop(datafile=carign_long,
           species=carign,
           length=longer,
           folder = "fastaredo/")
carign_stats_long<-statstable(carign_long,carign,longer)
MMD_previous(fasta = carign_long,
             species = carign,
             bp = longer)
previous_haplotypes_pairphi(carign_long,carign,longer,mhiorder)

#nexusfromfasta<-function(folderinput,folderoutput)
nexusfromfasta(folderinput = "fastaredo/",
               folderoutput = "nexusredo/")
##################################

#phist big table
phist_table_previous(folderinput="/Users/taylorely/Documents/Grad_Work/Chapter1/existing_data_species/short_previous_haplotypes/HaplotypesPairPhi/",
                     folderoutput="/Users/taylorely/Documents/Grad_Work/Chapter1/existing_data_species/short_previous_haplotypes",
                     runname="previous_tissue")



