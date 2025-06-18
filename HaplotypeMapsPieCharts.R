library("ape")
library("pegas")
library("ggplot2")
library("haplotypes")
library("seqinr") 
library("plyr")
library("dplyr")
library("stringr")
library("RColorBrewer")
library("graph4lg")

setwd("/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/Hapmaps")

#functions
reformat<-function(filepath,speciesname,removefirst,sangerlength){
  prev_species<-read.fasta(filepath, as.string = TRUE)
  prev_species<-as.data.frame(do.call(rbind, prev_species))
  prev_species$samples<-row.names(prev_species)
  colnames(prev_species)<-c("sequences","samples")
  split<-data.frame(do.call("rbind", strsplit(as.character(prev_species$samples), "_", fixed = TRUE)))
  combined<-cbind(prev_species,split)
  colnames(combined)<-c("sequences","samples","site","haplo_number","frequency")
  sequences<-rep(combined$sequences, times=combined$frequency)
  site<-rep(combined$site, times=combined$frequency)
  species_finished<-data.frame(cbind(sequences,site))
  species_finished$datatype<-paste("Sanger",sangerlength)
  species_finished$sequences <- as.character(species_finished$sequences)
  if(removefirst==TRUE){
  species_finished$sequences <- gsub("^.{0,1}", "", species_finished$sequences)
  }
  assign(  paste("sanger", speciesname, sep = ""), species_finished, envir = parent.frame())
}
#creating pie charts and haplotype maps for specific species
hap_pie<-function(eDNAdata, sangerdata, OTUnum, speciesname,folder) {
  folder_name2<-paste(folder,"/",speciesname,"2/",sep = "")
  dir.create(folder_name2)
  
  # take just the species of interest from eDNA data
  eDNAdataspecies<-eDNAdata[eDNAdata$OTU==OTUnum,]
  
  #make sanger sequences in upper case to be similar to eDNA format
  sangerdata$sequences<-toupper(sangerdata$sequences)
  
  #order haplotypes so the most abundant haplotypes are the same colors
  #get summary so we know how to order
  counte<-data.frame(table(eDNAdataspecies$sequences))
  eDNAdataspecies$counts<-counte$Freq[match(eDNAdataspecies$sequences,counte$Var1)]
  eDNAdataspecies<-eDNAdataspecies[order(eDNAdataspecies$counts, decreasing = TRUE),] 
  countsa<-data.frame(table(sangerdata$sequences))
  sangerdata$counts<-countsa$Freq[match(sangerdata$sequences,countsa$Var1)]
  sangerdata<-sangerdata[order(sangerdata$counts, decreasing = TRUE),]  
  
  #merge sanger and eDNA sequences to give unique names to sequences and add these to main datasets
  metadata<-data.frame(sequences=c(sangerdata$sequences,eDNAdataspecies$sequences))
  metadata2 <- transform(metadata, hap_id=paste("Haplotype ",match(sequences, unique(sequences)),sep = ""))
  #metadata<- metadata %>% select(sequences) %>% mutate(hap_id = paste("Haplotype",(dense_rank(sequences))))
  eDNAdataspecies$hap_id <- metadata2$hap_id[match(eDNAdataspecies$sequences,metadata2$sequences)]
  sangerdata$hap_id <- metadata2$hap_id[match(sangerdata$sequences,metadata2$sequences)]
  orderedhap<-paste("Haplotype ",1:length(unique(eDNAdataspecies$hap_id)),sep="")
  eDNAdataspecies$hap_id<-factor(eDNAdataspecies$hap_id,levels=orderedhap)
  orderedhap2<-paste("Haplotype ",1:length(unique(sangerdata$hap_id)),sep="")
  sangerdata$hap_id<-factor(sangerdata$hap_id,levels=orderedhap2)
  
  #set colors
  #make the first 10 colors the same 
  firstcol<-c("#90B9A1","#35775F","#52B9A6","#3F97CE","#C2E9EB","#F28170","#D3445B","#E8A9D1","#DCC8FE")
  n<-length(unique(metadata2$hap_id))
  n<-if(n>10) {n-9}
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual'|brewer.pal.info$category == "div",]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  haplo_col<-sample(col_vector, n)
  haplocol2<-c(firstcol,haplo_col)
  names(haplocol2)<-unique(metadata2$hap_id)
  #Pie charts for eDNA vs Sanger
  #eDNA
  tempr<-eDNAdataspecies
  tempr_counts <- ddply(eDNAdataspecies, c("hap_id"), summarise, nrows = length(hap_id))
  tempr_counts$total<-nrow(tempr)
  tempr_counts$percent<- tempr_counts$nrows / tempr_counts$total *100
  blank_theme <- theme_minimal()+ 
    theme(panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          plot.title=element_text(size=14, face="bold"),legend.background = element_rect(fill='transparent'), 
          legend.box.background = element_rect(fill='transparent') )
  
  png(file=paste(folder_name2, speciesname,"_piechart_eDNA.png", sep=""), width=10, height=10, units = "in", res=400,bg="transparent")
  piechart<- ggplot(tempr_counts, aes(x="", y=percent, fill=hap_id))+
    geom_bar(width = 1, stat = "identity")+
    scale_fill_manual(values = haplocol2)+
    coord_polar("y", start=0)+ 
    blank_theme +
    theme_void()
  print(piechart)
  dev.off()
  #Sanger
  tempr<-sangerdata
  tempr_counts <- ddply(sangerdata, c("hap_id"), summarise, nrows = length(hap_id))
  tempr_counts$total<-nrow(tempr)
  tempr_counts$percent<- tempr_counts$nrows / tempr_counts$total *100
  blank_theme <- theme_minimal()+ 
    theme(panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          plot.title=element_text(size=14, face="bold"),legend.background = element_rect(fill='transparent'), 
          legend.box.background = element_rect(fill='transparent') )
  
  png(file=paste(folder_name2, speciesname,"_piechart_sanger.png", sep=""), width=10, height=10, units = "in", res=400,bg="transparent")
  piechart<- ggplot(tempr_counts, aes(x="", y=percent, fill=hap_id))+
    geom_bar(width = 1, stat = "identity")+
    scale_fill_manual(values = haplocol2)+
    coord_polar("y", start=0)+ 
    blank_theme +
    theme_void()
  print(piechart)
  dev.off()
  
  #By Island make pie charts for eDNA data
  eIsland<-unique(eDNAdataspecies$site)
  for (i in 1:length(eIsland)) {
    tempr<-eDNAdataspecies[eDNAdataspecies$site==eIsland[i],]
    tempr_counts <- ddply(tempr, c("hap_id"), summarise, nrows = length(hap_id))
    tempr_counts$total<-nrow(tempr)
    tempr_counts$percent<- tempr_counts$nrows / tempr_counts$total *100
    
    #pie charts
    blank_theme <- theme_minimal()+ 
      theme(panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent', color=NA),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.ticks = element_blank(),
            plot.title=element_text(size=14, face="bold"),legend.background = element_rect(fill='transparent'), 
            legend.box.background = element_rect(fill='transparent') )
    
    png(file=paste(folder_name2, speciesname,"_",eIsland[i], "_piechart_eDNA_island.png", sep=""), width=10, height=10, units = "in", res=400,bg="transparent")
    piechart<- ggplot(tempr_counts, aes(x="", y=percent, fill=hap_id))+
      geom_bar(width = 1, stat = "identity")+
      scale_fill_manual(values = haplocol2)+
      coord_polar("y", start=0)+ 
      blank_theme +
      theme_void()
    print(piechart)
    dev.off()
  }
  
  #By region make pie charts for eDNA data
  Region<-unique(eDNAdataspecies$region)
  for (i in 1:length(Region)) {
    tempr<-eDNAdataspecies[eDNAdataspecies$region==Region[i],]
    tempr_counts <- ddply(tempr, c("hap_id"), summarise, nrows = length(hap_id))
    tempr_counts$total<-nrow(tempr)
    tempr_counts$percent<- tempr_counts$nrows / tempr_counts$total *100
    
    #pie charts
    blank_theme <- theme_minimal()+ 
      theme(panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent', color=NA),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.ticks = element_blank(),
            plot.title=element_text(size=14, face="bold"),legend.background = element_rect(fill='transparent'), 
            legend.box.background = element_rect(fill='transparent') )
    
    png(file=paste(folder_name2, speciesname,"_",Region[i], "_piechart_eDNA.png", sep=""), width=10, height=10, units = "in", res=400,bg="transparent")
    piechart<- ggplot(tempr_counts, aes(x="", y=percent, fill=hap_id))+
      geom_bar(width = 1, stat = "identity")+
      scale_fill_manual(values = haplocol2)+
      coord_polar("y", start=0)+ 
      blank_theme +
      theme_void()
    print(piechart)
    dev.off()
  }
  
  #By location/sample site make pie charts for eDNA data
  Locations<-unique(eDNAdataspecies$sample)
  for (i in 1:length(Locations)) {
    templ<-eDNAdataspecies[eDNAdataspecies$sample==Locations[i],]
    templ_counts <- ddply(templ, c("hap_id"), summarise, nrows = length(hap_id))
    templ_counts$total<-nrow(templ)
    templ_counts$percent<- templ_counts$nrows / templ_counts$total *100
    
    #pie charts
    blank_theme <- theme_minimal()+ 
      theme(panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent', color=NA),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.ticks = element_blank(),
            plot.title=element_text(size=14, face="bold"),legend.background = element_rect(fill='transparent'), 
            legend.box.background = element_rect(fill='transparent') )
    
    png(file=paste(folder_name2, speciesname,"_",Locations[i], "_piechart_eDNA.png", sep=""), width=10, height=10, units = "in", res=400,bg="transparent")
    piechart<- ggplot(templ_counts, aes(x="", y=percent, fill=hap_id))+
      geom_bar(width = 1, stat = "identity")+
      scale_fill_manual(values = haplocol2)+
      coord_polar("y", start=0)+ 
      blank_theme +
      theme_void()
    print(piechart)
    dev.off()
  }
  
  #By Island make pie charts for sanger data
  Islands<-unique(sangerdata$site)
  for (i in 1:length(Islands)) {
    templ<-sangerdata[sangerdata$site==Islands[i],]
    templ_counts <- ddply(templ, c("hap_id"), summarise, nrows = length(hap_id))
    templ_counts$total<-nrow(templ)
    templ_counts$percent<- templ_counts$nrows / templ_counts$total *100
    
    #pie charts
    blank_theme <- theme_minimal()+ 
      theme(panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent', color=NA),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.ticks = element_blank(),
            plot.title=element_text(size=14, face="bold"),legend.background = element_rect(fill='transparent'), 
            legend.box.background = element_rect(fill='transparent') )
    
    png(file=paste(folder_name2, speciesname,"_",Islands[i], "_piechart_sanger.png", sep=""), width=10, height=10, units = "in", res=400,bg="transparent")
    piechart<- ggplot(templ_counts, aes(x="", y=percent, fill=hap_id))+
      geom_bar(width = 1, stat = "identity")+
      scale_fill_manual(values = haplocol2)+
      coord_polar("y", start=0)+ 
      blank_theme +
      theme_void()
    print(piechart)
    dev.off()
  }
  
  #haplotype map for eDNA data
  temp <- eDNAdataspecies
  temp2<-temp[,c("sequences","region")]
  y <- as.list(temp2$sequences)
  y <- as.DNAbin(ape::as.alignment(y))
  h <- pegas::haplotype(y)
  net <- haploNet(h)
  net[,3] <- net[,3] #-1
  #abundances
  #remove unused haplos from color list
  temp3<-temp[,c("sequences","hap_id")]
  haplos<-unique(temp3$hap_id)
  haplo_col2<-haplocol2[which(names(haplocol2) %in% haplos)]
  #make haplo names for pie
  y <- as.list(temp3$sequences)
  y <- as.DNAbin(ape::as.alignment(y))
  R<-haploFreq(y, fac = c(temp3$hap_id), haplo = h)
  #order so same colors as the rest
  R2 <- R[, names(haplo_col2), drop = FALSE]
  #order R so correct colors
  png(file=paste(folder_name2, speciesname, "_hapmap_eDNA.jpg", sep=""), width=10, height=10, units = "in", res=400,bg="transparent")
  setHaploNetOptions(labels = FALSE,link.width = 5, mutations.sequence.length = 5,link.width.alt=0.0001)
  plot(net, size= attr(net, "freq"), scale.ratio = 800, cex = 0.8, pie = R2, bg=haplo_col2, fast = TRUE)
  #legend("topleft", title="Haplotype",legend=colnames(R2), col=haplo_col2, cex = 1,pch = 20)
  dev.off()
  
  #haplotype map for sanger data
  y2 <- as.list(sangerdata$sequences)
  y2 <- as.DNAbin(ape::as.alignment(y2))
  h2 <- pegas::haplotype(y2)
  net2 <- haploNet(h2)
  net2[,3] <- net2[,3] #-1
  #abundances
  #remove unused haplos from color list
  temp4<-sangerdata[,c("sequences","hap_id")]
  haplos<-unique(temp4$hap_id)
  haplo_col2<-haplocol2[which(names(haplocol2) %in% haplos)]
  #make haplo names for pie
  y2 <- as.list(temp4$sequences)
  y2 <- as.DNAbin(ape::as.alignment(y2))
  R<-haploFreq(y2, fac = c(temp4$hap_id), haplo = h2)
  #order so same colors as the rest
  R2 <- R[, names(haplo_col2), drop = FALSE]
  #order R so correct colors
  png(file=paste(folder_name2, speciesname, "_hapmap_sanger.jpg", sep=""), width=10, height=10, units = "in", res=400,bg="transparent")
  setHaploNetOptions(labels = FALSE,link.width = 5, mutations.sequence.length = 50,link.width.alt=0.0001)
  plot(net2, size= attr(net2, "freq"), scale.ratio = 10, cex = 0.8, pie = R2, bg=haplo_col2)
  #legend("topleft", title="Haplotype",legend=colnames(R2), col=haplo_col2, cex = 1,pch = 20)
  dev.off()
}
#bar plots
barplotcomp<-function(eDNAdata, sangerdata, OTUnum, speciesname,folder,eDNAproxy){
  folder_name2<-paste(folder,"/",speciesname,"_barplots","/",sep = "")
  dir.create(folder_name2)
  
  # take just the species of interest from eDNA data
  eDNAdataspecies<-eDNAdata[eDNAdata$OTU==OTUnum,]
  
  #make sanger sequences in upper case to be similar to eDNA format
  sangerdata$sequences<-toupper(sangerdata$sequences)
  
  #order haplotypes so the most abundant haplotypes are the same colors
  #get summary so we know how to order
  counte<-data.frame(table(eDNAdataspecies$sequences))
  eDNAdataspecies$counts<-counte$Freq[match(eDNAdataspecies$sequences,counte$Var1)]
  eDNAdataspecies<-eDNAdataspecies[order(eDNAdataspecies$counts, decreasing = TRUE),] 
  countsa<-data.frame(table(sangerdata$sequences))
  sangerdata$counts<-countsa$Freq[match(sangerdata$sequences,countsa$Var1)]
  sangerdata<-sangerdata[order(sangerdata$counts, decreasing = TRUE),]  
  
  #merge sanger and eDNA sequences to give unique names to sequences and add these to main datasets
  metadata<-data.frame(sequences=c(sangerdata$sequences,eDNAdataspecies$sequences))
  metadata2 <- transform(metadata, hap_id=paste("Haplotype ",match(sequences, unique(sequences)),sep = ""))
  #metadata<- metadata %>% select(sequences) %>% mutate(hap_id = paste("Haplotype",(dense_rank(sequences))))
  eDNAdataspecies$hap_id <- metadata2$hap_id[match(eDNAdataspecies$sequences,metadata2$sequences)]
  sangerdata$hap_id <- metadata2$hap_id[match(sangerdata$sequences,metadata2$sequences)]
  
  #if tissue differentiated by hilo/kona, rename
  sangerdata$site[sangerdata$site == "Kona"] <-"Hawaii"
  sangerdata$site[sangerdata$site == "Hilo"] <-"Hawaii"
  
  #only keep data when both islands were sampled
  eislands<-unique(eDNAdataspecies$site)
  sislands<-unique(sangerdata$site)
  toremove<-setdiff(eislands,sislands)
  if(length(toremove)>0){
  sangerdata<-sangerdata[sangerdata$site!=toremove,]
  eDNAdataspecies<-eDNAdataspecies[eDNAdataspecies$site!=toremove,]
  }
  #get islands for colors
  metadata<-c(eDNAdataspecies$hap_id,sangerdata$hap_id)
  #set colors
  firstcol<-c("lightblue","darkgreen","tan","darkblue","lightgreen","darkorchid4","rosybrown1","darkred","aquamarine","darkorange4")
  n<-length(unique(metadata))
  if(n<10){
    firstcol2<-firstcol[1:n]
    haplocol2<-firstcol2
  }
  if(n>10) {
    n<-n-10
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual'|brewer.pal.info$category == "div",]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    haplo_col<-sample(col_vector, n)
    haplocol2<-c(firstcol,haplo_col)
  }
  if(n==10){
    haplocol2<-firstcol
  }
  names(haplocol2)<-unique(metadata)
  coldf<-as.data.frame(haplocol2)
  coldf$hap_id<-rownames(coldf)
  rownames(coldf)<-c()
  #Pie charts for eDNA vs Sanger
  #eDNA
  tempr<-eDNAdataspecies
  tempr_counts <- ddply(eDNAdataspecies, c("hap_id"), summarise, nrows = length(hap_id))
  tempr_counts$total<-nrow(tempr)
  tempr_counts$percent<- tempr_counts$nrows / tempr_counts$total *100
  tempr_edna<-tempr_counts
  tempr_edna$Species<-speciesname
  tempr_edna$method<-"eDNA"
  #Sanger
  tempr2<-sangerdata
  tempr2_counts <- ddply(sangerdata, c("hap_id"), summarise, nrows = length(hap_id))
  tempr2_counts$total<-nrow(tempr2)
  tempr2_counts$percent<- tempr2_counts$nrows / tempr2_counts$total *100
  tempr_sanger<-tempr2_counts
  tempr_sanger$Species<-speciesname
  tempr_sanger$method<-"Tissue"
  tempr_combined<-rbind(tempr_edna,tempr_sanger)
  #reorder so that 
  tempr_combined$hap_id <- factor(tempr_combined$hap_id, levels = unique(tempr_combined$hap_id[order(tempr_combined$percent)]))
  tempr_combined<-merge(tempr_combined, coldf, by="hap_id",all = TRUE)
  #remove nas
  tempr_combined<- na.omit(tempr_combined)
  #write
  write.csv(tempr_combined,file = paste(folder,"/",speciesname,"_",eDNAproxy,"_haplofreq_byspecies.csv", sep=""))
  
  #make barplot of overall species relative haplotype frequencies
  png(file=paste(folder_name2, speciesname,"_",eDNAproxy,"_barplot_freq_byspecies.png", sep=""), width=10, height=8, units = "in", res=400,bg="transparent")
  piechart<-ggplot(tempr_combined, aes(fill=hap_id, y=percent, x=method)) + 
    geom_bar(position="fill", stat="identity", color="black")+
    ylab(label = "Relative Frequency") +
    xlab(label = "")+
    labs(fill="Haplotypes")+
    scale_fill_manual(values = haplocol2)+
    theme(panel.grid.major=element_blank(),
          panel.border = element_rect(color = "white",fill=NA),
          panel.background = element_rect(fill="grey95"),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=20), 
          axis.title.y = element_text(size=25), 
          legend.title = element_text(size=20), 
          legend.text = element_text(size=18)) 
  print(piechart)
  dev.off()
  
  #By Island barplots and data
  eIsland<-unique(eDNAdataspecies$site)
  island_tempr<-data.frame()
  for (i in 1:length(eIsland)) {
    tempr<-eDNAdataspecies[eDNAdataspecies$site==eIsland[i],]
    tempr_counts <- ddply(tempr, c("hap_id"), summarise, nrows = length(hap_id))
    tempr_counts$total<-nrow(tempr)
    tempr_counts$percent<- tempr_counts$nrows / tempr_counts$total *100
    tempr_counts$site<-eIsland[i]
    tempr_counts$method<-"eDNA"
    island_tempr<-rbind(island_tempr,tempr_counts)
  }
  Islands<-unique(sangerdata$site)
  for (i in 1:length(Islands)) {
    templ<-sangerdata[sangerdata$site==Islands[i],]
    templ_counts <- ddply(templ, c("hap_id"), summarise, nrows = length(hap_id))
    templ_counts$total<-nrow(templ)
    templ_counts$percent<- templ_counts$nrows / templ_counts$total *100
    templ_counts$site<-Islands[i]
    templ_counts$method<-"Tissue"
    island_tempr<-rbind(island_tempr,templ_counts)
  }
    #reorder
    island_tempr$hap_id <- factor(island_tempr$hap_id, levels = unique(island_tempr$hap_id[order(island_tempr$percent)]))
    island_tempr$species<-speciesname
    island_tempr<-merge(island_tempr, coldf, by="hap_id",all = TRUE)
    #remove nas
    island_tempr<- na.omit(island_tempr)
    write.csv(island_tempr,file = paste(folder,"/",speciesname,"_",eDNAproxy,"_haplofreq_byisland.csv", sep=""))
    
    #make barplot of overall species relative haplotype frequencies
    png(file=paste(folder_name2, speciesname,"_",eDNAproxy,"_barplot_freq_byisland.png", sep=""), width=10, height=8, units = "in", res=400,bg="transparent")
    piechart<-ggplot(island_tempr, aes(fill=hap_id, y=percent, x=method)) + 
      geom_bar(position="fill", stat="identity", color="black")+
      ylab(label = "Relative Frequency") +
      xlab(label = "")+
      labs(fill="Haplotypes")+
      scale_fill_manual(values = haplocol2)+
      theme(panel.grid.major=element_blank(),
            panel.border = element_rect(color = "white",fill=NA),
            panel.background = element_rect(fill="grey95"),
            axis.text.x = element_text(size=25),
            axis.text.y = element_text(size=20), 
            axis.title.y = element_text(size=25), 
            legend.title = element_text(size=20), 
            legend.text = element_text(size=18)) +
      facet_wrap(~site)
    print(piechart)
    dev.off()
}
#get shared haplotyoes
sharedhaplotable<-function(eDNAdata, sangerdata, OTUnum, speciesname,eDNAproxy){ 
  # take just the species of interest from eDNA data
  eDNAdataspecies<-eDNAdata[eDNAdata$OTU==OTUnum,]
  
  #make sanger sequences in upper case to be similar to eDNA format
  sangerdata$sequences<-toupper(sangerdata$sequences)
  
  #merge sanger and eDNA sequences to give unique names to sequences and add these to main datasets
  metadata<-data.frame(sequences=c(sangerdata$sequences,eDNAdataspecies$sequences))
  metadata<- metadata %>% select(sequences) %>% mutate(hap_id = paste("Haplotype",(dense_rank(sequences))))
  eDNAdataspecies$hap_id <- metadata$hap_id[match(eDNAdataspecies$sequences,metadata$sequences)]
  sangerdata$hap_id <- metadata$hap_id[match(sangerdata$sequences,metadata$sequences)]
  
  #get unique haplotype names from sanger and eDNA
  sangerhaplo<-unique(sangerdata$hap_id)
  eDNAhaplo<-unique(eDNAdataspecies$hap_id)
  sharedhaplo<-unlist(unique(sangerhaplo[sangerhaplo %in% eDNAhaplo]))
  
  #only keep haplos shared
  eDNAshared<-eDNAdataspecies[eDNAdataspecies$hap_id %in% sharedhaplo, ]
  sangershared<-sangerdata[sangerdata$hap_id %in% sharedhaplo, ]
  
  #for eDNA
  #convert DNA seqs to format for pairphist
  y <- as.list(eDNAshared$sequences)
  y <- as.DNAbin(ape::as.alignment(y))
  z <- as.dna(y)
  #identify populations for testing
  Island<-eDNAshared$site
  #Island<-eDNAshared$region
  #calculate phist and p-values
  pst<-pairPhiST(z, Island, nperm=500, negatives=TRUE, showprogbar=FALSE) #calculates distance in code
  #save data as table
  #How to create a prettier table
  #remove pops that it did not occur in
  i_order<-c("Kauai","Oahu","Maui","Hawaii")
  #i_order<-c("Kauai","Oahu","Maui","Hilo","Kona")
  pops<-unique(Island)
  i_order<-i_order[i_order %in% pops]
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
  write.csv(phist_reordered, file = paste(speciesname,"_eDNA_",eDNAproxy,"_phistat_onlyshared.csv", sep=""))
  
  #for sanger
  #convert DNA seqs to format for pairphist
  y <- as.list(sangershared$sequences)
  y <- as.DNAbin(ape::as.alignment(y))
  z <- as.dna(y)
  #identify populations for testing
  Island<-sangershared$site
  #calculate phist and p-values
  pst<-pairPhiST(z, Island, nperm=1000, negatives=TRUE, showprogbar=FALSE) #calculates distance in code
  #save data as table
  #How to create a prettier table
  #remove pops that it did not occur in
  i_order<-c("Kauai","Oahu","Maui","Hawaii")
  #i_order<-c("Kauai","Oahu","Maui","Hilo","Kona")
  pops<-unique(Island)
  i_order<-i_order[i_order %in% pops]
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
  write.csv(phist_reordered, file = paste(speciesname,"_sanger_","_phistat_onlyshared.csv", sep=""))
  
}


#load data
popgen_pa<-read.csv("cytb_ee025_r10_m90_nSC_popgen_pa.csv", stringsAsFactors = F)
popgen_r<-read.csv("cytb_ee025_r10_m90_nSC_popgen_r.csv", stringsAsFactors = F)
popgen_n<-read.csv("cytb_ee025_r10_m90_nSC_popgen_n.csv", stringsAsFactors = F)
popgen_n_reduced<-read.csv("cytb_ee025_rawpopgen_noSC_reduced.csv", stringsAsFactors = F)

reformat(filepath = "Abuabd_short.fasta",speciesname = "Abudefduf_abdominalis",removefirst = TRUE,sangerlength = "short")
reformat(filepath = "Abuvai_short.fasta",speciesname = "Abudefduf_vaigiensis",removefirst = TRUE,sangerlength = "short")
reformat(filepath = "Acanig_short.fasta",speciesname = "Acanthurus_nigrofuscus",removefirst = TRUE,sangerlength = "short")
reformat(filepath = "Acanthurus_nigroris_short.fasta",speciesname = "Acanthurus_nigroris",removefirst = FALSE,sangerlength = "short")
#reformat(filepath = "Cepharg_short.fasta",speciesname = "Cephalopholis_argus",removefirst = FALSE,sangerlength = "short")
#reformat(filepath = "Chae_lun_shorter.fasta",speciesname = "Chaetodon_lunulatus",removefirst = FALSE)
reformat(filepath = "Chafre_short.fasta",speciesname = "Chaetodon_fremblii",removefirst = FALSE,sangerlength = "short")
reformat(filepath = "Chamul_short.fasta",speciesname = "Chaetodon_multicinctus",removefirst = FALSE,sangerlength = "short")
reformat(filepath = "Chamil_short.fasta",speciesname = "Chaetodon_miliaris",removefirst = FALSE,sangerlength = "short")
reformat(filepath = "Chrvan_shorter.fasta",speciesname = "Chromis_vanderbilti",removefirst = FALSE,sangerlength = "short")
reformat(filepath = "Ctestr_short.fasta",speciesname = "Ctenochaetus_strigosus",removefirst = FALSE,sangerlength = "short")
reformat(filepath = "Mulflav_short.fasta",speciesname = "Mulloidichthys_flavolineatus",removefirst = FALSE,sangerlength = "short")
reformat(filepath = "Zebrasoma_flavescens_short.fasta",speciesname = "Zebrasoma_flavescens",removefirst = FALSE,sangerlength = "short")


## Make haplotype pie chart and barplots
#Abudefduf vaigiensis
hap_pie(eDNAdata = popgen_r,
        sangerdata = sangerAbudefduf_vaigiensis,
        OTUnum = "OTU_37",
        speciesname = "Abudefduf_vaigiensis",
        folder = "hapmaps_w_pie")
barplotcomp(eDNAdata = popgen_n_reduced,
            sangerdata = sangerAbudefduf_vaigiensis, 
            OTUnum = "OTU_37", 
            speciesname ="Abudefduf_vaigiensis",
            folder ="barplots",
            eDNAproxy ="raw")
sharedhaplotable(eDNAdata = popgen_n_reduced,
                 sangerdata=sangerAbudefduf_vaigiensis,
                 OTUnum="OTU_37",
                 speciesname="Abudefduf_vaigiensis",
                 eDNAproxy="raw")
#Abudefduf abdominalis
hap_pie(eDNAdata = popgen_r,
        sangerdata = sangerAbudefduf_abdominalis,
        OTUnum = "OTU_60",
        speciesname = "Abudefduf_abdominalis",
        folder = "hapmaps_w_pie")
barplotcomp(eDNAdata = popgen_n_reduced,
            sangerdata = sangerAbudefduf_abdominalis,
            OTUnum = "OTU_60",
            speciesname = "Abudefduf_abdominalis",
            folder ="barplots",
            eDNAproxy ="raw")
sharedhaplotable(eDNAdata = popgen_n_reduced,
                 sangerdata=sangerAbudefduf_abdominalis,
                 OTUnum="OTU_60",
                 speciesname="Abudefduf_abdominalis",
                 eDNAproxy="raw")
#Acanthurus nigrofuscus
hap_pie(eDNAdata = popgen_n_reduced,
        sangerdata = sangerAcanthurus_nigrofuscus,
        OTUnum = "OTU_33",
        speciesname = "Acanthurus_nigrofuscus",
        folder = "hapmaps_w_pie")
barplotcomp(eDNAdata = popgen_n_reduced,
            sangerdata = sangerAcanthurus_nigrofuscus,
            OTUnum = "OTU_33",
            speciesname = "Acanthurus_nigrofuscus",
            folder ="barplots",
            eDNAproxy ="raw")
sharedhaplotable(eDNAdata = popgen_n_reduced,
                 sangerdata=sangerAcanthurus_nigrofuscus,
                 OTUnum="OTU_33",
                 speciesname="Acanthurus_nigrofuscus",
                 eDNAproxy="raw")
temp <- popgen_r[popgen_r$OTU== "OTU_33",]
temp2<-temp[,c("sequences","region")]
y <- as.list(temp2$sequences)
y <- as.DNAbin(ape::as.alignment(y))
h <- pegas::haplotype(y)
net <- haploNet(h)
net[,3] <- net[,3] #-1
# print.default(net)
R<-haploFreq(y, fac = c(temps2$region),haplo = h)
png(file="hapmaps_w_pie_rawreads/AnigF_nopops_hapmap_eDNA2.jpg", width=10, height=10, units = "in", res=400,bg="transparent")
setHaploNetOptions(labels = FALSE,haplotype.inner.color = "lightgrey",link.width = 3, mutations.sequence.length = 50,link.width.alt=0.0001)
plot(net, size= attr(net, "freq"),scale.ratio = 5, fast = TRUE) #scale.ratio = 2,
dev.off()

y2 <- as.list(sangerAcanthurus_nigrofuscus$sequences)
y2 <- as.DNAbin(ape::as.alignment(y2))
h2 <- pegas::haplotype(y2)
net2 <- haploNet(h2)
net2[,3] <- net2[,3] #-1
# print.default(net)
R2<-haploFreq(y2, fac = c(sangerAcanthurus_nigrofuscus$site), haplo = h2)
png(file="hapmaps_w_pie_rawreads/AnigF_nopops_hapmap_sanger2.jpg", width=10, height=10, units = "in", res=400,bg="transparent")
setHaploNetOptions(labels = FALSE,haplotype.inner.color = "lightgrey",link.width = 2, mutations.sequence.length = 50,link.width.alt=0.0001)
plot(net2, size= attr(net2, "freq"), scale.ratio = 5, cex = 0.8, fast = TRUE)
dev.off()
#Acanthurus nigroris
hap_pie(eDNAdata = popgen_r,
        sangerdata = sangerAcanthurus_nigroris,
        OTUnum = "OTU_80",
        speciesname = "Acanthurus_nigroris",
        folder = "hapmaps_w_pie")
barplotcomp(eDNAdata = popgen_n_reduced,
            sangerdata = sangerAcanthurus_nigroris,
            OTUnum = "OTU_80",
            speciesname = "Acanthurus_nigroris",
            folder ="barplots",
            eDNAproxy ="raw")
sharedhaplotable(eDNAdata = popgen_n_reduced,
                 sangerdata=sangerAcanthurus_nigroris,
                 OTUnum="OTU_80",
                 speciesname="Acanthurus_nigroris",
                 eDNAproxy="raw")
#Cephalopholis argus
hap_pie(eDNAdata = popgen_r,
        sangerdata = sangerCephalopholis_argus,
        OTUnum = "OTU_554",
        speciesname = "Cephalopholis_argus",
        folder = "hapmaps_w_pie")
#Chaetodon fremblii
hap_pie(eDNAdata = popgen_r,
        sangerdata = sangerChaetodon_fremblii,
        OTUnum = "OTU_1173",
        speciesname = "Chaetodon_fremblii",
        folder = "hapmaps_w_pie")
barplotcomp(eDNAdata = popgen_n,
            sangerdata = sangerChaetodon_fremblii,
            OTUnum = "OTU_1173",
            speciesname = "Chaetodon_fremblii",
            folder ="barplots",
            eDNAproxy ="raw")
#Chaetodon multicinctus
hap_pie(eDNAdata = popgen_r,
        sangerdata = sangerChaetodon_multicinctus,
        OTUnum = "OTU_246",
        speciesname = "Chaetodon_multicinctus",
        folder = "hapmaps_w_pie")
barplotcomp(eDNAdata = popgen_n,
            sangerdata = sangerChaetodon_multicinctus,
            OTUnum = "OTU_246",
            speciesname = "Chaetodon_multicinctus",
            folder ="barplots",
            eDNAproxy ="raw")
sharedhaplotable(eDNAdata = popgen_n_reduced,
                 sangerdata=sangerChaetodon_multicinctus,
                 OTUnum="OTU_246",
                 speciesname="Chaetodon_multicinctus",
                 eDNAproxy="raw")
#Chaetodon miliaris
hap_pie(eDNAdata = popgen_n,
        sangerdata = sangerChaetodon_miliaris,
        OTUnum = "OTU_568",
        speciesname = "Chaetodon_miliaris",
        folder = "hapmaps_w_pie")
barplotcomp(eDNAdata = popgen_n,
            sangerdata = sangerChaetodon_miliaris,
            OTUnum = "OTU_568",
            speciesname = "Chaetodon_miliaris",
            folder ="barplots",
            eDNAproxy ="raw")
sharedhaplotable(eDNAdata = popgen_n_reduced,
                 sangerdata=sangerChaetodon_miliaris,
                 OTUnum="OTU_568",
                 speciesname="Chaetodon_miliaris",
                 eDNAproxy="raw")
temp <- popgen_r[popgen_r$OTU== "OTU_568",]
temp2<-temp[,c("sequences","region")]
y <- as.list(temp2$sequences)
y <- as.DNAbin(ape::as.alignment(y))
h <- pegas::haplotype(y)
net <- haploNet(h)
net[,3] <- net[,3] #-1
# print.default(net)
png(file="hapmaps_w_pie_rawreads/Cmili_nopops_hapmap_eDNA.jpg", width=10, height=10, units = "in", res=400,bg="transparent")
setHaploNetOptions(labels = FALSE,haplotype.inner.color = "lightgrey",link.width = 3, mutations.sequence.length = 50,link.width.alt=0.0001)
plot(net, size= attr(net, "freq"),scale.ratio = 5, fast = TRUE) #scale.ratio = 2,
dev.off()

y2 <- as.list(sangerChaetodon_miliaris$sequences)
y2 <- as.DNAbin(ape::as.alignment(y2))
h2 <- pegas::haplotype(y2)
net2 <- haploNet(h2)
net2[,3] <- net2[,3] #-1
# print.default(net)
R2<-haploFreq(y2, fac = c(sangerChaetodon_miliaris$site), haplo = h2)
png(file="hapmaps_w_pie_rawreads/Cmili_nopops_hapmap_sanger2.jpg", width=10, height=10, units = "in", res=400,bg="transparent")
setHaploNetOptions(labels = FALSE,haplotype.inner.color = "lightgrey",link.width = 2, mutations.sequence.length = 50,link.width.alt=0.0001)
plot(net2, size= attr(net2, "freq"), scale.ratio = 5, cex = 0.8, fast = TRUE)
dev.off()

#Chromis vanderbilti
hap_pie(eDNAdata = popgen_r,
        sangerdata = sangerChromis_vanderbilti,
        OTUnum = "OTU_68",
        speciesname = "Chromis_vanderbilti",
        folder = "hapmaps_w_pie")
barplotcomp(eDNAdata = popgen_n,
            sangerdata = sangerChromis_vanderbilti,
            OTUnum = "OTU_68",
            speciesname = "Chromis_vanderbilti",
            folder ="barplots",
            eDNAproxy ="raw")
sharedhaplotable(eDNAdata = popgen_n_reduced,
                 sangerdata=sangerChromis_vanderbilti,
                 OTUnum="OTU_68",
                 speciesname="Chromis_vanderbilti",
                 eDNAproxy="raw")
#Ctenochaetus strigosus
hap_pie(eDNAdata = popgen_n,
        sangerdata = sangerCtenochaetus_strigosus,
        OTUnum = "OTU_1089",
        speciesname = "Ctenochaetus_strigosus",
        folder = "hapmaps_w_pie_rawreads")
barplotcomp(eDNAdata = popgen_n,
            sangerdata = sangerCtenochaetus_strigosus,
            OTUnum = "OTU_1089",
            speciesname = "Ctenochaetus_strigosus",
            folder ="barplots",
            eDNAproxy ="raw")
#Mulloidichthys flavolineatus
hap_pie(eDNAdata = popgen_r,
        sangerdata = sangerMulloidichthys_flavolineatus,
        OTUnum = "OTU_149",
        speciesname = "Mulloidichthys_flavolineatus",
        folder = "hapmaps_w_pie")
barplotcomp(eDNAdata = popgen_n,
            sangerdata = sangerMulloidichthys_flavolineatus,
            OTUnum = "OTU_149",
            speciesname = "Mulloidichthys_flavolineatus",
            folder ="barplots",
            eDNAproxy ="raw")
sharedhaplotable(eDNAdata = popgen_n_reduced,
                 sangerdata=sangerMulloidichthys_flavolineatus,
                 OTUnum="OTU_149",
                 speciesname="Mulloidichthys_flavolineatus",
                 eDNAproxy="raw")
#Zebrasoma flavescens
hap_pie(eDNAdata = popgen_n,
        sangerdata = sangerZebrasoma_flavescens,
        OTUnum = "OTU_6",
        speciesname = "Zebrasoma_flavescens",
        folder = "hapmaps_w_pie_rawreads")
barplotcomp(eDNAdata = popgen_n,
            sangerdata = sangerZebrasoma_flavescens,
            OTUnum = "OTU_6",
            speciesname = "Zebrasoma_flavescens",
            folder ="barplots",
            eDNAproxy ="raw")
sharedhaplotable(eDNAdata = popgen_n_reduced,
                 sangerdata=sangerZebrasoma_flavescens,
                 OTUnum="OTU_6",
                 speciesname="Zebrasoma_flavescens",
                 eDNAproxy="raw")
#all species together
files<-list.files("barplots/tables_raw_redo")
setwd("/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/Hapmaps/barplots/tables_raw_redo")
allspecies<-do.call(rbind,lapply(files,read.csv))

#rename species
allspecies$Species<-sub("_"," ",allspecies$Species)

#get shared and not haplotypes
shared <- ddply(allspecies, c("Species","hap_id"), summarise, nrows = length(hap_id))
allspecies2<-merge(allspecies,shared,by=c("Species","hap_id"),all = TRUE)

#color all not shared haplotypes the same
#get count of not shared haplotypes per species
uniquehaplo <- ddply(allspecies2, c("Species","method","nrows.y"), summarise, nrows = length(nrows.y))
colnames(uniquehaplo)<-c("Species","method","unique1 or shared2","haplotype counts")
write.csv(uniquehaplo,"countsofuniquehaplotypes.csv")
#pull out the shared haplotypes
allspeciesshared<-allspecies2[allspecies2$nrows.y==2,]
#pull out all unique haplotypes
allspeciesunique<-allspecies2[allspecies2$nrows.y==1,]
#sum uniques by species and by methods
allspeciesunique2<-allspeciesunique %>%
  group_by(Species,method) %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
allspeciesunique2$hap_id<-"not shared"
allspeciesuniquefinal<-data.frame("Species"=allspeciesunique2$Species,
                                  "method"=allspeciesunique2$method,
                                  "hap_id"=allspeciesunique2$hap_id,
                                  "percent"=allspeciesunique2$percent)
allspeciessharedfinal<-data.frame("Species"=allspeciesshared$Species,
                                  "method"=allspeciesshared$method,
                                  "hap_id"=allspeciesshared$hap_id,
                                  "percent"=allspeciesshared$percent)
allspecies3<-rbind(allspeciessharedfinal,allspeciesuniquefinal)
write.csv(allspecies3,"needtorenamehapid.csv")
allspeciesfixed<-read.csv("needtorenamehapid.csv")


#make sure order of x axis is what I want
allspeciesfixed$method <- factor(allspeciesfixed$method, levels=c("Tissue", "eDNA"))
orderedhap<-paste("Haplotype ",1:9,sep="")
orderedhap<-c(orderedhap,"not shared")
allspeciesfixed$hap_id<-factor(allspeciesfixed$hap_id,levels=orderedhap)

#get colors the same
haplocol_redo <- c("lightsteelblue1","lightcyan2","cadetblue3","cadetblue4","darkslategray","seagreen","seagreen2","palegreen","aquamarine","turquoise1")
haplocol_redo<-c("#90B9A1","#35775F","#52B9A6","#3F97CE","#C2E9EB","#F28170","#D3445B","#E8A9D1","#DCC8FE","grey")
library(viridis)
colors_redo2<-viridis(n = 10)
#library(RColorBrewer)
# Get 10 colors from the YlGnBu palette
#ylgnbu_colors <- rev(brewer.pal(9, "YlGnBu"))
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual'|brewer.pal.info$category == "div",]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#haplo_col<-sample(col_vector, 69)
haplocol_redo<-c(colors_redo2[2:10],"aquamarine")
namcol<-c(paste("Haplotype ",1:9,sep=""),"not shared")
names(haplocol_redo)<-namcol

#make all species barplot
png(file="all_species_barplot_eDNA_raw_redo4.png", width=14.75, height=13, units = "in", res=400,bg="transparent")
piechart<-ggplot(allspeciesfixed, aes(fill=hap_id, y=percent, x=method)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity", color=NA)+
  ylab(label = "Relative Frequency") +
  xlab(label = "Method")+
  coord_cartesian(ylim = c(0, 1.09), expand = FALSE)+
  scale_fill_manual(values = haplocol_redo)+
  facet_wrap(~Species)+
  theme(strip.text = element_text(face = "italic",size = 16,margin = margin(4.4,4.4,40,4.4,"pt")),
        strip.background=element_rect(colour="black", fill="white"))+
  theme(panel.grid.major=element_blank(),
        panel.border = element_rect(color = "white",fill=NA),
        panel.background = element_rect(fill="grey95"),
        panel.spacing.x = unit(2,"lines"),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=18), 
        axis.title.y = element_text(size=25), 
        axis.title.x = element_text(size=25),
        legend.position="right") 
print(piechart)
dev.off()



