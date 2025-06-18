library("ggplot2") # load library
library("rfishbase")

setwd("/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/ComparativePopgen") #set working directory


################################ get shared data ################################
#load
#phist<-read.csv("cytb_ee025_m90_nSC_site_r_n20_phist_table_manyspecies.csv", stringsAsFactors = F) #load in phist values
reads<-read.csv("cytb_ee025_r10_m90_nSC_summmary_rarefied_haplotypes.csv", stringsAsFactors = F) #load in read count and sample numbers
taxonomy<-read.csv("Taxonomy.csv", stringsAsFactors = F) #load in taxonomy
#samples<-read.csv("cytb_ee025_site_sample_stats.csv", stringsAsFactors = F) #load in the number of samples each OTU was detected in 
geneticdiv<-read.csv("cytb_pa10_m90_nSC_ee25_site_r_basic_seq_stats.csv", stringsAsFactors = F)

#subset relevant columns from each df
#phist2<-phist[c(4,9:13,16)]
#reads2<-reads[c(2,7:19)]
#samples2<-samples[c(4:8)]
geneticdiv2<-geneticdiv[c(4:7,10:13,15:16)]
geneticdiv2<-geneticdiv2[geneticdiv2$Population=="overall",]

#match data into one df
#colnames(phist2)<-c("OTU","Population","Population2","Pairwise Phist","p-value","fdr p-value","significant with fdr") #rename so can merge
#combined_df<-merge(phist2,samples2, by=c("OTU","Population"),all.x=TRUE) #merge
#colnames(combined_df)<-c("OTU","Population 1","Population","Pairwise Phist","p-value","fdr p-value","significant with fdr","Pop 1 Read Depth","Pop 1 Samples Detected In","Pop 1 Number of Unique Haplos and Samples") #rename so can merge
#combined_df2<-merge(combined_df,samples2, by=c("OTU","Population"),all.x=TRUE) #merge
#colnames(combined_df2)<-c("OTU","Population 1","Population","Pairwise Phist","p-value","fdr p-value","significant with fdr","Pop 1 Read Depth","Pop 1 Samples Detected In","Pop 1 Number of Unique Haplos and Samples","Pop 2 Read Depth","Pop 2 Samples Detected In","Pop 2 Number of Unique Haplos and Samples") #rename so can merge
#combined_final<-merge(combined_df2,taxonomy, by="OTU",all.x=TRUE) #merge
combined_final<-merge(geneticdiv2,taxonomy, by="OTU",all.x=TRUE) #merge

#only marine species
combined_marine<-combined_final[combined_final$marine.freshwater != "freshwater",]
mspecieslist<-unique(combined_marine$OTU) #get list of marine species

#remove low sampling
combined_marine<-combined_marine[combined_marine$Sample.Size > 19,]

#remove non species level matches
combined_marine<-combined_marine[combined_marine$species.level.match.=="yes",]

#remove NAs
combined_marine_final<-na.omit(combined_marine)

#add info from fishbase
num_marinetaxa<-length(unique(combined_marine_final$OTU)) #get number of taxa
marspecies<-unique(combined_marine_final$Species)
larvalinfo<-larvae(species_list = marspecies)
larvalinfoIcareabout<-larvalinfo[c(1,9,36,40)]
speciesinfo<-species(species_list = marspecies)
speciesinfoIcareabout<-speciesinfo[c(1,15,23,28,29,40,41,58,76,79,81,92)]
combined_marine_lifehistory<-merge(combined_marine,speciesinfoIcareabout, by="Species",all.x=TRUE) #merge
combined_marine_lifehistory<-merge(combined_marine_lifehistory,larvalinfoIcareabout, by="Species",all.x=TRUE) #merge


#need to write this file so I can manually fill in the missing info from fishbase
write.csv(combined_marine_lifehistory,file = "cytb_r_ee025_geneticdiv_lifehistory_comp_N20.csv")





################################ paired t-test and boxplot per island ################################
library("PairedData")
Taxonomy<-read.csv("Taxonomy.csv", stringsAsFactors = F) #load in taxonomy
geneticdiv<-read.csv("cytb_pa10_m90_nSC_ee25_site_r_basic_seq_stats.csv", stringsAsFactors = F) #load in genetic diversity
popdiv<-read.csv("cytb_r10_m90_nSC_ee025_site_r_n10_rep1000_subsetmaui_overallphist.csv", stringsAsFactors = F) #load in pop gen stats
#only keep relevant columns
popdiv<-popdiv[c(4:5,8:11)] 
geneticdiv2<-geneticdiv[c(4:7,10:13,15:16)] 
# only keep by island comparisons
geneticdiv2<-geneticdiv2[geneticdiv2$Population!="overall",] 
popdiv<-popdiv[popdiv$Population!="overall",] 
# merge
combined_final<-merge(geneticdiv2,taxonomy, by="OTU",all.x=TRUE) 
combined_final2<-merge(combined_final,popdiv, by=c("OTU","Population"),all.x=TRUE) 
# only keep marine
unique(combined_final2$marine.freshwater) # check what options are to avoid mistakes like spaces
combined_marine<-combined_final2[combined_final2$marine.freshwater != "freshwater",] # only marine
#only keep species level matches
#unique(combined_marine$species.level.match.) # check what options are to avoid mistakes like spaces
#combined_mspecies<-combined_marine[combined_marine$species.level.match.!="no",] # only species level matches
#num_marinetaxa<-length(unique(combined_mspecies$OTU)) #get number of taxa
#remove low sample sizes
combined_marine<-combined_marine[combined_marine$Sample.Size > 19,]
#remove where only 1 haplotype was detected
combined_marine<-combined_marine[combined_marine$Number.of.Haplotypes > 1,]
#remove NAs
combined_mspecies<-na.omit(combined_mspecies)
#only shared species
all4<-table(combined_mspecies$OTU)
listall4<-names(all4[all4 == 4])
combined_mspecies_all<-combined_mspecies[combined_mspecies$OTU %in% listall4,]
length(unique(combined_mspecies_all$OTU))
#remove nas
combined_mspecies_all<-na.omit(combined_mspecies_all)
#change negative fst and phist to zeros
combined_mspecies_all$WC84.Fst[combined_mspecies_all$WC84.Fst<0] <- 0
combined_mspecies_all$phist[combined_mspecies_all$phist<0] <- 0
combined_mspecies$WC84.Fst[combined_mspecies$WC84.Fst<0] <- 0
combined_mspecies$phist[combined_mspecies$phist<0] <- 0

#function for paired t-test of 2 islands
statpairedttest<-function(datafile,genstat,island1,island2){
#only keep species that were detected in both islands
  islands2<-datafile[datafile$Population == island1 | datafile$Population == island2, ]
  tokeep<-names(table(islands2$OTU))[table(islands2$OTU) == 2]
  islands2<-islands2[islands2$OTU %in% tokeep,]
  # Subset
  is1 <- subset(islands2,  Population == island1, genstat, drop = TRUE)
  is2 <- subset(islands2,  Population == island2, genstat, drop = TRUE)
  jpeg(file=paste(genstat,"_",island1,"_",island2,"_pairedttestresults.jpg",sep = ""), width=8, height=6, units="in",res=600)
  pd <- paired(is1, is2)
  resfig<-plot(pd, type = "profile") + theme_bw()
  print(resfig)
  dev.off()
  # compute the difference
  d <- islands2[[genstat]][islands2$Population == island1] - islands2[[genstat]][islands2$Population == island2]
  #test normality
  normtest<-shapiro.test(d) 
  #paired t-test
  tresult<-t.test(is1, is2, paired = TRUE)
  #df to print
  results<-data.frame(normalityp=normtest$p.value[1],tstat=tresult$statistic[1],tdf=tresult$parameter[1],tpval=tresult$p.value[1],meandiff=tresult$estimate[1])
  write.csv(results,file = paste(genstat,"_",island1,"_",island2,"_pairedttestresults.csv",sep = ""))
}

statpairedttest(datafile=combined_mspecies,
                genstat = "R2",
                island1 = "Oahu",
                island2 = "Maui")

combined_mspecies_all$Population <- factor(combined_mspecies_all$Population, levels = c("Kauai","Oahu","Maui","Hawaii"))
jpeg("R2_Shared_r_n20_marine_islands.jpeg", width=8, height=6, units="in",res=600)
ggplot(combined_mspecies_all, aes(x = Population, y =  R2, fill = Population)) + 
  geom_boxplot()+
  ylab(label = "R2") +
  xlab(label = "Island") +
  #ylim(0,0.18)+ 
  scale_fill_manual(values=c("tan","lightblue","lightgreen","darkgrey"))+
  theme (panel.grid.major=element_blank(),
         panel.border = element_rect(color = "white",fill=NA),
         panel.background = element_rect(fill="grey95"),
         axis.text.x = element_text(size=18),
         axis.text.y = element_text(size=18), 
         axis.title.x = element_text(size=22),
         axis.title.y = element_text(size=22)) 
dev.off()




################################ paired t-test and boxplot per taxon ################################
Taxonomy<-read.csv("Taxonomy.csv", stringsAsFactors = F) #load in taxonomy
geneticdiv<-read.csv("cytb_pa10_m90_nSC_ee25_site_r_basic_seq_stats.csv", stringsAsFactors = F) #load in genetic diversity
popdiv<-read.csv("cytb_r10_m90_nSC_ee025_site_r_n20_rep1000_overallphist.csv", stringsAsFactors = F) #load in pop gen stats
#only keep relevant columns
popdiv<-popdiv[c(4:5,8:11)] 
geneticdiv2<-geneticdiv[c(4:7,10:13,15:16)] 
# only keep by species comparisons
geneticdiv2<-geneticdiv2[geneticdiv2$Population=="overall",] 
popdiv<-popdiv[popdiv$Population=="overall",] 
# merge
combined_final<-merge(geneticdiv2,taxonomy, by="OTU",all.x=TRUE) 
combined_final2<-merge(combined_final,popdiv, by="OTU",all.x=TRUE) 
# only keep marine
unique(combined_final2$marine.freshwater) # check what options are to avoid mistakes like spaces
combined_marine<-combined_final2[combined_final2$marine.freshwater != "freshwater",] # only marine
#only keep species level matches
unique(combined_marine$species.level.match.) # check what options are to avoid mistakes like spaces
combined_mspecies<-combined_marine[combined_marine$species.level.match.!="no",] # only species level matches
num_marinetaxa<-length(unique(combined_mspecies$OTU)) #get number of taxa
#remove low sample sizes
combined_mspecies<-combined_mspecies[combined_mspecies$Sample.Size > 19,]
#remove where only 1 haplotype was detected
combined_mspecies<-combined_mspecies[combined_mspecies$Number.of.Haplotypes > 1,]
#remove nas
combined_mspecies<-na.omit(combined_mspecies)
#change negative fst and phist to zeros
combined_mspecies$WC84.Fst[combined_mspecies$WC84.Fst<0] <- 0
combined_mspecies$phist[combined_mspecies$phist<0] <- 0

combined_mspecies <- combined_mspecies[order(combined_mspecies$class),]
combined_mspecies$order <- factor(combined_mspecies$order, levels = unique(combined_mspecies$order))

jpeg("haplodiv_r_n20_byorder.jpeg", width=12, height=8, units="in",res=600)
ggplot(combined_mspecies, aes(x = order, y = haplotype.diversity, fill = class)) + 
  geom_boxplot()+
  ylab(label = "Haplotype Diversity") +
  xlab(label = "Order") +
  labs(fill="Class") +
  scale_fill_manual(values=c("turquoise","darkblue","darkgreen"))+
  theme (panel.grid.major=element_blank(),
         panel.border = element_rect(color = "white",fill=NA),
         panel.background = element_rect(fill="grey95"),
         axis.text.x = element_text(size=18),
         axis.text.y = element_text(size=18), 
         axis.title.x = element_text(size=22),
         axis.title.y = element_text(size=22), 
         legend.title = element_text(size=18), 
         legend.text = element_text(size=18)) +
  coord_flip()
dev.off()
jpeg("nucleodiv_r_n20_byorder.jpeg", width=12, height=8, units="in",res=600)
ggplot(combined_mspecies, aes(x = order, y = nucleotide.diversity, fill = class)) + 
  geom_boxplot() +
  ylab(label = "Nucleotide Diversity") +
  xlab(label = "Order") +
  labs(fill="Class") +
  scale_fill_manual(values=c("turquoise","darkblue","darkgreen"))+
  theme (panel.grid.major=element_blank(),
         panel.border = element_rect(color = "white",fill=NA),
         panel.background = element_rect(fill="grey95"),
         axis.text.x = element_text(size=18),
         axis.text.y = element_text(size=18), 
         axis.title.x = element_text(size=22),
         axis.title.y = element_text(size=22), 
         legend.title = element_text(size=18), 
         legend.text = element_text(size=18)) +
  coord_flip()
dev.off()
jpeg("fst_r_n20_byorder.jpeg", width=12, height=8, units="in",res=600)
ggplot(combined_mspecies, aes(x = order, y = WC84.Fst, fill = class)) + 
  geom_boxplot() +
  ylab(label = "Fst") +
  xlab(label = "Order") +
  labs(fill="Class") +
  scale_fill_manual(values=c("turquoise","darkblue","darkgreen"))+
  theme (panel.grid.major=element_blank(),
         panel.border = element_rect(color = "white",fill=NA),
         panel.background = element_rect(fill="grey95"),
         axis.text.x = element_text(size=18),
         axis.text.y = element_text(size=18), 
         axis.title.x = element_text(size=22),
         axis.title.y = element_text(size=22), 
         legend.title = element_text(size=18), 
         legend.text = element_text(size=18)) +
  coord_flip()
dev.off()
jpeg("phist_r_n20_byorder.jpeg", width=12, height=8, units="in",res=600)
ggplot(combined_mspecies, aes(x = order, y = phist, fill = class)) + 
  geom_boxplot() +
  ylab(label = "Fst") +
  xlab(label = "Order") +
  labs(fill="Class") +
  scale_fill_manual(values=c("turquoise","darkblue","darkgreen"))+
  theme (panel.grid.major=element_blank(),
         panel.border = element_rect(color = "white",fill=NA),
         panel.background = element_rect(fill="grey95"),
         axis.text.x = element_text(size=18),
         axis.text.y = element_text(size=18), 
         axis.title.x = element_text(size=22),
         axis.title.y = element_text(size=22), 
         legend.title = element_text(size=18), 
         legend.text = element_text(size=18)) +
  coord_flip()
dev.off()

#remove low sampled orders
toremove<-names(table(combined_mspecies$order))[table(combined_mspecies$order) < 5]
for (i in 1:length(toremove)) {
  combined_mspecies<-combined_mspecies[combined_mspecies$order != toremove[i],]
}


#test normality
shapiro.test(combined_mspecies$haplotype.diversity)
shapiro.test(combined_mspecies$nucleotide.diversity)#not normal
shapiro.test(combined_mspecies$phist)#not normal
shapiro.test(combined_mspecies$WC84.Fst)#not normal

#kruskall-wallis test
kruskal.test(phist ~ family, data = combined_mspecies)




################################ figure out the ave num haplotypes ################################
mhaplo<-mean(reads$Total.Haplotypes)
mhaplo
maxhaplo<-max(reads$Total.Haplotypes)
maxhaplo
minhaplo<-min(reads$Total.Haplotypes)
minhaplo
sdhaplo<-sd(reads$Total.Haplotypes)
sdhaplo
sehaplo<-sdhaplo/sqrt(length(reads$OTU))
sehaplo
# figure out range of taxonomy
unique(Taxonomy$class)




################################ NMDS for haplotypes ################################
library("vegan")
library(tidyr)
library(dplyr)
library(pairwiseAdonis)
#load data
otudf<-read.csv("cytb_ee025_r10_m90_nSCspecies_haplotype_table_filtered_noSC.csv", stringsAsFactors = F)
metadata<-read.csv("Metadata_MHI.csv",stringsAsFactors = F)
#only keep marine species
motudf<-otudf[otudf$OTU %in% mspecieslist,]
#make matrix
otumatrix<-as.matrix(motudf[c(6:261)])
otumatrix<-t(otumatrix)
#sum replicates into sites
site<- metadata$Location[match(row.names(otumatrix),metadata$Samples)]
otumatrix2<-as.data.frame(cbind(site,otumatrix))
otumatrix2[2:2508] <- sapply(otumatrix2[2:2508], as.numeric)
row.names(otumatrix2)<-NULL
x<-otumatrix2
x <- aggregate(x[,2:ncol(x)], by = x[1], function(x) sum(x, na.rm = TRUE))

#only include species that are found on all islands
ix<-x
ix$site<- metadata$Island[match(x$site,metadata$Location)]
ix<-aggregate(ix[,2:ncol(ix)], by = ix[1], function(ix) sum(ix, na.rm = TRUE))
ix2<-as.data.frame(t(ix[,-1]))
colnames(ix2)<-ix$site
ix2$OTU<-motudf$OTU
ix3<-aggregate(ix2[,1:4], by = ix2[5], function(ix2) sum(ix2, na.rm = TRUE))
ix4<-ix3[ix3$Hawaii>0 & ix3$Kauai>0 & ix3$Maui>0 & ix3$Oahu>0,]
otustokeep<-ix4$OTU
motudfallislands<-motudf[motudf$OTU %in% otustokeep,]
motumatrixall<-as.matrix(motudfallislands[c(6:261)])
motumatrixall<-t(motumatrixall)
site<- metadata$Location[match(row.names(motumatrixall),metadata$Samples)]
otumatrix2<-as.data.frame(cbind(site,motumatrixall))
otumatrix2[2:2095] <- sapply(otumatrix2[2:2095], as.numeric)
row.names(otumatrix2)<-NULL
x<-otumatrix2
x <- aggregate(x[,2:ncol(x)], by = x[1], function(x) sum(x, na.rm = TRUE))

#sum haplotypes to OTUs
motudf2<-motudf[c(5:261)]
motudf2[2:257] <- sapply(motudf2[2:257], as.numeric)
row.names(motudf2)<-NULL
x<-motudf2
x <- aggregate(x[,2:ncol(x)], by = x[1], function(x) sum(x, na.rm = TRUE))
x2<-as.matrix(x[c(2:257)])
x2<-t(x2)

#eDNA index
indexmatrix<-wisconsin(otumatrix)

#run nmds
set.seed(1235)
nmds <- vegan::metaMDS(x[2:2508], distance = "bray", binary=TRUE)
#permanova
islands<- metadata$Island[match(x$site,metadata$Location)]
region<-metadata$Side[match(x$site,metadata$Location)]
site<-metadata$Location[match(row.names(otumatrix),metadata$Samples)]

#because betadisper was significant
ano = anosim(x[2:2508], islands, distance = "bray", permutations = 9999)
ano

permisland<-adonis2(x[2:2508]~islands+region, method = "bray",binary=TRUE)
permisland
islands<-metadata$Island[match(x$site,metadata$Location)]
site<-x$site
permisland<-adonis2(x[2:2508]~islands, method = "bray")
permisland
x2<-vegdist(x[2:2508], method="bray")
dispperm<-betadisper(x2, islands)
anova(dispperm)
adonis2(dist(dispperm$distances) ~ islands)
pairwise.adonis(x[2:2095],islands,sim.method = "jaccard")

#dataframe for the plot
data.scores = as.data.frame(scores(nmds)$sites)
data.scores$Sample <- x$site
#data.scores$Sample = row.names(otumatrix)
data.scores$Island <- metadata$Island[match(data.scores$Sample,metadata$Location)]
#data.scores$Location <- metadata$Location[match(data.scores$Sample,metadata$Samples)]
data.scores$Region <- metadata$Side[match(data.scores$Sample,metadata$Location)]
head(data.scores)
#plot
jpeg("nmds_bray_allOTU_haploinput_region_reads_marine.jpeg", width=12, height=8, units="in",res=600)
xx <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2, colour = Region)) + 
  geom_point(size = 4) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Island", y = "NMDS2")  + 
  stat_ellipse(size=1)
#+scale_colour_manual(values = c("#009E73", "#E69F00","purple","darkblue")) 
print(xx)
dev.off()





################################ NMDS for paireise phist ################################
phist<-read.csv("cytb_ee025_m90_nSC_sideofisland_withMauisub_n20_r_phist_table_manyspecies.csv", stringsAsFactors = F) #load in phist values
phist2<-phist[4:7]
phist2<-na.omit(phist2)
#remove south oahu, east kauai
phist2<-phist2[phist2$Population.1 != "East Kauai",]
phist2<-phist2[phist2$Population.2 != "East Kauai",]
phist2<-phist2[phist2$Population.1 != "South Oahu",]
phist2<-phist2[phist2$Population.2 != "South Oahu",]
phist2<-phist2[phist2$Population.1 != "Olowalu",]
phist2<-phist2[phist2$Population.2 != "Olowalu",]
phist2<-phist2[phist2$Population.1 != "South Makena",]
phist2<-phist2[phist2$Population.2 != "South Makena",]
#remove replicate rows
phist3<-phist2[!duplicated(phist2), ]
#only keep when all islands are present
otucount<-as.data.frame(table(phist3$OTU))
keepotu<-otucount$Var1[otucount$Freq >35]
phistall<-phist2[phist2$OTU %in% keepotu,]
phistall$PhiST.value[phistall$PhiST.value < 0] <- 0


#boxplot for pairwise phist
phist_island<-read.csv("cytb_ee025_m90_nSC_site_r_n20_phist_table_manyspecies_onlyallislands.csv", stringsAsFactors = F) #load in phist values
phistpairbox<-data.frame(populationcomp = paste(phist_island$Population.1," to ", phist_island$Population.2, sep = ""),phistval = phist_island$PhiST.value)
#phistpairbox$populationcomp <- factor(phistpairbox$populationcomp, levels = c("","","",""))
jpeg("pairwisephist_boxplot_shared_sideofisland_n20_r.jpeg", width=8, height=6, units="in",res=600)
ggplot(phistpairbox, aes(x = populationcomp, y =  phistval)) + 
  geom_boxplot()+
  ylab(label = "Pairwise ΦST") +
  xlab(label = "Island Comparisons") +
  #ylim(0,0.18)+ 
  #scale_fill_manual(values=c("tan","lightblue","lightgreen","darkgrey"))+
  theme (panel.grid.major=element_blank(),
         panel.border = element_rect(color = "white",fill=NA),
         panel.background = element_rect(fill="grey95"),
         axis.text.x = element_text(size=18,angle=90, hjust=1),
         axis.text.y = element_text(size=18), 
         axis.title.x = element_text(size=22),
         axis.title.y = element_text(size=22)) 
dev.off()

phistave<-aggregate(phistall$PhiST.value, by=list(phistall$Population.1, phistall$Population.2), mean)
phistave_table<-phistave %>% pivot_wider(names_from = Group.2, values_from = x)
write.csv(phistave_table,"phistave_maui_n20.csv")
phistave_table2<-read.csv("phistave.csv", stringsAsFactors = F,header = TRUE)
row.names(phistave_table2)<-phistave_table2$X
phistave_table2<-phistave_table2[2:9]
phistave_table3<-t(phistave_table2)
phistdist2<-as.dist(phistave_table3)
set.seed(689)
nmds <- vegan::metaMDS(phistdist2,distance = phistdist2)

#adonis
islands2<-c("Hawaii","Hawaii","Maui","Oahu","Oahu","Kauai","Kauai","Oahu")
pairwise.adonis(phistdist2,islands2)

#dataframe for the plot
data.scores = as.data.frame(scores(nmds))
data.scores$Sample <- gsub("\\.", " ", row.names(data.scores))
data.scores$Island <- metadata$Island[match(data.scores$Sample,metadata$Side)]
data.scores$Island[is.na(data.scores$Island)]<-"Maui"
head(data.scores)
#plot
jpeg("nmds_pairwisephist_r_n20_allislands_nonames.jpeg", width=12, height=8, units="in",res=600)
xx <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2, colour = Island,label = Sample)) + 
  geom_point(size = 10) + 
  #geom_text(hjust=0.8, vjust=-0.7,size = 6) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 18, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 20, colour = "black", face = "bold"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Island", y = "NMDS2")  + 
  scale_colour_manual(values = c("#009E73", "#E69F00","purple","darkblue")) 
print(xx)
dev.off()






################################ Random Forest ################################
library(randomForest)
library(datasets)
library(caret)
set.seed(222)
# can each replicate be assigned to an island?? site?? region??
#load data
otudf<-read.csv("cytb_ee025_r10_m90_nSCspecies_haplotype_table_filtered_noSC.csv", stringsAsFactors = F)
metadata<-read.csv("Metadata_MHI.csv",stringsAsFactors = F)
#only keep marine species
motudf<-otudf[otudf$OTU %in% mspecieslist,]
#sum haplotypes to OTUs
#motudf2<-motudf[c(5:261)]
#motudf2[2:257] <- sapply(motudf2[2:257], as.numeric)
#row.names(motudf2)<-NULL
#x<-motudf2
#x <- aggregate(x[,2:ncol(x)], by = x[1], function(x) sum(x, na.rm = TRUE))
#x2<-as.matrix(x[c(2:257)])
#x2<-t(x2)
#make matrix
otumatrix<-as.matrix(motudf[c(6:261)])
otumatrix<-t(otumatrix)
colnames(otumatrix)<-motudf$haplotype
#sum by site
#site<- metadata$Location[match(row.names(x2),metadata$Samples)]
#otumatrix2<-as.data.frame(cbind(site,x2))
#otumatrix2[2:161] <- sapply(otumatrix2[2:161], as.numeric)
#row.names(otumatrix2)<-NULL
#x<-otumatrix2
#x <- aggregate(x[,2:ncol(x)], by = x[1], function(x) sum(x, na.rm = TRUE))


# add classification variable
#island<- metadata$Island[match(x$site,metadata$Location)]
island<-metadata$Island[match(row.names(otumatrix),metadata$Samples)]
otumatrix2<-as.data.frame(cbind(island,otumatrix))
otumatrix2[2:ncol(otumatrix2)] <- sapply(otumatrix2[2:ncol(otumatrix2)], as.numeric)
otumatrix2$island<-as.factor(otumatrix2$island)
train <- otumatrix2

#cross-validation
cv<-trainControl(method = "cv",number =10) #resample 10 times
tin<-createDataPartition(train$island,p=0.8,list=FALSE) #train with 80% of the data
traindata<-train[tin, ] #subset the train data
testdata<-train[-tin, ] #subset the rest

rf_model<-train(island~., data = train, method = "rf",trControl = cv, importance = TRUE) # train
rf_model
p1 <- predict(rf_model, testdata) # predict
p1
confusionMatrix(p1, testdata$island) #get comparison
rf <- randomForest(island~., data=train, ntree=5000)


