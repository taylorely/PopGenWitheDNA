#This file is for all the correlations
##########

#load working directory and libraries
setwd("/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/Comparing_filtering") #set working directory
library("ggplot2") # load library
library("lmodel2") # load library

#######################################################################################################
######################## Different maxee threshold for genetic diversity comparisons ##################
#######################################################################################################
# load data
filterstatsJAMP<-read.csv("QualityFilterComp.csv", stringsAsFactors = F)
# for this we don't need the overall rows since it is not independent of the population values
filterstatsJAMP<-filterstatsJAMP[!(filterstatsJAMP$'Population'=="all"),] 
#subset one comparison at a time: either maxee = 0.25 or maxee = 0.5, one of presence-absence, scaled-by-reads, or raw reads proxies, and either short or long tissue sequences 
filterstatsJAMP<-filterstatsJAMP[,c('Scientific_Name',"Population",'n_short','H_short','h_short','π_short',"n_ee025_pa","H_ee025_pa","h_ee025_pa","π_ee025_pa")]
#remove incomplete comparisons
filterstatsJAMP<-na.omit(filterstatsJAMP)
#if needed remove low sample sizes
filterstatsJAMP<-filterstatsJAMP[(filterstatsJAMP$n_short>24 & filterstatsJAMP$n_ee025_pa>9),] # Remove any low samples

#test normality for each statistic
shapiro.test(filterstatsJAMP$π_short)
shapiro.test(filterstatsJAMP$π_ee025_pa)

# for each statistics: calculate spearman rank correlation, p-value, and degrees of freedom
statscor<-cor.test(filterstatsJAMP$π_short,filterstatsJAMP$π_ee025_pa, method="spearman",exact = FALSE) #run correlation
p<-round(as.numeric(statscor$p.value),digits = 6) # define and round p-value
rho<-round(as.numeric(statscor$estimate),digits = 3) # define and round spearman's rank
dfreedom<-nrow(filterstatsJAMP)-2
p
rho
dfreedom

# found maxee = 0.25 and 0.5 were very similar so we just used 0.25 from here out
##################################################################################################################################
########################## Comparing statistics, just change which columns for other statistics  #################################
##################################################################################################################################
# load data
phist_raw<-read.csv("PopgenComp.csv", stringsAsFactors = F) 
#subset one comparison at a time: one of presence-absence, scaled-by-reads, or raw reads proxies, and either short or long tissue sequences 
phist_raw<-phist_raw[c(1,3,4,6,36,38)] #keep columns for test, col 4 is pop size tissue, 6 is tissue long haplo diversity, 20 pop size pa, 22 is haplo diversity eDNA pa, 28 pop size scaled, 30 haplo diversity scaled, 36 pop size raw, 38 haplo diversity raw
#remove incomplete comparisons
phist_raw<-na.omit(phist_raw) # remove NAs
#remove the overall species values because that is not independent
phist_raw<-phist_raw[phist_raw$Population != "all",]

#For paired t-test
# Subset
tissue <- phist_raw$h_long
eDNA<- phist_raw$h_raw
jpeg(file="h_long_raw_pairedttestresults.jpg", width=8, height=6, units="in",res=600)
pd <- paired(eDNA, tissue)
resfig<-plot(pd, type = "profile") + theme_bw()
print(resfig)
dev.off()
# compute the difference
d <- eDNA - tissue
#test normality
normtest<-shapiro.test(d)
normtest
#paired t-test
wilcox.test(eDNA, tissue, paired = TRUE)
tresult<-t.test(eDNA, tissue, paired = TRUE)
#df to print
results<-data.frame(normalityp=normtest$p.value[1],tstat=tresult$statistic[1],tdf=tresult$parameter[1],tpval=tresult$p.value[1],meandiff=tresult$estimate[1])
write.csv(results,file ="h_long_raw_pairedttestresults.csv")

#remove low abundance
phist_raw<-phist_raw[(phist_raw$n_long>24 & phist_raw$n_long>24),] # Remove any low samples
phist_raw<-phist_raw[(phist_raw$n_raw>599 & phist_raw$n_raw>599),] # Remove any low samples
#test normality
shapiro.test(phist_raw$h_scaled)
#correlation
phicor<-cor.test(phist_raw$h_short, phist_raw$h_scaled, method="spearman",exact = FALSE) #run correlation
p<-round(as.numeric(phicor$p.value),digits = 3) # define and round p-value
rho<-round(as.numeric(phicor$estimate),digits = 2) # define and round spearman's rank
dfreedom<-nrow(phist_raw)-2
rho
p
dfreedom
#create the figure
jpeg(file=("h_originalshort_maxee025_scaled_e20ind_t25ind.jpeg"), width=14, height=10, units = "in", res=400) # pick settings for the figure
ggplot(data = phist_raw, aes(x=phist_raw$h_short,y=phist_raw$h_scaled))+ #pick data to use
  geom_point(aes(color = phist_raw$Scientific_Name),size = 10)+ # design points color and size
  #geom_text(label=phist_raw$Species) + # add in code if you want each point labeled by species
  xlab("h Tissue Short Sequences") + ylab("h eDNA Scaled by Reads") + #no x and y axis labels
  xlim(0,1) + ylim(0,0.7)+ #if you want to set the boundary of the axes
  geom_abline(intercept = 0, slope = 1,color="darkgrey",size=2)+ #add identity line
  annotate("text", x=0.2, y=c(0.65,0.60), label=c(paste("ρ = ",rho,sep = ""),paste("p-value = ",p,sep = "")), #add spearman rank correlation information
           color="black" ,size=12)+
  guides(color = guide_legend(title = "Species")) + #add legend for the species
  theme_bw()+ 
  theme(axis.text=element_text(size=40),axis.title=element_text(size=35),legend.text = element_text(size=15,face = "italic"),legend.title = element_text(size = 19)) #more graph design choices
dev.off()


##########################################################################################
#################################### Comparing Phist #####################################
##########################################################################################
# load data # Phist_ee025_pa.csv, PhiST_raw.csv, PhiST_ee025_scaled.csv
phist_raw<-read.csv("Phist_raw.csv", stringsAsFactors = F) 
phist_raw<-phist_raw[c(1:3,6:8,15,21:22)] #18 long #15 long # for raw
phist_raw<-phist_raw[c(1:3,7:9,16,22:23)] #19 long #16 long # for pa and scaled
#remove incomplete comparisons
phist_raw<-na.omit(phist_raw) # remove NAs
#make a column of pop comparisons
phist_raw$Islandcomb<-paste(phist_raw$Population.1,phist_raw$Population.2,sep = " and ")
# remove negatives
phist_raw$PhiST_sanger[phist_raw$PhiST_sanger<0] <- 0
phist_raw$PhiST_eDNA_pa[phist_raw$PhiST_eDNA_pa<0] <- 0
phist_raw$PhiST_eDNA_scaled[phist_raw$PhiST_eDNA_scaled<0] <- 0
phist_raw$PhiST_eDNA_raw[phist_raw$PhiST_eDNA_raw<0] <- 0

############# make bargraph #############
phist_raw<-read.csv("Phist_raw_bar_lowsampling.csv", stringsAsFactors = F) # load data # Phist_ee025_pa.csv, PhiST_raw.csv, PhiST_ee025_scaled.csv
phist_raw<-na.omit(phist_raw) # remove NAs
#make a column of pop comparisons
phist_raw$Islandcomb<-paste(phist_raw$Population.1,phist_raw$Population.2,sep = " and ")
phist_raw$IslandSpecies<-paste(phist_raw$Species,phist_raw$Islandcomb,sep = ": ")

#only long
phist_raw<-phist_raw[phist_raw$method!="Tissue Short",]
phist_raw<-phist_raw[phist_raw$method!="eDNA",]

#order the species
orderedspecies<-c("Mulloidichthys flavolineatus: Hawaii and Oahu","Mulloidichthys flavolineatus: Maui and Oahu","Mulloidichthys flavolineatus: Hawaii and Maui","Zebrasoma flavescens: Kona and Maui","Chromis vanderbilti: Maui and Kauai","Chromis vanderbilti: Oahu and Kauai","Chromis vanderbilti: Maui and Oahu","Chromis vanderbilti: Hawaii and Kauai","Chromis vanderbilti: Hawaii and Oahu","Chromis vanderbilti: Hawaii and Maui","Abudefduf abdominalis: Oahu and Kauai","Acanthurus nigrofuscus: Maui and Oahu","Acanthurus nigrofuscus: Hawaii and Oahu","Acanthurus nigrofuscus: Hawaii and Maui","Acanthurus nigrofuscus: Hawaii and Kauai","Acanthurus nigrofuscus: Maui and Kauai","Acanthurus nigrofuscus: Oahu and Kauai")
phist_raw$IslandSpecies<-factor(phist_raw$IslandSpecies,levels=orderedspecies)

#make graph
#"#90B9A1","#35775F","#52B9A6","#3F97CE","#C2E9EB","#F28170","#D3445B","#E8A9D1","#DCC8FE"
jpeg(file="pairwisephist_long_scaled_removelow_bargraph.jpg", width=10, height=11, units="in",res=600)
phibar<-ggplot(data = phist_raw, aes(x=IslandSpecies, y=PhiST.value,fill=method)) +  
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values = c("#52B9A6","#F28170") ) +
  theme(legend.position="right") +
  coord_flip()+
  theme_classic()
print(phibar)
dev.off()

############# paired t-test #############
# Subset
eDNA <- phist_raw$PhiST.value
tissue <- phist_raw$previous.long.phist
jpeg(file="pairwisephist_long_removelow_scaled_pairedttestresults.jpg", width=8, height=6, units="in",res=600)
pd <- paired(eDNA, tissue)
resfig<-plot(pd, type = "profile") + theme_bw()
print(resfig)
dev.off()
# compute the difference
d <- tissue - eDNA
#test normality
normtest<-shapiro.test(d)
normtest
#paired t-test
wresults<-wilcox.test(tissue, eDNA, paired = TRUE)
wresults
tresult<-t.test(tissue, eDNA, paired = TRUE)
tresult
#df to print
results<-data.frame(V = wresults$statistic[1],pval = wresults$p.value[1], meandiff = mean(d))
#results<-data.frame(normalityp=normtest$p.value[1],tstat=tresult$statistic[1],tdf=tresult$parameter[1],tpval=tresult$p.value[1],meandiff=tresult$estimate[1])
write.csv(results,file ="pairwisephist_long_removelow_scaled_pairedttestresults.csv")

############# correlation #############
#test normality
shapiro.test(phist_raw$PhiST_sanger)
#correlation
phicor<-cor.test(phist_raw$PhiST_sanger, phist_raw$PhiST_eDNA_scaled, method="spearman",exact = FALSE) #run correlation
p<-round(as.numeric(phicor$p.value),digits = 2) # define and round p-value
rho<-round(as.numeric(phicor$estimate),digits = 2) # define and round spearman's rank
dfreedom<-nrow(phist_raw)-2
dfreedom
p
#create the figure
jpeg(file=("phist_short_maxee025_scaled_nonegatives_rmlow_onlysharedhaplot_nomull.jpeg"), width=14, height=10, units = "in", res=400) # pick settings for the figure
ggplot(data = phist_raw, aes(x=phist_raw$PhiST_sanger,y=phist_raw$PhiST_eDNA_scaled))+ #pick data to use
  geom_point(aes(color = phist_raw$Species),size = 10)+ # design points color and size
  #geom_text(label=phist_raw$Species) + # add in code if you want each point labeled by species
  xlab("ΦST Tissue Short Sequences") + ylab("ΦST eDNA Scaled by Reads") + #no x and y axis labels
  #xlim(0,0.13) + ylim(0,0.08)+ #if you want to set the boundary of the axes
  geom_abline(intercept = 0, slope = 1,color="darkgrey",size=2)+ #add identity line
  annotate("text", x=0.02, y=c(0.05,0.045), label=c(paste("ρ = ",rho,sep = ""),paste("p-value = ",p,sep = "")), #add spearman rank correlation information
           color="black" ,size=12)+
  guides(color = guide_legend(title = "Species")) + #add legend for the species
  theme_bw()+ 
  theme(axis.text=element_text(size=40),axis.title=element_text(size=35),legend.text = element_text(size=15,face = "italic"),legend.title = element_text(size = 19)) #more graph design choices
dev.off()

##########################################################################################
################################### Chisquared test ######################################
##########################################################################################
# load data # Phist_ee025_pa.csv, PhiST_raw.csv, PhiST_ee025_scaled.csv
phist_raw<-read.csv("PhiST_ee025_scaled.csv", stringsAsFactors = F) 

#reformat data to only be island channels
phist_raw<-phist_raw[phist_raw$Population.1 != "Hilo",] # remove populations that only occur for 1 species
phist_raw<-phist_raw[phist_raw$Population.2 != "Hilo",] # remove populations that only occur for 1 species
phist_raw<-phist_raw[phist_raw$Population.1 != "Kona",] # remove populations that only occur for 1 species
phist_raw<-phist_raw[phist_raw$Population.2 != "Kona",] # remove populations that only occur for 1 species
#make a column of pop comparisons
phist_raw$Islandcomb<-paste(phist_raw$Population.1,phist_raw$Population.2,sep = " to ") 
#only keep island channels
phist_raw<-phist_raw[phist_raw$Islandcomb != "Hawaii to Kauai",]
phist_raw<-phist_raw[phist_raw$Islandcomb != "Hawaii to Oahu",]
phist_raw<-phist_raw[phist_raw$Islandcomb != "Maui to Kauai",]

#remove low abundance
#phist_raw<-phist_raw[(phist_raw$Pop.1.size.tissue>24 & phist_raw$Pop.2.size.tissue>24),] # Remove any low samples for tissue
phist_raw<-phist_raw[(phist_raw$Pop.1.size.eDNA>19 & phist_raw$Pop.2.size.eDNA>19),] # Remove any low samples for eDNA

phist_raw$fdr.tissue<-with(phist_raw, ifelse(phist_raw$previous.long.fdr.pvalue < 0.05, "yes", "no"))

#subset one island chain at a time
#phist_raw2<-phist_raw[phist_raw$Islandcomb == "Hawaii to Maui",]
phist_raw2<-phist_raw[phist_raw$Islandcomb == "Maui to Oahu",]
#phist_raw2<-phist_raw[phist_raw$Islandcomb == "Oahu to Kauai",]

#make contigency table
chitest<-as.matrix(table(phist_raw2$sig.fdr,phist_raw2$fdr.tissue))
chitest
chitest2<-matrix(c(0,1, 0, 5),nrow = 2,dimnames = list("eDNA" = c("yes", "no"),"tissue" = c("yes", "no")))
chitest2
#run mcnemars chisquared test
mcnemar.test(chitest)

##################################################################################
################################### IBD ##########################################
##################################################################################

phist_raw<-read.csv("PhiST_ee025_pa.csv", stringsAsFactors = F) # load data
phist_raw$Islandcomb<-paste(phist_raw$Population.1,phist_raw$Population.2,sep = " and ")

# remove rows that only have one population per species making it impossible for IBD comparison
phist_raw<-phist_raw[phist_raw$Species != "Caranx ignobilis",]
#phist_raw<-phist_raw[phist_raw$Species != "Caranx melampygus",]
phist_raw<-phist_raw[phist_raw$Species != "Lutjanus kasmira",]
#phist_raw<-phist_raw[phist_raw$Species != "Triaenodon obesus",]
phist_raw<-phist_raw[phist_raw$Species != "Stenella longirostris",]
phist_raw<-phist_raw[phist_raw$Species != "Chaetodon miliaris",]
#phist_raw<-phist_raw[phist_raw$Species != "Lutjanus fulvus",]
#phist_raw<-phist_raw[phist_raw$Species != "Zebrasoma flavescens",]
#phist_raw<-phist_raw[(phist_raw$Pop.1.size.eDNA>499 & phist_raw$Pop.2.size.eDNA>499),] # Remove any low samples
#phist_raw<-phist_raw[(phist_raw$Pop.1.size.tissue>24 & phist_raw$Pop.2.size.tissue>24),] # Remove any low samples

#order the graph
#for eDNA pa
phist_raw$Species <- factor(phist_raw$Species, levels=c("Abudefduf vaigiensis",  "Caranx melampygus","Mulloidichthys flavolineatus","Lutjanus fulvus","Chaetodon multicinctus","Abudefduf abdominalis","Acanthurus nigrofuscus","Cephalopholis argus","Chromis vanderbilti","Acanthurus nigroris","Zebrasoma flavescens"))

#make negatives zeros
phist_raw$PhiST.value[phist_raw$PhiST.value<0] <- 0

#plot all faceted by species
jpeg(file=("IBD_eDNA_pa_nonegative.jpeg"), width=32, height=23, units = "in", res=400)
ggplot(data = phist_raw, aes(x=phist_raw$short_dist,y=phist_raw$PhiST.value))+
  geom_point(aes(color = phist_raw$Islandcomb),size = 8)+
  #geom_text(label=phist_raw$Islandcomb) +
  xlab("Distance (km)") + ylab("ΦST from eDNA PA data") +
  #xlim(0,500) + ylim(0,0.5)+
  geom_smooth(method = "lm")+
  guides(color = guide_legend(title = "Pairwise Populations")) + 
  theme_bw()+
  theme(axis.text=element_text(size=30),axis.title=element_text(size=50),legend.text = element_text(size=30),legend.title = element_text(size = 35),strip.text = element_text(size = 30))+
  facet_wrap(~Species, scales = 'free',  nrow = 4)
dev.off()

### mantel tests for IBD
#load geographic distances
distislands4<-matrix(c(0,46,250,429,46,0,107,297,250,107,0,116,429,297,116,0), nrow = 4, ncol = 4)
distislands3<-matrix(c(0,46,250,46,0,107,250,107,0), nrow = 3, ncol = 3)
dist3HOK<-matrix(c(0,250,429,250,0,116,429,116,0), nrow = 3, ncol = 3)
distislandshilokona<-matrix(c(0,123,183,344,544,123,0,90,250,429,183,90,0,107,297,344,250,107,0,116,544,429,297,116,0), nrow = 5, ncol = 5)
distHMO<-matrix(c(0,46,250,46,0,107,250,107,0), nrow = 3, ncol = 3)
#make distance file
distislands3<-as.dist(distislands3)
distislands4<-as.dist(distislands4)
distislandshilokona<-as.dist(distislandshilokona)
dist3HOK<-as.dist(dist3HOK)
distHMO<-as.dist(distHMO)

#load raw phist values
ababdraw<-matrix(c(0,0.008038126,0.021182090,0.010980541,0.008038126,0,0.016562093,0.011004028,0.021182090,0.016562093,0,0.021174845,0.010980541,0.011004028,0.021174845,0),nrow = 4,ncol = 4)
abvairaw<-matrix(c(0,0.014709860,0.020398622,0.065576352,0.014709860,0,0.011260290,0.051643978,0.020398622,0.011260290,0,0.064814023,0.065576352,0.051643978,0.064814023,0),nrow = 4,ncol = 4)
chrvanrawa<-matrix(c(0,0.018100471,0.017779463,0.037048377,0.018100471,0,0.023734917,0.044309463,0.017779463,0.023734917,0,0.043843142,0.037048377,0.044309463,0.043843142,0),nrow = 4,ncol = 4)
mullflavraw<-matrix(c(0,0.046286554,0.064250527,0.046286554,0,0.064695121,0.064250527,0.064695121,0), nrow = 3, ncol = 3)
chaemulraw<-matrix(c(0,0.005878585,0.138244938,0.482167769,0.005878585,0,0.057344360,0.302996431,0.138244938,0.057344360,0,0.351011700,0.482167769,0.302996431,0.351011700,0), nrow = 4, ncol = 4)
zebraraw<-matrix(c(0,0.049576833,0.147438401,0.113430605,0.107735657,0.049576833,0,0.046136415,0.080584402,0.05964227,0.147438401,0.046136415,0,0.0332755,0.005947795,0.113430605,0.080584402,0.0332755,0,0.030309853,0.107735657,0.05964227,0.005947795,0.030309853,0), nrow = 5, ncol = 5)
trioberaw<-matrix(c(0,0,0.029129792,0,0,0.027015137,0.029129792,0.027015137,0), nrow = 3, ncol = 3)
#make distance files
chaemulraw<-as.dist(chaemul)
mullflavraw<-as.dist(mullflav)
ababdraw<-as.dist(ababd)
abvairaw<-as.dist(abvai)
chrvanraw<-as.dist(chrvan)
zebraraw<-as.dist(zebraraw)
trioberaw<-as.dist(trioberaw)

#load scaled phist values
abdvaiscaled<-matrix(c(0,0,0.025538211,0.0476638,0,0,0.011331394,0.019333611,0.025538211,0.011331394,0,0.025426924,0.0476638,0.019333611,0.025426924,0),nrow = 4,ncol = 4)
chrvanscaled<-matrix(c(0,0.007068351,0.009026535,0.016423909,0.007068351,0,0.013793415,0.020348161,0.009026535,0.013793415,0,0.021083333,0.016423909,0.020348161,0.021083333,0),nrow = 4,ncol = 4)
mullflavscaled<-matrix(c(0,0.02250766,0.030241935,0.033273201,0.02250766,0,0.030887223,0.050131926,0.030241935,0.030887223,0,0.047789054,0.033273201,0.050131926,0.047789054,0),nrow = 4,ncol = 4)
chaemulscaled<-matrix(c(0,0.111948332,0.454316129,0.111948332,0,0.304367922,0.454316129,0.304367922,0),nrow = 3,ncol = 3)
#make distance files
abdvaiscaled<-as.dist(abdvaiscaled)
chrvanscaled<-as.dist(chrvanscaled)
mullflavscaled<-as.dist(mullflavscaled)
chaemulscaled<-as.dist(chaemulscaled)

#load pa phist values
abdvaipa<-matrix(c(0,0,0.011256968,0.024909108,0,0,0,0,0.011256968,0,0,0.003785524,0.024909108,0,0.003785524,0),nrow = 4,ncol = 4)
chrvanpa<-matrix(c(0,0.002232202,0.003436925,0.002601887,0.002232202,0,0.005604234,0.004877605,0.003436925,0.005604234,0,0.005265246,0.002601887,0.004877605,0.005265246,0),nrow = 4,ncol = 4)
mullflavpa<-matrix(c(0,0.007738437,0.01623433,0.058175945,0.007738437,0,0.018056554,0.111662203,0.01623433,0.018056554,0,0.081669322,0.058175945,0.111662203,0.081669322,0),nrow = 4,ncol = 4)
lutfulpa<-matrix(c(0,0.019480135,0.027777778,0.019480135,0,0,0.027777778,0,0),nrow = 3,ncol = 3)
carmelpa<-matrix(c(0,5.30982849048188e-05,0,0,5.30982849048188e-05,0,0,0.000666218,0,0,0,0,0,0.000666218,0,0),nrow = 4,ncol = 4)
chaemulpa<-matrix(c(0,0.163051085,0.449131029,0.163051085,0,0.332061069,0.449131029,0.332061069,0),nrow = 3,ncol = 3)
#make distance files
abdvaipa<-as.dist(abdvaipa)
chrvanpa<-as.dist(chrvanpa)
mullflavpa<-as.dist(mullflavpa)
lutfulpa<-as.dist(lutfulpa)
carmelpa<-as.dist(carmelpa)
chaemulpa<-as.dist(chaemulpa)

#load tissue long phist values
abdvailong<-matrix(c(0,0,0,0,0,0,0,0.008719226,0,0,0,0,0,0.008719226,0,0),nrow = 4,ncol = 4)
stenlonglong<-matrix(c(0,0.018,0.024,0.018,0,0.003,0.024,0.003,0),nrow = 3,ncol = 3)
#make distance files
abdvailong<-as.dist(abdvailong)
stenlonglong<-as.dist(stenlonglong)

#load tissue short phist values
acanigRshort<-matrix(c(0,0,0.001923475,0,0,0,0.001923475,0,0),nrow = 3,ncol = 3)
#make distance files
acanigRshort<-as.dist(acanigRshort)

#run mantel test
mantel<-mantel.rtest(distHMO,trioberaw, nrepet = 9999)
#save output
capture.output(mantel,file = "mantel_chaemul_pa_IBD_9999rep.txt")


