#load working directory
setwd("/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/Comparing_filtering") #set working directory
library("ggplot2") # load library
library("lmodel2") # load library
#########################################################################
######################## genetic diversity comparisons ##################
#########################################################################
filterstatsJAMP<-read.csv("QualityFilterComp.csv", stringsAsFactors = F) # load data

filterstatsJAMP<-filterstatsJAMP[!(filterstatsJAMP$'Population'=="all"),] # for this we don't need the overall rows
#subset short and pa
filterstatsJAMP<-filterstatsJAMP[,c('Scientific_Name',"Population",'n_long','H_long','h_long','π_long',"n_ee025_raw","H_ee025_raw","h_ee025_raw","π_ee025_raw")]
#remove incomplete comparisons
filterstatsJAMP<-na.omit(filterstatsJAMP)

#function for paired t-test of 2 islands
#only keep species that were detected in both islands

#if needed remove low sample sizes
filterstatsJAMP<-filterstatsJAMP[(filterstatsJAMP$n_long>24 & filterstatsJAMP$n_ee025_raw>599),] # Remove any low samples

#test normality
shapiro.test(filterstatsJAMP$π_long)
shapiro.test(filterstatsJAMP$π_ee025_raw)

# calculate spearman rank correlation, p-value, and degrees of freedom
statscor<-cor.test(filterstatsJAMP$π_long,filterstatsJAMP$π_ee025_raw, method="spearman",exact = FALSE) #run correlation
p<-round(as.numeric(statscor$p.value),digits = 6) # define and round p-value
rho<-round(as.numeric(statscor$estimate),digits = 3) # define and round spearman's rank
dfreedom<-nrow(filterstatsJAMP)-2
p
rho
dfreedom

##########################################################################################
#################### same type of analysis but comparing eDNA to itself ##################
##########################################################################################

filterstatsJAMP<-read.csv("QualityFilterComp.csv", stringsAsFactors = F) # load data
filterstatsJAMP<-filterstatsJAMP[!(filterstatsJAMP$'Population'=="all"),] # for this we don't need the overall rows
#subset short and pa
filterstatsJAMP<-filterstatsJAMP[,c('Scientific_Name',"Population","n_ee025_raw","H_ee025_raw","h_ee025_raw","π_ee025_raw","n_ee05_raw","H_ee05_raw","h_ee05_raw","π_ee05_raw")]
#remove incomplete comparisons
filterstatsJAMP<-na.omit(filterstatsJAMP)


#if needed remove low sample sizes
filterstatsJAMP<-filterstatsJAMP[(filterstatsJAMP$n_ee025_raw>599 & filterstatsJAMP$n_ee025_raw>599),] # Remove any low samples

#test normality
shapiro.test(filterstatsJAMP$H_ee025_raw)
shapiro.test(filterstatsJAMP$H_ee025_raw)

# calculate spearman rank correlation, p-value, and degrees of freedom
statscor<-cor.test(filterstatsJAMP$π_ee025_raw,filterstatsJAMP$π_ee05_raw, method="spearman",exact = FALSE) #run correlation
p<-statscor$p.value # define and round p-value
rho<-round(as.numeric(statscor$estimate),digits = 3) # define and round spearman's rank
dfreedom<-nrow(filterstatsJAMP)-2
p
rho
dfreedom
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
#for eDNA raw
#phist_raw$Species <- factor(phist_raw$Species, levels=c("Abudefduf vaigiensis", "Chaetodon multicinctus", "Chromis vanderbilti","Mulloidichthys flavolineatus","Zebrasoma flavescens","Triaenodon obesus","Lutjanus fulvus","Chaetodon miliaris","Caranx melampygus","Abudefduf abdominalis","Acanthurus nigrofuscus","Cephalopholis argus","Acanthurus nigroris","Stenella longirostris"))
#for eDNA scaled
#phist_raw$Species <- factor(phist_raw$Species, levels=c("Abudefduf vaigiensis", "Chromis vanderbilti","Mulloidichthys flavolineatus","Chaetodon multicinctus","Lutjanus fulvus","Chaetodon miliaris","Abudefduf abdominalis","Acanthurus nigrofuscus","Cephalopholis argus","Acanthurus nigroris","Caranx melampygus","Zebrasoma flavescens","Stenella longirostris"))
#for eDNA pa
phist_raw$Species <- factor(phist_raw$Species, levels=c("Abudefduf vaigiensis",  "Caranx melampygus","Mulloidichthys flavolineatus","Lutjanus fulvus","Chaetodon multicinctus","Abudefduf abdominalis","Acanthurus nigrofuscus","Cephalopholis argus","Chromis vanderbilti","Acanthurus nigroris","Zebrasoma flavescens"))
#phist_raw$Species <- factor(phist_raw$Species, levels=c("Abudefduf vaigiensis", "Chromis vanderbilti", "Mulloidichthys flavolineatus",  "Caranx melampygus","Abudefduf abdominalis","Acanthurus nigrofuscus"))
#for tissue long
#phist_raw$Species <- factor(phist_raw$Species, levels=c("Abudefduf vaigiensis", "Stenella longirostris","Chromis vanderbilti","Caranx melampygus","Chaetodon miliaris", "Lutjanus fulvus", "Triaenodon obesus","Zebrasoma flavescens","Chaetodon multicinctus","Acanthurus nigroris","Cephalopholis argus","Abudefduf abdominalis","Mulloidichthys flavolineatus","Acanthurus nigrofuscus"))
#for tissue short
#phist_raw$Species <- factor(phist_raw$Species, levels=c("Acanthurus nigroris","Chaetodon miliaris","Zebrasoma flavescens","Acanthurus nigrofuscus","Abudefduf abdominalis","Chromis vanderbilti","Chaetodon multicinctus","Cephalopholis argus","Mulloidichthys flavolineatus","Abudefduf vaigiensis"))

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

##########################################################################################
#################################### Comparing Phist #####################################
##########################################################################################
phist_raw<-read.csv("Phist_raw.csv", stringsAsFactors = F) # load data # Phist_ee025_pa.csv, PhiST_raw.csv, PhiST_ee025_scaled.csv
phist_raw<-phist_raw[c(1:3,6:8,15,21:22)] #18 long #15 long # for raw
phist_raw<-phist_raw[c(1:3,7:9,16,22:23)] #19 long #16 long # for pa and scaled

#only shared haplotypes
#phist_raw<-read.csv("cytb_ee025_onlysharedhaplo_phist_rmlow.csv", stringsAsFactors = F) # load data # Phist_ee025_pa.csv, PhiST_raw.csv, PhiST_ee025_scaled.csv

phist_raw<-na.omit(phist_raw) # remove NAs
#make a column of pop comparisons
phist_raw$Islandcomb<-paste(phist_raw$Population.1,phist_raw$Population.2,sep = " and ")
# remove negatives
phist_raw$PhiST_sanger[phist_raw$PhiST_sanger<0] <- 0
phist_raw$PhiST_eDNA_pa[phist_raw$PhiST_eDNA_pa<0] <- 0
phist_raw$PhiST_eDNA_scaled[phist_raw$PhiST_eDNA_scaled<0] <- 0
phist_raw$PhiST_eDNA_raw[phist_raw$PhiST_eDNA_raw<0] <- 0

#make bargraph
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

#remove low abundance
#phist_raw<-subset(phist_raw, !(method == "Tissue Long" &  Pop.1.size.eDNA < 25 ))
#phist_raw<-subset(phist_raw, !(method == "Tissue Long" &  Pop.2.size.eDNA < 25))
#phist_raw<-subset(phist_raw, !(method == "Tissue Short" &  Pop.1.size.eDNA < 25))
#phist_raw<-subset(phist_raw, !(method == "Tissue Short" &  Pop.2.size.eDNA < 25))
#phist_raw<-subset(phist_raw, !(method == "eDNA" &  Pop.1.size.eDNA < 600))
#phist_raw<-subset(phist_raw, !(method == "eDNA" &  Pop.2.size.eDNA < 600))
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

#remove low abundance
#phist_raw<-phist_raw[(phist_raw$Pop.1.size.eDNA>19 & phist_raw$Pop.2.size.eDNA>19),] # Remove any low samples
#phist_raw<-phist_raw[(phist_raw$Pop.1.size.tissue>24 & phist_raw$Pop.2.size.tissue>24),] # Remove any low samples

#remove mulloidichthys
phist_raw<-phist_raw[phist_raw$Species!="Mulloidichthys_flavolineatus",]


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
################################# Comparing Tajima's D ###################################
##########################################################################################

phist_raw<-read.csv("PopgenComp.csv", stringsAsFactors = F) # load data
# long and eDNA
phist_raw<-phist_raw[c(1,3,4,8,28,32)] #keep columns for test, col 4 is pop size tissue, 8 is tajimas long, 16 is tajimas short, 20 pop size pa, 24 tajimas pa, 28 pop size scaled, 32 tajimas scaled, 36 pop size raw, 40 tajimas D raw
# short and eDNA
#phist_raw<-phist_raw[c(1,3,4,16,36,40)] #keep columns for test, col 4 is pop size tissue, 8 is tajimas long, 16 is tajimas short, 20 pop size pa, 24 tajimas pa, 28 pop size scaled, 32 tajimas scaled, 36 pop size raw, 40 tajimas D raw
phist_raw<-na.omit(phist_raw) # remove NAs
#remove the all column because that is not independent
phist_raw<-phist_raw[phist_raw$Population != "all",]

# Subset
eDNA <- phist_raw$tajima_scaled
tissue <- phist_raw$tajima_long
jpeg(file="tajimas_long_scaled_pairedttestresults.jpg", width=8, height=6, units="in",res=600)
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
wilcox.test(tissue, eDNA, paired = TRUE)
tresult<-t.test(tissue, eDNA, paired = TRUE)
#df to print
results<-data.frame(normalityp=normtest$p.value[1],tstat=tresult$statistic[1],tdf=tresult$parameter[1],tpval=tresult$p.value[1],meandiff=tresult$estimate[1])
write.csv(results,file ="tajimas_long_scaled_pairedttestresults.csv")


#remove low abundance
phist_raw<-phist_raw[(phist_raw$n_long>24 & phist_raw$n_long>24),] # Remove any low samples
phist_raw<-phist_raw[(phist_raw$n_scaled>19 & phist_raw$n_scaled>19),] # Remove any low samples
#test normality
shapiro.test(phist_raw$tajima_raw)
#correlation
phicor<-cor.test(phist_raw$tajima_short, phist_raw$tajima_raw, method="spearman",exact = FALSE) #run correlation
p<-round(as.numeric(phicor$p.value),digits = 3) # define and round p-value
rho<-round(as.numeric(phicor$estimate),digits = 2) # define and round spearman's rank
dfreedom<-nrow(phist_raw)-2
rho
p
dfreedom
#create the figure
jpeg(file=("tajimasD_originalshort_maxee025_raw_e0ind_t0ind.jpeg"), width=16, height=10, units = "in", res=400) # pick settings for the figure
ggplot(data = phist_raw, aes(x=phist_raw$tajima_short,y=phist_raw$tajima_raw))+ #pick data to use
  geom_point(aes(color = phist_raw$Scientific_Name),size = 8)+ # design points color and size
  #geom_text(label=phist_raw$Species) + # add in code if you want each point labeled by species
  xlab("Tajima's D Tissue Full Length") + ylab("Tajima's D eDNA raw reads") + #no x and y axis labels
  xlim(-2,2) + ylim(-2,2)+ #if you want to set the boundary of the axes
  geom_abline(intercept = 0, slope = 1,color="grey")+ #add identity line
  annotate("text", x=0, y=c(1.5,1), label=c(paste("rho = ",rho,sep = ""),paste("p = ",p,sep = "")), #add spearman rank correlation information
           color="black" ,size=20)+
  guides(color = guide_legend(title = "Species")) + #add legend for the species
  theme_bw()+ 
  theme(axis.text=element_text(size=40),axis.title=element_text(size=40),legend.text = element_text(size=30,face = "italic"),legend.title = element_text(size = 35)) #more graph design choices
dev.off()
##########################################################################################
########################## Comparing Haplotype Diversity #################################
##########################################################################################

phist_raw<-read.csv("PopgenComp.csv", stringsAsFactors = F) # load data
# long and eDNA
phist_raw<-phist_raw[c(1,3,4,6,36,38)] #keep columns for test, col 4 is pop size tissue, 6 is tissue long haplo diversity, 20 pop size pa, 22 is haplo diversity eDNA pa, 28 pop size scaled, 30 haplo diversity scaled, 36 pop size raw, 38 haplo diversity raw
# short and eDNA
#phist_raw<-phist_raw[c(1,3,4,14,36,38)] #14 tissue short haplo diversity
phist_raw<-na.omit(phist_raw) # remove NAs
#remove the all column because that is not independent
phist_raw<-phist_raw[phist_raw$Population != "all",]

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
########################## Comparing Nucleotide Diversity #################################
##########################################################################################

phist_raw<-read.csv("PopgenCompNucDiv6rounding.csv", stringsAsFactors = F) # load data
# long and eDNA
phist_raw<-phist_raw[c(1,3,4,15,36,39)] #keep columns for test, col 4 is pop size tissue, 7 is tissue long nucleo diversity, 20 pop size pa, 23 is diversity eDNA pa, 28 pop size scaled, 31 diversity scaled, 36 pop size raw, 39  diversity raw
# short and eDNA
#phist_raw<-phist_raw[c(1,3,4,15,20,23)] #15 tissue short nucleo diversity
phist_raw<-na.omit(phist_raw) # remove NAs
#remove the all column because that is not independent
phist_raw<-phist_raw[phist_raw$Population != "all",]
#remove low abundance
phist_raw<-phist_raw[(phist_raw$n_long>24 & phist_raw$n_long>24),] # Remove any low samples
phist_raw<-phist_raw[(phist_raw$n_raw>599 & phist_raw$n_raw>599),] # Remove any low samples

#remove lutkas
phist_raw<-phist_raw[(phist_raw$Scientific_Name != "Lutjanus kasmira"),]

# Subset
tissue <- phist_raw$π_short
eDNA<- phist_raw$π_raw
jpeg(file="nucleotidediv_short_raw_pairedttestresults.jpg", width=8, height=6, units="in",res=600)
pd <- paired(eDNA, tissue)
resfig<-plot(pd, type = "profile") + theme_bw()
print(resfig)
dev.off()
# compute the difference
d <-  tissue - eDNA 
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
write.csv(results,file ="nucleotidediv_short_raw_pairedttestresults.csv")

#test normality
shapiro.test(phist_raw$π_scaled)
#correlation
phicor<-cor.test(phist_raw$π_long, phist_raw$π_scaled, method="spearman",exact = FALSE) #run correlation
p<-round(as.numeric(phicor$p.value),digits = 3) # define and round p-value
rho<-round(as.numeric(phicor$estimate),digits = 2) # define and round spearman's rank
dfreedom<-nrow(phist_raw)-2
rho
p
p<-format(p, scientific = FALSE)
dfreedom
#create the figure
jpeg(file=("π_originallong_maxee025_scaled_e20ind_t25ind.jpeg"), width=14, height=10, units = "in", res=400) # pick settings for the figure
ggplot(data = phist_raw, aes(x=phist_raw$π_long,y=phist_raw$π_scaled))+ #pick data to use
  geom_point(aes(color = phist_raw$Scientific_Name),size = 10)+ # design points color and size
  #geom_text(label=phist_raw$Species) + # add in code if you want each point labeled by species
  xlab("π Tissue Long Sequences") + ylab("π eDNA Scaled by Reads") + #no x and y axis labels
  xlim(0,0.01) + ylim(0,0.007)+ #if you want to set the boundary of the axes
  geom_abline(intercept = 0, slope = 1,color="darkgrey",size=2)+ #add identity line
  annotate("text", x=0.008, y=c(0.006,0.0055), label=c(paste("ρ = ",rho,sep = ""),paste("p-value = ",p,sep = "")), #add spearman rank correlation information
           color="black" ,size=12)+
  guides(color = guide_legend(title = "Species")) + #add legend for the species
  theme_bw()+ 
  theme(axis.text=element_text(size=40),axis.title=element_text(size=35),legend.text = element_text(size=15,face = "italic"),legend.title = element_text(size = 19)) #more graph design choices
dev.off()

##########################################################################################
################################## Comparing R2 ##########################################
##########################################################################################

phist_raw<-read.csv("PopgenComp.csv", stringsAsFactors = F) # load data
# long and eDNA
#phist_raw<-phist_raw[c(1,3,4,10,28,34)] #keep columns for test, col 4 is pop size tissue, 10 is tissue long R2, 20 pop size pa, 26 is R2 pa, 28 pop size scaled, 34 R2 scaled, 36 pop size raw, 42  R2 raw
# short and eDNA
phist_raw<-phist_raw[c(1,3,4,18,28,34)] #18 tissue short R2
phist_raw<-na.omit(phist_raw) # remove NAs
#remove the all column because that is not independent
phist_raw<-phist_raw[phist_raw$Population != "all",]
#remove low abundance
phist_raw<-phist_raw[(phist_raw$n_long>24 & phist_raw$n_long>24),] # Remove any low samples
phist_raw<-phist_raw[(phist_raw$n_scaled>19 & phist_raw$n_scaled>19),] # Remove any low samples

# Subset
tissue <- phist_raw$R2_short
eDNA<- phist_raw$R2_scaled
jpeg(file="R2_short_removelow_scaled_pairedttestresults.jpg", width=8, height=6, units="in",res=600)
pd <- paired(eDNA, tissue)
resfig<-plot(pd, type = "profile") + theme_bw()
print(resfig)
dev.off()
# compute the difference
d <-  tissue - eDNA 
#test normality
normtest<-shapiro.test(d)
normtest
#paired t-test
wresults<-wilcox.test(tissue, eDNA, paired = TRUE)
wresults
tresult<-t.test(tissue, eDNA, paired = TRUE)
tresult
#df to print
#results<-data.frame(V = wresults$statistic[1],pval = wresults$p.value[1], meandiff = mean(d))
results<-data.frame(normalityp=normtest$p.value[1],tstat=tresult$statistic[1],tdf=tresult$parameter[1],tpval=tresult$p.value[1],meandiff=tresult$estimate[1])
write.csv(results,file ="R2_short_removelow_scaled_pairedttestresults.csv")


#test normality
shapiro.test(phist_raw$R2_raw)
#correlation
phicor<-cor.test(phist_raw$R2_long, phist_raw$R2_raw, method="spearman",exact = FALSE) #run correlation
p<-round(as.numeric(phicor$p.value),digits = 3) # define and round p-value
rho<-round(as.numeric(phicor$estimate),digits = 2) # define and round spearman's rank
dfreedom<-nrow(phist_raw)-2
rho
p
dfreedom
#create the figure
jpeg(file=("R2_originallong_maxee025_raw_e1500ind_t25ind.jpeg"), width=16, height=10, units = "in", res=400) # pick settings for the figure
ggplot(data = phist_raw, aes(x=phist_raw$R2_long,y=phist_raw$R2_raw))+ #pick data to use
  geom_point(aes(color = phist_raw$Scientific_Name),size = 8)+ # design points color and size
  #geom_text(label=phist_raw$Species) + # add in code if you want each point labeled by species
  xlab("R2 Tissue Full Length") + ylab("R2 eDNA Raw Reads") + #no x and y axis labels
  #xlim(0,0.18) + ylim(0,0.18)+ #if you want to set the boundary of the axes
  geom_abline(intercept = 0, slope = 1,color="grey")+ #add identity line
  annotate("text", x=0.14, y=c(0.07,0.06), label=c(paste("rho = ",rho,sep = ""),paste("p = ",p,sep = "")), #add spearman rank correlation information
           color="black" ,size=20)+
  guides(color = guide_legend(title = "Species")) + #add legend for the species
  theme_bw()+ 
  theme(axis.text=element_text(size=40),axis.title=element_text(size=40),legend.text = element_text(size=30,face = "italic"),legend.title = element_text(size = 35)) #more graph design choices
dev.off()

##########################################################################################
################################## Comparing r ##########################################
##########################################################################################

phist_raw<-read.csv("PopgenComp.csv", stringsAsFactors = F) # load data
# long and eDNA
#phist_raw<-phist_raw[c(1,3,4,11,28,35)] #keep columns for test, col 4 is pop size tissue, 11 is tissue long R2, 20 pop size pa, 27 is r pa, 28 pop size scaled, 35 r scaled, 36 pop size raw, 43  r raw
# short and eDNA
phist_raw<-phist_raw[c(1,3,4,19,28,35)] #19 tissue short r
phist_raw<-na.omit(phist_raw) # remove NAs
#remove the all column because that is not independent
phist_raw<-phist_raw[phist_raw$Population != "all",]
#remove low abundance
phist_raw<-phist_raw[(phist_raw$n_long>24 & phist_raw$n_long>24),] # Remove any low samples
phist_raw<-phist_raw[(phist_raw$n_scaled>19 & phist_raw$n_scaled>19),] # Remove any low samples

# Subset
tissue <- phist_raw$r_short
eDNA<- phist_raw$r_scaled
jpeg(file="r_short_removelow_scaled_pairedttestresults.jpg", width=8, height=6, units="in",res=600)
pd <- paired(eDNA, tissue)
resfig<-plot(pd, type = "profile") + theme_bw()
print(resfig)
dev.off()
# compute the difference
d <-  tissue - eDNA 
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
write.csv(results,file ="r_short_removelow_scaled_pairedttestresults.csv")


#test normality
shapiro.test(phist_raw$r_raw)
#correlation
phicor<-cor.test(phist_raw$r_long, phist_raw$r_raw, method="spearman",exact = FALSE) #run correlation
p<-round(as.numeric(phicor$p.value),digits = 3) # define and round p-value
rho<-round(as.numeric(phicor$estimate),digits = 2) # define and round spearman's rank
dfreedom<-nrow(phist_raw)-2
rho
p
dfreedom
#create the figure
jpeg(file=("r_originallong_maxee025_raw_e600ind_t25ind.jpeg"), width=16, height=10, units = "in", res=400) # pick settings for the figure
ggplot(data = phist_raw, aes(x=phist_raw$r_long,y=phist_raw$r_raw))+ #pick data to use
  geom_point(aes(color = phist_raw$Scientific_Name),size = 8)+ # design points color and size
  #geom_text(label=phist_raw$Species) + # add in code if you want each point labeled by species
  xlab("r Tissue Full Length") + ylab("r eDNA Raw Reads") + #no x and y axis labels
  #xlim(0,0.18) + ylim(0,0.18)+ #if you want to set the boundary of the axes
  geom_abline(intercept = 0, slope = 1,color="grey")+ #add identity line
  annotate("text", x=0.13 , y=c(0.7,0.6), label=c(paste("rho = ",rho,sep = ""),paste("p = ",p,sep = "")), #add spearman rank correlation information
           color="black" ,size=20)+
  guides(color = guide_legend(title = "Species")) + #add legend for the species
  theme_bw()+ 
  theme(axis.text=element_text(size=40),axis.title=element_text(size=40),legend.text = element_text(size=30,face = "italic"),legend.title = element_text(size = 35)) #more graph design choices
dev.off()
##########################################################################################
################################## Comparing Fu's Fs ##########################################
##########################################################################################

phist_raw<-read.csv("PopgenComp.csv", stringsAsFactors = F) # load data
# tissue and eDNA
phist_raw<-phist_raw[c(1,3,4,9,28,33)] #keep columns for test, col 4 is pop size tissue, 9 is tissue long Fs, 17 tissue short Fs, 21 pop size pa, 25 is fs pa, 28 pop size scaled, 33 fs scaled
#phist_raw<-phist_raw[c(1,3,4,17,28,33)] #keep columns for test, col 4 is pop size tissue, 9 is tissue long Fs, 17 tissue short Fs, 21 pop size pa, 25 is fs pa, 28 pop size scaled, 33 fs scaled
phist_raw<-na.omit(phist_raw) # remove NAs
#remove the all column because that is not independent
phist_raw<-phist_raw[phist_raw$Population != "all",]
#remove low abundance
phist_raw<-phist_raw[(phist_raw$n_long>24 & phist_raw$n_long>24),] # Remove any low samples
phist_raw<-phist_raw[(phist_raw$n_scaled>19 & phist_raw$n_scaled>19),] # Remove any low samples

# Subset
tissue <- phist_raw$fs_long
eDNA<- phist_raw$fs_scaled
jpeg(file="fusfs_long_removelow_scaled_pairedttestresults.jpg", width=8, height=6, units="in",res=600)
pd <- paired(eDNA, tissue)
resfig<-plot(pd, type = "profile") + theme_bw()
print(resfig)
dev.off()
# compute the difference
d <-  tissue - eDNA 
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
results<-data.frame(normalityp=normtest$p.value[1],tstat=tresult$statistic[1],tdf=tresult$parameter[1],tpval=tresult$p.value[1],meandiff=tresult$estimate[1])
write.csv(results,file ="fusfs_long_removelow_scaled_pairedttestresults.csv")

#test normality
shapiro.test(phist_raw$fs_pa)
#correlation
phicor<-cor.test(phist_raw$fs_long, phist_raw$fs_pa, method="spearman",exact = FALSE) #run correlation
p<-round(as.numeric(phicor$p.value),digits = 3) # define and round p-value
rho<-round(as.numeric(phicor$estimate),digits = 2) # define and round spearman's rank
dfreedom<-nrow(phist_raw)-2
rho
p
dfreedom
#create the figure
jpeg(file=("fs_originallong_maxee025_pa_e25ind_t25ind.jpeg"), width=16, height=10, units = "in", res=400) # pick settings for the figure
ggplot(data = phist_raw, aes(x=phist_raw$fs_long,y=phist_raw$fs_pa))+ #pick data to use
  geom_point(aes(color = phist_raw$Scientific_Name),size = 8)+ # design points color and size
  #geom_text(label=phist_raw$Species) + # add in code if you want each point labeled by species
  xlab("Fu's Fs Tissue Full Length") + ylab("Fu's Fs eDNA PA") + #no x and y axis labels
  #xlim(0,0.18) + ylim(0,0.18)+ #if you want to set the boundary of the axes
  geom_abline(intercept = 0, slope = 1,color="grey")+ #add identity line
  annotate("text", x=-4 , y=c(-20,-23), label=c(paste("rho = ",rho,sep = ""),paste("p = ",p,sep = "")), #add spearman rank correlation information
           color="black" ,size=20)+
  guides(color = guide_legend(title = "Species")) + #add legend for the species
  theme_bw()+ 
  theme(axis.text=element_text(size=40),axis.title=element_text(size=40),legend.text = element_text(size=30,face = "italic"),legend.title = element_text(size = 35)) #more graph design choices
dev.off()


##########################################################################################
################################### Chisquared test ######################################
##########################################################################################
phist_raw<-read.csv("PhiST_ee025_scaled.csv", stringsAsFactors = F) # load data

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

#get expected from tissue
#islchan<-c("O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hilo and Maui","Kona and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","O'ahu and Kaua'i","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu")
#table(islchan)
#tfexpok<-1/13
#tfexpom<-1/10
#tfexpmh<-1/10
#get the yes/no for tissue
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

############## compare to all eDNA data ############
phist_raw_all<-read.csv("cytb_ee025_m90_nSC_site_r_phist_table_manyspecies.csv", stringsAsFactors = F) 
#make a column of pop comparisons
phist_raw_all$Islandcomb<-paste(phist_raw_all$Population.1,phist_raw_all$Population.2,sep = " to ") 
#only keep island channels
phist_raw_all<-phist_raw_all[phist_raw_all$Islandcomb != "Hawaii to Kauai",]
phist_raw_all<-phist_raw_all[phist_raw_all$Islandcomb != "Hawaii to Oahu",]
phist_raw_all<-phist_raw_all[phist_raw_all$Islandcomb != "Maui to Kauai",]

#subset one island chain at a time
#phist_raw_all2<-phist_raw_all[phist_raw_all$Islandcomb == "Hawaii to Maui",]
#phist_raw_all2<-phist_raw_all[phist_raw_all$Islandcomb == "Maui to Oahu",]
phist_raw_all2<-phist_raw_all[phist_raw_all$Islandcomb == "Oahu to Kauai",]

chitest<-as.data.frame(table(phist_raw_all$sig.fdr,phist_raw_all$Islandcomb)) #sum info
colnames(chitest)<-c("sig","Islandcomb","freq") #rename columns
chitest$totalsamplesize<-sum(chitest$freq) #get total sample size
samplesize<-aggregate(chitest$freq, by=list(Category=chitest$Islandcomb), FUN=sum) #get by break sample size
colnames(samplesize)<-c("Islandcomb","SampleSize")
chitest<-chitest[chitest$sig=="yes",] #only keep significant rows
chitest2<-merge(chitest,samplesize,by="Islandcomb") #merge together
#expected if same ratios as subset of 18 species
#scaled by reads expected all: HM:0/8, MO:1/9, OK:3/14
expectedn<-c(0,1,3)
expectedt<-c(8,9,14)
chitest2$expected<-expectedn/expectedt
chitest2$expectedcount<-chitest2$expected*chitest2$SampleSize
#expected if equal across island channels
#chitest2$expected<-chitest2$SampleSize/chitest2$totalsamplesize #calculate weight of sample size
#chitest2$expectedcount<-(chitest2$SampleSize/chitest2$totalsamplesize)*sum(chitest2$freq) #calculate expected sig breaks
chiresults<-chisq.test(chitest2$freq, p = chitest2$expected,rescale.p=TRUE) #chisquared test
chitest3<-data.frame(chitest2$freq,chitest2$expectedcount)
chiresults2<-chisq.test(chitest3)
chiresults$expected

####################### for zeros do "nos" ##########################
phist_raw<-read.csv("PhiST_raw.csv", stringsAsFactors = F) # load data
#phist_raw<-read.csv("cytb_ee025_m90_nSC_site_r_phist_table_manyspecies.csv", stringsAsFactors = F) 

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
phist_raw<-phist_raw[(phist_raw$Pop.1.size.eDNA>599 & phist_raw$Pop.2.size.eDNA>599),] # Remove any low samples for eDNA

#get expected from tissue
#islchan<-c("O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hilo and Maui","Kona and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","O'ahu and Kaua'i","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu","O'ahu and Kaua'i","O'ahu and Kaua'i","Hawai'i Island and Maui","Maui and O'ahu")
#table(islchan)
tfexpok<-12/13
tfexpom<-9/10
tfexpmh<-9/10


#sum info
chitest<-as.data.frame(table(phist_raw$sig.fdr,phist_raw$Islandcomb)) 
colnames(chitest)<-c("sig","Islandcomb","freq") #rename columns
chitest$totalsamplesize<-sum(chitest$freq) #get total sample size
samplesize<-aggregate(chitest$freq, by=list(Category=chitest$Islandcomb), FUN=sum) #get by break sample size
colnames(samplesize)<-c("Islandcomb","SampleSize")
chitest<-chitest[chitest$sig=="yes",] #only keep significant rows
chitest2<-merge(chitest,samplesize,by="Islandcomb") #merge together
#get expected
chitest2$expected_long<-c(tfexpmh,tfexpom,tfexpok)
chitest2$expected_long_count<-chitest2$expected_long*chitest2$SampleSize
#chitest2$expected<-chitest2$SampleSize/chitest2$totalsamplesize #calculate weight of sample size
#chitest2$expectedcount<-(chitest2$SampleSize/chitest2$totalsamplesize)*sum(chitest2$freq) #calculate expected sig breaks
chitest2
#run chisquared
chiresults<-chisq.test(chitest2$freq, p = chitest2$expected_long_count, rescale.p = TRUE) 
chiresults



####### chisquared for everything
phist_raw<-read.csv("cytb_ee025_m90_nSC_site_pa_phist_table_manyspecies.csv", stringsAsFactors = F) # load data
#make a column of pop comparisons
phist_raw$Islandcomb<-paste(phist_raw$Population.1,phist_raw$Population.2,sep = " to ") 
#only keep island channels
phist_raw<-phist_raw[phist_raw$Islandcomb != "Hawaii to Kauai",]
phist_raw<-phist_raw[phist_raw$Islandcomb != "Hawaii to Oahu",]
phist_raw<-phist_raw[phist_raw$Islandcomb != "Maui to Kauai",]

chitest<-as.data.frame(table(phist_raw$sig.fdr,phist_raw$Islandcomb)) #sum info
colnames(chitest)<-c("sig","Islandcomb","freq") #rename columns
chitest$totalsamplesize<-sum(chitest$freq) #get total sample size
samplesize<-aggregate(chitest$freq, by=list(Category=chitest$Islandcomb), FUN=sum) #get by break sample size
colnames(samplesize)<-c("Islandcomb","SampleSize")
chitest<-chitest[chitest$sig=="yes",] #only keep significant rows
chitest2<-merge(chitest,samplesize,by="Islandcomb") #merge together
chitest2$expected<-chitest2$SampleSize/chitest2$totalsamplesize #calculate weight of sample size
chitest2$expectedcount<-(chitest2$SampleSize/chitest2$totalsamplesize)*sum(chitest2$freq) #calculate expected sig breaks
chiresults<-chisq.test(chitest2$freq, p = chitest2$expected) #chisquared test
chiresults$expected


