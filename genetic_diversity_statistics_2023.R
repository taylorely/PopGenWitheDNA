################## Notes/Instructions ############################
########## PLEASE READ SO YOU DONT HAVE ISSUES ################

#This code uses the following outputs from bioinformatics: 
#1) the haplo table output from the denoising step of JAMP 
#2) a tab file from ecotag/obitools
#3) a file created from taking the top match from blastn search using blastn -remote and awk code (awk 'c&&!--c;/Query/{printf "%s	%s", $0, c=5}' input.txt >> output.txt). Make sure to use find and replace to remove "# Query: " from each line
#4) an alignment file created after the first function load_species(). You can also just check to see if they are already aligned. In that case just load the fasta file created by load_species()

#With these outputs you also need more info about your data:  ###Side note: Make sure these are correct. Often there is an error because the sample name in the metadata file did not match the sample names in the JAMP output so double check this is correct. 
#1) a metadata file that links the sample names to Site, Island, Region, Side of the island, and year


#For functions you often need to include the following inputs (and the files mentioned above)
#1) runname or whatever you want to be a name modifier for your run so you don't overwrite different runs
#2) right now the code only works for using OTUs but code easily modify for species (choose from "OTU" or "SpeciesName") 
#3) what is the lowest taxonomy match you want to keep, this is to remove things that don't match well to anything and are likely bacteria or other non-target taxa
#4) lowest number of reads you want to keep, usually pick something from 1-10, this is to remove potential false positives
#5) number of blanks
#6) number of samples
#7) the folder that you want all files saved to. If you aren't doing multiple analyses just set to the same as working directory
#9) What frame the codon starts on 
#10) What your limit is for low sample size per population
#11) How you want to scale the reads for converting to a 0-4 scale (e.g. Turon et al. 2020 used <0.5,0.5<x<0.75, 0.75<x<0.9, >0.9  to divide into 0-4 scale)
#12) if you want to remove first nucleotide (TRUE or FALSE). Because of sequencing error this was often a change never observed in haplotypes extracted from tissue
#13) this code is meant for when you have multiple levels of population (Island vs side of the Island) so you sometimes have to set what population you are using

### note if having issues ##### There are some package issues with phyloseq and these ones so clear your environment and History or restart your R session if there are issues ########

################### Here is where the code starts ######################
# load libraries (not sure if all of these are still needed)
library("dplyr")
library("ape")
library("pegas")
library("reshape2")
library("iNEXT")
library("devtools")
library("ggplot2")
library("haplotypes")
library("stringr")
library("ade4")
library("taxize")
library("vegan")
library("seqinr")
library("mmod")
library("hierfstat")
library("adegenet")
library("graph4lg")
library("vegan")
library("phyloseq")
library("ape")
library("ggplot2")
library("treemapify")
library("plyr")
library("poppr")
library("tidyr")

################## here are the many functions #######################

#### functions for formating data
#load data and create dataframes and create proxies for individuals
load_species<-function(metadata,table,obitools,blast,lowestmatch,lowestreads,OTUvsspecies="OTU",folder,nblanks,nsamples,runname,scaled1,scaled2,scaled3,removefirstbp){
  ########### reformatting ############
  #remove blanks and create a subset file to see what had contamination
  if (nblanks > 0) {
  blanks<-metadata$Samples[metadata$Side=="Blank"]
  if (nblanks > 1) {contam<-table[rowSums(table[,c(blanks[1:nblanks])]!=0)>0,] #keep rows with something in the blanks
  } else if (nblanks==1) {contam<-table[table[blanks]>0,]}
  write.csv(contam,file=paste(folder,"/",runname,"_reads_in_blank.csv",sep = "")) #write these rows to a table
  table1 <- table[,!(names(table) %in% blanks)] #remove blanks from table
  }
  #get a list of total reads per sample which will later be used to rarefaction
  if (nblanks > 0) {
  samplecolumns<-metadata$Samples[metadata$Side!="Blank"]
  } else {samplecolumns<-metadata$Samples}
  sample_totals<-colSums(table1[,c(samplecolumns[1:nsamples])])
  
  # remove the line that says how many reads were removed
  lastrow<-(nrow(table1)-1)
  table1<-table1[1:lastrow,] 
  
  #remove first nucleotide (I decided to do this because it seems there is an error that some of the sequences have a G instead of A which is never seen in the tissue derived samples of longer sequence length)
  if(removefirstbp == TRUE) {
  table1$sequences <- as.character(table1$sequences)
  table1$sequences <- gsub("^.{1}", "", table1$sequences)
  }
  
  #remove sequences of different lengths
  otus <- unique(table1$OTU) #load a list of OTU numbers to refer to later
  i<-1
  id_all<-c()
  for (i in 1:length(otus)){
    temp <- table1[table1$OTU==otus[i],]
    x <- as.list(temp$sequences)
    bp<-nchar(x)
    temp$bp<-bp
    temp<-temp[(temp$bp != bp[1]),]
    id<-temp$haplotype
    id_all<-append(id_all,id)
  }
  if (length(id_all) != 0){
  for (i in 1:length(id_all)){
    table1<-table1[(table1$haplotype != id_all[i]),]
  }}
  
  #remove samples that had less than x reads and rows with all 0 reads
  table1[,c(samplecolumns[1:nsamples])][table1[,c(samplecolumns[1:nsamples])] <= lowestreads] <- 0 #replace low values with zero
  table1<-table1[rowSums(table1[,c(samplecolumns[1:nsamples])]!=0)>0,] #remove rows with all zeros
  
  # create combined file
  table1$id<-paste(table1$OTU,table1$haplotype,sep = "_")
  obitools<-obitools[1:lastrow,c(1,3:4,7,9:13,15:17)]
  colnames(obitools)<-c("id","best_identity","best_match","Family","Genus","Match_Count","Order","Rank","Scientific_name","SpeciesList","SpeciesName","TaxID")
  blast<-blast[c(1,4:8)]
  colnames(blast)<-c("id","blastspecies","blastID","evalue","querycover","blastpercent") # rename so not V1,V2...
  combined_first<-merge(obitools,blast, by="id",all.x=TRUE)
  combined<-merge(table1,combined_first,by="id",all.x=TRUE)
  combinedraw<-combined
  write.csv(combined,file=paste(folder,"/",runname,"_haplotypes_raw.csv",sep = "")) # this will be used for later comparisons
  
  ########### Remove non target taxa ############
  #remove non target taxa including non marine metazoans: humans, bacteria, algae, chicken, dog, ...
  combined<-combined[!grepl("9606", combined$blastID),] # remove any rows with ...
  combined<-combined[!grepl("Bilateria|Eumetazoa|Eukaryota|Bacteria|root|Bacilli|Philasterida|Flustrina|Tetrahymena|Lactobacillales|Ciliophora", combined$Scientific_name),]
  combined<-combined[!grepl("Galliformes", combined$Order),]
  combined<-combined[!grepl("Canidae|Bovidae|Suidae|Muridae", combined$Family),]
  combined<-combined[!grepl("strain",combined$Rank),]
  combined<-combined[!grepl("33317",combined$TaxID),] 
  
  ########### Remove haplotypes that didn't match too well ###########
  combined$blastpercent<-as.numeric(combined$blastpercent) # set percentages to numerics so we can use the following code
  combined3<-subset(combined, combined$best_identity>=lowestmatch & combined$blastpercent>=(lowestmatch*100))
  
  # remove the same ones from the table file not just the combined file
  ids<-combined3$id
  idremove<-unique(table1$id)
  for (i in 1:length(ids)) { # get a list of OTUs to remove
    idremove<-idremove[(idremove != ids[i])];
  }
  if (length(idremove)!=0){
  for (i in 1:length(idremove)) { # remove non target taxa from table
    table1<-table1[(table1$id != idremove[i]),]
  }}
  
  
  #get lowest shared taxonomy per OTU for blast searches, have this already for the ecopcr stuff 
  otus2<-unique(combined3$OTU) #get list of OTUs
  blast_lowest<-data.frame(matrix(ncol = 5, nrow = 0)) #create dataframe
  phylum<-c()
  i<-1 
  for (i in 1:length(otus2)) { # make a loop for getting lowest shared taxonomy from blast results
    temp <- combined3[combined3$OTU==otus2[i],]
    lowest<-c()
    try(lowest<-lowest_common(temp$blastID,db="ncbi"),  silent = TRUE)
    lowest<-append(lowest,c(0,0,0))
    lowest[1]<-as.character(lowest[1])
    try(phylum<-tax_name(sci = lowest[1], get = "phylum", db = "ncbi"), silent=TRUE)
    phylum<-append(phylum,c(0,0,0))
    blastlowest2<-c(otus2[i],lowest[1],lowest[2],lowest[3],phylum[3])
    blastlowest2<-unname(blastlowest2)
    blast_lowest<- rbind(blast_lowest,blastlowest2)
  }
  colnames(blast_lowest)<-c("OTU","name","rank","ID","phylum")
  combined3$blastlowest<-blast_lowest$name[match(combined3$OTU,blast_lowest$OTU)] 
  combined3$blastlowestid<-blast_lowest$ID[match(combined3$OTU,blast_lowest$OTU)] 
  combined3$blastlowestrank<-blast_lowest$rank[match(combined3$OTU,blast_lowest$OTU)] 
  combined3$phylum<-blast_lowest$phylum[match(combined3$OTU,blast_lowest$OTU)] 

  #save dataset
  final_combined<-combined3
  #write dataset to csv file
  write.csv(final_combined,file = paste(folder,"/",runname,"species_haplotype_table_filtered.csv",sep = ""))
  
  #write sequences as fasta file for the remove stop codon function
  seqs1<-as.list(combined3$sequences)
  names1<- as.list(combined3$haplotype)
  write.fasta(sequences=seqs1,names=names1,file = paste(folder,"/",runname,"haplotypes.fasta",sep = ""))

  ########### create a dataset with rarefied data ############
  endcolumn<-nsamples+3
  ttestrar<-t(combined3[5:(endcolumn+1)])
  S <- specnumber(ttestrar) # observed number of species
  Srare <- t(rrarefy(ttestrar, min(sample_totals))) #rarefy with the min read count
  rarefy_reads<-combined3[c(2:(nsamples+6),(nsamples+11),(nsamples+8),(nsamples+9),(nsamples+13),(nsamples+22),(nsamples+21),(nsamples+25))] 
  rarefy_reads[4:endcolumn]<-Srare
  rarefy_reads<-rarefy_reads[rowSums(rarefy_reads[,c(samplecolumns[1:nsamples])]!=0)>0,]
  combined_rare<-rarefy_reads
  write.csv(rarefy_reads,file=paste(folder,"/",runname,"_reads_rarefied.csv",sep = ""))
  
  # determine quantile scale
  rank_quantile<-data.frame(matrix(ncol=4, nrow = 0))
  otus<-unique(rarefy_reads[[OTUvsspecies]])
  for (i in 1:length(otus)) {
    raretable<-rarefy_reads[rarefy_reads[[OTUvsspecies]]==otus[i],]
    raretable<-raretable[4:endcolumn]
    raretable<-unlist(raretable)
    raretable<-raretable[raretable != 0]
    raretable<-as.vector(raretable)
    rank_test<-quantile(raretable,prob=c(scaled1,scaled2,scaled3))
    ranks<-c(otus[i],rank_test)
    rank_quantile<-rbind(rank_quantile,ranks)
  }
  colnames(rank_quantile)<-c("OTU","q1","q2","q3")
  rarefy_reads$q1 <- rank_quantile$q1[match(rarefy_reads$OTU,rank_quantile$OTU)] 
  rarefy_reads$q2 <- rank_quantile$q2[match(rarefy_reads$OTU,rank_quantile$OTU)] 
  rarefy_reads$q3 <- rank_quantile$q3[match(rarefy_reads$OTU,rank_quantile$OTU)] 
  
  #replace read counts with the 0-4 scale 
  readsmelt<-melt(rarefy_reads, id.vars=c("sort","haplotype","OTU","sequences","q1","q2","q3","best_identity","phylum","Order","Family","Genus","Scientific_name","blastlowest","blastpercent"))
  rarereadsmelt1<- readsmelt[readsmelt$value>0,]# remove rows without the haplotypes present
  rarereadsmelt1<-rarereadsmelt1[!is.na(rarereadsmelt1$value),]
  rarereadsmelt1$value<-as.numeric(rarereadsmelt1$value)
  rarereadsmelt1$q1<-as.numeric(rarereadsmelt1$q1)
  rarereadsmelt1$q2<-as.numeric(rarereadsmelt1$q2)
  rarereadsmelt1$q3<-as.numeric(rarereadsmelt1$q3)
  rarereadsmelt1<- rarereadsmelt1 %>% mutate(ranking =
                       case_when(rarereadsmelt1$value <= rarereadsmelt1$q1 ~ 1, 
                                 rarereadsmelt1$value <= rarereadsmelt1$q2 & rarereadsmelt1$value > rarereadsmelt1$q1 ~ 2,
                                 rarereadsmelt1$value <= rarereadsmelt1$q3 & rarereadsmelt1$value > rarereadsmelt1$q2 ~ 3,
                                 rarereadsmelt1$value > rarereadsmelt1$q3 ~ 4))
  
  
  #convert rarefied dataset to the right format for further use
  rarereadsmelt<- as.data.frame(lapply(rarereadsmelt1, rep, rarereadsmelt1$ranking))
  rarereadsmelt$site <- metadata$Island[match(rarereadsmelt$variable,metadata$Samples)] # add island sampled from
  rarereadsmelt$site<-as.factor(rarereadsmelt$site) # make sure as factor
  rarereadsmelt$region <- metadata$Region[match(rarereadsmelt$variable,metadata$Samples)] # add hilo/kona sampled from
  rarereadsmelt$region<-as.factor(rarereadsmelt$region) # make sure as factor
  rarereadsmelt$sample <- metadata$Location[match(rarereadsmelt$variable,metadata$Samples)] # add hilo/kona sampled from
  rarereadsmelt$sample <-as.factor(rarereadsmelt$sample)
  rarefied<-rarereadsmelt
  
  ###### create a normalized (rarefied) "raw" read dataset ##########
  #convert rarefied dataset to the right format for further use
  normreads<- as.data.frame(lapply(rarereadsmelt1, rep, rarereadsmelt1$value))
  normreads$site <- metadata$Island[match(normreads$variable,metadata$Samples)] # add island sampled from
  normreads$site<-as.factor(normreads$site) # make sure as factor
  normreads$region <- metadata$Region[match(normreads$variable,metadata$Samples)] # add hilo/kona sampled from
  normreads$region<-as.factor(normreads$region) # make sure as factor
  normreads$sample <- metadata$Location[match(normreads$variable,metadata$Samples)] # add hilo/kona sampled from
  normreads$sample <-as.factor(normreads$sample)
  normalized<-normreads
  
  ########### create a presence absence dataset ############
  data_presabsen<-combined3[5:(endcolumn+1)] #just sample columns
  data_presabsen[data_presabsen>0]<-1 
  presabsen<-combined3[c(2:(nsamples+6),(nsamples+11),(nsamples+8),(nsamples+9),(nsamples+13),(nsamples+22),(nsamples+21),(nsamples+25))] 
  presabsen[4:endcolumn]<-data_presabsen
  presabsen$total_pres_absen<-rowSums(presabsen[,4:endcolumn]) # add totals
  
  #change data to correct format and keeping location
  test3<-melt(presabsen, id.vars=c("sort","haplotype","OTU","sequences","total_pres_absen","best_identity","phylum","Order","Family","Genus","Scientific_name","blastlowest","blastpercent"))
  test4<- test3[test3$value>0,] # remove rows without the haplotypes present
  test4$site <- metadata$Island[match(test4$variable,metadata$Samples)] # add island sampled from
  test4$site<-as.factor(test4$site) # make sure as factor
  test4$region <- metadata$Region[match(test4$variable,metadata$Samples)] # add hilo/kona sampled from
  test4$region<-as.factor(test4$region) # make sure as factor
  test4$sample <- metadata$Location[match(test4$variable,metadata$Samples)] # add hilo/kona sampled from
  test4$sample <-as.factor(test4$sample) # make sure as factor
  presenceabsence<-test4

  #name the datasets so you don't write over if repeating this function
  assign(  paste("final_combined", runname, sep = ""), final_combined, envir = parent.frame())
  assign(  paste("presenceabsence", runname, sep = ""), presenceabsence, envir = parent.frame())
  assign(  paste("rarefied", runname, sep = ""), rarefied,envir = parent.frame() )
  assign(  paste("normalized", runname, sep = ""), normalized,envir = parent.frame() )
  assign(  paste("combinedraw", runname, sep = ""), combinedraw, envir = parent.frame())
  assign(  paste("combined_rare", runname, sep = ""), combined_rare, envir = parent.frame())
}
#for raw reads I need to reduce some species so that the functions don't run out of memory
reducereads<-function(datafile,nsamp,alignment,runname,folder,metadata){
  otus<-unique(datafile$OTU)
  i<-1
  popgen_new<-data.frame(matrix(ncol = 22, nrow = 0))
  for (i in 1:length(otus)){
    temp2 <- datafile[datafile$OTU==otus[i],] #subet one OTU at a time
    finalrow<-nsamp+4
    reads<-temp2[5:finalrow]
    totalreads<-sum(reads)
    if(totalreads>12000){
      divideby<-totalreads/12000
      temp2[5:finalrow]<-round(temp2[5:finalrow]/divideby)
    }
    readsmelt<-melt(temp2, id.vars=c("X","sort","haplotype","OTU","sequences","best_identity","phylum","Order","Family","Genus","Scientific_name","blastlowest","blastpercent"))
    rarereadsmelt1<- readsmelt[readsmelt$value>0,]# remove rows without the haplotypes present
    rarereadsmelt1<-rarereadsmelt1[!is.na(rarereadsmelt1$value),]
    normreads<- as.data.frame(lapply(rarereadsmelt1,rep,rarereadsmelt1$value))
    normreads$site <- metadata$Island[match(normreads$variable,metadata$Samples)] # add island sampled from
    normreads$site<-as.factor(normreads$site) # make sure as factor
    normreads$region <- metadata$Region[match(normreads$variable,metadata$Samples)] # add hilo/kona sampled from
    normreads$region<-as.factor(normreads$region) # make sure as factor
    normreads$sample <- metadata$Location[match(normreads$variable,metadata$Samples)] # add hilo/kona sampled from
    normreads$sample <-as.factor(normreads$sample)
    popgen_new<-rbind(popgen_new,normreads)
  }
  #REMOVE SEQs WITH STOP CODONS
  #remove sequences with stop codons as they are errors, psuedogenes, or numts
  i<-1
  id_stopcodons<-c()
  for (i in 1:length(alignment)) {
    temp<-alignment[[i]]
    temp2<-paste(temp,collapse = "")
    temp3<-s2c(temp2)
    codons<-seqinr::translate(temp3, frame=0, numcode = 2)
    stopc<-grepl("[[:punct:]]",codons)
    if(any(stopc=="TRUE")){
      id<-names(alignment)[i] 
      id_stopcodons<-append(id_stopcodons,id)
    }
  }
  if(length(id_stopcodons)!= 0){
    for (i in 1:length(id_stopcodons)){
      popgen_new<-popgen_new[(popgen_new$haplotype != id_stopcodons[i]),]
    }}
  popgen_reduced<<-popgen_new
  write.csv(popgen_new,file = paste(folder,"/",runname,"popgen_noSC_reduced.csv",sep = ""))
}
#remove stop codons
remove_stop_codons<-function(alignment,presenceabsence,rarefied,normalized,alldata,runname,frame1,numcode1,folder){
  #remove sequences with stop codons as they are errors, psuedogenes, or numts
  i<-1
  id_stopcodons<-c()
  for (i in 1:length(alignment)) {
    temp<-alignment[[i]]
    temp2<-paste(temp,collapse = "")
    temp3<-s2c(temp2)
    codons<-seqinr::translate(temp3, frame=frame1, numcode = numcode1)
    stopc<-grepl("[[:punct:]]",codons)
    if(any(stopc=="TRUE")){
      id<-names(alignment)[i] 
      id_stopcodons<-append(id_stopcodons,id)
    }
  }
  if(length(id_stopcodons)!= 0){
  for (i in 1:length(id_stopcodons)){
    presenceabsence<-presenceabsence[(presenceabsence$haplotype != id_stopcodons[i]),]
  }
  for (i in 1:length(id_stopcodons)){
    rarefied<-rarefied[(rarefied$haplotype != id_stopcodons[i]),]
  }
  for (i in 1:length(id_stopcodons)){
    alldata<-alldata[(alldata$haplotype != id_stopcodons[i]),]
  }
  for (i in 1:length(id_stopcodons)){
    normalized<-normalized[(normalized$haplotype != id_stopcodons[i]),]
  }}

  assign(  paste("presenceabsence", runname, sep = ""), presenceabsence, envir = parent.frame())
  assign(  paste("rarefied", runname, sep = ""), rarefied, envir = parent.frame())
  assign(  paste("normalized", runname, sep = ""), normalized, envir = parent.frame())
  assign(  paste("final_combined", runname, sep = ""), alldata, envir = parent.frame())
  write.csv(alldata,file = paste(folder,"/",runname,"species_haplotype_table_filtered_noSC.csv",sep = ""))
}
# Summarize how many haplotypes and haplotype counts per species, population, and with different proxies; This also creates files needed for popgen formulas
summary_eDNA<-function(OTUvsspecies="OTU",presenceabsence,rarefied,normalized,runname,nhap,folder,haplorestrict){
  ########### Summary for presence absence dataset ##########
  otus<-unique(presenceabsence[[OTUvsspecies]])
  summary_pa<-data.frame(matrix(ncol = 18, nrow = 0))
  i<-1
  
  for (i in 1:length(otus)){
    temp <- presenceabsence[presenceabsence[[OTUvsspecies]]==otus[i],] #pull out one otu at a time
    #only keep haplotypes with at least 2 occurances
    if(haplorestrict==TRUE){
      temp3<-count(temp$haplotype)
      temp3<-temp3[(temp3$freq==1),]
      i<-1
      lowhaplo<-temp3$x
      for (i in 1:length(lowhaplo)){
        temp<-temp[(temp$haplotype != lowhaplo[i]),]
      }
    }
    total_hap<-length(unique(temp$sequences)) # get number of haplotypes
    # This is where you need to change based on what populations you use
    Oahu<-temp[temp$site=="Oahu",] #pull info for haplotypes on Oahu
    Hawaii<-temp[temp$site=="Hawaii",] #pull info for haplotypes on Big Island
    Kona<-temp[temp$region=="Kona",] #pull info for haplotypes on Kona side
    Hilo<-temp[temp$region=="Hilo",] #pull info for haplotypes on Hilo side
    Kauai<-temp[temp$site=="Kauai",]
    Maui<-temp[temp$site=="Maui",]
    # count haplotypes per population
    hap_Oahu<-length(unique(Oahu$sequences)) # get number of haplotypes on Oahu
    hap_Hawaii<-length(unique(Hawaii$sequences)) # get number of haplotypes on Big Island
    hap_Kona<-length(unique(Kona$sequences)) # get number of haplotypes on Kona side
    hap_Hilo<-length(unique(Hilo$sequences)) # get number of haplotypes on Hilo side
    hap_Kauai<-length(unique(Kauai$sequences))
    hap_Maui<-length(unique(Maui$sequences))
    # count number of samples haplotypes were found in
    ind_Oahu<-nrow(Oahu) # get number of "individuals" detected on Oahu
    ind_Hawaii<-nrow(Hawaii) # get number of "individuals" detected on Big Island
    ind_Kona<-nrow(Kona) # get number of "individuals" detected on Kona side
    ind_Hilo<-nrow(Hilo) # get number of "individuals" detected on Hilo side
    ind_Kauai<-nrow(Kauai) 
    ind_Maui<-nrow(Maui) 
    counts<-c(otus[i],temp$Scientific_name[1],temp$best_identity[1],temp$blastlowest[1],temp$phylum[1],total_hap,hap_Oahu,ind_Oahu,hap_Hawaii,ind_Hawaii,hap_Kauai,ind_Kauai,hap_Maui,ind_Maui,hap_Kona,ind_Kona,hap_Hilo,ind_Hilo) # merge data
    summary_pa<-rbind(summary_pa, counts) # add data dataframe
  }
  colnames(summary_pa)<-c("OTU","Species","Percent Matched","BlastMatch","Phylum","Total Haplotypes", "Oahu Haplotypes","Oahu Individuals","Hawaii Haplotypes","Hawaii Individuals","Kauai Haplotypes","Kauai Individuals","Maui Haplotypes","Maui Individuals","Kona Haplotypes","Kona Individuals","Hilo Haplotypes","Hilo Individuals")      
  summary_pa$`Total Haplotypes`<-as.numeric(summary_pa$`Total Haplotypes`)
  summary_pa<-summary_pa[order(summary_pa$"Total Haplotypes", decreasing = TRUE),]
  write.csv(summary_pa, file = paste(folder,"/",runname, "_summmary_presenceabsence_haplotypes.csv", sep=""))
  
  #remove samples that have low sample sizes (less than 5 samples per island)
  cols.num <- c("Oahu Individuals","Hawaii Individuals","Kona Individuals","Hilo Individuals","Maui Individuals","Kauai Individuals")
  summary_pa[cols.num] <- sapply( summary_pa[cols.num], as.numeric )
  low_samples_df2_pa<-summary_pa[(summary_pa$`Total Haplotypes`==1 | summary_pa$`Oahu Individuals`< nhap & summary_pa$`Hawaii Individuals`< nhap & summary_pa$`Kauai Individuals`< nhap & summary_pa$`Maui Individuals`< nhap),] # get list of OTUs with less than 5 individuals on either Oahu or Big Island
  low_samples_pa<-unique(low_samples_df2_pa$OTU) # get list of OTUs
  low_region_df_pa<-summary_pa[(summary_pa$`Total Haplotypes`==1 | summary_pa$`Hilo Individuals`< nhap & summary_pa$`Kona Individuals`< nhap & summary_pa$`North Oahu Individuals`< nhap & summary_pa$`South Oahu Individuals`< nhap & summary_pa$`West Oahu Individuals`< nhap & summary_pa$`East Oahu Individuals`< nhap & summary_pa$`North Kauai Individuals`< nhap & summary_pa$`South Kauai Individuals`< nhap & summary_pa$`East Kauai Individuals`< nhap),]
  low_region_pa<-unique(low_region_df_pa$OTU)
  
  ########### Summary for rarefied dataset ##########
  otusr<-unique(rarefied[[OTUvsspecies]])
  summary_r<-data.frame(matrix(ncol = 18, nrow = 0))
  i<-1
  for (i in 1:length(otus)){
    temp <- rarefied[rarefied[[OTUvsspecies]]==otus[i],] #pull out one otu at a time
    #only keep haplotypes with at least 2 occurances
    if(haplorestrict==TRUE){
      temp3<-count(temp$haplotype)
      temp3<-temp3[(temp3$freq==1),]
      i<-1
      lowhaplo<-temp3$x
      for (i in 1:length(lowhaplo)){
        temp<-temp[(temp$haplotype != lowhaplo[i]),]
      }
    }
    total_hap<-length(unique(temp$sequences)) # get number of haplotypes
    # This is where you need to change based on what populations you use
    Oahu<-temp[temp$site=="Oahu",] #pull info for haplotypes on Oahu
    Hawaii<-temp[temp$site=="Hawaii",] #pull info for haplotypes on Big Island
    Kona<-temp[temp$region=="Kona",] #pull info for haplotypes on Kona side
    Hilo<-temp[temp$region=="Hilo",] #pull info for haplotypes on Hilo side
    Kauai<-temp[temp$site=="Kauai",]
    Maui<-temp[temp$site=="Maui",]
    # count haplotypes per population
    hap_Oahu<-length(unique(Oahu$sequences)) # get number of haplotypes on Oahu
    hap_Hawaii<-length(unique(Hawaii$sequences)) # get number of haplotypes on Big Island
    hap_Kona<-length(unique(Kona$sequences)) # get number of haplotypes on Kona side
    hap_Hilo<-length(unique(Hilo$sequences)) # get number of haplotypes on Hilo side
    hap_Kauai<-length(unique(Kauai$sequences))
    hap_Maui<-length(unique(Maui$sequences))
    # count number of samples haplotypes were found in
    ind_Oahu<-nrow(Oahu) # get number of "individuals" detected on Oahu
    ind_Hawaii<-nrow(Hawaii) # get number of "individuals" detected on Big Island
    ind_Kona<-nrow(Kona) # get number of "individuals" detected on Kona side
    ind_Hilo<-nrow(Hilo) # get number of "individuals" detected on Hilo side
    ind_Kauai<-nrow(Kauai) 
    ind_Maui<-nrow(Maui) 
    counts<-c(otus[i],temp$Scientific_name[1],temp$best_identity[1],temp$blastlowest[1],temp$phylum[1],total_hap,hap_Oahu,ind_Oahu,hap_Hawaii,ind_Hawaii,hap_Kauai,ind_Kauai,hap_Maui,ind_Maui,hap_Kona,ind_Kona,hap_Hilo,ind_Hilo) # merge data
    summary_r<-rbind(summary_r, counts) # add data dataframe
  }
  colnames(summary_r)<-c("OTU","Species","Percent Matched","BlastMatch","Phylum","Total Haplotypes", "Oahu Haplotypes","Oahu Individuals","Hawaii Haplotypes","Hawaii Individuals","Kauai Haplotypes","Kauai Individuals","Maui Haplotypes","Maui Individuals","Kona Haplotypes","Kona Individuals","Hilo Haplotypes","Hilo Individuals")      
  summary_r$`Total Haplotypes`<-as.numeric(summary_r$`Total Haplotypes`)
  summary_r<-summary_r[order(summary_r$"Total Haplotypes", decreasing = TRUE),]
  write.csv(summary_r, file = paste(folder,"/",runname, "_summmary_rarefied_haplotypes.csv", sep=""))
  
  #remove samples that have low sample sizes (less than 5 samples per island)
  cols.num <- c("Oahu Individuals","Hawaii Individuals","Kona Individuals","Hilo Individuals","Maui Individuals","Kauai Individuals")
  summary_r[cols.num] <- sapply( summary_r[cols.num], as.numeric )
  low_samples_df2_r<-summary_r[(summary_r$`Total Haplotypes`==1 | summary_r$`Oahu Individuals`< nhap & summary_r$`Hawaii Individuals`< nhap & summary_r$`Kauai Individuals`< nhap & summary_r$`Maui Individuals`< nhap),] # get list of OTUs with less than 5 individuals on either Oahu or Big Island
  low_samples_r<-unique(low_samples_df2_r$OTU) # get list of OTUs
  low_region_df_r<-summary_r[(summary_r$`Total Haplotypes`==1 | summary_r$`Hilo Individuals`< nhap | summary_r$`Kona Individuals`< nhap),]
  low_region_r<-unique(low_region_df_r$OTU)
  
  ########### Summary for normalized reads ##########
  otus<-unique(normalized[[OTUvsspecies]])
  summary_n<-data.frame(matrix(ncol = 18, nrow = 0))
  i<-1
  
  for (i in 1:length(otus)){
    temp <- normalized[normalized[[OTUvsspecies]]==otus[i],] #pull out one otu at a time
    #only keep haplotypes with at least 2 occurances
    if(haplorestrict==TRUE){
      temp3<-count(temp$haplotype)
      temp3<-temp3[(temp3$freq==1),]
      i<-1
      lowhaplo<-temp3$x
      for (i in 1:length(lowhaplo)){
        temp<-temp[(temp$haplotype != lowhaplo[i]),]
      }
    }
    total_hap<-length(unique(temp$sequences)) # get number of haplotypes
    # This is where you need to change based on what populations you use
    Oahu<-temp[temp$site=="Oahu",] #pull info for haplotypes on Oahu
    Hawaii<-temp[temp$site=="Hawaii",] #pull info for haplotypes on Big Island
    Kona<-temp[temp$region=="Kona",] #pull info for haplotypes on Kona side
    Hilo<-temp[temp$region=="Hilo",] #pull info for haplotypes on Hilo side
    Kauai<-temp[temp$site=="Kauai",]
    Maui<-temp[temp$site=="Maui",]
    # count haplotypes per population
    hap_Oahu<-length(unique(Oahu$sequences)) # get number of haplotypes on Oahu
    hap_Hawaii<-length(unique(Hawaii$sequences)) # get number of haplotypes on Big Island
    hap_Kona<-length(unique(Kona$sequences)) # get number of haplotypes on Kona side
    hap_Hilo<-length(unique(Hilo$sequences)) # get number of haplotypes on Hilo side
    hap_Kauai<-length(unique(Kauai$sequences))
    hap_Maui<-length(unique(Maui$sequences))
    # count number of samples haplotypes were found in
    ind_Oahu<-nrow(Oahu) # get number of "individuals" detected on Oahu
    ind_Hawaii<-nrow(Hawaii) # get number of "individuals" detected on Big Island
    ind_Kona<-nrow(Kona) # get number of "individuals" detected on Kona side
    ind_Hilo<-nrow(Hilo) # get number of "individuals" detected on Hilo side
    ind_Kauai<-nrow(Kauai) 
    ind_Maui<-nrow(Maui) 
    counts<-c(otus[i],temp$Scientific_name[1],temp$best_identity[1],temp$blastlowest[1],temp$phylum[1],total_hap,hap_Oahu,ind_Oahu,hap_Hawaii,ind_Hawaii,hap_Kauai,ind_Kauai,hap_Maui,ind_Maui,hap_Kona,ind_Kona,hap_Hilo,ind_Hilo) # merge data
    summary_n<-rbind(summary_n, counts) # add data dataframe
  }
  colnames(summary_n)<-c("OTU","Species","Percent Matched","BlastMatch","Phylum","Total Haplotypes", "Oahu Haplotypes","Oahu Individuals","Hawaii Haplotypes","Hawaii Individuals","Kauai Haplotypes","Kauai Individuals","Maui Haplotypes","Maui Individuals","Kona Haplotypes","Kona Individuals","Hilo Haplotypes","Hilo Individuals")      
  summary_n$`Total Haplotypes`<-as.numeric(summary_n$`Total Haplotypes`)
  summary_n<-summary_n[order(summary_n$"Total Haplotypes", decreasing = TRUE),]
  write.csv(summary_n, file = paste(folder,"/",runname, "_summmary_normalized_haplotypes.csv", sep=""))
  
  #remove samples that have low sample sizes (less than 5 samples per island)
  cols.num <- c("Oahu Individuals","Hawaii Individuals","Kona Individuals","Hilo Individuals","Maui Individuals","Kauai Individuals")
  summary_n[cols.num] <- sapply( summary_n[cols.num], as.numeric )
  low_samples_df2_n<-summary_n[(summary_n$`Total Haplotypes`==1 | summary_n$`Oahu Individuals`< nhap & summary_n$`Hawaii Individuals`< nhap & summary_n$`Kauai Individuals`< nhap & summary_n$`Maui Individuals`< nhap),] # get list of OTUs with less than 5 individuals on either Oahu or Big Island
  low_samples_n<-unique(low_samples_df2_n$OTU) # get list of OTUs
  low_region_df_n<-summary_n[(summary_n$`Total Haplotypes`==1 | summary_n$`Hilo Individuals`< nhap & summary_n$`Kona Individuals`< nhap & summary_n$`North Oahu Individuals`< nhap & summary_n$`South Oahu Individuals`< nhap & summary_n$`West Oahu Individuals`< nhap & summary_n$`East Oahu Individuals`< nhap & summary_n$`North Kauai Individuals`< nhap & summary_n$`South Kauai Individuals`< nhap & summary_n$`East Kauai Individuals`< nhap),]
  low_region_n<-unique(low_region_df_n$OTU)
  
  popgen_pa<-presenceabsence
  if(length(low_samples_pa)>0){
    i<-1
  for (i in 1:length(low_samples_pa)){
    popgen_pa<-popgen_pa[(popgen_pa$OTU != low_samples_pa[i]),]
  }}
  popgen_pa_region<-presenceabsence
  if(length(low_region_pa)>0){
    i<-1
  for (i in 1:length(low_region_pa)){
    popgen_pa_region<-popgen_pa_region[(popgen_pa_region$OTU != low_region_pa[i]),]
  }}
  popgen_r<-rarefied
  if(length(low_samples_r)>0){
    i<-1
  for (i in 1:length(low_samples_r)){
    popgen_r<-popgen_r[(popgen_r$OTU != low_samples_r[i]),]
  }}
  popgen_r_region<-rarefied
  if(length(low_region_r)>0){
    i<-1
  for (i in 1:length(low_region_r)){
    popgen_r_region<-popgen_r_region[(popgen_r_region$OTU != low_region_r[i]),]
  }}
  popgen_n<-normalized
  if(length(low_samples_n)>0){
    i<-1
  for (i in 1:length(low_samples_n)){
    popgen_n<-popgen_n[(popgen_n$OTU != low_samples_n[i]),]
  }}
  popgen_n_region<-normalized
  if(length(low_region_n)>0){
    i<-1
  for (i in 1:length(low_region_n)){
    popgen_n_region<-popgen_n_region[(popgen_n_region$OTU != low_region_n[i]),]
  }}
  #name the datasets so you don't write over if repeating this function
  assign(  paste("popgen_pa", runname, sep = ""), popgen_pa, envir = parent.frame())
  write.csv(popgen_pa, file = paste(folder,"/",runname, "_popgen_pa.csv", sep=""))
  assign(  paste("popgen_r", runname, sep = ""), popgen_r, envir = parent.frame())
  write.csv(popgen_r, file = paste(folder,"/",runname, "_popgen_r.csv", sep=""))
  assign(  paste("popgen_n", runname, sep = ""), popgen_n, envir = parent.frame())
  write.csv(popgen_n, file = paste(folder,"/",runname, "_popgen_n.csv", sep=""))
  assign(  paste("popgen_pa_region", runname, sep = ""), popgen_pa_region, envir = parent.frame())
  assign(  paste("popgen_r_region", runname, sep = ""), popgen_r_region, envir = parent.frame())
  assign(  paste("popgen_n_region", runname, sep = ""), popgen_n_region, envir = parent.frame())
}
# to create a large fasta file
fastafilealldata<-function(metadata,table,obitools,blast,lowestmatch,OTUvsspecies,folder,nblanks,nsamples,runname,removefirstbp,frame1,numcode1){
  ########### reformatting ############
  #remove blanks and create a subset file to see what had contamination
  if (nblanks > 0) {
    blanks<-metadata$Samples[metadata$Side=="Blank"]
    if (nblanks > 1) {contam<-table[rowSums(table[,c(blanks[1:nblanks])]!=0)>0,] #keep rows with something in the blanks
    } else if (nblanks==1) {contam<-table[table[blanks]>0,]}
    table1 <- table[,!(names(table) %in% blanks)] #remove blanks from table
  }
  
  # remove the line that says how many reads were removed
  lastrow<-(nrow(table1)-1)
  table1<-table1[1:lastrow,] 
  
  #get a list of total reads per sample which will later be used to rarefaction
  if (nblanks > 0) {
    samplecolumns<-metadata$Samples[metadata$Side!="Blank"]
  } else {samplecolumns<-metadata$Samples}
  sample_totals<-colSums(table1[,c(samplecolumns[1:nsamples])])
  
  
  #remove first nucleotide (I decided to do this because it seems there is an error that some of the sequences have a G instead of A which is never seen in the tissue derived samples of longer sequence length)
  if(removefirstbp == TRUE) {
    table1$sequences <- as.character(table1$sequences)
    table1$sequences <- gsub("^.{1}", "", table1$sequences)
  }
  
  #remove sequences of different lengths
  otus <- unique(table1$OTU) #load a list of OTU numbers to refer to later
  i<-1
  id_all<-c()
  for (i in 1:length(otus)){
    temp <- table1[table1$OTU==otus[i],]
    x <- as.list(temp$sequences)
    bp<-nchar(x)
    temp$bp<-bp
    temp<-temp[(temp$bp != bp[1]),]
    id<-temp$haplotype
    id_all<-append(id_all,id)
  }
  if (length(id_all) != 0){
    for (i in 1:length(id_all)){
      table1<-table1[(table1$haplotype != id_all[i]),]
    }}
  
  #remove samples that had all 0 reads
  table1<-table1[rowSums(table1[,c(samplecolumns[1:nsamples])]!=0)>0,] #remove rows with all zeros
  
  # create combined file
  table1$id<-paste(table1$OTU,table1$haplotype,sep = "_")
  obitools<-obitools[1:lastrow,c(1,3:4,7,9:13,15:17)]
  colnames(obitools)<-c("id","best_identity","best_match","Family","Genus","Match_Count","Order","Rank","Scientific_name","SpeciesList","SpeciesName","TaxID")
  blast<-blast[c(1,4:8)]
  colnames(blast)<-c("id","blastspecies","blastID","evalue","querycover","blastpercent") # rename so not V1,V2...
  combined_first<-merge(obitools,blast, by="id",all.x=TRUE)
  combined<-merge(table1,combined_first,by="id",all.x=TRUE)
  combinedraw<-combined
  
  ########### Remove non target taxa ############
  #remove non target taxa including non marine metazoans: humans, bacteria, algae, chicken, dog, ...
  combined<-combined[!grepl("9606", combined$blastID),] # remove any rows with ...
  combined<-combined[!grepl("Bilateria|Eumetazoa|Eukaryota|Bacteria|root|Bacilli|Philasterida|Flustrina|Tetrahymena|Lactobacillales|Ciliophora", combined$Scientific_name),]
  combined<-combined[!grepl("Galliformes", combined$Order),]
  combined<-combined[!grepl("Canidae|Bovidae|Suidae|Muridae", combined$Family),]
  combined<-combined[!grepl("strain",combined$Rank),]
  combined<-combined[!grepl("33317",combined$TaxID),] 
  
  ########### Remove haplotypes that didn't match too well ###########
  combined$blastpercent<-as.numeric(combined$blastpercent) # set percentages to numerics so we can use the following code
  combined3<-subset(combined, combined$best_identity>=lowestmatch & combined$blastpercent>=(lowestmatch*100))
  
  # remove the same ones from the table file not just the combined file
  ids<-combined3$id
  idremove<-unique(table1$id)
  for (i in 1:length(ids)) { # get a list of OTUs to remove
    idremove<-idremove[(idremove != ids[i])];
  }
  if (length(idremove)!=0){
    for (i in 1:length(idremove)) { # remove non target taxa from table
      table1<-table1[(table1$id != idremove[i]),]
    }}
  
  #remove stop codons
  seqs1<-as.list(combined3$sequences)
  names1<- as.list(combined3$haplotype)
  write.fasta(sequences=seqs1,names=names1,file = paste(folder,"/",runname,"haplotypes.fasta",sep = ""))
  alignment<-read.fasta(file = paste(folder,"/",runname,"haplotypes.fasta",sep = ""))
  i<-1
  id_stopcodons<-c()
  for (i in 1:length(alignment)) {
    temp<-alignment[[i]]
    temp2<-paste(temp,collapse = "")
    temp3<-s2c(temp2)
    codons<-seqinr::translate(temp3, frame=frame1, numcode = numcode1)
    stopc<-grepl("[[:punct:]]",codons)
    if(any(stopc=="TRUE")){
      id<-names(alignment)[i] 
      id_stopcodons<-append(id_stopcodons,id)
    }
  }
  if(length(id_stopcodons)!= 0){
    for (i in 1:length(id_stopcodons)){
      combined3<-combined3[(combined3$haplotype != id_stopcodons[i]),]
    }
  }
  #sum up reads to add to fasta name so I know how abundant
  combined3$totalabundance<-rowSums(combined3[,c(samplecolumns[1:nsamples])])
  
  #write fasta file per species
  folder_name2<-paste(folder,"/fastafiles_",runname,"_singletonsincluded/",sep = "")
  dir.create(folder_name2)
  otus<-unique(combined3$OTU)
  #set unique number so that each name is unique
  i<-1
  for (i in 1:length(otus)){
    temp <- combined3[combined3$OTU==otus[i],]
    #set sequences
    seqs_fasta<-as.list(temp$sequences)
    #set names ex: >OTU1 haplo_1 abundance:
    names_fasta<- as.list(paste(temp$OTU," ",temp$haplotype, " abundance: ",temp$totalabundance, sep=""))
    #write as fasta files
    write.fasta(sequences=seqs_fasta,names=names_fasta,file.out=paste(folder_name2,temp$Scientific_name[1],"_",temp$best_identity[1],"_",temp$blastspecies[1],"_",temp$blastpercent[1],"_",temp$OTU[1],"_",runname,".fasta",sep = ""))
  }
}
# to create fasta files per species (these were used for running other programs like arlequin and DNASP)
fastafilesperspecies<-function(datafile,popcolname,runname,folder){
  folder_name2<-paste(folder,"/fastafiles",runname,"/",sep = "")
  dir.create(folder_name2)
  otus<-unique(datafile$OTU)
  #set unique number so that each name is unique
  datafile$ID<-row.names(datafile)
  i<-1
  for (i in 1:length(otus)){
    temp <- datafile[datafile$OTU==otus[i],]
    #set sequences
    seqs_fasta<-as.list(temp$sequences)
    #set names ex: >OTU1 population:Oahu haplo_1
    names_fasta<- as.list(paste(temp$OTU,"_",temp$ID," population:",temp[[popcolname]]," ",temp$haplotype, sep=""))
    #write as fasta files
    write.fasta(sequences=seqs_fasta,names=names_fasta,file.out=paste(folder_name2,temp$Scientific_name[1],"_",temp$best_identity[1],"_",temp$blastlowest[1],"_",temp$OTU[1],"_",runname,".fasta",sep = ""))
  }
}
# to create fasta files per species and per population (these were used for running other programs like arlequin and DNASP)
fastafilesperspeciesperpopulation<-function(datafile,popcolname,runname,folder){
  folder_name2<-paste(folder,"/fastafiles",runname,"/",sep = "")
  dir.create(folder_name2)
  otus<-unique(datafile$OTU)
  #set unique number so that each name is unique
  datafile$ID<-row.names(datafile)
  i<-1
  for (i in 1:length(otus)){
    temp <- datafile[datafile$OTU==otus[i],]
    pops<-unique(temp[[popcolname]])
    for (i in 1:length(pops)){
      temp2 <- temp[temp[[popcolname]]==pops[i],]
      #set sequences
      seqs_fasta<-as.list(temp2$sequences)
      #set names ex: >OTU1 population:Oahu haplo_1
      names_fasta<- as.list(paste(temp2$OTU,"_",temp2$ID," population:",temp2[[popcolname]]," ",temp2$haplotype, sep=""))
      #write as fasta files
      write.fasta(sequences=seqs_fasta,names=names_fasta,file.out=paste(folder_name2,temp2$Scientific_name[1],"_",temp2$best_identity[1],"_",temp2$blastlowest[1],"_",temp2$OTU[1],"_",temp2[[popcolname]][1],"_",runname,".fasta",sep = ""))
    }
  }
}
#converts to nexus files for DNASP
fastatonexus<-function(folderinput,folderoutput,runname){
  folder_name2<-paste(folderoutput,"/nexusfiles_",runname,"/",sep = "")
  dir.create(folder_name2)
  files<-list.files(folderinput)
  i<-1
  for(i in 1:length(files)){
    fastaspecies<-read.fasta(file= paste0(folderinput,"/",files[i]))
    write.nexus.data(fastaspecies, file=paste(folder_name2,files[i],"_nexus.nex",sep = ""))
  }
}

##### functions for population genetics and genetic diversity
# for pairwise phist
eDNA_haplotypes_pairphi<-function(OTUvsspecies="OTU",folder,datafile,runname,order,numperm,nind,popcolname,haplorestrict){
  folder_name2<-paste(folder,"/PhiStats",runname,"/",sep = "")
  dir.create(folder_name2)
  otus<-unique(datafile[[OTUvsspecies]])
  i<-1
  for (i in 1:length(otus)){
    temp2 <- datafile[datafile[[OTUvsspecies]]==otus[i],] #subet one OTU at a time
    #identify any populations with less than nind individuals
    populations<-popcolname
    ind<-data.frame(table(temp2[[popcolname]])) # count number of occurances for each of the populations
    colnames(ind)<-c("site","ind") # rename the columns
    ind$ind<-as.numeric(ind$ind) # make number of individuals into a numeric
    lowpop<-as.vector(ind$site[ind$ind<nind])
    #remove any populations with less than nind individuals
    if(length(lowpop)!=0){
    i<-1
    for (i in 1:length(lowpop)){
      temp2<-temp2[(temp2[[popcolname]] != lowpop[i]),]
    }}
    #only keep haplotypes with at least 2 occurances
    if(haplorestrict==TRUE){
      temp3<-count(temp2$haplotype)
      temp3<-temp3[(temp3$freq==1),]
      i<-1
      lowhaplo<-temp3$x
      for (i in 1:length(lowhaplo)){
        temp2<-temp2[(temp2$haplotype != lowhaplo[i]),]
      }
    }
    if(length(unique(temp2[[popcolname]]))>1){
    #if subsample=TRUE then subsample if one population is more than double the other
    #if(subsample==TRUE){
      
    #}
    #convert DNA seqs to format for pairphist
    y <- as.list(temp2$sequences)
    y <- as.DNAbin(ape::as.alignment(y))
    z <- as.dna(y)
    #identify populations for testing
    Island<-temp2[[popcolname]]
    #calculate phist and p-values
    pst<-pairPhiST(z, Island, nperm=numperm, negatives=TRUE, showprogbar=FALSE) #calculates distance in code
    #How to create a prettier table
    #remove pops that it did not occur in
    i_order<-order
    pops<-unique(temp2[[popcolname]])
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
    write.csv(phist_reordered, file = paste(folder_name2, temp2$Scientific_name[1],"_",temp2$best_identity[1],"_",temp2$blastlowest[1],"_",temp2$OTU[1],"_",runname, "_phistat.csv", sep=""))
      } }
}
# for sequence diversity statistics
statstable<-function(folder,datafile,runname,popcolname){
  otus<-unique(datafile$OTU)
  summary_stats<-data.frame(matrix(ncol = 15, nrow = 0))
  i<-1
  for (i in 1:length(otus)){
  temp <- datafile[datafile$OTU==otus[i],]
  pops<-unique(temp[[popcolname]])
  x <- as.list(temp$sequences)
  x <- as.DNAbin(ape::as.alignment(x))
  h <- pegas::haplotype(x)
  N<-nrow(temp)
  H<-length(unique(temp$sequences))
  nd<-nuc.div(h,variance=TRUE) 
  nd_var<-paste(round(nd[1],digits = 4),format(round(nd[2],digits = 6),scientific=FALSE),sep = " ± ")
  hd<-hap.div(h,variance=TRUE)
  hd_var<-paste(round(hd[1],digits = 3),round(hd[2],digits = 3),sep = " ± ")
  tajima<-tajima.test(x) 
  population<-"overall"
  m<-x <- as.matrix(x)
  if(H>1){
    theta<-1
    B<-1000
    if (is.list(x)) x <- as.matrix(x)
    n <- dim(x)[1]
    k <- mean(dist.dna(x, "N"))
    ss <- seg.sites(x)
    U <- numeric(n)
    for (i in 1:n) for (j in ss)
      if (all(x[i, j, drop = TRUE] != x[-i, j])) U[i] <- U[i] + 1
    U <- (U - k/2)^2
    R2.obs <- sqrt(sum(U)/n)/length(ss)
    B<-1000
    R <- numeric(B)
    for (b in 1:B) {
      tr <- rcoal(n, rep("", n))
      tr$edge.length <- rpois(2*n - 2, theta * tr$edge.length)
      d <- cophenetic(tr)
      k <- mean(d)
      U <- tr$edge.length[tr$edge[, 2] <= n]
      U <- (U - k/2)^2
      R[b] <- sqrt(sum(U)/n)/sum(tr$edge.length)
    }
    R <- na.omit(R)
    R2_value<-round(R2.obs,digits = 3)
    R2_P = sum(R < R2.obs)/length(R)
    R2_P = round(R2_P,digits = 3)
  }else{
    R2_value<-"NA"
    R2_p<-"NA"
  }
  basicstats<-c(temp$Scientific_name[1],temp$blastlowest[1],temp$OTU[1],population,N,H,hd_var,nd_var,round(tajima$D,digits = 3),R2_value,tajima$Pval.normal,R2_P,tajima$Pval.beta,hd[1],nd[1])
  summary_stats<-rbind(summary_stats,basicstats)
  i<-1
  for (i in 1:length(pops)){
    temp2 <- temp[temp[[popcolname]]==pops[i],]
    x <- as.list(temp2$sequences)
    x <- as.DNAbin(ape::as.alignment(x))
    h <- pegas::haplotype(x)
    N<-nrow(temp2)
    H<-length(unique(temp2$sequences))
    nd<-nuc.div(h,variance=TRUE) 
    nd_var<-paste(round(nd[1],digits = 4),format(round(nd[2],digits = 6),scientific=F),sep = " ± ")
    hd<-hap.div(h,variance=TRUE)
    hd_var<-paste(round(hd[1],digits = 3),round(hd[2],digits = 3),sep = " ± ")
    tajima<-tajima.test(x) 
    population<-as.character(temp2[[popcolname]][1])
    m<-x <- as.matrix(x)
    if(H>1){
      theta<-1
      B<-1000
      if (is.list(x)) x <- as.matrix(x)
      n <- dim(x)[1]
      k <- mean(dist.dna(x, "N"))
      ss <- seg.sites(x)
      U <- numeric(n)
      for (i in 1:n) for (j in ss)
        if (all(x[i, j, drop = TRUE] != x[-i, j])) U[i] <- U[i] + 1
      U <- (U - k/2)^2
      R2.obs <- sqrt(sum(U)/n)/length(ss)
      B<-1000
      R <- numeric(B)
      for (b in 1:B) {
        tr <- rcoal(n, rep("", n))
        tr$edge.length <- rpois(2*n - 2, theta * tr$edge.length)
        d <- cophenetic(tr)
        k <- mean(d)
        U <- tr$edge.length[tr$edge[, 2] <= n]
        U <- (U - k/2)^2
        R[b] <- sqrt(sum(U)/n)/sum(tr$edge.length)
      }
      R <- na.omit(R)
      R2_value<-round(R2.obs,digits = 3)
      R2_P = sum(R < R2.obs)/length(R)
      R2_P = round(R2_P,digits = 3)
    }else{
      R2_value<-"NA"
      R2_p<-"NA"
    }
    basicstats<-c(temp2$Scientific_name[1],temp2$blastlowest[1],temp2$OTU[1],population,N,H,hd_var,nd_var,round(tajima$D,digits = 3),R2_value,tajima$Pval.normal,R2_P,tajima$Pval.beta,hd[1],nd[1])
    summary_stats<-rbind(summary_stats,basicstats)
  }
  colnames(summary_stats)<-c("Species","Blast Species","OTU","Population","Sample Size","Number of Haplotypes","hap div var","nuc div var","Tajima's D","R2","Tajima's D p-value","R2 p-value","Tajima's D p-value beta","haplotype diversity","nucleotide diversity")
  }
  write.csv(summary_stats, file = paste(folder,"/",runname,"_basic_seq_stats.csv", sep=""))
}
# for mismatch graphs
MMD_all<-function(folder,table,runname,popcolname){
  folder_name<-paste(folder,"MMD_",runname,"/",sep = "")
  dir.create(folder_name)
  otus<-unique(table$OTU)
  for (i in 1:length(otus)){
    temp <- table[table$OTU==otus[i],]
    y <- as.list(temp$sequences)
    y <- as.DNAbin(ape::as.alignment(y))
    d <- dist.dna(y, "N")
    maxdist<-max(d)
    if (maxdist > 0){
    if(maxdist < 4) {maxdist<-4}
    lcol = c("lightblue", "darkred")
    lty = c(1, 1)
    jpeg(file=paste(folder_name,temp$Scientific_name[1],"_",temp$best_identity[1],"_",temp$blastlowest[1],"_",temp$OTU[1],"_",runname,"_all_islands_MMD.jpeg",sep = ""), width=11, height=10, units = "in", res=400)
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
    }
    pops<-unique(temp[[popcolname]])
    for (i in 1:length(pops)){
      temppop <- temp[temp[[popcolname]]==pops[i],]
      y <- as.list(temppop$sequences)
      y <- as.DNAbin(ape::as.alignment(y))
      d <- dist.dna(y, "N")
      maxdist<-max(d)
      if(maxdist>0){
      if(maxdist < 4) {maxdist<-4}
      lcol = c("lightblue", "darkred")
      lty = c(1, 1)
      jpeg(file=paste(folder_name,temppop$Scientific_name[1],"_",temppop$best_identity[1],"_",temppop$blastlowest[1],"_",temppop$OTU[1],"_",runname,"_",temppop[[popcolname]][1],"_MMD.jpeg",sep = ""), width=11, height=10, units = "in", res=400)
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
#makes individual phist results into one nice table for analysis
convert_phist_to_table<-function(folderinput,folderoutput,runname){
  files<-list.files(folderinput)
  phisttable<-data.frame(matrix(ncol = 10, nrow = 0))
  i<-1
  for(i in 1:length(files)){
    #load a single phist table
    #grab phist values
    phistspecies<-read.csv(file= paste0(folderinput,files[i]), stringsAsFactors = F)
    phistspecies2<-phistspecies
    phistspecies2[upper.tri(phistspecies2)]<-"NA"
    phistspecies2<-melt(phistspecies2, id.vars="X")
    phistspecies2<-phistspecies2[phistspecies2$value!="NA",]
    phistspecies2$variable<-gsub("\\.", " ", phistspecies2$variable)
    #grab p-values
    phistspecies3<-phistspecies
    rownames(phistspecies3)<-phistspecies3$X
    phistspecies3<-phistspecies3[,-1]
    phistspecies3[lower.tri(phistspecies3)]<-"NA"
    phistspecies3$X<-rownames(phistspecies3)
    phistspecies3<-melt(phistspecies3, id.vars="X")
    phistspecies3<-na.omit(phistspecies3)
    phistspecies3<-phistspecies3[phistspecies3$value!="NA",]
    phistspecies3$variable<-gsub("\\.", " ", phistspecies3$variable)
    #correct p-values FDR
    phistspecies3$fdr<-p.adjust(p=phistspecies3$value, method = "fdr")
    #correct p-values bonforroni
    phistspecies3$bon<-p.adjust(p=phistspecies3$value, method = "bonferroni")
    #merge together
    colnames(phistspecies3)<-c("variable","X","p value","fdr p","bon p")
    phistall <- merge(phistspecies2,phistspecies3,by=c("X","variable"), all = TRUE)
    Species<-files[i]
    phistall<-cbind(Species,phistall)
    #if significant
    phistall$sigp <- with(phistall, ifelse(phistall$`p value` < 0.05, "yes", "no"))
    phistall$sigfdr <- with(phistall, ifelse(phistall$`fdr p` < 0.05, "yes", "no"))
    phistall$sigbon <- with(phistall, ifelse(phistall$`bon p`< 0.05, "yes", "no"))
    #bind to table
    phisttable<-rbind(phisttable,phistall)
      }
  colnames(phisttable)<-c("Species","Population 1","Population 2","PhiST value","p-value","fdr p-value","bonferroni p-value","sig p-value","sig fdr","sig bon")
  write.csv(phisttable, file = paste(folderoutput,"/",runname,"_phist_table_manyspecies.csv", sep=""))
}
# this is not used in the paper but was investigated
nmdsbyspecies<-function(folder,datafile,runname,popcolname,nsamples){
  foldername<-paste(folder,runname,sep = "")
  dir.create(foldername)
  otus<-datafile$OTU
  for (i in 1:length(otus)) {
  table1<-datafile[datafile$OTU=="OTU_515",]
  test<-table1[6:(5+nsamples)] # get samples
  #remove samples that don't have any detections
  columnsums<-colSums(test)
  df3<-data.frame(colnames(test),columnsums)
  emptysamples<-as.list(df3$colnames.test.[df3$columnsums==0])
  test<-select(test, -all_of(unlist(emptysamples)))
  rnam<-paste(table1$OTU,table1$blastspecies,table1$haplotype,sep="_")
  rownames(test)<-rnam # name rows on OTU + haplotype # because needs to be unique
  test<-as.matrix(test) # format data for phyloseq
  OTU <- otu_table(test, taxa_are_rows = TRUE) # format data for phyloseq
  combined_raw<-table1
  combined_raw[is.na(combined_raw)] <- "None" # some rows are blank so fill with info
  taxmat<-data.frame(Domain = combined_raw$Order, 
                     Phylum = combined_raw$Order,
                     Class = combined_raw$Order,
                     Order = combined_raw$Order,
                     Family = combined_raw$Family,
                     Genus = combined_raw$Genus,
                     Species = combined_raw$Scientific_name) # create tax file for phyloseq object
  rownames(taxmat) <- rownames(test) # name rows the same as OTU file
  taxmat<-as.matrix(taxmat) # format for phyloseq object
  TAX <- tax_table(taxmat) # format for phyloseq object
  physeq <- phyloseq(OTU, TAX) # create phyloseq object
  pmetadata<-metadata[!(metadata$Island=="Blank"),] # remove blanks
  emptysamples2<-unlist(emptysamples)
  if(length(emptysamples2)>0){
    i<-1
    for (i in 1:length(emptysamples2)){
      pmetadata<-pmetadata[(pmetadata$Samples != emptysamples2[i]),]
    }}
  pmetadata$Year<-as.factor(pmetadata$Year)
  rownames(pmetadata)<- sample_names(physeq) # name the rows
  pmetadata<-sample_data(pmetadata) # convert to phyloseq format
  random_tree <- rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq)) # not sure why this is needed but it is part of the phyloseq object
  physeq_raw <- merge_phyloseq(physeq, pmetadata, random_tree) # create file phyloseq object
  raw_ord<- ordinate(physeq_raw, "NMDS", "bray") # an ordination is needed to plot the nmds
  col.2<-c("darkgreen","lightblue") # pick colors for plot
  otuname<-table1$OTU[1]
  speciesname<-table1$Scientific_name[1]
  blastname<-table1$blastspecies[1]
  jpeg(paste(foldername,"/",otuname,"_",speciesname,"_",blastname,"_pca_island.jpg",sep = ""), width=12, height=10, units= "in", res=600)
  nmdsplot_IbyY_nf <- plot_ordination(physeq_raw, raw_ord, "samples", color="Island",shape = "Island") + 
    geom_polygon(aes(alpha = 0.8,fill=Island)) + geom_point(size=5) + ggtitle("NMDS of raw reads")
  print(nmdsplot_IbyY_nf)
  dev.off()
  jpeg(paste(foldername,"/",otuname,"_",speciesname,"_",blastname,"_pca_location.jpg",sep = ""), width=12, height=10, units= "in", res=600)
  nmdsplot_IbyY_nf2 <- plot_ordination(physeq_raw, raw_ord, "samples", color="Location") + 
    facet_wrap(~Island) + geom_polygon(aes(alpha = 0.8,fill=Location)) + geom_point(size=5) + ggtitle("NMDS of raw reads")
  print(nmdsplot_IbyY_nf2)
  dev.off()
  }
}
# not used in the paper but was investigated
pairwiseFst<-function(folder, datafile, runname, popcolname,order1,numperm,nind,haplorestrict){
  folder_name2<-paste(folder,"/PairwiseFst",runname,"/",sep = "")
  dir.create(folder_name2)
  otus<-unique(datafile$OTU)
  i<-1
  for (i in 1:length(otus)){
    temp2 <- datafile[datafile$OTU==otus[i],] #subet one OTU at a time
    #identify any populations with less than nind individuals
    populations<-popcolname
    ind<-data.frame(table(temp2[[popcolname]])) # count number of occurances for each of the populations
    colnames(ind)<-c("site","ind") # rename the columns
    ind$ind<-as.numeric(ind$ind) # make number of individuals into a numeric
    lowpop<-as.vector(ind$site[ind$ind<nind])
    #remove any populations with less than nind individuals
    if(length(lowpop)!=0){
      i<-1
      for (i in 1:length(lowpop)){
        temp2<-temp2[(temp2[[popcolname]] != lowpop[i]),]
      }}
    #only keep haplotypes with at least 2 occurances
    if(haplorestrict==TRUE){
      temp3<-count(temp2$haplotype)
      temp3<-temp3[(temp3$freq==1),]
      i<-1
      lowhaplo<-temp3$x
      for (i in 1:length(lowhaplo)){
        temp2<-temp2[(temp2$haplotype != lowhaplo[i]),]
      }
    }
    if(length(unique(temp2[[popcolname]]))>1){
      #if subsample=TRUE then subsample if one population is more than double the other
      #if(subsample==TRUE){
      #}
      #identify populations for testing
      Island<-temp2[[popcolname]]
      #convert to dataframe
      dfloci<-data.frame(population = Island,
                         seqs = temp2$sequences)
      #remove pops that it did not occur in
      i_order<-order1
      pops<-unique(temp2[[popcolname]])
      nopops<-setdiff(i_order, pops)
      if(length(nopops)!=0){
        i<-1
        for (i in 1:length(nopops)){
          i_order<-i_order[(i_order != nopops[i])]
        }}
      #convert to genind
      geninddata<-df2genind(dfloci[-1],pop = Island,sep = NULL,ploidy = 1)
      #convert to heirfstat
      #hdata<-genind2hierfstat(geninddata,pop=NULL)
      # run heirfstat
      pairwise.fst(dat = geninddata, pop=Island)
      resultsWC84<-pairwise.WCfst(locus_hierfstat,diploid=FALSE)
      resultsWC84<-as.matrix(resultsWC84) #convert to matrix
      resultsWC84[upper.tri(resultsWC84)] <- t(resultsWC84)[upper.tri(resultsWC84)] #copy so symettrical 
      resultsWC84_reordered <- reorder_mat(mat = resultsWC84, order = i_order) #reorder
      resultsWC84_reordered[upper.tri(resultsWC84_reordered)]<-NA # clear top
      #get p-values
      WC84_observed<-melt(resultsWC84_reordered) # melt so we can make a table of all perms
      WC84_observed$pairwisecomp<-paste(WC84_observed$Var1,WC84_observed$Var2,sep = " and ") # combine 2 populations into one column
      WC84_observed<-na.omit(WC84_observed) # remove duplicate info
      WC84_observed<-WC84_observed[WC84_observed$Var1!=WC84_observed$Var2,] # remove comparisons to itself 
      #Now start permutations
        #permuted_fsts<-data.frame(pairwisecomp = WC84_observed$pairwisecomp) #make dataset to fill with permutations
      #for (perm in seq_len(numperm)) {
        # Show progress every 10 permutations
        #if (perm %% 10 == 0 || perm == numperm) {
          #cat(sprintf("Permutations completed: %d / %d (%.1f%%)\n", perm, numperm, (perm / numperm) * 100))
        #}
        #permuted_populations <- sample(Island)  # Permute population labels
        #geninddata<-df2genind(dfloci[-1],pop = permuted_populations,sep = NULL,ploidy = 1) #convert to genind
        #hdata<-genind2hierfstat(geninddata,pop=NULL) #convert to heirfstat
        #resultsWC84<-genet.dist(hdata,diploid=FALSE,method="WC84") # run heirfstat
        #resultsWC84<-as.matrix(resultsWC84) #convert to matrix
        #resultsWC84[upper.tri(resultsWC84)] <- t(resultsWC84)[upper.tri(resultsWC84)] #copy so symettrical 
        ##remove pops that it did not occur in
        #i_order2<-order1
        #pops<-unique(permuted_populations)
        #nopops<-setdiff(i_order2, pops)
        #if(length(nopops)!=0){
          #i<-1
          #for (i in 1:length(nopops)){
            #i_order2<-i_order2[(i_order2 != nopops[i])]
          #}}
        #resultsWC84_reordered <- reorder_mat(mat = resultsWC84, order = i_order2) #reorder
        #resultsWC84_reordered[upper.tri(resultsWC84_reordered)]<-NA # clear top
        #WC84perms<-melt(resultsWC84_reordered) # melt so we can make a table of all perms
        #WC84perms$pairwisecomp<-paste(WC84perms$Var1,WC84perms$Var2,sep = " and ") # combine 2 populations into one column
        #WC84perms<-na.omit(WC84perms) # remove duplicate info
        #WC84perms<-WC84perms[WC84perms$Var1!=WC84perms$Var2,] # remove comparisons to itself 
        #print(WC84perms)
        #WC84perms2<-data.frame(pairwisecomp = WC84perms$pairwisecomp,perm = WC84perms$value)
        #permuted_fsts<-merge(permuted_fsts,WC84perms2, by="pairwisecomp",all=TRUE)
     #}
      # Calculate p-value as the proportion of permuted FST >= observed FST (but as a loop for each pairwise comp)
      #WC84_observed2<-WC84_observed[3:4]
      #obs_permuted_fsts<-merge(WC84_observed2,permuted_fsts, by="pairwisecomp",all = TRUE)
      #pvals<-c()
      #i<-1
      #for (i in 1:length(permuted_fsts$pairwisecomp)) {
      #  p_value <- mean(obs_permuted_fsts[i,-c(1:2)] >= obs_permuted_fsts$value[i], na.rm = TRUE)
      #  pvals<-c(pvals,p_value)
      #}
      #merge into 1 dataframe with fst and p-value
      #WC84_pvals<-data.frame(Pairwise_Comparisons = obs_permuted_fsts$pairwisecomp, Fst = obs_permuted_fsts$value, p.values = pvals)
      #write.csv(WC84_pvals,file = paste(folder_name2, temp2$Scientific_name[1],"_",temp2$best_identity[1],"_",temp2$blastlowest[1],"_",temp2$OTU[1],"_",runname, "_pairwisefst.csv", sep="")) #save
     }
  }
  files<-list.files(folder_name2)
  phisttable<-data.frame(matrix(ncol = 10, nrow = 0))
  i<-1
  for(i in 1:length(files)){
    #load a single phist table
    #grab phist values
    phistspecies<-read.csv(file= paste0(folderinput,files[i]), stringsAsFactors = F)
    phistspecies2<-phistspecies
    #correct p-values FDR
    #phistspecies2$fdr<-p.adjust(p=phistspecies2$value, method = "fdr")
    ##correct p-values bonforroni
    #phistspecies2$bon<-p.adjust(p=phistspecies2$value, method = "bonferroni")
    ##merge together
    #colnames(phistspecies2)<-c("Populations","Fst","p value","fdr p","bon p")
    #Species<-files[i]
    #phistall<-cbind(Species,phistspecies2)
    ##if significant
    #phistall$sigp <- with(phistall, ifelse(phistall$`p value` < 0.05, "yes", "no"))
    #phistall$sigfdr <- with(phistall, ifelse(phistall$`fdr p` < 0.05, "yes", "no"))
    #phistall$sigbon <- with(phistall, ifelse(phistall$`bon p`< 0.05, "yes", "no"))
    #bind to table
    phisttable<-rbind(phisttable,phistall)
  }
  colnames(phisttable)<-c("Species","Populations","FST value") #"p-value","fdr p-value","bonferroni p-value","sig p-value","sig fdr","sig bon"
  write.csv(phisttable, file = paste(folderoutput,"/",runname,"_fst_table_manyspecies.csv", sep=""))
}
#get info about taxonomy
taxonomy<-function(folder, datafile, runname){
  otus<-unique(datafile$OTU)
  taxa<-data.frame(matrix(ncol = 15, nrow = 0)) #create dataframe
  i<-1
  for (i in 1:length(otus)){
    temp2 <- datafile[datafile$OTU==otus[i],] #subet one OTU at a time
      #get obitools class
      try(lowest<-id2name(temp2$TaxID[1],db="ncbi"),  silent = TRUE)
      dflow<-as.data.frame(lowest[[1]])
      try(class<-tax_name(sci = dflow$name, get = "class", db = "ncbi"), silent=TRUE)
      if(length(class)==0){class<-data.frame(db="ncbi",query="NA",class="NA")}
      #get blast genus,family,order,class
      try(lowest<-id2name(temp2$blastlowestid[1],db="ncbi"),  silent = TRUE)
      dflow<-as.data.frame(lowest[[1]])
      try(bclass<-tax_name(sci = dflow$name, get = "class", db = "ncbi"), silent=TRUE)
      if(length(bclass)==0){bclass<-data.frame(db="ncbi",query="NA",class="NA")}
      try(border<-tax_name(sci = dflow$name, get = "order", db = "ncbi"), silent=TRUE)
      if(length(border)==0){border<-data.frame(db="ncbi",query="NA",class="NA")}
      try(bfamily<-tax_name(sci = dflow$name, get = "family", db = "ncbi"), silent=TRUE)
      if(length(bfamily)==0){bfamily<-data.frame(db="ncbi",query="NA",class="NA")}
      try(bgenus<-tax_name(sci = dflow$name, get = "genus", db = "ncbi"), silent=TRUE)
      if(length(bgenus)==0){bgenus<-data.frame(db="ncbi",query="NA",class="NA")}
      blastlowest2<-c(otus[i],temp2$best_identity[1],temp2$Scientific_name[1],temp2$Genus[1],temp2$Family[1],temp2$Order[1],class[3],temp2$blastpercent[1],temp2$blastlowest[1],bgenus[3],bfamily[3],border[3],bclass[3],temp2$phylum[1])
      blastlowest2<-unname(blastlowest2)
      taxa<- rbind(taxa,blastlowest2)
  }
  colnames(taxa)<-c("OTU","obitools_percent","obitools_species","obitools_genus","obitools_family","obitools_order","obitools_class","blast_percent","blast_species","blast_genus","blast_family","blast_order","blast_class","blast_phylum")
  write.csv(taxa,paste(folder,"/",runname,"_taxonomy_full.csv", sep=""))
  sumtaxa<-data.frame(matrix(ncol = 15, nrow = 0)) #create dataframe
  uotu<-length(unique(taxa$OTU))
  uobispecies<-length(unique(taxa$obitools_species))
  uobigen<-length(unique(taxa$obitools_genus))
  uobifam<-length(unique(taxa$obitools_family))
  uobiorder<-length(unique(taxa$obitools_order))
  uobiclass<-length(unique(taxa$obitools_class))
  ublspecies<-length(unique(taxa$blast_species))
  ublgen<-length(unique(taxa$blast_genus))
  ublfamily<-length(unique(taxa$blast_family))
  ublorder<-length(unique(taxa$blast_order))
  ublclass<-length(unique(taxa$blast_class))
  ublphy<-length(unique(taxa$blast_phylum))
  totals<-c(uotu,uobispecies,uobigen,uobifam,uobiorder,uobiclass,ublspecies,ublgen,ublfamily,ublorder,ublclass,ublphy)
  sumtaxa<-rbind(sumtaxa,totals)
  colnames(sumtaxa)<-c("otu_number","obitools_species","obitools_genus","obitools_family","obitools_order","obitools_class","blast_species","blast_genus","blast_family","blast_order","blast_class","blast_phylum")
  write.csv(sumtaxa,paste(folder,"/",runname,"_taxonomy_summary.csv", sep=""))
}
#part 1 of how to solve the memory issue for phist raw reads, have a function inside the function so it gets rid of files after running
epairphi<-function(temp2,runname,order,numperm,nind,popcolname,haplorestrict,folder_name2){
    #identify any populations with less than nind individuals
    populations<-popcolname
    ind<-data.frame(table(temp2[[popcolname]])) # count number of occurrences for each of the populations
    colnames(ind)<-c("site","ind") # rename the columns
    ind$ind<-as.numeric(ind$ind) # make number of individuals into a numeric
    lowpop<-as.vector(ind$site[ind$ind<nind])
    #remove any populations with less than nind individuals
    if(length(lowpop)!=0){
      i<-1
      for (i in 1:length(lowpop)){
        temp2<-temp2[(temp2[[popcolname]] != lowpop[i]),]
      }}
    #only keep haplotypes with at least 2 occurances
    if(haplorestrict==TRUE){
      temp3<-count(temp2$haplotype)
      temp3<-temp3[(temp3$freq==1),]
      i<-1
      lowhaplo<-temp3$x
      for (i in 1:length(lowhaplo)){
        temp2<-temp2[(temp2$haplotype != lowhaplo[i]),]
      }
    }
    if(length(unique(temp2[[popcolname]]))>1){
      #if subsample=TRUE then subsample if one population is more than double the other
      #if(subsample==TRUE){
      
      #}
      #convert DNA seqs to format for pairphist
      y <- as.list(temp2$sequences)
      y <- as.DNAbin(ape::as.alignment(y))
      z <- as.dna(y)
      #identify populations for testing
      Island<-temp2[[popcolname]]
      #calculate phist and p-values
      pst<-pairPhiST(z, Island, nperm=numperm, negatives=TRUE, showprogbar=FALSE) #calculates distance in code
      #How to create a prettier table
      #remove pops that it did not occur in
      i_order<-order
      pops<-unique(temp2[[popcolname]])
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
      write.csv(phist_reordered, file = paste(folder_name2, temp2$Scientific_name[1],"_",temp2$best_identity[1],"_",temp2$blastlowest[1],"_",temp2$OTU[1],"_",runname, "_phistat.csv", sep=""))
    } 
    }
#part 2 of how to solve the memory issue for phist raw reads, running phist
eDNA_haplotypes_pairphi_raw<-function(OTUvsspecies,folder,datafile,runname,order,numperm,nind,popcolname,haplorestrict){
  folder_name2<-paste(folder,"PhiStats",runname,"/",sep = "")
  dir.create(folder_name2)
  otus<-unique(datafile[[OTUvsspecies]])
  i<-1
  for (i in 1:length(otus)){
    popcolname<-popcolname
    order<-order
    numperm<-numperm
    nind<-nind
    haplorestrict<-haplorestrict
    folder<-folder
    temp2 <- datafile[datafile[[OTUvsspecies]]==otus[i],] #subset one OTU at a time
    epairphi(temp2,runname,order,numperm,nind,popcolname,haplorestrict,folder_name2)
  }
  }
#phist for the species and populations, not pairwise phist
overallphist<-function(folder,datafile,runname,popcolname,n_permutations,nind){
  otus<-unique(datafile$OTU) #get list for running the loop
  summary_stats<-data.frame(matrix(ncol = 10, nrow = 0)) #make empty df to fill
  i<-1
  for (i in 1:length(otus)){
    temp <- datafile[datafile$OTU==otus[i],] #subset each species
    pops<-unique(temp[[popcolname]]) #get list of unique populations
    N<-nrow(temp) #get sample size
    H<-length(unique(temp$sequences)) #get number of haplotypes
    nislands<-length(unique(temp[[popcolname]])) #get number of populations
    if(N>nind & H>1 & nislands>1){ #only run when it wont error because there aren't multiple samples, haplotypes, or populations
    population<-"overall" #name this run
    haplodf<-data.frame(populations = temp[[popcolname]],seqs = temp$sequences) #get df for seqs and populations for wc function
    ind<-data.frame(table(haplodf$populations)) # count number of occurances for each of the populations
    colnames(ind)<-c("site","ind") # rename the columns
    ind$ind<-as.numeric(ind$ind) # make number of individuals into a numeric
    notlowpop<-as.vector(ind$site[ind$ind>nind])
    #remove low pops
    haplodf<-haplodf[haplodf$populations %in% notlowpop,]
    countpop<-length(unique(haplodf$populations))
    if(countpop>1){
    haplodf2<-haplodf
    haplodf2$seqs <- as.numeric(as.factor(haplodf2$seqs)) #make seqs a factor
    test<-wc(haplodf2,diploid=FALSE) #calculate Fst
    fstwc<-test$FST[1] # save Fst value
    observedvalue<-fstwc
    temp2<-temp[temp[[popcolname]] %in% notlowpop,]
    y <- as.list(temp2$sequences) # get data ready for dist function
    y <- as.DNAbin(ape::as.alignment(y)) # get data ready for dist function
    #testdist2<-dist(y, method = "euclidean") #calculate euclidean distances
    testdist_dna<-dist.dna(y)
    testdist_dna_corrected<-cailliez(testdist_dna)
    Islands<-temp2[[popcolname]] #get list of pops
    Islands<-as.factor(Islands) #as factor
    #pamova<-pegas::amova(testdist2 ~  Islands, nperm = 1000) #run amova
    pamova2<-pegas::amova(testdist_dna_corrected ~  Islands, nperm = 1000) #run amova
    #sig2 <- setNames(pamova$varcomp$sigma2, rownames(pamova$varcomp))
    sig3 <- setNames(pamova2$varcomp$sigma2, rownames(pamova2$varcomp))
    #phitable<-getPhi(sig2) # Phi table
    phitable<-getPhi(sig3) # Phi table
    phist<-phitable[1]
    phist_pval<-pamova2$varcomp$P.value[1]
    #ttable<-table(temp$site,temp$sequences)
    #ttable2<-as.data.frame(ttable)
    #ttable3<-cast(ttable2,Var2~Var1,sum)
    #structure<-unique(temp$site)
    #structure2<-as.data.frame(structure)
    #aamova<-ade4::amova(ttable3,testdist2,structure2)
    # Permutation test for FST wc
    permuted_fsts <- numeric(n_permutations)
     for (perm in seq_len(n_permutations)) {
     permuted_populations <- sample(temp2[[popcolname]])  # Permute population labels
     haplodf<-data.frame(populations = permuted_populations, seqs = temp2$sequences) #get df for seqs and populations for wc function
     haplodf$seqs <- as.numeric(as.factor(haplodf$seqs)) #make seqs a factor
     test<-wc(haplodf,diploid=FALSE) #calculate Fst
     fstwc<-test$FST[1] # save Fst value
      permuted_fsts[perm] <- if (!is.null(fstwc)) fstwc else NA
    }
    ## Calculate p-value as the proportion of permuted FST >= observed FST
    wc_p_value <- mean(permuted_fsts >= observedvalue, na.rm = TRUE)
    basicstats<-c(temp$Scientific_name[1],temp$blastlowest[1],temp$OTU[1],population,N,H,fstwc,wc_p_value,phist,phist_pval)
    summary_stats<-rbind(summary_stats,basicstats)
    }}
    i<-1
    for (i in 1:length(notlowpop)){
      temp3 <- temp[temp[[popcolname]]==notlowpop[i],]
      N<-nrow(temp3)
      H<-length(unique(temp3$sequences))
      population<-notlowpop[i]
      samplength<-length(unique(temp3$region))
      if(N > nind & H > 1 & samplength > 1) {
      haplodf<-data.frame(populations = temp3$region,seqs = temp3$sequences)
      ind<-data.frame(table(haplodf$populations)) # count number of occurances for each of the populations
      colnames(ind)<-c("site","ind") # rename the columns
      ind$ind<-as.numeric(ind$ind) # make number of individuals into a numeric
      notlowpop2<-as.vector(ind$site[ind$ind>nind])
      #remove low pops
      haplodf<-haplodf[haplodf$populations %in% notlowpop2,]
      countpop<-length(unique(haplodf$populations))
      if(countpop>1){
      haplodf$seqs <- as.numeric(as.factor(haplodf$seqs))
     test<-wc(haplodf,diploid=FALSE)
     fstwc<-test$FST[1]
     observedvalue<-fstwc
     temp3<-temp3[temp3$region %in% notlowpop2,]
     y <- as.list(temp3$sequences) # get data ready for dist function
     y <- as.DNAbin(ape::as.alignment(y)) # get data ready for dist function
     #testdist2<-dist(y, method = "euclidean") #calculate euclidean distances
     Islands<-temp3$region #get list of pops
     Islands<-as.factor(Islands) #as factor
     testdist_dna<-dist.dna(y)
     testdist_dna_corrected<-cailliez(testdist_dna)
     #pamova<-pegas::amova(testdist2 ~  Islands, nperm = 1000) #run amova
     pamova2<-pegas::amova(testdist_dna_corrected ~  Islands, nperm = 1000)
     #sig2 <- setNames(pamova$varcomp$sigma2, rownames(pamova$varcomp))
     sig3 <- setNames(pamova2$varcomp$sigma2, rownames(pamova2$varcomp))
     #phitable<-getPhi(sig2) # Phi table
     phitable2<-getPhi(sig3) # Phi table
     phist<-phitable2[1]
     phist_pval<-pamova2$varcomp$P.value[1]
     permuted_fsts <- numeric(n_permutations)
     for (perm in seq_len(n_permutations)) {
       permuted_populations <- sample(temp3$region)  # Permute population labels
       haplodf<-data.frame(populations = permuted_populations, seqs = temp3$sequences) #get df for seqs and populations for wc function
       haplodf$seqs <- as.numeric(as.factor(haplodf$seqs)) #make seqs a factor
       test<-wc(haplodf,diploid=FALSE) #calculate Fst
       fstwc<-test$FST[1] # save Fst value
       permuted_fsts[perm] <- if (!is.null(fstwc)) fstwc else NA
     }
     # Calculate p-value as the proportion of permuted FST >= observed FST
     wc_p_value <- mean(permuted_fsts >= observedvalue, na.rm = TRUE)
     basicstats<-c(temp3$Scientific_name[1],temp3$blastlowest[1],temp3$OTU[1],population,N,H,fstwc,wc_p_value,phist,phist_pval)
     summary_stats<-rbind(summary_stats,basicstats)
     }}}
    }
  colnames(summary_stats)<-c("Species","blast lowest","OTU","Population","Sample Size","Number of Haplotypes","WC84 Fst","Fst p-value","phist","phist p-value")
  write.csv(summary_stats, file = paste(folder,"/",runname,"_overallphist.csv", sep=""))
}  
#phist accumulation curves
phistacccurv<-function(folder,datafile,runname,popcolname,numperm,sampleperm){
  folder_name2<-paste(folder,"PhiSTAccumulationCurves",runname,"/",sep = "")
  dir.create(folder_name2)
  otus<-unique(datafile$OTU) #get list for running the loop
  gooddata<-data.frame() #empty dataset
  i<-1
  for (i in 1:length(otus)){
    temp <- datafile[datafile$OTU==otus[i],] #subset each species
    ind<-data.frame(table(temp[[popcolname]])) # count number of occurances for each of the populations
    colnames(ind)<-c("site","ind") # rename the columns
    ind$ind<-as.numeric(ind$ind) # make number of individuals into a numeric
    minind<-ind$site[ind$ind > 1] #get populations that can't calculate phist
    temp2<-temp[temp$site %in% minind,] #remove populations that can't calculate phist
    maxind<-max(ind$ind[ind$ind!=max(ind$ind)]) #get the second highest sample size
    maxind<- maxind-1 #minus 1 because you can't resample the full amount
    if(length(unique(temp2[[popcolname]]))>1 & maxind > 2){
      folder_name3<-paste(folder_name2,temp$OTU[1],"_",temp$Scientific_name[1],"_",temp$best_identity[1],"_",temp$blastlowest[1],"_",temp$blastpercent[1],"/",sep = "") #name a folder for each species
      dir.create(folder_name3) #create the folder
      forfigure<-data.frame() #make an empty df to save stuff to
      sample_sizes <- seq(2, maxind, by = 1) # get sample sizes increasing from 2 to maxind by 1 increment
      # Loop through increasing sample sizes
      for (n in sample_sizes) {
        #remove any pops below sample size of n
        minind2<-ind$site[ind$ind > n] #get populations that can't calculate phist
        temp3<-temp2[temp2$site %in% minind2,] #remove populations that can't calculate phist
        #repeat so we have average and range
        for (i in 1:sampleperm) { 
          #subsample by n
          tempsample <- temp3 %>%
            group_by(temp3[[popcolname]]) %>%
            sample_n(n)
          #convert DNA seqs to format for pairphist
          y <- as.list(tempsample$sequences)
          y <- as.DNAbin(ape::as.alignment(y))
          z <- as.dna(y)
          #identify populations for testing
          Island<-tempsample[[popcolname]]
          #calculate phist and p-values
          pst<-pairPhiST(z, Island, nperm=numperm, negatives=TRUE, showprogbar=FALSE) 
          pstmelt<-melt(pst$PhiST) #format data for df
          pstmelt<-na.omit(pstmelt) #remove NAs due to matrix
          pstmelt$samplesize<-n # record the sample size
          forfigure<-rbind(forfigure,pstmelt) #save info
        }}
      forfigure$comp<-paste(forfigure$Var1, " and ",forfigure$Var2,sep = "") # get island pairwise labels
      popstoplot<-unique(forfigure$comp) # get list of each pairwise comparison
      #make a plot for each pairwise comparison
      for (i in 1:length(popstoplot)) {
        tempplot<-forfigure[forfigure$comp == popstoplot[i],] #get 1 pairwise comparison
        tempplot$samplesize<-as.numeric(tempplot$samplesize) 
        tempplot$value<-as.numeric(tempplot$value)
        #calculate mean sd, se, and se percent of mean for each sample size
        df.summary <- tempplot %>%
          dplyr::group_by(tempplot$samplesize) %>%
          dplyr::summarize(
            sd = sd(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE)/sqrt(length(value)),
            len = mean(value),
            percentofmean = (se/len)*100
          )
        #add if se is less than 5% of mean
        df.summary<-within(df.summary,{
          lessthan5per=NA
          lessthan5per[ percentofmean < 5] = "yes"
          lessthan5per[percentofmean > 5 | percentofmean == 5] = "no"
        })
        #add if se is less than 0.01
        df.summary<-within(df.summary,{
          lessthan001=NA
          lessthan001[ se < 0.01] = "yes"
          lessthan001[se > 0.01 | percentofmean == 0.01] = "no"
        })
        df.summary$species<-temp$Scientific_name[1]
        df.summary$OTU<-temp$OTU[1]
        df.summary$pairwisecomp<-tempplot$comp[1]
        colnames(df.summary)<-c("samplesize","sd","se","mean","se percent of the mean","se less than 5% of the mean","se less than 0.01","Species","OTU","pairwisecomp")
        write.csv(df.summary,file = paste(folder_name3,temp$OTU[1],"_",temp$Scientific_name[1],"_",temp$blastlowest[1],"_",tempplot$comp[1],"_",runname,"_pairwisephist_accumulation_data.csv"))
        #make figure
        jpeg(paste(folder_name3,temp$OTU[1],"_",temp$Scientific_name[1],"_",temp$blastlowest[1],"_",tempplot$comp[1],"_",runname,"_pairwisephist_accumulation_curve.jpeg",sep = ""), width=12, height=12, units="in",res=600)
        tempplot2<-ggplot(df.summary, aes(x = samplesize, y = mean)) +
          geom_point(data = tempplot,aes(x=samplesize,y=value))+
          geom_line(lwd=1.5,color="blue") +  
          geom_errorbar(data = df.summary, aes(ymin = mean - se, ymax = mean + se), width = 1,color = "blue") +  # Adds error bars
          labs(
            title = tempplot$comp[1],
            x = "Sample Size",
            y = "Pairwise ΦST")+
          theme (panel.grid.major=element_blank(),
                 panel.border = element_rect(color = "white",fill=NA),
                 panel.background = element_rect(fill="grey95"),
                 axis.text.x = element_text(size=18),
                 axis.text.y = element_text(size=18), 
                 axis.title.x = element_text(size=22),
                 axis.title.y = element_text(size=22)) 
        print(tempplot2)
        dev.off()
        last_row <- tail(df.summary, n=1) # get last line of df
        #get final info and save to df if it is well sampled
        if(last_row$`se less than 5% of the mean`=="yes" & last_row$`se less than 0.01` =="yes"){
          locations<-separate(last_row, pairwisecomp, into = c("island1", "island2"), sep = " and ") # get populations to calculate phist
          is1<-locations$island1[1] # island 1
          is2<-locations$island2[1] # island 2
          n_final<-locations$samplesize[1] # number of subsample
          temp4<-temp2[temp2[[popcolname]] == is1 | temp2[[popcolname]] == is2,] #subset for the right 2 islands
          # repeat for the number that sampleperm is
          tempstorage<-data.frame()
          for (i in 1:sampleperm) { 
            #subsample by n_final
            tempsample <- temp4 %>%
              group_by(temp4[[popcolname]]) %>%
              sample_n(n_final)
            #convert DNA seqs to format for pairphist
            y <- as.list(tempsample$sequences)
            y <- as.DNAbin(ape::as.alignment(y))
            z <- as.dna(y)
            #identify populations for testing
            Island<-tempsample[[popcolname]]
            #calculate phist and p-values
            pst<-pairPhiST(z, Island, nperm=1000, negatives=TRUE, showprogbar=FALSE) 
            pstmelt<-melt(pst$PhiST) #format data for df
            pstmelt<-na.omit(pstmelt) #remove NAs due to matrix
            pmelt<-melt(pst$p) #format data for df
            pmelt<-na.omit(pstmelt) #remove NAs due to matrix
            pstmelt$samplesize<-n # record the sample size
            pstmelt$pvalue<-pmelt$value
            tempstorage<-rbind(tempstorage,pstmelt)
          }
          avphist<-mean(tempstorage$value)
          sdphist<-sd(tempstorage$value)
          sephist<-sdphist/sqrt(length(tempstorage$value))
          avpval<-mean(tempstorage$pvalue)
          sdpval<-sd(tempstorage$pvalue)
          sepval<-sdpval/sqrt(length(tempstorage$pvalue))
          suminfo<-data.frame(tempstorage$Var1[1],tempstorage$Var2[1],tempstorage$samplesize,avphist,sdphist,sephist,avpval,sdpval,sepval)
          speciesinfo<-data.frame("Species" = temp$Scientific_name[1],
                                  "match" = temp$best_identity[1],
                                  "BlastSpecies" = temp$blastlowest[1],
                                  "OTU" = temp$OTU[1])
          longsuminfo<-cbind(speciesinfo,suminfo)
          gooddata<-rbind(gooddata,longsuminfo) #save info
          }
        }
      }
  }
  gooddata$percentsemeanphist<-gooddata$sephist/gooddata$avphist
  gooddata$percentsemeanpval<-gooddata$sepval/gooddata$avpval
  colnames(gooddata)<-c("Species","match","Blast Species","OTU","Island1","Island2","SampleSizeSubset","Mean Phist","sd Phist","se Phist","Mean pvalue","sd pvalue","se pvalue","se percent mean phist","se percent mean pval")
  gooddata$IslandComp<-paste(gooddata$Island1, " and ",gooddata$Island2,sep = "") # get island pairwise labels
  write.csv(gooddata,file = paste(folder_name2,runname,"_wellsampleddatafromse_phistlist.csv",sep = ""))
}



############################### my inputs ###############################
#set working directory
setwd("/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI")
#folder where you want things
folder1<-"/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/obitools/maxee025/"
#load metadata
metadata1<-read.csv("Metadata_MHI.csv", stringsAsFactors = F)
#load haplo table from JAMP output
table1<- read.csv("obitools/maxee025/MHI_ee025_haplo_table.csv", stringsAsFactors=F) # JAMP output
#load output from obitools
obitools1<-read.table("obitools/maxee025/MHI_ee025_tax_assigned.tab",header = TRUE, sep = "\t",fill=TRUE, quote="",blank.lines.skip=FALSE) # ecotag/obitools output
#load modified blastn output see notes at the top of this file
blast1<-read.table("obitools/maxee025/blast_MHI_ee025_zizka_awk.txt", header = FALSE, sep = "\t",fill=TRUE, blank.lines.skip=FALSE,quote = "") # blastn and awk output
#code for the translation of codons
numcode_vert_mt<-2

######################################################################################
#NOW RUN IT
load_species(metadata = metadata1,
             table = table1,
             obitools = obitools1,
             blast = blast1,
             lowestmatch = 0.90, #has to be 90% match
             lowestreads = 10, #sequences have to have at least 10 reads
             OTUvsspecies = "OTU",
             folder = folder1, #folder where everything will be saved to
             nblanks = 11,     # remember to change per run
             nsamples = 256,    # remember to change per run
             runname = "cytb_ee025_r10_m90",   # remember to change per run
             scaled1 = 0.5, #below this is smallest scale for scaled by reads, this is <50%
             scaled2 = 0.75, #this is 2nd  scale for scaled by reads, this is >50% but < 75%
             scaled3 = 0.9, #this is 3rd and 4th scale for scaled by reads, 3rd is >75% but < 90%, 4th is >90%
             removefirstbp = TRUE) # see why to do this in notes at the top of file

#read in alignment file created after the load_species function
alignment1<-read.fasta("obitools/maxee025/cytb_ee025_r10_m90haplotypes.fasta")

#REMOVE SEQs WITH STOP CODONS
remove_stop_codons(alignment = alignment1, # remember to read this in
                   presenceabsence = presenceabsencecytb_ee025_r10_m90, # presence-absence proxies, this df was created from load_species function
                   rarefied = rarefiedcytb_ee025_r10_m90, # scaled-by-reads proxies, this df was created from load_species function
                   normalized = normalizedcytb_ee025_r10_m90, # raw reads proxies, this df was created from load_species function
                   alldata = final_combinedcytb_ee025_r10_m90, # all data, this df was created from load_species function
                   runname = "cytb_ee025_r10_m90_nSC", #run name
                   frame1 = 0, #what frame to start codon on
                   numcode1 = numcode_vert_mt, #code for translation
                   folder = folder1) #folder where everything will be saved to

#reduce reads that are too large
reducereads(datafile=normalizedcytb_ee025_r10_m90, # this df was created from load_species function
            nsamp=256, #number of samples
            alignment=alignment1, #alignment fileloaded previously
            runname="cytb_ee025_raw", # run name
            folder=folder1, #folder where everything will be saved to
            metadata=metadata1) 

#get info about haplotypes and counts per species and population within a species
summary_eDNA(OTUvsspecies = "OTU",
             presenceabsence = presenceabsencecytb_ee025_r10_m90_nSC, # presence-absence proxies, this df was created from load_species function and modified by remove stop codon
             rarefied = rarefiedcytb_ee025_r10_m90_nSC, # scaled-by-reads proxies, this df was created from load_species function and modified by remove stop codon
             normalized = normalizedcytb_ee025_r10_m90_nSC, # raw reads proxies, this df was created from load_species function and modified by remove stop codon
             runname = "cytb_ee025_r10_m90_nSC", # run name
             nhap = 1, # minimum number of haplotypes
             folder = folder1, #folder where everything will be saved to
             haplorestrict = FALSE #TRUE: haplotypes need to be in multiple samples, FALSE: keep all data
             )

#get pairwise phist
#repeat this for popgen_rcytb_r10_m95_nSC, popgen_pacytb_r10_m95_nSC, popgen_ncytb_r10_m95_nSC
eDNA_haplotypes_pairphi(OTUvsspecies = "OTU",
                        folder = folder1, #folder where everything will be saved to
                        datafile = popgen_rcytb_ee025_r10_m90_nSC, # scaled-by-reads proxies, this was created from summary_eDNA function
                        runname = "cytb_r10_m90_nSC_ee025_r", #run name
                        numperm = 1000, #number of permutations for pairwise phist to calculate p-value
                        order = c("Kauai","Oahu","Maui","Hawaii"), #order for you final pairwise phist table
                        nind=19, #this means that there needs to be 20 or greater individuals per population to calculate pairwise phist
                        popcolname = "site", #which population level to use
                        haplorestrict = FALSE #I never used this function but I made it to test out
                        )

phistacccurv(folder=folder1, #folder where everything will be saved to
             datafile= popgen_rcytb_ee025_r10_m90_nSC, # scaled-by-reads proxies, this was created from summary_eDNA function
             runname="cytb_ee025_r10_m90_nSC_SE_rep100", #run name
             popcolname="site", #population level
             numperm=1, #this is for calculating p-value for phist which is not necessary so I just have it at 1
             sampleperm=100) #number of permutations for each sample size

eDNA_haplotypes_pairphi_raw(OTUvsspecies = "OTU",
                            folder = folder1, #folder where everything will be saved to
                            datafile = popgen_rcytb_ee025_r10_m90_nSC, # scaled-by-reads proxies, this was created from summary_eDNA function
                            runname = "cytb_r10_m90_nSC_ee025_site_raw_reduced_n600", # run name
                            numperm = 1000, #number of permutations for pairwise phist to calculate p-value
                            order = c("Kauai","Oahu","Maui","Hawaii"), #order for you final pairwise phist table
                            nind=599, #minimum number of proxies needed to run comparison
                            popcolname = "site", #which level of population to test
                            haplorestrict = FALSE)

#get table for phist and p-values stats
convert_phist_to_table(folderinput="/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/obitools/maxee025/PhiStatscytb_r10_m90_nSC_ee025_r/",
                       folderoutput="/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/obitools/maxee025",
                       runname="cytb_r10_m90_nSC_ee025_r")

#GET OVERALL STATS PER SPECIES
statstable(folder = "/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/obitools/maxee025",
           datafile = popgen_rcytb_ee025_r10_m90_nSC, # scaled-by-reads proxies, this was created from summary_eDNA function
           runname = "cytb_ee025_r10_m90_nSC_r_site", # run name
           popcolname = "site") #which level of population to test

#get fasta files fastafilesperspecies<-function(datafile,popcolname,runname,folder)
fastafilesperspecies(datafile=popgen_rcytb_ee025_r10_m90_nSC,
                     popcolname="site",
                     runname="cytb_r10_m90_nSC_ee05_site_r",
                     folder="/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/obitools/maxee025"
                       )
fastafilesperspeciesperpopulation(datafile=popgen_rcytb_ee025_r10_m90_nSC,
                                  popcolname="site",
                                  runname="cytb_r10_m90_nSC_ee025_n_bypop_region",
                                  folder=folder1)
fastatonexus(folderinput="/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/obitools/maxee025/fastafilescytb_r10_m90_nSC_ee025_r_bypop_site",
             folderoutput="/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/obitools/maxee025",
             runname="cytb_ee025_r_bypop_site_nexus")
  
MMD_all(folder=folder1,
        table = popgen_rcytb_ee025_r10_m90_nSC,
        runname = "cytb_r10_m90_nSC_ee025_r_site",
        popcolname = "site")

######## figure out if resample A. nigrofuscus if there would be a lot of haplotypes still
resample_haplo<-function(n_permutations,datafile,num_resample,comp_name,proxytype){
  permuted_fsts <- numeric(n_permutations)
  for (perm in seq_len(n_permutations)) {
    AcanigF_stiss<-datafile[sample(nrow(datafile), num_resample), ]
    nhaplo<-length(unique(AcanigF_stiss$haplotype))
    permuted_fsts[perm] <- if (!is.null(nhaplo)) nhaplo else NA
  }
  ave_haplo <- mean(permuted_fsts, na.rm = TRUE)
  sd_haplo<-sd(permuted_fsts, na.rm = TRUE)
  se_haplo<-sd(permuted_fsts, na.rm = TRUE)/sqrt(length(permuted_fsts))
  max_haplo<-max(permuted_fsts)
  min_haplo<-min(permuted_fsts)
  resampletable<<-c(comp_name,proxytype,num_resample,ave_haplo,sd_haplo,se_haplo,max_haplo,min_haplo)
}
#subsample only the species of interest
Anigrofuscus<-rawreads[rawreads$OTU=="OTU_33",]
Anigrofuscus_pa<-presenceabsence[presenceabsence$OTU=="OTU_33",]
Anigrofuscus_r<-rarefied[rarefied$OTU=="OTU_33",]
#make datatable to save to
AcanigF_resample<-data.frame()
#run function
resample_haplo(n_permutations=1000,
               datafile=Anigrofuscus,
               num_resample=232,
               comp_name = "Chaetodon fremblii Sample Size n232",
               proxytype="raw reads")
#rbind to datatable
AcanigF_resample<-rbind(AcanigF_resample,resampletable)
#rename columns
colnames(AcanigF_resample)<-c("Compared to","eDNA proxy","resampled N","Mean Num. Haplo","SD","SE","Max","Min")
write.csv(AcanigF_resample,file = "ResamplingAcanigrofuscusToSeeIfLargeHaploNumDueToReads.csv")

####################### not used for the paper ############################
# pairwiseFst<-function(folder, datafile, runname, popcolname,order1,numperm,nind,haplorestrict)
pairwiseFst(folder = folder1,
            datafile = presenceabsence,
            runname = "cytb_r10_m90_nSC_ee025_site_pa",
            numperm = 1000,
            order1 = c("Kauai","Oahu","Maui","Hawaii"),
            nind=2,
            popcolname = "site",
            haplorestrict = FALSE)

#overallphist(folder,datafile,runname,popcolname,n_permutations)
overallphist(folder=folder1,
             datafile=rarefied2,
             runname="cytb_r10_m90_nSC_ee025_site_r_n10_rep1000_subsetmaui",
             popcolname="site",
             n_permutations=1000,
             nind = 9)


convert_phist_to_table(folderinput="/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/obitools/maxee025/PhiStatscytb_r10_m90_nSC_ee03_site_normalized_reduced_test_n10/",
                       folderoutput="/Users/taylorely/Documents/Grad_Work/ProcessingSequences/MHI/obitools/maxee025",
                       runname="cytb_ee025_m90_nSC_site_raw_n10_reduced")

nmdsbyspecies(folder=folder1,
              datafile=ee025_raw,
              runname="cytb_ee025_raw",
              popcolname = "Region",
              nsamples=256)

cols1<-c("lightblue","lightblue","orange","darkblue","pink","lightblue","lightgreen","darkgreen","lightblue","orange","lightgreen","lightgreen","darkblue","darkgreen","orange","pink","darkgreen","tan","darkblue","tan","darkblue","purple","orange","tan","pink","lightgreen","pink","purple","darkgreen","tan","purple","purple")

jpeg("obitools/maxee025/allfish_ee025_raw_pca_island_location.jpg", width=12, height=10, units= "in", res=600)
nmdsplot_IbyY_nf <- plot_ordination(physeq_raw, raw_ord, "samples", color="Location") + 
  facet_wrap(~Island) + geom_polygon(aes(alpha = 0.8,fill=Location)) + geom_point(size=5) + ggtitle("NMDS of raw reads")+
  scale_color_manual(values=cols1)+
  scale_fill_manual(values=cols1)
print(nmdsplot_IbyY_nf)
dev.off()

sampledf_raw <- data.frame(sample_data(physeq_raw)) # get samples data
vOTU <- otu_table(physeq_raw) # convert phyloseq to vegan format
vOTU <- t(vOTU)
vOTU<-as(vOTU, "matrix")
#Bray curtis dissimilarity matrix
d_carn_raw <- vegdist(vOTU, method="bray") 
#Permanova
app_var_raw_location<-adonis(vOTU ~ sampledf_raw$Island+sampledf_raw$Location+sampledf_raw$Samples, method = "bray")$aov.tab
write.csv(app_var_raw_location,file="obitools/maxee025/PERMANOVA_raw.csv")

