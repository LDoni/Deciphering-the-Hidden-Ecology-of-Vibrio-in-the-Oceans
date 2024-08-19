# Deciphering-the-Hidden-Ecology-of-Vibrio-in-the-Oceans

# Metagenomic READS Data Processing Pipeline

This pipeline processes metagenomic data by performing the following steps: trimming reads, classifying them with Kraken2 using the RefSeq database, extracting Vibrio reads, classifying these reads with Enterobase, and finally refining the classifications using Bracken.

## Steps Overview

1. **Trimming**: Pre-process the raw reads to remove low-quality bases and adapters.
2. **Kraken2 Classification with RefSeq**: Classify the trimmed reads at the genus level using the RefSeq database.
3. **Bracken Genus-Level Refinement (RefSeq)**: Refine the genus-level classifications using Bracken.
4. **Extract Vibrio Reads (RefSeq)**: Extract reads classified as Vibrio from the RefSeq output.
5. **Kraken2 Classification with Enterobase**: Classify the extracted Vibrio reads at the species level using the Enterobase database.
6. **Bracken Species-Level Refinement (Enterobase)**: Refine the species-level classifications using Bracken.
7. **Extract Vibrio Reads (Enterobase)**: Extract Vibrio reads from the Enterobase output.

## Script

```bash

# Set the current path variable
currpath=$(pwd)

# Set the sample name file
SAMPLE="samples_TARA.txt"

# RefSeq Kraken2 database directory
REFSEQ_KRA_DB="/home/LAPO/krakenDB_bacteria_refseq"

# Enterobase Kraken2 database directory
ENT_KRA_DB="/home/LAPO/Desktop/TARA/enterobase_vibrio_KRAKEN2_db"

# Output directories
TRIM_OUT="${currpath}/trimming"
REFSEQ_OUT="${TRIM_OUT}/kraken_refseq"
ENT_OUT="${TRIM_OUT}/kraken_entero/entero_extraction"

#######################
# 1. Trimming of reads
#######################

cat $SAMPLE | parallel -j 10 'trim_galore \
--cores 5 --paired {}_1.fastq.gz {}_2.fastq.gz \
--trim-n --illumina --fastqc_args "-outdir ${TRIM_OUT}/quality_control" -o ${TRIM_OUT}'

###############################################################
# 2. Kraken2 classification with RefSeq database after trimming
###############################################################

cat $SAMPLE | parallel -j 10 \
"kraken2 --db ${REFSEQ_KRA_DB} \
--paired ${TRIM_OUT}/{}_1_val_1.fq.gz ${TRIM_OUT}/{}_2_val_2.fq.gz \
--threads 20 \
--use-names \
--gzip-compressed \
--report ${REFSEQ_OUT}/{}_report-kraken_refseq.txt \
--output ${REFSEQ_OUT}/{}_output_refseq.kraken"

#########################################################
# 3. Bracken genus-level refinement after Kraken2 (RefSeq)
#########################################################

cat $SAMPLE | parallel -j 10 \
"bracken -d ${REFSEQ_KRA_DB} \
-i ${REFSEQ_OUT}/{}_report-kraken_refseq.txt \
-o ${REFSEQ_OUT}/{}_bracken_genus_refseq.txt \
-l G"

combine_bracken_outputs.py --files *.txt -o ../braken_all_REFSEQ_prokEprot_merged/braken_all_REFSEQ_prokEprot_merged.csv ```


# FIG1 A




library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(readr)
library(plyr)
library(dplyr)
library(metagMisc)
library(forcats)
library(ggh4x)
#####                      importing data                           ####



abund_table<-read.csv("input/braken_all_REFSEQ_prokEprot_merged.csv",row.names=NULL, check.names=FALSE,sep = ";")

## aggregare colonne ripetute facendo la media ( per dereplicare le repliche dei campioni)

abund_table<-aggregate(.~Genus,abund_table,mean)
abund_table<-t(abund_table)

head(abund_table)
ncol(abund_table)
nrow(abund_table)
#TAXONOMY table
OTU_taxonomy<-read.csv("input/braken_all_REFSEQ_prokEprot_merged_fract_TAX.csv",row.names=1,check.names=FALSE,sep = ";")
nrow(OTU_taxonomy)

##  Metatable daereplicato su excel, filtro avanzato, non mostrare duplicati!
meta_table<-read.csv("input/braken_all_REFSEQ_prokEprot_merged_fract_METAWcooRd.csv",row.names=1,check.names=FALSE,sep = ";")
colnames(meta_table)
nrow(meta_table)
ncol(meta_table)




#Convert the data to phyloseq format
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = T)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)
physeq<-merge_phyloseq(phyloseq(OTU, TAX, SAM))



#Top 15 taxa

genus.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Genus"], sum, na.rm=TRUE)
top5phyla = names(sort(genus.sum, TRUE))[1:15]
GP1 = prune_taxa((tax_table(physeq)[, "Genus"] %in% top5phyla), physeq)



merge_samples_mean <- function(physeq, group){
  group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
  merged <- merge_samples(physeq, group)
  x <- as.matrix(otu_table(merged))
  if(taxa_are_rows(merged)){ x<-t(x) }
  out <- t(x/group_sums)
  out <- otu_table(out, taxa_are_rows = TRUE)
  otu_table(merged) <- out
  return(merged)
}

#merge by zone
phyloseq_obj_css_Zone<- merge_samples_mean(GP1,"Zone")


data_glom<- psmelt(phyloseq_obj_css_Zone) # create dataframe from phyloseq object
data_glom$Genus <- as.character(data_glom$Genus)
data_glom$Abundance<-(data_glom$Abundance*100)

 

newSTorder =c( "ANE","ANW","ASE",
               "ASW","ION" , "IOS","PON","PSE","PSW","MED" ,"RED","ARC","SOC")
data_glom$Sample<- as.character(data_glom$Sample)
data_glom$Sample <- factor(data_glom$Sample, levels=newSTorder)








ggplot(data_glom, aes(x = Sample, y = fct_reorder(Genus, NEG_TOT_ABUNDANCE), size=Abundance, fill=Ocean)) + 
  geom_point(alpha=0.5, shape=21, color="black") +
  scale_size_continuous(name = "Counts ", 
                        range = c(0, 16)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = NULL, y = NULL)+
  scale_size_continuous(name = "Abundance %",
                        breaks = bb,
                        labels = ll,
                        range = c(0.1, 20) )+
  scale_fill_viridis(discrete=TRUE, guide=FALSE, option="D")+ theme(axis.text=element_text(colour="black"))+
  scale_x_discrete(position = "top")+
   theme(text = element_text(size = 20))+scale_y_discrete(limits=rev)


###################################################################################### 
### HEATMAPS FIG 1B
############################################################################################################################# 
VIBRIO <- subset_taxa(GP1, Genus=="Vibrio" )



###################################################################################### 
### fraction
############################################################################################################################# 

variable1 = as.character(get_variable(VIBRIO, "Zone"))
variable2 = as.character(get_variable(VIBRIO, "Fraction3"))
# variable3 = as.character(get_variable(physeq, "Distance1"))

sample_data(VIBRIO)$NewPastedVar <- mapply(paste, variable1, variable2  
                                           , sep = "_")
# sample_data(VIBRIO)

physeq_MERGED_zone_FRACTIONS<-merge_samples_mean(VIBRIO, "NewPastedVar")

sample_variables(physeq_MERGED_zone_FRACTIONS)
sample_data(physeq_MERGED_zone_FRACTIONS)
sample_names(physeq_MERGED_zone_FRACTIONS)


#per far quello che facevo su excel
rm(df)
ttt<-nrow(as.data.frame(sample_names(physeq_MERGED_zone_FRACTIONS)))
df <- data.frame(matrix(ncol = 1, nrow = ttt))
df$names<-(as.data.frame(sample_names(physeq_MERGED_zone_FRACTIONS)))

df<-df[,2]
colnames(df)[1] <- "Samples"

foo <- data.frame(do.call('rbind', strsplit(as.character(df$Samples),'_',fixed=TRUE)))
df$Zone <-foo$X1
df$Size <-foo$X2
colnames(df)

rownames(df) <- df[,1]

df


physeq_MERGED_zone_FRACTIONS
#cambio il metadata

sample_data(physeq_MERGED_zone_FRACTIONS)<-sample_data(df)
physeq_MERGED_zone_FRACTIONS

sample_variables(physeq_MERGED_zone_FRACTIONS)
sample_data(physeq_MERGED_zone_FRACTIONS)



physeq_MERGED_zone_FRACTIONS<-subset_samples(physeq_MERGED_zone_FRACTIONS, !(  Size=="na" ))

data_glom_fract<- psmelt(physeq_MERGED_zone_FRACTIONS) # create dataframe from phyloseq object
data_glom_fract$Genus <- as.character(data_glom_fract$Genus) #convert to character
data_glom_fract$Abundance<-(data_glom_fract$Abundance*100)


newSTorder =c( "ANE","ANW","ASE",
               "ASW","ION" , "IOS","PON","PSE","PSW","MED" ,"RED","ARC","SOC")
data_glom_fract$Zone<- as.character(data_glom_fract$Zone)
data_glom_fract$Zone <- factor(data_glom_fract$Zone, levels=newSTorder)


data_glom_fract<-data_glom_fract %>% 
  mutate(Size = str_replace(Size, "Prokaryotes Fractions", "FLB"))%>% 
  mutate(Size = str_replace(Size, "Protist Fraction", "PAB"))



ggplot(data_glom_fract, aes(Zone, Size, fill= Abundance)) + 
  geom_tile()+ theme(axis.text.x = element_text(angle = 90))+scale_fill_gradient( na.value="black")+
  theme(panel.background = element_rect(fill = 'black'),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +scale_fill_viridis_c(option = "H",limits = c(0,4))+ 
  theme(legend.position = "none" )+ 
  theme(axis.text=element_text(colour="black"),axis.title.x = element_blank(),
        axis.text.x=element_blank(),axis.ticks.x =element_blank())+scale_y_discrete(expand = c(0, 0))+
  scale_x_discrete(expand = c(0, 0))+labs(y= "Fractions")+ theme(text = element_text(size = 20))



###################################################################################### 
### Depth
############################################################################################################################# 
variable1 = as.character(get_variable(VIBRIO, "Zone"))
variable2 = as.character(get_variable(VIBRIO, "Depth"))
# variable3 = as.character(get_variable(physeq, "Distance1"))

sample_data(VIBRIO)$NewPastedVar <- mapply(paste, variable1, variable2  
                                           , sep = "_")
sample_data(VIBRIO)

physeq_MERGED_zone_DEPTH<-merge_samples_mean(VIBRIO, "NewPastedVar")
sample_data(physeq_MERGED_zone_DEPTH)

#per far quello che facevo su excel
rm(df)
# (as.data.frame(sample_names(physeq_MERGED_zone_DEPTH)))
ttt<-nrow(as.data.frame(sample_names(physeq_MERGED_zone_DEPTH)))
df <- data.frame(matrix(ncol = 1, nrow = ttt))
df$names<-(as.data.frame(sample_names(physeq_MERGED_zone_DEPTH)))
df<-df[,2]
colnames(df)[1] <- "Samples"
foo <- data.frame(do.call('rbind', strsplit(as.character(df$Samples),'_',fixed=TRUE)))
df$Zone <-foo$X1
df$W_Layer <-foo$X2
colnames(df)

rownames(df) <- df[,1]

df


# physeq_MERGED_zone_DEPTH
#cambio il metadata

sample_data(physeq_MERGED_zone_DEPTH)<-sample_data(df)
physeq_MERGED_zone_DEPTH

head(sample_variables(physeq_MERGED_zone_DEPTH))
head(sample_data(physeq_MERGED_zone_DEPTH))

physeq_MERGED_zone_DEPTH<-subset_samples(physeq_MERGED_zone_DEPTH, !( ( W_Layer=="MIX" | W_Layer=="ZZZ"| W_Layer=="NA")))

data_glom_depth<- psmelt(physeq_MERGED_zone_DEPTH) # create dataframe from phyloseq object
data_glom_depth$Genus <- as.character(data_glom_depth$Genus) #convert to character
data_glom_depth$Abundance<-(data_glom_depth$Abundance*100)

depth_order = c(
  "SRF",
  "DCM",
  "MES")



data_glom_depth$W_Layer<- as.character(data_glom_depth$W_Layer)
data_glom_depth$W_Layer <- factor(data_glom_depth$W_Layer, levels=rev(depth_order) )

newSTorder =c( "ANE","ANW","ASE",
               "ASW","ION" , "IOS","PON","PSE","PSW","MED" ,"RED","ARC","SOC")
data_glom_depth$Zone<- as.character(data_glom_depth$Zone)
data_glom_depth$Zone <- factor(data_glom_depth$Zone, levels=newSTorder)

min(data_glom_depth$Abundance)
max(data_glom_depth$Abundance)
mean(data_glom_depth$Abundance)

bb_depth <- c(0, 0.5,1,1.5,2,2.5,3) # define breaks.
ll_depth <- c("0%","0.5%","1%","1.5%","2%","2.5%","3%") # labels.


#########################  plot 
library(ggh4x)

HMDEPTH<-  ggplot(data_glom_depth, aes(Zone, W_Layer, fill= Abundance)) + 
  geom_tile()+ theme(axis.text.x = element_text(angle = 90))+scale_fill_gradient( na.value="black")+
  theme(panel.background = element_rect(fill = 'black'),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +scale_fill_viridis_c(option = "H",limits = c(0,4))+ 
    guides(fill = guide_colourbar(barwidth = 40, barheight = 3,direction = "horizontal",title = "Abundance %",title.position = "top",title.hjust = 0.5))+
     theme(legend.position = "bottom" ) 
    theme(axis.text=element_text(colour="black"),axis.title.x = element_blank(),
                                               axis.text.x=element_blank(),axis.ticks.x =element_blank())+scale_y_discrete(expand = c(0, 0))+
    scale_x_discrete(expand = c(0, 0))+labs(y= "Depth")+ theme(text = element_text(size = 20)) 




######################################################################
# 4. Extract reads classified as Vibrio using RefSeq classified reads
######################################################################

cat $SAMPLE | parallel -j 10 \
"python3 ~/extract_kraken_reads.py \
-k ${REFSEQ_OUT}/{}_output_refseq.kraken \
-s1 ${TRIM_OUT}/{}_1_val_1.fq.gz \
-s2 ${TRIM_OUT}/{}_2_val_2.fq.gz \
--taxid 662 --include-children \
-o ${REFSEQ_OUT}/{}_extracted_reads_vibrio_1.fa \
-o2 ${REFSEQ_OUT}/{}_extracted_reads_vibrio_2.fa"

#####################################################################################
# 5. Kraken2 classification with Enterobase database using Vibrio extracted reads
#####################################################################################

cat $SAMPLE | parallel -j 10 \
"kraken2 --db ${ENT_KRA_DB} \
--paired ${REFSEQ_OUT}/{}_extracted_reads_vibrio_1.fa ${REFSEQ_OUT}/{}_extracted_reads_vibrio_2.fa \
--threads 20 \
--use-names \
--report ${ENT_OUT}/{}_report-kraken_entero.txt \
--output ${ENT_OUT}/{}_output_entero.kraken"

###############################################################
# 6. Bracken species-level refinement after Kraken2 (Enterobase)
###############################################################

cat $SAMPLE | parallel -j 10 \
"bracken -d ${ENT_KRA_DB} \
-i ${ENT_OUT}/{}_report-kraken_entero.txt \
-o ${ENT_OUT}/{}_bracken_species_entero.txt \
-l S"
combine_bracken_outputs.py --files *.txt -o ../braken_all_ENTERO_prokEprot_merged/braken_all_ENTERO_prokEprot_merged.txt



library(phyloseq)
library(ggside)
library(vegan)
library(ggplot2)
library(readr)
library(plyr)
library(dplyr)
library(hrbrthemes)
library(ggpubr)
library(viridis)
library(ggforce)
library(concaveman)
library(ggpubr)
library(ggside)
library(ggdist)
#####                      importing data                           ####


abund_table<-read.csv("input/OTUtable_ALL_ENTERO_derep.csv",row.names=1, check.names=FALSE,sep = ";")
abund_table<-t(abund_table)

#TAXONOMY table
OTU_taxonomy<-read.csv("input/braken_all_ENTERO_prokEprot_merged_fract_INPUT_R_OTU_TAX.csv",row.names=1,check.names=FALSE,sep = ";")
nrow(OTU_taxonomy)


meta_table<-read.csv("input/braken_all_REFSEQ_prokEprot_merged_fract_METAWcooRd.csv",row.names=1,check.names=FALSE,sep = ";")

#Convert the data to phyloseq format
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = T)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)
physeq<-merge_phyloseq(phyloseq(OTU, TAX, SAM))

#####################################################################################################################

#####   alfadiversity

#####################################################################################################################

###############   #fig 1C

OTU_ceiling = otu_table(as.matrix(ceiling(abund_table)), taxa_are_rows = T)

physeqCEILING<-merge_phyloseq(phyloseq(OTU_ceiling, TAX, SAM))

physeqCEILING<-subset_samples(physeqCEILING, !(  Fraction=="NA" ))

physeqCEILING<-subset_samples(physeqCEILING, !( ( Depth=="MIX" | Depth=="ZZZ"| Depth=="NA")))

p1<-plot_richness(physeqCEILING,x="Latitude",color = "Fraction3",measures=c("Observed"),title = "Vibrio Species Richness")

p1$data$Depth<- as.character(p1$data$Depth)
p1$data$Depth <- factor(p1$data$Depth, levels=newSTorder)


#latitudine 

min(p1$data$Latitude)
max(p1$data$Latitude)

ggplot(p1$data,aes(Latitude,value))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+  geom_point(aes(colour = Fraction3),alpha =0.45)+
  scale_color_manual(values=c("#020bf5", "#b80009"))+facet_grid(~Depth)+
  geom_ysideboxplot(aes(x=Fraction3,y=value,colour = Fraction3), orientation = "x") +
  theme(        ggside.panel.scale.y = .4)+scale_ysidex_discrete()+
  geom_smooth(aes(colour = Fraction3),size=1.5)+ ylab('Richness') + labs(color='Fraction')  +
  scale_x_continuous(breaks = seq(-65, 80, by = 20))#, expand = c(0, 0)  





#####################################################################################################################
###    FIG S2 t-PCoA
#####################################################################################################################

############merge average per zona


variable1 = as.character(get_variable(physeq, "Zone"))
variable2 = as.character(get_variable(physeq, "Fraction3"))
# variable3 = as.character(get_variable(physeq, "Distance1"))

sample_data(physeq)$NewPastedVar <- mapply(paste, variable1, variable2  
                                           , sep = "_")
sample_data(physeq)
# write.csv(sample_data(physeq),"sample_data_physeq.csv")#ho messo i nomi come quella dei kmers
# sample_data(physeq)<-sample_data(read.csv("sample_data_physeq.csv",row.names=1, check.names=FALSE,sep = ";"))


physeq_MERGED_zone_FRACTION<-merge_samples_mean(physeq, "NewPastedVar")#cambiato er fare la mantel, per fare b div ricalcolarlo DHN 







sample_variables(physeq_MERGED_zone_FRACTION)
sample_data(physeq_MERGED_zone_FRACTION)
sample_names(physeq_MERGED_zone_FRACTION)


#per far quello che facevo su excel
rm(df)
RR<-nrow(as.data.frame(sample_names(physeq_MERGED_zone_FRACTION)))
df <- data.frame(matrix(ncol = 1, nrow = RR))
df$names<-(as.data.frame(sample_names(physeq_MERGED_zone_FRACTION)))

df<-df[,2]
colnames(df)[1] <- "Samples"

foo <- data.frame(do.call('rbind', strsplit(as.character(df$Samples),'_',fixed=TRUE)))
df$Zone <-foo$X1
df$Fraction <-foo$X2
colnames(df)

rownames(df) <- df[,1]

# df$Ocean<-c(rep("Atlantic",9),rep("Polar",5),rep("Atlantic",8),rep("Indian",9),rep("Sea",5),rep("Pacific",12),
#            rep("Sea",4),rep("Polar",4))
df$Ocean<-c(rep("Atlantic",4),rep("Polar",2),rep("Atlantic",4),rep("Indian",4),rep("Sea",2),rep("Pacific",6),
            rep("Sea",2),rep("Polar",2))

physeq_MERGED_zone_FRACTION
#cambio il metadata




sample_data(physeq_MERGED_zone_FRACTION)<-sample_data(df)
physeq_MERGED_zone_FRACTION

sample_variables(physeq_MERGED_zone_FRACTION)
sample_data(physeq_MERGED_zone_FRACTION)



physeq_MERGED_zone_FRACTION1<-subset_samples(physeq_MERGED_zone_FRACTION, !(  Fraction=="NA" ))
# physeq_MERGED_zone_FRACTION1<-subset_samples(physeq_MERGED_zone_FRACTION, !(  Samples=="IOS_20-180" ))



##beta diversity  

otu.ord <- ordinate(physeq = physeq_MERGED_zone_FRACTION1, "PCoA")
sample_variables(physeq_MERGED_zone_FRACTION1)
#asse 1-2  
library(ggforce)
a<-plot_ordination(physeq = physeq_MERGED_zone_FRACTION1, otu.ord,color = "Fraction",shape = "Ocean",
                   axes =c(1,2))+
  geom_text(aes(label=Zone), Fraction = 3, vjust = 0,hjust=0, show.legend = FALSE)+
  theme(plot.title = element_text(hjust = 0.0))+ geom_point(Fraction = 2)+ theme_void()+theme_bw()
a1<-a+  geom_mark_ellipse(aes(color = Fraction), show.legend = FALSE)+ theme_void()+theme_bw() + geom_xsidedensity(aes(y=stat(density),fill=Fraction), alpha = 0.5, show.legend = FALSE) +
  geom_ysidedensity(aes(x=stat(density),fill=Fraction), alpha = 0.5, show.legend = FALSE) +
  theme_bw() +  scale_xsidey_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
  scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +scale_ysidex_discrete()+
  ggside::theme_ggside_void() 
a1

## PERMANOVA

OTU_matrix = as(otu_table(physeq_MERGED_zone_FRACTION), "matrix") 
OTUdf = as.data.frame(OTU_matrix)
BC_TAX<-vegdist(t(OTUdf),method = "bray")
env_data = as.data.frame(sample_data(physeq_MERGED_zone_FRACTION))
adonis2(BC_TAX~env_data$Fraction*env_data$Ocean, permutations = 9999)

## Mantel between bray taxonomy and kmers matrices 
BC_TAX<-vegdist(t(OTUdf),method = "bray")
BC_kmers=as.matrix(read.table("mat_abundance_braycurtis.csv",sep=";", header=TRUE, row.names=1))
BC_kmers[upper.tri(BC_kmers)] <- 0
BC_kmersDist <- as.dist(BC_kmers, diag = TRUE)
df <- data.frame( Kmers=BC_kmersDist[lower.tri(BC_kmersDist)], TAX=as.dist(BC_TAX)[lower.tri(as.dist(BC_TAX))])
set.seed(1234)
mantel(xdis =BC_TAX ,ydis = BC_kmersDist,method = "pearson",permutations = 9999)



## Plot Fig 1D k-PCoA
iMDS  <- ordinate(physeq_MERGED_zone_FRACTION, "PCoA", distance=BC_kmersDist)  
 physeq_MERGED_zone_FRACTION_kmears<-physeq_MERGED_zone_FRACTION
 sample_names(physeq_MERGED_zone_FRACTION_kmears)
 
 
 ### cambio i nomi
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("*_FLB$", "_BACT",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("*_PAB$", "_PROT", sample_names(physeq_MERGED_zone_FRACTION_kmears))
 
 
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("ANE_(\\w+)", "\\1_ANE",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("ANW_(\\w+)", "\\1_ANW",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("ARC_(\\w+)", "\\1_ARC",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("ASE_(\\w+)", "\\1_ASE",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("ASW_(\\w+)", "\\1_ASW",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("ION_(\\w+)", "\\1_ION",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("IOS_(\\w+)", "\\1_IOS",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("MED_(\\w+)", "\\1_MED",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("PON_(\\w+)", "\\1_PON",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("PSE_(\\w+)", "\\1_PSE",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("PSW_(\\w+)", "\\1_PSW",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("RED_(\\w+)", "\\1_RED",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
  sample_names(physeq_MERGED_zone_FRACTION_kmears) <- gsub("SOC_(\\w+)", "\\1_SOC",  sample_names(physeq_MERGED_zone_FRACTION_kmears))
 
  sample_data(physeq_MERGED_zone_FRACTION_kmears)<-sample_data(read.csv("sample_data_physeq_MERGED_zone_FRACTION.csv",row.names=1, check.names=FALSE,sep = ";"))
  
  
 
 
plot_ordination(physeq_MERGED_zone_FRACTION_kmears, iMDS,color = "Ocean",shape = "Ocean",
                   axes =c(1,2))+
  geom_text(aes(label=Zone), size = 3, vjust = 0,hjust=0, show.legend = FALSE)+
  theme(plot.title = element_text(hjust = 0.0))+ geom_point(size = 2)+ theme_void()+theme_bw()+
  geom_mark_ellipse(aes(color = Fraction ), show.legend = FALSE)+ theme_void()+theme_bw() + 
  geom_xsidedensity(aes(y=stat(density),fill=Fraction), alpha = 0.5, show.legend = FALSE) +
  geom_ysidedensity(aes(x=stat(density),fill=Ocean), alpha = 0.5, show.legend = FALSE) +
  theme_bw() +  scale_xsidey_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
  scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +scale_ysidex_discrete()+
  ggside::theme_ggside_void() 

 





####################################################################
# 7. Extract reads classified as Vibrio using Enterobase classified reads
####################################################################

cat $SAMPLE | parallel -j 10 \
"python3 ~/extract_kraken_reads.py \
-k ${ENT_OUT}/{}_output_entero.kraken \
-s1 ${REFSEQ_OUT}/{}_extracted_reads_vibrio_1.fa \
-s2 ${REFSEQ_OUT}/{}_extracted_reads_vibrio_2.fa \
--taxid 662 --include-children \
-o ${ENT_OUT}/{}_extracted_reads_vibrio_entero_1.fa \
-o2 ${ENT_OUT}/{}_extracted_reads_vibrio_entero_2.fa"


##  simka all samples 

simka -in imput_simka_PAV.txt -out results_PAV -kmer-size 31  -max-merge 4 -max-reads 0 -min-shannon-index 1.5
simka -in imput_simka_FLV.txt -out results_FLV -kmer-size 31  -max-merge 4 -max-reads 0 -min-shannon-index 1.5



## simka only superficial samples 
simka -in input_sinka_derep_SRF.txt -out simka_SRF -kmer-size 31 -max-reads 0 -min-shannon-index 1.5







#co-assembly per zone

cat sets.txt
ARC
ANE
ANW
ASE
ASW
ION
IOS
MED
PON
PSE
PSW
RED
SOC

for SET in `cat sets.txt`
do
R1s=`ls *_1.fa | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
R2s=`ls *_2.fa | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
 megahit -1 $R1s -2 $R2s -o $SET-co-assembly.fa 
done



# CAT TAXONOMY sui contigs 

for i in *.fa
do
CAT contigs -c $i -d /home/userbio/CAT_prepare_20210107/2021-01-07_CAT_database -t /home/userbio/CAT_prepare_20210107/2021-01-07_taxonomy -o BAT_output/${i}_BAT
CAT add_names -i BAT_output/${i}_BAT.contig2classification.txt -o BAT_off/${i}_.names_off.txt -t /home/userbio/CAT_prepare_20210107/2021-01-07_taxonomy --only_official
CAT summarise -c $i -i BAT_off/${i}_.names_off.txt -o BAT_class/${i}_BAT_summ.txt
done


#quantification on contigs with salmon

for i in *.fa
do
    prefix=$(basename $i .fa)
    salmon index -t $i -i ${prefix}_index 
    mv ${prefix}_index salmon
done

for i in *BACT*_index
do
prefix=$(basename $i _index)
sample=$(echo "$prefix" | cut -d "_" -f2)
echo $sample
salmon quant --meta -l A --index $i \
-1 ${sample}_ALL_READS_1.fa \
-2 ${sample}_ALL_READS_2.fa \
 -o quantification/${prefix}
done



for i in *PROT*_index
do
prefix=$(basename $i _index)
sample=$(echo "$prefix" | cut -d "_" -f2)
echo $sample
salmon quant --meta -l A --index $i  \
-1 ${sample}_ALL_READS_1.fa \
-2 ${sample}_ALL_READS_2.fa \
-o quantification/${prefix}
done





