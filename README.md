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


##  simka all samples --> used for the k-PCoA (fig 1D see line 599)

simka -in imput_simka.txt -out results -kmer-size 31  -max-merge 4 -max-reads 0 -min-shannon-index 1.5




## simka only superficial samples -> used for bioregions (Figure 2 A) and Vibrio surface dispersion through Oceans (fig 3B to G) 
simka -in input_sinka_derep_SRF.txt -out simka_SRF -kmer-size 31 -max-reads 0 -min-shannon-index 1.5



#####
####   BIOREGIONS
####
library(scales)

file_paths <- c(
  "SRF_0.22-3_mat_abundance_braycurtis.csv",
  "SRF_3-2_mat_abundance_braycurtis.csv",
  "SRF_20-180_mat_abundance_braycurtis.csv",
  "SRF_180_2000_mat_abundance_braycurtis.csv"
)

# Lettura dei file CSV in matrici
matrices <- lapply(file_paths, read.csv, row.names = 1, sep = ";")
matrices






# Esecuzione della PCoA
pcoa_results <- lapply(matrices, function(matrix) {
  dist_matrix <- as.dist(matrix) # Assicurati che sia una matrice di distanza
  cmdscale(dist_matrix, eig = TRUE, k = 3, add = TRUE)  
})

# Nomi per i risultati PCoA
names(pcoa_results) <- c("0.22-3", "3-20", "20-180", "180-2000")

# Calcolo dei valori RGB per ogni set di risultati PCoA
rgb_results <- lapply(pcoa_results, function(result) {
  coordinates <- result$points[, 1:3] # Primi tre assi principali
  eigenvalues <- result$eig[1:3] # Primi tre autovalori
  
  lambda_r <- 1 # Rapporto per il canale rosso
  lambda_g <- eigenvalues[2] / eigenvalues[1] # Rapporto per il canale verde
  lambda_b <- eigenvalues[3] / eigenvalues[1] # Rapporto per il canale blu
  
  # Conversione delle coordinate in valori RGB
  convert_to_rgb <- function(coord, lambda) {
    128 * (1 + lambda * coord / max(abs(coord)))
  }
  
  r <- convert_to_rgb(coordinates[, 1], lambda_r)
  g <- convert_to_rgb(coordinates[, 2], lambda_g)
  b <- convert_to_rgb(coordinates[, 3], lambda_b)
  
  r <- pmin(pmax(r, 0), 255)
  g <- pmin(pmax(g, 0), 255)
  b <- pmin(pmax(b, 0), 255)
  
  rgb_values <- data.frame(r = round(r), g = round(g), b = round(b))
  return(rgb_values)
})

# Lettura delle coordinate delle stazioni
coordinate <- read.csv2("cytoscape_samples_coordinates.csv", sep = ";")

# Loop per ciascuno dei set di dati
for (set_name in names(rgb_results)) {
  
  # Unione dei risultati RGB con le coordinate
  rgb_df <- data.frame(Sample = rownames(rgb_results[[set_name]]), rgb_results[[set_name]], stringsAsFactors = FALSE)
  merged_data <- merge(coordinate, rgb_df, by = "Sample")
  
  # Assicurati che r, g, b siano numerici e non fattori o caratteri
  merged_data$r <- as.numeric(as.character(merged_data$r))
  merged_data$g <- as.numeric(as.character(merged_data$g))
  merged_data$b <- as.numeric(as.character(merged_data$b))
  
  # Funzione per riscalare i valori in una gamma specifica (0-255)
  rescale_to_255 <- function(x) {
    (x - min(x)) / (max(x) - min(x)) * 255
  }
  
  # Riscalare i valori RGB
  normalize_rgb <- function(rgb_df) {
    rgb_df$r <- rescale_to_255(rgb_df$r)
    rgb_df$g <- rescale_to_255(rgb_df$g)
    rgb_df$b <- rescale_to_255(rgb_df$b)
    return(rgb_df)
  }
  
  # Applicare la normalizzazione ai valori RGB
  normalized_rgb_df <- normalize_rgb(rgb_df)
  
  # Unire i risultati normalizzati con le coordinate
  merged_data <- merge(coordinate, normalized_rgb_df, by = "Sample")
  
  # Conversione dei valori RGB in colori esadecimali
  merged_data$color <- apply(merged_data[, c("r", "g", "b")], 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
  
  #FIG2A
  ggplot() +
    borders("world", colour = "gray50", fill = "gray50") + # Aggiunge i confini del mondo
    geom_point(data = merged_data, aes(x = as.numeric(Longitude), y = as.numeric(Latitude), color = color), size = 3) +
    scale_color_identity() + # Usa i colori come forniti
    theme_minimal() +
    labs(x = "Longitudine", y = "Latitudine") +
    coord_fixed(1.3) + # Mantiene le proporzioni corrette
    ggtitle(paste("Mappa per il set", set_name))
  
  # Visualizzazione del grafico PCoA con colori normalizzati
  pcoa_coords <- as.data.frame(pcoa_results[[set_name]]$points[, 1:2])
  colnames(pcoa_coords) <- c("PC1", "PC2")
  
  merged_data_pcoa_rgb <- cbind(pcoa_coords, normalized_rgb_df)
  merged_data_pcoa_rgb$color <- apply(merged_data_pcoa_rgb[, c("r", "g", "b")], 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
  merged_data_pcoa_rgb$sample <- rownames(merged_data_pcoa_rgb)
  
  ggplot(merged_data_pcoa_rgb, aes(x = PC1, y = PC2, color = color)) +
    geom_point() +
    geom_text(aes(label = sample), vjust = -1, hjust = 0.5) +  
    scale_color_identity() +
    theme_minimal() +
    labs(x = "PC1", y = "PC2") +
    ggtitle(paste("PCoA per il set", set_name))

  c <- nlcor(merged_data_pcoa_rgb$PC1, merged_data_pcoa_rgb$PC2)
  
  # Creazione di un data frame e aggiunta alla lista
  df_new <- data.frame(Fraction = set_name, Value = c$cor.estimate)
  df_list[[set_name]] <- df_new
}


df_combined <- bind_rows(df_list)
## fig 2B
# Creazione del grafico unico con tutti i range
ggplot(df_combined, aes(x = Fraction, y = Value, fill = Fraction)) +
  geom_point(aes(color = Fraction)) +
  geom_line(aes(group = Fraction)) +
  theme_minimal() +
  labs(title = "NLCC across all fractions", x = "Fraction", y = "NLCC") +
  scale_fill_brewer(palette = "Set3")


















##################################################
############# cumulative correlarions############# 
##################################################
library(ggplot2)
library(plyr )
files <- c("cytoscape/file_conx/SRF_0.22-3_TimTrav.csv",
           "cytoscape/file_conx/SRF_5-20_TimTrav.csv",
           "cytoscape/file_conx/SRF_20-180_TimTrav.csv",
           "cytoscape/file_conx/SRF_180-2000_TimTrav.csv")

computeConfidence <- function(rho,N){
  upper <- tanh(atanh((rho+1.96)/sqrt(N-3)+0i))
  lower <- tanh(atanh((rho-1.96)/sqrt(N-3)+0i))
  return(c(Re(upper),Re(lower)))
}



cumCorr <- function(file_){
  CC<-read.csv(file_,sep = ",")
  CC<- CC[CC$tempo.medio.anni<=10,]
  CC2 <- CC[order(CC$tempo.medio.anni),]
  DB <- data.frame()
  for(i in 4:dim(CC2)[1]) {
    #print(i)
    CC_1<-CC2[1:i,]
    #print(CC_1)
    cc <- cor.test(as.numeric(CC_1$tempo.medio.anni),as.numeric(CC_1$dissimilarity),method = "spearman")$estimate
    ci <- computeConfidence(cc,i)
    DB <- rbind(DB,c(CC2[i,11],cc))
  }
  return(DB)
}

DBs <- sapply(files,cumCorr)
a <- data.frame(DBs[[1]], DBs[[2]], "0.22-3")
b <- data.frame(DBs[[3]], DBs[[4]], "5-20")
c <- data.frame(DBs[[5]], DBs[[6]], "20-180")
d <- data.frame(DBs[[7]], DBs[[8]], "180-2000")

# Creazione di una lista contenente i data frame
dataframes <- list(a, b, c, d)

# Loop per applicare le operazioni su tutti i data frame
for (i in 1:length(dataframes)) {
  colnames(dataframes[[i]])[1] <- "years.trav.time"
  colnames(dataframes[[i]])[2] <- "R"
  colnames(dataframes[[i]])[3] <- "Fraction"
}

# Assegnazione dei data frame modificati ad a, b, c, d
a <- dataframes[[1]]
b <- dataframes[[2]]
c <- dataframes[[3]]
d <- dataframes[[4]]
allDF <- rbind(a, b,c,d)
colnames(allDF)[1]<-"years.trav.time"
newSTorder = c( "0.22-3","5-20","20-180" ,"180-2000")

allDF$Fraction<- as.character(allDF$Fraction)
allDF$Fraction <- factor(allDF$Fraction, levels=newSTorder) 

#rangesssss
allDF$years.trav.time_ranges<-round(allDF$years.trav.time+0.5)
 
allDF$years.trav.time_ranges2<-round_any(allDF$years.trav.time, 0.5, f=ceiling)
 

#### fig 3B
ggplot(allDF, aes(x = years.trav.time, y = R, color = Fraction)) +
    geom_point(alpha = 0.1) +#, color = "black"
    geom_smooth(method = "gam", se = FALSE, formula = y ~ s(x, bs = "cs")) +  # Metodo GAM
    scale_y_continuous(limit = c(-0.25, 0.75)) +
    scale_x_continuous(limit = c(0.75, 5), breaks = round(seq(min(allDF$years.trav.time),
                                                               max(allDF$years.trav.time),
                                                               by = 0.5), 1)) +
    labs(color = "")+
  geom_smooth(method = "gam", se = FALSE, formula = y ~ s(x, bs = "cs"), color = "darkred", linetype = "dashed", size = 1.5) +
    theme(panel.background = element_blank(),  # Background bianco
          panel.grid.major = element_blank(),  # Griglia grigia
          panel.grid.minor = element_blank())  # Nascondi griglia minore











#########################################################
###############correlazioni simi and trav time
#######################################################





library(reshape2)
library(geodist)
library(ggplot2)
library(ggpubr)
library(diffcor)


coordinate<-read.csv2("cytoscape_samples_coordinates.csv",sep = ";")
color_palette <- c("SRF_0.22-3" = "red", "SRF_5-20" = "blue", "SRF_20-180" = "green", "SRF_180-2000" = "purple")
fraction_levels <- c("SRF_0.22-3", "SRF_5-20", "SRF_20-180", "SRF_180-2000")

list2<-data.frame(coordinate$Longitude,coordinate$Latitude,coordinate$Sample)

colnames(list2)
colnames(list2)[1] ="longitude"
colnames(list2)[2] ="latitude"
colnames(list2)[3] ="name"

distance_matrix <- geodist(list2, measure = 'geodesic' )/1000 #converting it to km

#also, check for other measures in the description

colnames(distance_matrix) <- list2$name
rownames(distance_matrix) <- list2$name

head.matrix(distance_matrix)



mat=as.matrix(distance_matrix)
melted<-reshape2::melt(mat)
newDF<-melted
 

newDF$merged<-paste(newDF$Var1,newDF$Var2,sep = "_")
colnames(newDF)[4] ="conx"
colnames(newDF)[3] ="dist_km"
head(newDF)

file_paths <- c(
  "cytoscape/file_conx/SRF_0.22-3_TimTrav.csv",
  "cytoscape/file_conx/SRF_5-20_TimTrav.csv",
  "cytoscape/file_conx/SRF_20-180_TimTrav.csv",
  "cytoscape/file_conx/SRF_180-2000_TimTrav.csv"
)

#head(read.csv("cytoscape/file_conx/SRF_0.22-3_TimTrav.csv", sep = ","))
##### loop per legegre i dati e filtarre <  1.5 !!!

correlation_data <- list()


for (file_path in file_paths) {
  cc <- read.csv(file_path, sep = ",")
  name <- gsub("cytoscape/file_conx/|_TimTrav.csv", "", file_path)
  cc$fraction <- rep(name, length(cc$conx))
  merged_data <- merge(cc, newDF, by = "conx", all.x = TRUE)
  merged_data <- merged_data[merged_data$tempo.medio.anni < 1.5,]
  correlation_data[[name]] <- merged_data
}


# fig 3C

plot <- ggplot()
for (fraction in names(correlation_data)) {
  df <- correlation_data[[fraction]]
  df$fraction <- factor(reorder(df$fraction, -df$tempo.medio.anni), levels = fraction_levels)  # Definizione dei livelli della variabile "fraction" con ordine inverso della dist_km
  plot <- plot + geom_smooth(data = df, aes(y = similarity, x = tempo.medio.anni, color = fraction), method = "lm") +
    geom_smooth(data = df, aes(y = similarity, x = tempo.medio.anni), method = "lm", color = color_palette[fraction], fill = color_palette[fraction], alpha = 0.2)
}

# Ordinamento della legenda
plot <- plot + scale_color_manual(values = color_palette, guide = guide_legend(order = 1))

# Visualizzazione del plot
plot
#Fisher's z-Tests Concerning Difference of Correlations


library(diffcor)

# Creazione del dataframe vuoto
df <- data.frame(fraction = character(), correlation = numeric(), p_value = numeric(), diff_corr = numeric(), stringsAsFactors = FALSE)

# Calcolo delle correlazioni e p-value per ogni frazione
for (fraction in names(correlation_data)) {
  data <- correlation_data[[fraction]]
  
  # Calcolo della correlazione
  r <- round(cor(data$similarity, data$tempo.medio.anni), 2)
  
  # Calcolo del p-value
  p <- cor.test(data$similarity, data$tempo.medio.anni)$p.value
  
  # Aggiunta delle informazioni al dataframe
  df <- rbind(df, data.frame(fraction = fraction, correlation = r, p_value = p, length=nrow(data)))
}

# Visualizzazione del dataframe
df
diff_df <- data.frame(fraction1 = character(), fraction2 = character(), diff_corr = numeric(), stringsAsFactors = FALSE)

# Calcolo delle differenze in correlazione per ogni combinazione di frazioni
for (i in 1:(nrow(df) - 1)) {
  for (j in (i + 1):nrow(df)) {
    fraction1 <- df$fraction[i]
    fraction2 <- df$fraction[j]
    correlation1 <- df$correlation[i]
    correlation2 <- df$correlation[j]
    length1 <- df$length[i]
    length2 <- df$length[j]
    
    # Calcolo della differenza in correlazione
    diff_corr <- diffcor.two(correlation1, correlation2, length1, length2, digit = 3)
    
    # Aggiunta delle informazioni al dataframe delle differenze in correlazione
    diff_df <- rbind(diff_df, data.frame(fraction1 = fraction1, fraction2 = fraction2, diff_corr))
  }
}

# Visualizzazione del dataframe delle differenze in correlazione
diff_df


#########################################################
###############correlazioni simi e km 
#######################################################


library(raster)
library(maptools)
library(ggplot2)
library(dplyr)
library(sf)
library(ggpubr)
library(reshape2)
library(geodist)
library(ggplot2)
library(ggpubr)
library(diffcor)








#####################  preparazione dei dati

coordinate<-read.csv2("cytoscape_samples_coordinates.csv",sep = ";")
head(coordinate)



rm(list2)
list2<-data.frame(coordinate$Longitude,coordinate$Latitude,coordinate$Sample)

colnames(list2)
colnames(list2)[1] ="longitude"
colnames(list2)[2] ="latitude"
colnames(list2)[3] ="name"

distance_matrix <- geodist(list2, measure = 'geodesic' )/1000 #converting it to km

#also, check for other measures in the description

colnames(distance_matrix) <- list2$name
rownames(distance_matrix) <- list2$name

 


mat=as.matrix(distance_matrix)
melted<-reshape2::melt(mat)

newDF<-melted
 
 

newDF$merged<-paste(newDF$Var1,newDF$Var2,sep = "_")
colnames(newDF)[4] ="conx"
colnames(newDF)[3] ="dist_km"
head(newDF)





# Leggi i dati delle correlazioni da diversi file CSV
file_paths <- c(
  "cytoscape/file_conx/SRF_0.22-3_TimTrav.csv",
  "cytoscape/file_conx/SRF_5-20_TimTrav.csv",
  "cytoscape/file_conx/SRF_20-180_TimTrav.csv",
  "cytoscape/file_conx/SRF_180-2000_TimTrav.csv"
)



color_palette <- c("SRF_0.22-3" = "red", "SRF_5-20" = "blue", "SRF_20-180" = "green", "SRF_180-2000" = "purple")


# Definizione dell'ordine dei livelli della variabile "fraction"
fraction_levels <- c("SRF_0.22-3", "SRF_5-20", "SRF_20-180", "SRF_180-2000")




correlation_data <- list()

for (file_path in file_paths) {
  cc <- read.csv(file_path, sep = ",")
  name <- gsub("cytoscape/file_conx/|_TimTrav.csv", "", file_path)
  cc$fraction <- rep(name, length(cc$conx))
  merged_data <- merge(cc, newDF, by = "conx", all.x = TRUE)
 # merged_data <- merged_data[merged_data$dist_km <= 5000,]
  correlation_data[[name]] <- merged_data
}



# Caricare la libreria dplyr
library(dplyr)

# Unire tutti i data frame in un unico data frame
all_data <- bind_rows(correlation_data)
dim(all_data)
# Selezionare solo le righe con 'conx' unici
unique_data <- all_data %>% distinct(conx, .keep_all = TRUE)
dim(unique_data)
head(unique_data)

 
unique_data<-unique_data[, !(names(unique_data) %in% c("X1", "X", "Var1.x" , "Var2.x","Var1.y"  ,"Var2.y"))]

# Salvare il data frame come CSV
write.csv(unique_data, "unique_conx_data.csv", row.names = FALSE)

import pandas as pd
import searoute as sr

# Leggi il file CSV
df = pd.read_csv('unique_conx_data.csv')

# Lista per memorizzare le distanze calcolate
distances = []

# Calcola la distanza tra le coppie di stazioni
for index, row in df.iterrows():
    origin = [row['Longitude_start'], row['Latitude_start']]
    destination = [row['Longitude_end'], row['Latitude_end']]
    try:
        route = sr.searoute(origin, destination)
        distance_km = route['properties']['length']
    except Exception as e:
        print(f"Errore nel calcolo della distanza per la riga {index}: {e}")
        distance_km = None
    distances.append(distance_km)

# Aggiungi la colonna con le distanze calcolate
df['searoutekm'] = distances

# Salva il risultato in un nuovo file CSV
df.to_csv('unique_conx_data_with_distances.csv', index=False)


### script fatto basandosi su searoutes nella cartella  km_onlywaters

conx_w_coord<-read.csv("km_onlywaters/unique_conx_data_with_distances.csv")

# Definisci una funzione per fare il merge per ogni data frame nella lista
merge_and_filter <- function(df, conx_w_coord) {
  df <- merge(df, conx_w_coord[, c("conx", "searoutekm")], by = "conx", all.x = TRUE)
  df <- df %>% filter(searoutekm <= 5000)
  return(df)
}

# Applica la funzione a ogni data frame nella lista e salva il risultato in correlation_data_SR
correlation_data_SR <- lapply(correlation_data, merge_and_filter, conx_w_coord)














#########################################################
###############correlazioni comparison log e zscore
#######################################################

# Funzione per calcolare lo Z-score
zscore <- function(x) {
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

# Lettura e filtro dei dati di correlazione per distanza e anni
correlation_data_combined <- data.frame()

for (file_path in file_paths) {
  cc <- read.csv(file_path, sep = ",")
  name <- gsub("cytoscape/file_conx/|_TimTrav.csv", "", file_path)
  cc$fraction <- rep(name, length(cc$conx))
  merged_data <- merge(cc, newDF, by = "conx", all.x = TRUE)
  merged_data <- merged_data[merged_data$dist_km <= 5000 & merged_data$tempo.medio.anni < 1.5,]
  merged_data$log_dist_km <- log1p(merged_data$dist_km)
  merged_data$log_tempo_anni <- log1p(merged_data$tempo.medio.anni)
  merged_data$z_log_dist_km <- zscore(merged_data$log_dist_km)
  merged_data$z_log_tempo_anni <- zscore(merged_data$log_tempo_anni)
  merged_data_long <- melt(merged_data, id.vars = c("conx", "similarity", "fraction"), measure.vars = c("z_log_dist_km", "z_log_tempo_anni"), variable.name = "type", value.name = "z_log_value")
  correlation_data_combined <- rbind(correlation_data_combined, merged_data_long)
}

# Creazione del plot combinato
plot <- ggplot(correlation_data_combined, aes(x = z_log_value, y = similarity, color = fraction, linetype = type)) +
  geom_smooth(method = "lm") +
  stat_cor(aes(label = paste(..r.label.., sep = "")), method = "pearson", geom = "text", position = position_jitter(width = 0.2, height = 0), show.legend = FALSE) +
  scale_linetype_manual(values = c("z_log_dist_km" = "dashed", "z_log_tempo_anni" = "solid"), labels = c("z_log_dist_km" = "Distance (z_log_km)", "z_log_tempo_anni" = "Time (z_log_years)")) +
  scale_color_manual(values = color_palette, guide = guide_legend(order = 1)) +
  labs(x = "Z-score Log-transformed Scale (z_log_km and z_log_years)", y = "Similarity") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(vjust = -0.5),
    legend.position = "top"
  )





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

#fig 4 A and B
library(phyloseq)
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
library(dplyr)
library(microbiomeMarker)
library("MiscMetabar")

abund_table<-read.csv("salmon/vibrio_tpm_frequency_table.csv",row.names=1, check.names=FALSE,sep = ",")
dfTAX<-read.csv("salmon/TAX_vibrio_tpm_frequency_table.csv",row.names=1, check.names=FALSE,sep = ",")
metaDF<-read.csv("salmon/metaTable_contigs.csv",row.names=1, check.names=FALSE,sep = ",")
 



physeq<-merge_phyloseq(phyloseq(otu_table(as.matrix(abund_table), taxa_are_rows = T), tax_table(as.matrix(dfTAX)), sample_data(metaDF)))
genus.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "specie"], sum, na.rm=TRUE)

top5phyla = names(sort(genus.sum, TRUE))[1:20]
physeq25 = prune_taxa((tax_table(physeq)[, "specie"] %in% top5phyla), physeq)

physeq_norm<-normalize(physeq25, method = "RLE")
# Define your phyloseq object
physeq <- physeq_norm   # replace with your phyloseq object
# Define your phyloseq object

# Ensure physeq is of class phyloseq
if (!inherits(physeq, "phyloseq")) {
  stop("physeq must be an object of class 'phyloseq'")
}

# Transpose OTU table if necessary
if (!physeq@otu_table@taxa_are_rows) {
  otu_tab <- t(physeq@otu_table)
} else {
  otu_tab <- physeq@otu_table
}

otu_tab <- as.data.frame(as(otu_tab, "matrix"))

# Summarize OTU table by taxa
tax_table <- as.data.frame(as(physeq@tax_table, "matrix"))
tax_table$Taxon <- rownames(tax_table)

# Merge OTU table with taxonomy table
otu_taxa <- merge(otu_tab, tax_table, by = "row.names")
rownames(otu_taxa) <- otu_taxa$Row.names
otu_taxa$Row.names <- NULL

# Specify the taxa level to use (e.g., "species")
taxa <- "species"  # replace with the desired taxonomic level, if different

# Summarize OTU counts by specified taxa rank and compartment (BACT or PROT)
taxa_counts <- otu_taxa %>%
  dplyr::group_by(across(all_of("species"))) %>%
  dplyr::summarise(across(starts_with("BACT_") | starts_with("PROT_"), sum)) %>%
  pivot_longer(cols = starts_with("BACT_") | starts_with("PROT_"),
               names_to = "Sample", values_to = "Abundance") %>%
  dplyr::mutate(Compartment = ifelse(grepl("^BACT_", Sample), "BACT", "PROT"))

 
taxa_counts <- taxa_counts %>%
  group_by(!!sym(taxa)) %>%
  mutate(Total = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Normalized_Abundance = Abundance / Total)



taxa_order <- taxa_counts %>%
  filter(Compartment == "PROT") %>%
  arrange(desc(Normalized_Abundance)) %>%
  distinct(!!sym(taxa)) %>%
  pull(!!sym(taxa))

# Convert taxa to a factor with levels ordered by PROT abundance
taxa_counts <- taxa_counts %>%
  mutate(!!sym(taxa) := factor(!!sym(taxa), levels = unique(taxa_order)))

ggplot(taxa_counts, aes(x = !!sym(taxa), y = Normalized_Abundance, fill = Compartment)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Make the bar plot horizontal
  labs(title = "Normalized Abundance of Taxa in BACT and PROT Compartments",
       x = "Taxa",
       y = "Normalized Abundance") +
  theme_minimal() +
  theme(legend.position = "right")


ggplot(taxa_counts, aes(x = !!sym(taxa), y = Normalized_Abundance, fill = Oceano   )) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Make the bar plot horizontal
  labs(title = "Normalized Abundance of Taxa in BACT and PROT Compartments by Ocean",
       x = "Taxa",
       y = "Normalized Abundance") +
  theme_minimal() +
  theme(legend.position = "right")


### distribution oceanic regions



# Ensure physeq is of class phyloseq
if (!inherits(physeq, "phyloseq")) {
  stop("physeq must be an object of class 'phyloseq'")
}

# Transpose OTU table if necessary
if (!physeq@otu_table@taxa_are_rows) {
  otu_tab <- t(physeq@otu_table)
} else {
  otu_tab <- physeq@otu_table
}

otu_tab <- as.data.frame(as(otu_tab, "matrix"))

# Summarize OTU table by taxa
tax_table <- as.data.frame(as(physeq@tax_table, "matrix"))
tax_table$Taxon <- rownames(tax_table)

# Merge OTU table with taxonomy table
otu_taxa <- merge(otu_tab, tax_table, by = "row.names")
rownames(otu_taxa) <- otu_taxa$Row.names
otu_taxa$Row.names <- NULL

# Specify the taxa level to use (e.g., "species")
taxa <- "species"  # replace with the desired taxonomic level, if different

# Add metadata information
metaDF <- data.frame(
  'Campione' = c(
    "BACT_ANE", "BACT_ARC", "BACT_ANW", "BACT_ASE", "BACT_ASW",
    "BACT_ION", "BACT_IOS", "BACT_MED", "BACT_PON", "BACT_PSE",
    "BACT_PSW", "BACT_RED", "BACT_SOC", "PROT_ANE", "PROT_ANW",
    "PROT_ARC", "PROT_ASE", "PROT_ASW", "PROT_ION", "PROT_IOS",
    "PROT_MED", "PROT_PON", "PROT_PSE", "PROT_PSW", "PROT_RED",
    "PROT_SOC"
  ),
  'Fraction' = c(
    rep("Bacteria", 13),
    rep("Protist", 13)
  ),
  'Oceano' = c(
    rep(c("Atlantic", "Polar", "Atlantic", "Atlantic", "Atlantic",
          "Indian", "Indian", "Seas", "Pacific", "Pacific",
          "Pacific", "Seas", "Polar"), 2)
  )
)

# Merge with taxa_counts
taxa_counts <- otu_taxa %>%
  group_by(across(all_of(taxa))) %>%
  summarise(across(starts_with("BACT_") | starts_with("PROT_"), sum)) %>%
  pivot_longer(cols = starts_with("BACT_") | starts_with("PROT_"),
               names_to = "Sample", values_to = "Abundance") %>%
  mutate(Compartment = ifelse(grepl("^BACT_", Sample), "BACT", "PROT")) %>%
  left_join(metaDF, by = c("Sample" = "Campione"))

# Normalize values so that each taxon has the same total abundance
taxa_counts <- taxa_counts %>%
  group_by(!!sym(taxa)) %>%
  mutate(Total = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Normalized_Abundance = Abundance / Total)

ggplot(taxa_counts, aes(x = !!sym(taxa), y = Normalized_Abundance, fill = Oceano   )) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Make the bar plot horizontal
  labs(title = "Normalized Abundance of Taxa in BACT and PROT Compartments by Ocean",
       x = "Taxa",
       y = "Normalized Abundance") +
  theme_minimal() +
  theme(legend.position = "right")
