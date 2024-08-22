# Deciphering-the-Hidden-Ecology-of-Vibrio-in-the-Oceans

This pipeline processes metagenomic TARA Ocean to analyze the *Vibrio* presence in the oceans

## Steps Overview
### Metagenomic Reads Analysis
1. **Trimming**
2. **Kraken2 Classification with RefSeq**
3. **Bracken Genus-Level Refinement (RefSeq)**  --> Fig 1A,1B
4. **Extract Vibrio Reads (RefSeq)**
5. **Kraken2 Classification with Enterobase**
6. **Bracken Species-Level Refinement (Enterobase)** --> fig 1C, t-PCoA
7. **Extract Vibrio Reads (Enterobase)**
8. **Travel time simulation**
9. **Simka kmers spectrum of Vibrio sequences** --> fig 1D (k-PCoA), 2A/B (Vibrio Bioregions), 3 B/C/D/E/F/G

### Contigs Analysis
10. **Megahit Co-Assembling using Vibrio reads**
11. **CAT Species Taxonomy**
12. **Salmon Quantification** --> fig 4A/B


## Script

```
# Set the current path variable
currpath=$(pwd)

# Set the sample name file
SAMPLE="samples_TARA.txt" #list of downloaded metagenomes

# RefSeq Kraken2 database directory
REFSEQ_KRA_DB="krakenDB_bacteria_refseq"

# Enterobase Kraken2 database directory
ENT_KRA_DB="enterobase_vibrio_KRAKEN2_db"

# Output directories
TRIM_OUT="${currpath}/trimming"
REFSEQ_OUT="${TRIM_OUT}/kraken_refseq"
ENT_OUT="${TRIM_OUT}/kraken_entero/entero_extraction"

``` 
 1. **Trimming**
```
cat $SAMPLE | parallel -j 10 'trim_galore \
--cores 5 --paired {}_1.fastq.gz {}_2.fastq.gz \
--trim-n --illumina --fastqc_args "-outdir ${TRIM_OUT}/quality_control" -o ${TRIM_OUT}'

```
2. **Kraken2 Classification with RefSeq**
```

cat $SAMPLE | parallel -j 10 \
"kraken2 --db ${REFSEQ_KRA_DB} \
--paired ${TRIM_OUT}/{}_1_val_1.fq.gz ${TRIM_OUT}/{}_2_val_2.fq.gz \
--threads 20 \
--use-names \
--gzip-compressed \
--report ${REFSEQ_OUT}/{}_report-kraken_refseq.txt \
--output ${REFSEQ_OUT}/{}_output_refseq.kraken"
```
3. **Bracken Genus-Level Refinement (RefSeq)**  --> Fig 1A,1B
```
cat $SAMPLE | parallel -j 10 \
"bracken -d ${REFSEQ_KRA_DB} \
-i ${REFSEQ_OUT}/{}_report-kraken_refseq.txt \
-o ${REFSEQ_OUT}/{}_bracken_genus_refseq.txt \
-l G"

combine_bracken_outputs.py --files *.txt -o ../braken_all_REFSEQ_prokEprot_merged/braken_all_REFSEQ_prokEprot_merged.csv ```
```

### Figure 1A (bubleplot)

```
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

abund_table<-read.csv("input/braken_all_REFSEQ_prokEprot_merged.csv",row.names=NULL, check.names=FALSE,sep = ";")
abund_table<-aggregate(.~Genus,abund_table,mean)
abund_table<-t(abund_table)
OTU_taxonomy<-read.csv("input/braken_all_REFSEQ_prokEprot_merged_fract_TAX.csv",row.names=1,check.names=FALSE,sep = ";")
meta_table<-read.csv("input/braken_all_REFSEQ_prokEprot_merged_fract_METAWcooRd.csv",row.names=1,check.names=FALSE,sep = ";")

OTU = otu_table(as.matrix(abund_table), taxa_are_rows = T)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)
physeq<-merge_phyloseq(phyloseq(OTU, TAX, SAM))

#Top 15 taxa

genus.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Genus"], sum, na.rm=TRUE)
top15phyla = names(sort(genus.sum, TRUE))[1:15]
GP1 = prune_taxa((tax_table(physeq)[, "Genus"] %in% top15phyla), physeq)

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
```

### Figure 1B (heatmap Vibrio)
```
VIBRIO <- subset_taxa(GP1, Genus=="Vibrio" )
 ```
####  Fraction
 ```

variable1 = as.character(get_variable(VIBRIO, "Zone"))
variable2 = as.character(get_variable(VIBRIO, "Fraction3"))
# variable3 = as.character(get_variable(physeq, "Distance1"))

sample_data(VIBRIO)$NewPastedVar <- mapply(paste, variable1, variable2  
                                           , sep = "_")
physeq_MERGED_zone_FRACTIONS<-merge_samples_mean(VIBRIO, "NewPastedVar")
rm(df)
ttt<-nrow(as.data.frame(sample_names(physeq_MERGED_zone_FRACTIONS)))
df <- data.frame(matrix(ncol = 1, nrow = ttt))
df$names<-(as.data.frame(sample_names(physeq_MERGED_zone_FRACTIONS)))

df<-df[,2]
colnames(df)[1] <- "Samples"

foox <- data.frame(do.call('rbind', strsplit(as.character(df$Samples),'_',fixed=TRUE)))
df$Zone <-foox$X1
df$Size <-foox$X2
rownames(df) <- df[,1]
sample_data(physeq_MERGED_zone_FRACTIONS)<-sample_data(df)
physeq_MERGED_zone_FRACTIONS<-subset_samples(physeq_MERGED_zone_FRACTIONS, !(  Size=="na" ))

data_glom_fract<- psmelt(physeq_MERGED_zone_FRACTIONS)  
data_glom_fract$Genus <- as.character(data_glom_fract$Genus) 
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



 ```
####  Depth
 ```
variable1 = as.character(get_variable(VIBRIO, "Zone"))
variable2 = as.character(get_variable(VIBRIO, "Depth"))
# variable3 = as.character(get_variable(physeq, "Distance1"))

sample_data(VIBRIO)$NewPastedVar <- mapply(paste, variable1, variable2  
                                           , sep = "_")

physeq_MERGED_zone_DEPTH<-merge_samples_mean(VIBRIO, "NewPastedVar")
 
 
rm(df)
 
ttt<-nrow(as.data.frame(sample_names(physeq_MERGED_zone_DEPTH)))
df <- data.frame(matrix(ncol = 1, nrow = ttt))
df$names<-(as.data.frame(sample_names(physeq_MERGED_zone_DEPTH)))
df<-df[,2]
colnames(df)[1] <- "Samples"
foo <- data.frame(do.call('rbind', strsplit(as.character(df$Samples),'_',fixed=TRUE)))
df$Zone <-foo$X1
df$W_Layer <-foo$X2
rownames(df) <- df[,1]
sample_data(physeq_MERGED_zone_DEPTH)<-sample_data(df)
 
#remove mix zzz and na depths
physeq_MERGED_zone_DEPTH<-subset_samples(physeq_MERGED_zone_DEPTH, !( ( W_Layer=="MIX" | W_Layer=="ZZZ"| W_Layer=="NA")))

data_glom_depth<- psmelt(physeq_MERGED_zone_DEPTH)  
data_glom_depth$Genus <- as.character(data_glom_depth$Genus)  
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


 ggplot(data_glom_depth, aes(Zone, W_Layer, fill= Abundance)) + 
  geom_tile()+ theme(axis.text.x = element_text(angle = 90))+scale_fill_gradient( na.value="black")+
  theme(panel.background = element_rect(fill = 'black'),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +scale_fill_viridis_c(option = "H",limits = c(0,4))+ 
    guides(fill = guide_colourbar(barwidth = 40, barheight = 3,direction = "horizontal",title = "Abundance %",title.position = "top",title.hjust = 0.5))+
     theme(legend.position = "bottom" ) 
    theme(axis.text=element_text(colour="black"),axis.title.x = element_blank(),
                                               axis.text.x=element_blank(),axis.ticks.x =element_blank())+scale_y_discrete(expand = c(0, 0))+
    scale_x_discrete(expand = c(0, 0))+labs(y= "Depth")+ theme(text = element_text(size = 20)) 

```


4. **Extract Vibrio Reads (RefSeq)**

```
cat $SAMPLE | parallel -j 10 \
"python3 ~/extract_kraken_reads.py \
-k ${REFSEQ_OUT}/{}_output_refseq.kraken \
-s1 ${TRIM_OUT}/{}_1_val_1.fq.gz \
-s2 ${TRIM_OUT}/{}_2_val_2.fq.gz \
--taxid 662 --include-children \
-o ${REFSEQ_OUT}/{}_extracted_reads_vibrio_1.fa \
-o2 ${REFSEQ_OUT}/{}_extracted_reads_vibrio_2.fa"
```


5. **Kraken2 Classification with Enterobase**
```
cat $SAMPLE | parallel -j 10 \
"kraken2 --db ${ENT_KRA_DB} \
--paired ${REFSEQ_OUT}/{}_extracted_reads_vibrio_1.fa ${REFSEQ_OUT}/{}_extracted_reads_vibrio_2.fa \
--threads 20 \
--use-names \
--report ${ENT_OUT}/{}_report-kraken_entero.txt \
--output ${ENT_OUT}/{}_output_entero.kraken"
```
6. **Bracken Species-Level Refinement (Enterobase)**
```
cat $SAMPLE | parallel -j 10 \
"bracken -d ${ENT_KRA_DB} \
-i ${ENT_OUT}/{}_report-kraken_entero.txt \
-o ${ENT_OUT}/{}_bracken_species_entero.txt \
-l S"
combine_bracken_outputs.py --files *.txt -o ../braken_all_ENTERO_prokEprot_merged/braken_all_ENTERO_prokEprot_merged.csv
```
### Figure 1C (Richness)
```
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


abund_table<-read.csv("input/braken_all_ENTERO_prokEprot_merged.csv",row.names=1, check.names=FALSE,sep = ";")
abund_table<-t(abund_table)
OTU_taxonomy<-read.csv("input/braken_all_ENTERO_prokEprot_merged_fract_INPUT_R_OTU_TAX.csv",row.names=1,check.names=FALSE,sep = ";")
meta_table<-read.csv("input/braken_all_REFSEQ_prokEprot_merged_fract_METAWcooRd.csv",row.names=1,check.names=FALSE,sep = ";")

#Convert the data to phyloseq format
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = T)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)
physeq<-merge_phyloseq(phyloseq(OTU, TAX, SAM))



OTU_ceiling = otu_table(as.matrix(ceiling(abund_table)), taxa_are_rows = T)

physeqCEILING<-merge_phyloseq(phyloseq(OTU_ceiling, TAX, SAM))

physeqCEILING<-subset_samples(physeqCEILING, !(  Fraction=="NA" ))

physeqCEILING<-subset_samples(physeqCEILING, !( ( Depth=="MIX" | Depth=="ZZZ"| Depth=="NA")))

p1<-plot_richness(physeqCEILING,x="Latitude",color = "Fraction3",measures=c("Observed"),title = "Vibrio Species Richness")

p1$data$Depth<- as.character(p1$data$Depth)
p1$data$Depth <- factor(p1$data$Depth, levels=newSTorder)


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
```
### Figure S2 (t-PCoA)
```
 
variable1 = as.character(get_variable(physeq, "Zone"))
variable2 = as.character(get_variable(physeq, "Fraction3"))

sample_data(physeq)$NewPastedVar <- mapply(paste, variable1, variable2  
                                           , sep = "_")

physeq_MERGED_zone_FRACTION<-merge_samples_mean(physeq, "NewPastedVar")#cambiato er fare la mantel, per fare b div ricalcolarlo DHN 
rm(df)
RR<-nrow(as.data.frame(sample_names(physeq_MERGED_zone_FRACTION)))
df <- data.frame(matrix(ncol = 1, nrow = RR))
df$names<-(as.data.frame(sample_names(physeq_MERGED_zone_FRACTION)))
df<-df[,2]
colnames(df)[1] <- "Samples"
foo <- data.frame(do.call('rbind', strsplit(as.character(df$Samples),'_',fixed=TRUE)))
df$Zone <-foo$X1
df$Fraction <-foo$X2
rownames(df) <- df[,1]
df$Ocean<-c(rep("Atlantic",4),rep("Polar",2),rep("Atlantic",4),rep("Indian",4),rep("Sea",2),rep("Pacific",6),
            rep("Sea",2),rep("Polar",2))

sample_data(physeq_MERGED_zone_FRACTION)<-sample_data(df)
physeq_MERGED_zone_FRACTION1<-subset_samples(physeq_MERGED_zone_FRACTION, !(  Fraction=="NA" ))

otu.ord <- ordinate(physeq = physeq_MERGED_zone_FRACTION1, "PCoA")


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
```
7. **Extract Vibrio Reads (Enterobase)**
```
cat $SAMPLE | parallel -j 10 \
"python3 ~/extract_kraken_reads.py \
-k ${ENT_OUT}/{}_output_entero.kraken \
-s1 ${REFSEQ_OUT}/{}_extracted_reads_vibrio_1.fa \
-s2 ${REFSEQ_OUT}/{}_extracted_reads_vibrio_2.fa \
--taxid 662 --include-children \
-o ${ENT_OUT}/{}_extracted_reads_vibrio_entero_1.fa \
-o2 ${ENT_OUT}/{}_extracted_reads_vibrio_entero_2.fa"
```
9. **Simka kmers spectrum of Vibrio sequences**

####  simka all samples --> used for the k-PCoA
```
simka -in imput_simka.txt -out results -kmer-size 31  -max-merge 4 -max-reads 0 -min-shannon-index 1.5
```
### Mantel between bray taxonomy and kmers matrices 
```
BC_TAX<-vegdist(t(OTUdf),method = "bray")
BC_kmers=as.matrix(read.table("mat_abundance_braycurtis.csv",sep=";", header=TRUE, row.names=1))
BC_kmers[upper.tri(BC_kmers)] <- 0
BC_kmersDist <- as.dist(BC_kmers, diag = TRUE)
df <- data.frame( Kmers=BC_kmersDist[lower.tri(BC_kmersDist)], TAX=as.dist(BC_TAX)[lower.tri(as.dist(BC_TAX))])
set.seed(1234)
mantel(xdis =BC_TAX ,ydis = BC_kmersDist,method = "pearson",permutations = 9999)
```
### Figure 1D (k-PCoA)
```
iMDS  <- ordinate(physeq_MERGED_zone_FRACTION, "PCoA", distance=BC_kmersDist)
physeq_MERGED_zone_FRACTION_kmears<-physeq_MERGED_zone_FRACTION
sample_names(physeq_MERGED_zone_FRACTION_kmears) 
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
```
####  simka only superficial samples -> used for bioregions (Figure 2 A/B) and Vibrio surface dispersion through Oceans (fig 3B to G) 
```
simka -in input_sinka_derep_SRF.txt -out simka_SRF -kmer-size 31 -max-reads 0 -min-shannon-index 1.5
```
### Figure2 A/B (bioregions)
```
library(scales)

file_paths <- c(
  "SRF_0.22-3_mat_abundance_braycurtis.csv",
  "SRF_3-2_mat_abundance_braycurtis.csv",
  "SRF_20-180_mat_abundance_braycurtis.csv",
  "SRF_180_2000_mat_abundance_braycurtis.csv"
)


matrices <- lapply(file_paths, read.csv, row.names = 1, sep = ";")
pcoa_results <- lapply(matrices, function(matrix) {
  dist_matrix <- as.dist(matrix) # Assicurati che sia una matrice di distanza
  cmdscale(dist_matrix, eig = TRUE, k = 3, add = TRUE)  
})


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

coordinate <- read.csv2("cytoscape_samples_coordinates.csv", sep = ";")
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
  rescale_rgb <- function(rgb_df) {
    rgb_df$r <- rescale_to_255(rgb_df$r)
    rgb_df$g <- rescale_to_255(rgb_df$g)
    rgb_df$b <- rescale_to_255(rgb_df$b)
    return(rgb_df)
  }
  
  rescaled_rgb_df <- rescale_rgb(rgb_df)
  rescaled_data <- merge(coordinate, rescaled_rgb_df, by = "Sample")
  
  # Conversione dei valori RGB in colori esadecimali
  merged_data$color <- apply(merged_data[, c("r", "g", "b")], 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
  
 
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


```


### Figure 3B (Cumulative Correlations)
```
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
dataframes <- list(a, b, c, d)

for (i in 1:length(dataframes)) {
  colnames(dataframes[[i]])[1] <- "years.trav.time"
  colnames(dataframes[[i]])[2] <- "R"
  colnames(dataframes[[i]])[3] <- "Fraction"
}


a <- dataframes[[1]]
b <- dataframes[[2]]
c <- dataframes[[3]]
d <- dataframes[[4]]
allDF <- rbind(a, b,c,d)
colnames(allDF)[1]<-"years.trav.time"
newSTorder = c( "0.22-3","5-20","20-180" ,"180-2000")

allDF$Fraction<- as.character(allDF$Fraction)
allDF$Fraction <- factor(allDF$Fraction, levels=newSTorder) 
allDF$years.trav.time_ranges<-round(allDF$years.trav.time+0.5)
allDF$years.trav.time_ranges2<-round_any(allDF$years.trav.time, 0.5, f=ceiling)
 

#### fig 3B
ggplot(allDF, aes(x = years.trav.time, y = R, color = Fraction)) +
    geom_point(alpha = 0.1) +#, color = "black"
    geom_smooth(method = "gam", se = FALSE, formula = y ~ s(x, bs = "cs")) +  
    scale_x_continuous(limit = c(0, 5), breaks = round(seq(min(allDF$years.trav.time),
                                                               max(allDF$years.trav.time),
                                                               by = 0.5), 1)) +
    labs(color = "")+
  geom_smooth(method = "gam", se = FALSE, formula = y ~ s(x, bs = "cs"), color = "darkred", linetype = "dashed", size = 1.5) +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank())  

```
### Figure 3C (Correlation similarity vs travel time)
```

library(reshape2)
library(geodist)
library(ggplot2)
library(ggpubr)
library(diffcor)


coordinate<-read.csv2("cytoscape_samples_coordinates.csv",sep = ";")
color_palette <- c("SRF_0.22-3" = "red", "SRF_5-20" = "blue", "SRF_20-180" = "green", "SRF_180-2000" = "purple")
fraction_levels <- c("SRF_0.22-3", "SRF_5-20", "SRF_20-180", "SRF_180-2000")

list2<-data.frame(coordinate$Longitude,coordinate$Latitude,coordinate$Sample)
colnames(list2)[1] ="longitude"
colnames(list2)[2] ="latitude"
colnames(list2)[3] ="name"

distance_matrix <- geodist(list2, measure = 'geodesic' )/1000  
colnames(distance_matrix) <- list2$name
rownames(distance_matrix) <- list2$name

mat=as.matrix(distance_matrix)
melted<-reshape2::melt(mat)
newDF<-melted
newDF$merged<-paste(newDF$Var1,newDF$Var2,sep = "_")
colnames(newDF)[4] ="conx"
colnames(newDF)[3] ="dist_km"
file_paths <- c(
  "cytoscape/file_conx/SRF_0.22-3_TimTrav.csv",
  "cytoscape/file_conx/SRF_5-20_TimTrav.csv",
  "cytoscape/file_conx/SRF_20-180_TimTrav.csv",
  "cytoscape/file_conx/SRF_180-2000_TimTrav.csv"
)

correlation_data <- list()
for (file_path in file_paths) {
  cc <- read.csv(file_path, sep = ",")
  name <- gsub("cytoscape/file_conx/|_TimTrav.csv", "", file_path)
  cc$fraction <- rep(name, length(cc$conx))
  merged_data <- merge(cc, newDF, by = "conx", all.x = TRUE)
  merged_data <- merged_data[merged_data$tempo.medio.anni < 1.5,]
  correlation_data[[name]] <- merged_data
}

plot <- ggplot()
for (fraction in names(correlation_data)) {
  df <- correlation_data[[fraction]]
  df$fraction <- factor(reorder(df$fraction, -df$tempo.medio.anni), levels = fraction_levels)  # Definizione dei livelli della variabile "fraction" con ordine inverso della dist_km
  plot <- plot + geom_smooth(data = df, aes(y = similarity, x = tempo.medio.anni, color = fraction), method = "lm") +
    geom_smooth(data = df, aes(y = similarity, x = tempo.medio.anni), method = "lm", color = color_palette[fraction], fill = color_palette[fraction], alpha = 0.2)
}
plot <- plot + scale_color_manual(values = color_palette, guide = guide_legend(order = 1))
plot


```
####Fisher's z-Tests
```
df <- data.frame(fraction = character(), correlation = numeric(), p_value = numeric(), diff_corr = numeric(), stringsAsFactors = FALSE)
for (fraction in names(correlation_data)) {
  data <- correlation_data[[fraction]]
  r <- round(cor(data$similarity, data$tempo.medio.anni), 2)
  p <- cor.test(data$similarity, data$tempo.medio.anni)$p.value
  df <- rbind(df, data.frame(fraction = fraction, correlation = r, p_value = p, length=nrow(data)))
}

df
diff_df <- data.frame(fraction1 = character(), fraction2 = character(), diff_corr = numeric(), stringsAsFactors = FALSE)
for (i in 1:(nrow(df) - 1)) {
  for (j in (i + 1):nrow(df)) {
    fraction1 <- df$fraction[i]
    fraction2 <- df$fraction[j]
    correlation1 <- df$correlation[i]
    correlation2 <- df$correlation[j]
    length1 <- df$length[i]
    length2 <- df$length[j]
    diff_corr <- diffcor.two(correlation1, correlation2, length1, length2, digit = 3)
    diff_df <- rbind(diff_df, data.frame(fraction1 = fraction1, fraction2 = fraction2, diff_corr))
  }
}

```
### Figure 3D (Correlation similarity vs km)
```

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

coordinate<-read.csv2("cytoscape_samples_coordinates.csv",sep = ";")

rm(list2)
list2<-data.frame(coordinate$Longitude,coordinate$Latitude,coordinate$Sample)
colnames(list2)[1] ="longitude"
colnames(list2)[2] ="latitude"
colnames(list2)[3] ="name"
distance_matrix <- geodist(list2, measure = 'geodesic' )/1000 
colnames(distance_matrix) <- list2$name
rownames(distance_matrix) <- list2$name
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

color_palette <- c("SRF_0.22-3" = "red", "SRF_5-20" = "blue", "SRF_20-180" = "green", "SRF_180-2000" = "purple")
fraction_levels <- c("SRF_0.22-3", "SRF_5-20", "SRF_20-180", "SRF_180-2000")
correlation_data <- list()

for (file_path in file_paths) {
  cc <- read.csv(file_path, sep = ",")
  name <- gsub("cytoscape/file_conx/|_TimTrav.csv", "", file_path)
  cc$fraction <- rep(name, length(cc$conx))
  merged_data <- merge(cc, newDF, by = "conx", all.x = TRUE)
  correlation_data[[name]] <- merged_data
}

all_data <- bind_rows(correlation_data)
unique_data <- all_data %>% distinct(conx, .keep_all = TRUE)

unique_data<-unique_data[, !(names(unique_data) %in% c("X1", "X", "Var1.x" , "Var2.x","Var1.y"  ,"Var2.y"))]


write.csv(unique_data, "unique_conx_data.csv", row.names = FALSE)

```
#### Searoute in python
```

import pandas as pd
import searoute as sr

df = pd.read_csv('unique_conx_data.csv')
distances = []
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
df['searoutekm'] = distances
df.to_csv('unique_conx_data_with_distances.csv', index=False)
```
```
conx_w_coord<-read.csv("km_onlywaters/unique_conx_data_with_distances.csv")
merge_and_filter <- function(df, conx_w_coord) {
  df <- merge(df, conx_w_coord[, c("conx", "searoutekm")], by = "conx", all.x = TRUE)
  df <- df %>% filter(searoutekm <= 5000)
  return(df)
}


correlation_data_SR <- lapply(correlation_data, merge_and_filter, conx_w_coord)


plot <- ggplot()
for (fraction in names(correlation_data_SR)) {
  df <- correlation_data_SR[[fraction]]
  df$fraction <- factor(reorder(df$fraction, -df$searoutekm), levels = fraction_levels)
  plot <- plot + 
    geom_smooth(data = df, aes(y = similarity, x = searoutekm, color = fraction), method = "lm") 
}
plot + scale_color_manual(values = color_palette, guide = guide_legend(order = 1))
```

####Fisher's z-Tests
```
df <- data.frame(fraction = character(), correlation = numeric(), p_value = numeric(), diff_corr = numeric(), stringsAsFactors = FALSE)
for (fraction in names(correlation_data_SR)) {
  data <- correlation_data_SR[[fraction]]
  r <- round(cor(data$similarity, data$tempo.medio.anni), 2)
  p <- cor.test(data$similarity, data$tempo.medio.anni)$p.value
  df <- rbind(df, data.frame(fraction = fraction, correlation = r, p_value = p, length=nrow(data)))
}

df
diff_df <- data.frame(fraction1 = character(), fraction2 = character(), diff_corr = numeric(), stringsAsFactors = FALSE)
for (i in 1:(nrow(df) - 1)) {
  for (j in (i + 1):nrow(df)) {
    fraction1 <- df$fraction[i]
    fraction2 <- df$fraction[j]
    correlation1 <- df$correlation[i]
    correlation2 <- df$correlation[j]
    length1 <- df$length[i]
    length2 <- df$length[j]
    diff_corr <- diffcor.two(correlation1, correlation2, length1, length2, digit = 3)
    diff_df <- rbind(diff_df, data.frame(fraction1 = fraction1, fraction2 = fraction2, diff_corr))
  }
}

```
### Figure S4 (Correlation comparison between travel time vs km)
```
zscore <- function(x) {
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

correlation_data_combined <- data.frame()
for (fraction in names(correlation_data_SR)) {
  df <- correlation_data_SR[[fraction]]
  df$fraction <- rep(fraction, nrow(df))
  
  df_km <- df[df$searoutekm <= 5000,]
  df_km$log_dist_km <- log1p(df_km$searoutekm)
  df_km$z_log_dist_km <- zscore(df_km$log_dist_km)
  df_km_long <- melt(df_km, id.vars = c("conx", "similarity", "fraction"), measure.vars = "z_log_dist_km", variable.name = "type", value.name = "z_log_value")
  
  df_time <- df[df$tempo.medio.anni < 1.5,]
  df_time$log_tempo_anni <- log1p(df_time$tempo.medio.anni)
  df_time$z_log_tempo_anni <- zscore(df_time$log_tempo_anni)
  df_time_long <- melt(df_time, id.vars = c("conx", "similarity", "fraction"), measure.vars = "z_log_tempo_anni", variable.name = "type", value.name = "z_log_value")
  correlation_data_combined <- rbind(correlation_data_combined, df_km_long, df_time_long)
}


plot <- ggplot(correlation_data_combined, aes(x = z_log_value, y = similarity, color = fraction, linetype = type)) +
  geom_smooth(method = "lm") +
  stat_cor(aes(label = paste(..r.label.., sep = "")), method = "pearson", geom = "text", position = position_jitter(width = 0.2, height = 0), show.legend = FALSE) +
  scale_linetype_manual(values = c("z_log_dist_km" = "dashed", "z_log_tempo_anni" = "solid"), labels = c("z_log_dist_km" = "Distance (z_log_km)", "z_log_tempo_anni" = "Time (z_log_years)")) +
  scale_color_manual(values = color_palette, guide = guide_legend(order = 1))  +
  theme_minimal() +
  theme(
    axis.title.x = element_text(vjust = -0.5),
    legend.position = "top"
  )

```
### Figure 3E (Frequency of edges and their strength of similarity)
```
library(dplyr)
library(reshape2)
library(geosphere)
library(ggplot2)

process_file <- function(file_path, time_travel_file, coordinates_file, matrix_name) {
  # Caricare la matrice dal file CSV
  mat <- as.matrix(read.table(file_path, sep=";", header=TRUE, row.names=1))

  df <- melt(mat)
  new_df <- df[!df$Var1 == df$Var2,]
  new_df$similarity <- 1 - new_df$value
  new_df$conx <- paste(new_df$Var1, new_df$Var2, sep = "_")
  B <- read.csv2(time_travel_file, sep = " ", header = FALSE)
  colnames(B) <- c("Stazione.partenza", "Stazione.arrivo", "tempo.medio.giorni", "devStand.giorni")
  B <- B[as.numeric(B$tempo.medio.giorni) >= 0, ]
  B <- B[B$Stazione.partenza != B$Stazione.arrivo, ]
  B$tempo.medio.anni <- as.numeric(B$tempo.medio.giorni) / 365
  B$conx <- paste(B$Stazione.partenza, B$Stazione.arrivo, sep = "_")
  C <- merge(new_df, B, by = "conx")
  my_network <- C[C$tempo.medio.anni <= 1.5,]
  coordinates <- read.csv(coordinates_file, sep = ";")
  points_list <- coordinates[coordinates$Sample %in% my_network$Var1,]
  dist_matrix <- distm(points_list[, c("Longitude", "Latitude")], fun = distVincentyEllipsoid)
  min_dist <- apply(dist_matrix, 2, function(x) {
    dists <- sort(x)
    dists[2] # La seconda distanza più piccola è la distanza dal punto più vicino
  })
  edges_df <- data.frame(
    from = rep(points_list$Sample, each = nrow(points_list)),
    to = rep(points_list$Sample, times = nrow(points_list)),
    distance = min_dist
  )
  edges_df <- edges_df[edges_df$from != edges_df$to, ]
  df2 <- merge(edges_df, my_network, by.x = c("from", "to"), by.y = c("Var1", "Var2"))
  df1 <- df2[df2$similarity > 0,]
  n <- 4
  breaks <- seq(0, 1, by = 1/n)
  breaks1 <- seq(0, 1.5, by = 0.5)
  df1$simi_ranges <- cut(as.numeric(df1$similarity), breaks = breaks, labels = FALSE)
  df1$anni_ranges <- cut(as.numeric(df1$tempo.medio.anni), breaks = breaks1, labels = FALSE)
  df1$matrix_name <- matrix_name
  df3<-table(df1$simi_ranges)
  df3$name<-gsub("_mat_abundance_braycurtis.csv", "", file_path)
  return(df3)
}
files <- c("SRF_0.22-3_mat_abundance_braycurtis.csv", 
           "SRF_3-2_mat_abundance_braycurtis.csv", 
           "SRF_20-180_mat_abundance_braycurtis.csv", 
           "SRF_180_2000_mat_abundance_braycurtis.csv")
matrix_names <- c("SRF_0.22-3", "SRF_3-2", "SRF_20-180", "SRF_180-2000")
time_travel_file <- "filetime_2023-09-29_17-02-30_tiempo_stdtiempo.txt"
coordinates_file <- "cytoscape_samples_coordinates.csv"
final_results <- do.call(rbind, lapply(1:length(files), function(i) {
  process_file(files[i], time_travel_file, coordinates_file, matrix_names[i])
}))
print(final_results)

final_df <- as.data.frame(final_results)
colnames(final_df) <- c("V1", "V2", "V3", "V4", "Size")
final_df$V1 <- as.numeric(final_df$V1)
final_df$V2 <- as.numeric(final_df$V2)
final_df$V3 <- as.numeric(final_df$V3)
final_df$V4 <- as.numeric(final_df$V4)

final_df_long <- melt(final_df, id.vars = "Size", variable.name = "Ranges", value.name = "Numero")
newSTorder = c("SRF_0.22-3", "SRF_3-2", "SRF_20-180", "SRF_180_2000")
final_df_long$Size <- factor(final_df_long$Size, levels = newSTorder)

newSTorder1 <- c("V4", "V3", "V2", "V1")
final_df_long$Ranges <- factor(final_df_long$Ranges, levels = newSTorder1)

ggplot(final_df_long, aes(x = Size, y = Numero, fill = Ranges)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("V1" = "green", "V2" = "yellow", "V3" = "darkorange", "V4" = "darkred")) +
  theme_bw() +
  labs(x = "Size", y = "Numero", fill = "Ranges")


```
### Figure 3F (Node centrality)
```
library(igraph)
library(dplyr)
library(reshape2)
library(ggpubr)
calcola_node_centrality <- function(file_list, file_time, file_coordinates) {
  risultati_lista <- list()
  
  for (file_abundance in file_list) {
    mat <- as.matrix(read.table(file_abundance, sep = ";", header = TRUE, row.names = 1))
    df <- melt(mat)
    new_df <- df[!df$Var1 == df$Var2, ]
    new_df$similarity <- 1 - new_df$value
    new_df$conx <- paste(new_df$Var1, new_df$Var2, sep = "_")
    B <- read.csv2(file_time)
    C <- merge(new_df, B, "conx")
    my_network <- C[C$tempo.medio.anni < 1.5, ]
    coordinates <- read.csv(file_coordinates, sep = ";")
    points_list <- coordinates[coordinates$Sample %in% my_network$Var1, ]
    dist_matrix <- distm(points_list[, c("Longitude", "Latitude")], fun = distVincentyEllipsoid)
    min_dist <- apply(dist_matrix, 2, function(x) {
      dists <- sort(x)
      dists[2]
    })
    edges_df <- data.frame(
      from = rep(points_list$Sample, each = nrow(points_list)),
      to = rep(points_list$Sample, times = nrow(points_list)),
      distance = min_dist
    )
    edges_df <- edges_df[edges_df$from != edges_df$to, ]
    edges_df <- merge(edges_df, points_list, by.x = "from", by.y = "Sample")
    edges_df <- merge(edges_df, points_list, by.x = "to", by.y = "Sample", suffixes = c("", "_to"))
    
    df2 <- merge(x = edges_df, y = my_network, by.x = c("from", "to"), by.y = c("Var1", "Var2"))
    df1 <- filter(df2, similarity > 0)
    edges <- df1[, c("from", "to", "similarity")]
    g <- graph.data.frame(edges, directed = TRUE)
    centr <- centr_degree(g)$res
    risultati <- data.frame(
      nodo = V(g)$name,
      node_centrality = centr
    )
    risultati <- risultati[order(-centr), ]
    risultati_lista[[file_abundance]] <- risultati
  }
  
  return(risultati_lista)
}
file_list <- c(
  "SRF_0.22-3_mat_abundance_braycurtis.csv",
  "SRF_180_2000_mat_abundance_braycurtis.csv",
  "SRF_20-180_mat_abundance_braycurtis.csv",
  "SRF_3-2_mat_abundance_braycurtis.csv"
)

file_time <- read.csv2("filetime_2023-09-29_17-02-30_tiempo_stdtiempo.txt",sep = " ", header = F)
colnames(file_time)<-c("Stazione.partenza"  ,"Stazione.arrivo","tempo.medio.giorni" ,"devStand.giorni")
file_time <- file_time[as.numeric(B$tempo.medio.giorni) >= 0, ]
file_time <- file_time[file_time$Stazione.partenza !=file_time$Stazione.arrivo, ]
file_time$tempo.medio.anni <-  as.numeric(file_time$tempo.medio.giorni) / 365
file_time$conx <- paste(file_time$Stazione.partenza, file_time$Stazione.arrivo, sep = "_")
file_coordinates <- "cytoscape_samples_coordinates.csv"

risultati_centralita <- calcola_node_centrality(file_list, file_time, file_coordinates)

risultati_file_1 <- risultati_centralita[["SRF_0.22-3_mat_abundance_braycurtis.csv"]]
risultati_file_2 <- risultati_centralita[["SRF_180_2000_mat_abundance_braycurtis.csv"]]
risultati_file_3 <- risultati_centralita[["SRF_20-180_mat_abundance_braycurtis.csv"]]
risultati_file_4 <- risultati_centralita[["SRF_3-2_mat_abundance_braycurtis.csv"]]

df_completo <- data.frame()
for (nome_file in names(risultati_centralita)) {
  risultati <- risultati_centralita[[nome_file]]
  risultati$nome_file <- gsub("_mat_abundance_braycurtis.csv", "", nome_file)
  df_completo <- bind_rows(df_completo, risultati)
}
print(df_completo)
ggplot(df_completo, aes(x=nome_file, y=node_centrality, fill=nome_file)) + 
  geom_boxplot(alpha=0.3,outlier.shape = NA) +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons)+   stat_compare_means(label.y = 70)  
df_completo$zona<-sub("_.*", "", df_completo$nodo)

```
### Figure  S5A (Distribution of node centrality boxplot)
```


ggplot(data = df_completo, aes(x = zona, y = node_centrality)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(color = nome_file), position = position_jitter(0.2), alpha = 0.7) +  
  labs(x = "Zona", y = "node_centrality") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_hline(yintercept = media_node_centrality, linetype = "dashed", color = "red")


```
### Figure  S5B (Top 30 stations node centrality)
```

  top_stazioni <- df_completo %>% 
     arrange(desc(node_centrality))
   
   totali_stazioni <- top_stazioni %>%
     group_by(nodo) %>%
     summarise(total_node_centrality = sum(node_centrality))
   top_stazioni <- top_stazioni %>%
     mutate(nodo = factor(nodo, levels = rev(totali_stazioni$nodo[order(-totali_stazioni$total_node_centrality)])))
   top_30_nodi <- top_stazioni %>%
     arrange(desc(node_centrality)) %>%  
     distinct(nodo, .keep_all = TRUE) %>% 
     slice(1:30) %>%  
     pull(nodo)  
   top_stazioni_filtrato <- top_stazioni %>%
     filter(nodo %in% top_30_nodi)
   print(top_stazioni_filtrato)
   ggplot(data = top_stazioni_filtrato, aes(x = nodo, y = node_centrality, fill = nome_file)) +
     geom_bar(stat = "identity") +
     coord_flip() +
     labs(x = "Stazione", y = "node_centrality", title = "Top 150 stazioni più rilevanti") +
     theme_minimal()
   
```
### Figure  S6A (Heatmap per mean travel time for 20-180)
```
 library(ggplot2)
 library(maps)
 library(reshape2)
 library(dplyr)
 library(geosphere)
 
 mat <- as.matrix(read.table("SRF_20-180_mat_abundance_braycurtis.csv", sep=";", header=TRUE, row.names=1))
 df <- melt(mat)
 new_df <- df[df$Var1 != df$Var2,]
 new_df$similarity <- 1 - new_df$value
 new_df$conx <- paste(new_df$Var1, new_df$Var2, sep="_")
 B <- read.csv2("filetime_2023-09-29_17-02-30_tiempo_stdtiempo.txt",sep = " ", header = F)
 colnames(B)<-c("Stazione.partenza"  ,"Stazione.arrivo","tempo.medio.giorni" ,"devStand.giorni")
 B <- B[as.numeric(B$tempo.medio.giorni) >= 0, ]
 B <- B[B$Stazione.partenza != B$Stazione.arrivo, ]
 B$tempo.medio.anni <-  as.numeric(B$tempo.medio.giorni) / 365
 B$conx <- paste(B$Stazione.partenza, B$Stazione.arrivo, sep = "_")
 C <- merge(new_df, B, by="conx")
 my_network <- C[C$tempo.medio.anni <= 1.5,]
 coordinates <- read.csv("cytoscape_samples_coordinates.csv", sep=";")
 points_list <- coordinates[coordinates$Sample %in% my_network$Var1,]
 dist_matrix <- distm(points_list[, c("Longitude", "Latitude")], fun = distVincentyEllipsoid)
 min_dist <- apply(dist_matrix, 2, function(x) sort(x)[2])
 edges_df <- data.frame(
   from = rep(points_list$Sample, each = nrow(points_list)),
   to = rep(points_list$Sample, times = nrow(points_list)),
   distance = min_dist
 )
 edges_df <- edges_df[edges_df$from != edges_df$to,]
 edges_df <- merge(edges_df, points_list, by.x="from", by.y="Sample")
 edges_df <- merge(edges_df, points_list, by.x="to", by.y="Sample", suffixes=c("", "_to"))
 df2 <- merge(edges_df, my_network, by.x=c("from", "to"), by.y=c("Var1", "Var2"))
 df1 <- filter(df2, similarity > 0)
 classify_ocean_region <- function(station) {
   if (grepl("^ANE", station)) {
     return("North Atlantic East")
   } else if (grepl("^ANW", station)) {
     return("North Atlantic West")
   } else if (grepl("^ASE", station)) {
     return("South Atlantic East")
   } else if (grepl("^ASW", station)) {
     return("South Atlantic West")
   } else if (grepl("^PON", station)) {
     return("North Pacific")
   } else if (grepl("^PSW", station)) {
     return("South Pacific West")
   } else if (grepl("^PSE", station)) {
     return("South Pacific East")
   } else if (grepl("^ION", station)) {
     return("North Indian")
   } else if (grepl("^IOS", station)) {
     return("South Indian")
   } else if (grepl("^SOC", station)) {
     return("Southern Ocean")
   } else if (grepl("^ARC", station)) {
     return("Arctic")
   } else if (grepl("^MED", station)) {
     return("Mediterranean")
   } else if (grepl("^RED", station)) {
     return("Red Sea")
   } else {
     return(NA)
   }
 }
 df1$Ocean_reg <- sapply(df1$from, classify_ocean_region)
 df1$Ocean_reg_to <- sapply(df1$to, classify_ocean_region)
 
 combined_links <- df1 %>%
    group_by(Ocean_reg, Ocean_reg_to) %>%
   summarise(
     count = n(),
     mean_travel_time = mean(as.numeric(tempo.medio.anni), na.rm = TRUE),
     mean_similarity = mean(similarity, na.rm = TRUE),
     .groups = 'drop'
   )
 ggplot(combined_links, aes(x = Ocean_reg, y = Ocean_reg_to, fill = mean_travel_time)) +
   geom_tile() +
   geom_text(aes(label = round(mean_travel_time, 2)), color = "white", size = 3) +
   scale_fill_gradient(low = "red", high = "blue") +
   labs(title = "Mean Travel Time Between and Within Ocean Regions 20-180",
        x = "Start Ocean Region",
        y = "Arrival Ocean Region",
        fill = "Mean Travel Time") +
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
```
### Figure  S6B (Heatmap per mean similarity for 20-180)
```
 ggplot(combined_links, aes(x = Ocean_reg, y = Ocean_reg_to, fill = mean_similarity)) +
   geom_tile() +
   geom_text(aes(label = round(mean_similarity, 2)), color = "white", size = 3) +
   scale_fill_gradientn(colors = c("blue",  "red")) +
   labs(title = "Mean Similarity Between and Within Ocean Regions 20-180",
        x = "Start Ocean Region",
        y = "Arrival Ocean Region",
        fill = "Mean Similarity") +
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

```
### Contigs Analysis
10. **Megahit Co-Assembling using Vibrio reads**
```
cat sets.txt
```
> ARC
> ANE
> ANW
> ASE
> ASW
> ION
> IOS
> MED
> PON
> PSE
> PSW
> RED
> SOC
```
for SET in `cat sets.txt`
do
R1s=`ls *_1.fa | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
R2s=`ls *_2.fa | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
 megahit -1 $R1s -2 $R2s -o $SET-co-assembly.fa 
done

```
11. **CAT Species Taxonomy**
```
for i in *.fa
do
CAT contigs -c $i -d CAT_prepare_20210107/2021-01-07_CAT_database -tCAT_prepare_20210107/2021-01-07_taxonomy -o BAT_output/${i}_BAT
CAT add_names -i BAT_output/${i}_BAT.contig2classification.txt -o BAT_off/${i}_.names_off.txt -t CAT_prepare_20210107/2021-01-07_taxonomy --only_official
CAT summarise -c $i -i BAT_off/${i}_.names_off.txt -o BAT_class/${i}_BAT_summ.txt
done
```
12. **Salmon Quantification**
```
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

```
### Figure  4A (Distribution of *Vibrio* species between FLV and PAV)
```
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

top20phyla = names(sort(genus.sum, TRUE))[1:20]
physeq20 = prune_taxa((tax_table(physeq)[, "specie"] %in% top20phyla), physeq)
physeq_norm<-normalize(physeq20, method = "RLE")
physeq <- physeq_norm   
if (!physeq@otu_table@taxa_are_rows) {
  otu_tab <- t(physeq@otu_table)
} else {
  otu_tab <- physeq@otu_table
}

otu_tab <- as.data.frame(as(otu_tab, "matrix"))
tax_table <- as.data.frame(as(physeq@tax_table, "matrix"))
tax_table$Taxon <- rownames(tax_table)
otu_taxa <- merge(otu_tab, tax_table, by = "row.names")
rownames(otu_taxa) <- otu_taxa$Row.names
otu_taxa$Row.names <- NULL
taxa <- "species"  
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
  mutate(scaled_Abundance = Abundance / Total)
taxa_order <- taxa_counts %>%
  filter(Compartment == "PROT") %>%
  arrange(desc(scaled_Abundance)) %>%
  distinct(!!sym(taxa)) %>%
  pull(!!sym(taxa))
taxa_counts <- taxa_counts %>%
  mutate(!!sym(taxa) := factor(!!sym(taxa), levels = unique(taxa_order)))

ggplot(taxa_counts, aes(x = !!sym(taxa), y = scaled_Abundance, fill = Compartment)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "right")

```
### Figure  4A and B (Distribution of *Vibrio* species among oceanic regions)
```

otu_tab <- as.data.frame(as(otu_tab, "matrix"))
tax_table <- as.data.frame(as(physeq@tax_table, "matrix"))
tax_table$Taxon <- rownames(tax_table)
otu_taxa <- merge(otu_tab, tax_table, by = "row.names")
rownames(otu_taxa) <- otu_taxa$Row.names
otu_taxa$Row.names <- NULL
taxa <- "species"  # replace with the desired taxonomic level, if different

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

taxa_counts <- otu_taxa %>%
  group_by(across(all_of(taxa))) %>%
  summarise(across(starts_with("BACT_") | starts_with("PROT_"), sum)) %>%
  pivot_longer(cols = starts_with("BACT_") | starts_with("PROT_"),
               names_to = "Sample", values_to = "Abundance") %>%
  mutate(Compartment = ifelse(grepl("^BACT_", Sample), "BACT", "PROT")) %>%
  left_join(metaDF, by = c("Sample" = "Campione"))

taxa_counts <- taxa_counts %>%
  group_by(!!sym(taxa)) %>%
  mutate(Total = sum(Abundance)) %>%
  ungroup() %>%
  mutate(scaled_Abundance = Abundance / Total)

ggplot(taxa_counts, aes(x = !!sym(taxa), y = scaled_Abundance, fill = Oceano   )) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "right")
