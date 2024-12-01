## script for ananlysing the gsSAR202 genome

# download SAR202 genome from: https://doi.org/10.6084/m9.figshare.26548255
############################################################################################
###############################SAR202 genome analysis##########################################
############################################################################################
#!/bin/bash

# Create a script for dRep genome deduplication, CheckM analysis, GTDB-Tk annotation, phylogenetic tree construction, gene prediction with Prokka, and functional annotation with KEGG and eggNOG

# Specify input and output directories
INPUT_DIR="path/to/your/genomes"  # Replace "path/to/your/genomes" with the directory containing your genome files
OUTPUT_DIR="path/to/output_directory"  # Replace "path/to/output_directory" with the output directory
CHECKM_OUTPUT_DIR="path/to/checkm_output_directory"  # Replace "path/to/checkm_output_directory" with the CheckM output directory
GTDBTK_OUTPUT_DIR="path/to/gtdbtk_output_directory"  # Replace "path/to/gtdbtk_output_directory" with the GTDB-Tk output directory
PHYLO_OUTPUT_DIR="path/to/phylogenetic_tree_directory"  # Replace "path/to/phylogenetic_tree_directory" with the phylogenetic tree output directory
PROKKA_OUTPUT_DIR="path/to/prokka_output_directory"  # Replace "path/to/prokka_output_directory" with the Prokka output directory
ANNOTATION_OUTPUT_DIR="path/to/annotation_output_directory"  # Replace "path/to/annotation_output_directory" with the functional annotation output directory

# Run dRep for deduplication
# -g specifies the input genome files, --S_algorithm selects the clustering algorithm (default here is fastANI),
# -pa 0.9 sets the primary ANI clustering threshold to 90%,
# -sa 0.95 sets the secondary ANI threshold to 95%.
dRep dereplicate $OUTPUT_DIR -g $INPUT_DIR/*.fasta --S_algorithm fastANI -pa 0.9 -sa 0.95

# Run CheckM to analyze the quality of the deduplicated genomes
# The analysis will focus on genomes with completeness greater than 50%
checkm lineage_wf $OUTPUT_DIR/dereplicated_genomes $CHECKM_OUTPUT_DIR -x fasta --tab_table
awk -F "\t" '$12 > 50' $CHECKM_OUTPUT_DIR/storage/bin_stats_ext.tsv > $CHECKM_OUTPUT_DIR/high_quality_genomes.tsv

# Extract genome file paths with completeness > 50%
awk -F "\t" '{print $1}' $CHECKM_OUTPUT_DIR/high_quality_genomes.tsv > $CHECKM_OUTPUT_DIR/high_quality_genomes_list.txt

# Run GTDB-Tk to obtain homologous genes from high-quality genomes
GTDBTK_INPUT_DIR="path/to/gtdbtk_input_directory"  # Replace with the directory where GTDB-Tk will use the filtered genomes
mkdir -p $GTDBTK_INPUT_DIR
while read genome; do
  cp $OUTPUT_DIR/dereplicated_genomes/$genome $GTDBTK_INPUT_DIR
done < $CHECKM_OUTPUT_DIR/high_quality_genomes_list.txt
gtdbtk classify_wf --genome_dir $GTDBTK_INPUT_DIR --out_dir $GTDBTK_OUTPUT_DIR --cpus 4

# Construct a phylogenetic tree using IQ-TREE
# Run IQ-TREE on aligned homologous sequences to build a phylogenetic tree
mkdir -p $PHYLO_OUTPUT_DIR
ALIGNMENT_FILE="$GTDBTK_OUTPUT_DIR/gtdbtk.bac120.msa.fasta"  # Replace if the file name differs
cd $PHYLO_OUTPUT_DIR
iqtree -s $ALIGNMENT_FILE -m MFP -bb 1000 -nt AUTO

# Run Prokka to predict genes from high-quality genomes
mkdir -p $PROKKA_OUTPUT_DIR
while read genome; do
  BASENAME=$(basename $genome .fasta)
  prokka --outdir $PROKKA_OUTPUT_DIR/$BASENAME --prefix $BASENAME $OUTPUT_DIR/dereplicated_genomes/$genome
done < $CHECKM_OUTPUT_DIR/high_quality_genomes_list.txt

# Perform KEGG and eggNOG annotation on predicted genes
mkdir -p $ANNOTATION_OUTPUT_DIR
while read genome; do
  BASENAME=$(basename $genome .fasta)
  PROKKA_GFF_FILE="$PROKKA_OUTPUT_DIR/$BASENAME/$BASENAME.gff"
  
  # KEGG annotation (using a hypothetical tool or script, replace with the actual command)
  kegg_annotation_tool -i $PROKKA_GFF_FILE -o $ANNOTATION_OUTPUT_DIR/$BASENAME"_kegg_annotation.txt"

  # eggNOG annotation (using eggNOG-mapper)
  emapper.py -i $PROKKA_GFF_FILE -o $ANNOTATION_OUTPUT_DIR/$BASENAME"_eggnog_annotation" --cpu 4

done < $CHECKM_OUTPUT_DIR/high_quality_genomes_list.txt

# Output information
echo "dRep deduplication completed, output directory: $OUTPUT_DIR"
echo "CheckM analysis completed, genomes with completeness > 50% are listed in: $CHECKM_OUTPUT_DIR/high_quality_genomes.tsv"
echo "GTDB-Tk annotation completed, output directory: $GTDBTK_OUTPUT_DIR"
echo "Phylogenetic tree construction completed, output directory: $PHYLO_OUTPUT_DIR"
echo "Prokka gene prediction completed, output directory: $PROKKA_OUTPUT_DIR"
echo "KEGG and eggNOG annotation completed, output directory: $ANNOTATION_OUTPUT_DIR"

######################################################################################
## Figures plot using R language

## setting up the working direction 

setwd("./data")

## adding the R packages

library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggalluvial)
library(RColorBrewer)
library(ggrepel)
library(ggthemes)
library(dplyr)
library(ggbreak)
library(patchwork)


## Plot the Figure 1A

data = read.table("Figure.1A.txt", header = T, sep = '\t')
ggplot(data, aes(EGS, GC content)) + geom_point(size = 1.5) + geom_smooth(method = 'loess', span = 0.5, se = TRUE, level = 0.95)+
  stat_cor(label.y =0.7, 
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+ 
  theme_bw()


## Plot the Figure 1B

data = read.table("Figure.1B.txt", header = T, sep = '\t')
ggplot(data, aes(value, Group, fill = variable)) + 
  geom_bar(position = 'stack', width = 0.7,stat = 'identity', color = 'white' )+
  scale_fill_manual(values=c( "#fa7f6f", "#82b0d2","#beb8dc", "#fbe7b3" ))+
  theme_bw()


## Plot the Figure 2


size = read.csv('Figure.2.EGS', header = T, sep = ',')
head(size)

p = ggplot(size, aes(abundance, depth, color = size)) + geom_point(aes(color = size), shape = 21)+scale_color_gradient(low = "#c6fedc",high = "#f7797d")+
    theme_classic()

pp = p+scale_y_break(c(1000,2000),scales = 1.5,space = 0.3)+
  scale_y_break(c(2500,3000),scales = 1.5,space = 0.3)


### plot GC%


fig.gc = read.csv('Figure2.GC.txt', header = T, sep = ',')

## plot
p1 = ggplot(fig.gc, aes(value, depth, color = GC)) + geom_point(aes(color = GC), shape = 21)+scale_color_gradient(low = "#108dc6",high = "#d68237")+
    theme_classic()

pp1 = p1+scale_y_break(c(1000,2000),scales = 1.5,space = 0.3)+
  scale_y_break(c(2500,3000),scales = 1.5,space = 0.3)

### combine fig

pp + pp1



## Plor Figure 3

fig = read.table('Figure.3A', header = T, sep = '\t')

ggplot(fig, aes(abundance, Deep)) + geom_boxplot(aes(color = group)) + 
  geom_point(aes(color = group), position = "jitter")+
  scale_color_manual(values = c("#1D91C0", "#41B6C4",  "#C7E9B4", "#FEE08B"))+ 
  theme_bw()


## Figure 3B

fig = read.csv('Figure.3B.txt', header = T, row.names = 1, sep = ',')
head(fig)

data = melt(fig, id = 'id')
head(data)

ggplot(data, aes(variable, value, fill = id)) + 
  geom_bar(position = 'stack', width = 0.7,stat = 'identity', color = 'white' )+
  scale_fill_manual(values=c( "#dbdbdb", "#a67d83", "#e7dad2","#f7e1ed",
   "#9acdde", "#82b0d2" ))+
  theme_classic()


# NETWORK FIGURE

library(igraph)
library(psych)
library(ggvegan)
## SRF water layer
sar202 = read.table("clipboard",head=T,row.names=1, sep = '\t')
head(sar202)

env = read.table("clipboard",head=T,row.names=1, sep = '\t')
head(env)
env.data.log = log1p(env)
fix(env.data.log)
env <- na.omit(env.data.log)


decorana(sar202)

ii = rda(sar202 ~ 1, env)

ii1 = rda(sar202 ~ ., env)

#choose environmental factors

vif.cca(ii1)

# < 10
ii1 <- rda(sar202 ~ env5 + env6  + env7  + env8 + env9 + env11 + env12 + env13 + env14+
             env16 + env18 + env19  + env24 + env26 + env29 + env31 + env33 + env34+env36+env38+
             env41+env47, env.data.log)

vif.cca(ii1)

ii1 <- rda(sar202 ~ env5  + env7  +  env13 + env18 + env19  + env24 +  env29 + env31 + env33 +env36+env38+
             env41+env47, env)

vif.cca(ii1)


#test again
ii1 <- rda(sar202 ~ env5  + env7  +  env13 + env18 + env19  + env24 +  env29 + env31 + env33 +env36+env38+
             env41+env47, env)

#all factoers < 10
vif.cca(ii1)

## build network
occor = corr.test(sar202, env, method = "spearman", adjust = "bonferroni")
occor.p = occor$p # p value
occor.r = occor$r
occor.r[occor.p>0.05|occor.r> -0.3] = 0

# write.csv(occor.r, 'SRF.occor.r.csv')
occor.r = read.csv('./SRF.occor.r.csv', header = T, row.names = 1, sep = ',')

igraph = graph_from_biadjacency_matrix(occor.r, mode = 'total',  weighted=TRUE)

igraph

set.seed(123)
plot(igraph,layout= layout_in_circle  ,vertex.frame.color=NA,vertex.label.cex=0.3,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))



























