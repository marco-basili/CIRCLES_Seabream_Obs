#######################
#  CIRCLES 16S all samples Observation
######################
 library(dada2)

#####
# CUTADAPT on the row seq

for i in *_R1_001.fastq.gz; do    
    SAMPLE=$(echo "$i" | sed 's/_R1_001.fastq.gz//');      
      cutadapt         -g "CCTACGGGNGGCNGCAG;max_error_rate=0.15"         
      -G "GACTACNNGGGTATCTAATCC;max_error_rate=0.15"         
      --discard-untrimmed         
      -o ${SAMPLE}_R1_001_noP.fastq.gz         
      -p ${SAMPLE}_R2_001_noP.fastq.gz         
      ${SAMPLE}_R1_001.fastq.gz         
      ${SAMPLE}_R2_001.fastq.gz;          
done
#####


path <- "CIRCLES_obs_seq/Obs/"
list(path)
fnFs <- sort(list.files(path, pattern="R1_001_noP.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_001_noP.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_L001_R"), `[`, 1)

plotQualityProfile(fnRs[1:2])

filtFs <- file.path(path, "filtered_2026", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered_2026", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out_CIRCLES <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,200),
                                maxN=0, maxEE=c(3,5), truncQ=2, rm.phix=TRUE,
                                compress=TRUE, multithread=FALSE, verbose = TRUE)


errF <- learnErrors(filtFs, multithread=TRUE, verbose =TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, verbose =TRUE)   

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, verbose =TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, verbose =TRUE)

mergers_IT <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab_CIRCLES <- makeSequenceTable(mergers_IT)
seqtab_CIRCLES_nochim <- removeBimeraDenovo(seqtab_CIRCLES, method="consensus", multithread=TRUE)
sum(seqtab_CIRCLES_nochim)
dim(seqtab_CIRCLES_nochim)
table(nchar(getSequences(seqtab_CIRCLES_nochim)))

seqtab_CIRCLES_prune <- seqtab_CIRCLES_nochim[,nchar(
  colnames(seqtab_CIRCLES_nochim)) %in% seq(360, 440)]
sum(seqtab_CIRCLES_prune)/sum(seqtab_CIRCLES_nochim)
# 0.9738982


taxa_seqtab_CIRCLES <- assignTaxonomy(seqtab_CIRCLES_prune, "~/Desktop/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out_CIRCLES, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers_IT, getN), 
               rowSums(seqtab_CIRCLES_prune))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

tab_Obs<-data.frame(sample_data(ps_CIRCLE_meta))
tab_Obs$ID <- rownames(tab_Obs)
track <- as.data.frame(track)
track$ID <- rownames(track)

merged <- merge(tab_Obs, track, by = "ID")
merged$keep <- merged$nonchim/merged$input

merged %>%  ggplot(aes(x=Timepoint, y=keep))+
  geom_boxplot(aes(fill=Country)) +
  # geom_smooth(aes(shannon, colour = Matrix, fill = Matrix))+
  #facet_grid(~Country)+
  #stat_pvalue_manual(pwc, hide.ns = TRUE) +
  theme_bw()

#################


##################
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("vegan")
library("microbiome")
library("viridis") 
library(microeco)
library(microbiome)
library(ggrepel)
# install.packages("BiocManager")
# BiocManager::install("phyloseq")
library(phyloseq)
library(gridExtra)
library(tidyr)
library(dplyr)
library(ggh4x)

#library(phangorn)
#library(tidyverse)
#library(DECIPHER)
#library(msa)


##########

# 2026

##############
# phyloseq form DADA Output
###############
ps_CIRCLE_Obs <- phyloseq(otu_table(seqtab_CIRCLES_prune, taxa_are_rows=FALSE), 
                          tax_table(taxa_seqtab_CIRCLES))
ps_CIRCLE_Obs
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 46883 taxa and 911 samples ]
# tax_table()   Taxonomy Table:    [ 46883 taxa by 6 taxonomic ranks ]


write.csv(sample_names(ps_CIRCLE_Obs_prune), "All_sample_Name.csv")

##############
#  filter the phyloseq object
################

Metadata_CIRCLES<-read.csv("Metadata_All_Obs_pool.csv")

row.names(Metadata_CIRCLES) <- Metadata_CIRCLES$Sample_name
ps_Obs_meta <-merge_phyloseq(ps_CIRCLE_Obs_prune, sample_data(Metadata_CIRCLES))


ps_Obs_meta
# ps_Obs_meta
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 46883 taxa and 619 samples ]
# sample_data() Sample Data:       [ 619 samples by 16 sample variables ]
# tax_table()   Taxonomy Table:    [ 46883 taxa by 6 taxonomic ranks ]

# using Metadata_allObs_pool che contiene i campioni filtrati 
# con criteri scelti insieme:
'
- Rimossi Time point non coerenti
- Rimossi Fecal sample [dopo valutazione -> qui di seguito]
- Rimossi Fillet
- Rimossi Waste water
- Rimossi Inlet water
- Rimossi Sed e Wat Italia (escluso water gabbia)
'

# extract the sequence
dna <- Biostrings::DNAStringSet(taxa_names(ps_Obs_meta))
names(dna) <- taxa_names(ps_Obs_meta)
ps_Obs_meta <- merge_phyloseq(ps_Obs_meta, dna)
taxa_names(ps_Obs_meta) <- paste0("ASV", seq(ntaxa(ps_Obs_meta)))

ps_Obs_meta
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 46883 taxa and 619 samples ]
# sample_data() Sample Data:       [ 619 samples by 16 sample variables ]
# tax_table()   Taxonomy Table:    [ 46883 taxa by 6 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 46883 reference sequences ]


################
# focus only on GUT for decide if remove the Feces
ps_Gut_all<-prune_samples(sample_data(ps_Obs_meta)$Sample_type == "Gut" ,  ps_Obs_meta)
ps_Gut_all <- prune_taxa(taxa_sums(ps_Gut_all) > 0, ps_Gut_all)
ps_Gut_all <- prune_samples(sample_sums(ps_Gut_all) > 500, ps_Gut_all)

ps_Gut_all <- transform_sample_counts(ps_Gut_all, function(x) (x / sum(x))*100)

dist_gut <- phyloseq::distance(ps_Gut_all, method = "bray")
NMDS_gut <- ordinate(ps_Gut_all, dist_gut, method = "NMDS", trymax=100)

plot_ordination(
  ps_Gut_all , # for PcOA
    #subset_samples(ps_Obs_rel, sample_names(ps_Obs_rel) != "T5_Adults_10_Skin_Spain_Farm"),
  NMDS_gut, title="Bray_NMDS") + 
   geom_point(aes(shape=Country, color =sample_type_add), alpha = 1, size = 5)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))



# Barplot Class on the ps pooled for timepoints
ps_Gut_all_cl <- tax_glom( ps_Gut_all, taxrank = "Class",NArm=FALSE)
ps_to_use <- ps_Gut_all_cl
# Identify top 10
top10_taxa <- names(
  sort(taxa_sums(ps_to_use), decreasing = TRUE)[1:10])
physeq_top10 <- prune_taxa(top10_taxa, ps_to_use)

plot_bar(physeq_top10, 
         #x = "sample_type_add", 
         fill = "Class")+
  facet_nested(~Matrix+Country+Timepoint,  scales = "free_x",
               space = "free_x") +
  scale_fill_brewer(palette = "Set3")  +
  theme(  axis.text.x = element_blank(),
    legend.position = "right",
    panel.spacing.x = unit(0.1, "lines")  # reduce horizontal spacing
  )

#################

#  Total reads
sum(sample_sums(ps_Obs_meta))
#   [1] 17946070

prune_samples(sample_data( ps_Obs_meta)$WP != "Wild" ,  ps_Obs_meta)
# Filter samples with less than 500 reads
ps_CIRCLE_meta <- prune_samples(sample_sums(ps_Obs_meta) > 500, ps_Obs_meta)
ps_CIRCLE_meta
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 46883 taxa and 467 samples ]
# sample_data() Sample Data:       [ 467 samples by 16 sample variables ]
# tax_table()   Taxonomy Table:    [ 46883 taxa by 6 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 46883 reference sequences ]
write.csv(sample_names(ps_CIRCLE_meta), "CIRCLE_used_sample.csv")

#############################
# pool samples
###############################

# Use dplyr to get the first occurrence of the character column for each group
detach("package:tidyr", unload = TRUE)
library(dplyr)

first_occurrence_df_Obs <- sample_data(ps_CIRCLE_meta) %>%
  group_by(Pool_code) %>%
  dplyr::slice(1) %>%
  select(Pool_code, c(Timepoint, Stage, Sample_type, Type_env, Country, WP,
                      food_chain, sample_site, sample_type_add, sample_type_number, phase, Matrix)) %>%
  ungroup()

ps_pooled <- merge_samples(ps_CIRCLE_meta, "Pool_code")

ps_pooled
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 46883 taxa and 325 samples ]
# sample_data() Sample Data:       [ 325 samples by 16 sample variables ]
# tax_table()   Taxonomy Table:    [ 46883 taxa by 6 taxonomic ranks ]

# Replace the sample_data in the phyloseq object
new_sample_data <- sample_data(as.data.frame(first_occurrence_df_Obs))
sample_names(new_sample_data) <- new_sample_data$Pool_code
sample_data(ps_pooled)<-new_sample_data

# check N samples fish and env (ci sono anche WILD) 
# Da capire come usarli e inserirle nella tabella



prune_samples(sample_data(ps_pooled)$Matrix == "Fish" , ps_Obs_rel)

# trasformare in percentuale
ps_Obs_rel <- transform_sample_counts(ps_pooled, function(x) (x / sum(x))*100)
ps_Obs_rel <- prune_taxa(taxa_sums(ps_Obs_rel) > 0,ps_Obs_rel)


##################
# removing chloroplast and Mitochondria 
##################
ps_Chloroplast <- ps_Obs_rel %>% subset_taxa( Order  == "Chloroplast" | Family  == "Mitochondria")

ps_Phylum_Cl <- tax_glom(ps_Chloroplast, taxrank="Phylum", NArm=FALSE)
plot(sample_sums(ps_Chloroplast))
# plot_bar(ps_Chloroplast, fill = "Phylum")
p1<-plot_bar(ps_Phylum_Cl, fill = "Phylum")+
  facet_grid(~Timepoint*Sample_type, scales = "free", space = "free")
p1

ps_no_Chloroplast <-ps_Obs_rel %>% subset_taxa( Order  != "Chloroplast" & Family  != "Mitochondria")
ps_no_Chloroplast <- prune_taxa(taxa_sums(ps_no_Chloroplast) > 0, ps_no_Chloroplast)
mean(sample_sums(ps_no_Chloroplast))
ps_no_Chloroplast

ps_Obs_rel <- transform_sample_counts(ps_no_Chloroplast, function(x) (x / sum(x))*100)
ps_Obs_rel
#############
# Variabili da usare
 
"Pool_code"          "Timepoint"          "Stage"             
"Sample_type"        "Type_env"           "Country"           
"WP"                 "food_chain"         "sample_site"       
"sample_type_add"    "sample_type_number" "phase"

################

ps_Phylum<- tax_glom(ps_Obs_rel, taxrank="Phylum", NArm=FALSE)
top.names = names(sort(taxa_sums(ps_Phylum), TRUE)[1:10])
top = prune_taxa(top.names,ps_Phylum)

p1<-plot_bar(top, fill = "Phylum")+
  facet_grid(~Timepoint*Sample_type, scales = "free", space = "free")
p1
########################
# Plot rarefaction curves
##############
otu <- as(otu_table(ps_Obs_meta), "matrix")
if(taxa_are_rows(ps_Obs_meta)){
  otu <- t(otu)
}

meta <- data.frame(sample_data(ps_Obs_meta))
site_colors <- as.factor(meta$Sample_type)
rarecurve(
  otu,
  step = 100,
  sample = min(rowSums(otu)),
  col = site_colors,
  label = FALSE
)




#######################
# Alpha div
######################

##################
# pool samples from the same fish
##################
# create new column for pool each fish sample
# Assuming 'physeq' is your phyloseq object
metadata_Farm <- sample_data(ps_Obs_rel)
# Merging 'Location' and 'Time' columns to create a new column 'Location_Time'
metadata_Farm$Code_pool_Alpha <- paste(metadata_Farm$Timepoint, metadata_Farm$Matrix, 
                                       metadata_Farm$Country, metadata_Farm$sample_type_number,
                                       sep = "_")
# Updating the phyloseq object with the modified metadata
sample_data(ps_Obs_rel) <- metadata_Farm
# Assuming 'physeq' is your phyloseq object
first_occurrence_df_Obs <- sample_data(ps_Obs_rel) %>%
  group_by(Code_pool_Alpha) %>%
  dplyr::slice(1) %>%
  select(Code_pool_Alpha, c(1:14)) %>%
  ungroup()

ps_Code_pool_Alpha<- merge_samples(ps_Obs_rel, "Code_pool_Alpha", fun=sum)
sample_sums(ps_Code_pool_Alpha)
ps_Code_pool_Alpha<- transform_sample_counts(ps_Code_pool_Alpha, function(x) ((x / sum(x))*100))

# Replace the sample_data in the phyloseq object
Obs_new_sample_data <- sample_data(as.data.frame(first_occurrence_df_Obs))
sample_names(Obs_new_sample_data) <- Obs_new_sample_data$Code_pool_Alpha
sample_data(ps_Code_pool_Alpha)<-Obs_new_sample_data


###################

#calculate different diversity indices
richness <- specnumber((otu_table(ps_Code_pool_Alpha))) #calculate different diversity indices
shannon <- microbiome::alpha(ps_Code_pool_Alpha, index = 'shannon' )
# anche estimateR funziona, ma in questo caso ha problemi con numeri non interi
shannon$richness <- richness

#convert cluster object to use with ggplot
# l'ordine va dato dopo così non si mescolano i valori 
meta_ps <- meta(ps_Code_pool_Alpha)
meta_ps$shannon <- shannon$diversity_shannon
meta_ps$richness <- shannon$richness


meta_ps$Timepoint <- factor(meta_ps$Timepoint, 
                       levels = c("T0","T1","T2","T3","T4","T5","T6","T7"), 
                       ordered = TRUE)

meta_ps %>%  ggplot(aes(x=as.numeric(Timepoint), y=richness))+
  geom_smooth(aes(color=Matrix), method = "lm", alpha = 0.2)+
  geom_point(aes(color=Matrix), size =2)+
  facet_grid(~Country)+
  #stat_pvalue_manual(pwc, hide.ns = TRUE) +
  theme_bw()

# If data are normal distributed: 
aov(Shannon ~ Group, data = alpha_df)
# se Wilcoxon test for 2 groups
wilcox.test(Shannon ~ Group, data = alpha_df)
# Use Kruskal–Wallis test per multiple groups
kruskal.test(Shannon ~ Group, data = alpha_df)





##################
# pool fish at each timepoint
##################
# create new column for pool each fish sample
# Assuming 'physeq' is your phyloseq object
metadata_Farm <- sample_data(ps_Obs_rel)
# Merging 'Location' and 'Time' columns to create a new column 'Location_Time'
metadata_Farm$Code_pool_Type <- paste(metadata_Farm$Timepoint, metadata_Farm$Matrix, 
                                       metadata_Farm$Country, metadata_Farm$Sample_type,
                                       sep = "_")
# Updating the phyloseq object with the modified metadata
sample_data(ps_Obs_rel) <- metadata_Farm
# Assuming 'physeq' is your phyloseq object
first_occurrence_df_Obs <- sample_data(ps_Obs_rel) %>%
  group_by(Code_pool_Type) %>%
  dplyr::slice(1) %>%
  select(Code_pool_Type, c(1:14)) %>%
  ungroup()

ps_Code_pool_Type<- merge_samples(ps_Obs_rel, "Code_pool_Type", fun=sum)
sample_sums(ps_Code_pool_Type)
ps_Code_pool_Type<- transform_sample_counts(ps_Code_pool_Type, function(x) ((x / sum(x))*100))

# Replace the sample_data in the phyloseq object
Obs_new_sample_data <- sample_data(as.data.frame(first_occurrence_df_Obs))
sample_names(Obs_new_sample_data) <- Obs_new_sample_data$Code_pool_Type
sample_data(ps_Code_pool_Type)<-Obs_new_sample_data


###################
# NMDS and PCoA
###############
ps_fish<-prune_samples(sample_data(ps_Obs_rel)$Matrix == "Fish" ,  ps_Obs_rel)
ps_env<-prune_samples(sample_data( ps_Obs_rel)$Matrix != "Fish" ,  ps_Obs_rel)
ps_It<-prune_samples(sample_data( ps_Obs_rel)$Country == "Italy" ,  ps_Obs_rel)
ps_Gr<-prune_samples(sample_data( ps_Obs_rel)$Country == "Greece" ,  ps_Obs_rel)
ps_Sp<-prune_samples(sample_data( ps_Obs_rel)$Country == "Spain" ,  ps_Obs_rel)
ps_Adult<-prune_samples(sample_data( ps_fish)$Stage == "Adults", ps_fish)

adonis2(phyloseq::distance(ps_Adult, method = "bray")~Country, 
        data = as(sample_data(ps_Adult), "data.frame"),
        permutations = 999)



'
dist_Obs_w <- phyloseq::distance(ps_Obs_rel, method = "bray")
dist_Obs_NMDS <- phyloseq::distance(
  subset_samples(ps_Obs_rel, sample_names(ps_Obs_rel) != "T5_Adults_10_Skin_Spain_Farm")
  , method = "bray")

NMDS_Obs <- ordinate(subset_samples(ps_Obs_rel, sample_names(ps_Obs_rel) != "T5_Adults_10_Skin_Spain_Farm"), 
                      dist_Obs_NMDS, method = "NMDS", trymax=100)
PcoA_Obs_w <- ordinate(ps_Obs_rel, dist_Obs_w, method = "PCoA", trymax=100)
'

 plot_ordination(
   ps_Obs_rel, # for PcOA
    #subset_samples(ps_Obs_rel, sample_names(ps_Obs_rel) != "T5_Adults_10_Skin_Spain_Farm"),
   PcoA_Obs_w, title="Bray_PCoA") + 
   #NMDS_Obs, title="Bray_PCoA")+
   geom_point(aes(shape=Country, color =Sample_type), alpha = 1, size = 3)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  scale_color_brewer(palette = "Accent")  

 
 adonis2(dist_Obs_w~Country, 
              data = as(metadata_Farm, "data.frame"),
              permutations = 999)

 
 Ad<-adonis2(dist_Obs_w~Matrix + Country  + Timepoint, 
                             data = as(metadata_Farm, "data.frame"),
                             permutations = 999)

 Ad$R2
 
 ################
# Agglomerate at Order level
#################

 


physeq_cl <- tax_glom(ps_Obs_rel, taxrank = "Class",NArm=FALSE)
physeq_fam <- tax_glom(ps_Obs_rel, taxrank = "Family",NArm=FALSE)
physeq_gen <- tax_glom(ps_Obs_rel, taxrank = "Genus",NArm=FALSE)

ps_Code_pool_Type
ps_Type_cl <- tax_glom(ps_Code_pool_Type, taxrank = "Class",NArm=FALSE)
ps_Type_fam <- tax_glom(ps_Code_pool_Type, taxrank = "Family",NArm=FALSE)
ps_Type_gen <- tax_glom(ps_Code_pool_Type, taxrank = "Genus",NArm=FALSE)

ps_Type_cl <- tax_glom(ps_Code_pool_Type, taxrank = "Class",NArm=FALSE)

########################
# Barplot Class on the ps pooled for timepoints
###################
ps_to_use <- ps_Type_cl
# Identify top 10
top10_taxa <- names(
  sort(taxa_sums(ps_to_use), decreasing = TRUE)[1:10])
physeq_top10 <- prune_taxa(top10_taxa, ps_to_use)

# Merge all other taxa as “Other”
physeq_other <- prune_taxa(!taxa_names(ps_to_use) %in% top10_taxa,
                           ps_to_use)

physeq_other <- merge_taxa(
  physeq_other,
  taxa_names(physeq_other),
  archetype = taxa_names(physeq_other)[1]
)
tax_table(physeq_other)[1, "Class"] <- "Other"
physeq_final <- merge_phyloseq(physeq_top10, physeq_other)  



p1<-plot_bar(physeq_final, x = "Sample_type", fill = "Class")+
  facet_nested(~Matrix+Country+Timepoint,  scales = "free_x",
             space = "free_x") +
  scale_fill_brewer(palette = "Set3")  +
  theme(
    legend.position = "right",
    panel.spacing.x = unit(0.1, "lines")  # reduce horizontal spacing
  )
p1
##############################################


#####################
# ASV variation in Time point
#################################

#  plot ASV summary
# create_empty list
n_ASV <- list()
n_newASV<- list()
n_lostASV<- list()
n_sharedASV<-list()
t_taxa_old<-list()

########################
# all samples
########################

ps_Fish <- prune_samples(sample_data(ps_Obs_rel)$Matrix == "Fish", ps_Obs_rel)
ps_Fish 

ps_Env <- prune_samples(sample_data(ps_Obs_rel)$Matrix != "Fish", ps_Obs_rel)

for (t in unique(ps_Env@sam_data$Timepoint)) {
  print(t)
  ps_t <- prune_samples(ps_Env@sam_data$Timepoint ==t, ps_Env)  # Subset phyloseq
  t_taxa = taxa_names(prune_taxa(taxa_sums(ps_t) > 0.00001, ps_t))                     # remove singleton
  n_ASV <- append(n_ASV, eval(length(t_taxa)))                                    # list of n ASV for timepoint
  n_sharedASV<- append(n_sharedASV, eval(length(Reduce(intersect, list(t_taxa_old,t_taxa)))))             # shared ASV
  nr_shared_ASV<-length(Reduce(intersect, list(t_taxa_old,t_taxa)))                 # create shared ASV for the 2 timepoint
  print(nr_shared_ASV)                                                                  
  nr_new_ASV <- eval(length(t_taxa))-nr_shared_ASV                                # create new ASV for the timepoint
  print(nr_new_ASV)                                                                    # print new ASV for the timepoint      
  n_newASV <- append(n_newASV, nr_new_ASV)                                            # add to the list new ASV for the timepoint
  nr_lost_ASV <- nr_shared_ASV - eval(length(t_taxa_old))                              # create new ASV for the timepoint
  print(nr_lost_ASV)                                                                    # print new ASV for the timepoint      
  n_lostASV<- append(n_lostASV,nr_lost_ASV)        
  t_taxa_old<-t_taxa                                                               # create phyloseq to compare
}

# Defining list
ls1 <- list(n_ASV, n_sharedASV, n_newASV, n_lostASV)

# Convert list to matrix
Tab_ASV <- matrix(unlist(ls1), nrow = 8, byrow = FALSE)
Tab_ASV_df <- as.data.frame(Tab_ASV)
colnames(Tab_ASV_df)<-c("n_ASV","sharedASV","newASV","lostASV")
Tab_ASV_df$Timepoint<- c("T0","T1","T2","T3","T4","T5","T6","T7")

# create the combined plot
p <- ggplot(Tab_ASV_df,aes(x=Timepoint)) +
  geom_line(aes( y=n_ASV,colour = "n_ASV",group = 1)) +
  geom_line(aes( y=sharedASV, colour = "sharedASV", group = 1)) +
  geom_col(aes(y = newASV, fill = "newASV"), 
           alpha = 0.5, width = 0.6) +
  geom_col(aes(y = lostASV, fill = "lostASV"), 
           alpha = 0.5, width = 0.6) +
  scale_colour_manual("", 
                      values = c("n_ASV"="black", "sharedASV"="red", "newASV" = "green", "lostASV"="red")) +
  labs(title = "ASVs trend along timepoint in Env samples", x = "Time point", y = "nr ASV in Environment") +
  theme_minimal()
p


##########

#  Abbondanza tassonomica per taxa e gruppo
# per PAPER

##############

ps_fish<-prune_samples(sample_data(ps_Obs_rel)$Matrix == "Fish" ,  ps_Obs_rel)
ps_env<-prune_samples(sample_data( ps_Obs_rel)$Matrix != "Fish" ,  ps_Obs_rel)
ps_It<-prune_samples(sample_data( ps_Obs_rel)$Country == "Italy" ,  ps_Obs_rel)
ps_Gr<-prune_samples(sample_data( ps_Obs_rel)$Country == "Greece" ,  ps_Obs_rel)
ps_Sp<-prune_samples(sample_data( ps_Obs_rel)$Country == "Spain" ,  ps_Obs_rel)
ps_Adult<-prune_samples(sample_data( ps_fish)$Stage == "Adults", ps_fish)
ps_early<-prune_samples(sample_data( ps_fish)$Stage != "Adults", ps_fish)
ps_wat<-prune_samples(sample_data( ps_env)$Sample_type != "Feed" ,  ps_env)

mean(sample_sums(subset_taxa(
  prune_samples(sample_data( ps_fish)$Type_env =="Hatchery" , ps_fish), 
 Family == "Flavobacteriaceae")))
sd(sample_sums(subset_taxa(
  prune_samples(sample_data( ps_fish)$Type_env =="Hatchery", ps_fish), 
  Family == "Flavobacteriaceae")))

mean(sample_sums(subset_taxa(ps_wat, Class == "Bacteroidia")))
sd(sample_sums(subset_taxa(ps_wat, Class == "Bacteroidia")))


mean(sample_sums(subset_taxa(
  prune_samples(sample_data( ps_fish)$Type_env =="Sea_cage" &
                  sample_data( ps_fish)$Sample_type =="Gut",
                ps_fish), 
  Class == "Brevinematia")))
sd(sample_sums(subset_taxa(
  prune_samples(sample_data( ps_fish)$Type_env =="Sea_cage" &
                  sample_data( ps_fish)$Sample_type =="Gut",
                ps_fish), 
  Class == "Brevinematia")))

####################
# circos per ASV-taxa sharing tra gruppi
######################

library(dplyr)
library(tidyr)
library(tibble)
library(circlize)
library(ggplot2)
library(ggalluvial)

################
# pool it

# Assuming 'physeq' is your phyloseq object
# create new phyloseq
###############
ps_Alluvial<-ps_Code_pool_Type
# substitute tissue with fish
sample_data(ps_Alluvial)$Sample_type <- sample_data(ps_Alluvial)$Sample_type %>%
  as.character() %>%
  replace(. %in% c("Gill","Skin","Gut"), "fish")
  
metadata_pool <- sample_data(ps_Alluvial)
# Merging 'Location' and 'Time' columns to create a new column 'Location_Time'
metadata_pool$Pool_shar <- paste(metadata_pool$Timepoint, 
        metadata_pool$Country, metadata_pool$Sample_type,
                                      sep = "_")
# Updating the phyloseq object with the modified metadata
sample_data(ps_Alluvial) <- metadata_pool
# Assuming 'physeq' is your phyloseq object
first_occurrence_df_Obs <- sample_data(ps_Alluvial) %>%
  group_by(Pool_shar) %>%
  dplyr::slice(1) %>%
  select(Pool_shar, c(1:14)) %>%
  ungroup()

ps_all<- merge_samples(ps_Alluvial, "Pool_shar", fun=sum)
sample_sums(ps_all)
ps_all<- transform_sample_counts(ps_all, function(x) ((x / sum(x))*100))

# Replace the sample_data in the phyloseq object
new_sample_data <- sample_data(as.data.frame(first_occurrence_df_Obs))
sample_names(new_sample_data) <- new_sample_data$Pool_shar
sample_data(ps_all)<-new_sample_data

##################

ps_all_gen<- tax_glom(ps_all, taxrank = "Genus", NArm = FALSE)

ps_Adult_all<-prune_samples(sample_data(ps_Code_pool_Type)$Stage == "Adults", ps_Code_pool_Type)
ps_A_It <-prune_samples(sample_data(ps_Adult_all)$Country == "Italy", ps_Adult_all)

ps_It <-prune_samples(sample_data(ps_all_gen)$Country == "Italy", ps_all_gen)
ps_Gr <-prune_samples(sample_data(ps_all_gen)$Country == "Greece", ps_all_gen)


######################
# set the phyloseq to use
physeq <- ps_Gr

df <- psmelt(physeq)
df <- df %>% filter(Abundance > 0)

df_sum <- df %>%
  group_by(Genus, Sample_type, Timepoint) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# remove the NAs
df_sum$Genus[is.na(df_sum$Genus)] <- "Unknown"


df_alluvial <- df_sum %>%
  mutate(
    Feed_val = ifelse(Sample_type == "Feed", Abundance, 0),
    #Fish_val = ifelse(Sample_type %in% c("Gill","Gut","Skin"), Abundance, 0),
    Fish_val = ifelse(Sample_type %in% c("Eggs","Larvae","Fry", 
                 "Juveniles", "fish"), Abundance, 0),
    Water_val = ifelse(Sample_type == "Rearing_water", Abundance, 0)
  ) %>%
  group_by(Timepoint, Genus) %>%
  summarise(
    Feed_val = sum(Feed_val, na.rm = TRUE),
    Fish_val = sum(Fish_val, na.rm = TRUE),
    Water_val = sum(Water_val, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Feed_val + Fish_val + Water_val > 0)

#df_filtered <- df_alluvial %>%
#  filter(rowSums(select(., Feed_val, Fish_val, Water_val) == 0) < 2)

# threshold = 1%
df_filtered2 <- df_alluvial %>%
  filter((Feed_val >= 1 & Fish_val >= 1 )|  (Fish_val >= 1 & Water_val >= 1))

sankey_df <- df_filtered2 %>%
  pivot_longer(cols = c(Feed_val, Fish_val, Water_val),
               names_to = "Stage",
               values_to = "Abundance") %>%
  mutate(Stage = factor(Stage, levels = c("Feed_val", "Fish_val", "Water_val")))

sankey_df_filtered <- sankey_df %>%
   filter(Abundance >= 1)

ggplot(sankey_df_filtered,
       aes(x = Stage, stratum = Genus, alluvium = Genus,
           y = Abundance, fill = Genus, label = Genus)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "white") +
  geom_stratum(alpha = 0.8, color = "white", size = 0) +
  facet_wrap(~Timepoint) +
  theme_minimal() +
  #scale_fill_brewer(palette = "Set3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "right")
#######################

###############
## Venn
############



dev.off()
library(VennDiagram)
library(venn)
library(ggvenn)

# Sankey for shared across feed and fish 
# farm_Fry <- tax_glom(farm_Fry, taxrank="Genus",  NArm = FALSE)

set1 <- c(colnames(otu_table(ps_Corr)[, colSums(otu_table(subset_samples(ps_Corr, sample_data(ps_Corr)$Type_env == "Wild"))) != 0]))
set2 <- c(colnames(otu_table(ps_Corr)[, colSums(otu_table(subset_samples(ps_Corr, sample_data(ps_Corr)$Type_env == "Sea_cage"))) != 0]))
#x <- list(Eggs = set1, Larvae = set2, Fry = set3, Juveiles = set4)
x <- list(Wild = set1, Farmed = set2)

set1 <- c(colnames(otu_table(ps_Corr_gut)[, colSums(otu_table(subset_samples(ps_Corr_gut, sample_data(ps_Corr_gut)$Type_env == "Wild"))) != 0]))
set2 <- c(colnames(otu_table(ps_Corr_gut)[, colSums(otu_table(subset_samples(ps_Corr_gut, sample_data(ps_Corr_gut)$Type_env == "Sea_cage"))) != 0]))
#x <- list(Eggs = set1, Larvae = set2, Fry = set3, Juveiles = set4)
x <- list(Wild = set1, Farmed = set2)

ggvenn(x, fill_color = c("lightgreen", "lightblue"),
       stroke_size = 0.5, set_name_size = 4)

##facet_grid()##########################
##   Lefse Analysis
############################

library(microeco)
library(dplyr)
library(magrittr)
library(tidyr)

otu_t<-as.data.frame(t(otu_table(ps_Corr)))
taxa_t<-as.data.frame(tax_table(ps_Corr))
sample_t <- as.data.frame(as.matrix(sample_data(ps_Corr)))
class(sample_t)

# make the taxonomic information unified, very important

taxa_t <- tidy_taxonomy(taxa_t)
dataset <- microtable$new(otu_table = otu_t, sample_table = sample_t, tax_table = taxa_t)
dataset
t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Type_env", 
                     taxa_level = "all", filter_thres = 0.001,  alpha = 0.01, lefse_subgroup = NULL, p_adjust_method = "none")
t1$plot_diff_bar(threshold = 4)
t1$res_diff %<>% subset(Significance %in% "***" )
# we show 20 taxa with the highest LDA (log10)
t1$plot_diff_bar(use_number = 1:30, width = 0.8)

# filter something not needed to show
t1$res_diff %<>% subset(Significance %in% "***" )
t1$plot_diff_abund(use_number = 1:30, add_sig = T, add_sig_label = "Significance")
t1$plot_diff_cladogram(use_taxa_num = 100, use_feature_num = 30, clade_label_level = 6)

#######################
# Netwotk with Core
##############
# Load the required libraries
library(phyloseq)
library(igraph)
library(ggplot2)
library(ggraph)


ps_Corr2 <- prune_taxa(taxa_sums(ps_Corr) > 1, ps_Corr)
ps_Corr_Gen = tax_glom(ps_Corr , taxrank="Genus", NArm=FALSE)
sum(sample_sums(ps_Corr2))/sum(sample_sums(ps_Corr))
#  [1] 0.9048743

ps_Corr_gut2 <- prune_taxa(taxa_sums(ps_Corr_gut) > 1, ps_Corr_gut)
ps_Corr_gut2_Gen = tax_glom(ps_Corr_gut2 , taxrank="Genus", NArm=FALSE)


# Sample data (metadata)
sample_data <- sample_data(ps_Corr_gut2)
# OTU/ASV (taxa) table
otu_table <- otu_table(ps_Corr_gut2)

# Convert OTU table to presence/absence data
otu_pa <- otu_table > 0 # TRUE if taxa is present in a sample

# Create a bipartite graph where samples and taxa are nodes, and edges represent presence of taxa in samples
network_graph <- graph_from_biadjacency_matrix(otu_pa)

# Add metadata to sample nodes
V(network_graph)$type <- c(rep("Sample", nsamples(ps_Corr_gut2)), rep("Taxa", ntaxa(ps_Corr_gut2)))

sample_data$WP <- as.factor(sample_data$WP)
# Add variable for coloring samples (e.g., from metadata)
V(network_graph)$color <- c(sample_data$WP, rep(NA, ntaxa(ps_Corr_gut2)))
#V(network_graph)$shape <- c(sample_data$Sample_type, rep("Taxa", ntaxa(ps_Corr2)))

# Size of taxa nodes based on abundance percentage
taxa_abundances <- taxa_sums(ps_Corr_gut2) / sum(taxa_sums(ps_Corr_gut2)) * 100
V(network_graph)$size <- c(rep(3, nsamples(ps_Corr_gut2)), taxa_abundances)  # Sample nodes get fixed size (e.g., 3)


# Create the network plot using ggraph
ggraph(network_graph, layout = "auto") +  # Choose a layout, e.g., Fruchterman-Reingold
  geom_edge_link(aes(), alpha = 0.6, color = "gray90") +  # Edges connecting samples and taxa
  geom_node_point(aes(color = factor(color), size = size), alpha = 0.8) +  # Nodes: samples and taxa, size and color as per metadata
  #scale_color_manual(values = c(Wild = "blue", "Farm" = "green", "taxa" = "red")) +  # Customize colors for nodes
  scale_shape_manual(values = c(15, 16,17,18))+
  theme_void() +  # Clean background
  labs(title = "Network of Samples and Taxa", color = "Node Type", size = "Abundance %") +
  theme(legend.position = "right")


##################
# core 
################
# Assume `physeq` is your phyloseq object
# Set a threshold for prevalence (e.g., present in at least 50% of the samples)
prevalence_threshold <- 0.75  # This means present in 50% of samples
# Set a relative abundance threshold (e.g., minimum relative abundance of 0.001 or 0.1%)
abundance_threshold <- 0.1
# Calculate core taxa based on prevalence and abundance
core_taxa_Obs_Gut <- core(ps_Corr_gut2_Gen, detection = abundance_threshold, prevalence = prevalence_threshold)
# View the core taxa in the core microbiome
core_taxa_Obs_Gut

plot_bar(core_taxa_Obs_Gut, fill = "Genus")+
  facet_grid(~Timepoint*Sample_type, scales = "free", space = "free")



#######################
# Taxa modeling check
#######################

ps_modeling <- subset_taxa(ps_Corr_gut,   Genus == "Pir4 lineage"|
                             Genus == "Planctomyces"| # NO
                             Genus == "Ochrobactrum"|
                             Genus == "Rhodopirellula"|
                             Genus == "Vibrio"|
                             Genus == "Endozoicomonas"|
                             Genus == "Legionella"|
                             Genus == "Exiguobacterium" |
                             Genus == "Clostridium sensu stricto 1"|
                             Genus == "Bacillus" |
                             Genus == "Cetobacterium")

ps_modeling_neg <- subset_taxa(ps_Corr_gut, Genus == "Bacteroides"|
                                 Genus == "Escherichia Shigella"|
                                 Genus == "Fusobacterium"|
                                 Genus == "Acinetobacter"|
                                 Genus == "Aeromonas"|
                                 Genus == "Flavobacterium")

ps_modeling_clos <- subset_taxa(ps_Corr_gut,  Family == "Clostridiaceae")
ps_modeling_vibrio <- subset_taxa(ps_Corr_gut,  Genus == "Vibrio")
ps_modeling_bac <- subset_taxa(ps_Corr_gut, Genus == "Bacteroides")

|Genus == "Endozoicomonas"
|Genus == "Exiguobacterium"
|Genus == "Idiomarina"
|Genus == "Sporosarcina"
|Genus == "Vibrio")

|
  |Genus == "Lactococcus"
|Genus == "Streptococcus"
|Genus == "Lactobacillus")
|Genus == "Bacillus"
|Genus == "Bacteroides"
|Genus == "Eubacteriumventriosumgroup"
|Genus == "Lactobacillus"
|Genus == "Streptococcus"
|Genus == "Blastopirellula"|
  Genus == "Mycobacterium"|
  Genus == "Pir4 lineage"|
  Genus == "Planctomyces"|
  Genus == "Rhodopirellula"




ps_modeling_gen <- tax_glom(ps_modeling, "Genus",NArm=FALSE)

ps_modeling_clos_gen <- tax_glom(ps_modeling_clos, "Genus")
ps_modeling_vibrio_gen <- tax_glom(ps_modeling_vibrio, "Genus")
ps_modeling_bac_gen <- tax_glom(ps_modeling_bac, "Genus")
ps_modeling_neg_gen <- tax_glom(ps_modeling_neg, "Genus")

plot_bar(ps_modeling_neg_gen, fill = "Genus")+
  facet_grid(~WP, scales = "free", space = "free")+
  theme(  axis.text.x = element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())


# psmelt
phyloseq::psmelt(ps_modeling_gen) %>%
  ggplot(data = ., aes(x = Pool, y = Abundance)) +
  geom_bar(aes(fill = OTU), outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_grid(~Diet*Time, scales = "free", space = "free")
#facet_wrap(~ OTU, scales = "free")


sample_alpha_env <- meta(ps_modeling_gen)
sample_alpha_env$Strepto <- sample_sums(subset_taxa(ps_modeling_gen, Genus == "Streptococcus"))
sample_alpha_env$Exi <- sample_sums(subset_taxa(ps_modeling_gen, Genus == "Exiguobacterium"))
sample_alpha_env$Flav <- sample_sums(subset_taxa(ps_modeling_gen, Genus == "Flavobacterium"))
sample_alpha_env$Bac <- sample_sums(subset_taxa(ps_modeling_gen, Genus == "Bacteroides"))
sample_alpha_env_T1<-subset(sample_alpha_env, Time == "T1")

fit <- aov(Flav~Diet, data=sample_alpha_env_T1)
summary(fit)

"Osservazione su tessuti e digesta di Gut
comparazione di abbondanza relativa tra Timepoint e Diete
Comparazione di abbondanza dei taxa rilevanti dal modeling (WP8)

Diet1-T2 riduzione di Streptococcus, presente invece in tutti gli altri trattamenti
Diet1 e Diet2 riduzione di Bacteroides e Flavobacterium, presente invece in CTRL-T1
Exiguobacterium presente In CTRL(tutti e tre i time point), ma non in Diet1 e Diet2


Streptococcus in Time2 - Anova tra CTRL-Diet1-Diet2 p = 0.02
Exiguobacterium in Time2 - Anova tra CTRL-Diet1-Diet2 p = 0.05
Bacteroides in T1-T2 Anova tra CTRL-Diet1-Diet2 p = 0.008
"
#############################


library(phyloseq)
library(ggplot2)

# Extract OTU table
otu_occ <- as(otu_table(ps_Adult), "matrix")

# Ensure taxa are rows
if (!taxa_are_rows(ps_Adult)) {
  otu_occ <- t(otu_occ)
}

# Number of samples where each taxon is present
occupancy <- rowSums((otu_occ > 0))
sum_ASV<-rowSums((otu_occ))
occ_abb<- sum_ASV/occupancy

# Proportion of samples occupied
occ_df <- data.frame(
  ASV = taxa_names(ps_Adult),
  Taxon = ps_fish@tax_table[,5],
  rel_Abb = occ_abb,
  n_sample = occupancy
)

# PLOT
library(ggrepel)
# Example: label ASVs with occupancy > 0.8
label_df <- subset(occ_df, rel_Abb > 10 | n_sample > 75)

ggplot(occ_df,
       aes(x = n_sample, y = rel_Abb)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(data = label_df, aes(label = Family), max.overlaps = Inf)+
  scale_y_log10() +
  theme_bw()
