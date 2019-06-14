library(tidyverse)
library(cowplot)

## Load the parameters from snakemake
input_dir <- snakemake@params$input_dir
output_figure <- snakemake@output$figure
output_table <- snakemake@output$table


## Read the metaquast summary for TM mixture
read_metaquast_tm <- function(input_dir=input_dir, file_name){
  file <- paste(input_dir, file_name, sep="/")
  TM.metaquast <- read.table(file, sep="\t", header=T, na.strings = "-", check.names = F) %>% 
    separate(Assemblies, c("sample", "assembler"), sep="\\.") %>% 
        filter(assembler!="velvet" & assembler!="disco")
  TM.metaquast$`TB40E` <- ifelse(TM.metaquast$sample=="TM_0_1", NA, TM.metaquast$`TB40E`)
  TM.metaquast$Merlin <- ifelse(TM.metaquast$sample=="TM_1_0", NA, TM.metaquast$Merlin)
  TM.metaquast$not_aligned <- NULL
  return(TM.metaquast)
}

## Read the metaquast summary for TM mixture
read_metaquast_ta <- function(input_dir=input_dir, file_name){
  file <- paste(input_dir, file_name, sep="/")
  TA.metaquast <- read.table(file, sep="\t", header=T, na.strings = "-", check.names = F) %>% 
    separate(Assemblies, c("sample", "assembler"), sep="\\.") %>% 
        filter(assembler!="velvet" & assembler!="disco")
  TA.metaquast$`TB40E` <- ifelse(TA.metaquast$sample=="TA_0_1", NA, TA.metaquast$`TB40E`)
  TA.metaquast$AD169 <- ifelse(TA.metaquast$sample=="TA_1_0", NA, TA.metaquast$AD169)
  TA.metaquast$not_aligned <- NULL
  return(TA.metaquast)
}


## Bigger font size theme
theme_big <- function(){
  p <- theme_bw(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "none")
  return(p)
}


theme_big_with_legend <- function(){
  p <- theme_bw(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle = 60, hjust = 1))
  return(p)
}

## The visualized criteria
metaquast_criteria = ["num_contigs", "Largest_contig", "Genome_fraction",
                   "Duplication_ratio", "Largest_alignment", "LGA50",
                   "NGA50", "num_misassemblies", "num_mismatches_per_100_kbp"]



## Genome fraction for TM and TA
TM.GenomeFraction <- read_metaquast_tm("TM.Genome_fraction.merged.tsv")
TM.GenomeFraction.long <- gather(TM.GenomeFraction, reference, fraction, -c(sample, assembler))
TA.GenomeFraction <- read_metaquast_ta("TA.Genome_fraction.merged.tsv")
TA.GenomeFraction.long <- gather(TA.GenomeFraction, reference, fraction, -c(sample, assembler))

## Duplicate ratio for TM and TA
TM.DuplicationRatio <- read_metaquast_tm("TM.Duplication_ratio.merged.tsv")
TM.DuplicationRatio.long <- gather(TM.DuplicationRatio, reference, DuplicationRatio, -c(sample, assembler))
TA.DuplicationRatio <- read_metaquast_ta("TA.Duplication_ratio.merged.tsv")
TA.DuplicationRatio.long <- gather(TA.DuplicationRatio, reference, DuplicationRatio, -c(sample, assembler))

## Largest alignment for TM and TA
TM.Largest_alignment <- read_metaquast_tm("TM.Largest_alignment.merged.tsv")
TM.Largest_alignment.long <- gather(TM.Largest_alignment, reference, LargestAlignment, -c(sample, assembler))
TA.Largest_alignment <- read_metaquast_ta("TA.Largest_alignment.merged.tsv")
TA.Largest_alignment.long <- gather(TA.Largest_alignment, reference, LargestAlignment, -c(sample, assembler))


## Largest contig for TM and TA
TM.Largest_contig <- read_metaquast_tm("TM.Largest_contig.merged.tsv")
TM.Largest_contig.long <- gather(TM.Largest_contig, reference, LargestContig, -c(sample, assembler))
TA.Largest_contig <- read_metaquast_ta("TA.Largest_contig.merged.tsv")
TA.Largest_contig.long <- gather(TA.Largest_contig, reference, LargestContig, -c(sample, assembler))


## Genome aligned L50 for TM and TA
TM.LGA50 <- read_metaquast_tm("TM.LGA50.merged.tsv")
TM.LGA50.long <- gather(TM.LGA50, reference, LGA50, -c(sample, assembler))
TA.LGA50 <- read_metaquast_ta("TA.LGA50.merged.tsv")
TA.LGA50.long <- gather(TA.LGA50, reference, LGA50, -c(sample, assembler))

## Genome aligned N50 for TM and TA
TM.NGA50 <- read_metaquast_tm("TM.NGA50.merged.tsv")
TM.NGA50.long <- gather(TM.NGA50, reference, NGA50, -c(sample, assembler))
TA.NGA50 <- read_metaquast_ta("TA.NGA50.merged.tsv")
TA.NGA50.long <- gather(TA.NGA50, reference, NGA50, -c(sample, assembler))


## Number of contigs for TM and TA
TM.contigs <- read_metaquast_tm("TM.num_contigs.merged.tsv")
TM.contigs.long <- gather(TM.contigs, reference, contigs, -c(sample, assembler))
TA.contigs <- read_metaquast_ta("TA.num_contigs.merged.tsv")
TA.contigs.long <- gather(TA.contigs, reference, contigs, -c(sample, assembler))

## Number of misassemblies for TM and TA
TM.misassemblies <- read_metaquast_tm("TM.num_misassemblies.merged.tsv")
TM.misassemblies.long <- gather(TM.misassemblies, reference, misassemblies, -c(sample, assembler))
TA.misassemblies <- read_metaquast_ta("TA.num_misassemblies.merged.tsv")
TA.misassemblies.long <- gather(TA.misassemblies, reference, misassemblies, -c(sample, assembler))

#TM.Total_aligned_length <- read_metaquast_tm("TM.Total_aligned_length.merged.tsv")
#TM.Total_aligned_length.long <- gather(TM.Total_aligned_length, reference, TotalAlignedLength, -c(sample, assembler))
#TA.Total_aligned_length <- read_metaquast_ta("TA.Total_aligned_length.merged.tsv")
#TA.Total_aligned_length.long <- gather(TA.Total_aligned_length, reference, TotalAlignedLength, -c(sample, assembler))


## Number of mismatches per 100 kbp for TM and TA
TM.num_mismatches_per_100_kbp <- read_metaquast_tm("TM.num_mismatches_per_100_kbp.merged.tsv")
TM.num_mismatches_per_100_kbp.long <- gather(TM.num_mismatches_per_100_kbp, reference, mismatches, 
                                             -c(sample, assembler))
TA.num_mismatches_per_100_kbp <- read_metaquast_ta("TA.num_mismatches_per_100_kbp.merged.tsv")
TA.num_mismatches_per_100_kbp.long <- gather(TA.num_mismatches_per_100_kbp, reference, mismatches, 
                                             -c(sample, assembler))

## Genome fraction visualization
p_gf_box <- ggplot(rbind(TM.GenomeFraction.long, TA.GenomeFraction.long), 
                   aes(assembler, fraction, fill=reference)) + 
  geom_boxplot() +
  ylab("% Genome Fraction") + 
  xlab("") +
  theme_big()


## Duplication ratio visualization
p_dr_box <- ggplot(rbind(TM.DuplicationRatio.long, TA.DuplicationRatio.long), 
                   aes(assembler, DuplicationRatio, fill=reference)) + 
  geom_boxplot() +
  ylab("Duplication Ratio") + 
  xlab("") +
  theme_big() 

## Largest alignment visualization
p_la_box <- ggplot(rbind(TM.Largest_alignment.long, TA.Largest_alignment.long), 
                   aes(assembler, LargestAlignment/1000, fill=reference)) + 
  geom_boxplot() +
  ylab("Largest Alignment (kb)") + 
  xlab("") +
  theme_big()

## Largest contig visualization
p_lc_box <- ggplot(rbind(TM.Largest_contig.long, TA.Largest_contig.long), 
                   aes(assembler, LargestContig/1000, fill=reference)) + 
  geom_boxplot() +
  ylab("Largest Contig (kb)") + 
  xlab("") +
  theme_big()

## LGA50 visualization
p_lga50_box <- ggplot(rbind(TM.LGA50.long, TA.LGA50.long), 
                      aes(assembler, LGA50, fill=reference)) + 
  geom_boxplot() +
  ylab("LGA50") + 
  xlab("") +
  theme_big()

## NGA50 visualization
p_nga50_box <- ggplot(rbind(TM.NGA50.long, TA.NGA50.long), 
                      aes(assembler, NGA50/1000, fill=reference)) + 
  geom_boxplot() +
  ylab("NGA50 (kb)") + 
  xlab("") +
  theme_big()


## Number of contigs visualization
p_c_box <- ggplot(rbind(TM.contigs.long, TA.contigs.long), 
                  aes(assembler, contigs, fill=reference)) + 
  geom_boxplot() +
  ylab("# Contigs") + 
  xlab("") +
  theme_big() 

## Only legend for all figures
p_c_box_legend <- ggplot(rbind(TM.contigs.long, TA.contigs.long), 
                         aes(assembler, contigs, fill=reference)) + 
  geom_boxplot() +
  ylab("# Contigs") + 
  xlab("") +
  theme_big_with_legend() 

## Number of misassemblies visualization
p_mis_box <- ggplot(rbind(TM.misassemblies.long, TA.misassemblies.long), 
                    aes(assembler, misassemblies, fill=reference)) + 
  geom_boxplot() +
  ylab("# Misassemblies") + 
  xlab("") +
  theme_big()

## Number of mismatches per 100kbp visualization
p_mism_box <- ggplot(rbind(TM.num_mismatches_per_100_kbp.long, TA.num_mismatches_per_100_kbp.long), 
                     aes(assembler, mismatches, fill=reference)) + 
  geom_boxplot() +
  ylab("# Mismatches/100kbp") + 
  xlab("") +
  theme_big() 


## Get the legend
legend <- get_legend(p_c_box_legend)

## Combine all visualization and legend
p <- plot_grid(p_c_box, p_lc_box, p_gf_box, p_dr_box, p_la_box, 
          p_lga50_box, p_nga50_box, p_mis_box, p_mism_box, 
          ncol = 3, labels = "AUTO"
)

assembly_evaluation_plot <- plot_grid(p, legend, rel_widths = c(3, 0.5))

ggsave(filename=output_figure, plot=assembly_evaluation_plot, width = 11, height = 11)



## Rank the assemblers
small_better <- function(df){
  
  criteria <- colnames(df)[4]
  final_rank <- df %>% filter(get(criteria) != "NA") %>%
    group_by(sample, reference) %>% mutate(rank = rank(!!as.name(criteria), ties.method = "first")) %>% 
    group_by(assembler) %>% summarise(rank=mean(rank)) %>% mutate(rank=rank(rank, ties.method = "first"))
}

big_better <- function(df){
  criteria <- colnames(df)[4]
  final_rank <- df %>% filter(get(criteria) != "NA") %>%
    group_by(sample, reference) %>% mutate(rank = rank(-!!as.name(criteria), ties.method = "first")) %>% 
    group_by(assembler) %>% summarise(rank=mean(rank)) %>% mutate(rank=rank(rank, ties.method = "first"))
}


GF.rank <- big_better(rbind(TM.GenomeFraction.long, TA.GenomeFraction.long))
GF.rank$criteria <- "Genome fraction"


LA.rank <- big_better(rbind(TM.Largest_alignment.long, TA.Largest_alignment.long))
LA.rank$criteria <- "Largest alignment"

LC.rank <- big_better(rbind(TM.Largest_contig.long, TA.Largest_contig.long))
LC.rank$criteria <- "Largest contig"


LGA50.rank <- small_better(rbind(TM.LGA50.long, TA.LGA50.long))
LGA50.rank$criteria <- "LGA50"

NGA50.rank <- big_better(rbind(TM.NGA50.long, TA.NGA50.long))
NGA50.rank$criteria <- "NGA50"


MA.rank <- small_better(rbind(TM.Misassembled_contigs_length.long, TA.Misassembled_contigs_length.long))
MA.rank$criteria <- "Misassembled contigs length"

NC.rank <- small_better(rbind(TM.contigs.long, TA.contigs.long))
NC.rank$criteria <- "Number of contigs"

DR.rank <- small_better(rbind(TM.DuplicationRatio.long, TA.DuplicationRatio.long))
DR.rank$criteria <- "Duplication ratio"


TAL.rank <- big_better(rbind(TM.Total_aligned_length.long, TA.Total_aligned_length.long))
TAL.rank$criteria <- "Total aligned length"

MM.rank <- small_better(rbind(TM.num_mismatches_per_100_kbp.long, TA.num_mismatches_per_100_kbp.long))
MM.rank$criteria <- "Mismatches per 100k"

summary.rank <- rbind(GF.rank, LA.rank, LGA50.rank, LC.rank, NGA50.rank, MA.rank, NC.rank, DR.rank, TAL.rank, MM.rank)


summary.rank.ordered <- summary.rank %>% group_by(criteria) %>% arrange(criteria, rank)

write.table(summary.rank.ordered, file=output_table, sep="\t")
