library(tidyverse)
library(cowplot)

## Load the parameters from snakemake
input_dir <- snakemake@params$input_dir
output_figure <- snakemake@output$figure
output_table <- snakemake@output$table


## Read the metaquast summary for TM mixture
read_metaquast_tm <- function(file_name, dir = input_dir) {
  file <- paste(dir, file_name, sep = "/")
  TM.metaquast <- read.table(file, sep = "\t", header = T, na.strings = "-", check.names = F) %>%
    separate(Assemblies, c("sample", "assembler"), sep = "\\.") 
  TM.metaquast$`TB40E` <- ifelse(TM.metaquast$sample == "TM_0_1", NA, TM.metaquast$`TB40E`)
  TM.metaquast$Merlin <- ifelse(TM.metaquast$sample == "TM_1_0", NA, TM.metaquast$Merlin)
  TM.metaquast$not_aligned <- NULL
  return(TM.metaquast)
}

## Read the metaquast summary for TM mixture
read_metaquast_ta <- function(file_name, dir = input_dir) {
  file <- paste(dir, file_name, sep = "/")
  TA.metaquast <- read.table(file, sep = "\t", header = T, na.strings = "-", check.names = F) %>%
    separate(Assemblies, c("sample", "assembler"), sep = "\\.") 
  TA.metaquast$`TB40E` <- ifelse(TA.metaquast$sample == "TA_0_1", NA, TA.metaquast$`TB40E`)
  TA.metaquast$AD169 <- ifelse(TA.metaquast$sample == "TA_1_0", NA, TA.metaquast$AD169)
  TA.metaquast$not_aligned <- NULL
  return(TA.metaquast)
}


## Nicer theme with and without legend
theme_with_legend <- function(){
  p<- guides(color = guide_legend(
    override.aes = list(size = 5),
    label.theme = element_text(size = 13)
  )) 
  return(p)
}

theme_wo_legend <- function(){
  p <- theme(axis.text.x = element_text(
      angle = 70,
      hjust = 1,
      size = 15
    ), legend.position = "none")
  return(p)
}


get_portion <- function(sample_summary) {
  sample_portion <- as.data.frame(do.call(
    rbind,
    strsplit(sample_summary$sample, split = "_")
  ))
  sample_portion$ref <- substr(sample_summary$reference, 1, 1)
  portion <- ifelse(substr(sample_portion[, 1], 1, 1) == sample_portion$ref,
                    paste(sample_portion[, 2], sample_portion[, 3], sep = "/"),
                    paste(sample_portion[, 3], sample_portion[, 2], sep = "/")
  )
  sample_summary$ratio <- portion
  sample_summary <- filter(
    sample_summary,
    ratio != "0",
    ratio != "0/1"
  )
}

color_dot_boxplot <- function(df, value, label) {
  g <- ggplot(get_portion(df), aes_string("assembler", value)) +
    geom_boxplot(show.legend = F, fill = "grey85") +
    geom_point(
      position = position_jitterdodge(
        jitter.width = 0,
        dodge.width = 0.3,
        seed = 1234
      ),
      aes(color = ratio),
      stroke = 0,
      size = 2
    ) + xlab("") + 
    theme_bw(base_size = 15) +
    ylab(label) + 
    scale_color_brewer(palette = "Set1")
  return(g)
}



## The visualized criteria
metaquast_criteria <- c(
  "num_contigs", "Largest_contig", "Genome_fraction",
  "Duplication_ratio", "Largest_alignment", "LGA50",
  "NGA50", "num_misassemblies", "num_mismatches_per_100_kbp"
)



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

## Number of mismatches per 100 kbp for TM and TA
TM.num_mismatches_per_100_kbp <- read_metaquast_tm("TM.num_mismatches_per_100_kbp.merged.tsv")
TM.num_mismatches_per_100_kbp.long <- gather(
  TM.num_mismatches_per_100_kbp, reference, mismatches,
  -c(sample, assembler)
)
TA.num_mismatches_per_100_kbp <- read_metaquast_ta("TA.num_mismatches_per_100_kbp.merged.tsv")
TA.num_mismatches_per_100_kbp.long <- gather(
  TA.num_mismatches_per_100_kbp, reference, mismatches,
  -c(sample, assembler)
)

## Genome fraction visualization
p_gf_box <- color_dot_boxplot(rbind(TM.GenomeFraction.long, 
                                    TA.GenomeFraction.long), 
                  "fraction", 
                  "% Genome fraction") +
  theme_wo_legend()
  



## Duplication ratio visualization
p_dr_box <- color_dot_boxplot(rbind(TM.DuplicationRatio.long, 
                                    TA.DuplicationRatio.long), 
                  "DuplicationRatio", 
                  "Duplication Ratio") +
  theme_wo_legend()


## Largest alignment visualization
p_la_box <- color_dot_boxplot(rbind(TM.Largest_alignment.long, 
                        TA.Largest_alignment.long), 
                  "LargestAlignment / 1000", 
                  "Largest Alignment (kb)") +
  theme_wo_legend()


## Largest contig visualization
p_lc_box <- color_dot_boxplot(rbind(TM.Largest_contig.long, 
                        TA.Largest_contig.long), 
                  "LargestContig / 1000", 
                  "Largest Contig (kb)") +
  theme_wo_legend()

## LGA50 visualization
p_lga50_box <- color_dot_boxplot(rbind(TM.LGA50.long, 
                        TA.LGA50.long), 
                  "LGA50", 
                  "LGA50") +
  theme_wo_legend()

## NGA50 visualization
p_nga50_box <- color_dot_boxplot(rbind(TM.NGA50.long, 
                        TA.NGA50.long), 
                  "NGA50 / 1000", 
                  "NGA50 (kb)") +
  theme_wo_legend()

## Number of contigs visualization
p_c_box <- color_dot_boxplot(rbind(TM.contigs.long, 
                        TA.contigs.long), 
                  "contigs", 
                  "# Contigs") +
  theme_wo_legend()

## Only legend for all figures
p_c_box_legend <- color_dot_boxplot(rbind(TM.contigs.long, 
                                          TA.contigs.long), 
                                    "contigs", 
                                    "# Contigs") +
  theme_with_legend()


## Number of misassemblies visualization
p_mis_box <- color_dot_boxplot(rbind(TM.misassemblies.long, 
                        TA.misassemblies.long), 
                  "misassemblies", 
                  "# Misassemblies") +
  theme_wo_legend()


## Number of mismatches per 100kbp visualization
p_mism_box <- color_dot_boxplot(rbind(TM.num_mismatches_per_100_kbp.long, 
                        TA.num_mismatches_per_100_kbp.long), 
                  "mismatches", 
                  "# Mismatches/100kbp") +
  theme_wo_legend()


## Get the legend
legend <- get_legend(p_c_box_legend)

## Combine all visualization and legend
p <- plot_grid(p_c_box, p_lc_box, p_gf_box, p_dr_box, p_la_box,
  p_lga50_box, p_nga50_box, p_mis_box, p_mism_box,
  ncol = 3, labels = "AUTO"
)

assembly_evaluation_plot <- plot_grid(p, legend, rel_widths = c(3, 0.5))

ggsave(filename = output_figure, plot = assembly_evaluation_plot, width = 11, height = 11)



## Rank the assemblers
small_better <- function(df) {
  criteria <- colnames(df)[4]
  final_rank <- df %>%
    filter(get(criteria) != "NA") %>%
    group_by(sample, reference) %>%
    mutate(rank = rank(!!as.name(criteria), ties.method = "first")) %>%
    group_by(assembler) %>%
    summarise(rank = mean(rank)) %>%
    mutate(rank = rank(rank, ties.method = "first"))
}

big_better <- function(df) {
  criteria <- colnames(df)[4]
  final_rank <- df %>%
    filter(get(criteria) != "NA") %>%
    group_by(sample, reference) %>%
    mutate(rank = rank(-!!as.name(criteria), ties.method = "first")) %>%
    group_by(assembler) %>%
    summarise(rank = mean(rank)) %>%
    mutate(rank = rank(rank, ties.method = "first"))
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


MA.rank <- small_better(rbind(TM.misassemblies.long, 
                        TA.misassemblies.long))
MA.rank$criteria <- "Number of misassembled"

NC.rank <- small_better(rbind(TM.contigs.long, TA.contigs.long))
NC.rank$criteria <- "Number of contigs"

DR.rank <- small_better(rbind(TM.DuplicationRatio.long, TA.DuplicationRatio.long))
DR.rank$criteria <- "Duplication ratio"


MM.rank <- small_better(rbind(TM.num_mismatches_per_100_kbp.long, TA.num_mismatches_per_100_kbp.long))
MM.rank$criteria <- "Mismatches per 100k"

summary.rank <- rbind(GF.rank, LA.rank, LGA50.rank, LC.rank, NGA50.rank, MA.rank, NC.rank, DR.rank, MM.rank)


summary.rank.ordered <- summary.rank %>%
  group_by(criteria) %>%
  arrange(criteria, rank)

write.table(summary.rank.ordered, file = output_table, sep = "\t")