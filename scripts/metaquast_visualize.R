library(tidyverse)
library(cowplot)

## Load the parameters from snakemake
input_dir <- snakemake@params$input_dir
output_figure <- snakemake@output$figure
output_table <- snakemake@output$table
output_score <- snakemake@output$score
read_metaquast <- function(file_name, dir = input_dir) {
  file <- paste(dir, file_name, sep = "/")

  metaquast <- read.table(file, sep = "\t", header = T, na.strings = "-", check.names = F) %>%
    separate(Assemblies, c("sample", "assembler"), sep = "\\.") 
  ref1 <- 3
  ref2 <- 4
  if (grepl("NGA50", file_name)){
      metaquast[,ref1] <- ifelse(endsWith(metaquast$sample, "_0_1"), -1, metaquast[,ref1])
      metaquast[,ref2] <- ifelse(endsWith(metaquast$sample, "_1_0"), -1, metaquast[,ref2])
      metaquast[is.na(metaquast)] <- 0
      metaquast[metaquast==-1] = NA
  }
  else{
      metaquast[,ref1] <- ifelse(endsWith(metaquast$sample, "_0_1"), NA, metaquast[,ref1])
      metaquast[,ref2] <- ifelse(endsWith(metaquast$sample, "_1_0"), NA, metaquast[,ref2])
  }
  metaquast$not_aligned <- NULL
  metaquast$assembler <- ifelse(metaquast$assembler=="metaspades", "metaspa", metaquast$assembler)
  return(metaquast)
}

metaquast_criteria <- c(
  "num_contigs", "Genome_fraction",
  "Duplication_ratio", "Largest_alignment", 
  "NGA50", "num_misassemblies", "num_mismatches_per_100_kbp"
)

files <- list.files(input_dir, pattern = "\\.merged\\.tsv$")
metaquast_list = list()
for (file in files){
    mix_criteria = gsub("\\.merged\\.tsv$", "", file)
    criteria = unlist(strsplit(mix_criteria, "\\."))[2]
    if (criteria %in% metaquast_criteria){
        df = read_metaquast(file)
        df$criteria = criteria
        df.long <- gather(df, reference, value, -c(sample, assembler, criteria))
        metaquast_list[[mix_criteria]] = df.long
    }
}

metaquast_all_criteria <- do.call(rbind, metaquast_list)

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

## Number of contigs visualization
p_c_box <- color_dot_boxplot(filter(metaquast_all_criteria, 
                            criteria=="num_contigs"), 
                  "value", 
                  "# Contigs") +
  scale_y_log10() +
  theme_wo_legend()

## Only legend for all figures
p_c_box_legend <- color_dot_boxplot(filter(metaquast_all_criteria, 
                            criteria=="num_contigs"), 
                                    "value", 
                                    "# Contigs") +
  theme_with_legend()


## Genome fraction visualization
p_gf_box <- color_dot_boxplot(filter(metaquast_all_criteria, 
                             criteria=="Genome_fraction"), 
                  "value", 
                  "% Genome fraction") +
  scale_y_log10() +
  theme_wo_legend() 
  


## Duplication ratio visualization
p_dr_box <- color_dot_boxplot(filter(metaquast_all_criteria, 
                               criteria=="Duplication_ratio"), 
                  "value", 
                  "Duplication ratio") +
  theme_wo_legend()


## Largest alignment visualization
p_la_box <- color_dot_boxplot(filter(metaquast_all_criteria, 
                              criteria=="Largest_alignment"), 
                  "value / 1000", 
                  "Largest Alignment (kb)") +
  scale_y_log10() +
  theme_wo_legend()


# ## Largest contig visualization
# p_lc_box <- color_dot_boxplot(rbind(TM.Largest_contig.long, 
#                         TA.Largest_contig.long), 
#                   "LargestContig / 1000", 
#                   "Largest Contig (kb)") +
#   scale_y_log10() +
#   theme_wo_legend()

# ## LGA50 visualization
# p_lga50_box <- color_dot_boxplot(rbind(TM.LGA50.long, 
#                         TA.LGA50.long), 
#                   "LGA50", 
#                   "LGA50") +
#   scale_y_log10() +
#   theme_wo_legend()

## NGA50 visualization
p_nga50_box <- color_dot_boxplot(filter(metaquast_all_criteria, 
                                criteria=="NGA50"), 
                  "value / 1000", 
                  "NGA50 (kb)") +
  theme_wo_legend()



## Number of misassemblies visualization
# p_mis_box <- color_dot_boxplot(rbind(TM.misassemblies.long, 
#                         TA.misassemblies.long), 
#                   "misassemblies", 
#                   "# Misassemblies") +
#   theme_wo_legend()


## Number of mismatches per 100kbp visualization
p_mism_box <- color_dot_boxplot(filter(metaquast_all_criteria, 
                            criteria=="num_mismatches_per_100_kbp"), 
                  "value", 
                  "# Mismatches/100kbp") +
  scale_y_log10() + 
  theme_wo_legend()


## Get the legend
legend <- get_legend(p_c_box_legend)


p6 <- plot_grid(p_c_box, p_la_box, p_gf_box, p_dr_box, 
  p_nga50_box, p_mism_box,
  ncol = 3, labels = "AUTO"
)

assembly_evaluation_plot6 <- plot_grid(p6, legend, rel_widths = c(3, 0.5))

ggsave(filename = output_figure, plot = assembly_evaluation_plot6, width = 11, height = 8)

metaquast_all_criteria %>% 
    group_by(assembler, criteria) %>% 
    summarize(average = mean(value, na.rm=T), sd = sd(value, na.rm=T)) %>% 
    group_by(criteria) %>% 
    mutate(rank = order(order(average, decreasing=ifelse(criteria %in% c(
    "Genome_fraction",
    "Largest_alignment", 
    "NGA50"), TRUE, FALSE)))) -> summary.rank.ordered

write.table(summary.rank.ordered, 
  file = output_table, 
  sep = "\t",
  quote=F,
  row.names=F)

summary.rank.for.score <- filter(summary.rank.ordered, criteria != "num_misassemblies")

# weights for 
# Duplication_ratio, 
# Genome_fraction
# Largest_alignment
# NGA50
# num_contigs
# num_mismatches_per_100_kbp
weights <- c(0.1, 0.25, 0.25, 0.2, 0.1, 0.1)

rank.max <- unique(summary.rank.for.score$assembler) %>% 
  length()
rank.min <- 1

score.max <- 10
score.min <- 1

mutate(summary.rank.for.score, score = (rank.max-rank)  * ((score.max-score.min)/(rank.max-rank.min)) + 1) %>%
    group_by(assembler) %>% 
    summarize(weighted.score=sum(score*weights)) ->
    summary.weighted.score


write.table(summary.weighted.score, 
  file = output_score, 
  sep = "\t",
  quote=F,
  row.names=F)