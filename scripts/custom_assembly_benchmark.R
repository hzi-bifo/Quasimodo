library(tidyverse)
library(cowplot)

## Load the parameters from snakemake
output_figure <- snakemake@output$assembly_benchmark_figure
#output_table <- snakemake@output$assembly_rank_table
criteria_files <- snakemake@input
## The visualized criteria
criteria <- c(
  "num_contigs", "Largest_contig", "Genome_fraction",
  "Duplication_ratio", "Largest_alignment", "LGA50",
  "NGA50", "num_misassemblies", "num_mismatches_per_100_kbp"
)

kb_scale_criteria <- c("Largest_alignment", "Largest_contig", "NGA50") 

## Read the metaquast summary for assemblers
read_metaquast <- function(metaquast_file) {
  filename <- sub("^([^.]*).*", "\\1", basename(metaquast_file)) 
  metaquast <- as.data.frame(t(read.table(metaquast_file, 
        sep = "\t", 
        header = F, 
        row.names = 1,
        na.strings = "-", 
        check.names = F)))
  
  metaquast$not_aligned <- NULL
  
  metaquast_long <- gather(metaquast, reference, value, -Assemblies)
  metaquast_long$value <- as.numeric(metaquast_long$value)
  if (filename %in% kb_scale_criteria){
      metaquast_long$value <- metaquast_long$value / 1000
  }
  return(metaquast_long)
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
  p <- theme(axis.text.x = element_blank(), 
  legend.position = "none")
  return(p)
}


color_jitter_boxplot <- function(df, label) {
  g <- ggplot(df, aes(Assemblies, value)) +
    geom_boxplot(show.legend = F, fill = "grey85") +
    geom_point(
      position = position_jitterdodge(
        jitter.width = 0,
        dodge.width = 0.3,
        seed = 1234
      ),
      aes(color = reference, shape=Assemblies),
      stroke = 0,
      size = 2
    ) + xlab("") + 
    theme_bw(base_size = 15) +
    ylab(label) + 
    scale_color_brewer(palette = "Set1")
  return(g)
}

color_jitterplot <- function(df, label) {
  g <- ggplot(df, aes(Assemblies, value, color=reference, shape=Assemblies)) +
    geom_jitter(stat = "identity", size=3) +
    xlab("") + 
    theme_bw(base_size = 15) +
    ylab(label) + 
    scale_color_brewer(palette = "Set1")
  return(g)
}

metaquast_criteria_list <- list()
metaquast_criteria_list[["table"]] <- list()
metaquast_criteria_list[["figure"]] <- list()

gg_func <- ifelse(snakemake@params$number_ref > 5, color_jitter_boxplot, color_jitterplot)


for (criterium_file in criteria_files){
    criterium_name <- sub("^([^.]*).*", "\\1", basename(criterium_file))
    criterium_table <- read_metaquast(criterium_file)
    metaquast_criteria_list[["table"]][[criterium_name]] <- criterium_table
    ylab <- gsub("_", " ", criterium_name)
    ylab <- ifelse(criterium_name %in% kb_scale_criteria, paste(ylab, "(Kb)", sep=" "), ylab)
    metaquast_criteria_list[["figure"]][[criterium_name]]  <- gg_func(criterium_table, ylab) + theme_wo_legend()
}

gg_legend <- get_legend(gg_func(metaquast_criteria_list[["table"]][[1]], "") + theme_with_legend())

## Combine all visualization and legend
p <- plot_grid(plotlist = metaquast_criteria_list[["figure"]],
  ncol = 3, labels = "AUTO"
)

assembly_evaluation_plot <- plot_grid(p, gg_legend, rel_widths = c(3, 0.8))

ggsave(filename = output_figure, plot = assembly_evaluation_plot, width = 12, height = 10)

