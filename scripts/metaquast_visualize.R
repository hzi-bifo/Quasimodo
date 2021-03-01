library(tidyverse)
library(cowplot)

## Load the parameters from snakemake
input_dir <- snakemake@params$input_dir
combined_ref_reports <- snakemake@input$combined_ref_reports

output_figure <- snakemake@output$figure
output_table <- snakemake@output$table
output_scaled <- snakemake@output$radarplot_table
output_score <- snakemake@output$score

assembler_map <- list(
  spades = "SPAdes", metaspa = "metaSPA",
  metaspades = "metaSPAdes",
  tadpole = "tadpole", abyss = "ABySS",
  megahit = "megahit", ray = "Ray",
  idba = "IDBA", vicuna = "Vicuna",
  iva = "IVA", savage = "Savage",
  virgena = "VirGenA"
)


read_metaquast <- function(file_name, dir = input_dir) {
  file <- paste(dir, file_name, sep = "/")

  metaquast <- read.table(file, sep = "\t", header = T, na.strings = "-", check.names = F) %>%
    separate(Assemblies, c("sample", "assembler"), sep = "\\.")
  ref1 <- "TB40E"
  ref2 <- ifelse(startsWith(basename(file_name), "TA"), "AD169", "Merlin")

  # ref1 <- 3
  # ref2 <- 4
  if (grepl("NGA50", file_name)) {
    metaquast[, ref1] <- ifelse(endsWith(metaquast$sample, "_0_1"), -1, metaquast[, ref1])
    metaquast[, ref2] <- ifelse(endsWith(metaquast$sample, "_1_0"), -1, metaquast[, ref2])
    metaquast[is.na(metaquast)] <- 0
    metaquast[metaquast == -1] <- NA
  }
  else {
    metaquast[, ref1] <- ifelse(endsWith(metaquast$sample, "_0_1"), NA, metaquast[, ref1])
    metaquast[, ref2] <- ifelse(endsWith(metaquast$sample, "_1_0"), NA, metaquast[, ref2])
  }
  metaquast$not_aligned <- NULL
  metaquast$assembler <- ifelse(metaquast$assembler == "metaspades", "metaspa", metaquast$assembler)
  return(metaquast)
}

metaquast_criteria <- c(
  "num_contigs", "Genome_fraction",
  "Duplication_ratio", "Largest_alignment",
  "NGA50", "num_mismatches_per_100_kbp"
)

files <- list.files(input_dir, pattern = "\\.merged\\.tsv$")
metaquast_list <- list()
for (file in files) {
  mix_criteria <- gsub("\\.merged\\.tsv$", "", file)
  criteria <- unlist(strsplit(mix_criteria, "\\."))[2]
  if (criteria %in% metaquast_criteria) {
    df <- read_metaquast(file)
    df$criteria <- criteria
    df.long <- gather(df, reference, value, -c(sample, assembler, criteria))
    metaquast_list[[mix_criteria]] <- df.long
  }
}

metaquast_all_criteria <- as.data.frame(do.call(rbind, metaquast_list))
# %>%
#   mutate(assembler=assembler_map[[assembler]]) ->
# print(unique(metaquast_all_criteria$assembler))

metaquast_all_criteria$assembler <- unlist(assembler_map[metaquast_all_criteria$assembler],
  use.names = FALSE
)

write.table(metaquast_all_criteria,
  file = snakemake@output$all_sample_table,
  sep = "\t",
  quote = F,
  row.names = F
)


## Nicer theme with and without legend
theme_with_legend <- function() {
  p <- guides(color = guide_legend(
    override.aes = list(size = 5),
    label.theme = element_text(size = 13)
  ))
  return(p)
}

theme_wo_legend <- function() {
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
  portion <- ifelse(portion == "1/0", "pure", portion)
  sample_summary$ratio <- portion
  sample_summary <- filter(
    sample_summary,
    ratio != "0",
    ratio != "0/1"
  )
}


metaquast_all_criteria_portion <- metaquast_all_criteria %>%
  get_portion()


color_dot_boxplot <- function(df, value, label) {
  g <- ggplot(df, aes_string("assembler", value)) +
    geom_boxplot(show.legend = F, fill = "grey85", outlier.shape = NA) +
    geom_point(
      position = position_jitterdodge(
        jitter.width = 0,
        dodge.width = 0.3,
        seed = 1234
      ),
      na.rm = TRUE,
      aes(color = ratio),
      stroke = 0,
      size = 2
    ) +
    xlab("") +
    theme_bw(base_size = 15) +
    ylab(label) +
    scale_color_brewer(palette = "Set2")
  return(g)
}

## Number of contigs visualization
p_c_box <- color_dot_boxplot(
  filter(
    metaquast_all_criteria_portion,
    criteria == "num_contigs"
  ),
  "value",
  "# Contigs"
) +
  scale_y_log10() +
  theme_wo_legend()

## Only legend for all figures
p_c_box_legend <- color_dot_boxplot(
  filter(
    metaquast_all_criteria_portion,
    criteria == "num_contigs"
  ),
  "value",
  "# Contigs"
) +
  theme_with_legend()


## Genome fraction visualization
p_gf_box <- color_dot_boxplot(
  filter(
    metaquast_all_criteria_portion,
    criteria == "Genome_fraction"
  ),
  "value",
  "% Genome fraction"
) +
  # scale_y_log10() +
  theme_wo_legend()



## Duplication ratio visualization
p_dr_box <- color_dot_boxplot(
  filter(
    metaquast_all_criteria_portion,
    criteria == "Duplication_ratio"
  ),
  "value",
  "Duplication ratio"
) +
  theme_wo_legend()


## Largest alignment visualization
p_la_box <- color_dot_boxplot(
  filter(
    metaquast_all_criteria_portion,
    criteria == "Largest_alignment"
  ),
  "value / 1000",
  "Largest Alignment (kb)"
) +
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
p_nga50_box <- color_dot_boxplot(
  filter(
    metaquast_all_criteria_portion,
    criteria == "NGA50"
  ),
  "value / 1000",
  "NGA50 (kb)"
) +
  theme_wo_legend()



## Number of misassemblies visualization
# p_mis_box <- color_dot_boxplot(rbind(TM.misassemblies.long,
#                         TA.misassemblies.long),
#                   "misassemblies",
#                   "# Misassemblies") +
#   theme_wo_legend()


## Number of mismatches per 100kbp visualization
p_mism_box <- color_dot_boxplot(
  filter(
    metaquast_all_criteria_portion,
    criteria == "num_mismatches_per_100_kbp"
  ),
  "value",
  "# Mismatches/100kbp"
) +
  scale_y_log10() +
  theme_wo_legend()


## Get the legend
legend <- get_legend(p_c_box_legend)


p6 <- plot_grid(p_c_box, p_la_box, p_gf_box, p_dr_box,
  p_nga50_box, p_mism_box,
  ncol = 3, labels = "AUTO"
)

assembly_evaluation_plot6 <- plot_grid(p6, legend, rel_widths = c(3, 0.5))

ggsave(filename = output_figure, plot = assembly_evaluation_plot6, width = 14, height = 8)

filter(metaquast_all_criteria, criteria == "NGA50") %>%
  select(criteria, sample, assembler, value) %>%
  group_by(criteria, sample, assembler) %>%
  summarize(value = mean(value, na.rm = T)) %>%
  as.data.frame() ->
av_individual_NGA50

# print(unique(av_individual_NGA50$assembler))

av_individual_NGA50[av_individual_NGA50$assembler == "metaSPA", "assembler"] <- "metaSPAdes"

big_better <- c("Genome fraction (%)", "Largest alignment")
small_better <- c("Duplication ratio", "# contigs", "# mismatches per 100 kbp")
combined_ref_criteria_names <- c(big_better, small_better)

print("Loading combined ref reports...")
combined_ref_report_list <- list()
for (report in combined_ref_reports) {
  df <- read.delim2(report, sep = "\t", header = T, row.names = 1, na.strings = "-", check.names = F, stringsAsFactors = F)
  sample <- gsub("\\..*$", "", colnames(df)[1])
  colnames(df) <- gsub("^[^\\.]+\\.|\\.scaffolds", "", colnames(df))
  # print(as.numeric(as.character(df[df$criteria == "# contigs",])))
  # print(as.numeric(gsub(" + .*", "", as.character(df[df$criteria == "# unaligned contigs",]))))
  num_aligned_contigs <- as.numeric(unname(df["# contigs", ])) - as.numeric(gsub(" \\+ .*", "", unname(df["# unaligned contigs", ])))

  if (endsWith(sample, "_0_1") || endsWith(sample, "_1_0")) {
    df["# contigs", ] <- num_aligned_contigs
  }
  else {
    df["# contigs", ] <- num_aligned_contigs / 2
  }

  df$criteria <- row.names(df)
  combined_ref_report_list[[sample]] <- df %>%
    filter(criteria %in% combined_ref_criteria_names) %>%
    mutate(sample = sample) %>%
    gather(assembler, value, -c(criteria, sample))
}

combined_ref_report_df <- as.data.frame(do.call(rbind, combined_ref_report_list))
combined_ref_report_df$value <- as.numeric(combined_ref_report_df$value)
combined_ref_report_df$criteria <- as.character(combined_ref_report_df$criteria)

# print(unique(combined_ref_report_df$assembler))

combined_ref_report_df$assembler <- unlist(assembler_map[combined_ref_report_df$assembler],
  use.names = FALSE
)
# print(combined_ref_report_df)

big_better <- c(big_better, "NGA50")

summary.rank.ordered <- as.data.frame(rbind(combined_ref_report_df, av_individual_NGA50)) %>%
  group_by(assembler, criteria) %>%
  summarize(average = mean(value, na.rm = T), sd = sd(value, na.rm = T)) %>%
  group_by(criteria) %>%
  mutate(rank = order(order(average, decreasing = ifelse(criteria %in% big_better, TRUE, FALSE))))

write.table(summary.rank.ordered,
  file = output_table,
  sep = "\t",
  quote = F,
  row.names = F
)

assembly_big_better <- summary.rank.ordered %>%
  filter(criteria %in% big_better) %>%
  group_by(criteria) %>%
  mutate(min = min(average, na.rm = T), max = (max(average, na.rm = T))) %>%
  mutate(scaled = (average - min) / (max - min))


assembly_small_better <- summary.rank.ordered %>%
  filter(criteria %in% small_better) %>%
  group_by(criteria) %>%
  mutate(min = min(average, na.rm = T), max = (max(average, na.rm = T))) %>%
  mutate(scaled = (max - average) / (max - min))


rbind(assembly_big_better, assembly_small_better) %>%
  select(assembler, criteria, scaled) %>%
  spread(criteria, scaled) ->
assembly.scaled.value


write.table(assembly.scaled.value,
  file = output_scaled,
  sep = "\t",
  quote = F,
  row.names = F
)

weights <- c(0.3, 0.3, 0.1, 0.1, 0.1, 0.1)

rbind(assembly_big_better, assembly_small_better) %>%
  select(assembler, criteria, scaled) %>%
  group_by(assembler) %>%
  summarize(weighted.score = 10 * sum(scaled * weights)) ->
summary.weighted.score

write.table(summary.weighted.score,
  file = output_score,
  sep = "\t",
  quote = F,
  row.names = F
)

## Calculate score and rank based on individual reference mapping in metaQUAST
# metaquast_all_criteria %>%
#     group_by(assembler, criteria) %>%
#     summarize(average = mean(value, na.rm=T), sd = sd(value, na.rm=T)) %>%
#     group_by(criteria) %>%
#     mutate(rank = order(order(average, decreasing=ifelse(criteria %in% c(
#     "Genome_fraction",
#     "Largest_alignment",
#     "NGA50"), TRUE, FALSE)))) -> summary.rank.ordered

# write.table(summary.rank.ordered,
#   file = output_table,
#   sep = "\t",
#   quote=F,
#   row.names=F)

# #summary.rank.for.score <- filter(summary.rank.ordered, criteria != "num_misassemblies")

# # weights for

# # Genome_fraction
# # Largest_alignment
# # NGA50
# # Duplication_ratio
# # num_contigs
# # num_mismatches_per_100_kbp
# #weights <- c(0.1, 0.25, 0.25, 0.2, 0.1, 0.1)
# #weights <- c(0.1, 0.3, 0.2, 0.1, 0.1, 0.2)
# #weights <- c(0.1, 0.3, 0.3, 0.1, 0.1, 0.1)
# weights <- c(0.3, 0.3, 0.1, 0.1, 0.1, 0.1)


# # rank.max <- unique(summary.rank.for.score$assembler) %>%
# #   length()
# # rank.min <- 1

# # score.max <- 10
# # score.min <- 1

# # mutate(summary.rank.for.score, score = (rank.max-rank)  * ((score.max-score.min)/(rank.max-rank.min)) + 1) %>%
# #     group_by(assembler) %>%
# #     summarize(weighted.score=sum(score*weights)) ->
# #     summary.weighted.score


# # write.table(summary.weighted.score,
# #   file = output_score,
# #   sep = "\t",
# #   quote=F,
# #   row.names=F)



# # Score based on average performance

# big_better <- c("Genome_fraction", "Largest_alignment", "NGA50")
# small_better <- c("Duplication_ratio", "num_contigs", "num_mismatches_per_100_kbp")


# summary.rank.ordered %>%
#   filter(criteria %in% big_better) %>%
#   group_by(criteria) %>%
#   mutate(min = min(average, na.rm=T), max = (max(average, na.rm=T))) %>%
#   mutate(scaled = (average - min)/(max-min)) ->
#   assembly_big_better

# summary.rank.ordered %>%
#   filter(criteria %in% small_better) %>%
#   group_by(criteria) %>%
#   mutate(min = min(average, na.rm=T), max = (max(average, na.rm=T))) %>%
#   mutate(scaled = (max - average)/(max-min)) ->
#   assembly_small_better

# rbind(assembly_big_better, assembly_small_better) %>%
#   select(assembler, criteria, scaled) %>%
#   group_by(assembler) %>%
#   summarize(weighted.score=10*sum(scaled*weights)) ->
#   summary.weighted.score

# write.table(summary.weighted.score,
#   file = output_score,
#   sep = "\t",
#   quote=F,
#   row.names=F)

# # if (!require(ggradar)) {
# #   devtools::install_github("ricardo-bion/ggradar",
# #                           dependencies = TRUE)
# # }
# # library(ggradar)

# # rbind(assembly_big_better, assembly_small_better) %>%
# #   select(assembler, criteria, scaled) %>%
# #   spread(criteria, scaled) %>%
# #   rename("Genome fraction" = "Genome_fraction",
# #          "Largest alignment" = "Largest_alignment",
# #          "Duplication ratio" = "Duplication_ratio",
# #          "number contigs" = "num_contigs",
# #          "mismatches"="num_mismatches_per_100_kbp") ->
# #   assembly_radar

# # radar_av_metric <- ggradar(assembly_radar)

# # ggsave(filename = output_radar, plot = radar_av_metric, width = 10, height = 8)