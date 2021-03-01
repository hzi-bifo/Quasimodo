library(tidyverse)
library(VennDiagram)
library(gridExtra)

## VCF files for mixtures
snp_vcf <- snakemake@input$snp

## Sample list
samples <- unlist(snakemake@params$mix_sample)

## SNP callers to compare using Venndigram
venn_snpcallers <- unlist(snakemake@params$venn_snpcallers)

## The genome difference by NUCmer
genome_diff_files <- snakemake@input$diff

## Make a list for the genome difference

rtg_summary <- snakemake@input$rtg_summary
rtg_snp_summary <- snakemake@input$rtg_snp_summary
rtg_roc <- snakemake@input$rtg_roc
rtg_snp_roc <- snakemake@input$rtg_snp_roc

caller_map <- list(
    bcftools = "BCFtools", clc = "CLC", freebayes = "FreeBayes",
    gatk = "GATK", lofreq = "LoFreq", varscan = "VarScan2"
)

make_snp_vector <- function(vcf) {
    snp_table <- tryCatch(
        {
            read.table(vcf,
                sep = "\t", comment.char = "#", header = F,
                colClasses = c("NULL", "character", "NULL", rep("character", 2), rep("NULL", 3))
            )
        },
        error = function(e) {
            snp_table <- data.frame(position = character(), ref = character(), alt = character())
        }
    )

    if (nrow(snp_table) > 0) {
        snp_table <- snp_table[, 1:3]
        colnames(snp_table) <- c("position", "ref", "alt")
        snp_vector <- do.call(paste, c(snp_table %>% filter(
            ref %in% c("A", "C", "G", "T"),
            alt %in% c("A", "C", "G", "T")
        ), sep = "-"))
    }
    else {
        snp_vector <- c()
    }

    snp_vector
}

genome_diff_list <- list()
for (genome_diff_file in genome_diff_files) {
    genome_diff_vector <- make_snp_vector(genome_diff_file)
    mixstrain <- gsub("\\..*$", "", basename(genome_diff_file))
    genome_diff_list[[mixstrain]] <- genome_diff_vector
}


## Empty dataframe for the performance of SNP callers
snpcaller_performance_summary <- data.frame(
    caller = character(), mixture = character(),
    genomediff = integer(), calleridentify = integer(),
    TP = integer(), FP = integer(), Precision = double(),
    Recall = double(), F1 = double()
)

snpcaller_performance_names <- colnames(snpcaller_performance_summary)
snp_list <- list()


for (file in snp_vcf) {
    snpcaller_ref_sample <- unlist(strsplit(basename(file), "\\."))
    sample <- snpcaller_ref_sample[1]
    ref_name <- snpcaller_ref_sample[2]
    caller_name_lower <- snpcaller_ref_sample[3]
    snpcaller <- caller_map[[caller_name_lower]]

    snp_vector <- make_snp_vector(file)
    snp_identify_count <- length(snp_vector)

    if (sample %in% samples) {
        mixstrain <- substr(sample, 1, 2)
        genome_diff_vector <- genome_diff_list[[mixstrain]]
        genome_diff_count <- length(genome_diff_vector)


        if (snp_identify_count > 0) {
            tp_count <- length(intersect(snp_vector, genome_diff_vector))
            fp_count <- length(setdiff(snp_vector, genome_diff_vector))
            fn_count <- length(setdiff(genome_diff_vector, snp_vector))
            precision <- round(tp_count / snp_identify_count, 3)
            recall <- round(tp_count / genome_diff_count, 3)
            f1 <- round(2 * (precision * recall) / (precision + recall), 3)
        }
        else {
            tp_count <- 0
            fp_count <- 0
            fn_count <- 0
            precision <- NA
            recall <- NA
            f1 <- NA
        }

        if (caller_name_lower %in% venn_snpcallers) {
            snpcaller <- caller_map[[caller_name_lower]]
            if (sample %in% names(snp_list)) {
                snp_list[[sample]][[snpcaller]] <- snp_vector
            }
            else {
                snp_list[[sample]] <- list(genome_diff_vector, snp_vector)
                names(snp_list[[sample]]) <- c("Genome", snpcaller)
            }
        }
    }
    else {
        genome_diff_count <- 0
        tp_count <- 0
        fp_count <- snp_identify_count
        precision <- 0
        recall <- NA
        f1 <- NA
    }

    snpcaller_performance <- data.frame(
        snpcaller, sample, genome_diff_count, snp_identify_count,
        tp_count, fp_count, precision, recall, f1
    )
    names(snpcaller_performance) <- snpcaller_performance_names
    snpcaller_performance_summary <- rbind(snpcaller_performance_summary, snpcaller_performance)
}


## Write the table for the performance summary of SNP callers
write.table(snpcaller_performance_summary,
    file = snakemake@output$performance_table,
    sep = "\t", row.names = F, quote = F
)


### RTG plots

make_rtg_df <- function(summary_files) {
    rtg_df_list <- list()

    header <- c(
        "threshold", "TP_baseline", "TP_call",
        "FP", "FN", "Precision", "Recall", "F1"
    )

    for (rtg in summary_files) {
        mixture <- unlist(strsplit(basename(dirname(rtg)), "\\."))[1]
        caller <- caller_map[[basename(dirname(dirname(rtg)))]]

        df <- read.table(rtg, sep = "", strip.white = T, stringsAsFactors = F, skip = 2)
        colnames(df) <- header
        df$mixture <- mixture
        df$caller <- caller
        rtg_df_list[[paste(mixture, caller, sep = ".")]] <- df [1, ]
    }

    rtg_df <- as.data.frame(do.call(rbind, rtg_df_list), stringsAsFactors = F)[, c(10, 9, 1:8)]

    callers <- unique(as.character(rtg_df$caller))
    mixtures <- unique(as.character(rtg_df$mixture))

    rtg_df$caller <- factor(rtg_df$caller, levels = callers)
    rtg_df$mixture <- factor(rtg_df$mixture, levels = sort(mixtures))
    return(rtg_df)
}


rtg_df <- make_rtg_df(rtg_summary)

rtg_snp_df <- make_rtg_df(rtg_snp_summary)

write.table(rtg_df,
    file = snakemake@output$rtg_performance_table,
    sep = "\t", row.names = F, quote = F
)

write.table(rtg_snp_df,
    file = snakemake@output$rtg_snp_performance_table,
    sep = "\t", row.names = F, quote = F
)

# rtg_df$caller=factor(rtg_df$caller, levels=callers)


point_plot <- function(df) {
    point_plot_whole <- ggplot(df, aes(Precision, Recall, color = caller, shape = mixture)) +
        geom_point(size = 4) +
        theme(
            legend.position = "right",
            legend.spacing.x = unit(0.3, "cm"),
            text = element_text(size = 15)
        ) +
        theme_bw(base_size = 15) +
        xlim(0, 1) +
        ylim(0, 1)
    #  guides(colour = guide_legend(override.aes = list(size=4)), shape = guide_legend(override.aes = list(size=4))) +


    print("Point plot inset zoom")
    inset_zoom <- ggplot(df, aes(Precision, Recall, color = caller, shape = mixture)) +
        geom_point(size = 4) +
        xlim(0.8, 1) +
        ylim(0.5, 1) +
        theme_bw(base_size = 16) +
        theme(legend.position = "none") +
        theme(text = element_text(size = 15)) +
        scale_colour_brewer(palette = "Set1") +
        xlab("") +
        ylab("")

    # Scatter plot
    print("Combine with inset plot")
    rect <- data.frame(x1 = 0.8, x2 = 1, y1 = 0.5, y2 = 1)
    point_plot <- point_plot_whole +
        annotation_custom(ggplotGrob(inset_zoom), xmin = -0.05, xmax = 0.7, ymin = 0.06, ymax = 0.92) +
        geom_rect(
            data = rect, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            alpha = 0.1, inherit.aes = F, color = "#A3A3A3"
        ) +
        geom_segment(aes(x = 0.678, y = 0.208, xend = 1, yend = 0.5), color = "#A3A3A3") +
        geom_segment(aes(x = 0.086, y = 0.892, xend = 0.8, yend = 1), color = "#A3A3A3") +
        scale_colour_brewer(palette = "Set1")
    return(point_plot)
}


# RTG Pointplot
rtg_point_plot <- point_plot(rtg_df)
# RTG Boxplot
print("Box plot")
rtg_box_plot <- ggplot(rtg_df, aes(caller, F1, fill = caller)) +
    geom_boxplot() +
    geom_jitter() +
    theme_bw(base_size = 16) +
    theme(
        legend.position = "none",
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 60, hjust = 1)
    ) +
    scale_fill_brewer(palette = "Set1")


# ggsave(file='rtg_threshold_precision_recall_f1.png', plot=rtg_threshold_recall_precision_f1_plot, width=12, height=6)

### SNP plots

### venndigram

callers <- unique(as.character(rtg_df$caller))
mixtures <- unique(as.character(rtg_df$mixture))

snpcaller_performance_summary$caller <- as.factor(snpcaller_performance_summary$caller)


mixture_snpcaller_performance_summary <- dplyr::filter(snpcaller_performance_summary, !grepl("-0", mixture))

mixture_snpcaller_performance_summary$mixture <- factor(mixture_snpcaller_performance_summary$mixture,
    levels = sort(mixtures)
)


# SNP score >=20 Pointplot

nortg_point_plot <- point_plot(mixture_snpcaller_performance_summary)

# Boxplot
# print("Box plot")
nortg_box_plot <- ggplot(mixture_snpcaller_performance_summary, aes(caller, F1, fill = caller)) +
    geom_boxplot() +
    geom_jitter() +
    theme_bw(base_size = 16) +
    theme(
        legend.position = "none",
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 60, hjust = 1)
    ) +
    scale_fill_brewer(palette = "Set1")

combined_plots <- list()

combined_plots[[1]] <- rtg_point_plot
combined_plots[[2]] <- rtg_box_plot
combined_plots[[3]] <- nortg_point_plot
combined_plots[[4]] <- nortg_box_plot


## Arrange the Venndigram
lay <- rbind(
    c(1, 1, 1, 1, 2, 2),
    c(1, 1, 1, 1, 2, 2),
    c(1, 1, 1, 1, 2, 2),
    c(3, 3, 3, 3, 4, 4),
    c(3, 3, 3, 3, 4, 4),
    c(3, 3, 3, 3, 4, 4)
)


print("Grid all point and box plots")
p <- grid.arrange(grobs = combined_plots, layout_matrix = lay)

ggsave(file = snakemake@output$performance_figure, plot = p, width = 12, height = 11) # ,


rtg_snp_point_plot <- point_plot(rtg_snp_df)
rtg_snp_box_plot <- ggplot(rtg_snp_df, aes(caller, F1, fill = caller)) +
    geom_boxplot() +
    geom_jitter() +
    theme_bw(base_size = 16) +
    theme(
        legend.position = "none",
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 60, hjust = 1)
    ) +
    scale_fill_brewer(palette = "Set1")

## Read the roc.tsv to extract precision, recall, F1 with score >= 20


make_rtg_roc_df <- function(rtg_roc_files) {
    roc_df_list <- list()

    for (roc_file in rtg_roc_files) {
        caller <- caller_map[[gsub("\\/.*$", "", gsub(".*rtg\\/", "", roc_file))]]
        mixture <- gsub("\\..*$", "", basename(gsub("\\/weighted_roc\\.tsv\\.gz", "", roc_file)))
        header <- c(
            "score", "TP_baseline", "FP",
            "TP_call", "FN", "Precision", "Recall", "F1"
        )

        roc_df <- read_tsv(roc_file, comment = "#", col_names = header)
        roc_df$mixture <- mixture
        roc_df$caller <- caller

        roc_df_list[[paste(caller, mixture, sep = ".")]] <- roc_df
    }
    roc_df_combined <- as.data.frame(do.call(rbind, roc_df_list))[, c(10, 9, 1:8)]
    roc_df_combined$caller <- factor(roc_df_combined$caller, levels = callers)
    roc_df_combined$mixture <- factor(roc_df_combined$mixture, levels = sort(mixtures))
    return(roc_df_combined)
}

# roc_df_list <- list()

# for (roc_file in rtg_roc) {
#     caller <- caller_map[[gsub("\\/.*$", "", gsub(".*rtg\\/", "", roc_file))]]
#     mixture <- tools::file_path_sans_ext(basename(gsub("\\/weighted_roc\\.tsv\\.gz", "", roc_file)))
#     header <- c(
#         "score", "TP_baseline", "FP",
#         "TP_call", "FN", "Precision", "Recall", "F1"
#     )

#     roc_df <- read_tsv(roc_file, comment = "#", col_names = header)
#     roc_df$mixture <- mixture
#     roc_df$caller <- caller

#     roc_df_list[[paste(caller, mixture, sep = ".")]] <- roc_df
# }

# roc_df_combined <- as.data.frame(do.call(rbind, roc_df_list))[, c(10, 9, 1:8)]

rtg_roc_df <- make_rtg_roc_df(rtg_roc)

write.table(rtg_roc_df,
    file = snakemake@output$rtg_roc_table,
    sep = "\t", row.names = F, quote = F
)

make_rtg_roc_df(rtg_roc) %>%
    dplyr::filter(score >= 20) %>%
    group_by(caller, mixture) %>%
    slice(which.min(score)) -> roc_score20_df


write.table(roc_score20_df,
    file = snakemake@output$rtg_score20_performance_table,
    sep = "\t", row.names = F, quote = F
)


rtg_score20_point_plot <- point_plot(roc_score20_df)
rtg_score20_box_plot <- ggplot(roc_score20_df, aes(caller, F1, fill = caller)) +
    geom_boxplot() +
    geom_jitter() +
    theme_bw(base_size = 16) +
    theme(
        legend.position = "none",
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 60, hjust = 1)
    ) +
    scale_fill_brewer(palette = "Set1")

rtg_snp_roc_df <- make_rtg_roc_df(rtg_snp_roc)

write.table(rtg_snp_roc_df,
    file = snakemake@output$rtg_snp_roc_table,
    sep = "\t", row.names = F, quote = F
)


rtg_snp_roc_df %>%
    dplyr::filter(score >= 20) %>%
    group_by(caller, mixture) %>%
    slice(which.min(score)) -> roc_score20_snp_df

write.table(roc_score20_snp_df,
    file = snakemake@output$rtg_snp_score20_performance_table,
    sep = "\t", row.names = F, quote = F
)

rtg_snp_score20_point_plot <- point_plot(roc_score20_snp_df)
rtg_snp_score20_box_plot <- ggplot(
    roc_score20_snp_df,
    aes(caller, F1, fill = caller)
) +
    geom_boxplot() +
    geom_jitter() +
    theme_bw(base_size = 16) +
    theme(
        legend.position = "none",
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 60, hjust = 1)
    ) +
    scale_fill_brewer(palette = "Set1")


rtg_combined_plots <- list()

rtg_combined_plots[[1]] <- rtg_snp_point_plot
rtg_combined_plots[[2]] <- rtg_snp_box_plot
rtg_combined_plots[[3]] <- rtg_score20_point_plot
rtg_combined_plots[[4]] <- rtg_score20_box_plot
rtg_combined_plots[[5]] <- rtg_snp_score20_point_plot
rtg_combined_plots[[6]] <- rtg_snp_score20_box_plot

## Arrange the Venndigram
lay2 <- rbind(
    c(1, 1, 1, 1, 2, 2),
    c(1, 1, 1, 1, 2, 2),
    c(1, 1, 1, 1, 2, 2),
    c(3, 3, 3, 3, 4, 4),
    c(3, 3, 3, 3, 4, 4),
    c(3, 3, 3, 3, 4, 4),
    c(5, 5, 5, 5, 6, 6),
    c(5, 5, 5, 5, 6, 6),
    c(5, 5, 5, 5, 6, 6)
)

print("Grid all point and box plots")
rtg_p <- grid.arrange(grobs = rtg_combined_plots, layout_matrix = lay2)

ggsave(
    file = snakemake@output$rtg_performance_figure,
    plot = rtg_p,
    width = 12,
    height = 16
)


## ROC curves

rtg_roc_df$mix <- gsub("-.*$", "", rtg_roc_df$mixture)
rtg_roc_df$ratio <- gsub("^[^-]+-", "", rtg_roc_df$mixture)

roc_p <- ggplot(rtg_roc_df, aes(Recall, Precision, color = caller)) +
    geom_point(size = 1.5, alpha = 0.4, stroke = 0) +
    theme_bw(base_size = 20) +
    facet_grid(ratio ~ mix) +
    scale_color_brewer(palette = "Set1") +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))


ggsave(
    file = snakemake@output$rtg_roc_figure,
    plot = roc_p,
    width = 12,
    height = 14
)



rtg_snp_roc_df$mix <- gsub("-.*$", "", rtg_snp_roc_df$mixture)
rtg_snp_roc_df$ratio <- gsub("^[^-]+-", "", rtg_snp_roc_df$mixture)

snp_roc_p <- ggplot(rtg_snp_roc_df, aes(Recall, Precision, color = caller)) +
    geom_point(size = 1.5, alpha = 0.4, stroke = 0) +
    theme_bw(base_size = 20) +
    facet_grid(ratio ~ mix) +
    scale_color_brewer(palette = "Set1") +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))


ggsave(
    file = snakemake@output$rtg_snp_roc_figure,
    plot = snp_roc_p,
    width = 12,
    height = 14
)


## venndigrams
combined_venndiagrams <- list()

n_panel <- 0

for (mix in samples) {
    n_panel <- n_panel + 1
    futile.logger::flog.threshold(futile.logger::ERROR,
        name = "VennDiagramLogger"
    )
    # Venndiagram
    plot <- venn.diagram(
        x = snp_list[[mix]],
        filename = NULL,
        main = mix,
        col = "transparent",
        fill = c("blueviolet", "#dc143c", "dodgerblue1", "green3"),
        alpha = 0.60,
        fontfamily = "sans",
        cat.fontfamily = "sans",
        main.fontfamily = "sans",
        cex = 1.2,
        main.cex = 1.3,
        cat.cex = 1.2
    )
    combined_venndiagrams[[n_panel]] <- grobTree(plot)
}

## Arrange the Venndigram
venn_lay <- rbind(
    c(1, 2, 3),
    c(4, 5, 6)
)


print("Grid all venn plots")
p2 <- grid.arrange(grobs = combined_venndiagrams, layout_matrix = venn_lay)

ggsave(
    file = snakemake@output$snp_venndiagram_figure,
    plot = p2,
    width = 12,
    height = 8
)