library(tidyverse)
library(VennDiagram)
library(gridExtra)

## Two mixtures
mix_vcf <- snakemake@input$mix
#all_caller_table <- data.frame(caller=character(), mixture=character(), 
#                               genome=character(), position=character(), 
#                               ref=character(), alt=character())


#samples <- gsub("^.*\\.", "", snakemake@params$mix_sample)
#samples <- gsub("\\..*$", "", snakemake@params$mix_sample)

## Sample list
samples <- unlist(snakemake@params$mix_sample)

## SNP callers to compare using Venndigram
venn_snpcallers <- unlist(snakemake@params$venn_snpcallers)

## The genome difference by NUCmer
genome_diff_files <- snakemake@input$diff

## Make a list for the genome difference
genome_diff_list <- list()
for (genome_diff_file in genome_diff_files){
    genome_diff_table <- read.table(genome_diff_file, sep="\t", comment.char="#", header=F,
                                    colClasses=c(rep("character", 3), rep("NULL", 9)))
    colnames(genome_diff_table) <- c("position", "ref", "alt")
    genome_diff_table_filtered <- genome_diff_table %>% filter(ref!=".", alt!=".")
    genome_diff_vector <- do.call(paste, c(genome_diff_table_filtered, sep="-"))
    mixstrain <- gsub("\\..*$", "", basename(genome_diff_file))
    genome_diff_list[[mixstrain]] <- genome_diff_vector
}

## Empty dataframe for the performance of SNP callers
snpcaller_performance_summary <- data.frame(caller=character(), mixture=character(),
                                            genomediff=integer(), calleridentify=integer(), 
                                            TP=integer(),  FP=integer(), precision=double(), 
                                            recall=double(), f1=double())

snpcaller_performance_names <- colnames(snpcaller_performance_summary)
snp_list = list()


for (file in mix_vcf){
    snpcaller_ref_sample <- unlist(strsplit(basename(file), "\\."))
    sample <- snpcaller_ref_sample[1]
    ref_name <- snpcaller_ref_sample[2]
    snpcaller <- snpcaller_ref_sample[3]
 
    mixstrain <- substr(sample, 1, 2)
    genome_diff_vector <- genome_diff_list[[mixstrain]]
    genome_diff_count <- length(genome_diff_vector)

    snp_table <- tryCatch({
        read.table(file, sep="\t", comment.char="#", header=F, 
               colClasses=c("NULL", "character", "NULL", rep("character", 2), rep("NULL", 3)))
        }, 
        error = function(e){
            snp_table <- data.frame(position=character(), ref=character(), alt=character())
    })

    #all_caller_table <- rbind(all_caller_table, snp_table)
        if (nrow(snp_table) > 0){
            snp_table <- snp_table[,1:3]
            colnames(snp_table) <- c("position", "ref", "alt")
            snp_vector <- do.call(paste, c(snp_table %>% filter(ref %in% c("A", "C", "G" ,"T"), 
                                                  alt %in% c("A", "C", "G" ,"T")), sep="-"))
            snp_identify_count <- length(snp_vector)
            tp_count <- length(intersect(snp_vector, genome_diff_vector))
            fp_count <- length(setdiff(snp_vector, genome_diff_vector))
            fn_count <- length(setdiff(genome_diff_vector, snp_vector))
            precision <- round(tp_count/snp_identify_count, 3)
            recall <- round(tp_count/genome_diff_count, 3)
            f1 <- round(2*(precision*recall)/(precision+recall), 3)
        }
        else{
            snp_vector <- c()
            snp_identify_count <- 0
            tp_count <- 0
            fp_count <- 0
            fn_count <- 0
            precision <- NA
            recall <- NA
            f1 <- NA
        }

        #snp_table$caller <- snpcaller
        #snp_table$mixture <- sample
        #snp_table$genome <- ref_name

    snpcaller_performance <- data.frame(snpcaller, sample, genome_diff_count, snp_identify_count, 
                                                 tp_count, fp_count, precision, recall, f1)
    names(snpcaller_performance) <- snpcaller_performance_names
    snpcaller_performance_summary <- rbind(snpcaller_performance_summary, snpcaller_performance)
    
    if (snpcaller %in% venn_snpcallers){
        if (sample %in% names(snp_list)){
            snp_list[[sample]][[snpcaller]] <- snp_vector      
        }
        else{
            snp_list[[sample]] <- list(genome_diff_vector, snp_vector)
            names(snp_list[[sample]]) <- c("Genome", snpcaller)
        }
    }
}


## Write the table for the performance summary of SNP callers
write.table(snpcaller_performance_summary, 
            file=snakemake@output$snpcaller_performance_summary, 
            sep="\t", row.names=F, quote=F)

combined_plots <- list()

point_plot_whole <- ggplot(snpcaller_performance_summary, aes(precision, recall, color=caller, shape=mixture)) + 
    geom_point(size=4) +
    theme_bw(base_size=13) +
    xlim(0,1) + 
    ylim(0,1)
#  guides(colour = guide_legend(override.aes = list(size=4)), shape = guide_legend(override.aes = list(size=4))) +


inset_zoom <- ggplot(snpcaller_performance_summary, aes(precision, recall, color=caller, shape=mixture)) + 
    geom_point(size=4) +
    xlim(0.8, 1) +
    ylim(0.5, 1) +
    theme_bw(base_size=12) +
    theme(legend.position="none") +
    scale_colour_brewer(palette="Set1") +
    xlab("") +
    ylab("")

rect <- data.frame(x1=0.8, x2=1, y1=0.5, y2=1)
point_plot  <- point_plot_whole + 
    annotation_custom(ggplotGrob(inset_zoom), xmin=-0.07, xmax=0.75, ymin=0.15, ymax=0.8) + 
    geom_rect(data=rect, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha=0.1, inherit.aes=F, color="#A3A3A3") + 
    geom_segment(aes(x=0.726, y=0.271, xend=1, yend=0.5), color="#A3A3A3") + 
    geom_segment(aes(x=0.083, y=0.780, xend=0.8, yend=1), color="#A3A3A3") + 
    scale_colour_brewer(palette="Set1")

box_plot <- ggplot(snpcaller_performance_summary, aes(caller, f1, fill=caller)) + 
  geom_boxplot() + geom_jitter() +
  theme_bw(base_size = 14) + 
  scale_fill_brewer(palette = "Set1")


combined_plots[[1]] <- point_plot
combined_plots[[2]] <- box_plot
n_panel <- 2
for (mix in samples){
    print(mix)
    n_panel <- n_panel + 1
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    plot <- venn.diagram(x=snp_list[[mix]],
            filename = NULL,
            main = mix,
            col="transparent",
            fill=c("blueviolet","#dc143c","dodgerblue1", "green3"),
            alpha = 0.60, 
            fontfamily = "sans",
            cat.fontfamily = "sans",
            main.fontfamily = "sans",
            cex=1.1, 
            main.cex = 1.2,
            cat.cex=1.1
    )
    combined_plots[[n_panel]] <- grobTree(plot)
}


## Arrange the Venndigram
lay <- rbind(c(1,1,1,2,2,2),
             c(1,1,1,2,2,2),
             c(1,1,1,2,2,2),
             c(3,3,4,4,5,5),
             c(3,3,4,4,5,5),
             c(6,6,7,7,8,8),
             c(6,6,7,7,8,8))

#lay <- rbind(c(1,1,2,2,3,3),
#             c(1,1,2,2,3,3),
#             c(4,4,4,5,5,5),
#             c(4,4,4,5,5,5),
#             c(4,4,4,5,5,5))
p <- grid.arrange(grobs=combined_plots, layout_matrix = lay)

ggsave(file=snakemake@output$snpcaller_performance_figure_pdf, plot=p, width=12, height=12)#, height=12, width=15)
#ggsave(plot=p, snakemake@output$snpcaller_performance_figure_png, width=12, height=12)#, height=12, width=15)

