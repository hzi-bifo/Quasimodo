library(tidyverse)
library(VennDiagram)
library(gridExtra)

## VCF files for mixtures
mix_vcf <- snakemake@input$vcfs

## SNP callers to compare using Venndigram
venn_snpcallers <- unlist(snakemake@params$snp_venn)

## The genome difference by NUCmer
genome_diff_file <- snakemake@input$genome_diff

## Make a list for the genome difference


genome_diff_table <- read.table(genome_diff_file, sep="\t", comment.char="#", header=F,
                                colClasses=c(rep("character", 3), rep("NULL", 9)))
colnames(genome_diff_table) <- c("position", "ref", "alt")
genome_diff_table_filtered <- genome_diff_table %>% filter(ref!=".", alt!=".")
genome_diff_vector <- do.call(paste, c(genome_diff_table_filtered, sep="-"))


## Empty dataframe for the performance of SNP callers
snpcaller_performance_summary <- data.frame(caller=character(),
                                            genomediff=integer(), calleridentify=integer(), 
                                            TP=integer(),  FP=integer(), precision=double(), 
                                            recall=double(), f1=double())

snpcaller_performance_names <- colnames(snpcaller_performance_summary)
snp_list = list()
snp_list[["Genome"]] = genome_diff_vector

for (file in mix_vcf){
    snpcaller <- gsub('\\.filtered\\.vcf$', '', basename(file), ignore.case=TRUE)
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

    snpcaller_performance <- data.frame(snpcaller, genome_diff_count, snp_identify_count, 
                                                 tp_count, fp_count, precision, recall, f1)
    names(snpcaller_performance) <- snpcaller_performance_names
    snpcaller_performance_summary <- rbind(snpcaller_performance_summary, snpcaller_performance)
    snp_list[[snpcaller]] <- snp_vector
}


## Write the table for the performance summary of SNP callers
write.table(snpcaller_performance_summary, 
            file=snakemake@output$snp_benchmark_table, 
            sep="\t", row.names=F, quote=F)

print("Point plot whole")
point_plot <- ggplot(snpcaller_performance_summary, aes(precision, recall, color=caller)) + 
    geom_point(size=4) +
    theme(legend.position = 'right',
        legend.spacing.x = unit(0.3, 'cm'),
        text = element_text(size=14)) +
    theme_bw(base_size=15) 


futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
# Venndiagram
venn <- venn.diagram(x=snp_list,
        filename = NULL,
        col="transparent",
        fill=c("blueviolet","#dc143c","dodgerblue1", "green3"),
        alpha = 0.60, 
        fontfamily = "sans",
        cat.fontfamily = "sans",
        main.fontfamily = "sans",
        cex=1.2, 
        main.cex = 1.3,
        cat.cex=1.2
)

ggsave(file=snakemake@output$snp_benchmark_figure, plot=point_plot, width=8, height=6)
ggsave(file=snakemake@output$snp_venn_figure, plot=grobTree(venn), width=6, height=6)