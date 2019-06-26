library(tidyverse)
library(VennDiagram)
library(gridExtra)

## FP SNP vcf
fp_vcf <- snakemake@input$fp

## Sample list
samples <- unlist(snakemake@params$mix_sample)

## SNP callers to compare using Venndigram
venn_snpcallers <- unlist(snakemake@params$fp_compared_snpcallers)

snp_list = list()


for (file in fp_vcf){
    snpcaller_ref_sample <- unlist(strsplit(basename(file), "\\."))
    sample <- snpcaller_ref_sample[1]
    ref_name <- snpcaller_ref_sample[2]
    snpcaller <- snpcaller_ref_sample[3]
 
    mixstrain <- substr(sample, 1, 2)

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
    }
    else{
        snp_vector <- c()
        snp_identify_count <- 0
    }

    if (snpcaller %in% venn_snpcallers){
        if (sample %in% names(snp_list)){
            snp_list[[sample]][[snpcaller]] <- snp_vector      
        }
        else{
            snp_list[[sample]] <- list(snp_vector)
            names(snp_list[[sample]]) <- c(snpcaller)
        }
    }
}

combined_plots <- list()
n_panel <- 0
for (mix in samples){
    print(mix)
    n_panel <- n_panel + 1
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    # Venndiagram
    plot <- venn.diagram(x=snp_list[[mix]],
            filename = NULL,
            main = mix,
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
    combined_plots[[n_panel]] <- grobTree(plot)
}

## Arrange the Venndigram
lay <- rbind(c(1,2,3),
             c(4,5,6))

print("Grid all plots")
p <- grid.arrange(grobs=combined_plots, layout_matrix = lay)
ggsave(file=snakemake@output$fp_compare_figure, plot=p, width=12, height=8)