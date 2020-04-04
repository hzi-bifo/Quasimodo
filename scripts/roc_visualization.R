library(tidyverse)

roc_file_list <- snakemake@input$roc_tsv
caller_map <- list(bcftools="BCFtools",clc="CLC",freebayes="FreeBayes",
        gatk="GATK",lofreq="LoFreq",varscan="VarScan2")

roc_df_list <- list()

for (roc_file in roc_file_list){
    caller <- caller_map[[gsub("\\/.*$", "", gsub(".*rtg\\/", "", roc_file))]]
    sample <- tools::file_path_sans_ext(basename(gsub("\\/weighted_roc\\.tsv\\.gz", "", roc_file)))
    header <- c("score", "TP_baseline", "FP",
     "TP_call", "FN", "Precision", "Recall", "F1")

    roc_df <- read_tsv(roc_file, comment = "#", col_names = header)
    roc_df$sample <- sample
    roc_df$caller <- caller
    
    roc_df_list[[paste(caller, sample, sep=".")]] <- roc_df
}


roc_df_combined <- as.data.frame(do.call(rbind, roc_df_list))[,c(10,9,1:8)]
write.table(roc_df_combined, 
            file=snakemake@output$rtg_score_performance_table, 
            sep="\t", row.names=F, quote=F)

# print(unique(roc_df_combined$caller))

callers <- unique(roc_df_combined$caller)

roc_df_combined$caller <- factor(roc_df_combined$caller, levels=callers)

roc_df_combined$mix = gsub("-.*$", "", roc_df_combined$sample)
roc_df_combined$ratio = gsub("^[^-]+-", "", roc_df_combined$sample)

p <- ggplot(roc_df_combined, aes(Recall, Precision, color=caller)) + 
geom_point(size=1.5, alpha=0.4, stroke=0) + 
theme_bw(base_size = 20) +
facet_grid(ratio~mix) +
scale_color_brewer(palette = "Set1") +
guides(color = guide_legend(override.aes = list(size = 5, alpha=1)))


ggsave(file=snakemake@output$recall_precision_pdf, plot=p, width=12, height=14)

p2 <- ggplot(roc_df_combined, aes(FP, TP_call, color=caller)) + 
geom_point(size=1.5, alpha=0.3, stroke=0) + 
theme_bw(base_size = 20) +
facet_grid(ratio~mix) +
scale_color_brewer(palette = "Set1") +
guides(color = guide_legend(override.aes = list(size = 5, alpha=1)))

ggsave(file=snakemake@output$fp_tp_pdf, plot=p2, width=12, height=14)