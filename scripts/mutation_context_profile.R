library(tidyverse)
library(cowplot)
#library(vcfR)
library(SomaticSignatures)
library(MutationalPatterns)
library(reshape2)
#library(SomaticCancerAlterations)
#library(ggplot2)


vcf_list <- snakemake@input$vcf
#tp_snp <- snakemake@input$tp
#fp_snp <- snakemake@input$fp
ref <- snakemake@params$ref

ref_fh <- FaFile(ref, index= paste(ref, '.fai', sep=""))


filterVCF <- function(vcf){
  	var <- read.table(vcf, sep="\t", comment.char="#", header=F, 
			colClasses=c("NULL", "numeric", "NULL", rep("character", 2), 
						"numeric", "NULL", "character"))
	colnames(var) <- c("Position", "ref", "alt", "qual", "Frequency")
  	var_filtered <-dplyr::filter(var, ref %in% c("A", "T", "C", "G"), alt %in% c("A", "T", "C", "G"))
	var_filtered$Frequency <- as.numeric(gsub(".*AF=([01]\\.[0-9]+);.*$", "\\1", var_filtered$Frequency, perl=T))
    return(var_filtered)
}

## Allele frequency profile plot
varPlot <- function(sample, fp_vcf, tp_vcf="None"){
    fp_snp <- filterVCF(fp_vcf)
    fp_snp$type <- "FP"
    if(tp_vcf == "None"){
        all_snp <- fp_snp
    }
    else{
        tp_snp <- filterVCF(tp_vcf)
        tp_snp$type <- "TP"
        all_snp <- rbind(fp_snp, tp_snp)
    }
  	jitter <- ggplot(all_snp, aes(Position, Frequency, color=type, size=I(1))) + 
    	geom_point() + 
    	theme_bw(base_size = 10) + 
    	xlab("Genome position (bp)") + 
    	ylab(paste("AF", sample, sep=" ")) 
}


## Calculate mutation context based on VCF files
vcfMutationContext <- function(vcfFile, genomeName, genomeFaFile, study){
  vrange <- readVcfAsVRanges(vcfFile, genome=genomeName)
  vrange$study <- study
  mc <- mutationContext(vrange, genomeFaFile)
}

## to make plotting nicer
theme_ss <- function() {
  t = theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.text.y = element_text(hjust = 0.5)
   
    )
  return(t)
}

theme_small_axis <- function(x = TRUE, y = TRUE, size = 10) {
  ## decrease the x/y-axis label size
  template = element_text(size = size)
  t = theme_ss()
  if (x)
    t = t + theme(axis.text.x = template)
  if (y)
    t = t + theme(axis.text.y = template)
  return(t)
}


## Mutation context plot based on mutation motif contribution
plotSpectrum <- function(x, colorby = c("sample", "alteration")) {
  colorby = match.arg(colorby)
  ## reused part of 'meltSignatures'
  w_df = melt(x, varnames = c("motif", "sample"))
  w_df$alteration = sub("([ACGTN])([ACGTN]) .+", "\\1>\\2", w_df$motif)
  w_df$context = sub("[ACGTN][ACGTN] (.+)", "\\1", w_df$motif)
  p <- ggplot(w_df)
  p <- p + geom_bar(aes_string(x = "context", y = "value", fill = colorby),
                    stat = "identity", position = "identity")
  p <- p + facet_grid(sample ~ alteration)
  p <- p + theme(legend.position = "none")
  p <- p + xlab("Motif") + ylab("Contribution")
  p <- p + theme_ss() + theme_small_axis()
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}


## Mutation context for each sample
for(vcf in vcf_list){
    pathname <- dirname(vcf)
    basename <- basename(vcf)
    snpcaller_ref_sample <- unlist(strsplit(basename, "\\."))
    #print(snpcaller_ref_sample)
    sample <- snpcaller_ref_sample[1]
    ref_name <- snpcaller_ref_sample[2]
    snpcaller <- snpcaller_ref_sample[3]
    #print(ref_fh)
    #print(paste(sample, "(TP)", sep=" "))
    study <- ifelse(endsWith(sample, "-0-1"), paste("unmixed", sample, sep=" "), sample)
    mc <- vcfMutationContext(vcf, ref_name, ref_fh, study)#sample)

    if (!endsWith(sample, "-1-0") & !endsWith(sample, "-0-1")){
	    tp_vcf <- file.path(pathname, "tp", gsub("\\.filtered", "\\.tp", basename, perl=T))
	    fp_vcf <- file.path(pathname, "fp", gsub("\\.filtered", "\\.fp", basename, perl=T))

        tp_mc <- vcfMutationContext(tp_vcf, ref_name, ref_fh, paste(sample, "(TP)", sep=" "))
        fp_mc <- vcfMutationContext(fp_vcf, ref_name, ref_fh, paste(sample, "(FP)", sep=" "))
    }
    if (endsWith(sample, "-0-1")){
        seq_error <- mc
		snp_profile_0_1 <- varPlot(sample, vcf)
    }
    else if (endsWith(sample, "-1-1")){
        mc_1_1 <- mc
        tp_mc_1_1 <- tp_mc
        fp_mc_1_1 <- fp_mc
		#tp_snp_profile_1_1 <- varPlot(tp_vcf, paste(sample, "(TP)", sep=" "))
		#fp_snp_profile_1_1 <- varPlot(fp_vcf, paste(sample, "(FP)", sep=" "))
        snp_profile_1_1 <- varPlot(sample, fp_vcf, tp_vcf)

    }
    else if (endsWith(sample, "-1-10")){
        mc_1_10 <- mc
        tp_mc_1_10 <- tp_mc
        fp_mc_1_10 <- fp_mc
		#tp_snp_profile_1_10 <- varPlot(tp_vcf, paste(sample, "(TP)", sep=" "))
		#fp_snp_profile_1_10 <- varPlot(fp_vcf, paste(sample, "(FP)", sep=" "))
        snp_profile_1_10 <- varPlot(sample, fp_vcf, tp_vcf)
    }
    else if (endsWith(sample, "-1-50")){
        mc_1_50 <- mc
        tp_mc_1_50 <- tp_mc
        fp_mc_1_50 <- fp_mc
		#tp_snp_profile_1_50 <- varPlot(tp_vcf, paste(sample, "(TP)", sep=" "))
		#fp_snp_profile_1_50 <- varPlot(fp_vcf, paste(sample, "(FP)", sep=" "))
        snp_profile_1_50 <- varPlot(sample, fp_vcf, tp_vcf)
    }
    else{
        next
    }
}

## Combine mutation context for all, TP and FP SNPs
vr <- c(mc_1_1, mc_1_10, mc_1_50)
tp_vr <- c(tp_mc_1_1, tp_mc_1_10, tp_mc_1_50)
fp_vr <- c(fp_mc_1_1, fp_mc_1_10, fp_mc_1_50, seq_error)

sca_mm_mc = motifMatrix(vr, group = "study", normalize = TRUE)
sca_mm_mc_tp = motifMatrix(tp_vr, group = "study", normalize = TRUE)
sca_mm_mc_fp = motifMatrix(fp_vr, group = "study", normalize = TRUE)

plot_vr <- plotMutationSpectrum(vr, "study") + 
  scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

plot_tp_vr <- plotMutationSpectrum(tp_vr, "study") + 
  scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

plot_fp_vr <- plotMutationSpectrum(fp_vr, "study") + 
  scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

snp_profile_plot <- plot_grid(snp_profile_0_1, snp_profile_1_1, 
			                  snp_profile_1_10, snp_profile_1_50, 
                              ncol=1, align="v", labels = "AUTO")

mutationcontext_plot <- plot_grid(plot_vr, plot_tp_vr, plot_fp_vr, ncol=1, align="v", labels = "AUTO")

means_mm <- data.frame(all=rowMeans(sca_mm_mc), TP=rowMeans(sca_mm_mc_tp), FP=rowMeans(sca_mm_mc_fp))

means_mm_plot <- plotSpectrum(as.matrix(means_mm)) +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position="none") +
  xlab("SNPs identified by LoFreq")


ggsave(snakemake@output$snp_profile_figure, plot=snp_profile_plot, width=12, height=16)
ggsave(snakemake@output$mutationcontext_figure, plot=mutationcontext_plot, width=12, height=18)
ggsave(snakemake@output$averaged_mc_figure, plot=means_mm_plot, width=12, height=6)
