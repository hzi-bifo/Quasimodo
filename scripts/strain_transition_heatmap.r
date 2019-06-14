library(tidyverse)
library(cowplot)
library(dplyr)
library(SomaticSignatures)
library(MutationalPatterns)

nucmer_vcf_path <- "../../results/nucmer/HIV/vcf"
fa_path <- "../../data/HIV/fasta"


## Function to generate the mutation context matrix
vcfMutationContext <- function(vcfFile, genomeName, genomeFaFile, study){
  vrange <- readVcfAsVRanges(vcfFile, genome=genomeName)
  ref_fh <- FaFile(genomeFaFile, index= paste(genomeFaFile, '.fai', sep=""))
  vrange$study <- study
  mc <- mutationContext(vrange, ref_fh)
}


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


#transition_ratio_df <- data.frame(ref1=character(), ref2=character(), ratio=double())
transition_ratio_list <- list()
multi_mm <- list()

for(filename in dir(nucmer_vcf_path, pattern = ".vvv$")){
  vcf_filepath <- paste(nucmer_vcf_path, filename, sep="/")
  vcf_filename_split <- unlist(strsplit(filename, "\\."))
  print(filename)
  ref1_name <- vcf_filename_split[1]
  ref2_name <- vcf_filename_split[2]
  study = paste(ref1_name, ref2_name, sep=".")
  
  ref_fa <- paste(fa_path, "/", ref1_name, ".fa", sep="")
  mc <- vcfMutationContext(vcf_filepath, ref1_name, ref_fa, study)
  
  transition <- as_tibble(mc) %>% dplyr::filter(alteration=="CT" | alteration=="TC") %>%nrow()
  all <- as_tibble(mc) %>% nrow()
  transition_ratio_df_single <- data.frame(ref1=ref1_name, ref2=ref2_name, ratio=transition/all)
  #transition_ratio_df <- rbind(transition_ratio_df, transition_ratio_df_single)
  #print(all)
  mm <- motifMatrix(mc, group = "study", normalize = TRUE)
  multi_mm[[study]] <- mm
  transition_ratio_list[[study]] <- transition_ratio_df_single
  #mc_vector <- c(mc_vector, mc)
}

transition_ratio_heatmap <- do.call(rbind, transition_ratio_list)
multi_mm_profile <- do.call(cbind, multi_mm)

print(transition_ratio_heatmap)

#plotSpectrum(multi_mm_profile) + 
#  scale_fill_brewer(palette = "Set1") + 
#  xlab("Motif( SNPs by NUCmer)")
