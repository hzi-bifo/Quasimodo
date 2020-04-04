#!/usr/bin/env Rscript
library(optparse)

library(ggplot2)
library(tidyverse)

option_list <- list(
    make_option(c("-i", "--indir"),
        type = "character", default = NULL,
        help = "The input dir where the snakemake benchmark outputs are",
        metavar = "character"
    ),
    make_option(c("-o", "--out"),
        type = "character", default = "assembly_resource_benchmark",
        help = "output file name of the resource use benchmark [default= %default]",
        metavar = "character"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$indir)) {
    print_help(opt_parser)
    stop("At least one argument must be supplied (input dir)", call. = FALSE)
}

assembler_map <- list(
  spades = "SPAdes", metaspa = "metaSPA",
  metaspades = "metaSPAdes",
  tadpole = "tadpole", abyss = "ABySS",
  megahit = "megahit", ray = "Ray",
  idba = "IDBA", vicuna = "Vicuna",
  iva = "IVA", savage = "Savage"
)

header <- c(
    "time (hours)", "hours", "memory (GB)", "max_vms",
    "max_uss", "max_pss", "io_in", "io_out (GB)", "mean_load"
)
resource_use_list <- list()
for (file in dir(opt$indir, pattern = ".benchmark.txt$")) {
    sample <- gsub("\\..*$", "", gsub("\\.benchmark\\.txt", "", file))
    assembler <- assembler_map[[gsub("^[^\\.]+\\.", "", gsub("\\.benchmark\\.txt", "", file))]]
    file_path <- file.path(opt$indir, file)
    read_tsv(file_path, comment = "#", col_names = header, skip = 1) %>%
        select(-c(hours, max_pss, max_uss, max_vms, mean_load, io_in)) %>%
        mutate(
            `time (hours)` = `time (hours)` / 3600,
            `memory (GB)` = `memory (GB)` / 1000,
            `io_out (GB)` = `io_out (GB)` / 1000
        ) -> resource_use_
    resource_use_$sample <- sample
    resource_use_$assembler <- assembler
    resource_use_list[[paste(sample, assembler, sep = ".")]] <- resource_use_
}

resource_use_df <- as.data.frame(do.call(rbind, resource_use_list))[, c(4, 5, 1:3)]


write.table(resource_use_df,
    file = paste0(opt$out, ".tsv"),
    sep = "\t", row.names = F, quote = F
)

resource_use_df %>%
    gather(key = "metrics", value = "value", -c(sample, assembler)) ->
resource_use_df.long
# %>%
# filter(!(metrics == "io_out (GB)" & value > 600000)) %>%
# filter(!(metrics == "time (hours)" & value > 1000000))



# assembler.resource <- read.table("./assembler.resource.benchmark.summary.txt",
#     sep = "\t",
#     header = T,
#     check.names = F
# ) %>%
#     separate(sample.assembler, c("sample", "assembler"), sep = "\\.") %>%
#     select(-c(`h:m:s`, max_pss, max_uss, max_vms, mean_load, io_in))

# colnames(assembler.resource) <- c(
#     "sample",
#     "assembler", "time (s)",  "memory (MB)", "io_out (MB)"
# )

# assembler.resource.long <- assembler.resource %>%
#     gather(key = "metrics", value = "value", -c(sample, assembler)) %>%
#     filter(!(metrics == "io_out (MB)" & value > 60000)) %>%
#     filter(!(metrics == "time (s)" & value > 60000))


g <- ggplot(resource_use_df.long, aes(assembler, value)) +
    geom_boxplot() +
    geom_point(
        position = position_jitterdodge(
            jitter.width = 0,
            dodge.width = 0.3,
            seed = 1234
        ),
        show.legend = F,
        aes(fill = sample),
        stroke = 0,
        size = 1
    ) +
    theme_bw(base_size = 15) +
    scale_color_brewer(palette = "Set1") +
    facet_wrap(~metrics, ncol = 2, scales = "free_y") +
    theme(axis.text.x = element_text(size = 15, angle = 60, hjust = 1)) +
    xlab("") +
    ylab("") +
    scale_y_log10()

ggsave(paste0(opt$out, ".pdf"), plot = g, width = 8, height = 7)
# ggsave("assembler.resource.use.log.nocolor.png", plot = g, width = 8, height = 7)