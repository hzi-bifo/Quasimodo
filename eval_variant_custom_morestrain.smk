include: "rules/load_config_custom.smk"

try:
    vcfs = list(map(str.strip, config["vcfs"].split(",")))
except AttributeError as e:
    raise PathNotGiven(
        "The VCF files from SNP calling are not specified.")

snp_dir = results_dir + "/snp"
snpcall_dir = snp_dir + "/callers"
ref = refs[0]
alt_genomes = refs[1:]

# remove all extension after the first '.' to get the genome names
genome_names = [os.path.basename(ref).split(".", 1)[0] for ref in refs]

genome_name_file_dict = dict(zip(genome_names, refs))


labels = config['labels']
novenn = config['novenn']
mix_ratio = config['ratio']


callers = [os.path.splitext(os.path.basename(vcf))[0] for vcf in vcfs] if labels is None else [label for label in labels.split(',')]

caller_vcf_dict = dict(zip(callers, vcfs))


rule all:
    input:

        expand(snp_dir + "/nucmer/{alt_genome_name}.{ref_genome_name}.maskrepeat.variants", 
            alt_genome_name=genome_names[1:], 
            ref_genome_name=genome_names[0]),
        expand(snp_dir + "/nucmer/mix_genomes.{ref_genome_name}.maskrepeat.variants.vcf",
            ref_genome_name=genome_names[0]),
        expand(snp_dir + "/rtg/{caller}/mix_genomes.{ref_genome_name}/weighted_roc.tsv.gz", 
            caller=callers, 
            ref_genome_name=genome_names[0])


rule generate_ref_sdf:
    input:
        ref
    output:
        directory(ref + '.sdf')
    conda:
        "config/conda_env.yaml"
    shell:
        """
        rtg format {input} -o {output}
        """

rule nucmer:
    input:
        ref_genome = lambda wc:genome_name_file_dict[wc.ref_genome_name],
        alt_genome = lambda wc:genome_name_file_dict[wc.alt_genome_name]
    output:
        delta = temp(snp_dir + "/nucmer/{alt_genome_name}.{ref_genome_name}.delta"),
        mask_repeat_variants = snp_dir + "/nucmer/{alt_genome_name}.{ref_genome_name}.maskrepeat.variants"
    conda:
        "config/conda_env.yaml"
    params:
        genome_diff_prefix = snp_dir + "/nucmer/{alt_genome_name}.{ref_genome_name}"
    shell:
        """
        nucmer --prefix={params.genome_diff_prefix} {input.ref_genome} {input.alt_genome}
        
        show-snps -CTHlr <(delta-filter -r -q {output.delta}) > {output.mask_repeat_variants}
        
        """

#! The minimum genome depth of each genome is 10, otherwise some positions may not be covered and causes FN

rule create_genome_diff_vcf:
    input:
        mask_repeat_variants = expand(snp_dir + "/nucmer/{alt_genome_name}.{ref_genome_name}.maskrepeat.variants", 
            alt_genome_name = genome_names[1:], 
            ref_genome_name = genome_names[0]),
        ref = ref
    output:
        variants_vcf = snp_dir + "/nucmer/mix_genomes.{ref_genome_name}.maskrepeat.variants.vcf",
        variants_vcf_bgz = snp_dir + "/nucmer/mix_genomes.{ref_genome_name}.maskrepeat.variants.vcf.gz"
    conda:
        "config/conda_env.yaml"
    params:
        genome_diff_prefix = lambda wc: snp_dir + "/nucmer/{alt_genome_name}.{ref_genome_name}",
        mix_ratio_param = "--ratio " + mix_ratio if mix_ratio is not None else ""
    shell:
        """
        python3 program/multi_mummer2vcf.py {input.mask_repeat_variants} {params.mix_ratio_param} --output-header -n -g {input.ref} > \
            {output.variants_vcf}
        bgzip -c {output.variants_vcf} > {output.variants_vcf_bgz}
        tabix -p vcf {output.variants_vcf_bgz}

        """


rule rtg:
    input:
        ref_genome = lambda wc:genome_name_file_dict[wc.ref_genome_name],
        ref_genome_sdf = rules.generate_ref_sdf.output,
        gs_vcf = rules.create_genome_diff_vcf.output.variants_vcf_bgz,
        vcf = lambda wc: caller_vcf_dict[wc.caller]
    output:
        sorted_vcf = snp_dir + "/rtg/{caller}.{ref_genome_name, [^\.\/]+}.vcf",
        vcf_bgz = snp_dir + "/rtg/{caller}.{ref_genome_name, [^\.\/]+}.vcf.gz",
        rtg_summary = snp_dir + "/rtg/{caller, [^\.\/]+}/mix_genomes.{ref_genome_name, [^\.\/]+}/summary.txt",
        rtg_roc = snp_dir + "/rtg/{caller}/mix_genomes.{ref_genome_name, [^\.\/]+}/weighted_roc.tsv.gz"
    conda:
        "config/conda_env.yaml"
    params:
        rtg_outdir = snp_dir + "/rtg/{caller}/mix_genomes.{ref_genome_name}"
    shell:
        """
        rm -rf {params.rtg_outdir}
        awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' {input.vcf} > {output.sorted_vcf}
        bgzip -c {output.sorted_vcf} > {output.vcf_bgz}
        tabix -p vcf {output.vcf_bgz}
        rtg vcfeval -b {input.gs_vcf} -c {output.vcf_bgz} -o {params.rtg_outdir} \
            -t {input.ref_genome_sdf} --squash-ploidy -f QUAL --sample ALT
        """
