rule mutationcontext:
    input:
        vcf = lambda wc: expand(snpcall_dir + "/{snpcaller}/{sample_ref}.{snpcaller}.filtered.vcf",
                                sample_ref=[sample for sample in sample_ref if sample.startswith(
                                    wc.mix) and not "-1-0" in sample],
                                snpcaller=selected_snpcaller),
        tp = lambda wc: expand(snpcall_dir + "/{snpcaller}/tp/{sample_ref}.{snpcaller}.tp.vcf",
                               snpcaller=selected_snpcaller,
                               sample_ref=[sample for sample in sample_ref_tp if sample.startswith(wc.mix)]),
        fp = lambda wc: expand(snpcall_dir + "/{snpcaller}/fp/{sample_ref}.{snpcaller}.fp.vcf",
                               snpcaller=selected_snpcaller,
                               sample_ref=[sample for sample in sample_ref_tp if sample.startswith(wc.mix)])
    output:
        snp_profile_figure = results_dir + \
            "/final_figures/{mix}.{selected_snpcaller}.snp.profile.pdf",
        mutationcontext_figure = results_dir + \
            "/final_figures/{mix}.{selected_snpcaller}.mutationcontext.total.tp.fp.pdf",
        averaged_mc_figure = results_dir + \
            "/final_figures/{mix}.{selected_snpcaller}.mutationcontext.total.tp.fp.averaged.pdf"
    conda:
        "../config/conda_env.yaml"
    params:
        ref = lambda wc: merlin_ref if wc.mix == "TM" else ad_ref
    script:
        "../scripts/mutation_context_profile.R"
