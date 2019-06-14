rule metaquast:
    input:
        scaffolds = lambda wc: expand(assembly_dir + "/{assembler}/{{sample}}.{assembler}.scaffolds.fa", assembler=assemblers)
    output:
        report = metaquast_dir + "/{mix}/{sample, [A-Z]+-[0-9\-]+}/report.html"
    params:
        metaquast_outdir = metaquast_dir + "/{mix}/{sample}",
        ref = lambda wc: ",".join([tb_ref, ad_ref]) if wc.mix == "TA" else ",".join([tb_ref, merlin_ref])
    shell:
        """
        metaquast.py --unique-mapping -o {params.metaquast_outdir} -R {params.ref} {input.scaffolds}
        """
