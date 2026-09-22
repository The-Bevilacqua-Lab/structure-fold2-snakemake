############################################################################
# Transcript counts and abundance from the -DMS channel
############################################################################

rule rtsc_total_stops:
    """Total -DMS RT-stop counts per transcript (raw-count proxy for expression)."""
    input:
        f"{TMP}/output/se/{{id}}/combined_minus.rtsc",
    output:
        f"{config['output_dir']}/{{id}}/counts_minus.csv",
    log:
        "logs/rtsc_total_stops/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/rtsc_tools.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python3 workflow/scripts/rtsc_total_stops.py --rtsc {input} --output {output} >{log} 2>&1
        """


rule calculate_transcript_abundance:
    """RPKM/TPM abundance from -DMS RT-stops (run in a scratch dir: rtsc_abundances.py can't take an output path)."""
    input:
        f"{TMP}/output/se/{{id}}/combined_minus.rtsc",
    output:
        f"{config['output_dir']}/{{id}}/abundance_{{mode}}.csv",
    log:
        "logs/calculate_transcript_abundance/{id}_{mode}.log",
    wildcard_constraints:
        id="[^/]+",
        mode="RPKM|TPM",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_abundances.py",
        tmpdir=f"{TMP}/structurefold2/{{id}}_abundance_{{mode}}",
    shell:
        """
        (
            mkdir -p {params.tmpdir} $(dirname {output}) \
                && cp {input} {params.tmpdir}/combined_minus.rtsc \
                && cd {params.tmpdir} \
                && python3 {params.workdir}/{params.script} {wildcards.mode} -f combined_minus.rtsc \
                && cd - >/dev/null \
                && mv {params.tmpdir}/combined_minus_{wildcards.mode}.csv {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """
