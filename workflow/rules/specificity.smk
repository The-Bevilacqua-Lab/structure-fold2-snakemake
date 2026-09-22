############################################################################
# A/C specificity of RT-stops
############################################################################

rule calculate_specificity:
    """A/C specificity of the RT-stops per ID and treatment."""
    input:
        f"{TMP}/output/se/{{id}}/combined_{{treatment}}.rtsc",
    output:
        f"{config['output_dir']}/{{id}}/specificity_{{treatment}}.csv",
    log:
        "logs/calculate_specificity/{id}_{treatment}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        output=f"{config['output_dir']}/{{id}}/specificity_{{treatment}}",
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_specificity.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 {params.workdir}/{params.script} -rtsc {input} -name {params.output} -index {params.transcriptome} >{log} 2>&1
        """


rule calculate_specificity_by_sample:
    """A/C specificity of the RT-stops per sample."""
    input:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.rtsc",
    output:
        f"{config['output_dir']}/specificity_{{sample}}.csv",
    log:
        "logs/calculate_specificity_by_sample/{sample}.log",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        output=f"{config['output_dir']}/specificity_{{sample}}",
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_specificity.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 {params.workdir}/{params.script} -rtsc {input} -name {params.output} -index {params.transcriptome} >{log} 2>&1
        """


rule plot_specificity:
    """Plot A/C specificity for every sample."""
    input:
        expand(
            "{output_dir}/specificity_{sample}.csv",
            output_dir=config["output_dir"],
            sample=SAMPLES,
        ),
    output:
        report(
            f"{config['output_dir']}/qc/specificity_plot.png",
            category="QC",
            subcategory="Specificity",
            caption="../report_captions/specificity_plot.rst",
        ),
    log:
        "logs/plot_specificity.log",
    conda:
        "../envs/plotting.yaml"
    params:
        samples=SAMPLES,
        conditions=samples["condition"].tolist(),
    shell:
        """
        python3 workflow/scripts/plot_specificity.py \
            --inputs {input} \
            --samples {params.samples} \
            --conditions {params.conditions} \
            --output {output} >{log} 2>&1
        """


# plot_specificity_comparison needs config key `specificity_comparison_runs`: a
# list of {name, output_dir, samplesheet} entries, one per pipeline run.

rule plot_specificity_comparison:
    """Combined specificity figure across several pipeline runs."""
    input:
        _specificity_comparison_inputs,
    output:
        plot=report(
            f"{config['output_dir']}/qc/specificity_comparison.png",
            category="QC",
            subcategory="Specificity",
            caption="../report_captions/specificity_comparison.rst",
        ),
        manifest=f"{config['output_dir']}/qc/specificity_comparison_manifest.tsv",
    run:
        import pandas as pd

        rows = []
        for run in config["specificity_comparison_runs"]:
            sheet = pd.read_csv(run["samplesheet"], sep="\t")
            for _, row in sheet.iterrows():
                rows.append(
                    {
                        "file": f"{run['output_dir']}/specificity_{row['sample']}.csv",
                        "group": run["name"],
                        "condition": row["condition"],
                        "replicate": row.get("ID", row["sample"]),
                    }
                )
        pd.DataFrame(rows).to_csv(output.manifest, sep="\t", index=False)
        shell(
            "Rscript workflow/scripts/plot_specificity_combined.R"
            " {output.manifest} {output.plot}"
        )
