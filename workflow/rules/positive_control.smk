############################################################################
# Positive-control structure plot: p4p6 ONLY.
#
# Colors p4p6 DMS reactivity onto the Tetrahymena group I intron P4-P6
# domain secondary structure (workflow/scripts/qc_plots/plot_p4p6_structure.py, matplotlib), using the XRNA drawing in
# resources/p4p6/p4p6.xrna (structure + nucleotide coordinates). The p4p6
# construct (resources/p4p6/p4p6.fa) carries 2 extra 5' nucleotides, so
# reactivity position p is drawn at XRNA nucleotide p - 2.
############################################################################

P4P6_XRNA = "resources/p4p6/p4p6.xrna"


rule plot_p4p6_react:
    """p4p6 reactivity on the P4-P6 structure, one ID."""
    input:
        react=f"{config['output_dir']}/{{id}}/reactivity.csv",
        xrna=P4P6_XRNA,
    output:
        report(
            f"{config['output_dir']}/{{id}}/p4p6_react.png",
            category="Positive control (p4p6)",
            subcategory="Per-replicate",
            labels={"id": "{id}"},
            caption="../report_captions/p4p6_react.rst",
        ),
    log:
        "logs/plot_p4p6_react/{id}.log",
    conda:
        "../envs/plotnine.yaml"
    shell:
        """
        python3 workflow/scripts/qc_plots/plot_p4p6_structure.py \
            --inputs {input.react} --labels {wildcards.id} --xrna {input.xrna} \
            --output {output} >{log} 2>&1
        """


rule plot_p4p6_all_samples:
    """p4p6 reactivity on the P4-P6 structure, every ID side by side."""
    input:
        react=expand(f"{config['output_dir']}/{{id}}/reactivity.csv", id=IDS),
        xrna=P4P6_XRNA,
    output:
        report(
            f"{config['output_dir']}/qc/p4p6_all_samples.png",
            category="Positive control (p4p6)",
            subcategory="All replicates",
            caption="../report_captions/p4p6_all_samples.rst",
        ),
    log:
        "logs/plot_p4p6_all_samples.log",
    conda:
        "../envs/plotnine.yaml"
    params:
        labels=IDS,
    shell:
        """
        python3 workflow/scripts/qc_plots/plot_p4p6_structure.py \
            --inputs {input.react} --labels {params.labels} --xrna {input.xrna} \
            --output {output} >{log} 2>&1
        """


rule plot_p4p6_roc:
    """ROC of p4p6 reactivity vs the P4-P6 structure: unpaired A/C (positive) vs WCF-paired A/C (negative)."""
    input:
        react=expand(f"{config['output_dir']}/{{id}}/reactivity.csv", id=IDS),
        xrna=P4P6_XRNA,
    output:
        plot=report(
            f"{config['output_dir']}/qc/p4p6_roc.png",
            category="Positive control (p4p6)",
            subcategory="ROC",
        ),
        auc=f"{config['output_dir']}/qc/p4p6_auc.csv",
    log:
        "logs/plot_p4p6_roc.log",
    conda:
        "../envs/plotnine.yaml"
    params:
        labels=IDS,
        trim3=P4P6_ROC_TRIM3,
    shell:
        """
        python3 workflow/scripts/qc_plots/plot_p4p6_roc.py \
            --inputs {input.react} --labels {params.labels} --xrna {input.xrna} --trim3 {params.trim3} \
            --plot {output.plot} --auc-csv {output.auc} >{log} 2>&1
        """
