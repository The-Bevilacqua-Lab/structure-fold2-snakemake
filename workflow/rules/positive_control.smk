############################################################################
# Positive-control structure plot: p4p6 ONLY.
#
# These rules color reactivity directly onto a hand-annotated diagram of
# the p4p6 (Tetrahymena group I intron P4-P6 domain) secondary structure
# (resources/p4p6/p4p6_structure.png + p4p6_pairing.txt), using a fixed
# pixel-coordinate map that is specific to that one construct and image.
############################################################################


rule plot_p4p6_react:
    input:
        f"{config['output_dir']}/{{id}}/reactivity.csv",
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
        "../envs/plotting.yaml"
    shell:
        "python3 workflow/scripts/plot_p4p6.py --input {input} --output {output} > {log} 2>&1"


rule plot_p4p6_all_samples:
    input:
        expand(f"{config['output_dir']}/{{id}}/reactivity.csv", id=IDS),
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
        "../envs/plotting.yaml"
    shell:
        "python3 workflow/scripts/plot_p4p6_all_samples.py "
        "--input {input} "
        "--samples {IDS} "
        "--output {output} > {log} 2>&1"
