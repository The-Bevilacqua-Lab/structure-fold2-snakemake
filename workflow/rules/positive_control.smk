############################################################################
# Positive-control structure plot: p4p6 ONLY.
#
# These rules color reactivity directly onto a hand-annotated diagram of
# the p4p6 (Tetrahymena group I intron P4-P6 domain) secondary structure
# (resources/p4p6/p4p6_structure.png + p4p6_pairing.txt), using a fixed
# pixel-coordinate map that is specific to that one construct and image.
# There is no way to generalize this to an arbitrary positive control or
# organism without a bespoke coordinate map for it, so this file is only
# included when config['positive_control_name'] == "p4p6" (see the main
# Snakefile). For any other (or no) positive control, only the generic
# alignment-% QC in workflow/rules/common.smk applies.
############################################################################

rule plot_p4p6_react:
    input:
        f"{config['output_dir']}/{{id}}/reactivity.csv"
    output:
        f"{config['output_dir']}/{{id}}/p4p6_react.png"
    conda:
        "../envs/plotting.yaml"
    shell:
        "python3 workflow/scripts/plot_p4p6.py --input {input} --output {output}"


rule plot_p4p6_all_samples:
    input:
        expand(f"{config['output_dir']}/{{id}}/reactivity.csv", id=IDS)
    output:
        f"{config['output_dir']}/qc/p4p6_all_samples.png"
    conda:
        "../envs/plotting.yaml"
    shell:
        "python3 workflow/scripts/plot_p4p6_all_samples.py "
        "--input {input} "
        "--samples {IDS} "
        "--output {output}"


rule plot_p4p6_react_plus_only:
    """
    Same as plot_p4p6_react, but for the +DMS-only (no -DMS subtraction)
    reactivity.
    """
    input:
        f"{config['output_dir']}/{{id}}/reactivity_plus_only.csv"
    output:
        f"{config['output_dir']}/{{id}}/p4p6_react_plus_only.png"
    conda:
        "../envs/plotting.yaml"
    shell:
        "python3 workflow/scripts/plot_p4p6.py --input {input} --output {output}"


rule plot_p4p6_all_samples_plus_only:
    """
    Same as plot_p4p6_all_samples, but for the +DMS-only (no -DMS
    subtraction) reactivity.
    """
    input:
        expand(f"{config['output_dir']}/{{id}}/reactivity_plus_only.csv", id=IDS)
    output:
        f"{config['output_dir']}/qc/p4p6_all_samples_plus_only.png"
    conda:
        "../envs/plotting.yaml"
    shell:
        "python3 workflow/scripts/plot_p4p6_all_samples.py "
        "--input {input} "
        "--samples {IDS} "
        "--output {output}"
