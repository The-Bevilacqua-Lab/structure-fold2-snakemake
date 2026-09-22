############################################################################
# Reactivity calculation and CSV conversion
############################################################################

rule rtsc_to_react:
    """Convert RT-stop counts to 2-8% normalized reactivities (-DMS subtracted); the trim3 tail is masked to NA."""
    input:
        plus=get_plus_rtsc_for_reactivity,
        minus=get_minus_rtsc_for_reactivity,
        transcriptome=TRANSCRIPTOME_UPPER,
    output:
        f"{config['output_dir']}/{{id}}/reactivity.react",
    log:
        "logs/rtsc_to_react/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_to_react.py",
        output_prefix=f"{config['output_dir']}/{{id}}/reactivity",
        trim3=TRIM3,
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python3 {params.workdir}/{params.script} {input.minus} {input.plus} {input.transcriptome} \
                -name {params.output_prefix} -trim3 {params.trim3} >{log} 2>&1 \
            && python3 workflow/scripts/mask_trim3_react.py \
                --input {output} --output {output}.masked --trim3 {params.trim3} >>{log} 2>&1 \
            && mv {output}.masked {output}
        """


rule rtsc_to_raw_react:
    """Convert RT-stop counts to raw reactivities (no 2-8% scaling, uncapped)."""
    input:
        plus=get_plus_rtsc_for_reactivity,
        minus=get_minus_rtsc_for_reactivity,
    output:
        f"{config['output_dir']}/{{id}}/raw_reactivity.react",
    log:
        "logs/rtsc_to_raw_react/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_to_react.py",
        transcriptome=TRANSCRIPTOME_UPPER,
        output_prefix=f"{config['output_dir']}/{{id}}/raw_reactivity",
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python3 {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} \
                -name {params.output_prefix} -nrm_off -threshold 1000000 >{log} 2>&1
        """


rule convert_react_to_csv:
    """Convert reactivity.react to CSV."""
    input:
        react=f"{config['output_dir']}/{{id}}/reactivity.react",
        transcriptome=TRANSCRIPTOME_UPPER,
    output:
        f"{config['output_dir']}/{{id}}/reactivity.csv",
    log:
        "logs/convert_react_to_csv/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python3 workflow/scripts/react_to_csv.py --react {input.react} --output {output} --fasta {input.transcriptome} >{log} 2>&1
        """


rule convert_raw_react_to_csv:
    """Convert raw_reactivity.react to CSV."""
    input:
        f"{config['output_dir']}/{{id}}/raw_reactivity.react",
    output:
        f"{config['output_dir']}/{{id}}/raw_reactivity.csv",
    log:
        "logs/convert_raw_react_to_csv/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome} >{log} 2>&1
        """


# Normalization comparison (experimental, on demand): every combination of
# trim3 (notrim|trim125), ln (noln|ln) and 2-8% scaling (nonrm|nrm), always
# on each replicate's own unpooled RT-stop counts.
NORM_COMBOS = expand(
    "{trim}_{nlog}_{norm}",
    trim=["notrim", "trim125"],
    nlog=["noln", "ln"],
    norm=["nonrm", "nrm"],
)



rule rtsc_to_react_comparison:
    """Reactivities for one trim/ln/nrm option combination (experimental, on demand)."""
    input:
        plus=f"{TMP}/output/se/{{id}}/combined_plus.rtsc",
        minus=f"{TMP}/output/se/{{id}}/combined_minus.rtsc",
    output:
        f"{config['output_dir']}/{{id}}/norm_comparison/{{trim}}_{{nlog}}_{{norm}}/reactivity.react",
    log:
        "logs/rtsc_to_react_comparison/{id}_{trim}_{nlog}_{norm}.log",
    wildcard_constraints:
        trim="notrim|trim125",
        nlog="noln|ln",
        norm="nonrm|nrm",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_to_react.py",
        transcriptome=TRANSCRIPTOME_UPPER,
        output_prefix=lambda wc: f"{config['output_dir']}/{wc.id}/norm_comparison/{wc.trim}_{wc.nlog}_{wc.norm}/reactivity",
        trim_flag=lambda wc: "-trim3 125" if wc.trim == "trim125" else "",
        trim3=lambda wc: 125 if wc.trim == "trim125" else 0,
        ln_flag=lambda wc: "-ln_off" if wc.nlog == "noln" else "",
        nrm_flag=lambda wc: "-nrm_off" if wc.norm == "nonrm" else "",
    shell:
        """
        python3 {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} \
            -name {params.output_prefix} \
            {params.trim_flag} {params.ln_flag} {params.nrm_flag} >{log} 2>&1 \
        && python3 workflow/scripts/mask_trim3_react.py \
            --input {output} --output {output}.masked --trim3 {params.trim3} >>{log} 2>&1 \
        && mv {output}.masked {output}
        """


rule convert_react_comparison_to_csv:
    """Convert a norm-comparison reactivity.react to CSV."""
    input:
        f"{config['output_dir']}/{{id}}/norm_comparison/{{trim}}_{{nlog}}_{{norm}}/reactivity.react",
    output:
        f"{config['output_dir']}/{{id}}/norm_comparison/{{trim}}_{{nlog}}_{{norm}}/reactivity.csv",
    log:
        "logs/convert_react_comparison_to_csv/{id}_{trim}_{nlog}_{norm}.log",
    wildcard_constraints:
        trim="notrim|trim125",
        nlog="noln|ln",
        norm="nonrm|nrm",
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome} >{log} 2>&1
        """
