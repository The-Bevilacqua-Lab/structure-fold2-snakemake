############################################################################
# Reactivity calculation and CSV conversion
############################################################################

rule rtsc_to_react:
    """Convert RT-stop counts to 2-8% normalized reactivities (-DMS subtracted)."""
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
        norm=REACTIVITY_NORM_ARGS,
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python3 {params.workdir}/{params.script} {input.minus} {input.plus} {input.transcriptome} \
                -name {params.output_prefix} {params.norm} >{log} 2>&1
        """


rule rtsc_to_raw_react:
    """Convert RT-stop counts to raw reactivities (no 2-8% scaling, uncapped)."""
    input:
        plus=get_plus_rtsc_for_reactivity,
        minus=get_minus_rtsc_for_reactivity,
        transcriptome=TRANSCRIPTOME_UPPER,
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
        output_prefix=f"{config['output_dir']}/{{id}}/raw_reactivity",
        norm=REACTIVITY_NORM_ARGS,
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python3 {params.workdir}/{params.script} {input.minus} {input.plus} {input.transcriptome} \
                -name {params.output_prefix} -nrm_off -threshold 1000000 {params.norm} >{log} 2>&1
        """


rule react_to_csv:
    """
    Any reactivity .react (standard, raw, heat-corrected, norm comparison) -> CSV with columns
    transcript, position, base, reactivity (NA where the .react has no value).
    """
    input:
        react=f"{config['output_dir']}/{{path}}.react",
        transcriptome=TRANSCRIPTOME_UPPER,
    output:
        f"{config['output_dir']}/{{path}}.csv",
    log:
        "logs/react_to_csv/{path}.log",
    wildcard_constraints:
        # <id>/.../reactivity[_<suffix>] or raw_reactivity; depth_matched/ writes its CSVs directly
        path=r"(?!depth_matched/)(?:[^/]+/)+(?:raw_)?reactivity(?:_[^/]+)?",
    conda:
        "../envs/biopython.yaml"
    shell:
        "python3 workflow/scripts/reactivity/react_to_csv.py --react {input.react} --output {output} "
        "--fasta {input.transcriptome} > {log} 2>&1"


# RNA Framework XML (config rnaframework_xml): the same reactivities in the rf-norm output format,
# one <transcript>.xml per transcript with at least one value, readable by rf-fold, rf-correlate,
# rf-compare etc. Attributes: Ding scoring (ln(count + 1)), 2-8% over the whole transcript, AC, cap 7.
RF_XML_SCRIPT = "workflow/scripts/reactivity/react_to_rnaframework_xml.py"


rule react_to_rnaframework_xml:
    """reactivity.react -> rnaframework_xml/<transcript>.xml (rf-norm format)."""
    input:
        react=f"{config['output_dir']}/{{id}}/reactivity.react",
        transcriptome=TRANSCRIPTOME_UPPER,
    output:
        directory(f"{config['output_dir']}/{{id}}/rnaframework_xml"),
    log:
        "logs/react_to_rnaframework_xml/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python3 {RF_XML_SCRIPT} --react {input.react} --fasta {input.transcriptome} --outdir {output} >{log} 2>&1
        """


# Normalization comparison (experimental, on demand): every combination of
# trim3 (notrim|trim125), ln (noln|ln) and 2-8% scaling (nonrm|nrm), always
# on each replicate's own unpooled RT-stop counts.



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
        && python3 workflow/scripts/reactivity/mask_trim3_react.py \
            --input {output} --output {output}.masked --trim3 {params.trim3} >>{log} 2>&1 \
        && mv {output}.masked {output}
        """


