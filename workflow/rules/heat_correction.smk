############################################################################
# Optional heat correction between two temperatures
############################################################################

# Heat correction between the two temperatures of the samplesheet's optional
# 'temperature' column (HEAT_CORRECTION etc. are set in the Snakefile).
# Transcripts are restricted to those with RT-stop coverage >= 1 in the pooled
# +DMS samples of both temperatures, and the 2-8% scale is generated once from
# the pooled lower-temperature data and reused for the higher temperature
# (Su et al. 2018 PNAS SI, "Determination of DMS reactivity" steps 3a/3b), so
# the higher temperature does not renormalize away its own heat-induced signal.

if HEAT_CORRECTION:

    rule heat_correction_plus_coverage:
        """RT-stop coverage of the pooled +DMS samples at both temperatures."""
        input:
            lower=f"{TMP}/output/se/pooled_{DISTINCT_TEMPERATURES[0]}/combined_plus.rtsc",
            higher=f"{TMP}/output/se/pooled_{DISTINCT_TEMPERATURES[1]}/combined_plus.rtsc",
        output:
            f"{config['output_dir']}/qc/heat_correction_plus_coverage.csv",
        log:
            "logs/heat_correction_plus_coverage.log",
        threads: 2
        conda:
            "../envs/rtsc_tools.yaml"
        params:
            workdir=f"{workflow.basedir}/workflow",
            script="scripts/StructureFold3/rtsc_coverage.py",
            transcriptome=os.path.abspath(TRANSCRIPTOME_UPPER),
            tmpdir=f"{TMP}/structurefold2/heat_correction_plus_coverage",
            lower_name=f"{DISTINCT_TEMPERATURES[0]}.rtsc",
            higher_name=f"{DISTINCT_TEMPERATURES[1]}.rtsc",
        shell:
            """
            mkdir -p {params.tmpdir} $(dirname {output}) $(dirname {log}) \
                && cp {input.lower} {params.tmpdir}/{params.lower_name} \
                && cp {input.higher} {params.tmpdir}/{params.higher_name} \
                && OUT_ABS=$(readlink -f {output}) \
                && LOG_ABS=$(readlink -f {log}) \
                && cd {params.tmpdir} \
                && python3 {params.workdir}/{params.script} \
                    {params.transcriptome} -f {params.lower_name} {params.higher_name} \
                    -threads {threads} -name "$OUT_ABS" >"$LOG_ABS" 2>&1 \
                && cd - >/dev/null \
                && rm -rf {params.tmpdir}
            """


    rule heat_correction_shared_transcripts:
        """List transcripts with coverage >= 1 at both temperatures."""
        input:
            f"{config['output_dir']}/qc/heat_correction_plus_coverage.csv",
        output:
            f"{config['output_dir']}/qc/heat_correction_shared_transcripts.txt",
        log:
            "logs/heat_correction_shared_transcripts.log",
        conda:
            "../envs/rtsc_tools.yaml"
        params:
            workdir=f"{workflow.basedir}/workflow",
            script="scripts/StructureFold3/coverage_overlap.py",
            tmpdir=f"{TMP}/structurefold2/heat_correction_shared_transcripts",
            threshold=1.0,
        shell:
            """
            mkdir -p {params.tmpdir} $(dirname {output}) $(dirname {log}) \
                && IN_ABS=$(readlink -f {input}) \
                && LOG_ABS=$(readlink -f {log}) \
                && cd {params.tmpdir} \
                && python3 {params.workdir}/{params.script} -f "$IN_ABS" -n {params.threshold} >"$LOG_ABS" 2>&1 \
                && cd - >/dev/null \
                && mv {params.tmpdir}/*_overlap_{params.threshold}.txt {output} \
                && rm -rf {params.tmpdir}
            """


    rule heat_correction_lower_temperature_scale:
        """Generate the 2-8% scale from the pooled lower-temperature data (Su et al. 2018, step 3a)."""
        input:
            plus=get_plus_rtsc_for_id(f"pooled_{DISTINCT_TEMPERATURES[0]}"),
            minus=f"{TMP}/output/se/pooled_{DISTINCT_TEMPERATURES[0]}/combined_minus.rtsc",
            restrict=f"{config['output_dir']}/qc/heat_correction_shared_transcripts.txt",
        output:
            react=f"{config['output_dir']}/qc/heat_correction_lower_temperature.react",
            scale=f"{config['output_dir']}/qc/heat_correction_lower_temperature.scale",
        log:
            "logs/heat_correction_lower_temperature_scale.log",
        conda:
            "../envs/rtsc_tools.yaml"
        params:
            workdir=f"{workflow.basedir}/workflow",
            script="scripts/StructureFold3/rtsc_to_react.py",
            transcriptome=TRANSCRIPTOME_UPPER,
            output_prefix=f"{config['output_dir']}/qc/heat_correction_lower_temperature",
        shell:
            """
            mkdir -p $(dirname {output.react}) \
                && python3 {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} \
                    -name {params.output_prefix} -restrict {input.restrict} >{log} 2>&1
            """


    rule rtsc_to_react_heat_shared:
        """Reactivity per ID using the shared lower-temperature scale, restricted to the shared transcripts."""
        input:
            plus=get_plus_rtsc_for_reactivity,
            minus=get_minus_rtsc_for_reactivity,
            restrict=f"{config['output_dir']}/qc/heat_correction_shared_transcripts.txt",
            scale=f"{config['output_dir']}/qc/heat_correction_lower_temperature.scale",
        output:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity_precorrection.react",
        log:
            "logs/rtsc_to_react_heat_shared/{id}.log",
        wildcard_constraints:
            id="[^/]+",
        conda:
            "../envs/rtsc_tools.yaml"
        params:
            workdir=f"{workflow.basedir}/workflow",
            script="scripts/StructureFold3/rtsc_to_react.py",
            transcriptome=TRANSCRIPTOME_UPPER,
            output_prefix=f"{config['output_dir']}/{{id}}/heat_correction/reactivity_precorrection",
        shell:
            """
            mkdir -p $(dirname {output}) \
                && python3 {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} \
                    -name {params.output_prefix} -restrict {input.restrict} -scale {input.scale} >{log} 2>&1
            """


    rule heat_correction_shared_react:
        """Restrict all pre-correction reactivities to the transcripts present in every one."""
        input:
            expand(
                f"{config['output_dir']}/{{id}}/heat_correction/reactivity_precorrection.react",
                id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
            ),
        output:
            expand(
                f"{config['output_dir']}/{{id}}/heat_correction/reactivity.react",
                id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
            ),
        log:
            "logs/heat_correction_shared_react.log",
        conda:
            "../envs/rtsc_tools.yaml"
        params:
            script=f"{workflow.basedir}/workflow/scripts/StructureFold3/react_intersect_transcripts.py",
        shell:
            """
            python3 {params.script} -in {input} -out {output} >{log} 2>&1
            """


    rule convert_pre_heat_correction_react_to_csv:
        """Convert pre-correction reactivities to CSV."""
        input:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity.react",
        output:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity.csv",
        log:
            "logs/convert_pre_heat_correction_react_to_csv/{id}.log",
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


    rule heat_correct_reactivity:
        """Rescale lower- and higher-temperature reactivities to a common overall signal."""
        input:
            lower=expand(
                f"{config['output_dir']}/{{id}}/heat_correction/reactivity.react",
                id=HEAT_CORRECTION_LOWER_IDS,
            ),
            higher=expand(
                f"{config['output_dir']}/{{id}}/heat_correction/reactivity.react",
                id=HEAT_CORRECTION_HIGHER_IDS,
            ),
        output:
            lower=expand(
                f"{config['output_dir']}/{{id}}/heat_correction/reactivity_{HEAT_CORRECTION_SUFFIX}.react",
                id=HEAT_CORRECTION_LOWER_IDS,
            ),
            higher=expand(
                f"{config['output_dir']}/{{id}}/heat_correction/reactivity_{HEAT_CORRECTION_SUFFIX}.react",
                id=HEAT_CORRECTION_HIGHER_IDS,
            ),
        log:
            f"{config['output_dir']}/qc/heat_correction_scale_factors.log",
        conda:
            "../envs/rtsc_tools.yaml"
        params:
            workdir=f"{workflow.basedir}/workflow",
            script="scripts/StructureFold3/react_heat_correct.py",
            suffix=HEAT_CORRECTION_SUFFIX,
        shell:
            """
            python3 {params.workdir}/{params.script} \
                -lower {input.lower} \
                -higher {input.higher} \
                -suffix {params.suffix} >{log} 2>&1
            """


    rule convert_heat_corrected_react_to_csv:
        """Convert heat-corrected reactivities to CSV."""
        input:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity_{HEAT_CORRECTION_SUFFIX}.react",
        output:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity_{HEAT_CORRECTION_SUFFIX}.csv",
        log:
            "logs/convert_heat_corrected_react_to_csv/{id}.log",
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


    # Comparison run: the same correction on every shared transcript, with no
    # coverage threshold (config['heat_correction_compare_all_transcripts']).
    if HEAT_CORRECTION_COMPARE_ALL_TRANSCRIPTS:

        rule heat_correction_all_transcripts_shared_react:
            """Restrict plain reactivities to the transcripts present in every ID (no coverage filter)."""
            input:
                expand(
                    f"{config['output_dir']}/{{id}}/reactivity.react",
                    id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
                ),
            output:
                expand(
                    f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity.react",
                    id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
                ),
            log:
                "logs/heat_correction_all_transcripts_shared_react.log",
            conda:
                "../envs/rtsc_tools.yaml"
            params:
                script=f"{workflow.basedir}/workflow/scripts/StructureFold3/react_intersect_transcripts.py",
            shell:
                """
                mkdir -p $(dirname {output[0]}) \
                    && python3 {params.script} -in {input} -out {output} >{log} 2>&1
                """


        rule heat_correct_reactivity_all_transcripts:
            """Heat-correct on all shared transcripts (no coverage filter)."""
            input:
                lower=expand(
                    f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity.react",
                    id=HEAT_CORRECTION_LOWER_IDS,
                ),
                higher=expand(
                    f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity.react",
                    id=HEAT_CORRECTION_HIGHER_IDS,
                ),
            output:
                lower=expand(
                    f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity_{HEAT_CORRECTION_SUFFIX}.react",
                    id=HEAT_CORRECTION_LOWER_IDS,
                ),
                higher=expand(
                    f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity_{HEAT_CORRECTION_SUFFIX}.react",
                    id=HEAT_CORRECTION_HIGHER_IDS,
                ),
            log:
                f"{config['output_dir']}/qc/heat_correction_all_transcripts_scale_factors.log",
            conda:
                "../envs/rtsc_tools.yaml"
            params:
                workdir=f"{workflow.basedir}/workflow",
                script="scripts/StructureFold3/react_heat_correct.py",
                suffix=HEAT_CORRECTION_SUFFIX,
            shell:
                """
                python3 {params.workdir}/{params.script} \
                    -lower {input.lower} \
                    -higher {input.higher} \
                    -suffix {params.suffix} >{log} 2>&1
                """


        rule convert_heat_correction_all_transcripts_to_csv:
            """Convert the all-transcripts heat-corrected reactivities to CSV."""
            input:
                f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity_{HEAT_CORRECTION_SUFFIX}.react",
            output:
                f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity_{HEAT_CORRECTION_SUFFIX}.csv",
            log:
                "logs/convert_heat_correction_all_transcripts_to_csv/{id}.log",
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
