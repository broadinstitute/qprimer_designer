version 1.0

import "tasks.wdl" as tasks

workflow AlignAndEvaluate {
    input {
        File primer_fasta
        File target_fasta
        String target_name
        File features
        File params_file
        Int seq_count
        Boolean is_on_target
        Boolean is_evaluate_mode
        File? prev_eval
        Int threads = 2
    }

    Int multi_map_raw = seq_count * 5
    Int multi_map = if multi_map_raw > 50000 then 50000 else multi_map_raw
    Int min_cov = if is_on_target then round(seq_count * 0.95) else 0
    Int build_threads = if seq_count > 1000 then 4 else 1

    call tasks.build_index {
        input:
            fasta = target_fasta,
            target_name = target_name,
            threads = build_threads
    }

    call tasks.align {
        input:
            fasta = primer_fasta,
            index_files = build_index.index_files,
            target_name = target_name,
            multi_map = multi_map,
            threads = threads
    }

    call tasks.parse_map {
        input:
            sam = align.sam
    }

    call tasks.process_map {
        input:
            sam = align.sam,
            parsed = parse_map.parsed
    }

    call tasks.check_coverage {
        input:
            mapped_temp = process_map.mapped_temp,
            is_on_target = is_on_target,
            min_coverage = min_cov
    }

    String ref_type = if is_on_target || is_evaluate_mode then "on" else "off"

    call tasks.prepare_input {
        input:
            mapped = check_coverage.mapped,
            ref_fasta = target_fasta,
            features = features,
            params_file = params_file,
            ref_type = ref_type,
            prev_eval = prev_eval,
            skip_length_filter = is_evaluate_mode
    }

    call tasks.evaluate_ml {
        input:
            ml_input = prepare_input.ml_input,
            ref_fasta = target_fasta,
            ref_type = ref_type,
            threads = threads
    }

    output {
        File eval_result = evaluate_ml.eval_result
    }
}
