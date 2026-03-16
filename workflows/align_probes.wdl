version 1.0

import "tasks.wdl" as tasks

workflow AlignProbes {
    input {
        File probe_fasta
        File target_fasta
        String target_name
        Int seq_count
        Boolean is_on_target
        Int probe_mismatch_penalty
        Int threads = 2
    }

    Int multi_map_raw = seq_count * 5
    Int multi_map = if multi_map_raw > 50000 then 50000 else multi_map_raw
    Int build_threads = if seq_count > 1000 then 4 else 1

    call tasks.build_index {
        input:
            fasta = target_fasta,
            target_name = target_name,
            threads = build_threads
    }

    call tasks.align_probes_task {
        input:
            fasta = probe_fasta,
            index_files = build_index.index_files,
            target_name = target_name,
            multi_map = multi_map,
            is_on_target = is_on_target,
            probe_mismatch_penalty = probe_mismatch_penalty,
            threads = threads
    }

    call tasks.parse_probe_mapping {
        input:
            sam = align_probes_task.sam,
            output_name = target_name
    }

    output {
        File probe_csv = parse_probe_mapping.probe_csv
    }
}
