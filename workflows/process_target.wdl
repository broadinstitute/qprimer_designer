version 1.0

import "tasks.wdl" as tasks
import "align_and_evaluate.wdl" as aae
import "align_probes.wdl" as ap

workflow ProcessTarget {
    input {
        File target_fasta
        String virus_name
        Array[File] off_target_fastas
        Array[String] off_target_names
        Array[Int] off_target_seq_counts
        File params_file
        Boolean probe_mode
        Int probe_mismatch_penalty
        Int threads = 2
    }

    # Step 1: Decide whether to use original or representative sequences
    call tasks.count_sequences {
        input:
            fasta = target_fasta
    }

    call tasks.choose_target_seq {
        input:
            seq_count = count_sequences.count
    }

    # Step 2: MSA + representative selection (only if needed)
    if (choose_target_seq.choice == "representative") {
        call tasks.make_msa {
            input:
                fasta = target_fasta,
                virus_name = virus_name,
                threads = threads
        }

        call tasks.pick_representative_seqs {
            input:
                alignment = make_msa.alignment,
                params_file = params_file,
                virus_name = virus_name
        }
    }

    File selected_fasta = select_first([pick_representative_seqs.representative_fa, target_fasta])

    # Step 3: Generate primers
    call tasks.generate_primers {
        input:
            fasta = selected_fasta,
            params_file = params_file,
            virus_name = virus_name
    }

    # Step 4: Generate probes (if probe mode)
    if (probe_mode) {
        call tasks.generate_probes {
            input:
                fasta = selected_fasta,
                params_file = params_file,
                virus_name = virus_name
        }
    }

    # Step 5: ON-target alignment and evaluation (using init primers)
    call aae.AlignAndEvaluate as on_target_eval {
        input:
            primer_fasta = generate_primers.primer_fasta,
            target_fasta = target_fasta,
            target_name = virus_name,
            features = generate_primers.primer_features,
            params_file = params_file,
            seq_count = count_sequences.count,
            is_on_target = true,
            is_evaluate_mode = false,
            threads = threads
    }

    # Step 5b: ON-target probe alignment (if probe mode)
    if (probe_mode) {
        call ap.AlignProbes as on_target_probe_align {
            input:
                probe_fasta = select_first([generate_probes.probe_fasta]),
                target_fasta = target_fasta,
                target_name = virus_name,
                seq_count = count_sequences.count,
                is_on_target = true,
                probe_mismatch_penalty = probe_mismatch_penalty,
                threads = threads
        }
    }

    # Step 6: Filter primers
    call tasks.filter_primer_list {
        input:
            eval_result = on_target_eval.eval_result,
            init_fasta = generate_primers.primer_fasta,
            params_file = params_file,
            virus_name = virus_name,
            probe_mapping_csv = on_target_probe_align.probe_csv,
            probe_seqs_fasta = generate_probes.probe_fasta,
            probe_mode = probe_mode
    }

    # Step 7: OFF-target alignments (using filtered primers)
    scatter (idx in range(length(off_target_names))) {
        call aae.AlignAndEvaluate as off_target_eval {
            input:
                primer_fasta = filter_primer_list.filt_fasta,
                target_fasta = off_target_fastas[idx],
                target_name = off_target_names[idx],
                features = generate_primers.primer_features,
                params_file = params_file,
                seq_count = off_target_seq_counts[idx],
                is_on_target = false,
                is_evaluate_mode = false,
                prev_eval = on_target_eval.eval_result,
                threads = threads
        }

        # OFF-target probe alignment (if probe mode)
        if (probe_mode) {
            call ap.AlignProbes as off_target_probe_align {
                input:
                    probe_fasta = filter_primer_list.filt_probes_fasta,
                    target_fasta = off_target_fastas[idx],
                    target_name = off_target_names[idx],
                    seq_count = off_target_seq_counts[idx],
                    is_on_target = false,
                    probe_mismatch_penalty = probe_mismatch_penalty,
                    threads = threads
            }
        }
    }

    # Step 8: Build final output
    Array[File] off_eval_results = off_target_eval.eval_result
    Array[File?] off_probe_csvs_optional = off_target_probe_align.probe_csv

    call tasks.build_final_output {
        input:
            eval_on = on_target_eval.eval_result,
            eval_off = off_eval_results,
            primer_fasta = filter_primer_list.filt_fasta,
            ref_fasta = target_fasta,
            params_file = params_file,
            target_name = virus_name,
            probe_mapping_on_csv = on_target_probe_align.probe_csv,
            probe_mapping_off_csvs = select_all(off_probe_csvs_optional),
            probe_seqs_fasta = filter_primer_list.filt_probes_fasta,
            probe_mode = probe_mode
    }

    output {
        File final_csv = build_final_output.final_csv
    }
}
