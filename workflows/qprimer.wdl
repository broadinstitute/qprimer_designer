version 1.0

import "tasks.wdl" as tasks
import "align_and_evaluate.wdl" as aae
import "process_target.wdl" as pt

# Helper task: resolve a single FASTA by name from directory
task resolve_single_fasta {
    input {
        String target_dir
        String name
    }

    command <<<
        cp "~{target_dir}/~{name}.fa" "~{name}.fa"
    >>>

    output {
        File fasta = "~{name}.fa"
    }

    runtime {
        docker: "ubuntu:22.04"
        memory: "1 GB"
        cpu: 1
    }
}

# Helper task: filter panel members to exclude a specific target
task filter_panel {
    input {
        Array[String] panel
        String exclude
    }

    command <<<
        for name in ~{sep=' ' panel}; do
            if [ "$name" != "~{exclude}" ]; then
                echo "$name"
            fi
        done
    >>>

    output {
        Array[String] filtered = read_lines(stdout())
    }

    runtime {
        docker: "ubuntu:22.04"
        memory: "1 GB"
        cpu: 1
    }
}

# Main workflow: orchestrates the full pipeline for singleplex, multiplex, or evaluate mode.
workflow qprimer_pipeline {
    input {
        # Mode selection
        Boolean multiplex = false
        Boolean evaluate = false
        Boolean probe_mode = false

        # Singleplex inputs
        Array[String] targets = []
        Array[String] cross = []
        Array[String] host = []

        # Multiplex inputs
        Array[String] panel = []

        # Evaluate inputs
        String? forward_seq
        String? reverse_seq
        File? pset_file

        # Directory containing target FASTA files (each named <target>.fa)
        String target_dir = "target_seqs/original"

        # Parameters
        File params_file

        # Runtime
        Int threads = 2

        # Probe mismatch penalty (derived from PROBE_MAX_MISMATCHES in params)
        Int probe_mismatch_penalty = 10
    }

    # =========================================================================
    # EVALUATE MODE
    # =========================================================================
    if (evaluate) {
        String eval_pset_name = if defined(pset_file)
            then basename(select_first([pset_file]), ".fa")
            else "eval_primers"

        String primary_target = targets[0]

        call resolve_single_fasta as eval_resolve_primary {
            input:
                target_dir = target_dir,
                name = primary_target
        }

        call tasks.prepare_pset_fasta {
            input:
                pset_path = pset_file,
                forward_seq = forward_seq,
                reverse_seq = reverse_seq,
                pset_name = eval_pset_name
        }

        call tasks.prepare_features {
            input:
                fasta = prepare_pset_fasta.pset_fasta,
                output_name = eval_pset_name
        }

        call tasks.count_sequences as eval_count_primary {
            input:
                fasta = eval_resolve_primary.fasta
        }

        call aae.AlignAndEvaluate as eval_on_target {
            input:
                primer_fasta = prepare_pset_fasta.pset_fasta,
                target_fasta = eval_resolve_primary.fasta,
                target_name = primary_target,
                features = prepare_features.features,
                params_file = params_file,
                seq_count = eval_count_primary.count,
                is_on_target = true,
                is_evaluate_mode = true,
                threads = threads
        }

        Array[String] eval_off_targets = flatten([cross, host])

        scatter (off_name in eval_off_targets) {
            call resolve_single_fasta as eval_resolve_off {
                input:
                    target_dir = target_dir,
                    name = off_name
            }

            call tasks.count_sequences as eval_count_off {
                input:
                    fasta = eval_resolve_off.fasta
            }

            call aae.AlignAndEvaluate as eval_off_target {
                input:
                    primer_fasta = prepare_pset_fasta.pset_fasta,
                    target_fasta = eval_resolve_off.fasta,
                    target_name = off_name,
                    features = prepare_features.features,
                    params_file = params_file,
                    seq_count = eval_count_off.count,
                    is_on_target = false,
                    is_evaluate_mode = true,
                    threads = threads
            }
        }

        call tasks.export_report {
            input:
                eval_on = eval_on_target.eval_result,
                eval_off = eval_off_target.eval_result,
                primer_ids = prepare_pset_fasta.primer_ids,
                output_dir_name = eval_pset_name
        }
    }

    # =========================================================================
    # SINGLEPLEX MODE
    # =========================================================================
    if (!evaluate && !multiplex) {
        Array[String] singleplex_off_names = flatten([cross, host])

        # Resolve off-target FASTAs and count sequences
        scatter (off_name in singleplex_off_names) {
            call resolve_single_fasta as singleplex_resolve_off {
                input:
                    target_dir = target_dir,
                    name = off_name
            }

            call tasks.count_sequences as singleplex_count_off {
                input:
                    fasta = singleplex_resolve_off.fasta
            }
        }

        scatter (virus in targets) {
            call resolve_single_fasta as singleplex_resolve_target {
                input:
                    target_dir = target_dir,
                    name = virus
            }

            call pt.ProcessTarget {
                input:
                    target_fasta = singleplex_resolve_target.fasta,
                    virus_name = virus,
                    off_target_fastas = singleplex_resolve_off.fasta,
                    off_target_names = singleplex_off_names,
                    off_target_seq_counts = singleplex_count_off.count,
                    params_file = params_file,
                    probe_mode = probe_mode,
                    probe_mismatch_penalty = probe_mismatch_penalty,
                    threads = threads
            }
        }
    }

    # =========================================================================
    # MULTIPLEX MODE
    # =========================================================================
    if (!evaluate && multiplex) {
        scatter (on_target in panel) {
            call resolve_single_fasta as multiplex_resolve_target {
                input:
                    target_dir = target_dir,
                    name = on_target
            }

            # Get other panel members as off-targets
            call filter_panel {
                input:
                    panel = panel,
                    exclude = on_target
            }

            Array[String] all_off_names = flatten([filter_panel.filtered, host])

            # Resolve all off-target FASTAs for this panel member
            scatter (off_name in all_off_names) {
                call resolve_single_fasta as multiplex_resolve_off {
                    input:
                        target_dir = target_dir,
                        name = off_name
                }

                call tasks.count_sequences as multiplex_count_off {
                    input:
                        fasta = multiplex_resolve_off.fasta
                }
            }

            call pt.ProcessTarget as multiplex_process {
                input:
                    target_fasta = multiplex_resolve_target.fasta,
                    virus_name = on_target,
                    off_target_fastas = multiplex_resolve_off.fasta,
                    off_target_names = all_off_names,
                    off_target_seq_counts = multiplex_count_off.count,
                    params_file = params_file,
                    probe_mode = probe_mode,
                    probe_mismatch_penalty = probe_mismatch_penalty,
                    threads = threads
            }
        }

        call tasks.select_best_multiplex {
            input:
                csvs = multiplex_process.final_csv,
                output_name = "multiplex_output",
                params_file = params_file,
                probe_mode = probe_mode
        }
    }

    # =========================================================================
    # OUTPUTS
    # =========================================================================
    output {
        # Evaluate mode outputs
        Array[File]? evaluation_reports = export_report.reports

        # Singleplex mode outputs
        Array[File]? singleplex_csvs = ProcessTarget.final_csv

        # Multiplex mode outputs
        Array[File]? multiplex_per_target_csvs = multiplex_process.final_csv
        File? multiplex_combined_csv = select_best_multiplex.multiplex_csv
    }
}
