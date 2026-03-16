version 1.0

# Individual task definitions for the qPrimer Designer pipeline.

task normalize_fasta_extensions {
    input {
        Array[File] fasta_files
    }

    command <<<
        set -euo pipefail
        mkdir -p target_seqs_normalized
        for f in ~{sep=' ' fasta_files}; do
            base=$(basename "$f")
            stem="${base%.*}"
            ext="${base##*.}"
            case "$ext" in
                fa|fasta|fna|FA|FASTA|FNA)
                    cp "$f" "target_seqs_normalized/${stem}.fa"
                    ;;
                *)
                    echo "WARNING: skipping file with unrecognized extension: $f" >&2
                    ;;
            esac
        done
    >>>

    output {
        Array[File] normalized = glob("target_seqs_normalized/*.fa")
    }

    runtime {
        docker: "ubuntu:22.04"
        memory: "1 GB"
        cpu: 1
    }
}

task count_sequences {
    input {
        File fasta
    }

    command <<<
        grep -c '^>' '~{fasta}'
    >>>

    output {
        Int count = read_int(stdout())
    }

    runtime {
        docker: "ubuntu:22.04"
        memory: "1 GB"
        cpu: 1
    }
}

task choose_target_seq {
    input {
        Int seq_count
        Int max_seqs_to_bypass_msa = 5
    }

    command <<<
        if [ ~{seq_count} -le ~{max_seqs_to_bypass_msa} ]; then
            echo "original"
        else
            echo "representative"
        fi
    >>>

    output {
        String choice = read_string(stdout())
    }

    runtime {
        docker: "ubuntu:22.04"
        memory: "1 GB"
        cpu: 1
    }
}

task make_msa {
    input {
        File fasta
        String virus_name
        Int threads = 2
    }

    command <<<
        mafft --auto --quiet --thread ~{threads} '~{fasta}' > '~{virus_name}.aln'
    >>>

    output {
        File alignment = "~{virus_name}.aln"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: threads
    }
}

task pick_representative_seqs {
    input {
        File alignment
        File params_file
        String virus_name
    }

    command <<<
        qprimer pick-representatives \
            --in '~{alignment}' \
            --out '~{virus_name}_representative.fa' \
            --params '~{params_file}' \
            --name '~{virus_name}'
    >>>

    output {
        File representative_fa = "~{virus_name}_representative.fa"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: 1
    }
}

task generate_primers {
    input {
        File fasta
        File params_file
        String virus_name
    }

    command <<<
        qprimer generate \
            --in '~{fasta}' \
            --out '~{virus_name}_init.fa' \
            --params '~{params_file}' \
            --name '~{virus_name}'
    >>>

    output {
        File primer_fasta = "~{virus_name}_init.fa"
        File primer_features = "~{virus_name}_init.feat"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: 1
    }
}

task generate_probes {
    input {
        File fasta
        File params_file
        String virus_name
    }

    command <<<
        qprimer generate-probe \
            --in '~{fasta}' \
            --out '~{virus_name}_probe.fa' \
            --params '~{params_file}' \
            --name '~{virus_name}'
    >>>

    output {
        File probe_fasta = "~{virus_name}_probe.fa"
        File probe_features = "~{virus_name}_probe.feat"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: 1
    }
}

task prepare_pset_fasta {
    input {
        File? pset_path
        String? forward_seq
        String? reverse_seq
        String pset_name
    }

    command <<<
        set -euo pipefail
        mkdir -p evaluate

        if [ -n "~{default='' pset_path}" ]; then
            cp '~{pset_path}' 'evaluate/~{pset_name}.fa'
        else
            FOR_SEQ='~{forward_seq}'
            REV_SEQ='~{reverse_seq}'
            # Reverse complement the reverse primer
            REV_RC=$(echo "$REV_SEQ" | tr 'ACGTacgt' 'TGCAtgca' | rev)

            # Build readable pair ID
            F_LEN=${#FOR_SEQ}
            R_LEN=${#REV_SEQ}
            if [ $F_LEN -ge 10 ]; then
                F_TAG=$(echo "${FOR_SEQ:0:5}${FOR_SEQ: -5}" | tr 'a-z' 'A-Z')
            else
                F_TAG=$(echo "$FOR_SEQ" | tr 'a-z' 'A-Z')
            fi
            if [ $R_LEN -ge 10 ]; then
                R_TAG=$(echo "${REV_SEQ:0:5}${REV_SEQ: -5}" | tr 'a-z' 'A-Z')
            else
                R_TAG=$(echo "$REV_SEQ" | tr 'a-z' 'A-Z')
            fi
            PAIR_ID="f_${F_TAG}_r_${R_TAG}"

            printf ">%s_for\n%s\n>%s_rev\n%s\n" "$PAIR_ID" "$FOR_SEQ" "$PAIR_ID" "$REV_RC" \
                > 'evaluate/~{pset_name}.fa'
        fi

        # Validate primer pairs
        python3 -c "
import sys
seen = {}
with open('evaluate/~{pset_name}.fa') as f:
    for line in f:
        if not line.startswith('>'):
            continue
        name = line[1:].strip()
        if name.endswith('_for'):
            seen.setdefault(name[:-4], set()).add('for')
        elif name.endswith('_rev'):
            seen.setdefault(name[:-4], set()).add('rev')
ids = sorted(pid for pid, tags in seen.items() if tags == {'for', 'rev'})
if not ids:
    print('ERROR: No valid primer pairs found (expected *_for and *_rev entries)', file=sys.stderr)
    sys.exit(1)
print(f'{len(ids)} primer set(s) will be evaluated')
# Write IDs to file for downstream use
with open('primer_ids.txt', 'w') as out:
    for i in ids:
        out.write(i + '\n')
"
    >>>

    output {
        File pset_fasta = "evaluate/~{pset_name}.fa"
        Array[String] primer_ids = read_lines("primer_ids.txt")
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "1 GB"
        cpu: 1
    }
}

task prepare_features {
    input {
        File fasta
        String output_name
    }

    command <<<
        qprimer prepare-features \
            --fa '~{fasta}' \
            --out '~{output_name}.feat'
    >>>

    output {
        File features = "~{output_name}.feat"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "2 GB"
        cpu: 1
    }
}

task build_index {
    input {
        File fasta
        String target_name
        Int threads = 1
    }

    command <<<
        mkdir -p bt2_index
        bowtie2-build --threads ~{threads} '~{fasta}' 'bt2_index/~{target_name}'
    >>>

    output {
        Array[File] index_files = glob("bt2_index/~{target_name}.*")
        File index_sentinel = "bt2_index/~{target_name}.1.bt2"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: threads
    }
}

task align {
    input {
        File fasta
        Array[File] index_files
        String target_name
        Int multi_map
        String map_option = "--mp 2,2 --rdg 4,4 --rfg 4,4 -L 8 -N 1 --score-min L,-0.6,-0.6"
        Int threads = 2
    }

    command <<<
        # Reconstruct index directory
        mkdir -p bt2_index
        for f in ~{sep=' ' index_files}; do
            ln -s "$f" "bt2_index/$(basename $f)"
        done

        bowtie2 -x 'bt2_index/~{target_name}' \
            -U '~{fasta}' -f \
            -p ~{threads} -k ~{multi_map} ~{map_option} \
            --no-hd --no-unal > output.sam
    >>>

    output {
        File sam = "output.sam"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: threads
    }
}

task align_probes_task {
    input {
        File fasta
        Array[File] index_files
        String target_name
        Int multi_map
        Boolean is_on_target
        Int probe_mismatch_penalty
        Int threads = 2
    }

    command <<<
        mkdir -p bt2_index
        for f in ~{sep=' ' index_files}; do
            ln -s "$f" "bt2_index/$(basename $f)"
        done

        if [ "~{is_on_target}" = "true" ]; then
            MAP_OPT="--end-to-end --mp ~{probe_mismatch_penalty},~{probe_mismatch_penalty} --rdg 20,20 --rfg 20,20 -N 1 --score-min L,-0.8,-0.8"
        else
            MAP_OPT="--mp 2,2 --rdg 4,4 --rfg 4,4 -L 8 -N 1 --score-min L,-0.6,-0.6"
        fi

        bowtie2 -x 'bt2_index/~{target_name}' \
            -U '~{fasta}' -f \
            -p ~{threads} -k ~{multi_map} $MAP_OPT \
            --no-hd --no-unal > output.sam
    >>>

    output {
        File sam = "output.sam"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: threads
    }
}

task parse_probe_mapping {
    input {
        File sam
        String output_name
    }

    command <<<
        qprimer parse-probe-mapping \
            --sam '~{sam}' --out '~{output_name}.probe.csv'
    >>>

    output {
        File probe_csv = "~{output_name}.probe.csv"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "2 GB"
        cpu: 1
    }
}

task parse_map {
    input {
        File sam
    }

    command <<<
        sam2pairwise < '~{sam}' > parsed.txt
    >>>

    output {
        File parsed = "parsed.txt"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: 1
    }
}

task process_map {
    input {
        File sam
        File parsed
    }

    command <<<
        set -euo pipefail
        awk 'NR % 4 == 2' '~{parsed}' > pseqs.txt
        awk 'NR % 4 == 3' '~{parsed}' > matches.txt
        awk 'NR % 4 == 0' '~{parsed}' > tseqs.txt
        cut -f1-4 '~{sam}' | paste - pseqs.txt | paste - tseqs.txt | paste - matches.txt > mapped.temp
    >>>

    output {
        File mapped_temp = "mapped.temp"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: 1
    }
}

task check_coverage {
    input {
        File mapped_temp
        Boolean is_on_target
        Int min_coverage
    }

    command <<<
        set -euo pipefail
        if [ "~{is_on_target}" = "false" ]; then
            cp '~{mapped_temp}' mapped.out
        else
            awk '{print $1, $3}' '~{mapped_temp}' | sort -u | \
                awk '{count[$1]++} END {for (i in count) print i, count[i]}' > cov.txt
            awk '$2 >= ~{min_coverage} {print $1}' cov.txt > keys.txt
            awk 'NR==FNR {keys[$1]; next} ($1 in keys)' keys.txt '~{mapped_temp}' > mapped.out
        fi
    >>>

    output {
        File mapped = "mapped.out"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: 1
    }
}

task prepare_input {
    input {
        File mapped
        File ref_fasta
        File features
        File params_file
        String ref_type
        File? prev_eval
        Boolean skip_length_filter = false
    }

    command <<<
        set -euo pipefail
        PREV_ARG=""
        if [ -n "~{default='' prev_eval}" ]; then
            PREV_ARG="--prev '~{prev_eval}'"
        fi
        SKIP_ARG=""
        if [ "~{skip_length_filter}" = "true" ]; then
            SKIP_ARG="--skip-length-filter"
        fi
        qprimer prepare-input \
            --in '~{mapped}' \
            $PREV_ARG \
            --out output.input \
            --ref '~{ref_fasta}' \
            --reftype ~{ref_type} \
            --features '~{features}' \
            --params '~{params_file}' \
            $SKIP_ARG
    >>>

    output {
        File ml_input = "output.input"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: 1
    }
}

task evaluate_ml {
    input {
        File ml_input
        File ref_fasta
        String ref_type
        Int threads = 2
    }

    command <<<
        qprimer evaluate \
            --in '~{ml_input}' --out output.eval \
            --ref '~{ref_fasta}' --reftype ~{ref_type} \
            --threads ~{threads}
    >>>

    output {
        File eval_result = "output.eval"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "8 GB"
        cpu: threads
        gpuCount: 1
        gpuType: "nvidia-tesla-t4"
    }
}

task filter_primer_list {
    input {
        File eval_result
        File init_fasta
        File params_file
        String virus_name
        File? probe_mapping_csv
        File? probe_seqs_fasta
        Boolean probe_mode
    }

    command <<<
        set -euo pipefail
        PROBE_ARGS=""
        if [ "~{probe_mode}" = "true" ]; then
            PROBE_ARGS="--probe-mapping '~{probe_mapping_csv}' --probe-seqs '~{probe_seqs_fasta}' --probe-out '~{virus_name}_filt_probes.csv'"
        fi
        qprimer filter \
            --init '~{init_fasta}' \
            --scores '~{eval_result}' \
            --out '~{virus_name}_filt.fa' \
            --params '~{params_file}' \
            $PROBE_ARGS
        touch '~{virus_name}_filt_probes.csv' '~{virus_name}_filt_probes.fa'
    >>>

    output {
        File filt_fasta = "~{virus_name}_filt.fa"
        File filt_csv = "~{virus_name}_filt.csv"
        File filt_probes_csv = "~{virus_name}_filt_probes.csv"
        File filt_probes_fasta = "~{virus_name}_filt_probes.fa"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: 1
    }
}

task build_final_output {
    input {
        File eval_on
        Array[File] eval_off
        File primer_fasta
        File ref_fasta
        File params_file
        String target_name
        File? probe_mapping_on_csv
        Array[File] probe_mapping_off_csvs = []
        File? probe_seqs_fasta
        Boolean probe_mode
    }

    command <<<
        set -euo pipefail

        OFF_ARG=""
        OFF_FILES="~{sep=' ' eval_off}"
        if [ -n "$OFF_FILES" ]; then
            OFF_ARG="--off $OFF_FILES"
        fi

        PROBE_ARGS=""
        if [ "~{probe_mode}" = "true" ]; then
            PROBE_OFF_FILES="~{sep=' ' probe_mapping_off_csvs}"
            PROBE_ARGS="--probe-mapping-on '~{probe_mapping_on_csv}'"
            if [ -n "$PROBE_OFF_FILES" ]; then
                PROBE_ARGS="$PROBE_ARGS --probe-mapping-off $PROBE_OFF_FILES"
            fi
            PROBE_ARGS="$PROBE_ARGS --probe-seqs '~{probe_seqs_fasta}'"
        fi

        qprimer build-output \
            --on '~{eval_on}' \
            $OFF_ARG \
            --fa '~{primer_fasta}' \
            --out '~{target_name}.csv' \
            --name '~{target_name}' \
            --params '~{params_file}' \
            --ref '~{ref_fasta}' \
            $PROBE_ARGS
    >>>

    output {
        File final_csv = "~{target_name}.csv"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: 1
    }
}

task select_best_multiplex {
    input {
        Array[File] csvs
        String output_name
        File? params_file
        Boolean probe_mode
    }

    command <<<
        set -euo pipefail
        PARAMS_ARG=""
        if [ "~{probe_mode}" = "true" ] && [ -n "~{default='' params_file}" ]; then
            PARAMS_ARG="--params '~{params_file}'"
        fi
        qprimer select-multiplex \
            ~{sep=' ' csvs} \
            --out '~{output_name}.csv' \
            $PARAMS_ARG
    >>>

    output {
        File multiplex_csv = "~{output_name}.csv"
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: 1
    }
}

task export_report {
    input {
        File eval_on
        Array[File] eval_off
        Array[String] primer_ids
        String output_dir_name
    }

    command <<<
        set -euo pipefail

        OFF_ARG=""
        OFF_FILES="~{sep=' ' eval_off}"
        if [ -n "$OFF_FILES" ]; then
            OFF_ARG="--off $OFF_FILES"
        fi

        qprimer export-report \
            --on '~{eval_on}' \
            $OFF_ARG \
            --out '~{output_dir_name}' \
            --names ~{sep=' ' primer_ids}
    >>>

    output {
        Array[File] reports = glob("~{output_dir_name}/*.xlsx")
    }

    runtime {
        docker: "ghcr.io/broadinstitute/qprimer_designer:latest"
        memory: "4 GB"
        cpu: 1
    }
}
