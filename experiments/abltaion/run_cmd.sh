#!/usr/bin/env bash
set -euo pipefail

ablation_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$ablation_root/../.." && pwd)"
main_result_root="$project_root/result"
ablation_result_root="$ablation_root/result"

if [[ -d "$project_root/preprocessing/scLightGAT_data" ]]; then
    data_root="$project_root/preprocessing/scLightGAT_data"
else
    data_root="$(cd "$project_root/.." && pwd)/preprocessing/scLightGAT_data"
fi

sync_seed_schedule() {
    local source_result="$1"
    local target_result="$2"
    local source_seed_dir="$source_result/.benchmark_state/seeds"
    local target_seed_dir="$target_result/.benchmark_state/seeds"

    for run_id in 1 2 3 4 5; do
        if [[ ! -s "$source_seed_dir/run_${run_id}.json" ]]; then
            echo "Error: missing main-experiment seed: $source_seed_dir/run_${run_id}.json" >&2
            echo "Run the main experiment before the ablation suite." >&2
            exit 2
        fi
    done

    mkdir -p "$target_seed_dir"
    for run_id in 1 2 3 4 5; do
        cp "$source_seed_dir/run_${run_id}.json" "$target_seed_dir/run_${run_id}.json"
    done
}

sync_seed_schedule "$main_result_root" "$ablation_result_root"

bash "$ablation_root/run_script.sh" \
    --train_path "$data_root/Integrated_training/train.h5ad" \
    --test_path "$data_root/Independent_testing" \
    --save_path "$ablation_result_root" \
    --train_celltype Celltype_training \
    --test_celltype Celltype_training \
    --train_batch Data_source \
    --test_batch None \
    --test_batch_map "GSE115978=samples,GSE123139=Sample.Short.Name,GSE153935=orig.ident,GSE166555=sample_id,Zhengsorted=None,lung_full=None,sapiens_full=donor" \
    --n_runs 5

sync_seed_schedule \
    "$main_result_root/CAF" \
    "$ablation_result_root/CAF"

bash "$ablation_root/run_script.sh" \
    --train_path "$data_root/caf.data/caf_train.h5ad" \
    --test_path "$data_root/caf.data/caf_test.h5ad" \
    --save_path "$ablation_result_root/CAF" \
    --train_celltype Celltype_training \
    --test_celltype Celltype_training \
    --train_batch SampleID \
    --test_batch None \
    --n_runs 5

subtype_data="$data_root/Subtype"
subtype_result="$ablation_result_root/Subtype"

for broad_type in B_cell CD4_T_cell CD8_T_cell DC Plasma_cell; do
    sync_seed_schedule \
        "$main_result_root/Subtype/$broad_type" \
        "$subtype_result/$broad_type"

    bash "$ablation_root/run_script.sh" \
        --train_path "$subtype_data/${broad_type}_train.h5ad" \
        --test_path "$subtype_data/${broad_type}_test.h5ad" \
        --save_path "$subtype_result/$broad_type" \
        --train_celltype Celltype_subtraining \
        --test_celltype Celltype_subtraining \
        --train_batch Data_source \
        --test_batch Data_source \
        --n_runs 5
done
