# scLightGAT

This directory contains the complete scLightGAT implementation. All
project-owned paths are derived from the location of the scripts, so the
project directory can be renamed without changing the code.

## Structure

```text
scLightGAT/
├── run_cmd.sh
├── run_script.sh
├── bar_chart.py
├── heatmap.py
├── plot_umap.py
├── result/
├── utils/
├── experiments/
│   └── abltaion/
│       ├── No_CDVAE/
│       ├── No_GAT/
│       ├── Only_LightGBM/
│       ├── result/
│       ├── plots/
│       ├── bar_chart.py
│       ├── heatmap.py
│       ├── plot_umap.py
│       ├── run_cmd.sh
│       └── run_script.sh
└── scLightGAT/
    ├── main.py
    └── scLightGAT/
        ├── CDVAE.py
        ├── LightGBM.py
        ├── GAT.py
        ├── data.py
        ├── graph.py
        ├── hparam.py
        ├── random.py
        └── scLightGAT.py
```

The production model always runs the complete C-DVAE, LightGBM and GAT
pipeline. Ablation options are not exposed by the production entry point.

## Run scLightGAT

Run the full benchmark suite, including the independent-test datasets, CAF
and the five Subtype datasets, with:

```bash
bash run_cmd.sh
```

The launcher writes below this project's `result/` directory and uses the
project's configured scLightGAT Python environment.

To run a custom train/test pair, use:

```bash
bash run_script.sh \
  --train_path /path/to/train.h5ad \
  --test_path /path/to/test.h5ad \
  --save_path ./result \
  --train_celltype Celltype_training \
  --test_celltype Celltype_training \
  --train_batch Data_source \
  --test_batch None \
  --n_runs 5
```

`--test_path` may be a directory. In that case the reference is trained once
for each run and is used to predict every `.h5ad` file in that directory.

Production runs use the single `CPU_WORKERS = 12` resource setting in
`scLightGAT/scLightGAT/hparam.py` for LightGBM, PyTorch and native numeric
libraries. The loader keeps every cell, gene and value in `X`,
but reads only the annotation columns and expression layers used by scLightGAT.

## Run ablations

The three ablations are independent source copies with fixed behavior:

- `No_CDVAE`: LightGBM and GAT without C-DVAE features.
- `No_GAT`: C-DVAE and LightGBM without GAT refinement.
- `Only_LightGBM`: DGE-selected expression directly classified by LightGBM;
  both C-DVAE and GAT are removed.

Each variant defines `CPU_WORKERS = 12` in its own `scLightGAT/hparam.py`.
Its launcher applies that value to LightGBM, PyTorch when present, and native
numeric libraries before the variant starts.

Run all three across the independent-test datasets, CAF and the five Subtype
datasets:

```bash
bash experiments/abltaion/run_cmd.sh
```

Run the main experiment first. The ablation launcher reuses its persisted
random run schedule so the full model, `No_CDVAE`, `No_GAT` and
`Only_LightGBM` are compared with paired seeds; it does not introduce a fixed
seed.

To provide paths manually, pass the normal annotation arguments to:

```bash
bash experiments/abltaion/run_script.sh \
  --train_path /path/to/train.h5ad \
  --test_path /path/to/test.h5ad \
  --save_path ./experiments/abltaion/result \
  --train_celltype Celltype_training \
  --test_celltype Celltype_training \
  --train_batch Data_source \
  --test_batch None \
  --n_runs 5
```

The ablation launcher only accepts output paths under
`experiments/abltaion/result/`. This keeps ablation runs separate from the
full-model results.

Plot the ablation results against the full model baseline with:

```bash
python3 experiments/abltaion/bar_chart.py
python3 experiments/abltaion/heatmap.py
python3 experiments/abltaion/plot_umap.py
```

These plots are written below `experiments/abltaion/plots/`. General
Subtype plots are written to `plots/bar_chart/Subtype/` and
`plots/umap/Subtype/` inside the ablation directory. Immune and non-immune
bar charts are also produced by the single `bar_chart.py` entry point under
`plots/bar_chart/immune/` and `plots/bar_chart/non-immune/`. Compartment plots
only cover datasets with valid cells in those metrics. The heatmap command
creates detailed and Overall figures for broad type, Subtype + CAF, immune,
and non-immune results.
Dataset blocks are arranged from left to right, with metrics on the x-axis,
the three ablations on the y-axis, and dataset labels above each block. Each
Overall cell is the equal-weight mean of all displayed datasets/tasks in its
family. All circles have the same size; circle color and the centered number
show the paired mean score difference. The full model is hidden and used only
as the paired score reference. Pass `--absolute` to display absolute ablation
scores instead.

## Plot results

The plotting programs read from `result/` and write under `plots/` by
default. Copy the other tools' result directories into this result tree before
plotting.

Create the main benchmark heatmaps with:

```bash
python3 heatmap.py
```

The heatmap command creates detailed and Overall figures for broad type,
Subtype + CAF, immune, and non-immune results. Heatmaps place tools on the
y-axis and dataset-grouped annotation metrics on the x-axis. Absolute scores
use a red-to-blue scale (red is low, blue is high). Circle area represents run
SD in detailed plots and between-dataset SD in Overall plots. Within each
dataset-metric column that has at least two available tools, the highest
unrounded mean score is marked with a thick black circle edge (all exact ties
are marked); missing completed results remain visible as `–`.

```bash
python3 bar_chart.py
python3 plot_umap.py
```

`bar_chart.py` is the single bar-chart entry point. It writes the standard
figures directly below `plots/bar_chart/`, Subtype figures below
`plots/bar_chart/Subtype/`, and compartment figures below
`plots/bar_chart/immune/` and `plots/bar_chart/non-immune/`.
It also writes `plots/bar_chart/Summary_Broad_Type_Average_Runtime.png`.
For each completed run, this is the equal-dataset mean of
`train_runtime_seconds + predict_runtime_seconds` across the seven broad-type
datasets. Bars show the mean across runs and error bars show one sample
standard deviation across runs.

Subtype results are read from `result/Subtype/<Broad type>/` and written to
`plots/bar_chart/Subtype/` and `plots/umap/Subtype/`.

Plotting uses only numeric run directories with `_RUN_COMPLETE.json`.
Datasets or methods that have not finished are reported as `Skipped`; they do
not stop the remaining plots.
