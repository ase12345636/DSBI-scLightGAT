#!/usr/bin/env python3
"""Create score heatmaps for the main scLightGAT benchmark.

Rows are benchmark tools and columns are the four annotation metrics grouped
by dataset. Circle color is the absolute mean score (blue high, red low) and
circle area is the sample SD. Detailed plots use SD across runs; Overall plots
use SD across dataset/task means. Missing completed results are shown as a
dash. Figures are created for broad type, Subtype + CAF, immune and non-immune
benchmark scopes.
"""

from __future__ import annotations

import argparse
import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "result"
DEFAULT_OUT_DIR = PROJECT_ROOT / "plots" / "heatmap"

ALL_METHODS = (
    "Celltypist",
    "CHETAH",
    "Seurat",
    "SingleR",
    "scGPT",
    "scLightGAT",
)

COMPARTMENTS = (("immune", "Immune"), ("non-immune", "Non-immune"))

METRICS = (
    ("accuracy", "Accuracy"),
    ("balanced_accuracy", "Balanced Acc."),
    ("f1_macro", "Macro-F1"),
    ("f1_weighted", "Weighted-F1"),
)

BROAD_DATASET_ORDER = (
    "GSE115978",
    "GSE123139",
    "GSE153935",
    "GSE166555",
    "Zhengsorted",
    "lung_full",
    "sapiens_full",
)

SUBTYPE_DATASET_ORDER = (
    "B_cell",
    "CD4_T_cell",
    "CD8_T_cell",
    "DC",
    "Plasma_cell",
)

DATASET_LABELS = {
    "B_cell": "B cell",
    "CD4_T_cell": "CD4 T cell",
    "CD8_T_cell": "CD8 T cell",
    "DC": "DC",
    "Plasma_cell": "Plasma cell",
    "CAF": "CAF",
    "lung_full": "Lung",
    "sapiens_full": "Sapiens",
}


@dataclass(frozen=True)
class DatasetSpec:
    """One displayed dataset and its path below a result root."""

    key: str
    parts: tuple[str, ...]
    label: str


def _natural_run_key(path: Path) -> tuple[int, int | str]:
    if path.name.isdigit():
        return 0, int(path.name)
    return 1, path.name


def _completion_marker_is_valid(marker_path: Path) -> bool:
    try:
        marker = json.loads(marker_path.read_text(encoding="utf-8"))
    except Exception as exc:
        warnings.warn(
            f"Skipping {marker_path.parent}: invalid completion marker ({exc})",
            stacklevel=2,
        )
        return False
    return str(marker.get("status", "complete")).casefold() == "complete"


def _read_method_metrics(
    method_dir: Path,
    method: str,
    compartment: str | None = None,
) -> dict[str, tuple[float, float]]:
    """Return mean and sample SD for each metric across completed runs."""
    metric_names = tuple(metric for metric, _ in METRICS)
    collected: dict[str, list[float]] = {
        metric: [] for metric in metric_names
    }
    if not method_dir.is_dir():
        return {metric: (np.nan, np.nan) for metric in metric_names}

    for run_dir in sorted(method_dir.iterdir(), key=_natural_run_key):
        if not run_dir.is_dir() or not run_dir.name.isdigit():
            continue
        marker_path = run_dir / "_RUN_COMPLETE.json"
        metrics_path = _metrics_path(run_dir, method, compartment)
        if not marker_path.is_file() or not metrics_path.is_file():
            continue
        if not _completion_marker_is_valid(marker_path):
            continue

        try:
            frame = pd.read_csv(metrics_path)
        except Exception as exc:
            warnings.warn(
                f"Skipping {metrics_path}: could not read CSV ({exc})",
                stacklevel=2,
            )
            continue

        lookup = {str(column).casefold(): column for column in frame.columns}
        if compartment is not None:
            count_source = lookup.get("n_test_cells")
            if count_source is not None:
                counts = pd.to_numeric(frame[count_source], errors="coerce")
                frame = frame.loc[counts > 0]
        for metric in metric_names:
            source = lookup.get(metric.casefold())
            if source is None:
                continue
            numeric = pd.to_numeric(frame[source], errors="coerce")
            numeric = numeric.replace([np.inf, -np.inf], np.nan).dropna()
            if not numeric.empty:
                # One run contributes one scalar even if a CSV unexpectedly
                # contains multiple rows; otherwise its SD would be overweighted.
                collected[metric].append(float(numeric.mean()))

    stats: dict[str, tuple[float, float]] = {}
    for metric, values in collected.items():
        if not values:
            stats[metric] = (np.nan, np.nan)
        else:
            stats[metric] = (
                float(np.mean(values)),
                float(np.std(values, ddof=1)) if len(values) > 1 else np.nan,
            )
    return stats


def _method_dir(root: Path, spec: DatasetSpec, method: str) -> Path:
    return root.joinpath(*spec.parts, method)


def _metrics_path(
    run_dir: Path,
    method: str,
    compartment: str | None = None,
) -> Path:
    if compartment is None:
        return run_dir / f"{method}_metrics.csv"
    return run_dir / compartment / f"{method}_metrics.csv"


def _candidate_roots(
    results_dir: Path,
    sclightgat_results_dir: Path,
    method: str,
) -> list[Path]:
    preferred = (
        (sclightgat_results_dir, results_dir)
        if method == "scLightGAT"
        else (results_dir, sclightgat_results_dir)
    )
    roots: list[Path] = []
    for root in preferred:
        if root not in roots:
            roots.append(root)
    return roots


def _read_preferred_metrics(
    spec: DatasetSpec,
    method: str,
    results_dir: Path,
    sclightgat_results_dir: Path,
    compartment: str | None = None,
) -> dict[str, tuple[float, float]]:
    metric_names = tuple(metric for metric, _ in METRICS)
    fallback = {metric: (np.nan, np.nan) for metric in metric_names}
    for root in _candidate_roots(
        results_dir,
        sclightgat_results_dir,
        method,
    ):
        values = _read_method_metrics(
            _method_dir(root, spec, method),
            method,
            compartment,
        )
        if any(np.isfinite(mean) for mean, _ in values.values()):
            return values
    return fallback


def _has_completed_run(
    root: Path,
    spec: DatasetSpec,
    method: str,
    compartment: str | None = None,
) -> bool:
    method_dir = _method_dir(root, spec, method)
    if not method_dir.is_dir():
        return False
    if compartment is not None:
        values = _read_method_metrics(method_dir, method, compartment)
        return any(np.isfinite(mean) for mean, _ in values.values())
    return any(
        run_dir.is_dir()
        and run_dir.name.isdigit()
        and (run_dir / "_RUN_COMPLETE.json").is_file()
        and _metrics_path(run_dir, method, compartment).is_file()
        for run_dir in method_dir.iterdir()
    )


def _spec_has_results(
    spec: DatasetSpec,
    results_dir: Path,
    sclightgat_results_dir: Path,
    compartment: str | None = None,
) -> bool:
    return any(
        _has_completed_run(root, spec, method, compartment)
        for method in ALL_METHODS
        for root in _candidate_roots(
            results_dir,
            sclightgat_results_dir,
            method,
        )
    )


def _discover_broad_specs(
    results_dir: Path,
    sclightgat_results_dir: Path,
    compartment: str | None = None,
) -> list[DatasetSpec]:
    excluded = {
        "Subtype",
        "CAF",
        "CAF_test",
        "CAF_train",
        "_shared_models",
    }
    names = set()
    for root in (results_dir, sclightgat_results_dir):
        if not root.is_dir():
            continue
        names.update(
            path.name
            for path in root.iterdir()
            if path.is_dir()
            and not path.name.startswith((".", "_"))
            and path.name not in excluded
        )

    ordered = [name for name in BROAD_DATASET_ORDER if name in names]
    ordered.extend(sorted(names - set(ordered)))
    specs = []
    for name in ordered:
        spec = DatasetSpec(
            key=name,
            parts=(name,),
            label=DATASET_LABELS.get(name, name.replace("_", " ")),
        )
        if _spec_has_results(
            spec,
            results_dir,
            sclightgat_results_dir,
            compartment,
        ):
            specs.append(spec)
    return specs


def _discover_subtype_caf_specs(
    results_dir: Path,
    sclightgat_results_dir: Path,
) -> list[DatasetSpec]:
    specs = []
    for name in SUBTYPE_DATASET_ORDER:
        spec = DatasetSpec(
            key=name,
            parts=("Subtype", name),
            label=DATASET_LABELS[name],
        )
        if _spec_has_results(spec, results_dir, sclightgat_results_dir):
            specs.append(spec)

    caf = DatasetSpec(key="CAF", parts=("CAF",), label="CAF")
    if _spec_has_results(caf, results_dir, sclightgat_results_dir):
        specs.append(caf)
    return specs


def collect_panel(
    specs: Sequence[DatasetSpec],
    results_dir: Path,
    sclightgat_results_dir: Path,
    compartment: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Collect method means and run SDs into aligned wide matrices."""
    metric_names = tuple(metric for metric, _ in METRICS)
    columns = pd.MultiIndex.from_tuples(
        [
            (spec.key, metric)
            for spec in specs
            for metric in metric_names
        ],
        names=("dataset", "metric"),
    )
    means = pd.DataFrame(
        np.nan,
        index=pd.Index(ALL_METHODS, name="method"),
        columns=columns,
        dtype=float,
    )
    stds = means.copy()
    available_counts: dict[str, int] = {}

    for spec in specs:
        count = 0
        for method in ALL_METHODS:
            values = _read_preferred_metrics(
                spec,
                method,
                results_dir,
                sclightgat_results_dir,
                compartment,
            )
            if any(np.isfinite(mean) for mean, _ in values.values()):
                count += 1
            for metric in metric_names:
                mean, std = values[metric]
                means.at[method, (spec.key, metric)] = mean
                stds.at[method, (spec.key, metric)] = std
        available_counts[spec.key] = count

    return means, stds, available_counts


def make_overall_panel(
    matrix: pd.DataFrame,
    specs: Sequence[DatasetSpec],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return equal-dataset mean and between-dataset sample SD."""
    metric_names = tuple(metric for metric, _ in METRICS)
    means = pd.DataFrame(
        np.nan,
        index=pd.Index(ALL_METHODS, name="method"),
        columns=pd.Index(metric_names, name="metric"),
        dtype=float,
    )
    stds = means.copy()
    if matrix.empty or not specs:
        empty = pd.DataFrame(
            index=pd.Index(ALL_METHODS, name="method"),
            columns=pd.MultiIndex.from_tuples(
                [],
                names=("dataset", "metric"),
            ),
            dtype=float,
        )
        return empty, empty.copy()

    dataset_keys = [spec.key for spec in specs]
    for metric in metric_names:
        columns = [(dataset, metric) for dataset in dataset_keys]
        values = matrix.loc[:, columns]
        complete = values.notna().all(axis="columns")
        means.loc[complete, metric] = values.loc[complete].mean(axis="columns")
        if len(columns) > 1:
            stds.loc[complete, metric] = values.loc[complete].std(
                axis="columns",
                ddof=1,
            )
        else:
            stds.loc[complete, metric] = np.nan

    overall_columns = pd.MultiIndex.from_tuples(
        [("Overall", metric) for metric in metric_names],
        names=("dataset", "metric"),
    )
    means.columns = overall_columns
    stds.columns = overall_columns
    return means, stds


MIN_BUBBLE_AREA = 500.0
MAX_BUBBLE_AREA = 1500.0


def _nice_std_limit(matrices: Sequence[pd.DataFrame]) -> float:
    """Choose one readable SD scale shared by every figure in this script."""
    finite = []
    for matrix in matrices:
        if matrix.empty:
            continue
        values = matrix.to_numpy(dtype=float)
        finite.extend(values[np.isfinite(values) & (values >= 0.0)])
    max_pp = max((float(value) * 100.0 for value in finite), default=0.0)
    max_pp = max(max_pp, 1.0)
    magnitude = 10.0 ** np.floor(np.log10(max_pp))
    scaled = max_pp / magnitude
    if scaled <= 1.0:
        nice = 1.0
    elif scaled <= 2.0:
        nice = 2.0
    elif scaled <= 5.0:
        nice = 5.0
    else:
        nice = 10.0
    return nice * magnitude / 100.0


def _bubble_area(std: float, std_limit: float) -> float:
    if not np.isfinite(std) or std <= 0.0:
        return MIN_BUBBLE_AREA
    fraction = float(np.clip(std / std_limit, 0.0, 1.0))
    return MIN_BUBBLE_AREA + (MAX_BUBBLE_AREA - MIN_BUBBLE_AREA) * fraction


def _format_score(value: float) -> str:
    text = f"{value:.2f}"
    return text[1:] if text.startswith("0.") else text


def _add_size_legend(
    legend_ax: plt.Axes,
    std_limit: float,
    title: str,
) -> None:
    levels = (0.0, std_limit / 2.0, std_limit)
    handles = [
        legend_ax.scatter(
            [],
            [],
            s=_bubble_area(level, std_limit),
            facecolor="#B8B8B8",
            edgecolor="#666666",
            linewidth=0.5,
        )
        for level in levels
    ]
    legend_ax.legend(
        handles,
        [f"{level * 100.0:.1f}" for level in levels],
        title=title,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.10),
        ncol=1,
        frameon=False,
        fontsize=9,
        title_fontsize=9.5,
        handlelength=3.0,
        handletextpad=0.9,
        labelspacing=1.45,
        borderaxespad=0.0,
    )


def plot_panel(
    matrix: pd.DataFrame,
    std_matrix: pd.DataFrame,
    specs: Sequence[DatasetSpec],
    save_path: Path,
    title: str,
    *,
    std_limit: float,
    std_legend_title: str,
) -> None:
    """Plot a bubble heatmap: color is mean and circle area is SD."""
    if matrix.empty or not specs:
        warnings.warn(f"Skipped {title}: no completed results", stacklevel=2)
        return
    if not matrix.index.equals(std_matrix.index) or not matrix.columns.equals(
        std_matrix.columns
    ):
        raise ValueError("Mean and SD matrices are not aligned")

    values = matrix.to_numpy(dtype=float)
    std_values = std_matrix.to_numpy(dtype=float)
    cmap = LinearSegmentedColormap.from_list(
        "benchmark_score",
        ("#D56D78", "#FAFAFA", "#689BC8"),
    )
    normalizer = Normalize(vmin=0.0, vmax=1.0)

    n_columns = values.shape[1]
    n_rows = values.shape[0]
    width = max(7.2, 0.78 * n_columns + 3.8)
    height = max(6.5, 0.65 * n_rows + 2.8)
    fig, ax = plt.subplots(figsize=(width, height))
    ax.set_facecolor("white")

    row_grid, column_grid = np.indices(values.shape)
    finite = np.isfinite(values)
    column_maxima = np.full(n_columns, np.nan, dtype=float)
    for column in range(n_columns):
        if np.any(finite[:, column]):
            column_maxima[column] = np.nanmax(values[:, column])
    is_best = finite & np.isclose(
        values,
        column_maxima[np.newaxis, :],
        rtol=1e-9,
        atol=1e-12,
    )
    # A lone available result is not framed as a "winner" because there was
    # no competing tool in that dataset-metric column.
    is_best &= (finite.sum(axis=0) >= 2)[np.newaxis, :]
    edge_colors = np.where(is_best, "#000000", "#B0B0B0")
    edge_widths = np.where(is_best, 1.5, 0.45)
    sizes = np.array(
        [
            _bubble_area(std, std_limit)
            for std in std_values.ravel()
        ],
        dtype=float,
    ).reshape(std_values.shape)
    ax.scatter(
        column_grid[finite],
        row_grid[finite],
        s=sizes[finite],
        c=values[finite],
        cmap=cmap,
        norm=normalizer,
        marker="o",
        edgecolors=edge_colors[finite],
        linewidths=edge_widths[finite],
        zorder=2,
    )
    ax.set_xlim(-0.5, n_columns - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)

    metric_labels = {metric: label for metric, label in METRICS}
    ax.set_xticks(np.arange(n_columns))
    ax.set_xticklabels(
        [metric_labels[metric] for _, metric in matrix.columns],
        rotation=42,
        ha="right",
        rotation_mode="anchor",
        fontsize=10,
    )
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(list(matrix.index), fontsize=11.5)

    ax.set_xticks(np.arange(-0.5, n_columns, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="#E1E1E1", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            value = values[row, column]
            if np.isfinite(value):
                std = std_values[row, column]
                std_label = f"{std * 100.0:.1f}" if np.isfinite(std) else "–"
                label = f"{_format_score(value)}\n± {std_label}"
            else:
                label = "–"
            ax.text(
                column,
                row,
                label,
                ha="center",
                va="center",
                fontsize=7.4,
                fontweight="semibold",
                color="#111111",
                linespacing=0.95,
                zorder=3,
            )

    block_width = len(METRICS)
    for index, spec in enumerate(specs):
        start = index * block_width
        end = start + block_width
        ax.text(
            (start + end - 1) / 2,
            1.045,
            spec.label,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=11.5,
            fontweight="bold",
        )
        if index:
            ax.axvline(start - 0.5, color="white", linewidth=5.0)
            ax.axvline(start - 0.5, color="#555555", linewidth=0.7)

    ax.set_title(title, fontsize=16, fontweight="bold", pad=48)
    mappable = ScalarMappable(norm=normalizer, cmap=cmap)
    mappable.set_array([])
    colorbar = fig.colorbar(
        mappable,
        ax=ax,
        pad=0.015,
        fraction=0.022,
        shrink=0.64,
        anchor=(0.0, 1.0),
    )
    colorbar.set_ticks((0.0, 0.5, 1.0))
    colorbar.set_ticklabels(("0.0 · low", "0.5", "1.0 · high"))
    colorbar.ax.tick_params(labelsize=9.5, colors="#111111")
    colorbar.outline.set_visible(False)

    _add_size_legend(
        colorbar.ax,
        std_limit,
        f"{std_legend_title}\nBest mean: black edge",
    )
    fig.subplots_adjust(left=0.12, right=0.945, bottom=0.27, top=0.78)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _print_coverage(
    panel_name: str,
    specs: Sequence[DatasetSpec],
    available_counts: Mapping[str, int],
) -> None:
    details = ", ".join(
        f"{spec.label}={available_counts.get(spec.key, 0)}/{len(ALL_METHODS)}"
        for spec in specs
    )
    print(f"{panel_name} available tools: {details}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results_dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help="Root containing benchmark tool results.",
    )
    parser.add_argument(
        "--sclightgat_results_dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help="Root containing scLightGAT results when stored separately.",
    )
    parser.add_argument("--out_dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    roots = {args.results_dir, args.sclightgat_results_dir}
    if not any(root.is_dir() for root in roots):
        print("Skipped: no result directory was found.")
        return

    broad_specs = _discover_broad_specs(
        args.results_dir,
        args.sclightgat_results_dir,
    )
    subtype_caf_specs = _discover_subtype_caf_specs(
        args.results_dir,
        args.sclightgat_results_dir,
    )
    compartment_specs = {
        compartment: _discover_broad_specs(
            args.results_dir,
            args.sclightgat_results_dir,
            compartment,
        )
        for compartment, _ in COMPARTMENTS
    }

    broad, broad_stds, broad_counts = collect_panel(
        broad_specs,
        args.results_dir,
        args.sclightgat_results_dir,
    )
    subtype_caf, subtype_caf_stds, subtype_caf_counts = collect_panel(
        subtype_caf_specs,
        args.results_dir,
        args.sclightgat_results_dir,
    )
    broad_overall, broad_overall_stds = make_overall_panel(
        broad,
        broad_specs,
    )
    subtype_caf_overall, subtype_caf_overall_stds = make_overall_panel(
        subtype_caf,
        subtype_caf_specs,
    )
    compartment_panels = {}
    for compartment, scope_label in COMPARTMENTS:
        specs = compartment_specs[compartment]
        matrix, stds, counts = collect_panel(
            specs,
            args.results_dir,
            args.sclightgat_results_dir,
            compartment,
        )
        overall, overall_stds = make_overall_panel(matrix, specs)
        compartment_panels[compartment] = {
            "label": scope_label,
            "specs": specs,
            "matrix": matrix,
            "stds": stds,
            "counts": counts,
            "overall": overall,
            "overall_stds": overall_stds,
        }
    detail_std_limit = _nice_std_limit((broad_stds, subtype_caf_stds))
    overall_std_limit = _nice_std_limit(
        (broad_overall_stds, subtype_caf_overall_stds)
    )
    compartment_detail_std_limit = _nice_std_limit(
        tuple(
            panel["stds"]
            for panel in compartment_panels.values()
        )
    )
    compartment_overall_std_limit = _nice_std_limit(
        tuple(
            panel["overall_stds"]
            for panel in compartment_panels.values()
        )
    )

    _print_coverage("Broad type", broad_specs, broad_counts)
    _print_coverage("Subtype + CAF", subtype_caf_specs, subtype_caf_counts)
    for compartment, scope_label in COMPARTMENTS:
        panel = compartment_panels[compartment]
        _print_coverage(
            scope_label,
            panel["specs"],
            panel["counts"],
        )

    broad_path = args.out_dir / "Heatmap_Broad_Type.png"
    subtype_caf_path = args.out_dir / "Heatmap_Subtype_CAF.png"
    broad_overall_path = args.out_dir / "Heatmap_Broad_Type_Overall.png"
    subtype_caf_overall_path = (
        args.out_dir / "Heatmap_Subtype_CAF_Overall.png"
    )
    compartment_paths = {
        compartment: (
            args.out_dir
            / f"Heatmap_{scope_label.title().replace('-', '_')}.png",
            args.out_dir
            / (
                f"Heatmap_{scope_label.title().replace('-', '_')}"
                "_Overall.png"
            ),
        )
        for compartment, scope_label in COMPARTMENTS
    }

    plot_panel(
        broad,
        broad_stds,
        broad_specs,
        broad_path,
        "Benchmark Heatmap · Broad type",
        std_limit=detail_std_limit,
        std_legend_title="SD (pp): runs",
    )
    plot_panel(
        subtype_caf,
        subtype_caf_stds,
        subtype_caf_specs,
        subtype_caf_path,
        "Benchmark Heatmap · Subtype + CAF",
        std_limit=detail_std_limit,
        std_legend_title="SD (pp): runs",
    )
    overall_spec = (DatasetSpec("Overall", ("Overall",), "Overall"),)
    plot_panel(
        broad_overall,
        broad_overall_stds,
        overall_spec,
        broad_overall_path,
        "Benchmark Heatmap · Broad type Overall",
        std_limit=overall_std_limit,
        std_legend_title=(
            "SD (pp): datasets"
        ),
    )
    plot_panel(
        subtype_caf_overall,
        subtype_caf_overall_stds,
        overall_spec,
        subtype_caf_overall_path,
        "Benchmark Heatmap · Subtype + CAF Overall",
        std_limit=overall_std_limit,
        std_legend_title=(
            "SD (pp): datasets"
        ),
    )

    for compartment, scope_label in COMPARTMENTS:
        panel = compartment_panels[compartment]
        detail_path, compartment_overall_path = compartment_paths[compartment]
        plot_panel(
            panel["matrix"],
            panel["stds"],
            panel["specs"],
            detail_path,
            f"Benchmark Heatmap · {scope_label}",
            std_limit=compartment_detail_std_limit,
            std_legend_title="SD (pp): runs",
        )
        plot_panel(
            panel["overall"],
            panel["overall_stds"],
            overall_spec,
            compartment_overall_path,
            f"Benchmark Heatmap · {scope_label} Overall",
            std_limit=compartment_overall_std_limit,
            std_legend_title="SD (pp): datasets",
        )

    output_paths = [
        broad_path,
        subtype_caf_path,
        broad_overall_path,
        subtype_caf_overall_path,
    ]
    for compartment, _ in COMPARTMENTS:
        output_paths.extend(compartment_paths[compartment])
    for path in output_paths:
        if path.is_file():
            print(f"Saved {path}")


if __name__ == "__main__":
    main()
