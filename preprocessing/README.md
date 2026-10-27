# scLightGAT Data Preprocessing

The datasets used by scLightGAT are assembled from multiple public sources, including the Gene Expression Omnibus (GEO) and data archives released with the corresponding publications. This pipeline automatically downloads the required files, performs standardized preprocessing, and harmonizes cell-type labels across sources. Training annotations are mapped to the predefined scLightGAT label set. Test annotations are mapped to the most specific supported class within the broad cell-type ontology; when a supported source annotation is coarser than this ontology, it is retained rather than inferred or subdivided. Coarse labels such as `T.cells` are not retained in the final training ground truth.

<table>
<thead>
  <tr>
    <th>Dataset Role</th>
    <th>Biological Context</th>
    <th>Source / Accession</th>
    <th>Cells (Filtered)</th>
    <th>Subtype</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td rowspan="5">Training Datasets</td>
    <td>NSCLC</td>
    <td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907">GSE131907</a></td>
    <td>73,593</td>
    <td>Follicular B cells<br>Germinal B cells<br>MALT B cells<br>CD4+Tfh/Th cells<br>CD4+exhausted T cells<br>CD4+naive T cells<br>CD4+reg T cells<br>CD8+Naive T cells<br>CD8+exhausted T cells<br>DC<br>cDC<br>Plasma cells</td>
  </tr>

  <tr>
      <td>Breast cancer</td>
      <td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078">GSE176078</a></td>
      <td>83,472</td>
      <td>Memory B cells<br>Naive B cells<br>CD4+Tfh/Th cells<br>CD4+memory T cells<br>CD4+naive T cells<br>CD4+reg T cells<br>CD8+exhausted T cells<br>CD8+memory T cells<br>DC<br>cDC<br>Plasmablasts</td>
  </tr>

  <tr>
      <td>Hepatocellular carcinoma</td>
      <td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140228">GSE140228</a></td>
      <td>17,182</td>
      <td>CD4+Tfh/Th cells<br>CD4+memory T cells<br>CD4+naive T cells<br>CD4+reg T cells<br>CD8+MAIT T cells<br>CD8+Naive T cells<br>CD8+exhausted T cells<br>CD8+memory T cells<br>DC<br>cDC<br>Plasma cells</td>
  </tr>

  <tr>
      <td>Colorectal cancer</td>
      <td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465">GSE132465</a></td>
      <td>33,092</td>
      <td>CD4+Tfh/Th cells<br>cDC<br>IgA+ Plasma<br>IgG+ Plasma</td>
  </tr>

  <tr>
      <td>Prostate cancer</td>
      <td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141445">GSE141445</a></td>
      <td>29,636</td>
      <td></td>
  </tr>

  <tr>
    <td rowspan="7">Test Datasets</td>
    <td>NSCLC</td>
    <td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153935">GSE153935</a></td>
    <td>4,844</td>
    <td></td>
  </tr>

  <tr>
    <td>Melanoma</td>
    <td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978">GSE115978</a></td>
    <td>4,634</td>
    <td></td>
  </tr>

  <tr>
    <td>Melanoma immune cells</td>
    <td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123139">GSE123139</a></td>
    <td>7,903</td>
    <td></td>
  </tr>

  <tr>
    <td>Colorectal cancer</td>
    <td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166555">GSE166555</a></td>
    <td>35,382</td>
    <td></td>
  </tr>

  <tr>
    <td>Sorted PBMC</td>
    <td><a href="https://doi.org/10.5281/zenodo.3357167">Zheng</a></td>
    <td>17,385</td>
    <td></td>
  </tr>

  <tr>
    <td>Pulmonary fibrosis lung</td>
    <td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893">GSE135893</a></td>
    <td>113,002</td>
    <td></td>
  </tr>

  <tr>
    <td>Tabula Sapiens</td>
    <td><a href="https://doi.org/10.6084/m9.figshare.25335436.v2">GSE201333</a><br>Download from scDiffusion</td>
    <td>49,357</td>
    <td></td>
  </tr>

  <tr>
    <td>CAF Dataset</td>
    <td>Pan-cancer CAF atlas</td>
    <td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210347">GSE210347</a></td>
    <td>73,402</td>
    <td>CAF_myo<br>CAF_infla<br>Fibroblasts<br>CAF_adi<br>CAF_EndoMT<br>CAF_PN<br>CAF_api</td>
  </tr>
  </tbody>
</table>

## Preprocessing

The five training sources are first restricted to their shared genes and combined, after which quality control is applied once to the integrated reference. The five standard-QC test datasets and the two CAF splits are processed separately using the same sequence: cells must contain at least 200 detected genes, genes must be detected in at least three cells, observation and gene names are made unique, and genes whose symbols begin with `RP` are removed. After recalculating quality-control metrics, cells are retained only when both total counts and the number of detected genes fall strictly between the dataset-specific 2nd and 98th percentiles and the mitochondrial fraction is below 10%.

After standard quality control, Scrublet is run on the integrated training reference and the five standard-QC test datasets with an expected doublet rate of 0.06. Predicted doublets are removed from the integrated training reference, whereas test-set doublet scores and predictions are retained without removing cells. Scrublet is not applied to the CAF data.

GSE135893 (`lung_full`) and GSE201333 (`sapiens_full`) follow the scDiffusion source-filtering policy instead: cells must contain at least 10 detected genes and genes must be detected in at least three cells. These two datasets do not undergo `RP`-gene removal, percentile filtering, mitochondrial-threshold filtering, or Scrublet. GSE135893 is reconstructed from the original-study count matrix, whereas GSE201333 uses the explicitly stored raw-count layer from the scDiffusion Figshare release.

For all final outputs, source counts are preserved before library-size normalization to 10,000 counts per cell followed by `log1p` transformation.

## Label Mapping

All cell-type labels used for scLightGAT training and testing are mapped from published source annotations; no model-predicted labels are used as ground truth. The table below lists the metadata fields used to map cell types and subtypes. For some datasets, a fallback field is used when the primary field is missing or cannot be mapped; all listed fields come from the same published annotation source.

<table>
<thead>
  <tr>
    <th>Dataset Role</th>
    <th>Source / Accession</th>
    <th>Cell-type Source</th>
    <th>Subtype Source</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td rowspan="5">Training Datasets</td>
    <td>GSE131907</td>
    <td>GEO: <code>obs["Cell_subtype"]</code> (fallback: <code>obs["Cell_type.refined"]</code> / <code>obs["Cell_type"]</code>)</td>
    <td>GEO: <code>obs["Cell_subtype"]</code></td>
  </tr>
  <tr>
    <td>GSE176078</td>
    <td>GEO: <code>obs["celltype_subset"]</code> (fallback: <code>obs["celltype_minor"]</code>)</td>
    <td>GEO: <code>obs["celltype_subset"]</code></td>
  </tr>
  <tr>
    <td>GSE140228</td>
    <td>GEO: <code>obs["celltype_sub"]</code></td>
    <td>GEO: <code>obs["celltype_sub"]</code></td>
  </tr>
  <tr>
    <td>GSE132465</td>
    <td><a href="https://doi.org/10.5281/zenodo.7419236">scATOMIC</a>: <code>obs["harmonized_original_annotation"]</code> (fallback: <code>obs["original_annotation.1"]</code>)</td>
    <td><a href="https://doi.org/10.5281/zenodo.7419236">scATOMIC</a>: <code>obs["original_annotation.2"]</code></td>
  </tr>
  <tr>
    <td>GSE141445</td>
    <td><a href="https://doi.org/10.5281/zenodo.7419236">scATOMIC</a>: <code>obs["harmonized_original_annotation"]</code> (fallback: <code>obs["original_annotation.1"]</code>)</td>
    <td></td>
  </tr>
  <tr>
    <td rowspan="7">Test Datasets</td>
    <td>GSE153935</td>
    <td><a href="https://doi.org/10.5281/zenodo.7400873">Hanley et al.</a>: <code>obs["Cell.type2"]</code> (fallback: <code>obs["Lineage"]</code>)</td>
    <td></td>
  </tr>
  <tr>
    <td>GSE115978</td>
    <td>GEO: <code>obs["cell.types"]</code></td>
    <td></td>
  </tr>
  <tr>
    <td>GSE123139</td>
    <td><a href="https://li-melanoma-2018.s3.eu-west-1.amazonaws.com/Li2018_umi_counts.tar.gz">Li2018</a>: <code>obs["t_nk_mc_group"]</code> with <code>obs["Cell.type"]</code> for shared-naive T cells; <code>obs["non_t_nk_mc_group"]</code> for non-T/NK cells</td>
    <td></td>
  </tr>
  <tr>
    <td>GSE166555</td>
    <td>GEO: <code>obs["main_cell_type"]</code> with <code>obs["cell_type_imm_simple"]</code> / <code>obs["cell_type_str_simple"]</code></td>
    <td></td>
  </tr>
  <tr>
    <td>Zheng</td>
    <td><a href="https://doi.org/10.5281/zenodo.3357167">Abdelaal</a>: <code>Labels.csv</code></td>
    <td></td>
  </tr>
  <tr>
    <td>GSE135893</td>
    <td>GEO: <code>obs["celltype"]</code></td>
    <td></td>
  </tr>
  <tr>
    <td>GSE201333</td>
    <td><a href="https://doi.org/10.6084/m9.figshare.25335436.v2">scDiffusion</a>: <code>obs["free_annotation"]</code></td>
    <td></td>
  </tr>
  <tr>
    <td>CAF Dataset</td>
    <td>GSE210347</td>
    <td>GEO: <code>obs["celltype"]</code> (fibroblast/CAF selection)</td>
    <td>GEO: <code>obs["cluster"]</code> (C1&ndash;C8 mapping)</td>
  </tr>
</tbody>
</table>

The scLightGAT training labels are harmonized to the cell-type hierarchy shown below. For test datasets, annotations are mapped to the most specific supported class within the broad cell-type ontology; a supported coarse source annotation is retained rather than assigned to a more specific class without supporting evidence. For example, `T.cells` remains at that level rather than being mapped to `CD4+T cells` or `CD8+T cells`.

<table>
  <thead>
    <tr>
      <th>Major Type</th>
      <th>Subtype</th>
      <th>Major Type</th>
      <th>Subtype</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td rowspan="5">CD4+T cells</td>
      <td>CD4+Tfh/Th cells</td>
      <td rowspan="4">Plasma cells</td>
      <td>IgA+ Plasma</td>
    </tr>
    <tr>
      <td>CD4+exhausted T cells</td>
      <td>IgG+ Plasma</td>
    </tr>
    <tr>
      <td>CD4+memory T cells</td>
      <td>Plasma cells</td>
    </tr>
    <tr>
      <td>CD4+naive T cells</td>
      <td>Plasmablasts</td>
    </tr>
    <tr>
      <td>CD4+reg T cells</td>
      <td rowspan="2">DC</td>
      <td>DC</td>
    </tr>
    <tr>
      <td rowspan="4">CD8+T cells</td>
      <td>CD8+MAIT T cells</td>
      <td>cDC</td>
    </tr>
    <tr>
      <td>CD8+Naive T cells</td>
      <td>NK cells</td>
      <td></td>
    </tr>
    <tr>
      <td>CD8+exhausted T cells</td>
      <td>pDC</td>
      <td></td>
    </tr>
    <tr>
      <td>CD8+memory T cells</td>
      <td>Macrophages</td>
      <td></td>
    </tr>
    <tr>
      <td rowspan="5">B cells</td>
      <td>Follicular B cells</td>
      <td>Monocytes</td>
      <td></td>
    </tr>
    <tr>
      <td>Germinal B cells</td>
      <td>Epithelial cells</td>
      <td></td>
    </tr>
    <tr>
      <td>MALT B cells</td>
      <td>Endothelial cells</td>
      <td></td>
    </tr>
    <tr>
      <td>Memory B cells</td>
      <td>Fibroblasts</td>
      <td></td>
    </tr>
    <tr>
      <td>Naive B cells</td>
      <td>Mast cells</td>
      <td></td>
    </tr>
  </tbody>
</table>

### CAF Target Mapping

The GSE210347 source clusters are mapped to the seven CAF targets as follows:

| Source Cluster | Output Target |
| --- | --- |
| `C1` | `CAF_myo` |
| `C2` | `CAF_infla` |
| `C3` | `Fibroblasts` |
| `C4` | `CAF_adi` |
| `C5` | `Fibroblasts` |
| `C6` | `CAF_EndoMT` |
| `C7` | `CAF_PN` |
| `C8` | `CAF_api` |

## Installation

The scLightGAT preprocessing pipeline is implemented in Python. Clone the repository and create the Conda environment as follows:

```bash
git clone https://github.com/chenh2lab/scLightGAT
cd scLightGAT/preprocessing
conda create --name scLightGAT_data \
  --override-channels \
  --strict-channel-priority \
  --channel conda-forge \
  --channel bioconda \
  --file requirements.txt
conda activate scLightGAT_data
```

## Usage

Run the pipeline from the `preprocessing` directory to download the required public data, preprocess the datasets, harmonize the annotations, and generate all model-ready `.h5ad` files:

```bash
python3 main.py
```

This is equivalent to `python3 main.py build`. By default, the pipeline processes all datasets, uses seed `0` for Scrublet, stores checkpoints in `./.cache`, and writes the outputs to `./scLightGAT_data`. Downloaded source files are removed after a successful full build; run `python3 main.py build --keep-raw` to retain them.

Outputs are organized under `./scLightGAT_data/Integrated_training`, `./scLightGAT_data/Independent_testing`, `./scLightGAT_data/Subtype`, and `./scLightGAT_data/caf.data`.

Subtype datasets use the fixed donor-held-out assignments in [`subtype_donor_split.json`](subtype_donor_split.json), while the CAF dataset uses the fixed tissue-held-out split defined in [`utils/config.py`](utils/config.py).

## Output Format

The generated AnnData files expose the following core fields:

| AnnData Field | Contents |
| --- | --- |
| `.X` | Library-size-normalized and `log1p`-transformed expression |
| `.layers["counts"]` / `.layers["raw_counts"]` | Preserved pre-normalization source count matrices |
| `.layers["log_transformed"]` | Copy of the normalized and `log1p`-transformed `.X` matrix |
| `.obs["Celltype_training"]` | Model-facing target label: a harmonized major or supported coarse cell type for broad datasets, or one of the seven CAF targets for CAF outputs |
| `.obs["Celltype_subtraining"]` | Source-backed subtype label; present only in the integrated-training and derived subtype outputs |
| `.obs["Data_source"]` | Source dataset identifier |

<!-- 

# scLightGAT 公開來源資料重建與標註說明

本目錄用公開、可追溯的作者釋出資料，重建可直接交給 scLightGAT 的
AnnData 輸入。重建的目的不是複製原始 scLightGAT 資料檔裡來源不明、可能混合
raw count 與 log expression 的 `.X`，而是從作者發布的 count matrix 重新完成：

1. 明確選定論文使用的 cohort；
2. 保留來源 count matrix；
3. 執行一致的 QC 與 `normalize_total(1e4) + log1p`；
4. 將原作者標籤保守地對應到 broad cell-type ontology；
5. 對作者只提供 coarse broad label 的 test cell，原樣保留可追溯的
   `Celltype_training`，不以 marker 人工拆分；
6. 輸出 broad benchmark 所需欄位與 expression layers；
7. 將來源可追溯的 training subtype 在每個 broad task 內做固定、donor-held-out
   的約 80/20 切分。

這裡的「raw counts」是作者在 GEO、Zenodo 或 Figshare 公開的計數矩陣，不一定是
由 FASTQ 在本專案重新比對所得。受控存取的原始 reads 不會被繞過，也不會把
normalized `.X` 猜測成 counts。

論文目標、資料集配置與標籤 ontology 參考本專案中的
[學位論文](../scLightGAT/基於自注意力圖神經網路與輕型梯度提升機的單細胞RNA定序細胞類型註釋之整合框架開發.docx)。
本 README 記錄的是目前程式碼實際執行的重建方法；若論文文字與公開資料的實際
schema 不同，以可稽核的來源欄位與本目錄程式碼為準。

## 一鍵重建

在 `preprocessing` 目錄中執行：

```bash
python3 -m pip install -r requirements.txt
python3 main.py build \
  --cache-dir ./.cache \
  --output ./scLightGAT_data \
  --only all \
  --keep-raw \
  --seed 0
python3 main.py status
```

If integrated training already finished QC/Scrublet but stopped during label
validation, resume from the saved QC checkpoint instead of repeating Scrublet:

```bash
python3 main.py build --only train --from-stage annotated --keep-raw --seed 0
```

This rebuilds only integrated annotations, normalization, and the final training
file. It does not rerun integration QC or Scrublet.

直接執行 `python3 main.py` 等同於使用預設路徑建立全部資料，但不保留下載的 raw
cache。若已經完成下載與 QC，可從 `annotated` checkpoint 重新套用目前的來源標籤
policy 並收尾：

```bash
python3 main.py build --only all --from-stage annotated --keep-raw --seed 0
python3 main.py status
```

只有在 `normalized.h5ad` 已由目前這版 source-label 程式產生、且只是最後寫檔被
中止時，才可保留 cache 並只重做 `written` 階段：

```bash
python3 main.py build --only all --from-stage written --keep-raw --seed 0
python3 main.py status
```

這會直接沿用 normalized checkpoints。若 cache 曾由舊的人工標註版本建立，請勿
使用 `--from-stage normalized` 或 `written`，必須從 `annotated` 重建。
Pipeline 會先寫出 test、extra 與 CAF files，
再執行較耗記憶體的 training integration；training concat 只載入一份 counts，不會把
parsed H5AD 中相同的 `X`、`counts`、`raw_counts` 三份矩陣一起合併。

成功的完整 release 中，五個 training source rows 會停在 `parsed`（因為它們只作為
integration components）；七個 independent/extra rows 與兩個 CAF rows 必須到達
`written`，而 `__TRAIN_INTEGRATION__` 必須到達 `qc -> annotated -> normalized ->
written`。實際 cell/gene 數請使用本次 build 產生的
`.cache/reports/final_outputs.csv`，不要引用舊 scLightGAT 檔案或 Word 中歷史版本的
過濾後數量。

## 輸出檔案

```text
scLightGAT_data/
├── Integrated_training/
│   └── train.h5ad
├── Subtype/
│   ├── B_cell_train.h5ad
│   ├── B_cell_test.h5ad
│   ├── CD4_T_cell_train.h5ad
│   ├── CD4_T_cell_test.h5ad
│   ├── CD8_T_cell_train.h5ad
│   ├── CD8_T_cell_test.h5ad
│   ├── DC_train.h5ad
│   ├── DC_test.h5ad
│   ├── Plasma_cell_train.h5ad
│   ├── Plasma_cell_test.h5ad
│   ├── subtype_donor_assignments.csv
│   ├── subtype_split_summary.csv
│   └── subtype_split_manifest.json
├── Independent_testing/
│   ├── GSE153935.h5ad
│   ├── GSE115978.h5ad
│   ├── GSE123139.h5ad
│   ├── GSE166555.h5ad
│   ├── Zhengsorted.h5ad
│   ├── lung_full.h5ad
│   └── sapiens_full.h5ad
└── caf.data/
    ├── caf_train.h5ad
    └── caf_test.h5ad
```

`Subtype/` 只使用 integrated training 中具有來源 subtype annotation 的 cells。
程式以 `Group_ID = Data_source + ":" + Donor_ID` 作為切分單位，並直接讀取
[`subtype_donor_split.json`](subtype_donor_split.json) 內已凍結的完整 254 筆
`broad task × donor group` train/test assignments。該 broad 內同一 donor 的所有
sample、平台與所有 subtype cells 只能全部在 train 或全部在 test；
不同 broad 是獨立分類 task，所以同一 donor 在不同 broad 可有不同的固定
assignment。正式重建不會重新搜尋或隨機產生切分，也沒有可調整的 subtype
split seed。

固定清單同時列出 train 與 test，不會把新增 donor 默認塞入 train。
若 integrated training 出現 JSON 未登記的 donor、JSON 中的 donor 在當次資料缺失、
同一 broad 內重複／衝突 assignment，或切分後不再滿足每個 subtype 兩側至少
兩位 donor 與各來源兩側均有 donor，preprocessing 會直接報錯，不會退回
cell-level split。

因每位 donor 的 cell 數不同，各 subtype 的 test cell 比例不會也不應被強迫恰為
20%。`subtype_donor_split.json` 是唯一具有決定力的 assignment 設定；
`Subtype/subtype_donor_assignments.csv` 只是每次重建後依 JSON 與實際資料產生的衍生
稽核表，不會被 preprocessing 當作輸入。實際 cell/group 數另記錄於
`subtype_split_summary.csv`、`subtype_split_manifest.json` 及每個 H5AD 的
`.uns["subtype_split"]`。

每個輸出檔的主要契約如下：

| AnnData 欄位 | 內容 |
|---|---|
| `.X` | 每個 cell library-size normalize 到 10,000 後的 `log1p` expression |
| `.layers["log_transformed"]` | 與 `.X` 相同，供 scLightGAT 直接讀取，避免再次 `log1p` |
| `.layers["counts"]` | 作者發布且通過檢查的原始計數值 |
| `.layers["raw_counts"]` | `counts` 的保留副本 |
| `.raw` | QC 後完整基因空間的 normalized/log1p AnnData |
| `.obs["Celltype_training"]` | broad ground truth；training 為 13 類 exact label，test 為 exact label 或該 dataset 明示的來源 coarse label |
| `.obs["Celltype_subtraining"]` | 僅 integrated training／Subtype 輸出使用；由原作者 fine annotation 精確 harmonize 的 subtype，沒有來源依據時為 NA |
| `.obs["eligible_subtype"]` | 是否具有可用的 `Celltype_subtraining` ground truth |
| `.obs["Subtype_source_label"]` | harmonization 前的原作者 subtype/fine label |
| `.obs["Subtype_provenance"]` / `.obs["Subtype_resolution_source"]` | 實際採用的原始 metadata 欄位與解析狀態 |
| `.obs["Data_source"]` | 每個 cell 的資料來源 |
| `.obs["Sample_ID"]` | 原作者定義的 specimen/library sample；只用於 provenance，不作 split unit |
| `.obs["Donor_ID"]` | 同一研究內的 biological individual；同 donor 的多 sample/platform 會合併 |
| `.obs["Group_ID"]` | `Data_source:Donor_ID`；Subtype train/test 的唯一 grouping unit |
| `.obs["Celltype_provenance"]` | exact 或 coarse ground truth 的原作者 metadata 路徑 |
| `.obs["Major_source_label"]` | 對應前的原始 broad/fine label |
| `.obs["eligible_broad"]` | 是否可進入 broad model |

各 test file 保留自己的完整 QC 後基因空間。HVG、DGE、train/test gene
intersection、Harmony、class balancing 與模型特徵工程均留給 scLightGAT，不在此處
提前截成舊資料的固定 gene set。

## 標籤 ontology

### Broad cell type

`Celltype_training` 僅允許下列 13 類：

| 類別 | 類別 | 類別 |
|---|---|---|
| CD4+T cells | CD8+T cells | NK cells |
| pDC | B cells | Plasma cells |
| DC | Macrophages | Monocytes |
| Epithelial cells | Endothelial cells | Fibroblasts |
| Mast cells |  |  |

Integrated training 的 `Celltype_training` 不接受 coarse label；test 的同一欄可保留
`T.cells`、`Mono/Mac`、`Immune` 或 `Mesenchymal`；下游依 dataset-specific policy
判斷相容預測，例如 `T.cells` 接受 CD4 或 CD8。`T/NK cells` 與 `Unknown` 不是評估
ground truth。Malignant/tumor、smooth muscle、PVL/perivascular、
pericyte、glia、osteoclast 與 CD34+ progenitor 等不屬於這 13 類；無充分證據時不會被
硬塞入相近類別。

### CAF targets

CAF 是獨立的 tissue-held-out task。`CAF_train` 與 `CAF_test` 直接使用七個 target：

`CAF_EndoMT`, `CAF_adi`, `CAF_api`, `CAF_myo`, `CAF_PN`, `CAF_infla`,
`Fibroblasts`。

GSE210347 作者 `cluster` 欄位的固定轉換為：

| 原始 cluster | 輸出 target |
|---|---|
| C1 | CAF_myo |
| C2 | CAF_infla |
| C3 | Fibroblasts |
| C4 | CAF_adi |
| C5 | Fibroblasts |
| C6 | CAF_EndoMT |
| C7 | CAF_PN |
| C8 | CAF_api |

CAF task 的 target 只作為獨立分類任務寫入 `Celltype_training`；標籤完全來自作者
cluster mapping。

## Label provenance 的判定順序

每個 cell 依下列順序處理：

1. **原始研究標籤**：使用 GEO/作者 metadata 中明確的 broad 或 fine label。
2. **保留的原始 reference metadata**：GSE132465 與 GSE141445 從 scATOMIC
   Zenodo 封存的 training/reference metadata 取回原研究 annotation。程式禁止使用
   `scATOMIC_pred` 或 `layer_*` prediction 欄位當 ground truth。
3. **明確且 branch-aware 的 mapping**：只在 dataset-specific mapping 中轉換已知
   同義詞；不以跨資料集的模糊 substring 猜測 lineage。
4. **來源 coarse evaluation label**：只有 mapper 明列的 coarse 值才寫入
   test 的 `Celltype_training`；不把它猜成較細的 13 類。
5. **未解決**：保留 NA。沒有可支援的 exact/coarse broad ground truth 的 cell 留在
   `annotated.h5ad` 與 cell-fate report 供稽核，但不進 final model file。

因此目前 model-ready release 的 broad ground truth 是「source-backed + conservative
harmonization」；build 不讀取任何人工標註。最準確的 per-cell 說法應讀取
`Celltype_training` 與 `Celltype_provenance`。

## Training datasets：expression、cohort 與標籤來源

五個 training dataset 先各自下載及配對 metadata，再以共同 gene intersection 合併。
QC 與 Scrublet 在合併後執行一次。

Donor identity 在 source metadata 尚未被 integration 丟棄前標準化：GSE131907 使用
官方 feature summary 的 `Sample -> Patient id`（因此 paired `LUNG_Nxx/Txx` 同組）；
GSE176078 使用 26 個官方完整 CID；GSE140228 使用 cell-info 的 `Sample` 與 `Donor`
（跨 Droplet/Smart-seq2 合併）；GSE132465 使用 `Sample` 與 `Patient`；GSE141445
依 GEO `Sample_description` 公布的 column suffix 1--13 到 patient ID 對照。所有
`Group_ID` 都加上 `Data_source` namespace，避免不同研究恰好同名的 donor 被誤合併。

| Dataset | Expression 與 cohort | Broad ground-truth 證據 | 主要 reference |
|---|---|---|---|
| GSE131907, LUAD/NSCLC | GEO `raw_UMI_matrix`; 只保留 `Sample` 符合 regex `LUNG_[NT]\d+`，即論文的 11 primary tumor + 11 distant-normal lung samples | 優先使用 GEO `Cell_subtype` 推導可追溯 broad label，再回退 `Cell_type.refined` / `Cell_type`；Pericyte/Smooth muscle 不猜成 Fibroblast | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907); [Kim et al., 2020](https://doi.org/10.1038/s41467-020-16164-1) |
| GSE176078, breast cancer | GEO 作者發布的 sparse count matrix；26 個 primary tumors（11 ER+, 5 HER2+, 10 TNBC） | 作者 `celltype_subset`、`celltype_minor`，必要時回退 `celltype_major`；PVL 不強制映射為 fibroblast | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078); [Wu et al., 2021](https://doi.org/10.1038/s41588-021-00911-1) |
| GSE140228, HCC immune atlas | GEO droplet UMI 與 Smart-seq2 `read_counts`; 兩個平台都依 `Tissue_sub` 只保留 tumor | 作者 `celltype_sub`；明確 CD4/CD8、B/plasma、NK、macrophage、monocyte、DC、mast；模糊 `CD4/CD8-*` 與 ILCs 保留 unresolved | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140228); [Zhang et al., 2019](https://doi.org/10.1016/j.cell.2019.10.003) |
| GSE132465, colorectal cancer | GEO `raw_UMI_count_matrix`; 依作者 cell annotation 保留 `-T`/tumor cells | scATOMIC Zenodo preserved source/reference metadata，以完整 cell ID 配對；不使用 `scATOMIC_pred` | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465); [Lee et al., 2020](https://doi.org/10.1038/s41588-020-0636-z); [scATOMIC data](https://doi.org/10.5281/zenodo.7419236) |
| GSE141445, prostate cancer | GEO `data.raw.matrix`; 只配對 preserved metadata 中的 GSE141445/Chen prostate rows | scATOMIC Zenodo 封存的原研究/reference metadata；prostate cancer/basal-intermediate 對應 epithelial | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141445); [Chen et al., 2021](https://doi.org/10.1038/s41556-020-00613-6); [scATOMIC data](https://doi.org/10.5281/zenodo.7419236) |

GSE140228 Smart-seq2 檔名雖為 `read_counts`，實際包含作者發布的非負 fractional read
assignments（例如 0.33、0.5）。本重建把它視為作者發布的 raw measurement 並原值
保留，不進行四捨五入；其餘來源預設要求非負整數 counts。

## Independent testing datasets

| Dataset | Expression 與 cohort | Broad evaluation ground truth | 主要 reference |
|---|---|---|---|
| GSE153935, NSCLC fibroblast study | GEO `TLDS_AllCells` count matrix；原研究整合 metadata 從 Zenodo `TLDS_AllCells_RPCAintegrated.Rdata` 依 cell ID 配對 | 原作者 `Cell.type2`, `Cell.type`, `Lineage`；精確類別直接映射，`T.cells`、`Mono/Mac`、`Immune`、`Mesenchymal` 保留為 coarse truth；不以 marker 補標 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153935); [Hanley et al., 2023](https://doi.org/10.1038/s41467-023-35832-6); [Zenodo metadata](https://doi.org/10.5281/zenodo.7400873) |
| GSE115978, melanoma | GEO `counts.csv` + `cell.annotations.csv` | 原作者 `cell.types` / `cell_type` / `CellType`；coarse T 寫為 `T.cells`，malignant 為 ontology 外；不人工拆 T 或重標 B/plasma | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978); [Jerby-Arnon et al., 2018](https://doi.org/10.1016/j.cell.2018.09.006) |
| GSE123139, melanoma immune cells | GEO per-sample count matrices + Li et al. 作者 per-cell metadata；限定 `type=MARS`, `treatment=N`, `location=(S)C or p(S)C` | `t_nk_mc_group`, `non_t_nk_mc_group`, `all_mc_group`, `Cell.type`；shared naive 沒有可靠 CD4/CD8 時以 `T.cells` 評估；不做 marker review | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123139); [Li et al., 2019](https://doi.org/10.1016/j.cell.2018.11.086); [作者資料包](https://li-melanoma-2018.s3.eu-west-1.amazonaws.com/Li2018_umi_counts.tar.gz) |
| GSE166555, colorectal cancer | GEO patient count matrices + metadata；若 condition/tissue/origin 可識別且有足夠 cells，保留 tumor | `main_cell_type` 搭配 compartment-specific `cell_type_str_simple` / `cell_type_imm_simple`；不以 `sct_cell_type` 或 marker review 覆寫 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166555); [Uhlitz et al., 2021](https://doi.org/10.15252/emmm.202114123) |
| Zheng_sorted, sorted PBMC benchmark | Abdelaal et al. Zenodo `Filtered_DownSampled_SortedPBMC_data` 與 `Labels.csv` | `benchmark_label` 明確映射 B、CD4/CD8、monocyte、NK、DC；CD34+ 不猜成成熟免疫類別 | [benchmark data](https://doi.org/10.5281/zenodo.3357167); [Abdelaal et al., 2019](https://doi.org/10.1186/s13059-019-1795-z); [Zheng et al., 2017](https://doi.org/10.1038/ncomms14049) |

Testing dataset 的 Scrublet score 與 `predicted_doublets` 會保留供分析，不像 training
integration 一樣直接移除 predicted doublets。這是為了維持論文 test-set 評估範圍；
使用者仍可在下游另外做敏感度分析。

## scDiffusion-derived extra benchmarks

| Dataset | Expression 與範圍 | Broad evaluation ground truth | 主要 reference |
|---|---|---|---|
| `lung_full` | 依 scDiffusion human-lung benchmark 的 31-label vocabulary 重建，但 expression 回到 GSE135893 的 raw MatrixMarket counts，metadata 使用 `IPF_metadata.csv`；不是沿用 scDiffusion 檔中可能已 normalized 的 `.X` | GSE135893 作者 `celltype` 直接映射；`T Cells` / `Proliferating T Cells` 原樣寫為 `T.cells` coarse truth，smooth muscle 為 ontology 外 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893); [Habermann et al., 2020](https://doi.org/10.1126/sciadv.aba1972); [scDiffusion](https://doi.org/10.1093/bioinformatics/btae518) |
| `sapiens_full` | scDiffusion 官方 Figshare v2 的 `sapiens` package；以精確的六類 vocabulary 找到目標 H5AD，只接受 `raw_counts` 或 `counts` layer | source label 為五個 organs 中的 fibroblast、macrophage、memory B、naive B、NK、plasma，保守映射為 exact broad truth | [Figshare v2](https://doi.org/10.6084/m9.figshare.25335436.v2); [scDiffusion](https://doi.org/10.1093/bioinformatics/btae518); [Tabula Sapiens](https://doi.org/10.1126/science.abl4896) |

`lung_full` 是「使用 scDiffusion benchmark 的資料範圍/label vocabulary，但從原始研究
counts 重建」；`sapiens_full` 則直接取 scDiffusion 官方資料包中明示的 raw-count
layer。兩者不可簡稱為完全相同的上游處理方式。

## CAF dataset 與 tissue split

`CAF_train` 和 `CAF_test` 都來自 GSE210347 pan-cancer CAF atlas 的作者 count RDS 與
metadata。程式先以作者 `celltype` 選出 fibroblast/CAF cells，再用作者 `tissue` 建立
固定且互斥的 tissue split：

| Split | Tissues |
|---|---|
| CAF_train | Pancreas, Bladder, Ovary, Colon |
| CAF_test | Stomach, Lung, Breast, Thyroid, Prostate, Intrahepatic, Liver |

標籤完全取自作者 `cluster` 的 C1–C8 mapping，未使用人工 marker annotation。
GSE210347 是對 12 個既有研究、10 種 solid cancer 的 pan-cancer reanalysis；因此此處的
「來源」應同時引用 GSE210347/PanCAF 論文，而不是把所有 cells 說成同一次新收案。

主要引用：[GSE210347](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210347)、
[Luo et al., 2022](https://doi.org/10.1038/s41467-022-34395-2)、
[PanCAF code archive](https://doi.org/10.5281/zenodo.7095147)。

## QC 與 normalization

### Main、test 與 CAF

處理順序固定為：

1. cell 至少偵測 200 genes；
2. gene 至少在 3 cells 中出現；
3. observation/variable names 唯一化；
4. 移除所有 gene symbol 以 `RP` 開頭的 genes；
5. 在此矩陣重算 QC metrics；
6. `total_counts` 與 `n_genes_by_counts` 都必須嚴格位於各資料的第 2–98 percentile；
7. mitochondrial fraction `< 10%`。

五個 training sources 合併後執行 Scrublet，expected doublet rate 0.06，predicted
doublets 會移除。五個 independent tests 使用相同 Scrublet 參數計分但保留 cells。
CAF 不執行 test-set Scrublet 分支。

### scDiffusion extra datasets

`lung_full` 與 `sapiens_full` 使用 scDiffusion source policy：cell 至少 10 genes、gene
至少 3 cells，不再套用 thesis percentile filters，以避免把 benchmark 範圍改造成另
一個 cohort。

### Expression integrity

- counts 必須非負、finite；除已明確記錄的 GSE140228 Smart-seq2 外必須為整數。
- 重複 gene rows 先相加，不丟失 counts。
- normalization 前保留 `counts` / `raw_counts`。
- `.X` 只由 counts 經 `normalize_total(1e4) + log1p` 產生。
- 缺少 raw-count layer 的上游 H5AD 直接報錯，不從數值分布猜測 `.X`。

## Audit reports 與 cell 去向

每個 logical dataset 都會產生：

- `.cache/reports/<dataset>_annotation_resolution.csv`：原始 source label、實際 mapping
  欄位、model label 與 cell count。
- `.cache/reports/<dataset>_cell_fate.csv`：每個 parsed cell 是否通過 QC、是否取得
  exact/coarse broad truth、annotation status 與 exclusion reason。
- `.cache/reports/final_outputs.csv`：本次實際寫出的 H5AD、cell/gene 數與 non-null
  label 數。

QC 失敗、training predicted doublet、ontology 外與 broad unresolved 是不同的
exclusion reason，不應在方法說明中合併稱為「低品質 cell」。

## 主要引用

### Dataset 與 annotation resources

1. Kim N. et al. *Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma.* Nature Communications (2020). [DOI](https://doi.org/10.1038/s41467-020-16164-1), [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907).
2. Wu S. Z. et al. *A single-cell and spatially resolved atlas of human breast cancers.* Nature Genetics (2021). [DOI](https://doi.org/10.1038/s41588-021-00911-1), [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078).
3. Zhang Q. et al. *Landscape and Dynamics of Single Immune Cells in Hepatocellular Carcinoma.* Cell (2019). [DOI](https://doi.org/10.1016/j.cell.2019.10.003), [GSE140228](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140228).
4. Lee H. O. et al. *Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer.* Nature Genetics (2020). [DOI](https://doi.org/10.1038/s41588-020-0636-z), [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465).
5. Chen S. et al. *Single-cell analysis reveals transcriptomic remodellings in distinct cell types that contribute to human prostate cancer progression.* Nature Cell Biology (2021). [DOI](https://doi.org/10.1038/s41556-020-00613-6), [GSE141445](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141445).
6. Hanley C. J. et al. *Single-cell analysis reveals prognostic fibroblast subpopulations linked to molecular and immunological subtypes of lung cancer.* Nature Communications (2023). [DOI](https://doi.org/10.1038/s41467-023-35832-6), [GSE153935](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153935), [Zenodo](https://doi.org/10.5281/zenodo.7400873).
7. Jerby-Arnon L. et al. *A Cancer Cell Program Promotes T Cell Exclusion and Resistance to Checkpoint Blockade.* Cell (2018). [DOI](https://doi.org/10.1016/j.cell.2018.09.006), [GSE115978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978).
8. Li H. et al. *Dysfunctional CD8 T Cells Form a Proliferative, Dynamically Regulated Compartment within Human Melanoma.* Cell (2019). [DOI](https://doi.org/10.1016/j.cell.2018.11.086), [GSE123139](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123139).
9. Uhlitz F. et al. *Mitogen-activated protein kinase activity drives cell trajectories in colorectal cancer.* EMBO Molecular Medicine (2021). [DOI](https://doi.org/10.15252/emmm.202114123), [GSE166555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166555).
10. Zheng G. X. Y. et al. *Massively parallel digital transcriptional profiling of single cells.* Nature Communications (2017). [DOI](https://doi.org/10.1038/ncomms14049).
11. Abdelaal T. et al. *A comparison of automatic cell identification methods for single-cell RNA sequencing data.* Genome Biology (2019). [DOI](https://doi.org/10.1186/s13059-019-1795-z), [benchmark archive](https://doi.org/10.5281/zenodo.3357167).
12. Habermann A. C. et al. *Single-cell RNA sequencing reveals profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary fibrosis.* Science Advances (2020). [DOI](https://doi.org/10.1126/sciadv.aba1972), [GSE135893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893).
13. Luo H. et al. *Pan-cancer single-cell analysis reveals the heterogeneity and plasticity of cancer-associated fibroblasts in the tumor microenvironment.* Nature Communications (2022). [DOI](https://doi.org/10.1038/s41467-022-34395-2), [GSE210347](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210347).
14. Tabula Sapiens Consortium. *The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans.* Science (2022). [DOI](https://doi.org/10.1126/science.abl4896).
15. Luo E. et al. *scDiffusion: conditional generation of high-quality single-cell data using diffusion model.* Bioinformatics (2024). [DOI](https://doi.org/10.1093/bioinformatics/btae518), [paper data v2](https://doi.org/10.6084/m9.figshare.25335436.v2).
16. Nofech-Mozes I. et al. *Pan-cancer classification of single cells in the tumour microenvironment.* Nature Communications (2023). [DOI](https://doi.org/10.1038/s41467-023-37353-8), [preserved data/code](https://doi.org/10.5281/zenodo.7419236).

## 可直接用於方法章的簡短敘述

> We reconstructed scLightGAT-compatible AnnData objects from publicly
> released count matrices and source metadata rather than reusing the
> expression matrix of the original model-data release. Source annotations
> were conservatively harmonized to 13 broad cell types. Author-provided
> coarse labels were retained as dataset-specific evaluation ground truth
> without marker-based reassignment; unsupported cells were excluded from the
> model-ready output. Raw counts, normalized log1p expression, per-cell source
> provenance, and exclusion reasons were retained for reproducibility. No
> manual marker annotation was used. Within each broad-type subtype task, one
> predefined donor-held-out assignment was shared by all subtypes; samples,
> platforms, and subtype cells from the same donor never crossed train and test.

引用本段時仍應依實際分析對象一併引用上方對應的原始 dataset paper、資料庫 accession
以及原始資料 references。 -->
