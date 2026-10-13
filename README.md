# Inclusive jets in Au+Au collisions at 200 GeV

This repository contains the STAR analysis workflow for inclusive jets in
Au+Au collisions at $\sqrt{s_{\mathrm{NN}}}=200$ GeV. It covers jet tree
production for data and embedding, output merging, QA and efficiency studies,
RooUnfold response/closure tests, data unfolding, and preliminary
$R_{CP}$/$R_{AA}$ checks.

## Repository layout

```text
data/                  STAR maker, steering macro, file lists, and job submission
trees/                 merge utilities and merged ROOT files
analysis/
  efficiencies/        matching, trigger, and purity efficiencies
  histograms/           embedding/data spectra, QA, JES/JER, and comparisons
  prescales/            run matching and prescale preparation
  qa/                   additional QA plots
  unfolding/            RooUnfold response, closure, and data unfolding
  results/              preliminary R_CP checks
  raa/                  preliminary R_AA checks
  pp_reference/         pp reference spectra used by the R_AA macro
```

The main implementation is in
`data/StRoot/StPicoHFJetMaker/StPicoHFJetMaker.{h,cxx}`. Tree production is
configured and launched through `data/StRoot/macros/runPicoHFJetMaker.C`.

## Environment

Tree production uses the STAR software stack and 64-bit architecture. On the
Alma9 submission hosts, enter the STAR SL7 container for compilation and local
tests:

```bash
singularity exec -e \
  -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 \
  /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif csh

setup 64b
```

`cons` and local `root4star` tests must run inside this environment. Batch
submission must run outside the container.

FastJet is normally supplied by the STAR environment. Its installation can be
checked with `fastjet-config --version` or by inspecting `LD_LIBRARY_PATH`. If a
private installation is needed, see `data/fastjet_install.sh` and update the
environment variables near the top of `data/runLocal.sh`.

The unfolding scripts use Apptainer and
`analysis/unfolding/roounfold.sif`, which provides ROOT 6 and RooUnfold.

## End-to-end workflow

Run commands from the directory shown in each section. Most wrapper scripts
accept optional input/output arguments; invoking them without arguments uses
the defaults documented below.

### 1. Produce trees

#### Local test

From inside the STAR container:

```bash
cd data
./runLocal.sh             # data and embedding
./runLocal.sh data        # data only
./runLocal.sh embedding   # embedding only
```

The script compiles with `cons` and processes the 10-file test lists in
`data/filelists/`. The resulting test ROOT files are written in `data/`.

#### Batch production

From a STAR submission host, outside the container:

```bash
cd data
./run.sh             # submit data and all embedding pThat bins
./run.sh data        # submit data only
./run.sh embedding   # submit embedding only
```

Embedding is submitted in 11 $\hat{p}_{T}$ bins, from `pt3_5` through
`pt50_-1`. Data uses `filelists/pico_low_14.list`. Jobs are submitted with
`star-submit-beta` and written under
`data/submit/<YYYY-MM-DD>/job_<file-list>/`.

Submission defaults, including the STAR version and number of files per job,
are set in `data/submit/submit.sh`. Physics cuts, triggers, jet radii, matching
settings, cross-section weights, and systematic switches are set in
`data/StRoot/macros/runPicoHFJetMaker.C`.

### 2. Merge production outputs

From `trees/`:

```bash
cd trees
./merge_data_all.sh
./merge_embedding_all.sh
./merge_pThatbins.sh
```

The scripts select the newest directory under `data/submit/` by default and
produce:

- `data_merged.root`
- `embedding_merged.root`
- `merged_all/embedding_pt*.root` (one file per $\hat{p}_{T}$ bin)

To use a specific production or output location:

```bash
./merge_data_all.sh ../data/submit/2025-11-28 /path/to/output data_custom.root
./merge_embedding_all.sh ../data/submit/2025-11-28 /path/to/output embedding_custom.root
./merge_pThatbins.sh ../data/submit/2025-11-28 /path/to/output
```

The environment variables `MIN_SIZE`, `BATCH`, and `CLEAN_SHARDS` control the
minimum accepted file size, number of files per merge shard, and temporary
shard cleanup.

### 3. Make histograms and QA plots

From `analysis/histograms/`:

```bash
cd analysis/histograms
./run_hists.sh                         # embedding_merged.root -> hists.root
./run_make_hists_data.sh               # data_merged.root -> hists_data.root
./run_make_qa_hists.sh                 # PDFs under pdf/QA and pdf/Jets
```

Input and output basenames can be overridden:

```bash
./run_hists.sh embedding_custom.root hists_custom.root
./run_make_hists_data.sh data_custom.root hists_data_custom.root
./run_make_qa_hists.sh data_custom.root /path/to/qa-output
```

Additional wrappers in this directory provide $\hat{p}_{T}$-bin diagnostics,
JES/JER studies, raw-calibration and tree-spectrum comparisons, leading-track
comparisons, and analysis-setting comparisons. Each script prints or documents
its accepted arguments; for example:

```bash
./run_make_hists_pthat_diagnostics.sh
./run_make_jes_jer.sh
./run_draw_jes_residuals.sh
./run_compare_raw_calibrations.sh
./run_compare_tree_spectra.sh
```

### 4. Calculate efficiencies

The efficiency macro derives matching efficiency versus truth jet $p_T$,
trigger efficiency versus reconstructed corrected $p_T$, and reconstruction
purity. Results are produced for radii 0.2, 0.3, and 0.4; the three centrality
classes; and leading-track thresholds of 0, 5, 7, and 9 GeV/$c$.

From `analysis/efficiencies/`:

```bash
cd analysis/efficiencies
./run_efficiencies.sh embedding_merged.root efficiencies.root
```

The input basename is resolved under `trees/`, while the output is written in
`analysis/efficiencies/`. With no arguments, the current wrapper expects
`trees/embedding_merged_MCReco2p0.root` and writes `efficiencies.root`.

To compare efficiency files against a reference $\hat{p}_{T}$ bin and create
ratio PDFs:

```bash
./run_make_efficiency_ratio_pdfs.sh [input-directory] [output-directory] [reference-pattern]
```

The defaults are the current directory, `PDFs_Ratios/`, and reference pattern
`50_-1`. `make_efficiency_overlays.C` can also be run directly with ROOT to
draw overlays from a single efficiency file.

### 5. Build responses and run closure tests

From `analysis/unfolding/`:

```bash
cd analysis/unfolding
./run_unfold_embedding.sh
```

The default command uses Bayesian unfolding, reads
`trees/embedding_merged.root`, and writes to `out_embedding_BAYES/`. The main
ROOT output is `responses_embedding.root`; diagnostic closure plots are written
alongside it. The embedding sample is split deterministically into training and
test subsets for closure.

The method, input, and output directory can be selected explicitly:

```bash
./run_unfold_embedding.sh BAYES embedding_merged.root out_embedding_BAYES
./run_unfold_embedding.sh SVD /absolute/path/embedding.root out_embedding_SVD
```

The response includes matching inefficiency and fake/purity information. The
analysis is evaluated for all configured radii, centralities, and leading-track
cuts.

### 6. Unfold data

Data unfolding combines high-tower and minimum-bias inputs, the embedding
response, equivalent-minimum-bias normalization, and a configurable number of
Bayesian iterations:

```bash
cd analysis/unfolding
./run_unfold_data.sh \
  data_merged_OrigCalib.root \
  dataMB_merged.root \
  out_embedding_BAYES_OrigCalib/responses_embedding.root \
  ../efficiencies/efficiencies.root \
  out_data \
  4
```

Relative data basenames are resolved under `trees/`. The response and
efficiency arguments may be relative to the current directory or absolute as
shown by the wrapper defaults. The output directory contains
`unfolded_data.root` and diagnostic plots.

The current `unfold_data_wEff.cxx` keeps the efficiency-file argument for
backward compatibility, but does not apply it separately: matching efficiency
and fake/purity information are already encoded in the response.

Each input directory must contain `unfolded_data.root`.

### 7. Preliminary nuclear-modification checks

After unfolding, the quick $R_{CP}$ check can be run with:

```bash
cd analysis/results
./run_rcp.sh ../unfolding/out_data/unfolded_data.root
```

By default it writes `out_rcp/rcp.root` and plots under `out_rcp/pdf/`. An
optional third argument selects the unfolded histogram name.

The preliminary $R_{AA}$ wrapper currently uses the paths configured near the
top of `analysis/raa/run_raa_check.sh`, including the unfolded Au+Au result and
the references in `analysis/pp_reference/`:

```bash
./analysis/raa/run_raa_check.sh
```

Outputs are written under `analysis/raa/out_raa/`.

## Tree organization and branches

Each production ROOT file stores one `JetTree` for every radius and centrality
combination:

```text
R0.2/
  CENT_0_10/JetTree
  MID_20_40/JetTree
  PERI_60_80/JetTree
R0.3/...
R0.4/...
```

All trees contain:

| Branch | Description |
| --- | --- |
| `runId`, `eventId` | Run and event identifiers |
| `centralityWeight` | Event centrality weight |
| `reco_pt` | Raw reconstructed jet transverse momentum |
| `reco_pt_corr` | Background-corrected reconstructed jet transverse momentum |
| `reco_eta`, `reco_phi` | Reconstructed jet direction |
| `reco_area`, `reco_rho` | Jet area and event background density |
| `reco_pt_lead` | Leading reconstructed constituent transverse momentum |
| `reco_n_constituents` | Number of reconstructed constituents |
| `reco_neutral_fraction` | Reconstructed neutral-energy fraction |
| `reco_trigger_match` | High-tower trigger-match flag |

Embedding trees additionally contain:

| Branch | Description |
| --- | --- |
| `xsecWeight` | Cross-section weight identifying the $\hat{p}_{T}$ bin |
| `deltaR` | Reconstructed-to-truth jet matching distance |
| `mc_pt`, `mc_eta`, `mc_phi`, `mc_area` | Truth-jet kinematics and area |
| `mc_pt_lead` | Leading truth constituent transverse momentum |
| `mc_n_constituents` | Number of truth constituents |
| `mc_neutral_fraction` | Truth-jet neutral fraction |
| `mc_sum_pt` | Sum of generated-particle transverse momentum used for diagnostics |

Centrality is encoded by the directory containing the tree rather than by a
`centrality` branch.

## Analysis conventions

The current efficiency and unfolding macros share the following defaults:

- jet radii: $R=0.2, 0.3, 0.4$;
- centralities: 0--10%, 20--40%, and 60--80%;
- leading-track thresholds: 0, 5, 7, and 9 GeV/$c$;
- minimum jet areas: 0.07, 0.20, and 0.40 for the three radii;
- maximum reconstructed neutral fraction: 0.95;
- matching requirement in the efficiency macro: $\Delta R < 0.6R$.

These values are repeated in several ROOT macros. When changing the analysis
configuration, keep tree production, histogramming, efficiencies, and
unfolding synchronized.

## Notes

- Large ROOT outputs, generated plots, submission directories, and container
  images are workflow products and should not generally be committed.
- The scripts assume the STAR/BNL filesystem layout and are not intended as a
  portable standalone analysis environment.
- Check script defaults before a production run: several filenames encode a
  particular calibration or reconstruction variation.
