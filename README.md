# supra_mic_growth_2026

**This repository accompanies a manuscript currently under peer review.**

This repository holds everything needed to (1) redraw every panel of Figures 1–6
from the deposited source data, and (2) regenerate that source data from the raw
reads in SRA.

# 1. Repository Structure

```text
<repo>/
├── README.md
├── SourceData/
│   ├── KB_Cohort_SourceData.xlsx        [Raw] Experimental data
│   ├── kbpap_SourceData.xlsx            [Raw] Experimental data
│   ├── smg_SourceData.xlsx              [Raw] Experimental data
│   ├── {parent,monoclonal,selected,passage}_pap_feature_results.csv   [Intermediate] PAP Result [1]
│   ├── gene_counts.txt                  [Intermediate] From HPC Pipeline
│   ├── master_wgs_table.tsv             [Intermediate] From HPC Pipeline
│   ├── het_chrom_summary.tsv            [Intermediate] From HPC Pipeline
│   ├── mosdepth_chr_normalized_depth_long.tsv   [Intermediate] From HPC Pipeline
│   ├── trackA_ploidy_summary.tsv        [Intermediate] From HPC Pipeline
│   ├── smallvar_4state_long.tsv         [Intermediate] From HPC Pipeline
│   ├── target_panel_callability.tsv     [Intermediate] From HPC Pipeline
│   ├── tree.core_SNPs.ML.tre            [Intermediate] From HPC Pipeline
│   ├── Cn_all_merged_results.txt        [Intermediate] From HPC Pipeline
│   └── refs/                            [Public] NCBI + FungiDB refs
│       ├── GCF_000149245.1_CNA3_genomic.fna / .gtf
│       └── FungiDB-68_CneoformansH99_GO.gaf.gz / _Curated_GO.gaf.gz
├── Plot_Codes/
│   ├── fig1_submit.R
│   ├── fig2_submit.R
│   ├── fig3_submit.R
│   ├── fig4_submit.R
│   ├── fig5_submit.R
│   └── fig6_submit.R
├── Derived/                               (created by fig3–fig5; intermediate DE and phenotype tables)
├── Figure1_panels/ … Figure6_panels/      (created by the figure scripts)
└── HPC_Pipeline/                           ← See the README.md in this folder for details.
```

[1] PAP results obtained from forked diskimageR repository  [![diskImageR GitHub Repository](https://img.shields.io/badge/GitHub-diskImageR-24292f?logo=github&logoColor=white)](https://github.com/tony27786/diskImageR)

# 2. Quick start to reproduce the figures

```bash
git clone https://github.com/tony27786/supra_mic_growth_2026.git
cd supra_mic_growth_2026
git lfs pull                 # one large table is stored with Git LFS
cd Plot_Codes
# Run in order: fig5 reads DE tables that fig3 and fig4 write to Derived/
for i in 1 2 3 4 5 6; do Rscript fig${i}_submit.R; done
```

# 3. Descriptions

1. The analyses were performed on high-performance computing (HPC) clusters. Detailed system information is provided below:

   ```shell
   OS: CentOS Linux 7 (Core)
   Architecture: x86_64
   CPU: Intel(R) Xeon(R) Gold 6240R CPU @ 2.40GHz
   Scheduler: slurm 22.05.6
   ```

2. `SourceData/smallvar_4state_long.tsv` is ~883 MB and is tracked with [Git LFS](https://git-lfs.com). Install it before cloning, or run `git lfs pull`
   afterwards. Without it the file arrives as a small pointer stub and `fig5_submit.R` will fail while reading it.

3. The disk-diffusion and population-analysis plates were quantified in Fiji/ImageJ driven by a fork of [`diskImageR`](https://github.com/tony27786/diskImageR) that adds six-well plate cropping, a fixed-ROI reader and the PAP wrappers used here. The fork descends from [acgerstein/diskImageR](https://github.com/acgerstein/diskImageR).
   
4. Figure scripts should be run in numerical order. `fig3_submit.R` and `fig4_submit.R` generate differential-expression and phenotype tables in `Derived/` (`fig3_DE_batch01_vs_02.csv`, `fig4_DE_batch03_vs_01.csv` and `fig4_master_out.csv`), which are subsequently read by `fig5_submit.R`. Running `fig5_submit.R` before these upstream scripts will therefore result in a missing-file error.

5. Figures 1, 2 and 6 depend only on files in `SourceData/` and can be generated independently. Each figure script writes its output panels as PDF files to the corresponding `FigureN_panels/` directory in the repository root.
