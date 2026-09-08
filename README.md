# supra_mic_growth_2026

This repository is provided exclusively for peer review under manuscript submission number `placeholder`.



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
└── HPC_Pipeline/                           ← HPC Pipelines
```

[1] PAP results obtained from forked diskimageR repository  <a href="https://github.com/tony27786/diskImageR.git" target="_blank" style="display:inline-flex; align-items:center; padding:5px 10px; border-radius:10px; background:#24292f; color:#fff; text-decoration:none; font-weight:600; box-shadow:0 6px 18px rgba(0,0,0,.12);"><img src="data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAxNiAxNiIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik04IDBDMy41OCAwIDAgMy41OCAwIDhjMCAzLjU0IDIuMjkgNi41MyA1LjQ3IDcuNTkuNC4wNy41NS0uMTcuNTUtLjM4IDAtLjE5LS4wMS0uODItLjAxLTEuNDktMi4wMS4zNy0yLjUzLS40OS0yLjY5LS45NC0uMDktLjIzLS40OC0uOTQtLjgyLTEuMTMtLjI4LS4xNS0uNjgtLjUyLS4wMS0uNTMuNjMtLjAxIDEuMDguNTggMS4yMy44Mi43MiAxLjIxIDEuODcuODcgMi4zMy42Ni4wNy0uNTIuMjgtLjg3LjUxLTEuMDctMS43OC0uMi0zLjY0LS44OS0zLjY0LTMuOTUgMC0uODcuMzEtMS41OS44Mi0yLjE1LS4wOC0uMi0uMzYtMS4wMi4wOC0yLjEyIDAgMCAuNjctLjIxIDIuMi44Mi42NC0uMTggMS4zMi0uMjcgMi0uMjcuNjggMCAxLjM2LjA5IDIgLjI3IDEuNTMtMS4wNCAyLjItLjgyIDIuMi0uODIuNDQgMS4xLjE2IDEuOTIuMDggMi4xMi41MS41Ni44MiAxLjI3LjgyIDIuMTUgMCAzLjA3LTEuODcgMy43NS0zLjY1IDMuOTUuMjkuMjUuNTQuNzMuNTQgMS40OCAwIDEuMDctLjAxIDEuOTMtLjAxIDIuMiAwIC4yMS4xNS40Ni41NS4zOUE4LjAxMyA4LjAxMyAwIDAwMTYgOGMwLTQuNDItMy41OC04LTgtOHoiLz48L3N2Zz4=" width="20" height="20" style="margin-right:8px;"><span>diskimageR GitHub Repository</span></a>

# 2. Descriptions

1. The analyses were performed on high-performance computing (HPC) clusters. Detailed system information is provided below:

   ```shell
   OS: CentOS Linux 7 (Core)
   Kernel: Linux 5.4.0-100-generic
   Architecture: x86_64
   CPU: Intel(R) Xeon(R) Gold 6240R CPU @ 2.40GHz
   Scheduler: slurm 22.05.6
   ```

2. HPC pipelines contains two part.