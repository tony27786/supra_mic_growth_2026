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

[1] PAP results obtained from forked diskimageR repository  [![diskImageR GitHub Repository](https://img.shields.io/badge/GitHub-diskImageR-24292f?logo=github&logoColor=white)](https://github.com/tony27786/diskImageR)

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
