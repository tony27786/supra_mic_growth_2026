# HPC_Pipeline

This directory contains the compute pipeline that turns the sequencing reads deposited in SRA into the intermediate tables in `SourceData/`, which the figure scripts in `Plot_Codes/` then read.

Because several steps are computationally intensive and require substantial memory, regeneration from the raw sequencing data is intended to be performed on high-performance computing (HPC) systems.

# 1. What runs here

```text
SRA (WGS, 80 runs)
  └─ wgs_pipeline/            SLURM, one submission chain
       submit_all.sh
         00_fastp  →  01_map  ─┬─  02_mosdepth
                               └─  03_vcf  →  04_extract_het
       submitted separately, all depend on 01_map or the clean reads:
         05_smallvar_4state · trackA_kmer_ploidy · cn_ksnp4 · cn_mlst
         00_make_mappability_mask   (run once per reference, before 04)
  └─ wgs_aftershell/          R, on the cluster
       merge_batches.sh  →  hetero_cnv.R  →  hetero_qc.R  →  wgs_aftershell.R
         └──────────────────────────────────► SourceData/*.tsv

SRA (RNA-seq, 60 runs)
  └─ rna_pipeline/            SLURM
       cn_rnaseq.slurm        fastp → STAR → featureCounts
         └──────────────────────────────────► SourceData/gene_counts.txt
         
Plate photographs[1]
  └─ image_analysis/          Using our forked diskimageR
       run_pap_analysis.R
         └──────────────────────────────────► SourceData/*_pap_feature_results.csv
```

[1] image_analysis codes will provide later.

**Note: **Reference for every step: NCBI RefSeq **GCF_000149245.1** — *Cryptococcus neoformans* var. *grubii* H99, assembly CNA3. WGS and RNA-seq share this one coordinate system.

# 2. Which script writes which intermediate file

| `SourceData/` file                                           | Written by                                   |
| ------------------------------------------------------------ | -------------------------------------------- |
| `mosdepth_chr_normalized_depth_long.tsv`                     | `wgs_aftershell/hetero_cnv.R`                |
| `het_chrom_summary.tsv`                                      | `wgs_aftershell/wgs_aftershell.R`            |
| `master_wgs_table.tsv`                                       | `wgs_aftershell/wgs_aftershell.R`            |
| `trackA_ploidy_summary.tsv`                                  | `wgs_pipeline/trackA_kmer_ploidy.slurm`      |
| `smallvar_4state_long.tsv`                                   | `wgs_pipeline/05_smallvar_4state.slurm`      |
| `target_panel_callability.tsv`                               | `wgs_pipeline/05_smallvar_4state.slurm`      |
| `tree.core_SNPs.ML.tre`                                      | `wgs_pipeline/cn_ksnp4.slurm` → `ksnp4_out/` |
| `Cn_all_merged_results.txt`                                  | `wgs_pipeline/cn_mlst.slurm`                 |
| `gene_counts.txt`                                            | `rna_pipeline/cn_rnaseq.slurm`               |
| `{parent,monoclonal,selected,passage}_pap_feature_results.csv` | `image_analysis/run_pap_analysis.R`[1]       |

[1] image_analysis codes will provide later.

# 3. Structure in this folder

```text
HPC_Pipeline/
├── README.md                     this file
├── manifest/
│   ├── wgs_samples.tsv           80 WGS runs
│   ├── rna_samples.tsv           60 RNA-seq runs
│   └── cohort_groups.csv         isolate → phenotype group
├── wgs_pipeline/
│   ├── 00_config.sh              all shared settings; sourced by every stage
│   ├── submit_all.sh             chains 00 → 01 → {02, 03 → 04}
│   ├── 00_fastp.slurm            adapter/quality trimming
│   ├── 00_make_mappability_mask.slurm   GenMap unique-region mask (once per reference)
│   ├── 01_map.slurm              bwa mem → sort → index → QC
│   ├── 02_mosdepth.slurm         per-chromosome depth
│   ├── 03_vcf.slurm + Snakefile  per-contig bcftools call (ploidy 2) → concat → filter
│   ├── 04_extract_het.slurm      biallelic heterozygous SNPs inside the mask
│   ├── mask_het_table.sh         apply the mask to an existing heterozygous-SNP table
│   ├── 05_smallvar_4state.slurm  four-state small-variant evidence + callability
│   ├── trackA_kmer_ploidy.slurm  FastK → GenomeScope2 reference-free ploidy
│   ├── cn_ksnp4.slurm            shovill → Kchooser4 → kSNP4 core-SNP tree
│   ├── cn_mlst.slurm             stringMLST, assembly-free
│   ├── fig6_panel25_cds.bed      frozen 25-gene panel, union CDS (BED4, col4 = CNAG)
│   └── fig6_panel25.tsv          the same panel, human readable
├── wgs_aftershell/
│   ├── wgs_analysis_config.R     paths, batch tag, heterozygosity baseline switch
│   ├── merge_batches.sh          combine run directories (see §7)
│   ├── hetero_cnv.R              mosdepth summaries → normalized depth table
│   ├── hetero_qc.R               mapping and depth QC
│   ├── wgs_aftershell.R          heterozygosity, clonality, genome-state calls
│   └── group.csv                 isolate → phenotype group
├── rna_pipeline/
│   ├── cn_rnaseq.slurm           fastp → STAR → featureCounts
│   └── make_sample_tsv.sh        build samples.tsv from a FASTQ directory
├── image_analysis/
│   ├── run_pap_analysis.R        Fiji + diskImageR → PAP feature tables [1]
│   └── run_kb_analysis.R         Fiji + diskImageR → disk-diffusion RAD/FoG [1]
└── environment/
    ├── software_versions.tsv     frozen tool versions
    ├── system_info.txt           kernel, SLURM and conda versions
    ├── 99_audit_software_versions.slurm   the audit that produced them
    └── conda_{lists,history,explicit}/    per-environment package manifests
```

[1] image_analysis codes will provide later.

# 4. Descriptions

1. The codes in `image_analysis` folder will be provided later.
2. Several steps are computationally intensive and require substantial memory, regeneration from the raw sequencing data is intended to be performed on high-performance computing (HPC) systems.
3. Our HPC systems are scheduled using Slurm. Details can be found in the `environment` folder.