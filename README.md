# Macrolide-resistant *Mycoplasma pneumoniae* ST3 subclades in Taiwan

This repository contains the code and processed data needed to reproduce the bioinformatic analyses in:

> **Emergence of Novel Macrolide-Resistant *Mycoplasma pneumoniae* Sequence Type 3 Subclades in Children after COVID-19 Pandemic, Taiwan.**
> *Emerging Infectious Diseases* (2026). <!-- TODO: authors, volume(issue), pages, DOI -->

## Repository structure

```
.
├── data/
│   ├── meta.csv                 # Metadata of all strains
│   ├── P1_type1/                # P1 type 1 strains, mapped to M129
│   │   ├── ref_M129.node_labelled.final_tree.tre
│   │   └── ref_M129.summary_of_snp_distribution.vcf
│   └── P1_type2/                # P1 type 2 strains, mapped to FH
│       ├── p2FH.node_labelled.final_tree.tre
│       └── p2FH.summary_of_snp_distribution.vcf
├── ref/
│   ├── NC_000912.fasta          # M129 reference genome
│   ├── NC_000912.gff3           # M129 genome annotation
│   └── Genetic_Codes.txt        # Genetic code table used for annotation
├── Fig3.R                       # Reproduces Figure 3
├── Fig4.R                       # Reproduces Figure 4
├── snp_stat.R                   # Clade-specific SNP testing -> clade_specific_snp.csv
├── clade_specific_snp.csv
├── variants_annotation.ipynb    # Variant annotation against M129 -> TableS1.csv
└── TableS1.csv
```

## Data

| Path | Description |
|---|---|
| `data/meta.csv` | Metadata of all *M. pneumoniae* strains included in the study. |
| `data/P1_type1/` | Phylogenomic tree and variant profile of all P1 type 1 strains, with reads mapped to the M129 reference genome (NC_000912). |
| `data/P1_type2/` | Phylogenomic tree and variant profile of all P1 type 2 strains, with reads mapped to the FH reference genome. <!-- TODO: FH accession --> |

Strains were analyzed separately by P1 type because mapping P1 type 2 strains to the type 1 reference (M129) produces spurious variants in the divergent P1 adhesin region.

Variants were called with [Snippy](https://github.com/tseemann/snippy) and combined into a core alignment with `snippy-core`. <!-- TODO: confirm the following step and add versions --> The tree and variant distribution files were generated with [Gubbins](https://github.com/nickjcroucher/gubbins).

Raw sequencing reads are available at <!-- TODO: e.g. NCBI BioProject PRJNAxxxxxx -->.

## Analyses

### Figures

- `Fig3.R` reproduces Figure 3 of the article.
- `Fig4.R` reproduces Figure 4 of the article.

### Clade-specific SNPs (`snp_stat.R`)

Identifies ST3 clade-specific single-nucleotide polymorphisms using a two-sided two-proportion test and writes the results to `clade_specific_snp.csv`.

### Variant annotation (`variants_annotation.ipynb`)

Annotates variants against the M129 reference genome (`ref/`) and produces `TableS1.csv` (Supplementary Table 1).

A variant was considered a clade-specific mutation when the alternative allele was present in more than 95% of strains within the target clade. Variants below this threshold were removed manually.

## Requirements

<!-- TODO: list R / Python versions and packages. To collect them:
     grep -hE "library\(|require\(" *.R | sort -u
-->

- R
- Python with Jupyter

## Usage

Clone the repository and run the scripts from the repository root, so that relative paths to `data/` and `ref/` resolve correctly:

```bash
git clone https://github.com/YinchengChen23/REPO.git   # TODO: repo name
cd REPO

Rscript Fig3.R
Rscript Fig4.R
Rscript snp_stat.R
jupyter notebook variants_annotation.ipynb
```

## Versions

| Release | Description |
|---|---|
| `v2.0-published` | Code and data corresponding to the published article. **Use this version to reproduce the results.** |
| `v1.0-submission` | Version at initial submission, kept for reference. |

## Citation

If you use this code or data, please cite:

> <!-- TODO: full citation with DOI -->

## License

See [LICENSE](LICENSE).
