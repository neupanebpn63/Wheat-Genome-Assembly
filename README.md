
# Wheat Genome Assembly and Annotation Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24323076.svg)](https://doi.org/10.5281/zenodo.24323076)

This repository documents the computational workflows and scripts used for whole-genome assembly, structural and functional annotation, and comparative genomic analyses of two Pacific Northwest (PNW) hexaploid wheat cultivars: **Madsen**, a soft white winter wheat (*Triticum aestivum* ssp. *aestivum*), and **Pritchett**, a soft white winter club wheat (*Triticum aestivum* ssp. *compactum*).

Both genomes were generated using PacBio HiFi sequencing and processed using the same computational workflow. The pipeline provides a reproducible record of the assembly, annotation, quality assessment, and downstream analyses performed in these studies.

## Computational Workflow

The workflow encompasses the major stages of genome assembly and analysis:

1. **Read processing and quality assessment:** FASTQ processing, k-mer analysis, and read quality evaluation.
2. **Genome assembly and scaffolding:** De novo contig assembly using hifiasm and reference-guided scaffolding using RagTag.
3. **Assembly evaluation:** Assessment of genome completeness and assembly statistics using BUSCO, QUAST, and GFAstats.
4. **Repeat identification and masking:** Characterization and masking of repetitive sequences using RepeatMasker and curated transposable element libraries.
5. **Structural gene annotation:** Transfer of reference gene models using Liftoff, evidence-guided gene prediction using AUGUSTUS, and integration of gene annotations.
6. **Functional annotation:** Protein similarity searches against UniProt Swiss-Prot, domain and Gene Ontology annotation using InterProScan, and integration of functional information into GFF3 files.
7. **Comparative genomics:** Whole-genome alignment, structural rearrangement analysis, synteny, and gene collinearity using MUMmer, SyRI, and MCScanX.
8. **Genome visualization:** Generation of chromosome-scale comparative plots and circular genome visualizations.
9. **Data preparation and submission:** Standardization of sequence identifiers, preparation of genome annotation resources, and submission of sequencing reads and assemblies to the European Nucleotide Archive (ENA).

The documentation includes SLURM job scripts, command-line workflows, custom Python scripts, software versions, and generalized file paths. Example identifiers and paths are used in place of user-specific and project-specific configurations.

### Pipeline Documentation

The complete workflow is available in two formats:

- **[WGS_Assembly_Pipeline.md](WGS_Assembly_Pipeline.md)** — Markdown documentation containing the computational workflow, commands, and scripts.
- **[WGS_Assembly_Pipeline.html](WGS_Assembly_Pipeline.html)** — HTML version for convenient viewing.

## Associated Genome Projects

### 1. Madsen — Soft White Winter Wheat

**Cultivar:** Madsen  
**Subspecies:** *Triticum aestivum* ssp. *aestivum*  
**Genome type:** Hexaploid wheat (AABBDD; 2n = 6x = 42)

Madsen is a Pacific Northwest-adapted soft white winter wheat cultivar. Its pseudomolecule-scale genome assembly provides a genomic resource for studying wheat genome organization, structural variation, functional gene content, and adaptation in regional breeding germplasm.

**Genome and Data Resources**

| Resource | Accession / Link |
|---|---|
| ENA BioProject | [PRJEB121356](https://www.ebi.ac.uk/ena/browser/view/PRJEB121356) |
| ENA BioSample | [ERS30837014](https://www.ebi.ac.uk/ena/browser/view/ERS30837014) |
| USDA-ARS Ag Data Commons | [10.15482/USDA.ADC/32592210.v1](https://doi.org/10.15482/USDA.ADC/32592210.v1) |

Raw PacBio HiFi sequencing reads and the chromosome-scale genome assembly are deposited in ENA. Additional resources, including structural and functional annotation files (GFF3), predicted proteins, transcript and coding sequences (CDS), repeat annotation files, and filtered InterProScan functional annotations, are available through the USDA National Agricultural Library's Ag Data Commons.

**Associated publication:** Neupane et al. (2026), *Journal of Plant Registrations* (accepted; full citation to be updated upon publication).

### 2. Pritchett — Soft White Winter Club Wheat

**Cultivar:** Pritchett  
**Subspecies:** *Triticum aestivum* ssp. *compactum*  
**Genome type:** Hexaploid wheat (AABBDD; 2n = 6x = 42)

Pritchett is a Pacific Northwest-adapted soft white winter club wheat cultivar characterized by its compact spike morphology and desirable end-use quality. Its pseudomolecule-scale genome assembly provides a genomic resource for investigating club wheat genome organization, structural variation, gene content, and the genetic architecture of the compactum phenotype.

**Genome and Data Resources**

ENA accession information, Ag Data Commons resources, and associated publication details will be updated upon public availability.

## Software and Reproducibility

The pipeline documents the software versions, parameters, computational procedures, and scripts used in the genome assembly and analysis workflows. Scripts were executed in a high-performance computing environment using the SLURM workload manager.

Paths, account identifiers, filenames, and submission credentials have been generalized for public documentation. Users should adapt these settings to their computing environment and verify software dependencies and reference-data compatibility before execution.

## Archiving and Citation

A versioned snapshot of this repository is permanently archived in Zenodo:

**Neupane, B. (2026). Wheat Genome Assembly and Annotation Pipeline. Zenodo.**

**DOI:** https://doi.org/10.5281/zenodo.24323076

Please cite the archived software release and the relevant genome assembly publication when using or adapting the workflows provided in this repository.

The GitHub repository may receive future updates, while the archived Zenodo release preserves a fixed version of the workflow for reproducibility.

## Contact

**Bipin Neupane**  
Department of Plant Pathology  
Washington State University, Pullman, WA, USA

GitHub: [@neupanebpn63](https://github.com/neupanebpn63)

Questions, suggestions, and feedback related to the pipeline are welcome through the repository's GitHub Issues.
