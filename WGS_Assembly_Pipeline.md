# Whole Genome Assembly & Annotation Pipeline for Wheat
This document contains the complete set of scripts and SLURM jobs used for whole genome assembly and annotation of hexaploid wheat (*Triticum aestivum*, 2n = 6x = 42). Steps are arranged serially from raw read processing through final assembly, repeat masking, gene annotation, and downstream analyses including synteny, collinearity, and genome visualization. For questions, clarifications, or updates, feel free to reach out.

---
## SBATCH Header Reference
```bash
#!/bin/bash
#SBATCH --partition=<general>                # specify the partition/queue
#SBATCH --account=<general_user>             # specify the account
#SBATCH --nodes=1                            # allocate nodes
#SBATCH --cpus-per-task=4                    # allocate CPU cores
#SBATCH --qos=normal                         # specify QoS
#SBATCH --time=1-00:00:00                    # set maximum run time
#SBATCH --job-name="JOB TITLE"               # Job name
#SBATCH --mail-user=your_email@example.com   # Notification email
#SBATCH --mail-type=BEGIN,END,FAIL           # send email notification
#SBATCH --output=job_output_%j.txt           # standard output file
#SBATCH --error=job_error_%j.txt             # error messages file
#SBATCH --export=ALL                         # export all user environment variables
```

---

## Step 00 — Merge FASTQ Files

Concatenates all raw FASTQ chunks into a single merged file.

```bash
cat sample*.fastq.gz > sample.fastq.gz
```
---
---

## Step 01 — GenomeScope 2.0 (v2.0)

Estimates genome size, heterozygosity, and ploidy from k-mer frequency histogram.

```bash
#!/bin/bash
#SBATCH --job-name=Gscope2
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-10:00:00
#SBATCH --output=genomescope_%j.out
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

GS_SCRIPT="/project/<account>/WGS_Data/genomescope2.0/genomescope.R"

# List all FASTQ files
ls *.fastq.gz > fastq_files.txt

# Create temp directory for KMC
mkdir -p kmc_tmp

# Count k-mers (k=21) with KMC
# -ci1: min count; -cs10000: max counter value
kmc -k21 -t16 -m64 -ci1 -cs10000 @fastq_files.txt reads_db kmc_tmp/

# Convert KMC database to histogram
kmc_tools transform reads_db histogram reads.histo -cx10000

# Run GenomeScope 2.0 (-p 6 for hexaploid)
Rscript $GS_SCRIPT -i reads.histo -o genomescope_output -k 21 -p 6

echo "GenomeScope analysis complete. Results are in genomescope_output/"
```

---
---

## Step 02 — FastQC (v0.12.1)

Generates quality control report for raw reads.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --qos=normal
#SBATCH --time=1-00:00:00
#SBATCH --job-name="fastqc_sample"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=fastqc_output.txt

module load fastqc

source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

fastqc sample.fastq.gz -t 4 -o /project/<account>/WGS_Data/new_data/sample/01-fastqc
```

---
---

## Step 03 — FastpLong (v0.3.0)

Quality filters and trims long reads; generates an HTML report.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --qos=normal
#SBATCH --time=1-00:00:00
#SBATCH --job-name="fastplong_sample"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=fastplong_output_%j.Log
#SBATCH --error=fastplong_error_%j.Log

module purge

source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

Report_HTML="/project/<account>/WGS_Data/new_data/sample/011-fastplong/fastplong_sample_report.html"

fastplong -i sample.fastq.gz -o sample_out.fq.gz -h $Report_HTML
```

---
---

## Step 04 — Hifiasm Assembly (v0.25.0)

Assembles HiFi reads into contigs.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=512000
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="sample_hifiasm"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=hifiasm_%j.out
#SBATCH --error=hifiasm_%j.err
#SBATCH --export=ALL

module purge
module load gcc
module load zlib

source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

INPUT_DIR="/project/<account>/WGS_Data/new_data/sample/01-genome"
FASTQ_FILES="${INPUT_DIR}/sample1.fastq.gz ${INPUT_DIR}/sample2.fastq.gz"
OUTPUT_PREFIX="/project/<account>/WGS_Data/new_data/sample/03-hifiasm"

hifiasm -o ${OUTPUT_PREFIX} -t 32 ${FASTQ_FILES}
```

---
---

## Step 05 — BUSCO on Contigs (v6.0.0)

Assesses contig assembly completeness using BUSCO with the Poales lineage.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256000
#SBATCH --qos=normal
#SBATCH --time=3-00:00:00
#SBATCH --job-name="sample_busco"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=busco_%j.out
#SBATCH --error=busco_%j.err
#SBATCH --export=ALL

module purge

INPUT_FILE="/project/<account>/WGS_Data/new_data/sample/03-hifiasm/sample_contig.fa"
LINEAGE_DATASET="poales_odb10"
OUT_NAME="busco_sample_contig_analysis"

busco -i $INPUT_FILE -m genome -l $LINEAGE_DATASET -c $SLURM_CPUS_PER_TASK -o $OUT_NAME
```

---
---

## Step 06 — GFAstats on Contigs (v1.3.11)

Computes assembly statistics (N50, total length, etc.) from the contig FASTA.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=2
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=40000
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="gfastats_sample_contig"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=gfastats_%j.out
#SBATCH --error=gfastats_%j.err
#SBATCH --export=ALL

module purge

EXECUTABLE="/project/<account>/WGS_Data/gfastats/build/bin/gfastats"
INPUT_FILE="/project/<account>/WGS_Data/new_data/sample/03-hifiasm/sample_contig.fa"
OUTPUT="/project/<account>/WGS_Data/new_data/sample/04-contig_quality/02-gfastat/gfastats_sample_contig_summary.txt"

genome_size="16000000000"
chmod +x $EXECUTABLE
set -x

$EXECUTABLE $INPUT_FILE $genome_size --discover-paths > $OUTPUT 2> debug.log
```

---
---

## Step 07 — QUAST on Contigs (v5.3.0)

Evaluates contig assembly quality metrics (NGA50, mismatches, etc.).

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256000
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="quast_sample"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=sample_quast%j.out
#SBATCH --error=sample_quast%j.err
#SBATCH --export=ALL

set -e
set -x

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

ASSEMBLY="/project/<account>/WGS_Data/new_data/sample/04-contig_quality/03-quast/sample_contig.fa"
OUTDIR="/project/<account>/WGS_Data/new_data/sample/04-contig_quality/03-quast"

quast -o "$OUTDIR" --large --threads 16 "$ASSEMBLY"
```

---
---

## Step 08 — Scaffolding with RagTag (v2.1.0)

Orders and orients contigs into scaffolds using a reference genome.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=768000
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="sample_scaffold"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=sample_%j.out
#SBATCH --error=sample_%j.err
#SBATCH --export=ALL

set -e
set -x

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

REF="/project/<account>/WGS_Data/reference/reference.fa"
QUERY="/project/<account>/WGS_Data/new_data/sample/03-hifiasm/sample_contig.fa"
OUT="/project/<account>/WGS_Data/new_data/sample/05-scaffolding"

ragtag.py scaffold "$REF" "$QUERY" -o "$OUT" -t 24 -u

conda deactivate
```

---
---

## Step 09 — BUSCO on Scaffolds (v6.0.0)

Assesses scaffold assembly completeness using BUSCO.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=768000
#SBATCH --qos=normal
#SBATCH --time=3-00:00:00
#SBATCH --job-name="sample_busco"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=busco_sample_%j.out
#SBATCH --error=busco_sample_%j.err
#SBATCH --export=ALL

module purge

INPUT_FILE="/project/<account>/WGS_Data/new_data/sample/05-scaffolding/ragtag.scaffold.fasta"
LINEAGE_DATASET="poales_odb10"
OUT_NAME="busco_sample_scaffold_analysis"

export _JAVA_OPTIONS="-Xmx700g"

busco -i $INPUT_FILE -m genome -l $LINEAGE_DATASET -c $SLURM_CPUS_PER_TASK -o $OUT_NAME
```

---
---

## Step 10 — GFAstats on Scaffolds (v1.3.11)

Computes assembly statistics from the scaffolded FASTA.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=2
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=40000
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="gfastats_scaff"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=gfastats_scaff_%j.out
#SBATCH --error=gfastats_scaff_%j.err
#SBATCH --export=ALL

module purge

EXECUTABLE="/project/<account>/WGS_Data/gfastats/build/bin/gfastats"
INPUT_FILE="/project/<account>/WGS_Data/new_data/sample/05-scaffolding/ragtag.scaffold.fasta"
OUTPUT="/project/<account>/WGS_Data/new_data/sample/06-scaff_quality/02-gfastat/gfastats_scaff_summary.txt"

genome_size="16000000000"
chmod +x $EXECUTABLE
set -x

$EXECUTABLE $INPUT_FILE $genome_size --discover-paths > $OUTPUT 2> debug.log
```

---
---

## Step 11 — QUAST on Scaffolds (v5.3.0)

Evaluates scaffold assembly quality metrics.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256000
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="sample_quast"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=sample_quast%j.out
#SBATCH --error=sample_quast%j.err
#SBATCH --export=ALL

set -e
set -x

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

ASSEMBLY="/project/<account>/WGS_Data/new_data/sample/05-scaffolding/ragtag.scaffold.fasta"
OUTDIR="/project/<account>/WGS_Data/new_data/sample/06-scaff_quality/03-quast"

quast -o "$OUTDIR" --large --threads 16 "$ASSEMBLY"
```

---
---

## Step 12 — RepeatMasker (v4.2.0)

Masks repetitive elements in the scaffold using a custom TE library.

### Prepare TE Library

Download and combine TREP and ClariTeRep databases:

```bash
wget https://github.com/jdaron/CLARI-TE/raw/master/clariTeRep.fna
wget http://botserv2.uzh.ch/kelldata/trep-db/downloads/trep-db_complete_Rel-16.fasta.gz
```

Reformat headers to RepeatMasker `#Class/Family` format using the script below, then combine:

```bash
python reformat_te_headers.py trep-db_complete_Rel-16.fasta clariTeRep.fna > combined_te_lib_clean.fa
```

**`reformat_te_headers.py`** — Normalizes TREP and ClariTeRep FASTA headers to `Name#Class/Family` format required by RepeatMasker:

```python
#!/usr/bin/env python3
"""
Reformat TREP and ClariTeRep headers to RepeatMasker #Class/Family format.
Usage:
    python reformat_te_headers.py trep-db_complete_Rel-16.fasta clariTeRep.fna > combined_te_lib_clean.fa
"""

import sys
import re

PREFIX_MAP = {
    "RLC": "LTR/Copia",    "RLG": "LTR/Gypsy",    "RLB": "LTR/Bel-Pao",
    "RLR": "LTR/Retrovirus","RLX": "LTR/Unknown",
    "RIL": "LINE/L1",       "RII": "LINE/I",        "RIX": "LINE/Unknown",
    "RST": "SINE/tRNA",     "RSX": "SINE/Unknown",  "RXX": "Unknown",
    "DTC": "DNA/CACTA",     "DTH": "DNA/Harbinger", "DTM": "DNA/MuDR",
    "DTT": "DNA/Mariner",   "DTA": "DNA/hAT",       "DTX": "DNA/Unknown",
    "DTP": "DNA/P",         "DTB": "DNA/PiggyBac",
    "DHH": "RC/Helitron",   "DHX": "RC/Helitron",
    "DMX": "DNA/Maverick",  "XXX": "Unknown",
}

def get_class_family(prefix):
    return PREFIX_MAP.get(prefix.upper(), "Unknown")

def reformat_trep(header):
    name = header[1:].split()[0]
    prefix = name[:3]
    cf = get_class_family(prefix)
    return f">{name}#{cf}"

def reformat_clarite(header):
    name = header[1:].split()[0]
    sf_match = re.search(r'sF=([A-Z]{3})', header)
    if sf_match:
        return f">{name}#{get_class_family(sf_match.group(1))}"
    cl_match = re.search(r'classif=([A-Z]{3})', header)
    if cl_match:
        return f">{name}#{get_class_family(cl_match.group(1))}"
    return f">{name}#Unknown"

def process_file(filepath, mode):
    opener = open
    if filepath.endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(filepath, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                print(reformat_trep(line) if mode == "trep" else reformat_clarite(line))
            else:
                print(line)

def main():
    if len(sys.argv) < 3:
        print("Usage: python reformat_te_headers.py trep.fasta clariTeRep.fna > combined_te_lib_clean.fa")
        sys.exit(1)
    process_file(sys.argv[1], "trep")
    process_file(sys.argv[2], "clarite")

if __name__ == "__main__":
    main()
```

### Run RepeatMasker

```bash
#!/bin/bash
#SBATCH --job-name=repeat_masker
#SBATCH --output=repeat_masker_%j.out
#SBATCH --error=repeat_masker_%j.err
#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=150G
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL

module purge
module load gcc perl
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

GENOME="/project/<account>/WGS_Data/new_data/sample/05-scaffolding/ragtag.scaffold.fasta"
OUTPUT_DIR="/project/<account>/WGS_Data/new_data/sample/07-repeatmasker"
TE_LIB="/project/<account>/WGS_Data/new_data/sample/07-repeatmasker/combined_te_lib.fasta"

RepeatMasker \
  -lib ${TE_LIB} \
  -dir ${OUTPUT_DIR} \
  -gff \
  -xsmall \
  -e ncbi \
  -pa $SLURM_NTASKS_PER_NODE \
  ${GENOME}
```
### Summarize RepeatMasker Output (Python script)

Parses RepeatMasker `.out` file and summarizes TE content by superfamily with genome coverage percentages.

```python
#!/usr/bin/env python3
"""
Usage: python summarize_repeats.py ragtag.scaffold.fasta.out
"""

import sys, collections

GENOME_SIZE = 14800150786  # update to match your .tbl total length

def classify(family):
    f = family.strip()
    if f.startswith("LTR/Gypsy"):       return "Gypsy (RLG)"
    if f.startswith("LTR/Copia"):       return "Copia (RLC)"
    if f in ("LTR/mixture","LTR") or f.startswith("LTR/mixture") or f in ("LTR/Unknown","LTR/unknown","pararetrovirus","mixture"): return "Unclassified LTR (RLX)"
    if f.startswith("LINE"):            return "LINE (RIX)"
    if f.startswith("SINE"):            return "SINE (SIX)"
    if f.startswith("TIR/EnSpm_CACTA") or f.startswith("DNA/CACTA"):     return "CACTA (DTC)"
    if f.startswith("TIR/MuDR_Mutator") or f.startswith("DNA/MuDR"):     return "Mutator (DTM)"
    if f.startswith("TIR/PIF_Harbinger") or f.startswith("DNA/Harbinger") or f == "DNA/PIF-Harbinger": return "Harbinger (DTH)"
    if f.startswith("TIR/Tc1_Mariner") or f.startswith("DNA/Mariner"):   return "Mariner (DTT)"
    if f.startswith("TIR/hAT") or f.startswith("DNA/hAT"):               return "hAT (DTA)"
    if f.startswith("DNA/Unknown") or f.startswith("DNA/unknown") or f.startswith("TIR"): return "Unclassified TIR (DTX)"
    if f.startswith("Helitron") or f.startswith("RC/Helitron"):           return "Helitron (DHH)"
    if f == "Simple_repeat":            return "_Simple_repeat"
    if f == "Low_complexity":           return "_Low_complexity"
    if f.startswith("Satellite"):       return "_Satellite"
    if f in ("Small_RNA",) or f.startswith(("rRNA","tRNA","snRNA")): return "_Small_RNA"
    return "Unclassified TE (XXX)"

def main(out_file):
    bp_counts = collections.defaultdict(int)
    elem_counts = collections.defaultdict(int)
    with open(out_file) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith(("SW","score")): continue
            if line.lstrip().startswith("*"): line = line.lstrip()[1:].lstrip()
            parts = line.split()
            if len(parts) < 11: continue
            try:
                begin, end = int(parts[5]), int(parts[6])
            except ValueError:
                continue
            cat = classify(parts[10])
            bp_counts[cat] += end - begin + 1
            elem_counts[cat] += 1

    G = GENOME_SIZE
    def pct(bp): return 100.0 * bp / G

    # Print summary table (see original for full formatted output)
    for cat, bp in sorted(bp_counts.items()):
        print(f"{cat:<40} {bp:>14,}  {pct(bp):>6.2f}%")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python summarize_repeats.py ragtag.scaffold.fasta.out")
        sys.exit(1)
    main(sys.argv[1])
```
---
---

## Step 13 — GMAP CDS Alignment (v2024.11.20)

Aligns reference CDS sequences to the repeat-masked assembly to generate hints for AUGUSTUS.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=500G
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="gmap-sample"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=gmap-sample_%j.out
#SBATCH --error=gmap-sample_%j.err
#SBATCH --export=ALL

set -euo pipefail
set -x

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

GENOME="/project/<account>/WGS_Data/new_data/sample/07-repeatmasker/ragtag.scaffold.fasta.masked"
CDS="/project/<account>/WGS_Data/genome.file/reference_cds.fasta"
OUT="/project/<account>/WGS_Data/new_data/sample/09-gmap"
DB_NAME="sample_masked"

[[ ! -f "$GENOME" || ! -f "$CDS" ]] && { echo "Error: Missing input files!"; exit 1; }
mkdir -p "$OUT"

echo "Building GMAP index..."
gmap_build -D "$OUT" -d "$DB_NAME" -k 14 "$GENOME"

echo "Running GMAP..."
gmapl -D "$OUT" -d "$DB_NAME" \
  -f gff3_gene \
  -n 1 \
  -t 24 \
  -B 5 \
  --min-intronlength=20 \
  --max-intronlength=500000 \
  "$CDS" > "$OUT/sample_gmap.gff3"

conda deactivate
```

Convert GMAP GFF to AUGUSTUS hints format:

```bash
awk 'BEGIN {OFS="\t"} {
  if ($3 == "exon" || $3 == "CDS") {
    hint_type = ($3 == "CDS") ? "CDSpart" : "exonpart"
    match($9, /Parent=([^;]+)/, arr)
    grp = (arr[1] != "") ? arr[1] : "g"NR
    print $1, "E", hint_type, $4, $5, ".", $7, ".", "src=E;grp=" grp
  }
}' sample_gmap.gff3 > sample_gmap.hints.gff3
```

---
---
## Step 14 — Liftoff Gene Transfer (v1.6.3)

Transfers gene annotations from IWGSC reference to the new assembly.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=500G
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="liftoff-sample"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=liftoff-sample_%j.out
#SBATCH --error=liftoff-sample_%j.err
#SBATCH --export=ALL

set -euo pipefail
set -x

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

TARGET_GENOME="/project/<account>/WGS_Data/new_data/sample/05-scaffolding/ragtag.scaffold.fasta"
REF_GENOME="/project/<account>/WGS_Data/genome.file/IWGSC2.1/iwgsc_refseqv2.1_assembly.fa"
REF_GFF="/project/<account>/WGS_Data/genome.file/IWGSC2.1/iwgsc2.1.HC.gff3"
OUT_GFF="/project/<account>/WGS_Data/new_data/sample/08-liftoff/sample_liftoff.gff3"
UNMAPPED_GFF="/project/<account>/WGS_Data/new_data/sample/08-liftoff/sample_liftoff_unmapped.gff3"

liftoff "$TARGET_GENOME" "$REF_GENOME" -g "$REF_GFF" -o "$OUT_GFF" -u "$UNMAPPED_GFF" -p 24

conda deactivate
```

> **Reference data sources:**
> - GFF: https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gff3/triticum_aestivum/
> - Annotations (HC): https://urgi.versailles.inrae.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v2.1/
> - Assembly: https://urgi.versailles.inrae.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v2.1/

---
---

## Step 15 — AUGUSTUS: Split Genome

Splits masked genome into per-chromosome FASTA files for parallel AUGUSTUS runs.

```bash
#!/bin/bash
#SBATCH --job-name=split_genome
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --output=split_genome.out
#SBATCH --error=split_genome.err
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
set -x

GENOME="/project/<account>/WGS_Data/new_data/sample/07-repeatmasker/ragtag.scaffold.fasta.masked"
CHR_LIST="/project/<account>/WGS_Data/new_data/sample/10-augustus/chr_headers.txt"
OUTDIR="/project/<account>/WGS_Data/new_data/sample/10-augustus/split"

cd "$(dirname "$GENOME")"
samtools faidx "$GENOME"

mkdir -p "$OUTDIR"

# Extract named chromosomes
while read -r chr; do
    samtools faidx "$GENOME" "$chr" > "$OUTDIR/${chr//|/_}.fasta"
done < "$CHR_LIST"

# Extract unplaced contigs
cut -f1 "$GENOME.fai" | grep -v -Ff "$CHR_LIST" > "$OUTDIR/unplaced.list"
xargs -a "$OUTDIR/unplaced.list" samtools faidx "$GENOME" > "$OUTDIR/unplaced.fasta"
```
---
---

## Step 16 — AUGUSTUS: Run Gene Prediction (v3.5.0)

Runs AUGUSTUS in array mode (one job per chromosome split) using GMAP hints.

```bash
#!/bin/bash
#SBATCH --job-name=augustus_array
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=7-00:00:00
#SBATCH --array=1-22
#SBATCH --output=augustus_%A_%a.out
#SBATCH --error=augustus_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com

set -euo pipefail
set -x

source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

SPLIT_DIR="/project/<account>/WGS_Data/new_data/sample/10-augustus/split"
HINTS="/project/<account>/WGS_Data/new_data/sample/09-gmap/sample_gmap.hints.gff3"
OUTDIR="/project/<account>/WGS_Data/new_data/sample/10-augustus/augustus_output"

export AUGUSTUS_CONFIG_PATH=/home/<username>/.conda/envs/<env_name>/config

mkdir -p "$OUTDIR"

FASTAS=($(ls "$SPLIT_DIR"/*.fasta | sort))
CHR_FASTA="${FASTAS[$SLURM_ARRAY_TASK_ID-1]}"
CHR_NAME=$(basename "$CHR_FASTA" .fasta)

augustus \
  --species=wheat \
  --hintsfile="$HINTS" \
  --extrinsicCfgFile="$AUGUSTUS_CONFIG_PATH/extrinsic/extrinsic.MPE.cfg" \
  --UTR=off --gff3=on --uniqueGeneId=true \
  --softmasking=1 \
  "$CHR_FASTA" > "$OUTDIR/$CHR_NAME.gff3"

conda deactivate
```
---
---
## Step 17 — Final Gene Annotation Merge (v1.7.0 — AGAT)

Filters AUGUSTUS predictions, removes TEs and redundant models, merges with Liftoff annotations, and outputs a final GFF3.

```bash
#!/bin/bash
#SBATCH --job-name="sample_gene_annotation"
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --output=sample_geneannot_%j.out
#SBATCH --error=sample_geneannot_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com

LIFTOFF_GFF="/project/<account>/WGS_Data/new_data/sample/08-liftoff/sample_liftoff.gff3"
GENOME="/project/<account>/WGS_Data/new_data/sample/07-repeatmasker/ragtag.scaffold.fasta.masked"
AUGUSTUS_GFF="/project/<account>/WGS_Data/new_data/sample/10-augustus/sample_augustus.gff3"
SPROT_FA="uniprot_sprot.fasta"
TE_FA="Viridiplantae_v4.0.fasta"

source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

# 0. Extract raw Augustus proteins
gffread "$AUGUSTUS_GFF" -g "$GENOME" -y sample_aug.raw.prot.faa -S
echo "Raw proteins: $(grep -c '>' sample_aug.raw.prot.faa)"

# 1. Length filter (500–15000 aa)
seqkit seq --min-len 500 sample_aug.raw.prot.faa | seqkit seq --max-len 15000 > aug.step1.faa
cp aug.step1.faa aug.complete.faa
echo "After length filter: $(grep -c '>' aug.complete.faa)"

# 2. Swiss-Prot filter (keep proteins with significant homology)
diamond makedb --in "$SPROT_FA" --db uniprot_sprot.dmnd || true
diamond blastp \
    --query aug.complete.faa --db uniprot_sprot.dmnd \
    --out aug.vs_sprot.tsv \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
    --max-target-seqs 1 --evalue 1e-15 --query-cover 60 --subject-cover 60 \
    --ultra-sensitive --threads 40
awk '$11 < 1e-25 && $12 / ($13 * 3) > 0.38 {print $1}' aug.vs_sprot.tsv | sort -u > keep_sprot.txt
seqkit grep -n -f keep_sprot.txt aug.complete.faa > sample_augustus_clean.prot.faa
echo "After Swiss-Prot filter: $(grep -c '>' sample_augustus_clean.prot.faa)"

# 3. TE filter (remove transposon-like proteins)
diamond makedb --in "$TE_FA" --db rex_te.dmnd || true
diamond blastp \
    --query sample_augustus_clean.prot.faa --db rex_te.dmnd \
    --out aug.vs_rex.tsv \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --max-target-seqs 1 --evalue 1e-25 --query-cover 40 --id 30 --threads 40
cut -f1 aug.vs_rex.tsv | sort -u > te_hits.txt
seqkit grep -v -n -f te_hits.txt sample_augustus_clean.prot.faa > sample_augustus_clean_noTE.prot.faa
echo "After TE filter: $(grep -c '>' sample_augustus_clean_noTE.prot.faa)"

# 4. Keep longest isoform and extract clean gene IDs
agat_sp_keep_longest_isoform.pl --gff "$AUGUSTUS_GFF" -o sample_augustus_clean.longest.gff3
grep ">" sample_augustus_clean_noTE.prot.faa \
    | sed 's/>//' | sed 's/ .*//' \
    | sed 's/\.t[0-9]*$//' | sed 's/\.mRNA.*$//' \
    | sort -u > gene_ids.txt
agat_sp_filter_feature_from_keep_list.pl \
    --gff sample_augustus_clean.longest.gff3 \
    --keep_list gene_ids.txt --attribute ID \
    --output sample_augustus_clean.gff3

# 5. Overlap removal vs Liftoff (keep Augustus genes with <50% overlap)
awk -F'\t' '$3=="gene" {if(match($9,/ID=([^;]+)/,a)) print $1"\t"($4-1)"\t"$5"\t"a[1]"\t0\t"$7}' \
    "$LIFTOFF_GFF" | sort -k1,1 -k2,2n > liftoff_genes.bed
awk -F'\t' '$3=="gene" {if(match($9,/ID=([^;]+)/,a)) print $1"\t"($4-1)"\t"$5"\t"a[1]"\t0\t"$7}' \
    sample_augustus_clean.gff3 | sort -k1,1 -k2,2n > augustus_genes.bed
bedtools intersect -a augustus_genes.bed -b liftoff_genes.bed -v -f 0.5 -r > augustus_novel.bed
awk '{print $4}' augustus_novel.bed > augustus_novel_ids.txt
agat_sp_filter_feature_from_keep_list.pl \
    --gff sample_augustus_clean.gff3 \
    --keep_list augustus_novel_ids.txt \
    --output sample_augustus_novel.gff3
sed -i 's/\ttranscript\t/\tmRNA\t/' sample_augustus_novel.gff3
agat_convert_sp_gxf2gxf.pl --gff sample_augustus_novel.gff3 --output sample_augustus_novel_final.gff3

# 6. Prefix novel Augustus IDs to avoid conflicts
agat_sp_manage_IDs.pl \
    --gff sample_augustus_novel_final.gff3 \
    --prefix "aug_novel_" \
    --output sample_augustus_novel_prefixed.gff3

# 7. Merge Liftoff + novel Augustus (cat only; AGAT for formatting)
cat "$LIFTOFF_GFF" sample_augustus_novel_prefixed.gff3 > sample_combined_raw.gff3
agat_convert_sp_gxf2gxf.pl --gff sample_combined_raw.gff3 --output sample_final_merged.gff3

# 8. Count verification
LIFTOFF_COUNT=$(grep -P "\tgene\t" "$LIFTOFF_GFF" | wc -l)
NOVEL_COUNT=$(grep -P "\tgene\t" sample_augustus_novel_prefixed.gff3 | wc -l)
FINAL_COUNT=$(grep -P "\tgene\t" sample_final_merged.gff3 | wc -l)
echo "Liftoff: $LIFTOFF_COUNT | Novel: $NOVEL_COUNT | Final: $FINAL_COUNT"

# 9. Verify no novel genes were dropped
grep -P "\tgene\t" sample_augustus_novel_prefixed.gff3 \
    | grep -oP 'ID=[^;]+' | sed 's/ID=//' | sort > expected_aug_ids.txt
grep -P "\tgene\t" sample_final_merged.gff3 \
    | grep -oP 'ID=[^;]+' | sed 's/ID=//' | grep "aug_novel_" | sort > actual_aug_ids.txt
comm -23 expected_aug_ids.txt actual_aug_ids.txt > truly_missing.txt
echo "Novel genes missing: $(wc -l < truly_missing.txt)"

# 10. Summary statistics
agat_sp_statistics.pl --gff sample_final_merged.gff3 -o sample_final_merged.stats.txt

echo "Pipeline complete. Final annotation: sample_final_merged.gff3"
conda deactivate
```
---
---

## Step 18 — Extract Sequences with GFFRead (v0.12.9)

Extracts CDS, protein, and transcript sequences from the merged annotation GFF.

```bash
# CDS sequences
gffread /project/<account>/WGS_Data/new_data/sample/11-gene_annotation/sample_final_merged.gff3 \
  -g /project/<account>/WGS_Data/new_data/sample/05-scaffolding/ragtag.scaffold.fasta \
  -x sample_final_merged_cds.fa

# Protein sequences
gffread /project/<account>/WGS_Data/new_data/sample/11-gene_annotation/sample_final_merged.gff3 \
  -g /project/<account>/WGS_Data/new_data/sample/05-scaffolding/ragtag.scaffold.fasta \
  -y sample_final_merged_proteins.fa

# Transcript sequences
gffread /project/<account>/WGS_Data/new_data/sample/11-gene_annotation/sample_final_merged.gff3 \
  -g /project/<account>/WGS_Data/new_data/sample/05-scaffolding/ragtag.scaffold.fasta \
  -w sample_final_merged_transcripts.fa
```

Clean protein FASTA (remove stops, short seqs, non-standard characters):

```bash
awk 'BEGIN{RS=">"; ORS=""} NR>1 {
  header=substr($0,1,index($0,"\n")-1);
  seq=substr($0,index($0,"\n")+1);
  gsub(/[ \t\n]/,"",seq); seq=toupper(seq);
  gsub(/\./,"",seq); gsub(/\*$/,"",seq);
  if (length(seq) >= 10) print ">"header"\n"seq"\n";
}' sample.proteins.fa > sample.proteins.clean.fa
```
---
---

## Step 19 — BLASTp vs UniProt Swiss-Prot (v2.14.1)

Assigns functional descriptions to predicted proteins via BLASTp against Swiss-Prot.

Prepare database first:

```bash
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -dbtype prot -parse_seqids -out uniprot_sprot
```

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256000
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="Sample_blastp"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=blastp_sample_%j.out
#SBATCH --error=blastp_sample_%j.err
#SBATCH --export=ALL

module purge
module load blast-plus/2.14.1

source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

QUERY="/project/<account>/WGS_Data/new_data/sample/12-inter_blastp/sample_final_merged.proteins.clean.fa"
DB="uniprot_sprot"
OUT="/project/<account>/WGS_Data/new_data/sample/12-inter_blastp/02-blastp/sample_combo_blastp.out"

blastp -query "$QUERY" -db "$DB" -out "$OUT" -evalue 1e-5 -outfmt 6 -max_target_seqs 5 -num_threads 16
```
---
---

## Step 20 — InterProScan (v5.75-106.0)

Assigns functional domains, GO terms, and pathway annotations to proteins via InterProScan array job.

Install:

```bash
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/interproscan-5.75-106.0-64-bit.tar.gz
tar -xvzf interproscan-5.75-106.0-64-bit.tar.gz
cd interproscan-5.75-106.0
chmod +x interproscan.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=sample_interproscan
#SBATCH --array=1-5
#SBATCH --output=iprscan_%j.out
#SBATCH --error=iprscan_%j.err
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --cpus-per-task=48
#SBATCH --mem=786G
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
set -x

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

INPUT="/project/<account>/WGS_Data/new_data/sample/12-inter_blastp/01-interproscan/split/split.fasta.part_00${SLURM_ARRAY_TASK_ID}.fa"
OUTDIR="/project/<account>/WGS_Data/new_data/sample/12-inter_blastp/01-interproscan/interproscan_out_${SLURM_ARRAY_TASK_ID}"

mkdir -p "$OUTDIR"

export JAVA_OPTS="-Xmx740g -XX:+UseG1GC -XX:MaxGCPauseMillis=200"

/project/<account>/WGS_Data/interproscan/interproscan-5.75-106.0/interproscan.sh \
  -i "$INPUT" \
  -f tsv,xml \
  -goterms \
  -pa \
  -cpu 16 \
  -dp --tempdir /project/<account>/WGS_Data/tmp \
  -b "$OUTDIR/output"

conda deactivate
```
---
---

## Step 21 — Dotplot (NUCmer per Chromosome) (v4.0.1 — MUMmer)

Generates per-chromosome dotplots between the new assembly and a reference (e.g., Chinese Spring).

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=768000
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="mummer_dotplot"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=sample_%j.out
#SBATCH --error=sample_%j.err
#SBATCH --export=ALL

set -e
set -x

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

CS="/project/<account>/WGS_Data/genome.file/IWGSC2.1/cs.21chr.fa"
MAD="/project/<account>/WGS_Data/new_data/sample/05-scaffolding/sample.21chr.fa"

mkdir -p per_chr nucmer_out

chroms=(Chr1A Chr1B Chr1D Chr2A Chr2B Chr2D Chr3A Chr3B Chr3D \
        Chr4A Chr4B Chr4D Chr5A Chr5B Chr5D Chr6A Chr6B Chr6D \
        Chr7A Chr7B Chr7D)

echo "Extracting chromosomes..."
for chr in "${chroms[@]}"; do
    samtools faidx $CS $chr > per_chr/"${chr}_CS.fa"
    samtools faidx $MAD $chr > per_chr/"${chr}_SAM.fa"
done

echo "Running NUCmer per chromosome..."
for chr in "${chroms[@]}"; do
    ref=per_chr/${chr}_CS.fa
    qry=per_chr/${chr}_SAM.fa
    out=nucmer_out/${chr}
    nucmer --maxmatch -l 100 -c 500 -p $out $ref $qry
    delta-filter -1 ${out}.delta > ${out}.1to1.delta
    show-coords -rcl ${out}.1to1.delta > ${out}.coords
done

for chr in Chr1A Chr1B Chr1D Chr2A Chr2B Chr2D Chr3A Chr3B Chr3D \
           Chr4A Chr4B Chr4D Chr5A Chr5B Chr5D Chr6A Chr6B Chr6D Chr7A Chr7B Chr7D; do
    delta-filter -i 98 -l 15000 ${chr}.delta > ${chr}.filtered.delta
    mummerplot --large --png --layout -p ${chr} ${chr}.filtered.delta
done

conda deactivate
```
---
---

## Step 22 — MUMmer / NUCmer Whole-Genome (v4.0.1)

Runs whole-genome NUCmer alignment for synteny/SyRI input.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=768000
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="mummer_sample"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=sample_%j.out
#SBATCH --error=sample_%j.err
#SBATCH --export=ALL

set -e
set -x

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

REF="/project/<account>/WGS_Data/genome.file/IWGSC2.1/cs.21chr.fa"
QUERY="/project/<account>/WGS_Data/new_data/sample/05-scaffolding/sample.21chr.fa"
OUT="/project/<account>/WGS_Data/new_data/sample/13-synteny/02-syri/"

nucmer --prefix="$OUT/sam_vs_cs" --maxmatch -l 100 -c 500 "$REF" "$QUERY"
delta-filter -i 95 -l 5000 "$OUT/sam_vs_cs.delta" > "$OUT/sam_vs_cs.filtered.delta"
show-coords -THrd sam_vs_cs.filtered.delta > sam_vs_cs.filtered.coords

conda deactivate
```
---
---

## Step 23 — SyRI (v1.7.1) + PlotSR (v1.1.0)

Identifies structural rearrangements (inversions, translocations) and visualizes synteny between assemblies.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=256000
#SBATCH --qos=normal
#SBATCH --time=8-00:00:00
#SBATCH --job-name="syri"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=syri_sample_%j.out
#SBATCH --error=syri_sample_%j.err
#SBATCH --export=ALL

set -e
set -x

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

COORDS="sam_vs_cs.filtered.coords"
REF="/project/<account>/WGS_Data/genome.file/IWGSC2.1/cs.21chr.fa"
QUERY="/project/<account>/WGS_Data/new_data/sample/05-scaffolding/sam.21chr.fa"

syri -c "$COORDS" -r "$REF" -q "$QUERY" -F T -k --nosnp --nc 40

# Filter to events ≥1 kb ref / ≥10 kb query and plot
awk '($3 - $2 >= 2000) && ($8 - $7 >= 2000)' syri.out > syri.filtered2k.out
plotsr --sr syri.filtered2k.out --genomes genomes.txt --nodup -o Output.png

conda deactivate
```
---
---

## Step 24 — MCScanX Synteny & Collinearity (v1.0.0)

Detects collinear gene blocks between the new assembly and a reference proteome using MCScanX.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=640G
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="mcscanx_sample"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=sample_%j.out
#SBATCH --error=sample_%j.err
#SBATCH --export=ALL

set -euo pipefail
set -e
set -x

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

SAMPLE_GFF="/project/<account>/WGS_Data/new_data/sample/11-gene_annotation/sample_final_merged.gff3"
SAMPLE_PEP="/project/<account>/WGS_Data/new_data/sample/12-inter_blastp/sample_final_merged.proteins.clean.fa"
CS_GFF="/project/<account>/WGS_Data/genome.file/IWGSC2.1/iwgsc2.1.HC.gff3"
CS_PEP="/project/<account>/WGS_Data/genome.file/IWGSC2.1/iwgsc2.1.HC.pep.fasta"
OUTDIR="data"
mkdir -p "$OUTDIR"

# Step 1: Extract GFF for MCScanX (chr, gene, start, end with strand handling)
extract_gff () {
    local gff=$1; local outgff=$2; local chr_prefix=$3; local gene_prefix=$4
    awk -v chr_pref="$chr_prefix" -v gene_pref="$gene_prefix" '
    BEGIN {OFS="\t"}
    $3 == "gene" {
        split($9, a, ";");
        for (i in a) {
            if (a[i] ~ /^ID=/) {
                gsub(/^ID=/, "", a[i]); gene_id = a[i];
                sub(/\.[0-9]+$/, "", gene_id); sub(/\.t[0-9]+$/, "", gene_id);
                chr_raw = $1;
                if (chr_raw ~ /^[Cc]hr/) { sub(/^[Cc]hr/, "", chr_raw); sub(/_RagTag$/, "", chr_raw); chr = chr_pref chr_raw; }
                else { chr = chr_pref chr_raw; }
                start=$4; end=$5; strand=$7;
                if (strand == "-") { print chr, gene_pref gene_id, end, start; }
                else               { print chr, gene_pref gene_id, start, end; }
            }
        }
    }' "$gff" > "$outgff"
}
extract_gff "$SAMPLE_GFF" sample.mc.gff "ma" "ma_"
extract_gff "$CS_GFF"    cs.mc.gff    "cs" "cs_"

# Step 2: Combine and sort GFF
cat sample.mc.gff cs.mc.gff > combined.mc.gff
sort -k1,1 -k3,3n combined.mc.gff > "$OUTDIR/ma_cs.gff"

# Step 3: Clean PEP headers
clean_pep () {
    local pep=$1; local out=$2; local prefix=$3
    awk -v pref="$prefix" '/^>/ { sub(/^>/, "", $0); id=$0; sub(/\.[0-9]+$/, "", id); sub(/\.t[0-9]+$/, "", id); print ">" pref id; next } {print}' "$pep" > "$out"
}
clean_pep "$SAMPLE_PEP" sample.cleaned.pep "ma_"
clean_pep "$CS_PEP"     cs.cleaned.pep    "cs_"

# Step 4: Combine PEP, keep longest per ID
cat sample.cleaned.pep cs.cleaned.pep > combined.pep.tmp
awk '/^>/ {
    if (seq_id && length(seq) > maxlen[seq_id]) { seqs[seq_id]=seq; maxlen[seq_id]=length(seq); }
    seq_id=substr($0,2); seq=""; next
} { seq=seq $0 }
END { if (length(seq) > maxlen[seq_id]) { seqs[seq_id]=seq; maxlen[seq_id]=length(seq); }
      for (id in seqs) { print ">" id "\n" seqs[id]; }
}' combined.pep.tmp > combined.pep.clean.fa
awk '/^>/ {print; next} {gsub(/[^ACDEFGHIKLMNPQRSTVWY]/, ""); print}' combined.pep.clean.fa > "$OUTDIR/ma_cs.pep"

# Step 5: BLAST all-vs-all
makeblastdb -in "$OUTDIR/ma_cs.pep" -dbtype prot -parse_seqids -out "$OUTDIR/ma_cs_db"
blastp -query "$OUTDIR/ma_cs.pep" -db "$OUTDIR/ma_cs_db" \
       -out "$OUTDIR/ma_cs.blast" -evalue 1e-5 -outfmt 6 -max_target_seqs 5 -num_threads 16

# Step 6: Filter BLAST to match GFF IDs
cut -f2 "$OUTDIR/ma_cs.gff" | sort -u > gff_ids.txt
awk 'FNR==NR {a[$1]=1; next} a[$1] && a[$2]' gff_ids.txt "$OUTDIR/ma_cs.blast" > "$OUTDIR/ma_cs.filtered.blast"
mv "$OUTDIR/ma_cs.filtered.blast" "$OUTDIR/ma_cs.blast"

# Step 7: Run MCScanX
MCSCANX_PATH="/project/<account>/WGS_Data/MCScanX"
cd "$OUTDIR"
$MCSCANX_PATH/MCScanX ma_cs
echo "Done! Check $OUTDIR/ma_cs.collinearity for results."

conda deactivate
```
The collinearity results from MCScanX (.collinearity and .gff files) were visualized using AccuSyn (https://accusyn.usask.ca/), a web-based tool for interactive synteny visualization.

---
---
## Step 25 — Mashmap (v3.1.3)

Computes approximate whole-genome alignments for large-scale structural comparisons.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=768000
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="mashmap_sample"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=sample_%j.out
#SBATCH --error=sample_%j.err
#SBATCH --export=ALL

set -e
set -x

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

REF="/project/<account>/WGS_Data/genome.file/ref_genome_A.fa"
QUERY="/project/<account>/WGS_Data/genome.file/ref_genome_B.fa"
OUT="/project/<account>/WGS_Data/new_data/sample/14-translocation/02-mashmap/"

mashmap -s 500000 --perc_identity 98 --filter_mode one-to-one -t 12 \
  -o "$OUT/sam_vs_ref_50k.paf" -r "$REF" -q "$QUERY"

# Filter specific chromosomes
awk '$1 == "Chr2A" {split($13, a, ":"); if (a[3] > 0.98) print $0}' sam_vs_ref_50k.paf > Chr2A_98_50k.paf
awk '$1 == "Chr7D" {split($13, a, ":"); if (a[3] > 0.98) print $0}' sam_vs_ref_50k.paf > Chr7D_98_50k.paf

conda deactivate
```
---
---

## Step 26 — BLASTn for Translocation Analysis (v2.16.0)

Aligns specific query sequences against a target chromosome to investigate translocation breakpoints.

```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name=blast_chr_batch
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=blast_chr_batch.out
#SBATCH --error=blast_chr_batch.err
#SBATCH --export=ALL

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

QUERY_FILE="/project/<account>/WGS_Data/new_data/sample/14-translocation/01-blast/query_sequences.fa"
DB_DIR="/project/<account>/WGS_Data/new_data/sample/14-translocation/01-blast/"
DB_NAME="sam.2a_db"
OUTPUT_DIR="/project/<account>/WGS_Data/new_data/sample/14-translocation/01-blast/output"
BASE_NAME=$(basename "$QUERY_FILE" .fasta)

blastn -query "$QUERY_FILE" -db ${DB_DIR}/${DB_NAME} \
    -out "${OUTPUT_DIR}/${BASE_NAME}_vs_chr.txt" \
    -evalue 1e-10 -num_threads 4 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# Create summary table (coordinates in base pairs)
awk -v base="$BASE_NAME" '
    BEGIN { OFS="\t"; print "Query","Subject","Identity(%)","Align_Len","Q_start","Q_end","Chr_Start(bp)","Chr_End(bp)","Evalue","BitScore" }
    { printf "%s\t%s\t%.2f\t%d\t%d\t%d\t%d\t%d\t%e\t%.1f\n", $1,$2,$3,$4,$7,$8,$9,$10,$11,$12 }
' "${OUTPUT_DIR}/${BASE_NAME}_vs_chr7.txt" > "${OUTPUT_DIR}/${BASE_NAME}_summary_bp.txt"
```
---
---

## Step 27 — Circos Plot (vx.x)

Generates a circular genome visualization showing GC content, gene density (main + novel), and TE density across 21 chromosomes.

### 27a — Prepare Data (`prep_data.slurm`)

Runs all data preparation scripts sequentially.

```bash
#!/bin/bash
#SBATCH --job-name=circos_prep
#SBATCH --output=logs/prep_%j.out
#SBATCH --error=logs/prep_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

mkdir -p logs results data

bash 01_karyotype.sh           && echo "01 done" \
&& bash 02_windows.sh          && echo "02 done" \
&& bash 03_gene_density_main.sh    && echo "03 done" \
&& bash 04_gene_density_novel.sh   && echo "04 done" \
&& bash 05_te_density.sh       && echo "05 done" \
&& bash 06_gc_content.sh       && echo "06 done"

echo "All prep complete"
```

### 27b — Individual Data Scripts

**`01_karyotype.sh`** — Extracts chromosome sizes and builds Circos karyotype file:

```bash
#!/bin/bash
FA=/project/<account>/WGS_Data/new_data/sample/05-scaffolding/sam.21chr.fa

grep -P "^Chr[1-7][ABD]\t" ${FA}.fai | cut -f1,2 | sort -k1,1V > chrom_sizes.txt
python3 make_karyotype.py chrom_sizes.txt > karyotype.txt
cat karyotype.txt
```

**`make_karyotype.py`** — Converts chromosome sizes to Circos karyotype format, coloring by subgenome (A/B/D):

```python
# Usage: python3 make_karyotype.py chrom_sizes.txt > karyotype.txt
import sys

colors = {'A': 'subA', 'B': 'subB', 'D': 'subD'}

with open(sys.argv[1]) as f:
    for line in f:
        chrom, size = line.strip().split('\t')
        subgenome = chrom[-1]
        label = chrom.replace('Chr', '')
        color = colors.get(subgenome, 'grey')
        print(f"chr - {chrom} {label} 0 {size} {color}")
```

**`02_windows.sh`** — Creates 5 Mb genomic windows for density calculations:

```bash
#!/bin/bash
bedtools makewindows -g chrom_sizes.txt -w 5000000 > windows_5mb.bed
echo "Windows: $(wc -l < windows_5mb.bed)"
```

**`03_gene_density_main.sh`** — Computes main annotation gene density per 5 Mb window (genes/Mb):

```bash
#!/bin/bash
GFF=/project/<account>/WGS_Data/new_data/sample/11-gene_annotation/sample_final_merged.gff3

grep -v "^#" $GFF | awk '$3=="gene"' \
    | grep -P "^chr[1-7][ABD]_RagTag\t" \
    | sed 's/chr\([1-7][ABD]\)_RagTag/Chr\1/g' \
    | awk '{print $1"\t"($4-1)"\t"$5}' | sort -k1,1 -k2,2n > genes_main.bed

bedtools intersect -a windows_5mb.bed -b genes_main.bed -c > gene_density_main_raw.txt

awk '{window_mb=($3-$2)/1000000; density=$4/window_mb; printf "%s\t%s\t%s\t%.4f\n",$1,$2,$3,density}' \
    gene_density_main_raw.txt > gene_density_main.txt
```

**`04_gene_density_novel.sh`** — Computes novel Augustus gene density per 5 Mb window:

```bash
#!/bin/bash
GFF=/project/<account>/WGS_Data/new_data/sample/11-gene_annotation/sample_augustus_novel_final.gff3

grep -v "^#" $GFF | awk '$3=="gene"' \
    | grep -P "^chr[1-7][ABD]_RagTag\t" \
    | sed 's/chr\([1-7][ABD]\)_RagTag/Chr\1/g' \
    | awk '{print $1"\t"($4-1)"\t"$5}' | sort -k1,1 -k2,2n > genes_novel.bed

bedtools intersect -a windows_5mb.bed -b genes_novel.bed -c > gene_density_novel_raw.txt

awk '{window_mb=($3-$2)/1000000; density=$4/window_mb; printf "%s\t%s\t%s\t%.4f\n",$1,$2,$3,density}' \
    gene_density_novel_raw.txt > gene_density_novel.txt
```

**`05_te_density.sh`** — Computes TE density (fraction 0–1) per 5 Mb window from RepeatMasker output:

```bash
#!/bin/bash
RM_OUT=/project/<account>/WGS_Data/new_data/sample/07-repeatmasker/ragtag.scaffold.fasta.out

awk 'NR>3 && NF>=9 {
    chrom=$5
    if (chrom ~ /^chr[1-7][ABD]_RagTag$/) {
        sub(/_RagTag$/,"",chrom); sub(/^chr/,"Chr",chrom)
        start=$6; end=$7
        if (start > end) { tmp=start; start=end; end=tmp }
        print chrom"\t"start"\t"end
    }
}' $RM_OUT | sort -k1,1 -k2,2n > repeats.bed

bedtools intersect -a windows_5mb.bed -b repeats.bed -wao \
    | awk '{print $1"\t"$2"\t"$3"\t"$7}' \
    | bedtools groupby -i - -g 1,2,3 -c 4 -o sum > te_density_raw.txt

awk '{window_size=$3-$2; fraction=$4/window_size; if(fraction>1) fraction=1; printf "%s\t%s\t%s\t%.4f\n",$1,$2,$3,fraction}' \
    te_density_raw.txt > te_density.txt
```

**`06_gc_content.sh`** — Computes GC content per 5 Mb window using bedtools nuc:

```bash
#!/bin/bash
FA=/project/<account>/WGS_Data/new_data/sample/05-scaffolding/sam.21chr.fa

bedtools nuc -fi $FA -bed windows_5mb.bed \
    | awk 'NR>1 { printf "%s\t%s\t%s\t%.4f\n", $1, $2, $3, $5 }' > gc_content.txt
```

**`rename_chrom.sh`** — Reusable helper to convert `chr1A_RagTag` → `Chr1A` in any file via stdin:

```bash
#!/bin/bash
grep -P "^chr[1-7][ABD]_RagTag\b" \
    | sed 's/chr\([1-7][ABD]\)_RagTag/Chr\1/g'
```

### 27c — Run Circos (`run_circos.slurm`)

Executes Circos with the assembled configuration.

```bash
#!/bin/bash
#SBATCH --job-name=circos_run
#SBATCH --output=logs/circos_%j.out
#SBATCH --error=logs/circos_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

module purge
source /home/<username>/WGS_Data/etc/profile.d/conda.sh
conda activate <env_name>

cd /project/<account>/WGS_Data/new_data/sample/15-circos

circos -conf circos.conf \
    && echo "Success: result/sample_circos.png" \
    || echo "FAILED - check logs"
```

### 27d — Circos Configuration Files

**`colors.conf`**

```
subA = 68,114,196
subB = 112,173,71
subD = 237,125,49
```

**`ideogram.conf`**

```
<ideogram>
  <spacing>
    default = 0.004r
    break   = 0.5r
  </spacing>
  radius           = 0.85r
  thickness        = 35p
  fill             = yes
  stroke_color     = dgrey
  stroke_thickness = 2p
  show_label       = yes
  label_font       = bold
  label_radius     = dims(ideogram,radius) + 0.12r
  label_size       = 48
  label_parallel   = yes
</ideogram>
```

**`ticks.conf`**

```
show_ticks       = yes
show_tick_labels = yes

<ticks>
  radius      = 1r
  color       = black
  thickness   = 2p
  multiplier  = 1e-6
  format      = %d

  <tick>
    spacing    = 100000000
    size       = 15p
    show_label = yes
    label_size = 18p
    format     = %d Mb
  </tick>

  <tick>
    spacing    = 50000000
    size       = 8p
    show_label = no
  </tick>
</ticks>
```

**`circos.conf`** (main configuration):

```
karyotype = karyotype.txt

<colors>
  <<include <circos_env_path>/etc/colors.conf>>
  subA = 68,114,196
  subB = 112,173,71
  subD = 237,125,49
</colors>

<fonts>
serif = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunrm.otf
serif_roman = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunrm.otf
serif_default = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunrm.otf
serif_normal = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunrm.otf
serif_light = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunrm.otf
serif_bold = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbx.otf
serif_italic = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunui.otf
serif_bolditalic = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbi.otf
serif_italicbold = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbi.otf
serif_condensed = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunssdc.otf
serif_condensedbold = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunssdc.otf
serif_boldcondensed = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunssdc.otf
sans_serif = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbmr.ttf
sans = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbmr.ttf
sansserif = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbmr.ttf
light = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbmr.ttf
normal = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbmr.ttf
default = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbmr.ttf
semibold = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbsr.otf
bold = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbbx.otf
italic = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbmo.otf
bolditalic = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbxo.otf
italicbold = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbxo.otf
condensed = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunssdc.otf
condensedbold = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunssdc.otf
boldcondensed = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunssdc.otf
mono = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmuntt.ttf
fixed = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmuntt.ttf
mono_light = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbtl.otf
fixed_light = /home/karol.marlowe/.conda/envs/circos_env/fonts/modern/cmunbtl.otf
symbols = /home/karol.marlowe/.conda/envs/circos_env/fonts/symbols/circossymbols.otf
symbol = /home/karol.marlowe/.conda/envs/circos_env/fonts/symbols/circossymbols.otf
glyph = /home/karol.marlowe/.conda/envs/circos_env/fonts/symbols/circossymbols.otf
wingding = /home/karol.marlowe/.conda/envs/circos_env/fonts/symbols/circossymbols.otf
wingding1 = /home/karol.marlowe/.conda/envs/circos_env/fonts/symbols/circossymbols.otf
wingding2 = /home/karol.marlowe/.conda/envs/circos_env/fonts/symbols/circossymbols.otf
wingding3 = /home/karol.marlowe/.conda/envs/circos_env/fonts/symbols/circossymbols.otf
glyphweb = /home/karol.marlowe/.conda/envs/circos_env/fonts/symbols/circossymbols.otf
webding = /home/karol.marlowe/.conda/envs/circos_env/fonts/symbols/circossymbols.otf
</fonts>

<<include ideogram.conf>>
<<include ticks.conf>>

<image>
  dir    = results
  file   = sample_circos.png
  png    = yes
  svg    = yes
  radius = 1800p
  angle_offset      = -90
  auto_alpha_colors = yes
  background        = white
</image>

chromosomes_units           = 1000000
chromosomes_display_default = yes

# Group by subgenome A → B → D
chromosomes_order = Chr1A,Chr2A,Chr3A,Chr4A,Chr5A,Chr6A,Chr7A,\
                    Chr1B,Chr2B,Chr3B,Chr4B,Chr5B,Chr6B,Chr7B,\
                    Chr1D,Chr2D,Chr3D,Chr4D,Chr5D,Chr6D,Chr7D

<plots>
  # Chromosome size labels
  <plot>
    type         = text
    file         = chrom_size_label.txt
    r0           = 0.90r
    r1           = 0.98r
    label_radius = 0.945r
    label_size   = 28p
    label_font   = bold
    label_color  = black
    label_snuggle    = no
    label_rotate     = no
    label_parallel   = yes
    <backgrounds>
      <background>
        color = vlgrey
      </background>
    </backgrounds>
  </plot>

  # Track 1: GC Content
  <plot>
    type  = line
    file  = gc_content.txt
    r0    = 0.83r
    r1    = 0.89r
    min   = 0.40   # update after checking your range
    max   = 0.52   # update after checking your range
    color = vdblue
    thickness = 2
    <axes>
      <axis>
        spacing   = 0.25r
        color     = dgrey
        thickness = 1
      </axis>
    </axes>
    <backgrounds>
      <background>
        color = vlgrey
      </background>
    </backgrounds>
  </plot>

  # Track 2: Main Gene Density
  <plot>
    type       = histogram
    file       = gene_density_main.txt
    r0         = 0.74r
    r1         = 0.82r
    min        = 0
    max        = 31   # update after checking your max
    fill_color = dgreen
    color      = dgreen
    thickness  = 1
    <axes>
      <axis>
        spacing   = 0.5r
        color     = lgrey
        thickness = 1
      </axis>
    </axes>
    <backgrounds>
      <background>
        color = vlgrey
      </background>
    </backgrounds>
  </plot>

  # Track 3: Novel Gene Density
  <plot>
    type       = histogram
    file       = gene_density_novel.txt
    r0         = 0.63r
    r1         = 0.72r
    min        = 0
    max        = 3    # update after checking your max
    fill_color = vdblue
    color      = vdblue
    thickness  = 1
    extend_bin = no
    <axes>
      <axis>
        spacing   = 0.5r
        color     = lgrey
        thickness = 1
      </axis>
    </axes>
    <backgrounds>
      <background>
        color = vlgrey
      </background>
    </backgrounds>
  </plot>

  # Track 4: TE Density (innermost)
  <plot>
    type       = histogram
    file       = te_density.txt
    r0         = 0.52r
    r1         = 0.62r
    min        = 0
    max        = 1
    fill_color = dorange
    color      = dorange
    thickness  = 1
    <axes>
      <axis>
        spacing   = 0.5r
        color     = lgrey
        thickness = 1
      </axis>
    </axes>
    <backgrounds>
      <background>
        color = vlgrey
      </background>
    </backgrounds>
  </plot>

</plots>

<<include <circos_env_path>/etc/housekeeping.conf>>
```
---
---