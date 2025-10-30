# Whole Genome Assembly of Wheat Cultivar “Madsen” and “Pritchett”

> Sequencing: PacBio Revio (UW PacBio Sequencing Services)
> Modules: bamtools, blast-plus, busco, fastqc, gfastats, pigz, python

## Steps with SLURM Script

1. Fastqc
For assessing the quality of high-throughput sequencing data. It is typically performed during the initial steps of a sequencing pipeline to ensure the data is of sufficient quality for downstream analyses, including genome assembly and variant calling.
```bash
#!/bin/bash
#SBATCH --partition=<general>			      specify the partition/queue
#SBATCH --account=<general_user>		      specify the account
#SBATCH --nodes=1 				     	      use a single node
#SBATCH --cpus-per-task=4			     	  allocate 24 CPU cores
#SBATCH --qos=normal			    	      specify QoS
#SBATCH --time=1-00:00:00				      set maximum run time
#SBATCH --job-name=”JOB TITLE”	    	      Job name
#SBATCH --mail-user=your_email@example.com	  Notification email
#SBATCH --mail-type=BEGIN,END,FAIL		      send email notification
#SBATCH --output=fastqc_output_%j.txt	      Specifies the file for the standard output
#SBATCH --error=fastqc_error_%j.txt		      Specifies the file for error messages
#SBATCH –output=fastqc_output.txt			  Export all user environment variables

#Load require modules
module load fastqc

# Run Fastqc 
fastqc input_file.fastq -t 4 -o /output/directory
```
(where, t = number of threads)

2. HiFiAdapterFilt
```bash
git clone https://github.com/sheinasim/HiFiAdapterFilt.git
cd HiFiAdapterFilt
chmod +x hifiadapterfiltFCS.sh  (making the script executable)
./hifiadapterfiltFCS.sh
```
```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --qos=normal
#SBATCH --time=1-00:00:00
#SBATCH --job-name=”JOB TITLE”
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL

#Load require modules
module load bamtools/2.5.2
module load blast-plus/2.14.1
module load pigz/2.8

#Activate the Conda environment
source /home/your_username/.conda/envs/your_environment/bin/activate

# Change to the project directory containing HiFiAdapterFilt
cd /project/directory/HiFiAdapterFilt/

# Run HiFiAdapterFilt with specified options
bash hifiadapterfilt.sh -p input_file -l 50 -m 95 -t 4 -o /project/directory/filtered_outputs
```
where -p specifies the prefix for the output files, -l 50 sets the minimum read length threshold, -m 95 sets the maximum mismatch rate threshold, -t allocate the 24 threads, and -o specifies the directory for filtered output files

3. Merging Fastq files, bamfiles, and converting bam files to fastq files


4. Hifiasm
```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name=”JOB TITLE”
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=hifiasm1_%j.out
#SBATCH -- error=hifiasm1_%j.err
#SBATCH --export=ALL

#Load required
module purge
module load gcc/14.1.0
module load zlib/1.3.1
```

5. gfastats
```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name=”JOB TITLE”
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=4
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=40000
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="JOB TITLE"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=gfastats_%j.out
#SBATCH --error=gfastats_%j.err
#SBATCH --export=ALL
 
module purge
 
EXECUTABLE="/path/to/gfastats/build/bin/gfastats"
INPUT_FILE="/path/to/assembled.bp.p_ctg.gfa"
OUTPUT="/path/to/gfastatsmd1_summary.txt"
genome_size="15000000000"

chmod +x $EXECUTABLE
set -x
$EXECUTABLE $INPUT_FILE $genome_size --discover-paths > $OUTPUT 2> debug.log
```

6. BUSCO
```bash
#Converting .gfa files from hifiasm into .fa files
awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa  

#cloning:
git clone https://gitlab.com/ezlab/busco.git
cd busco/

#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256000
#SBATCH --qos=normal
#SBATCH --time=3-00:00:00
#SBATCH --job-name="JOB TITLE"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=busco_pt2_%j.out
#SBATCH --error=busco_pt2_%j.err
#SBATCH --export=ALL
 
module purge

INPUT_FILE="/path/to/filename.fa"
LINEAGE_DATASET="poales_odb10"
OUT_NAME="busco_analysis"

export _JAVA_OPTIONS="-Xmx64G" (Used if the Java memory gets less and detailed results is not published else not needed) 
busco -i $INPUT_FILE -m genome -l $LINEAGE_DATASET -c $SLURM_CPUS_PER_TASK -o OUT_NAME
```

6.1. Busco_Visualize results
```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name=”JOB TITLE”
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

module purge
module load python/3.12.5
module load r/4.4.0

python -m venv busco_env
source busco_env/bin/activate

pip install .
mkdir -p BUSCO_filename

cp /path/to/short_summary.specific.poales_odb10.OUT_NAME.txt BUSCO_filename
export PYTHONPATH=/path/to/busco:$PYTHONPATH

python3 scripts/generate_plot.py –working_directory BUSCO_filename
deactivate
```

7. RagTag
```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=768000
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="JOB TITLE"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=scaffold_%j.out
#SBATCH --error=scaffold_%j.err
#SBATCH --export=ALL

set -e
set -x

module purge
source/home/user.name/etc/profile.d/conda.sh
conda activate <environment>

REF="/path/to/attraction.fasta "
QUERY="/path/to/example.fasta"
OUT="/path/to/ragtag_example"

ragtag.py scaffold “$REF” “$QUERY” -o “$OUT” -t -24 -u
conda deactivate
```

8. Quast
```bash

```

9. RepeatMasker
```bash
conda create -n repeatmask_env_new python=3.8 -y
conda activate repeatmask_env_new

conda install -y -c bioconda -c conda-forge repeatmodeler=2.0.3 repeatmasker rmblast=2.11.0 trf=4.09 recon=1.08 repeatscout=1.0.6 cd-hit=4.8.1 mafft=7.487 ltr_retriever=2.9.0 perl=5.32.1 perl-json perl-file-which=1.24 perl-uri

mkdir -p ~/repeat-annotation
cd ~/repeat-annotation

REPEATMASKER_DIR=$(find "$CONDA_PREFIX" -type d -name "RepeatMasker" | head -n 1)
cd "$REPEATMASKER_DIR"
perl ./configure -rmblast_dir "$CONDA_PREFIX/bin" -trf_prgm "$CONDA_PREFIX/bin/trf" -default_search_engine rmblast

wget https://github.com/jdaron/CLARI-TE/raw/master/clariTeRep.fna
wget http://botserv2.uzh.ch/kelldata/trep-db/downloads/trep-db_complete_Rel-16.fasta.gz
gunzip/project/wrsggl/WGS_Marlowe/md/07-Gene annotation/trep-db_complete_Rel-16.fasta.gz

#!/bin/bash
#SBATCH --job-name=<JOB TITLE>
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

# Load required modules
module load gcc perl RepeatMasker

# Set paths
GENOME="/path/to/ragtag.scaffold.fasta"
OUTPUT_DIR="/path/to/output"
TE_LIB="${OUTPUT_DIR}/combined_te_lib.fasta"

# Combine TE libraries
cat ${OUTPUT_DIR}/clari_te.fasta ${OUTPUT_DIR}/trep-db_complete_Rel-16.fasta > ${TE_LIB}

# Run RepeatMasker with custom library
RepeatMasker \
  -lib ${TE_LIB} \
  -dir ${OUTPUT_DIR} \
  -gff \
  -xsmall \
  -e ncbi \
  -pa $SLURM_NTASKS_PER_NODE \
  ${GENOME}

# Rename output files
mv ${OUTPUT_DIR}/ragtag.scaffold.fasta.masked ${OUTPUT_DIR}/example_masked.fasta
mv ${OUTPUT_DIR}/ragtag.scaffold.fasta.out.gff ${OUTPUT_DIR}/example_te.gff
```
It got failed because the headers was greater than 50 characters. So do the following steps:

```bash
cp /path/to/ragtag.scaffold.fasta example_ragtag_original.fasta

#Generate shortened FASTA with header mapping (using seqkit)
awk '/^>/ {
  old=$0;
  gsub("^>", "", old);
  new=">seq" ++i;
  print new >> "example_ragtag_short.fasta";
  print old "\t" substr(new, 2) >> "header_mapping.txt";
  next;
}
{
  print >> "example_ragtag_short.fasta";
}' example_ragtag_original.fasta 

#Extract headers (without ">")
grep '^>' example_ragtag_original.fasta | sed 's/^>//' > original_headers.txt
grep '^>' example_ragtag_short.fasta | sed 's/^>//' > short_headers.txt

# Combine both into a mapping file
paste original_headers.txt short_headers.txt > header_mapping.txt

#Update your SLURM script
Change:
GENOME="/path/to/ragtag.scaffold.fasta" To: GENOME="/path/to/example_ragtag_short.fasta"
```

10. GMAP-GSNAP
```bash
#!/bin/bash
#SBATCH --partition=<general>
#SBATCH --account=<general_user>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=500G
#SBATCH --qos=normal
#SBATCH --time=14-00:00:00
#SBATCH --job-name="JOB TITLE"
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=gmap_%j.out
#SBATCH --error=gmap_%j.err
#SBATCH --export=ALL

set -euo pipefail
set -x

module purge
source /home/user.name/etc/profile.d/conda.sh
conda activate <environment>

GENOME="/path/to/example_ragtag_short.fasta.masked"
CDS="/path/to/Attraktion_6103_618.fasta"
OUT="/path/to/gmap_madsen"
DB_NAME="original_masked"

[[ ! -f "$GENOME" || ! -f "$CDS" ]] && { echo "Error: Missing input files!"; exit 1; }
mkdir -p "$OUT"

echo "Building GMAP index..."
gmap_build -D "$OUT" -d "$DB_NAME" -k 14 "$GENOME"

echo "Running GMAP..."
gmap -D "$OUT" -d "$DB_NAME" \
  -f gff3_gene \
  -n 1 \
  -t 24 \
  -B 5 \
  --min-intronlength=20 \
  --max-intronlength=500000 \
  --intron-density=0.1 \
  "$CDS" > "$OUT/original_gmap.gff3"

if [[ -f "$OUT/original_gmap.gff3" ]]; then
    echo "✅ GMAP completed."
    echo "📊 GFF stats:"
    grep -c "^##" "$OUT/original_gmap.gff3" && grep -c -v "^##" "$OUT/original_gmap.gff3"
    echo "🧾 Gene count:"
    grep -c "gene" "$OUT/original_gma p.gff3"
else
    echo "❌ GMAP failed. Check $OUT/*.log and gmap_mj.err"
    exit 1
fi
conda deactivate
```
We reformatted the GMAP output into AUGUSTUS-compatible hints because AUGUSTUS requires a specific GFF structure with standardized fields—such as source=hints, feature types like CDSpart, and a src= tag (e.g., src=E)—to correctly interpret external evidence during gene prediction. Without this format, AUGUSTUS will either ignore the hints or throw parsing errors.

```bash
awk 'BEGIN {OFS="\t"} {
  if ($3 == "exon" || $3 == "CDS" || $3 == "intron") {
    hint_type = ($3 == "CDS") ? "CDSpart" : $3;
    print $1, "hints", hint_type, $4, $5, ".", $7, ".", "src=E;grp=gene"NR
  }
}' /path/to/original_gmap_hints.gff > /path/to/original_augustus_hints.gff

```
11. 
===================================

# Dillinger
## _The Last Markdown Editor, Ever_
[![N|Solid](https://cldup.com/dTxpPi9lDf.thumb.png)](https://nodesource.com/products/nsolid)
Dillinger is a cloud-enabled, mobile-ready, offline-storage compatible,
AngularJS-powered HTML5 Markdown editor.
- Type some Markdown on the left

## Features
Markdown is a lightweight markup language based on the formatting conventions
that people naturally use in email.
> The overriding design goal for Markdown's
> formatting syntax is to make it as readable
For production environments...
```sh
npm install --production
NODE_ENV=production node app
```
| Plugin | README |
| ------ | ------ |
| Dropbox | [plugins/dropbox/README.md][PlDb] |
**Free Software, Hell Yeah!**

