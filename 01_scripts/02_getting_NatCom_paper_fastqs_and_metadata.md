# Geeting the fastq files form the 687 samples available in SRA archive (code PRJNA1283282)

We will also generate a metadata file from it

```bash
#################################################################
# Step 1: Download and prepare FASTQ files and sample metadata  #
#################################################################

# Create working directory structure
cd ${WORK_DIR}/ITS_trichuris
mkdir 01_raw
mkdir -p 01_raw/PRJNA1283282/{runinfo,fastq,sra}
cd 01_raw/PRJNA1283282

# Load required tools
module load sra-tools/3.0.7--h9f5acd7_0

# -----------------------------------------------------------
# Download RunInfo metadata file manually from NCBI:
# 1. Go to: https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1283282
# 2. Click "Send to → Run Selector"
# 3. Click "Metadata → Download RunInfo"
# 4. Save as: PRJNA1283282.runinfo.csv in the runinfo/ folder
# -----------------------------------------------------------

# Create accession list
cut -d',' -f1 runinfo/PRJNA1283282.runinfo.csv | tail -n +2 > runinfo/SRR_AccList.txt

# Download SRA files
prefetch --option-file runinfo/SRR_AccList.txt --output-directory sra

# Convert SRA to FASTQ
fasterq-dump sra/*/*.sra \
  --threads 8 \
  --split-files \
  --outdir fastq \
  --temp tmp

# Extract relevant columns from metadata (Run, BioSample, geo_loc_name, Sample Name)
gawk '
BEGIN {
  FPAT = "([^,]+)|(\"[^\"]+\")";
  OFS = "\t"
}
NR==1 {
  for (i=1; i<=NF; i++) {
    f=$i; gsub(/^\"|\"$/, "", f)
    if (f=="Run") run=i
    if (f=="BioSample") bios=i
    if (f=="geo_loc_name") geo=i
    if (f=="Sample Name") sname=i
  }
  print "Run","BioSample","country","Sample_Name"
  next
}
{
  r=$run; b=$bios; c=$geo; s=$sname
  gsub(/^\"|\"$/, "", r); gsub(/^\"|\"$/, "", b);
  gsub(/^\"|\"$/, "", c); gsub(/^\"|\"$/, "", s)
  print r,b,c,s
}' runinfo/PRJNA1283282.runinfo.csv > runinfo/run_metadata_core.tsv

# Get list of FASTQ runs
ls fastq/*.fastq* \
  | sed 's#.*/##' \
  | sed -E 's/\.fastq(\.gz)?$//' \
  | sort -u > runinfo/run_list_from_fastq.txt

# Map runs to FASTQ paths
awk '
{
  r=$1
  fqgz="fastq/"r".fastq.gz"
  fq="fastq/"r".fastq"
  if (system("[ -f "fqgz" ]")==0)
    print r "\t" fqgz "\t"
  else if (system("[ -f "fq" ]")==0)
    print r "\t" fq "\t"
  else
    print r "\t\t"
}' runinfo/run_list_from_fastq.txt > runinfo/run_fastq_paths.tsv

# Sort and join metadata tables
sort runinfo/run_metadata_core.tsv > runinfo/run_metadata_core.sorted.tsv
sort runinfo/run_fastq_paths.tsv > runinfo/run_fastq_paths.sorted.tsv

join -t $'\t' \
  runinfo/run_metadata_core.sorted.tsv \
  runinfo/run_fastq_paths.sorted.tsv > runinfo/metadata_full.tsv

#################################################################
# Step 2: Fix and normalize metadata table (path, spacing, country)
#################################################################

BASEDIR="${WORK_DIR}/ITS_trichuris/01_raw"
RUNINFO="${BASEDIR}/runinfo/PRJNA1283282.runinfo.csv"
FASTQDIR="${BASEDIR}/fastq"
OUTDIR="${BASEDIR}/runinfo"

mkdir -p "${OUTDIR}"

# Absolute FASTQ paths
ls "${FASTQDIR}"/*.fastq* 2>/dev/null \
  | awk -v fqdir="${FASTQDIR}" '
    BEGIN{OFS="\t"}
    {
      fn=$0; sub(/^.*\//,"",fn); run=fn
      sub(/\.fastq(\.gz)?$/,"",run)
      print run, $0
    }' | sort -k1,1 > "${OUTDIR}/run_fastq_abs.tsv"

# Parse and normalize metadata columns
gawk '
BEGIN {
  FPAT = "([^,]+)|(\"[^\"]+\")"
  OFS  = "\t"
}
NR==1 {
  for (i=1; i<=NF; i++) {
    f=$i; gsub(/^\"|\"$/, "", f)
    if (f=="Run") run=i
    if (f=="BioSample") bios=i
    if (f=="geo_loc_name") geo=i
    if (f=="Sample Name") sname=i
  }
  print "Run","BioSample","geo_loc_name","Sample_Name"
  next
}
{
  r=$run; b=$bios; g=$geo; s=$sname
  gsub(/^\"|\"$/, "", r)
  gsub(/^\"|\"$/, "", b)
  gsub(/^\"|\"$/, "", g)
  gsub(/^\"|\"$/, "", s)
  print r,b,g,s
}' "${RUNINFO}" | sort -k1,1 > "${OUTDIR}/run_metadata_core.tsv"

# Merge with FASTQ paths
join -t $'\t' -1 1 -2 1 \
  "${OUTDIR}/run_metadata_core.tsv" \
  "${OUTDIR}/run_fastq_abs.tsv" > "${OUTDIR}/metadata_fixed.tsv"

# Normalize country names
awk -F'\t' 'BEGIN{OFS="\t"}
NR==1 {print $0,"country_norm"; next}
{
  country=$3; norm=country
  gsub(/'\''/,"",norm)
  gsub(/[[:space:]]+/,"_",norm)
  gsub(/[^A-Za-z0-9_]/,"",norm)
  print $0,norm
}' "${OUTDIR}/metadata_fixed.tsv" > "${OUTDIR}/metadata_fixed.wide.tsv"

# Write CSV version (escaped for Excel)
awk -F'\t' '
BEGIN{OFS=","}
function csvq(x,   y){
  y=x; gsub(/"/,"\"\"",y)
  return (y ~ /[,
]/) ? "\""y"\"" : y
}
NR==1 {print "Run","BioSample","geo_loc_name","Sample_Name","FASTQ","country_norm"; next}
{
  print csvq($1),csvq($2),csvq($3),csvq($4),csvq($5),csvq($6)
}' "${OUTDIR}/metadata_fixed.wide.tsv" > "${OUTDIR}/metadata_fixed.csv"

# Export simplified tab/CSV versions (e.g., for plotting, clustering)
cut -f1,2,4,5,6,7 "${OUTDIR}/metadata_fixed.wide.tsv" > "${OUTDIR}/metadata_fixed.wide.nospace.tsv"
cut -d',' -f1,2,4,5,6,7 "${OUTDIR}/metadata_fixed.csv" > "${OUTDIR}/metadata_fixed.nospace.csv"
```
In this code, we generate metadata_fixed.wide.nospace.tsv, which will be used in late scripts
