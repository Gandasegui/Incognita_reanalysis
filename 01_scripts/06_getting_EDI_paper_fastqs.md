# Geeting the necessary fastq files form the Emerging Infectious Diseases (EDI) paper

The paper by Venkatesan et al. indlude fastqfile form the ITSs as well as from mitochondrial seqence

They analyised samples from more than 80 participants

### We need to:
- Keep only the fastq files that correspond to ITS-1 and ITS-2 sequences
- Use the sample ID form the Supplementary Data 6 in Nat Com paper to keep only those fastqs from the same participant 

```bash
# ===============================
# 0. Setup workspace
# ===============================
cd ${WORK_DIR}/ITS_trichuris
mkdir 03_validation
mkdir 03_validation/01_illumina_PRJNA1131306
cd 03_validation/01_illumina_PRJNA1131306

# ===============================
# Step 1: Prepare RunInfo metadata
# ===============================
module load sra-tools/3.0.7--h9f5acd7_0

# parse_runinfo.py extracts Run, Marker, PatientID and normalizes country field
cat > parse_runinfo.py <<'PY'
#!/usr/bin/env python3
import csv
import sys
import re

infile = sys.argv[1]
outfile = sys.argv[2]

def norm_country_from_geo(geo_loc: str) -> str:
    """
    geo_loc_name examples:
      - "Cote d'Ivoire"
      - "Tanzania"
      - "Laos"
      - sometimes "Country: something" or includes extra text
    We normalize to: Tanzania, Cote_dIvoire, Laos (as requested).
    """
    if geo_loc is None:
        return "NA"
    s = geo_loc.strip()

    # common normalizations
    s_low = s.lower()
    if "tanzania" in s_low:
        return "Tanzania"
    if "laos" in s_low:
        return "Laos"

    # handle Cote d'Ivoire variations
    # e.g., "Cote d'Ivoire", "Côte d’Ivoire", etc.
    if "cote" in s_low and "ivoire" in s_low:
        return "Cote_dIvoire"
    if "côte" in s_low and "ivoire" in s_low:
        return "Cote_dIvoire"

    # fallback: sanitize spaces/punctuation into underscore-friendly token
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^\w]+", "", s)  # remove non-word chars
    return s if s else "NA"

with open(infile, newline='') as f, open(outfile, "w", newline='') as out:
    r = csv.DictReader(f)

    out.write("\t".join(["Run","PatientID","Marker","LibraryName","SampleName","Country"]) + "\n")

    for row in r:
        run = (row.get("Run") or "").strip()
        lib = (row.get("Library Name") or "").strip()
        sam = (row.get("Sample Name") or "").strip()

        # IMPORTANT: use geo_loc_name (not geo_loc_name_country)
        geo = (row.get("geo_loc_name") or "").strip()
        country = norm_country_from_geo(geo)

        # Parse Library Name like BL-B8823-ITS-2 or BL-7374-ITS-1, etc.
        patient = ""
        marker = ""
        parts = lib.split("-")
        if len(parts) >= 3:
            patient = "-".join(parts[0:2])   # BL-B8823
            marker  = "-".join(parts[2:])    # ITS-2 / ITS-1 / cox-1 ...
        else:
            m = re.search(r"(BL-[A-Za-z0-9]+)", lib)
            patient = m.group(1) if m else ""
            marker = ""

        out.write("\t".join([run, patient, marker, lib, sam, country]) + "\n")
PY
chmod +x parse_runinfo.py

IN=PRJNA1131306.RunInfo.csv
OUT=PRJNA1131306.parsed.tsv
./parse_runinfo.py "$IN" "$OUT"

# Filter ITS-1 and ITS-2 entries
grep -P '^Run\t|\tITS-1\t|\tITS-2\t' "$OUT" > PRJNA1131306.ITS1_ITS2.tsv
grep -P '^Run\t|\tITS-1\t' "$OUT" > PRJNA1131306.ITS1.tsv
grep -P '^Run\t|\tITS-2\t' "$OUT" > PRJNA1131306.ITS2.tsv

cut -f1 PRJNA1131306.ITS1.tsv | tail -n +2 > ITS1_SRRs.txt
cut -f1 PRJNA1131306.ITS2.tsv | tail -n +2 > ITS2_SRRs.txt

# ===============================
# 2. Find overlapping samples with Nature paper
# ===============================
cd ..
# Have to extract Nautre_paper_metadata.csv from Supplementary Data 6: Subset metadata for study participants aged 6–18 years
cut -d',' -f2 Nautre_paper_metadata.csv | tail -n +2 | sort -u > nature.patientIDs.txt
cut -f2 PRJNA1131306.ITS1_ITS2.tsv | tail -n +2 | sort -u > ITS1_ITS2.patientIDs.txt
comm -12 nature.patientIDs.txt ITS1_ITS2.patientIDs.txt > shared.patientIDs.txt

# Filter PRJNA table to shared samples
awk -F'\t' 'NR==FNR{keep[$1]=1; next} NR==1 || ($2 in keep && ($3=="ITS-1" || $3=="ITS-2"))' \
  shared.patientIDs.txt PRJNA1131306.ITS1_ITS2.tsv > PRJNA1131306.ITS1_ITS2_shared.tsv

cut -f1 PRJNA1131306.ITS1_ITS2_shared.tsv | tail -n +2 | sort -u > ITS1_ITS2_shared_SRRs.txt

# ===============================
# 3. Download & convert SRA to FASTQ
# ===============================
mkdir -p 02_sra_shared40_EDI/{sra,fastq,logs}
prefetch --progress -O 02_sra_shared40_EDI/sra --option-file ITS1_ITS2_shared_SRRs.txt

find 02_sra_shared40_EDI/sra -name "*.sra" > 02_sra_shared40_EDI/sra/sra_files.list

while read -r SRA; do
  echo "Converting $SRA"
  fasterq-dump -e 8 -O 07_sra_shared40/fastq --split-files "$SRA"
done < 02_sra_shared40_EDI/sra/sra_files.list

# ===============================
# 4. Match Sample IDs to FASTQ paths
# ===============================
META_TSV="01_raw/runinfo/metadata_fixed.wide.nospace.tsv"
NATURE_CSV="Nautre_paper_metadata.csv"

awk -F',' 'NR>1{print $2"\t"$3}' "$NATURE_CSV" | sed 's/\r\$//' > nature.PatientID_to_SampleID.tsv

awk 'NR==FNR{keep[$1]=1; next} ($1 in keep)' \
  shared.patientIDs.txt nature.PatientID_to_SampleID.tsv > shared40.PatientID_to_SampleID.tsv

awk -F'\t' '$1!="BL-B8823"' shared40.PatientID_to_SampleID.tsv > shared39.PatientID_to_SampleID.tsv

awk -F'\t' 'NF>=4{print $3"\t"$4"\t"$1}' "$META_TSV" | sed 's/\r\$//' > nature.SampleID_to_fastq_and_SRR.tsv

awk 'NR==FNR{fq[$1]=$2"\t"$3; next} ($2 in fq){print $1"\t"$2"\t"fq[$2]}' \
  nature.SampleID_to_fastq_and_SRR.tsv shared39.PatientID_to_SampleID.tsv > shared40.Nature_fastqs.tsv

# ===============================
# 5. Organize overlapping FASTQ files
# ===============================

OUTDIR="08_shared40_Nature_fastq"
mkdir -p "$OUTDIR"

cut -f3 shared40.Nature_fastqs.tsv | while read -r FQ; do
  if [ -s "$FQ" ]; then
    cp -v "$FQ" "$OUTDIR"/
  else
    echo "MISSING_OR_EMPTY: $FQ" >&2
  fi
done

#let's remove the dupiucate:
cat 01_illumina_PRJNA1131306/PRJNA1131306.ITS1_ITS2.tsv | grep '8823'
#SRR31228438     BL-B8823        ITS-2   BL-B8823-ITS-2  BL-B8823        Cote_dIvoire
#SRR31228722     BL-B8823        ITS-1   BL-B8823-ITS-1  BL-B8823        Cote_dIvoire
rm 02_sra_shared40_EDI/fastq/SRR31228438*
rm 02_sra_shared40_EDI/fastq/SRR31228722*
```

