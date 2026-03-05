# ITS-1 and ITS-2 amplicon sequences fro the reference genomes

```bash
#Generating the working environment
WORK_DIR='/data/pam/team333/jg34/scratch/'
cd ${WORK_DIR}
mkdir -p ITS_trichuris
cd ${WORK_DIR}/ITS_trichuris

# ===================================
# STEP 1: Download and Prepare Reference Genomes
# ==========================================
mkdir -p 00_ref
cd 00_ref

# Uncompress Trichuris trichiura genome
# Source: https://github.com/stephenrdoyle/ancient_trichuris/blob/master/02_data/trichuris_trichiura.fa.gz
gunzip trichuris_trichiura.fa.gz

# Download Trichuris suis genome (not used for amplicons, but available if needed)
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/trichuris_suis/PRJNA179528/trichuris_suis.PRJNA179528.WBPS19.genomic.fa.gz

# ==========================================
# STEP 2: Define Primers for ITS Regions
# ==========================================

# Define ITS2 primers
cat > primers.tab <<'EOF'
ITS2	ATGTCGACGCTACGCCTGTC	TAGCCTCGTCTGATCTGAGG
EOF

# ==========================================
# STEP 3: Extract ITS2 from Trichuris trichiura Genome
# ==========================================

mkdir -p 01_its2

# Extract amplicons from both strands using primers
seqkit -j 8 amplicon -m 2 -p primers.tab -r 1:-1 trichuris_trichiura.fa > 01_its2/trichiura.ITS2.raw.fa
seqkit -j 8 amplicon -m 2 -p primers.tab -r 1:-1 --bed trichuris_trichiura.fa > 01_its2/trichiura.ITS2.raw.bed

# Filter amplicons to expected length
seqkit seq -m 450 -M 800 01_its2/trichiura.ITS2.raw.fa > 01_its2/trichiura.ITS2.fa

# ==========================================
# STEP 4: Download ITS Reference Sequences from NCBI
# ==========================================

mkdir -p 02_ncbi
cat > 02_ncbi/trichuris_its_accessions.txt << 'EOF'
AY439019
FN543136
FR849690
FR849691
FR849693
FR849695
GQ301551
GQ301553
GQ301554
GQ301555
HQ844233
JF680987
JF690952
KJ588133
KJ588134
KJ588135
KJ588142
KJ588147
KJ588152
KJ588156
KJ588160
KJ588164
KJ588165
KJ588166
KJ588167
KL363401
KL363419
KL367758
KP336476
KP336477
KP336483
KP336484
KT186234
KT581252
KX961642
KX961643
KX961644
KX961653
LC320148
LC377164
LC500223
LC512871
LC512872
LR535751
LR743565
LR743566
LR743567
LR743568
LR743569
LR743572
LR743571
LR743570
LR743573
LR743574
MG007679
MG007680
MG007681
MG007682
MG007683
MG007711
MG007715
MG656438
MG656439
MG656440
MG656441
MG656442
MG656443
MG656444
MN428343
MN428344
MN428345
MN428346
MN428347
MN428348
MN428349
MN428350
MN428351
MN428352
MN428353
MN428398
MN428406
MN447323
MN447326
MN447327
MN447328
OP824836
KC006432
PRJEB126
PRJNA179528
PRJNA208415
PRJNA208416
ERP128004
EOF

OUT=02_ncbi/trichuris_its_genbank.fasta
: > "$OUT"

while read -r ACC; do
  [[ -z "$ACC" ]] && continue
  curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=text" \
    >> "$OUT"
done < 02_ncbi/trichuris_its_accessions.txt

# Extract amplicons from GenBank sequences and deduplicate
seqkit -j 8 amplicon -m 2 -p primers.tab -r 1:-1 "$OUT" > 02_ncbi/trichuris_its_genbank.amplicon.fasta
seqkit rmdup -s 02_ncbi/trichuris_its_genbank.amplicon.fasta -o 02_ncbi/trichuris_its_genbank.amplicon.derep.fasta

# ==========================================
# STEP 5: Identify ITS2 in Trichuris incognita Genome
# ==========================================

mkdir -p 03_map_incognita

# Align ITS2 references to incognita genome
minimap2 -x asm5 -t 8 trichuris_incognita.fa 02_ncbi/trichuris_its_genbank.amplicon.derep.fasta \
  > 03_map_incognita/incognita_vs_genebank.paf

awk 'BEGIN{OFS="\t"} {id=$11; if(id>=400) print $0, id}' \
  03_map_incognita/incognita_vs_genebank.paf \
  > 03_map_incognita/incognita_vs_genebank.filtered.paf

awk 'BEGIN{OFS="\t"} {q=$1; t=$6; s=$8; e=$9; if(s>e){tmp=s;s=e;e=tmp} name=q"__"t":"s"-"e; print t, s, e, name}' \
  03_map_incognita/incognita_vs_genebank.filtered.paf \
  > 03_map_incognita/incognita_ITS2_candidates.bed

# Extract candidate ITS2 regions
seqkit subseq --bed 03_map_incognita/incognita_ITS2_candidates.bed \
  trichuris_incognita.fa \
  > 03_map_incognita/incognita_ITS2_candidates.fa

# Deduplicate candidate sequences
seqkit rmdup -s 03_map_incognita/incognita_ITS2_candidates.fa \
  -o 03_map_incognita/incognita_ITS2_candidates.derep.fa

# Extract best match (optional alignment checks skipped for clarity)
seqkit grep -n -p "TTRE_chr2_scaffold39_13832528-13833110" \
  03_map_incognita/incognita_ITS2_candidates.derep.fa \
  > 03_map_incognita/incognita_ITS2.primary.fa

# Final primer-bounded region
seqkit subseq -r 194:791 \
  03_map_incognita/incognita_ITS2.primary.fa \
  > 03_map_incognita/incognita_ITS2.primerbounded.fa

# ==========================================
# STEP 6: Combine ITS2 References for Downstream Analyses
# ==========================================

cat 01_its2/trichiura.ITS2.fa \
    03_map_incognita/incognita_ITS2.primerbounded.fa \
    > 00_ref/01_its2/trichuris_tri_inc_ITS2_refs.fa

# Label sequences
cat > 00_ref/01_its2/trichuris_tri_inc_ITS2_refs.labeled.fa << 'EOF'
>Trichuris_trichiura_ITS2
ATGTCGACGCTACGCCTGTCTGAGGGTCGTTAAGCATAATAGCGAATGCGCCGCTCAGGC
TACAGGTTGAGGTTGGTGGCGAGCACCGGACAAACCTGCATCCGCGCGCGAGCGAGCGTG
ACGCCGCAGCTCCGTTGCCAGCGAGCCGGGATGGCAACTGGTAGCCGGAGCAGCGGAGAG
CGGCCAAGTCAGCGTAGCGCGAAGACTACCCGACTTGGCTACCGGCCGCGCCGTCGGCGT
ACAGCAGTTGAGCAGGGAGCGGTGACGTTGCCGGGCTCGTAGTAGCCAAGTGTTCGTCGT
CGTTGCACGGCAGCAGCAGCAGCAGTCGACGACGACGAAACCGTTCGCCTTGCTGCGGCG
GCGGCGGGCACCGCTCGACCGACGACGGCGGCAAACGTTGCTGACACACTGCTACCGGTC
AGGCGAGCTGCCGGTAGCGGCGTCTTTGCCGTTCGTCAACTAGCGGTTCGCCGCGTGGCT
CCACCCGTACGTCTCCTCGATGTTGACCTCAGATCAGACGAGGCTA
>Trichuris_incognita_ITS2
ATGTCGACGCTACGCCTGTCTGAGGGTCGTTACGAAATAAAGCAAATGCGCCGCTCAGCA
GGCTGCTACTGCCGCTGGGAACTAGCGGTAGCAAGCAGCGCGGGTACGGCTGCCCGTTGG
TTGGTCTCAGCGAGCGCGACGCCGCAGCTGCTCCTGCTGCTACTGGCAGCGACGGCAGGT
GGCCGTCATCGCTGACAGGCAGCCGCAGCTTCTGCGGAGAGCGGCTAACTCAGCGCAGTA
CGGAAGCTGCCCGAGTTGGCTATGTCGCTACATCGTCAGCGTAAAGCCGGCGAACGACCG
TTGACCACCGAGCGACCACTGCGGGCGCGAGCGCAGTCGTCCGTCTTCGTCGCCGCCCCC
TAGATCGACGGCAGCGGCAGTCGACTAGAAGACGACGACGGCGTCTGCTGCTACGCGCTC
TGCAGCGGCGATCTCGGCTGCGACGGCAAACGTTGCTGACACACCGCTGCTACCACTTGC
TGCTGCTGCTGCTGCTGGTAGCGGCGGCGTCTTTGCCGTACGCTGGTCAGTCACGGCCGC
AACGCCACTGCGCTCCCGTTGTTGTCTCTTCGTGTTGACCTCAGATCAGACGAGGCTA
EOF

# ==========================================
# STEP 7: Extract ITS1 Amplicons from Reference Genomes
# ==========================================

# Set up paths and primers
OUT=00_ref/ITS1_ref_build
GENOME_TRI=00_ref/trichuris_trichiura.fa
GENOME_INC=00_ref/trichuris_incognita.fa
mkdir -p "$OUT/01_genbank" "$OUT/02_trichiura" "$OUT/03_incognita"

ITS1_F="ACTGGCCGAACCAAGCCATC"
ITS1_R="TCTACGAGCCAAGTGATCCAC"

cat > "$OUT/primers_ITS1.tab" <<EOF
ITS1	${ITS1_F}	${ITS1_R}
EOF

# Download GenBank ITS1 refs
cat > "$OUT/01_genbank/ITS1_accessions.txt" << 'EOF'
GQ301555
AM992984
AM992997
KJ588071
KP336484
KT344830
MH390366
AM992999
MG656442
EOF

GBRAW="$OUT/01_genbank/ITS1_genbank_raw.fa"
: > "$GBRAW"
while read -r ACC; do
  [[ -z "$ACC" ]] && continue
  curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=text" \
    >> "$GBRAW"
done < "$OUT/01_genbank/ITS1_accessions.txt"

# Deduplicate
GBAMP_DEREP="$OUT/01_genbank/ITS1_genbank.amplicon.derep.fa"
seqkit rmdup -s "$GBRAW" -o "$GBAMP_DEREP"

# Map to T. trichiura genome
minimap2 -t 8 -x asm5 --secondary=no "$GENOME_TRI" "$GBAMP_DEREP" \
  > "$OUT/02_trichiura/genbankITS1_vs_trichiura.paf"

# Extract best match BED
awk 'BEGIN{OFS="\t"}{t=$6; s=$8; e=$9; if(s>e){tmp=s;s=e;e=tmp} ps=s-500; if(ps<0)ps=0; pe=e+500; print t, ps, pe, "trichiura_ITS1_pad500"}' \
  "$OUT/02_trichiura/genbankITS1_vs_trichiura.paf" > "$OUT/02_trichiura/trichiura_ITS1.pad500.bed"

# Extract sequence and trim with primers
seqkit subseq --bed "$OUT/02_trichiura/trichiura_ITS1.pad500.bed" "$GENOME_TRI" \
  > "$OUT/02_trichiura/trichiura_ITS1.pad500.fa"

seqkit -j 8 amplicon -m 2 -p "$OUT/primers_ITS1.tab" -r 1:-1 "$OUT/02_trichiura/trichiura_ITS1.pad500.fa" \
  > "$OUT/02_trichiura/ITS1_trichiura.fa"

# Repeat for T. incognita
minimap2 -t 8 -x asm5 --secondary=no "$GENOME_INC" "$OUT/02_trichiura/ITS1_trichiura.fa" \
  > "$OUT/03_incognita/trichiuraITS1_vs_incognita.paf"

awk 'BEGIN{OFS="\t"}{t=$6; s=$8; e=$9; if(s>e){tmp=s;s=e;e=tmp} ps=s-500; if(ps<0)ps=0; pe=e+500; print t, ps, pe, "incognita_ITS1_pad500"}' \
  "$OUT/03_incognita/trichiuraITS1_vs_incognita.paf" > "$OUT/03_incognita/incognita_ITS1.pad500.bed"

seqkit subseq --bed "$OUT/03_incognita/incognita_ITS1.pad500.bed" "$GENOME_INC" \
  > "$OUT/03_incognita/incognita_ITS1.pad500.fa"

seqkit -j 8 amplicon -m 2 -p "$OUT/primers_ITS1.tab" -r 1:-1 "$OUT/03_incognita/incognita_ITS1.pad500.fa" \
  > "$OUT/03_incognita/ITS1_incognita.fa"

# ==========================================
# STEP 8: Combine and Label ITS1 Sequences
# ==========================================

cat "$OUT/02_trichiura/ITS1_trichiura.fa" \
    "$OUT/03_incognita/ITS1_incognita.fa" \
    > "$OUT/trichuris_tri_inc_ITS1_refs.labeled.fa"

# Merge ITS1 and ITS2 for downstream analysis
cat "$OUT/trichuris_tri_inc_ITS1_refs.labeled.fa" \
    00_ref/01_its2/trichuris_tri_inc_ITS2_refs.labeled.fa \
    > 00_ref/trichuris_tri_inc_ITS1_ITS2_refs.labeled.fa

# Split by marker
mkdir -p 00_ref/by_marker

seqkit grep -r -p "ITS1$" 00_ref/trichuris_tri_inc_ITS1_ITS2_refs.labeled.fa \
  > 00_ref/by_marker/trichuris_tri_inc_ITS1_refs.fa

seqkit grep -r -p "ITS2$" 00_ref/trichuris_tri_inc_ITS1_ITS2_refs.labeled.fa \
  > 00_ref/by_marker/trichuris_tri_inc_ITS2_refs.fa
```
