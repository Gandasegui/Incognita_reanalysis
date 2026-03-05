# Haplotyoe calling, comparison with ASV and seq depth on ONT data from NatCom paper

```bash
####################################
# Step 0: Setting directories and loading required modules
####################################
cd ${WORK_DIR}/ITS_trichuris
mkdir 02_hap
cd 02_hap
#modules
module load minimap2/2.28--h577a1d6_4
module load samtools/1.9--h91753b0_8
module load racon/1.3.3--he860b03_1
module load medaka/2.1.1
module load cutadapt/1.15--py36_0
module load seqkit/v2.9.0

####################################
# Step 1: Configuration
####################################
THREADS=8
BASEDIR="/data/pam/team333/jg34/scratch/ITS_trichuris/01_raw"
META="${BASEDIR}/runinfo/metadata_fixed.wide.nospace.tsv"
FWD_PRIMER="ATGTCGACGCTACGCCTGTC"
REV_PRIMER="CCTCAGATCAGACGAGGCTA"
MINLEN_TRIM=480
MAXLEN_TRIM=740
ID=0.985
MIN_CLUSTER_READS=100
MEDAKA_MODEL="r1041_e82_400bps_sup_v5.2.0"

####################################
# Step 2: Primer Trimming and Length Filtering
####################################
mkdir -p 02_qc
awk -F'\t' '{print $1"\t"$4}' "$META" | while IFS=$'\t' read -r RUN FQ; do
  cutadapt -j $THREADS \
    -g "^${FWD_PRIMER}...${REV_PRIMER}$" \
    --discard-untrimmed \
    -o "02_qc/${RUN}.cutadapt.fastq" "$FQ"

  seqkit seq -g -m $MINLEN_TRIM -M $MAXLEN_TRIM "02_qc/${RUN}.cutadapt.fastq" > "02_qc/${RUN}.qc.fastq"
  #rm "02_qc/${RUN}.cutadapt.fastq"
done

seqkit stats 02_qc/*.qc.fastq | tail -n +2 > 02_qc/qc_report.txt

####################################
# Step 3: Pool Filtered Reads and Dereplicate
####################################
mkdir -p 03_pool 04_cluster
cat 02_qc/*.qc.fastq > 03_pool/all.qc.fastq
seqkit fq2fa 03_pool/all.qc.fastq > 03_pool/all.qc.fasta

vsearch --derep_fulllength 03_pool/all.qc.fasta \
  --output 04_cluster/all.derep.vsearch.fasta \
  --sizeout \
  --relabel uniq_

vsearch --fastx_filter 04_cluster/all.derep.vsearch.fasta \
  --minsize $MIN_CLUSTER_READS \
  --fastaout 04_cluster/all.derep.vsearch.min${MIN_CLUSTER_READS}.fasta

vsearch --cluster_fast 04_cluster/all.derep.vsearch.min${MIN_CLUSTER_READS}.fasta \
  --id $ID \
  --centroids 04_cluster/haplotypes.reps.vsearch.min${MIN_CLUSTER_READS}.fasta \
  --uc 04_cluster/clusters.vsearch.min${MIN_CLUSTER_READS}.uc \
  --threads $THREADS

####################################
# Step 4: Polishing Haplotypes (Racon + Medaka), then clean to remove chimeras
####################################
mkdir -p 05_polish
minimap2 -t $THREADS -x map-ont --secondary=no \
  04_cluster/haplotypes.reps.vsearch.min${MIN_CLUSTER_READS}.fasta \
  03_pool/all.qc.fastq > 05_polish/all_to_reps.paf

racon -t $THREADS \
  03_pool/all.qc.fastq \
  05_polish/all_to_reps.paf \
  04_cluster/haplotypes.reps.vsearch.min${MIN_CLUSTER_READS}.fasta \
  > 05_polish/haplotypes.racon.fasta

medaka_consensus -t $THREADS \
  -i 03_pool/all.qc.fastq \
  -d 05_polish/haplotypes.racon.fasta \
  -o 05_polish/medaka_vsearch_out \
  -m $MEDAKA_MODEL

cp 05_polish/medaka_vsearch_out/consensus.fasta 05_polish/haplotypes.polished.fasta

#Let's chek the presence of the primers
seqkit locate -p ATGTCGACGCTACGCCTGTC 05_polish/haplotypes.polished.fasta
'seqID   patternName     pattern strand  start   end     matched
uniq_6251;size=155      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       18      37      ATGTCGACGCTACGCCTGTC
uniq_9686;size=101      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       1       20      ATGTCGACGCTACGCCTGTC
uniq_8544;size=114      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       1       20      ATGTCGACGCTACGCCTGTC
uniq_8544;size=114      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       57      76      ATGTCGACGCTACGCCTGTC
uniq_8688;size=113      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       14      33      ATGTCGACGCTACGCCTGTC
uniq_5324;size=183      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       13      32      ATGTCGACGCTACGCCTGTC
uniq_8672;size=113      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       1       20      ATGTCGACGCTACGCCTGTC
uniq_4989;size=195      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       14      33      ATGTCGACGCTACGCCTGTC
uniq_4380;size=224      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       13      32      ATGTCGACGCTACGCCTGTC
uniq_5253;size=185      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       13      32      ATGTCGACGCTACGCCTGTC
uniq_4101;size=240      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       1       20      ATGTCGACGCTACGCCTGTC
uniq_4101;size=240      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       21      40      ATGTCGACGCTACGCCTGTC
uniq_4101;size=240      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       41      60      ATGTCGACGCTACGCCTGTC
uniq_5227;size=186      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       13      32      ATGTCGACGCTACGCCTGTC
uniq_5227;size=186      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       33      52      ATGTCGACGCTACGCCTGTC
uniq_7514;size=129      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       14      33      ATGTCGACGCTACGCCTGTC
uniq_7665;size=127      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       1       20      ATGTCGACGCTACGCCTGTC
uniq_5336;size=182      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       13      32      ATGTCGACGCTACGCCTGTC
uniq_368;size=1475      ATGTCGACGCTACGCCTGTC    ATGTCGACGCTACGCCTGTC    +       13      32      ATGTCGACGCTACGCCTGTC'

seqkit locate -p CCTCAGATCAGACGAGGCTA 05_polish/haplotypes.polished.fasta
#Nothing here, so we have to work on the FWD sequence

#After doing a clustalW in EMBL, I can see that all samples begins with the same pattern 'AGGGTCGTTA'
#Before that all are chimeras and repeas fo the fwd primer, lo let's remove it
seqkit locate -p AGGGTCGTTA 05_polish/haplotypes.polished.fasta
#IT look like all have in record, except for
seqkit locate -p AGGGTCGTTA 05_polish/haplotypes.polished.fasta | grep '6314'
'uniq_6314;size=154      AGGGTCGTTA      AGGGTCGTTA      +       1       10      AGGGTCGTTA
uniq_6314;size=154      AGGGTCGTTA      AGGGTCGTTA      +       26      35      AGGGTCGTTA'

MOTIF="AGGGTCGTTA"

#first run
cutadapt -e 0.2 \
  -g "${MOTIF}" \
  --discard-untrimmed \
  -o 05_polish/haplotypes.polished.trim.fasta 05_polish/haplotypes.polished.fastaf

cutadapt -e 0.2 \
  -g "${MOTIF}" \
  --prefix "${MOTIF}" \
  -o 05_polish/haplotypes.polished.trim2.fasta 05_polish/haplotypes.polished.trim.fasta

#Now, the clustalW shows that the haplotype uniq_4537 has a wired begining, let's remove it
cutadapt \
  -g "AGCAGGGTAATA" \
  -o 05_polish/haplotypes.polished.clean.fasta 05_polish/haplotypes.polished.trim2.fasta

seqkit replace -p '^AGGGTCGTTA' -r '' 05_polish/haplotypes.polished.clean.fasta \
| seqkit replace -p '^([^;]+);.*' -r '$1' \
> 05_polish/haplotypes.polished.clean.renamed.fasta

seqkit rmdup -n 05_polish/haplotypes.polished.clean.renamed.fasta \
> 05_polish/haplotypes.polished.clean.renamed.dedup.fasta # no duplicates

####################################
# Step 5: Count Haplotypes per Sample
####################################
mkdir -p 06_counts

for FQ in 02_qc/*.qc.fastq; do
  RUN=$(basename "$FQ" .qc.fastq)
  minimap2 -t $THREADS -x map-ont --secondary=no \
    05_polish/haplotypes.polished.clean.renamed.fasta "$FQ" \
    > 06_counts/${RUN}.paf
done

bsub.py 9 hap_map -q long "bash map.sh"

{
  echo -e "Run\tHaplotype\tReads"
  for P in 06_counts/*.paf; do
    RUN=$(basename "$P" .paf)
    awk 'BEGIN{FS="\t"; OFS="\t"}
         {key=$6"\t"$1; if(!(key in seen)) { seen[key]=1; c[$6]++ }}
         END{for(h in c) print "'"$RUN"'", h, c[h]}' "$P"
  done
} > 06_counts/run_haplotype_counts.tsv

####################################
# Step 6: Species Assignment (BLAST against ITS2 DB)
####################################
module load blast/2.14.1--pl5321h6f7f691_0
makeblastdb -in 00_ref/01_its2/trichuris_tri_inc_ITS2_refs.labeled.fa -dbtype nucl -out 06_counts/db_ITS2refs

blastn -query 05_polish/haplotypes.polished.clean.renamed.fasta \
  -db 06_counts/db_ITS2refs \
  -task megablast \
  -max_target_seqs 5 \
  -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' \
  -out 05_polish/haplotypes_vs_ITS2refs.b6

awk 'BEGIN{FS=OFS="\t"} {
  cov = ($5 > 0) ? $4/$5 : 0;
  if(cov >= 0.70 && (!($1 in best) || $13 > best[$1])) {
    best[$1]=$13; ref[$1]=$2; pid[$1]=$3; clen[$1]=$4; qlen[$1]=$5
  }} END {
  print "Haplotype","Species","pident","aln_len","qlen","bitscore"
  for(h in best) print h, ref[h], pid[h], clen[h], qlen[h], best[h]
}' 05_polish/haplotypes_vs_ITS2refs.b6 > 06_counts/haplotype_species.tsv

# Merge counts and species
awk 'NR>1{map[$1]=$2} END{for(k in map) print k,map[k]}' 06_counts/haplotype_species.tsv | sort > 06_counts/hap2sp.tsv
tail -n +2 06_counts/run_haplotype_counts.tsv | sort -k2,2 > 06_counts/counts.sorted.tsv

{
  echo -e "Run\tHaplotype\tReads\tSpecies"
  join -t $'\t' -1 2 -2 1 06_counts/counts.sorted.tsv 06_counts/hap2sp.tsv \
    | awk 'BEGIN{OFS="\t"} {print $2,$1,$3,$4}'
} > 06_counts/run_haplotype_counts.with_species.tsv

{
  head -n 1 06_counts/run_haplotype_counts.with_species.tsv
  tail -n +2 06_counts/run_haplotype_counts.with_species.tsv \
    | sort -t $'\t' -k1,1 -k3,3nr
} > 06_counts/run_haplotype_counts.with_species.sorted.tsv

####################################
# Step 7: Assigment ASV - haplotyope at 100% identity
####################################
makdir -p 07_ASVs
mamba activate vclust
vsearch --usearch_global 05_polish/haplotypes.polished.clean.renamed.fasta \
  --db 07_ASVs/ASVs_clean.fasta \
  --id 1 \
  --strand plus \
  --maxaccepts 50 \
  --maxrejects 0 \
  --blast6out 07_ASVs/hap_vs_asv.b6

{
  echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"
  cat 07_ASVs/hap_vs_asv.b6
} > 07_ASVs/hap_vs_asv.with_header.b6

seqkit rmdup -s 07_ASVs/ASVs_clean.fasta -o ASVs_clean.derep.fasta
#[INFO] 4 duplicated records removed, but let's keep it as originally

####################################
# Step 8: Sequencing depth ???
####################################
mkdir -p 08_depth
# ONT: reads + bases per qc FASTQ
{
  echo -e "Run\treads\tbases\tfile"
  for fq in 02_qc/*.qc.fastq; do
    run=$(basename "$fq" .qc.fastq)
    awk -v RUN="$run" -v FQ="$fq" '
      NR%4==2 { bases += length($0) }
      END { reads = NR/4; printf "%s\t%d\t%d\t%s\n", RUN, reads, bases, FQ }
    ' "$fq"
  done
} > 08_depth/seqdepth_ONT_qc.tsv
```
### With this code, we are obtianing two key files:
- run_haplotype_counts.with_species.sorted.tsv; the file with the haplotype counts and species assigment
- hap_vs_asv.with_header.b6; the file that lins haplotypes and ASVs with 100% identities
