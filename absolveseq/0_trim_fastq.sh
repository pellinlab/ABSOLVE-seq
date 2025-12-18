read1=$1
read2=$2
threads=$3
TRIMMOMATIC='/Users/jieconglin/Documents/jiecong/ABSOLVE_WTCas9/Trimmomatic-0.36'
FASTQ_DIR='/Users/jieconglin/Documents/jiecong/ABSOLVE_WTCas9'
OUT_DIR='/Users/jieconglin/Documents/jiecong/ABSOLVE_WTCas9/wtCas9_trimmed'
java -jar $TRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 \
${FASTQ_DIR}/${read1}_R1.fastq.gz \
${FASTQ_DIR}/${read1}_R2.fastq.gz \
${OUT_DIR}/trimmed_${read1}_R1_paired.fastq.gz ${OUT_DIR}/trimmed_${read1}_R1_unpaired.fastq.gz \
${OUT_DIR}/trimmed_${read1}_R2_paired.fastq.gz ${OUT_DIR}/trimmed_${read1}_R2_unpaired.fastq.gz \
ILLUMINACLIP:/Users/jieconglin/Documents/jiecong/ABSOLVE_WTCas9/scripts/adapters.fa:2:30:10 LEADING:3 TRAILING:3 AVGQUAL:20 MINLEN:48 \
-threads $threads