#Processing raw reads running fastp v0.21.0 and FastQC v0.11.4 and multiqc v1.13.dev0
cat samples.txt |  parallel -n 1 --jobs 1 --joblog fastp_run.txt -I % fastp -i Raw_reads/%_R1_001.fastq.gz -I Raw_reads/%_R2_001.fastq.gz -o fastp_processed/%_fastp_R1_001.fastq.gz -O fastp_processed/%_fastp_R2_001.fastq.gz
fastqc -f fastq -out ../FASTQC_fastp_processed *fastq.gz -t 30
#Redoing but cutting off first 14nts based on fastqc results
cat samples.txt |  parallel -n 1 --jobs 1 --joblog fastp_run.txt -I % fastp -i Raw_reads/%_R1_001.fastq.gz -I Raw_reads/%_R2_001.fastq.gz -o fastp_processed/%_fastp_R1_001.fastq.gz -O fastp_processed/%_fastp_R2_001.fastq.gz -f 14 -F 14
fastqc -f fastq -out ../FASTQC_fastp_processed_trimmed *fastq.gz -t 30

#Merging samples as everything was run on 2 lanes
ls -1 | awk -F "_" '{print $1"_"$2"}' > ../merged_samples.txt
cat ../merged_samples.txt |  parallel -n 1 --jobs 1 --joblog ../merging.txt -I % "cat %_L001_fastp_R1_001.fastq.gz %_L002_fastp_R1_001.fastq.gz > ../merged_fastp_processed_trimmed/%_L001_L002_fastp_R1_001.fastq.gz"
cat ../merged_samples.txt |  parallel -n 1 --jobs 1 --joblog ../merging.txt -I % "cat %_L001_fastp_R2_001.fastq.gz %_L002_fastp_R2_001.fastq.gz > ../merged_fastp_processed_trimmed/%_L001_L002_fastp_R2_001.fastq.gz"

#Making salmon index
grep "^>" GCF_000325575.1_P_alecto_ASM32557v1_genomic.fna | cut -d " " -f 1 | tr -d ">" > decoys.txt
cat GCF_000325575.1_P_alecto_ASM32557v1_rna.fna GCF_000325575.1_P_alecto_ASM32557v1_genomic.fna > P_alecto_gentrome.fa
salmon index -t P_alecto_gentrome.fa -d decoys.txt -p 30 -i P_alecto__salmon_index --gencode

#Running salmon v1.4.0
cat P_alecto_samples.txt | parallel -n 1 --jobs 1 -I % --joblog P_alecto_salmon_run.txt salmon quant -i Salmon_indexes/P_alecto/P_alecto_salmon_index -l A -1 merged_fastp_processed_trimmed/%_L001_L002_fastp_R1_001.fastq.gz -2 merged_fastp_processed_trimmed/%_L001_L002_fastp_R2_001.fastq.gz --validateMappings -o Salmon_quants/P_alecto/%_transcripts_quant
