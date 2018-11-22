#!/bin/bash

cd ~/Jenna/rawdata

#### Combine all left and right reads together
cat *_1.fq.gz > reads_1.fq.gz
cat *_2.fq.gz > reads_2.fq.gz


#### Quality trim sequences
java -jar /home/heath-lab/Software/trinity/trinity-plugins/Trimmomatic-0.32/trimmomatic.jar PE ./raw_data/reads_1.fq.gz ./raw_data/reads_2.fq.gz ./raw_data/reads_P1.fq.gz ./raw_data/reads_U1.fq.gz ./raw_data/reads_P2.fq.gz ./raw_data/reads_U2.fq.gz ILLUMINACLIP:/home/heath-lab/Software/trinityrnaseq-2.0.6/trinity-plugins/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25

##################################
#### Code for denovo assembly ####
##################################
Trinity --seqType fq --max_memory 30G --left reads_P1.fq.gz --right reads_P2.fq.gz --CPU 2 --normalize_reads > run_1.log 2>&1 &


#########################################
#### Code for genome-guided assembly ####
#########################################
#### Map to Lamprey Genome Reference
tophat2 -p 8 -o ./th_out ./Lamprey_Genome/lamprey_genome ./raw_data/reads_P1.fq.gz ./raw_data/reads_P2.fq.gz


#### Code for genome guided Trinity
Trinity --genome_guided_bam ./th_out/accepted_hits.bam --genome_guided_max_intron 20000 --CPU 2 --max_memory 30G


cd ~/Jenna/raw_data

#### Build bowtie and RSEM references
~/Software/trinity/util/align_and_estimate_abundance.pl --transcripts ../Lamprey_Trinity_GG/Trinity-GG.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference


#### Bowtie align reads to transciptome and counts with RSEM
Files=`ls *_1.fq.gz | sed 's/_1.fq.gz//g'`
for F in $Files; do {
	~/Software/trinity/util/align_and_estimate_abundance.pl --transcripts ../Lamprey_Trinity_GG/Trinity-GG.fasta --seqType fq --left ${F}_1.fq.gz --right ${F}_2.fq.gz --est_method RSEM --aln_method bowtie --out_dir ./${F} --trinity_mode
};done


#### Make isoforms level counts matrix from RSEM files
~/Software/trinity/util/abundance_estimates_to_matrix.pl --est_method RSEM --name_sample_by_basedir --out_prefix Trinity_trans \
	./CC_01/RSEM.isoforms.results \
	./CC_02/RSEM.isoforms.results \
	./CC_03/RSEM.isoforms.results \
	./CC_04/RSEM.isoforms.results \
	./C5_01/RSEM.isoforms.results \
	./C5_02/RSEM.isoforms.results \
	./C5_03/RSEM.isoforms.results \
	./C5_04/RSEM.isoforms.results \
	./C10_01/RSEM.isoforms.results \
	./C10_02/RSEM.isoforms.results \
	./C10_03/RSEM.isoforms.results \
	./C10_04/RSEM.isoforms.results \
	./C30_01/RSEM.isoforms.results \
	./C30_02/RSEM.isoforms.results \
	./C30_03/RSEM.isoforms.results \
	./C30_04/RSEM.isoforms.results

#### Make genes level counts matrix from RSEM files
~/Software/trinity/util/abundance_estimates_to_matrix.pl --est_method RSEM --name_sample_by_basedir --out_prefix Trinity_genes \
	./CC_01/RSEM.genes.results \
	./CC_02/RSEM.genes.results \
	./CC_03/RSEM.genes.results \
	./CC_04/RSEM.genes.results \
	./C5_01/RSEM.genes.results \
	./C5_02/RSEM.genes.results \
	./C5_03/RSEM.genes.results \
	./C5_04/RSEM.genes.results \
	./C10_01/RSEM.genes.results \
	./C10_02/RSEM.genes.results \
	./C10_03/RSEM.genes.results \
	./C10_04/RSEM.genes.results \
	./C30_01/RSEM.genes.results \
	./C30_02/RSEM.genes.results \
	./C30_03/RSEM.genes.results \
	./C30_04/RSEM.genes.results

