#main pipeline to align FASTqs
for sample in $(ls -d */ | cut -f1 -d'/')
do
    if [ ! -e ${sample}/${sample}.dedup.numAligned.txt ]
    then
        cd ${sample}
        
        genefuse -r hg19.fa -f cancer.hg19.csv -1 ${sample}_1.fastq -2 ${sample}_2.fastq -h ${sample}_report.html -t 30 > ${sample}_genefususion_result.txt
        genefuse -r hg19.fa -f cancer.hg19.csv -u 1 -1 ${sample}_1.fastq -2 ${sample}_2.fastq -h ${sample}_uset1_report.html -t 30 > ${sample}_uset1_genefususion_result.txt

        trim_galore --fastqc --paired -j 25 ${sample}_1.fastq ${sample}_2.fastq

        bwa mem -R '@RG\tID:1\tSM:'${sample}'\tPL:illumina\tLB:lib1\tPU:unit1' -t 30 hg19.fa ${sample}_1_val_1.fq ${sample}_2_val_2.fq | samtools view -bS -F 0x900 - > ${sample}.aligned.bam
        samtools view -c -F 0x4 ${sample}.aligned.bam > ${sample}.numAligned.txt

        rm -R ${sample}_1_val_1.fq ${sample}_2_val_2.fq

        samtools sort -@ 25 -m 4G -o ${sample}.aligned.bam ${sample}.aligned.bam
        samtools index ${sample}.aligned.bam

        picard MarkDuplicates --INPUT ${sample}.aligned.bam --OUTPUT ${sample}.dedup.aligned.bam --METRICS_FILE ${sample}.metrics.txt --ASSUME_SORTED True --REMOVE_DUPLICATES True
        samtools sort -@ 25 -m 4G -o ${sample}.dedup.aligned.bam ${sample}.dedup.aligned.bam
        samtools index ${sample}.dedup.aligned.bam

        samtools view -c -F 0x4 ${sample}.dedup.aligned.bam > ${sample}.dedup.numAligned.txt
        cd ..
    fi
done


for sample in $(ls -d */ | cut -f1 -d'/')
do
	cd ${sample}
	bam-readcount/build/bin/bam-readcount -f hg19.fa  ${sample}.dedup.aligned.bam  chr17:7577547-7577547 -w 0 > ${sample}_TP53_VAF.txt
	bam-readcount/build/bin/bam-readcount -f hg19.fa  ${sample}.dedup.aligned.bam  chr17:41226461-41226461 -w 0 > ${sample}_BRCA1_VAF.txt
	bam-readcount/build/bin/bam-readcount -f hg19.fa  ${sample}.dedup.aligned.bam  chr7:140481402-140481402 -w 0 > ${sample}_BRAF_VAF.txt
	bam-readcount/build/bin/bam-readcount -f hg19.fa  ${sample}.dedup.aligned.bam  chr10:123279605-123279605 -w 0 > ${sample}_FGFR2_VAF.txt
	bam-readcount/build/bin/bam-readcount -f hg19.fa  ${sample}.dedup.aligned.bam  chr1:115256529-115256529 -w 0 > ${sample}_NRAS_VAF.txt
	bam-readcount/build/bin/bam-readcount -f hg19.fa  ${sample}.dedup.aligned.bam  chr4:1808911-1808911 -w 0 > ${sample}_fusion_denominator_VAF.txt
	cd ..
done

#TACC3 coordinates chr4:1741279-1742970 


#running gatk now
#run only once
gatk CreateSequenceDictionary -R hg19.fa

for sample in $(ls -d */ | cut -f1 -d'/')
do
	if [ ! -e ${sample}.GATK.somatic_filtered.dedup.vcf ]
    	then
		cd ${sample}
		normal_file=$(cat normal_file.txt | head -n 1 | tail -1)
		
		gatk Mutect2 -R hg19.fa -I ${sample}.dedup.aligned.bam -I ${normal_file}/${normal_file}.dedup.aligned.bam -normal ${normal_file} -O ${sample}.GATK.somatic_unfiltered.dedup.vcf --native-pair-hmm-threads 30
		gatk FilterMutectCalls -R hg19.fa -V ${sample}.GATK.somatic_unfiltered.dedup.vcf -O ${sample}.GATK.somatic_filtered.dedup.vcf
		mkdir COSMIC_analyses/${sample}
     		bcftools filter -O z -o COSMIC_analyses/${sample}/final_${sample}.GATK.somatic_filtered.dedup.vcf -i '%FILTER="PASS"' ${sample}.GATK.somatic_filtered.dedup.vcf
		cd ..	
	fi
done





