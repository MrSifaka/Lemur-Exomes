import os

configfile: "LemurExomes.json"

# User paths
fastq_directory = "2017_Lemur_Exomes/fastq/"
temp_directory = "tmp/"
gatk = "" # path to jar file for GATK 3.7

# Set fastq and sample lists
all_fastq_prefixes = config["lemur_fastq_prefixes"]
prop_samples = config["Prop_males"] + config["Prop_females"]
eul_samples = config["Eulemur_males"] + config["Eulemur_females"]
all_samples = config["Eulemur_males"] + config["Eulemur_females"] + config["Prop_males"] + config["Prop_females"]

rule all:
	input:
		expand(
			"2017_Lemur_Exomes/fastq/fastqc/{fq_prefix}_fastqc.html",
			fq_prefix=all_fastq_prefixes),
		expand(
			"2017_Lemur_Exomes/bams/{sample}.Pcoq_1_0_genomic.mkdup.sorted.bam",
			sample=prop_samples),
		expand(
			"2017_Lemur_Exomes/bams/{sample}.Eflav_k33QCA_genomic.mkdup.sorted.bam",
			sample=eul_samples),
		expand("2017_Lemur_Exomes/vcf/genotyped_vcfs/filtered_vcfs/{genome}.filtered_allvariant_MP{MQ}.vcf.gz",
			MQ=["20", "30", "40"],
			genome=["Eflav_k33QCA_genomic", "Pcoq_1_0_genomic"])

rule prepare_reference:
	input:
		ref = lambda wildcards: config["genome_paths"][wildcards.assembly]
	output:
		new = "Seq_Data/References/{assembly}.fasta",
		fai = "Seq_Data/References/{assembly}.fasta.fai",
		amb = "Seq_Data/References/{assembly}.fasta.amb",
		dict = "Seq_Data/References/{assembly}.dict"
	run:
		shell(
			"ln -s ../{} {{output.new}} && touch -h {{output.new}}".format(input.ref))
		# faidx
		shell(
			"samtools faidx {output.new}")
		# .dict
		shell(
			"samtools dict -o {output.dict} {output.new}")
		# bwa
		shell(
			"bwa index {output.new}")


rule fastqc_analysis:
	input:
		os.path.join(fastq_directory, "{fq_prefix}_001.fastq.gz")
	output:
		"2017_Lemur_Exomes/fastq/fastqc/{fq_prefix}_fastqc.html"
	shell:
		"fastqc -o fastqc {input}"


rule map_and_process_trimmed_reads:
	input:
		fq1 = os.path.join(fastq_directory, "{sample}_R1_001.fastq.gz"),
		fq2 = os.path.join(fastq_directory, "{sample}_R2_001.fastq.gz"),
		fai = "Seq_Data/References/{genome}.fasta.fai",
		ref = "Seq_Data/References/{genome}.fasta"
	output:
		"2017_Lemur_Exomes/bams/{sample}.{genome}.mkdup.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"]
	threads: 4
	shell:
		" bwa mem -t {threads} -R "
		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| samblaster "
		"| samtools fixmate -O bam - - | samtools sort "
		"-O bam -o {output}"


rule index_bam:
	input:
		"2017_Lemur_Exomes/bams/{sample}.{genome}.mkdup.sorted.bam"
	output:
		"2017_Lemur_Exomes/bams/{sample}.{genome}.mkdup.sorted.bam.bai"
	shell:
		"samtools index {input}"


rule bam_stats:
	input:
		"2017_Lemur_Exomes/bams/{sample}.{genome}.mkdup.sorted.bam"
	output:
		"2017_Lemur_Exomes/stats/{sample}.{genome}.mkdup.sorted.bam.stats"
	shell:
		"samtools stats {input} | grep ^SN | cut -f 2- > {output}"


rule generate_callable_sites:
	input:
		ref = "Seq_Data/References/{genome}.fasta",
		bam = "2017_Lemur_Exomes/bams/{sample}.{genome}.mkdup.sorted.bam",
		bai = "2017_Lemur_Exomes/bams/{sample}.{genome}.mkdup.sorted.bam.bai"
	output:
		"2017_Lemur_Exomes/callable_sites/{sample}.{genome}.callablesites"
	params:
		temp_dir = "tmp/",
		gatk_path = gatk,
		summary = "2017_Lemur_Exomes/stats/{sample}.{genome}.callable.summary"
	shell:
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T CallableLoci -R {input.ref} -I {input.bam} --minDepth 4 --minMappingQuality 20 --minBaseQuality 30 --summary {params.summary} -o {output}"


rule extract_callable_sites:
	input:
		"2017_Lemur_Exomes/callable_sites/{sample}.{genome}.callablesites"
	output:
		"2017_Lemur_Exomes/callable_sites/{sample}.{genome}.ONLYcallablesites.bed"
	shell:
		"sed -e '/CALLABLE/!d' {input} > {output}"


rule combine_callable_sites_prop:
	input:
		expand(
			"2017_Lemur_Exomes/callable_sites/{sample}.Pcoq_1_0_genomic.ONLYcallablesites.bed",
			sample=prop_samples)
	output:
		"2017_Lemur_Exomes/callable_sites/combined.Pcoq_1_0_genomic.COMBINEDcallablesites.bed"
	shell:
		"cat {input} | sort -k1,1 -k2,2n | bedtools merge -i stdin > {output}"


rule combine_callable_sites_eul:
	input:
		expand(
			"2017_Lemur_Exomes/callable_sites/{sample}.Eflav_k33QCA_genomic.ONLYcallablesites.bed",
			sample=eul_samples)
	output:
		"2017_Lemur_Exomes/callable_sites/combined.Eflav_k33QCA_genomic.COMBINEDcallablesites.bed"
	shell:
		"cat {input} | sort -k1,1 -k2,2n | bedtools merge -i stdin > {output}"


rule gatk_gvcf_eul:
	input:
		ref = "Seq_Data/References/Eflav_k33QCA_genomic.fasta",
		bam = "2017_Lemur_Exomes/bams/{sample}.Eflav_k33QCA_genomic.mkdup.sorted.bam",
		bai = "2017_Lemur_Exomes/bams/{sample}.Eflav_k33QCA_genomic.mkdup.sorted.bam.bai",
		callable = "2017_Lemur_Exomes/callable_sites/combined.Eflav_k33QCA_genomic.COMBINEDcallablesites.bed"
	output:
		"2017_Lemur_Exomes/vcf/{sample}.Eflav_k33QCA_genomic.g.vcf.gz"
	params:
		temp_dir = "tmp/"
	threads:
		4
	shell:
		"""gatk --java-options "-Xmx15g" HaplotypeCaller"""
		""" -R {input.ref} -I {input.bam} -L {input.callable} """
		"""-ERC GVCF -O {output}"""


rule gatk_gvcf_prop:
	input:
		ref = "Seq_Data/References/Pcoq_1_0_genomic.fasta",
		bam = "2017_Lemur_Exomes/bams/{sample}.Pcoq_1_0_genomic.mkdup.sorted.bam",
		bai = "2017_Lemur_Exomes/bams/{sample}.Pcoq_1_0_genomic.mkdup.sorted.bam.bai",
		callable = "2017_Lemur_Exomes/callable_sites/combined.Pcoq_1_0_genomic.COMBINEDcallablesites.bed"
	output:
		"2017_Lemur_Exomes/vcf/{sample}.Pcoq_1_0_genomic.g.vcf.gz"
	params:
		temp_dir = "tmp/"
	threads:
		4
	shell:
		"""gatk --java-options "-Xmx15g" HaplotypeCaller"""
		""" -R {input.ref} -I {input.bam} -L {input.callable} """
		"""-ERC GVCF -O {output}"""


rule combine_gvcfs_prop:
    input:
        gvcfs_list = expand("2017_Lemur_Exomes/vcf/{sample}.Pcoq_1_0_genomic.g.vcf.gz", sample=prop_samples),
        ref = "Seq_Data/References/Pcoq_1_0_genomic.fna"
    output:
        "2017_Lemur_Exomes/vcf/combined_gvcfs/Pcoq_1_0_genomic.combined.g.vcf.gz"
    params:
        files = lambda wildcards, input: " -V ".join(input.gvcfs_list),
        temp_dir = "tmp/"
    log:
        "logs/gatk/combinegvcfs.log"
    shell:
        """ gatk --java-options "-Xmx16g" """
        """ CombineGVCFs -V {params.files} """
        """ -R {input.ref} """
        """ -O {output} """


rule combine_gvcfs_eul:
    input:
        gvcfs_list = expand("2017_Lemur_Exomes/vcf/{sample}.Eflav_k33QCA_genomic.g.vcf.gz", sample=eul_samples),
        ref = "Seq_Data/References/Pcoq_1_0_genomic.fna"
    output:
        "2017_Lemur_Exomes/vcf/combined_gvcfs/Eflav_k33QCA_genomic.combined.g.vcf.gz"
    params:
        files = lambda wildcards, input: " -V ".join(input.gvcfs_list),
        temp_dir = "tmp/"
    log:
        "logs/gatk/combinegvcfs.log"
    shell:
        """ gatk --java-options "-Xmx16g" """
        """ CombineGVCFs -V {params.files} """
        """ -R {input.ref} """
        """ -O {output} """


rule genotype_gvcfs:
    input:
        ref = "Seq_Data/References/{genome}.fna",
        gvcf = "2017_Lemur_Exomes/vcf/combined_gvcfs/{genome}.combined.g.vcf.gz",
        callable = "2017_Lemur_Exomes/callable_sites/combined.{genome}.COMBINEDcallablesites.bed"
    output:
        "2017_Lemur_Exomes/vcf/genotyped_vcfs/{genome}.genotyped.vcf.gz"
    params:
        temp_dir = "tmp/"
    shell:
        """ gatk --java-options "-Xmx16g" """
        """ GenotypeGVCFs -V {input.gvcf} -L {input.callable} """
        """ -R {input.ref} --include-non-variant-sites """
        """ -O {output} """


rule filter_vcfs:
	input:
		expand("2017_Lemur_Exomes/vcf/genotyped_vcfs/{genome}.genotyped.vcf.gz", genome=["Eflav_k33QCA_genomic", "Pcoq_1_0_genomic"])
	output:
		expand("2017_Lemur_Exomes/vcf/genotyped_vcfs/filtered_vcfs/{genome}.filtered_allvariant_MP{MQ}.vcf.gz",
			MQ=["20", "30", "40"],
			genome=["Eflav_k33QCA_genomic", "Pcoq_1_0_genomic"])
	params:
		bgzip = "~/anaconda3/envs/LemurEx17/bin/bgzip",
		bcftools = "~/anaconda3/envs/LemurEx17/bin/bcftools"
	shell:
		"{params.bcftools} filter -i "
		"'QUAL >= 30 && MQ >= {MQ} && QD > 2' {input} | "
		"{params.bcftools} filter -i 'FMT/DP >= 10 & FMT/GQ >= 30' -S . - | "
		"{params.bgzip} > {output}"
