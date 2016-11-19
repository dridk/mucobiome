import glob 
import re 

rule merge_all :
	input: 
		[m.group(0)+".rename.fasta" for m in (re.search("\d+",l) for l in glob.glob("raw/*.fastq.gz")) if m is not None]

	output:
		"all.merge.fasta"
	shell:
		"cat {input} > {output}"

rule remove_chimer : 
	input:
		"all.merge.fasta"
	output:
		"all.unchimera.fasta"
	shell:
		"vsearch --uchime_denovo {input} --nonchimeras {output}"



rule assign_taxonomy : 
	input:
		reads = "all.merge.fasta",
		greengene = "greengene.trim.fasta"
	output:
		biom = "raw.biom",
		otu  = "otu.txt"
	shell:
		"vsearch --usearch_global {input.reads} --db {input.greengene} --id 0.90 --sizein --threads 40 --biomout {output.biom} --otutabout {output.otu}"


rule add_metadata:
	input:
		biom="raw.biom",
		col="greengene.taxonomy.txt"
	output:
		"final.biom"
	shell:
		"biom add-metadata -i {input.biom} --observation-metadata-fp {input.col} -m mapping.txt -o {output} --sc-separated taxonomy"


rule create_taxonomy_for_biom:
	input:
		"greengene.trim.fasta"
	output:
		"greengene.taxonomy.txt"
	shell:
		"echo -e '#OTU_ID\ttaxonomy'>{output};"
		"cat {input}|grep '>'|sed 's/>//g'|cut -f1,4 >>{output}"





rule merge :
	input:
		forward = "raw/{sample}_1.fastq.gz",	
		reverse = "raw/{sample}_2.fastq.gz"
	output:
		"{sample}.merged.fastq"
	log:
		"{sample}.merged.log"
	shell: 
		"vsearch --fastq_mergepairs {input.forward} --reverse {input.reverse} --fastqout {output} 2> {log}"





rule clean : 
	input:
		"{sample}.merged.fastq"
	output:
		"{sample}.clean.fastq"
	log:
		"{sample}.clean.log"
	shell:
		"sickle se -q 30 -l 500  -f {input} -t sanger -o {output}  > {log}"


rule reverse:
	input:
		"{sample}.clean.fastq"
	output:
		"{sample}.fasta"
	shell:
		"seqtk seq -A -r {input} > {output}"



rule trim_primers:
	input:
		"{sample}.fasta"
	output:
		"{sample}.trim.fasta"
	log:
		forward = "{sample}.cutadapt_forward",
		reverse = "{sample}.cutadapt_reverse"
		
	shell:
		"cutadapt --discard-untrimmed -g {config[primer_forward]} {input} 2> {log.forward}|"
		"cutadapt --discard-untrimmed -a $(echo {config[primer_reverse]}|rev|tr ATGCRM TACGYK) - 2> {log.reverse} > {output}"



rule dereplicate : 
	input:
		"{sample}.trim.fasta"
	output:
		"{sample}.dereplicate.fasta"
	log:
		"{sample}.dereplicate.log"

	shell:
		"vsearch --derep_fulllength {input} --output {output} --sizeout --relabel_keep 2> {log}"

rule rename:
	input:
		"{sample}.dereplicate.fasta"
	output:
		"{sample}.rename.fasta"
	shell:
		"cat {input}|sed -e 's/>.*/&sample={wildcards.sample};/g' > {output}"
