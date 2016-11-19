import glob 
import re 

rule merge_all :
	input: 
		set([m.group(0)+".rename.fasta" for m in (re.search("\d+",l) for l in glob.glob("raw/*.fastq.gz")) if m is not None])

	output:
		"all.merge.fasta"
	shell:
		"cat {input} > {output}"

rule remove_chimer : 
	input:
		"all.merge.fasta"
	output:
		"all.unchimera.fasta"
	log:
		"all.unchimera.log"
	shell:
		"vsearch --uchime_denovo {input} --nonchimeras {output} 2> {log}"



rule assign_taxonomy : 
	input:
		reads = "all.unchimera.fasta",
		reference = config['reference']
	output:
		otu  = "otu.txt"
	shell:
		"vsearch --usearch_global {input.reads} --db {input.reference} --id 0.97 --sizein --threads 40 --otutabout {output.otu}"




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
