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
	threads:
		60
	shell:
		"vsearch --uchime_denovo {input} --nonchimeras {output} --threads {threads} 2> {log}"



rule assign_taxonomy : 
	input:
		reads = "all.merge.fasta",
		greengene = "greengene.trim.fasta"
	output:
		biom = "raw.biom",
		otu  = "otu.txt"

	log:
		"assign_taxonomy.log"
	threads:
		60

	shell:
		"vsearch --usearch_global {input.reads} --db {input.greengene} --id {config[threshold]} --sizein --threads {threads} --biomout {output.biom} --otutabout {output.otu} 2> {log}"


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
		file="{sample}.merged.fastq",
		log="{sample}.merged.log"
	shell: 
		"vsearch --fastq_mergepairs {input.forward} --reverse {input.reverse} --fastqout {output.file} 2> {output.log}"





rule clean : 
	input:
		"{sample}.merged.fastq"
	output:
		file="{sample}.clean.fastq",
		log="{sample}.clean.log"

	message:
		"Clean .. "

	shell:
		"sickle se -q {config[quality]} -l 500  -f {input} -t sanger -o {output.file} > {output.log}"



rule merge_clean_report:
	input:
		merge  ="{sample}.merged.log",
		clean  ="{sample}.clean.log",
		trim_f ="{sample}.cutadapt_forward.log",
		trim_r ="{sample}.cutadapt_reverse.log"
	output:
		final="{sample}.merge_clean.report",
		mid=temp("{sample}.merge_clean.temp"),
		mid2=temp("{sample}.merge_clean.temp2"),
		

	message:
		"parse log for {wildcards.sample}"
	shell:
		"cat {input.merge}|grep -A2 'Merging'|tail -n 2|grep -oE '[0-9]+ '|tr -d ' ' > {output.mid};" 
		"cat {input.clean}|grep 'kept'|grep -Eo '[0-9]+' >> {output.mid};" 
		"cat {input.trim_f}|grep 'Reads written'|grep -Eo '[0-9,]+ '|tr -d ',' >> {output.mid};"
		"cat {input.trim_r}|grep 'Reads written'|grep -Eo '[0-9,]+ '|tr -d ',' >> {output.mid};"
		"cat {output.mid}|awk -v OFS='\t' 'NR==1{{FINAL=$1; print $0, 100}} NR!=1{{print $1,$1*100/FINAL}}' > {output.mid2};"
		"echo 'raw\nmerge\nclean\nforward\nreverse'|paste - {output.mid2} > {output.final}" 	

	



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
		file    ="{sample}.trim.fasta",
		forward ="{sample}.cutadapt_forward.log",
		reverse ="{sample}.cutadapt_reverse.log"
		
	shell:
		"cutadapt --discard-untrimmed -g {config[primer_forward]} {input} 2> {output.forward}|"
		"cutadapt --discard-untrimmed -a $(echo {config[primer_reverse]}|rev|tr ATGCRM TACGYK) - 2> {output.reverse} > {output.file}"



rule dereplicate : 
	input:
		"{sample}.trim.fasta"
	output:
		"{sample}.dereplicate.fasta"
	shell:
		"vsearch --derep_fulllength {input} --output {output} --sizeout --relabel_keep"

rule rename:
	input:
		"{sample}.dereplicate.fasta"
	output:
		"{sample}.rename.fasta"
	shell:
		"cat {input}|sed -e 's/>.*/&sample={wildcards.sample};/g' > {output}"
