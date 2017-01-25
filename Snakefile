import glob 
import re 


# rule downlad_greengene:
#     input:
#         FTP.remote("greengenes.microbio.me/greengenes_release/current/gg_13_5.fasta.gz",keep_local=True,username="anonymous", password="")
#     output:
#     	"greengene.fasta"
#     shell:
#     	"mv {input} {output}.gz; gzip -d {output}.gz"


# Generate Biom file . Final point of the pipeline
rule final:
	input: "final.biom"


# Merge all fasta file after all preprocess has been done
rule merge_all :
	input: 
		set([m.group(0)+".rename.fasta" for m in (re.search("\d+",l) for l in glob.glob("{}/*.fastq.gz".format(config["raw_folder"]))) if m is not None])

	output:
		"all.merge.fasta"
	shell:
		"cat {input} > {output}"

# Remove chimera
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


# Taxonomy assignement 
rule assign_taxonomy : 
	input:
		reads           = "all.merge.fasta",
		target_database = "target_database.trim.fasta"
	output:
		biom = "raw.biom",
		otu  = "otu.txt"

	log:
		"assign_taxonomy.log"
	threads:
		60

	shell:
		"vsearch --usearch_global {input.reads} --db {input.target_database} --id {config[threshold]} --sizein --threads {threads} --biomout {output.biom} --otutabout {output.otu} 2> {log}"

# Add metadata : Column are replaced by taxon name and row by sample data 
rule add_metadata:
	input:
		biom="raw.biom",
		col ="target_taxonomy.txt",
		row =config["sample_data"]
	output:
		"final.biom"
	shell:
		"biom add-metadata -i {input.biom} --observation-metadata-fp {input.col} -m {input.row} -o {output} --sc-separated taxonomy"

# Create 2 columns file with Taxon ID and Taxon Name as header
rule create_taxonomy_for_biom:
	input:
		config["database_taxonomy"]
	output:
		"target_taxonomy.txt"
	shell:
		"zcat {input}|sed '1i#OTU_ID\ttaxonomy' > {output}"


# Merge pair end file with vsearch
rule merge_with_vsearch :
	input:
		forward = config["raw_folder"]+"/{sample}_1.fastq.gz",	
		reverse = config["raw_folder"]+"/{sample}_2.fastq.gz"

	output:
		file ="{sample}.merged_vsearch.fastq",  #Vsearch log
		log1 ="{sample}.merged_vsearch.log",    # Standard log 
		log2 ="{sample}.merged.log" 
	shell: 
		"vsearch --fastq_mergepairs {input.forward} --reverse {input.reverse} --fastqout {output.file} 2> {output.log1};"
		"cat {output.log1}|grep -A2 'Merging'|tail -n 2|grep -oE '[0-9]+ '|tr -d ' ' > {output.log2};"



# Merge pair end file with flash
rule merge_with_flash :
	input:
		forward = config["raw_folder"]+"/{sample}_1.fastq.gz",	
		reverse = config["raw_folder"]+"/{sample}_2.fastq.gz"

	output:
		file  ="{sample}.merged_flash.fastq",
		log1  ="{sample}.merged_flash.log",
		log2  ="{sample}.merged.log"

	shell: 
		"flash {input.forward} {input.reverse} -o {wildcards.sample}  >{output.log1}; mv {wildcards.sample}.extendedFrags.fastq {output.file};"
		"cat {output.log1}|grep -A2 \"Read combination\"|grep -oE '[0-9]*' > {output.log2}"

rule clean : 
	input:
		"{sample}.merged_"+config["merge_tool"] +".fastq"
	output:
		file ="{sample}.clean.fastq",
		log  ="{sample}.clean.log"

	message:
		"Clean .. "

	shell:
		"sickle se -q {config[quality]} -l 500  -f {input} -t sanger -o {output.file} > {output.log}"



rule merge_clean_report:
	input:
		merge  ="{sample}.merged_"+config["merge_tool"]+".log",
		clean  ="{sample}.clean.log",
		trim_f ="{sample}.cutadapt_forward.log",
		trim_r ="{sample}.cutadapt_reverse.log"
	output:
		final="{sample}.resume.log",
		mid=temp("{sample}.merge_clean.temp"),
		mid2=temp("{sample}.merge_clean.temp2"),
		

	message:
		"parse log for {wildcards.sample}"
	shell:
		"cat {wildcards.sample}.merged.log > {output.mid};" 
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

rule uncompress_database:
	input:
		config["database_fasta"]
	output:
		"target_database.fasta"
	shell:
		"gzip -d < {input} > {output}"
	message:
		"Uncompress target database"

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
