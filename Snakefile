import glob 
import re 
# merge all fasta file into 1 
rule combine_all:
	input:  [m.group(0)+".fasta" for m in (re.search("\d+",l) for l in glob.glob("raw/*.fastq.gz")) if m is not None]

	output:
		"all.merge.fasta" 
	shell:
		"cat {input} > {output}"



# Execute Flash to merge pair end 
rule merge:
	input:
		forward="raw/{sample}_1.fastq.gz",
		reverse="raw/{sample}_2.fastq.gz"
	output:
		"{sample}.extendedFrags.fastq.gz"
	shell:
		"flash {input.forward} {input.reverse} -z -o {wildcards.sample}"


# Clean sequences 
rule clean_seq:
	input:
		"{sample}.extendedFrags.fastq.gz"
	output:
		"{sample}.clean.fastq.gz"
	shell:
		"sickle se -f {input} -t sanger -o {output} -g -q 35"


# Convert fastq to fasta and use sample name as fasta header 
rule to_fasta:
	input:
		"{sample}.clean.fastq.gz"
	output:
		"{sample}.fasta" 
	shell:
		"seqtk rename {input} '{wildcards.sample}_'|seqtk seq -A > {output}"

# Clusting fasta file using vsearch 
rule clustering:
	input:
		"all.merge.fasta"
	output:
		cluster= "cluster.uc",
		centroid= "centroid.fasta"	
	
	threads : 64 
	
	shell:
		"vsearch -threads {threads} --cluster_fast {input} --id 0.97 --centroids {output.centroid} --uc {output.cluster} --relabel_keep --relabel_sha1 --sizeout" 


#rule create_otu:
#	input:
#		"centroid.fasta"
#	output:
#		"otu.txt" 
#	shell:
#		"cat centroid.fasta|grep -oE '>.+\s'|sed 's/>//g'|awk -f genotu.awk > {output}"



rule create_biom:
	input:
		cluster="cluster.uc",
		centroid="centroid.fasta"
	output:
		"raw.biom"
	shell:
		"biom from-uc -i {input.cluster} -o {output} --rep-set-fp {input.centroid}"


rule add_metadata:
	input:
		row = "map.txt",
		col = "otu.txt",
		raw = "raw.biom"
	output:
		"final.biom" 
	shell:
		"biom add-metadata -i {input.raw} --observation-metadata-fp {input.col} -m {input.row} -o {output} --sc-separated taxonomy"



rule assign_taxonomy:
	input :
		 "centroid.fasta"
	output:
		"classified.txt"

	shell:
		"java -Xmx10g -jar /PROGS/EXTERN/RDPTools/classifier.jar  classify -c 0.5 -o {output} -h soil.txt  centroid.fasta"

rule parse_classifed:
	input :
		"classified.txt"
	output:
		"otu.txt"
	shell:
		"cat classified.txt|awk 'BEGIN{{OFS=\"\t\";print(\"#OTUSID\",\"taxonomy\")}} {{print($1,$3\";k__\"$6\";p__\"$9\";c__\"$12\";o__\"$15\";f__\"$18)}}' > {output}"



	

