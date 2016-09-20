# merge all fasta file into 1 
rule combine_all:
	input: 
		"1001.fasta","1080.fasta", "1126.fasta", "1015.fasta", "1029.fasta", "1036.fasta", "1042.fasta" 
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

# Convert fastq to fasta and use sample name as fasta header 
rule to_fasta:
	input:
		"{sample}.extendedFrags.fastq.gz"
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
	shell:
		"vsearch -threads {threads} --cluster_fast {input} --id 0.97 --centroids {output.centroid} --uc {output.cluster} --relabel_keep --relabel_sha1" 


rule create_otu:
	input:
		"centroid.fasta"
	output:
		"otu.txt" 
	shell:
		"cat centroid.fasta|grep -oE '>.+\s'|sed 's/>//g'|awk -f genotu.awk > {output}"



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


