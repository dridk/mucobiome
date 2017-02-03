from Bio import SeqIO
import random
import yaml 

print("generate sim reads ")

def make_mutation(seq , proba = 0.1):
	count     = int(len(seq) * proba)
	nucleotid = ["A","C","G","T"] 
	seq_list  = list(seq)
	for index in random.sample(range(len(seq)), count):
		seq_list[index] = random.choice(nucleotid)

	return "".join(seq_list) 



CONFIG_PATH = "sim_reads.yaml"

with open(CONFIG_PATH, "r") as config:
	yaml  = yaml.load(config)
	TOTAL = int(yaml["total"])
	SOURCE = yaml["source"]

	for sample in yaml["samples"]:
		print("simulate ", sample["name"])
		OUTPUT    = "sim_"+sample["name"] + ".fastq"
		MUT       = sample["proba"]
		ABUNDANCE = sample["species"]
		FACTOR    = sample["factor"]

		sum_abundunce = sum(ABUNDANCE.values())

		for key in ABUNDANCE.keys():
			ABUNDANCE[key] = ABUNDANCE[key] / sum_abundunce

		with open(OUTPUT,"w") as fastq:
			with open(SOURCE,"r") as fasta:
				records = SeqIO.parse(fasta, "fasta")
				for record in records:
					if int(record.name) in ABUNDANCE.keys():
						count = ABUNDANCE[int(record.name)] * TOTAL;
						record.letter_annotations["phred_quality"]  = [40] * len(record)
						record.seq = make_mutation(record.seq, MUT)
						for i in range(int(count) * FACTOR ):
							SeqIO.write(record, fastq, "fastq")	




			
	




