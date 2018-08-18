import pandas as pd
import numpy as np
import biom
import yaml
import re
from Bio import SeqIO
import sys

PATH      = sys.argv[1]
BIOM_FILE = sys.argv[2]

biom_file = open(PATH+"/"+BIOM_FILE, "w")



# Load assignement as a dict 
# contains : 16sRNA Hash sequence => Taxon ID 
assignement = {}
with open(PATH+"/all.assignement.txt") as file:
    for line in file:
        row = line.rstrip().split("\t")
        assignement[row[0]] = row[1]


# Load samples from config file 
config = yaml.load(open(PATH+"/config.yaml"))


# Create OTU Matrix as a pandas data frame 
df = pd.DataFrame(columns=np.unique(list(assignement.values())))

# Loop over all sample, extract Sequence Hash and read count (size=)
for sample in config["samples"]:
    file = PATH+"/{}.dereplicate.fasta".format(sample)
    
    sample_dict = {}

    
    for record in SeqIO.parse(file, "fasta"):
        # Extract size from id 
        # 7316bc425bf317aa495fc6d1ec2a69bc25440ba0;size=1768;
        
        size_regexp = re.search('(\w+);size=(\d+)', record.id, re.IGNORECASE)
        uuid   = size_regexp.group(1)
        size   = int(size_regexp.group(2))
        
        try:
            taxaid = assignement[uuid]
            sample_dict[taxaid] = [size]
        except:
            sample_dict[taxaid] = [np.NaN]
        
        
    df = df.append(pd.DataFrame(sample_dict, index=[str(sample)]))

# Replace NaN by 0       
df = df.replace(np.nan, 0)

df = df.transpose()

# Create a biom data Table and write as a file 
table = biom.table.Table(df.values, observation_ids=df.index, sample_ids=df.columns)

biom_file.write(table.to_json("mucobiome"))
biom_file.close()