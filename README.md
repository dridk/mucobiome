# Mucobiome
Mucobiome is a simple pipeline which aims to analyse 16S RNA genomics data from high throughput sequencer.
Neither QIIME nor MOTHUR are required. Few dependency like vsearch are necessary and are listed bellow.
Mucobiome use [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home) as the backbone. This tools make it possible to run
the pipeline optimally using multithreading. 

## How it works ? **Fastq** -> **biom**   
  
the pipeline takes several pair-end fastq as input ( one per sample) and generate one biom file which contains OTU table with taxonomy and sample 
meta data. 
 - For each reads pairs :
  - Merge fastq file using vsearch or flash 
  - Clean fastq file using sickle 
  - Reverse reads with seqtk  
  - Trim adaptators from reads using cutadapts 
  - Dereplicate reads with vsearch 
  - merge all reads pairs into one fasta file 
 - With Greengene database 
  - Extract interesting region from greengene 16S database using cutadapts and user's adaptator 
 - Make a taxonomy assignement using vsearch --usearch_global
  - Compare sequence from the merged file and greengene database
  - Create a Biom file
  - Add taxonomy and sample metadata into the biom file


## Installation 
### Python depedencies 
Mucobiome has been written with Python 3.4. 

    pip install -r requirements.txt 

### Install vsearch
```
wget https://github.com/torognes/vsearch/archive/v2.3.4.tar.gz
tar xzf v2.3.4.tar.gz
cd vsearch-2.3.4
./autogen.sh
./configure
make
make install  # as root or sudo make install
```

### Install seqtk
```
git clone https://github.com/lh3/seqtk.git;
cd seqtk;
make
sudo make install
```

### Install sickle 
```
https://github.com/najoshi/sickle.git
cd sickle
make
sudo make install
```

## Download Database 
Mucobiome works with greengene. But you can use another database if you respect the same format. 
Run download_greengene.sh from the database folder to download greengene data. 

     cd database; sh download_greengene.sh 
     
## usage 
### Test your installation 
The actual repositories contains a simple dataset. Try the following commands which do nothing but display all pipeline command.

```
# you are in the main directory 
snakemake -d working_directory -np --configfile config.yaml
```
<table>
<tr><th>Options</th><th>Description</th></tr>
<tr><th>-d</th><td>The working directory. All generated files will be drop here</td></tr>
<tr><th>-n</th><td>Don't execute any commands.</td></tr>
<tr><th>-p</th><td>Print commands</td></tr>
<tr><th>--configfile</th><td>Tell which config file to use</td></tr>
</table>

### Input file 
input data are paired fastq.gz files. You must put all your datas into the data/raw folder. Both paired files must respect the following syntax.
*{Sample}* is your samplename.

     {SAMPLE}_1.fastq.gz 
     {SAMPLE}_2.fastq.gz

### config.yaml 
This file contains all parameters required to perform an analysis. 

<table>
<tr><th>Options</th><th>Description</th></tr>
<tr><th>raw_folder</th><td>the directory which contains fastq input files</td></tr>
<tr><th>primer_forward</th><td>The forward primer. By defaut primers select the V3-V5 region</td></tr>
<tr><th>reverse_primer</th><td>The reverse primer. By defaut primers select the V3-V5 region</td></tr>
<tr><th>database_fasta</th><td>This is a fasta file which contains complete 16S sequence. By default it use greengene</td></tr>
<tr><th>database_taxonomy</th><td>This is a two columns file. Fasta header ID from database_fasta and the taxonomy</td></tr>
<tr><th>sample_data</th><td>Sample meta data</td></tr>
<tr><th>threshold</th><td>Taxonomy assignement threshold. By default 97%</td></tr>
<tr><th>quality</th><td>Remove all reads bellow this threshold. You should set this value to 20</td></tr>
<tr><th>min_len</th><td>Minimum reads length accepted</td></tr>
<tr><th>max_len</th><td>Maximum reads length accepted</td></tr>
<tr><th>merge_tool</th><td>Use Flash or vsearch to perform merging. By default it use vsearch</td></tr>
</table>

### sample_data.tsv
This files contains metadata for each sample. Fill it with your own data. 


### Run your experiment 
When all requirements are done, you can run one of the following commands

       # Run your pipeline using 60 threads
       snakemake -d working_directory -p --configfile config.yaml --cores 60
       
       # Force the pipeline to rebuild the last step
       snakemake -d working_directory -fp --configfile config.yaml --cores 60
       
       # Force the pipeline to rebuild everything
       snakemake -d working_directory -Fp --configfile config.yaml --cores 60


     







