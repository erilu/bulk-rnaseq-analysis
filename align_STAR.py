import subprocess
import os

# modify these variables based on where you store your data / genome index
output_directory = "STAR_output/"
genome_directory = "/data/genomes/h38/STAR/"
fastq_directory = "/data/analysis/hypoxia/fastq/"

os.mkdir(output_directory)
for fastq in os.listdir(fastq_directory):
    # only process files that end in fastq.gz
    if fastq.endswith('.fastq.gz'):
        prefix=fastq.strip(".fastq.gz") + "_output"
        # make an output folder for the current fastq file
        os.mkdir(output_directory + prefix)
        print ("Currently mapping: " + fastq)
        # run STAR on the current fastq file
        subprocess.call("STAR --runThreadN 64 --genomeDir " + genome_directory + " --readFilesCommand zcat --outFilterType BySJout --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterMultimapNmax 20 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesIn "+ fastq_directory + fastq + " --clip3pAdapterSeq GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix " + output_directory + prefix + "/", shell=True)
