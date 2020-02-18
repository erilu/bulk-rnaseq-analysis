
import os
import csv

# modify desired_column based on strandedness of library kit
desired_column = 3
output_directory = "STAR_output/"
sample_dict = {}
sample_names = []

# loop through each folder and store data from ReadsPerGene.out.tab in sample_dict
for folder in os.listdir(output_directory):
    if folder.endswith("_output"):
        file_name=folder.strip('_output')
        sample_names.append(file_name)
        gene_dict = {}
        with open(output_directory + folder + "/ReadsPerGene.out.tab") as tabfile:
            print("Column " + str(desired_column) + " of " + file_name + " stored")
            reader = csv.reader(tabfile, delimiter = "\t")
            for row in reader:
                gene_dict[row[0]]=row[desired_column]
        # store gene_dict for the current file in sample_dict
        sample_dict[file_name] = gene_dict

# write sample_dict to output files
# the qc metrics are stored in qc.csv, and the gene counts are stored in raw_counts.csv
with open("raw_counts.csv", "wt") as counts_file, open("qc.csv", "wt") as qc_file:
    counts_writer = csv.writer(counts_file)
    qc_writer = csv.writer(qc_file)
    counts_writer.writerow( ['ensembl_id']+ sample_names )
    qc_writer.writerow( ['qc_metric']+ sample_names )
    sorted_genes = sorted(list(sample_dict[sample_names[0]].keys()))
    for gene in sorted_genes:
        output=[gene]
        for sample in sample_names:
            output.append(sample_dict[sample][gene])
        if gene.startswith("N_"):
            qc_writer.writerow(output)
        else:
            counts_writer.writerow(output)
