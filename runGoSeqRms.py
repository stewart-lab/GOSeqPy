import json
from collections import defaultdict

import goseqpy as gsp

import sys
import os
import rmstime
import rmscmdline
import rmslogging
COMMAND_LINE_DEF_FILE = "./runGoSeqCommandLine.txt"

def main():
    (start_time_secs, pretty_start_time) = rmstime.get_time_and_pretty_time()
    print("pretty_start:", pretty_start_time)
    
    my_args = rmscmdline.get_args(start_time_secs, pretty_start_time, COMMAND_LINE_DEF_FILE)
    logfile = rmslogging.open_log_file(my_args["log_file"])
    rmslogging.write_args_and_files(my_args, sys.argv[0], COMMAND_LINE_DEF_FILE, logfile)
    
    in_file = my_args["in_file"]
    out_file = os.path.join(my_args["out_dir"], "output.tsv")
    final_json = os.path.join(my_args["out_dir"], "json_final.json")
    has_no_header = my_args["has_no_header_line"]
    print(has_no_header)
    p_value_col = my_args["p_value_col"]
    p_val_cutoff = my_args["p_value"]
    log2_fc_col = my_args["log2_fc_col"]
    log2_fc_cutoff = my_args["log2_fc"]
    # Load example data
    #with open('example_data.json', 'r') as f:
    genes_to_write = []
    json_file = "./json_file.json"
    if (in_file.endswith(".json")):
        json_file = in_file
    else:
        print("here need to process input file and write new json file")
        with open(in_file, "r") as in_gene_file: #, open('example_data_template.json') as json_template:  
            line_ctr = 0
            num_lines = len(in_gene_file.readlines())
            in_gene_file.seek(0) #reset the file back to the beginning
            for line in in_gene_file:
                line_ctr = line_ctr + 1
                if (has_no_header == False and line_ctr == 1):
                    continue 
                parts = line.strip().split("\t")
                gene = parts[0]
                if (p_value_col):
                    p_value_col = int(p_value_col)
                    p_val = float(parts[p_value_col])
                    p_val_cutoff = float(p_val_cutoff)
                    if (p_val >= p_val_cutoff):
                        print ("continuingP:", p_val,  " cutoff:", p_val_cutoff)
                        continue
                if (log2_fc_col):
                    log2_fc_col = int(log2_fc_col)
                    log2_fc = float(parts[log2_fc_col])
                    log2_fc_cutoff = float(log2_fc_cutoff)
                    if (abs(log2_fc) < log2_fc_cutoff):
                        print ("continuingFC:", log2_fc,  " cutoff:", log2_fc_cutoff)
                        continue        
                #maybe_comma = ","
                #if line_ctr == num_lines:
                #    maybe_comma = ""
                genes_to_write.append(gene)
                #out_json_file.write('"' + gene + '"' + maybe_comma + '\n')
            #out_json_file.write(" ]\n}\n")    
        num_genes_to_write = len(genes_to_write)
        with open('example_data_template.json') as json_template, open(json_file, 'w') as out_json_file:
            for line in json_template:
                out_json_file.write(line)
            gene_ctr = 0
            for gene in genes_to_write:
                gene_ctr = gene_ctr + 1
                maybe_comma = ","
                if gene_ctr == num_genes_to_write:
                    maybe_comma = ""
                out_json_file.write('"' + gene + '"' + maybe_comma + '\n')
            out_json_file.write(" ]\n}\n")       
    with open(json_file, 'r') as f, open(final_json, 'w') as final_json_out_file:
        example_data = json.load(f)
        f.seek(0)
        for line in f:
            final_json_out_file.write(line)
    all_genes = example_data['all_genes']
    gene_lens = example_data['gene_lengths']
    de_genes = example_data['de_genes']

    # Load gene set to genes mapping
    gs_to_genes = gsp.parse_gmt('./gene_sets/c5.bp.v7.1.symbols.gmt')

    # Map each gene to its gene set
    gene_to_gene_sets = defaultdict(lambda: [])
    for gene_set, genes in gs_to_genes.items():
        for gene in genes:
            gene_to_gene_sets[gene].append(gene_set)    
    gene_to_gene_sets = dict(gene_to_gene_sets)

    # Run GOSeq
    results = gsp.run_GOseq(
        de_genes, 
        all_genes, 
        gene_to_gene_sets, 
        gene_lens
    )

    # Show the results
    #print(results)
    results.to_csv(out_file, sep='\t')

    rmslogging.close_log_file(logfile)  
    (end_time_secs, x) = rmstime.get_time_and_pretty_time()
    total_elapsed_time = end_time_secs - start_time_secs
    print("All done. Total elapsed time: " + str(total_elapsed_time) + " seconds.\n")      

if __name__ == "__main__":
    main()