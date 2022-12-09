import json
from collections import defaultdict
import goseqpy as gsp
import sys
import os
import cmdlogtime

COMMAND_LINE_DEF_FILE = "./runGoSeqCommandLine.txt"
def main():
    (start_time_secs, pretty_start_time, my_args, logfile) = cmdlogtime.begin(COMMAND_LINE_DEF_FILE)   
    in_file = my_args["in_file"]
    out_file = os.path.join(my_args["out_dir"], "output.tsv")
    final_json = os.path.join(my_args["out_dir"], "json_final.json")
    has_no_header = my_args["has_no_header_line"]
    p_value_col = my_args["p_value_col"]
    p_val_cutoff = my_args["p_value"]
    log2_fc_col = my_args["log2_fc_col"]
    log2_fc_cutoff = my_args["log2_fc"]
    species = my_args["species"]
    geneset_file = my_args["geneset_file"]
    id_col = int(my_args["id_col"])
    # Load example data
    genes_to_write = []
    json_file = "./json_file.json"
    if (in_file.endswith(".json")):
        json_file = in_file
    else:
        with open(in_file, "r") as in_gene_file: #, open('example_data_template.json') as json_template:  
            line_ctr = 0
            num_lines = len(in_gene_file.readlines())
            in_gene_file.seek(0) #reset the file back to the beginning
            for line in in_gene_file:
                line_ctr = line_ctr + 1
                if (has_no_header == False and line_ctr == 1):
                    continue 
                parts = line.strip().split("\t")
                gene = parts[id_col]
                if (p_value_col):
                    p_value_col = int(p_value_col)
                    p_val = float(parts[p_value_col])
                    p_val_cutoff = float(p_val_cutoff)
                    if (p_val >= p_val_cutoff):
                        continue
                if (log2_fc_col):
                    log2_fc_col = int(log2_fc_col)
                    log2_fc = float(parts[log2_fc_col])
                    log2_fc_cutoff = float(log2_fc_cutoff)
                    if (abs(log2_fc) < log2_fc_cutoff):
                        continue        
                genes_to_write.append(gene)
        num_genes_to_write = len(genes_to_write)
        
        if (species == "human"):
            json_template_file = './data/example_data_template.json'
        elif (species == "drosophila_flydb"):
            #json_template_file = 'genes_and_gene_lengths_symbol_drosophila.json'  #RMS make this file a parameter.
            json_template_file = '../data/genes_and_gene_lengths_flydb_drosophila.json'
        elif (species == "drosophila_symbol"):
            json_template_file = '../data/genes_and_gene_lengths_symbol_drosophila.json'
        else:
            print("Species: ", species, " not yet supported.")
            sys.exit(1)
            
        with open(json_template_file) as json_template, open(json_file, 'w') as out_json_file:
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
    if (species == "human"):
        gs_to_genes = gsp.parse_gmt('../data/gene_sets/c5.bp.v7.1.symbols.gmt')  #all #RMS. make this file a parameter
    elif (species == "drosophila_flydb"):
        #gs_to_genes = gsp.parse_gmt('./gene_sets/drosophila_genesets_symbol4.txt')  #RMS. make this file a parameter
        gs_to_genes = gsp.parse_gmt('../data/gene_sets/drosophila_genesets_flydb.txt')  #all
    elif (species == "drosophila_symbol"):
        #gs_to_genes = gsp.parse_gmt('./gene_sets/drosophila_genesets_symbol4.txt')  #RMS make this file a parameter
        gs_to_genes = gsp.parse_gmt('../data/gene_sets/drosophila_genesets_symbol.txt')  #all
    else:
        print("Species: ", species, " not yet supported.")
        sys.exit(1)
    if (not geneset_file.endswith("ALLZZZ")): 
        gs_to_genes = gsp.parse_gmt(geneset_file)
    #gs_to_genes = gsp.parse_gmt('./gene_sets/c5.bp.v7.1.symbols_leuk_only.gmt')  #leuk only

    # Map each gene to its gene set
    gene_to_gene_sets = defaultdict(lambda: [])
    for gene_set, genes in gs_to_genes.items():
        for gene in genes:
            gene_to_gene_sets[gene].append(gene_set)    
    gene_to_gene_sets = dict(gene_to_gene_sets)
    
    try: 
        # Run GOSeq
        results = gsp.run_GOseq(
            de_genes, 
            all_genes, 
            gene_to_gene_sets, 
            gene_lens,
            my_args["out_dir"]
        )
        results.to_csv(out_file, sep='\t')
    except Exception as e:
        print("Error: ", e)
        if (e.__str__().endswith("incorrect number of dimensions\n")):
            print("Might  be a problem with using -spec=drosophila_flydb and the ids used in the Gene Set file are symbols. Check parameters")
        else:
            print("Some error. Not the dimension error.")
    
    cmdlogtime.end(logfile, start_time_secs) 

if __name__ == "__main__":
    main()
