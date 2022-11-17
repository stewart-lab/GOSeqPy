import cmdlogtime
#import stew_util as su
import os
from collections import defaultdict
from pathlib import Path as ph

COMMAND_LINE_DEF_FILE = "./fix_go_descriptions_from_goseq_commandLine.txt"
def main():
    (start_time_secs, pretty_start_time, my_args, addl_logfile)= cmdlogtime.begin(COMMAND_LINE_DEF_FILE)   
    out_goseq_fixed = os.path.join(my_args["out_dir"], "out_goseq_fixed.tsv")    
    go_descs, go_namespaces = read_go_desc_files(my_args['in_go_desc1'], my_args['in_go_desc_suppl'], my_args['go_subset'])
    # write out fixed goseq file
    print_goseq_fixed(out_goseq_fixed, my_args['in_goseq'], go_descs, go_namespaces)    
    #import pdb
    #pdb.set_trace()
    cmdlogtime.end(addl_logfile, start_time_secs) 

# ----------------------------------- FUNCTIONS  -------------------------------------
def read_go_desc_files(in_go_desc1, in_go_desc_suppl, go_subset):
    go_descs = defaultdict()
    go_namespaces = defaultdict()
    go_descs, go_namespaces = process_desc_lines(in_go_desc1, go_descs, go_subset, go_namespaces)
           
    if (not in_go_desc_suppl.endswith("ALLZZZ")):  # RMS, probably dont' need to do this as the supplemental file is usually a subset of the big file
        go_descs, go_namespaces = process_desc_lines(in_go_desc_suppl, go_descs, go_subset, go_namespaces) 
     
    return go_descs, go_namespaces

def process_desc_lines(in_desc_file, go_descs, go_subset, go_namespaces):
    desc = "" 
    print("indescfile:", in_desc_file)   
    with open(in_desc_file, 'r') as in_descs:   
        lines = in_descs.readlines()
        for line in lines:            
            cols = line.strip().split(": ")
            if (cols[0] == 'id'):  #flybase specific. will need to parameterize for other orgs.
                go_term = cols[1]
                continue
            if (cols[0] == 'name'):
                desc = cols[1]
            if (cols[0] == "namespace"):
                current_go_subset = cols[1]
                if (current_go_subset  == "biological_process"):
                    go_subset_abbrev = "BP"
                elif (current_go_subset  == "molecular_function"):
                    go_subset_abbrev = "MF"
                elif (current_go_subset  == "cellular_component"):
                    go_subset_abbrev = "CC"
                else: 
                    go_subset_abbrev = current_go_subset
                if (go_subset == "ALLZZZ"):  # don't restrict by subset
                    #if (go_term == "GO:0000278"):
                    #    print(go_term)
                    #    print(desc)
                    #    import pdb
                    #    pdb.set_trace()
                    go_descs[go_term] = desc
                    go_namespaces[go_term] = go_subset_abbrev
                else:
                    if (go_subset == current_go_subset):
                        go_descs[go_term] = desc
                        go_namespaces[go_term] = go_subset_abbrev
    return go_descs, go_namespaces

def print_goseq_fixed(out_goseq_fixed, in_go_seq, go_descs, go_namespaces):
    print(out_goseq_fixed)
    with open(out_goseq_fixed, "w") as outf:
        with open(in_go_seq, 'r') as inf:
            ctr = 0
            for line in inf.readlines():
                ctr = ctr + 1
                if (ctr == 1):
                    outf.write(line)
                    continue
                parts = line.strip().split("\t")
                if (parts[5] == ""):
                    outf.write(parts[0] + "\t" + parts[1] + "\t" + parts[2] + "\t" +
                               parts[3] + "\t" + parts[4] + "\t" +
                               go_descs[parts[0]] + "\t" + 
                               go_namespaces[parts[0]] + "\t"  + parts[7] + "\n")
                else:
                    outf.write(line) 
            
def make_dir(dir):
    x = ph(dir).mkdir(exist_ok=True)  
   
if __name__ == "__main__":
    main()