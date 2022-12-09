##############################################################################
#   Run the GOSeq R package from within Python by using r2py
##############################################################################

from optparse import OptionParser
import json
import numpy as np
import pandas as pd
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import sys
import os
from collections import defaultdict

def run_GOseq(de_genes, all_genes, gene_to_gene_sets, gene_lens, out_dir):
    out_pwf = os.path.join(out_dir, "pwfplot.pdf") 
    out_rData = os.path.join(out_dir, "environment.RData") 
    out_r_msgs = os.path.join(out_dir, "r_messages.txt") 
    de_genes = set(de_genes)
    de_status = [int(g in de_genes) for g in all_genes]
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_de_status = ro.conversion.py2rpy(de_status)
        r_all_genes = ro.conversion.py2rpy(all_genes)
        r_gene_sets = ro.vectors.ListVector(gene_to_gene_sets)
        r_gene_lens = ro.conversion.py2rpy(gene_lens)
    rstring="""
        function(de_status, all_genes, gene_sets, gene_lens, out_r_msgs, out_pwf, out_rData) {
            out_msgs_f <- file(out_r_msgs, open="wt")
            sink(out_msgs_f, type="message")
            library(goseq)
            gene_lens <- unlist(gene_lens)
            de_status <- unlist(de_status)
            all_genes <- unlist(all_genes)
            names(de_status) <- all_genes
            gene_sets <- lapply(gene_sets, as.character)
            pdf(out_pwf, width = 4, height = 4)
            pwf <- nullp(de_status, bias.data = gene_lens)
            
            dev.off()
            res <- goseq(pwf, gene2cat=gene_sets)
            sink(type="message")
            close(out_msgs_f)
            fdr <- p.adjust(res$over_represented_pvalue, method="BH")
            res$fdr <- fdr
            #res <- res[p.adjust(res$over_represented_pvalue, method="BH")<.05]
            df <- data.frame(res)
            save(list = ls(all.names = TRUE), file = out_rData)
            
            df
        }
    """
    # original lines for human analysis. Not sure why this worked for drosophila!
    #pwf <- nullp(de_status, "hg19", "geneSymbol", bias.data = gene_lens)
    #res <- goseq(pwf, "hg19","geneSymbol", gene2cat=gene_sets)
    #pwf <- nullp(de_status, "dm3", "geneSymbol", bias.data = gene_lens)
    #res <- goseq(pwf, "dm3","geneSymbol", gene2cat=gene_sets)
    # NOTE FOR POSTERITY-- We don't need the above stuff th at specifies a species (e.g. "hg19")
    # and an identifier ("geneSymbol") because we are always passing in our own genes and gene lengths.
    base = importr('base')
    
    r_func = ro.r(rstring)
    r_res = r_func(r_de_status, r_all_genes, r_gene_sets, r_gene_lens, out_r_msgs, out_pwf, out_rData)
    base.warnings()
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_res = ro.conversion.rpy2py(r_res)
    df_res = df_res.set_index('category')
    return df_res


def parse_gmt(gene_sets_f):
    gene_set_to_genes = {}
    with open(gene_sets_f, 'r') as f:
        for l in f:
            toks = l.split('\t')
            gene_set = toks[0]
            genes = [x.strip() for x in toks[2:]]
            gene_set_to_genes[gene_set] = genes
    return gene_set_to_genes


if __name__ == "__main__":
    main()
