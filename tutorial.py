import json
from collections import defaultdict

import goseqpy as gsp

# Load example data
with open('example_data.json', 'r') as f:
    example_data = json.load(f)
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
print(results)
