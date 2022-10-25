# GOSeqPy

 ## A Python wrapper around [GOSeq](https://bioconductor.org/packages/release/bioc/html/goseq.html), a method for performing GO set enrichment on differentially expressed gene sets from RNA-seq experiments.
 
Steps To Run:
- pip install GOSeqPy
- Somehow get fish out the runGoSeqRms.py executable, either from the assets subdirectory from your installed GOSeqPackages or directly from github: https://raw.githubusercontent.com/stewart-lab/GOSeqPy/main/runGoSeqRms.py
- run it with a command line that looks something like:

```bash
python runGoSeqRms.py <my_outdir_path> <my_de_genes_file> -pvalcol=5 -pval=0.05 -log2fc=1.5 -log2fccol=2
```