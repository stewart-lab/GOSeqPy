# GOSeqPy

A Python wrapper around [GOSeq](https://bioconductor.org/packages/release/bioc/html/goseq.html), a method for performing GO set enrichment on differentially expressed gene sets from RNA-seq experiments.

Check requirements.txt for python requirements. Note that you'll need rpy2, which is non-standard.

On the R Side of the world, you'll need the goseq R package. 
See:  https://bioconductor.org/packages/release/bioc/html/goseq.html
for instructions on how to install goseq into your R.

-------------- Running the example command line ----------------------
Make sure you are in the src/ directory then:

For running on Ron's old Trashcan:
pwd
/Users/Ron/Desktop/GeneSetEnrichment_General/GOSeqPy/src

python runGoSeqRms.py ../goseq_out/ ../ActiveCollaborations/Tattersall/GOSeq_Tattersall/inputData_goseq/TattersallExp2//vaper.post_v_pre/group/all.Up.deSeq2.tsv  -pvalcol=5 -pval=0.05 -log2fc=1.5 -log2fccol=2

For running anywhere:

be in the src directory then:
python runGoSeqRms.py ../goseq_out/ ../data/all.Up.deSeq2.tsv  -pvalcol=5 -pval=0.05 -log2fc=1.5 -log2fccol=2
