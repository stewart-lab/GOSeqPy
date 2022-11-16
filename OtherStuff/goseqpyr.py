# -*- coding: utf-8 -*-
"""GoSeqPYR.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/11UAx2nNftWhkRt-UU4Bmdm5A14oR3wKd
"""

from google.colab import drive
drive.mount('/content/drive')

!ln -s  '/content/drive/My Drive/' /content/drive/MyDrive/

!ln -s "/usr/local/lib/R/site-library" "/content/drive/MyDrive/GoSeqPy/usr/local/lib/R/site-library"

!ls "/content/drive/MyDrive/"

# Commented out IPython magic to ensure Python compatibility.
# %cd "/content/drive/MyDrive/GoSeqPy"

# Commented out IPython magic to ensure Python compatibility.
# %env PYTHONPATH=/content/drive/MyDrive/cmdlogtime

# Commented out IPython magic to ensure Python compatibility.
# %load_ext rpy2.ipython

# Commented out IPython magic to ensure Python compatibility.
# %%R
# .libPaths('usr/local/lib/R/site-library')

# Commented out IPython magic to ensure Python compatibility.
# %%R
# .libPaths()

# Commented out IPython magic to ensure Python compatibility.
# %%R
# #.libPaths("usr/local/lib/R/site-library/")  # this doesn't seem to work. The package gets installed, but rpy2 cannot find it.
# .libPaths("/usr/local/lib/R/site-library/")  # this might work?  Yes. But slow to load. Takes 15-20 minutes to do the process.
# #.libPaths("library")  # this works. Current theory is rpy2 knows how to look in "library", but not other R places. probably could configure this..
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("goseq")

# Commented out IPython magic to ensure Python compatibility.
# %%R
# library("goseq")

!python /content/drive/MyDrive/GoSeqPy/runGoSeqRms.py ./goseq_out/ ./all.Up.deSeq2.tsv  -pvalcol=5 -pval=0.05 -log2fc=1.5 -log2fccol=2

"""ONE TIME STUFF BELOW HERE -----------------------------------------------------"""

# Commented out IPython magic to ensure Python compatibility.
# %%R
# .libPaths("library")

!echo $PYTHONPATH

# Commented out IPython magic to ensure Python compatibility.
# %ls "/usr/local/lib/R/site-library"

!ls usr/local/lib/R/site-library/

# Commented out IPython magic to ensure Python compatibility.
# %%R
# tar("library.tar.gz", "/usr/local/lib/R/site-library")
#

# Commented out IPython magic to ensure Python compatibility.
# %%R
# #NOTE-- THIS DOESN"T WORK
# install.packages("goseq")
# tar("library.tar.gz", "/usr/local/lib/R/site-library")

cp library.tar.gz /content/drive/MyDrive/src/

# Commented out IPython magic to ensure Python compatibility.
# %%R
# remove.packages("goseq")

# Commented out IPython magic to ensure Python compatibility.
# %%R
# .libPaths("library")
#

!tar xf /content/drive/MyDrive/src/library.tar.gz

# Commented out IPython magic to ensure Python compatibility.
# %%R
# .libPaths()