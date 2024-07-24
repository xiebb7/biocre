## BioCRE
<img src="figure/BioCRE_readme.png" width="600" align='center'>
BioCRE, an innovative computational framework, employs a sophisticated bi-orientation regression model to analyze multi-omics datasets. This approach deciphers the complex interactions within Gene Regulatory Networks (GRNs) at the chromosomal scale, pinpointing crucial Cis-Regulatory Elements (CREs) that play pivotal roles in gene expression and regulation. By integrating diverse genomic information, BioCRE enhances our understanding of cellular processes and disease mechanisms, offering new avenues for therapeutic intervention and personalized medicine strategies.

## Installation
Installation with virtual environment are suggested:
```
conda create -n biocre python=3.8
```
To install BioCRE, make sure you have [PyTorch](https://pytorch.org/) installed.
```
pip install torch torchvision torchaudio -i https://mirrors.cloud.tencent.com/pypi/simple
```
Then install BioCRE by pip:
```
pip install biocre
```

## Usage of BioCRE
Executing the BioCRE computational framework necessitates the provision of three primary inputs: `rna_adata`, `atac_adata`, and `meta_data`. Here's a refined explanation of these components:
`rna_adata` - This represents the processed single-cell RNA sequencing (sncRNA-seq) data encapsulated in an AnnData object. 
`atac_adata` - Similarly structured as rna_adata, but this dataset comprises single-cell Assay for Transposase-Accessible Chromatin using sequencing (snATAC-seq) information. 
`meta_data` - This data serves as a genomic annotation resource. It includes precise genomic locations for genes and peaks. Typically, the cellranger output file named 'features.tsv.gz' serves as the metadata.


