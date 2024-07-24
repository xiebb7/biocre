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
Executing the BioCRE pipeline necessitates the provision of three primary inputs: `rna_adata`, `atac_adata`, and `meta_data`:
* `rna_adata` - This represents the processed single-cell RNA sequencing (sncRNA-seq) data encapsulated in an AnnData object. 
* `atac_adata` - Similarly structured as rna_adata, but this dataset comprises single-cell Assay for Transposase-Accessible Chromatin using sequencing (snATAC-seq) information. The cells in these two adata object needs to be identifcal and the pro-precessed cells are suggested as input.
* `meta_data` - This data serves as a genomic annotation resource. It includes precise genomic locations for genes and peaks. Typically, the cellranger output file named 'features.tsv.gz' serves as the metadata.

To ensure accurate integration and comparison across multi-omics data, it is imperative that the cells represented in both AnnData objects (rna_adata and atac_adata) are identical. This means that each cell in one dataset should correspond directly to a cell in the other dataset, maintaining consistency in cell identity and order.Pre-processing steps should be applied to the cells prior to inputting them into the analysis pipeline. These pre-processing measures typically include quality control checks, normalization, filtering out low-quality cells or features, and batch effect correction, if necessary. By performing these steps, you enhance the reliability and interpretability of downstream multi-omics analyses, ensuring that any observed correlations or differences are biologically meaningful rather than artifacts of technical variability.

