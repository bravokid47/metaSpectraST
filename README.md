![logo](metaSpectraST_logo.png)

# metaSpectraST

![version](https://img.shields.io/badge/metaSpectraST-v0.0-green) ![license](https://img.shields.io/badge/license-MIT-green) ![commit](https://img.shields.io/github/last-commit/bravokid47/metaSpectraST?color=orange)

metaSpectraST is an unsupervised and database-independent analysis tools for metaproteomic MS/MS data using spectrum clustering. It clusters all experimentally observed MS/MS spectra based on their spectral similarity and create a representative consensus spectrum for each cluster by using the spectrum clustering algorithm implemented in the spectral library search engine, [SpectraST](http://tools.proteomecenter.org/wiki/index.php?title=Software:SpectraST). 

Spectrally similar MS/MS spectra that are grouped in one spectral cluster are presumed to originate from the sampe peptide sequence, and therefore metaSpecraST treats them as replicate spectra and quantitatively profiles samples by counting the number (spectral count, SC) or intentisity of replicate spectra (spectral index, SI<sub>N</sub>) in each spectral cluster.

The metaSpectraST spectral clusters also offer a portal to integrate and reconcile multiple peptide identification approacheds, including database search, open modification search, and *de novo* sequencing. For each spectral cluster, sequences of raw spectra and their cosensus spectrum assigned by different indentification methods vote for the consensus peptide sequence of the spectral cluster through a heuristic reconciliation scheme and the majority rule.

With metaSpectraST you can,

1. Fast profile and compare the microbial communities of your sample;
2. Classify your metaproteomic (or proteomic) samples;
3. Validate biological/technical replicates;
4. Integrate and reconcile multiple peptide/protein identification approaches for further taxonomic or functional studies.

# Contents
- [Installation](https://github.com/bravokid47/metaSpectraST/edit/main/README.md#installation)
  - [Dependencies](https://github.com/bravokid47/metaSpectraST/edit/main/README.md#dependencies)
  - [Installing metaSpectraST](https://github.com/bravokid47/metaSpectraST/edit/main/README.md#installing-metaspectrast)
- [Quick start]()
- [About]()

# Installation
## Dependencies
- **Python version >= 3.7, R version 4.1.3**
- **SpectraST (v5.0)**

[SpectraST](http://tools.proteomecenter.org/wiki/index.php?title=Software:SpectraST) is an integral component of the [Trans Proteomic Pipeline suite (TPP)](http://tools.proteomecenter.org/wiki/index.php?title=Software:TPP) of software. A compiled executable file is included here, which can be used alone without other TPP components.

We encourage useres to download and install the entire TPP suite, which provides other useful functionaliteissuch as raw data importation, automatic validation of search results, protein inference, and quantification and visualization. Please refer to the guides for [TPP Linux installation](http://tools.proteomecenter.org/wiki/index.php?title=Linux_Installation_Guides), and the official download site for [Windows installer](http://tools.proteomecenter.org/wiki/index.php?title=TPP:5.2_Installation).

- **edgeR (v3.34.0)**

metaSpectraST normalizes the data using the trimmed mean of M-values (TMM) normalization method implemented in the edgeR package. edgeR is not necessary if you would like to normalize the data with other methods. Please refer to [Bioconductor-edgR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html) for further information.

To install the edgeR package, start R and enter:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
```

You may also need to install the **limma** (v3.48.1) package

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
```

## Installing metaSpectraST
1. Download or clone the repository:

```shell
git clone https://github.com/bravokid47/metaSpectraST.git
```

2. Make metaSpectraST executable by adding the directory ```yourpath/metaspectrast/``` to the environment variable ```$PATH```, or just copy the following line to the ```~/.bashrc``` or ```~/.bash_profile``` file.

```shell
export PATH="$PATH:yourpath/metaspectrast";
```




