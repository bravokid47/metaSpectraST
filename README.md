# metaSpectraST

![version](https://img.shields.io/badge/metaSpectraST-v0.0-green) ![license](https://img.shields.io/badge/license-MIT-green) ![commit](https://img.shields.io/github/last-commit/bravokid47/metaSpectraST?color=orange)

metaSpectraST is an unsupervised and database-independent analysis tools for metaproteomic MS/MS data using spectrum clustering. It clusters all experimentally observed MS/MS spectra based on their spectral similarity and create a representative consensus spectrum for each cluster by using the spectrum clustering algorithm implemented in the spectral library search engine, [SpectraST](http://tools.proteomecenter.org/wiki/index.php?title=Software:SpectraST). 

Spectrally similar MS/MS spectra that are grouped in one spectral cluster are presumed to originate from the sampe peptide sequence, and therefore metaSpecraST treats them as replicate spectra and quantitatively profiles samples by counting the number (spectral count, SC) or intentisity of replicate spectra (spectral index, \\[SI_N\\]) in each spectral cluster.
