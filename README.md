<img width="263" alt="image" src="https://github.com/bioinformatic-seragnoli/BOBaFIT/assets/52926104/66a79d80-690e-4955-9f8b-9a5603e4efbd">

# **BOBaFIT: Refitting diploid region profiles using a clustering procedure**

This is the repository of the code used in R package "BOBaFIT", published in Bioconductor Release (3.18). 

See BIOCONDUCTOR page for installation, documentation and details: https://www.bioconductor.org/packages/release/bioc/html/BOBaFIT.html 



This package provides a method to refit and correct the diploid region in copy number profiles. It uses a clustering algorithm to identify pathology-specific normal (diploid) chromosomes and then use their copy number signal to refit the whole profile. 
The package is composed by three functions: DRrefit (the main function), ComputeNormalChromosome and PlotCluster.

## **Please cite:**
_Mazzocchetti, G., Poletti, A., Solli, V., Borsi, E., Martello, M., Vigliotta, I., ... & Terragna, C. (2022). 
BoBafit: A copy number clustering tool designed to refit and recalibrate the baseline region of tumors’ profiles. Computational and Structural Biotechnology Journal, 20, 3718-3728._ DOI: https://doi.org/10.1016/j.csbj.2022.06.062
