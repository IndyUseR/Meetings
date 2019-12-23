
<style>
.small-code pre code {
  font-size: 1em;
}
</style>
Analysis Of RNA-seq Data Using R
========================================================
author: Nadia Atallah Lanman
date: November 19, 2019
autosize: true
css: custom.css

Central Dogma of Molecular Biology
=========================================================
<div class="midcenter" style="margin-left:-500px; margin-top:-100px;">
<img src="./AnalysisOfGeneExpressionDataUsingR-figure/centralDogma.png"></img>
</div>

- By measuring how many mRNA molecules there are for each gene, we can quantify which genes are "on" or "off" 

Why do gene expression analyses?
========================================================

- Now feasible and simple to measure the expression levels of thousands of genes simultaneously
- Gene Annotation
- Functional studies
- Some molecular features are only seen at RNA level

Several types of gene expression analyses
========================================================

- qRT-PCR
- microarray
- RNA-seq
  - bulk
  - single-cell

RNA-Seq
========================================================

- High-throughput sequencing of RNA
- Allows for quantification of gene expression and differential expression analyses
- Characterization of alternative splicing
- de novo transcriptome assembly
  - no genome sequence necessary!

Basic RNA-seq Protocol
========================================================
<div class="midcenter" style="margin-left:-300px; margin-top:-300px;">
<img src="./AnalysisOfGeneExpressionDataUsingR-figure/rnaseqProtocol.png"></img>
</div>

Terminology
========================================================
- Reads = the sequenced portion of cDNA fragments
- Depth= (read length)(number of reads) / (haploid genome length)
- Single-end= cDNA fragments are sequenced from only one end (1x100)
- Paired-end= cDNA fragments are sequenced from both ends (2x100)
- Strand-specific= you know whether the read originated from the + or - strand

FASTQ file format
=======================================================
- Raw sequence data generally is in FASTQ format
  - Text files with both sequence and quality information
  - Integer scores converted to ASCII characters

$$ Q = -10log_{10}P $$
<div class="midcenter" style="margin-left:-375px; margin-top:100px;">
<img src="./AnalysisOfGeneExpressionDataUsingR-figure/fastqRead.png"></img>
</div>

Basic RNA-seq analysis workflow
========================================================
<div class="midcenter" style="margin-left:-350px; margin-top:-300px;">
<img src="./AnalysisOfGeneExpressionDataUsingR-figure/workflow.png"></img>
</div>

Install TCGAbiolinks
=========================================================
class: small-code

```r
#install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
#citation("TCGAbiolinks")
```
- TCGAbiolinks is a package for accessing and doing basic analyses on TCGA (The Cancer Genome Atlas) data.
- TCGA data surveys over 33 different tumor types
  - Allows comparisons to be made between cancer and normal tumors


Install SummarizedExperiment
==========================================================
class: small-code

```r
#BiocManager::install("SummarizedExperiment")

library(SummarizedExperiment)
```
- The SummarizedExperiment Class is a container to store matrices of experimental results and is often used to store gene expression data (as well as other sequencing-related datasets)
- This container is very helpful in keeping experimental data and meta data in sync

Download TCGA Count Data
========================================================
class: small-code
- Use the Genomic Data Commons to search for TCGA data for prostate cancer

```r
# Gene expression aligned against hg38
query.counts <- GDCquery(project = "TCGA-PRAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type="HTSeq - Counts")
GDCdownload(query.counts)
#read the downloaded data and prepare in into an R object
data.counts <- GDCprepare(query.counts)
```

```
|                                                    |  0%                      |                                                    |0.1814882% ~2 m remaining |                                                    |0.3629764% ~2 m remaining |                                                    |0.5444646% ~2 m remaining |                                                    |0.7259528% ~1 m remaining |                                                    |0.907441% ~3 m remaining  |                                                    |1.088929% ~3 m remaining  |                                                    |1.270417% ~3 m remaining  |                                                    |1.451906% ~2 m remaining  |                                                    |1.633394% ~2 m remaining  |                                                    |1.814882% ~2 m remaining  |=                                                   |1.99637% ~2 m remaining   |=                                                   |2.177858% ~2 m remaining  |=                                                   |2.359347% ~2 m remaining  |=                                                   |2.540835% ~2 m remaining  |=                                                   |2.722323% ~2 m remaining  |=                                                   |2.903811% ~2 m remaining  |=                                                   |3.085299% ~2 m remaining  |=                                                   |3.266788% ~2 m remaining  |=                                                   |3.448276% ~2 m remaining  |=                                                   |3.629764% ~2 m remaining  |=                                                   |3.811252% ~2 m remaining  |==                                                  |3.99274% ~1 m remaining   |==                                                  |4.174229% ~1 m remaining  |==                                                  |4.355717% ~1 m remaining  |==                                                  |4.537205% ~1 m remaining  |==                                                  |4.718693% ~1 m remaining  |==                                                  |4.900181% ~1 m remaining  |==                                                  |5.08167% ~1 m remaining   |==                                                  |5.263158% ~1 m remaining  |==                                                  |5.444646% ~1 m remaining  |==                                                  |5.626134% ~1 m remaining  |===                                                 |5.807623% ~1 m remaining  |===                                                 |5.989111% ~1 m remaining  |===                                                 |6.170599% ~1 m remaining  |===                                                 |6.352087% ~1 m remaining  |===                                                 |6.533575% ~1 m remaining  |===                                                 |6.715064% ~1 m remaining  |===                                                 |6.896552% ~1 m remaining  |===                                                 |7.07804% ~1 m remaining   |===                                                 |7.259528% ~1 m remaining  |===                                                 |7.441016% ~1 m remaining  |===                                                 |7.622505% ~1 m remaining  |====                                                |7.803993% ~1 m remaining  |====                                                |7.985481% ~1 m remaining  |====                                                |8.166969% ~1 m remaining  |====                                                |8.348457% ~1 m remaining  |====                                                |8.529946% ~1 m remaining  |====                                                |8.711434% ~1 m remaining  |====                                                |8.892922% ~1 m remaining  |====                                                |9.07441% ~1 m remaining   |====                                                |9.255898% ~1 m remaining  |====                                                |9.437387% ~1 m remaining  |=====                                               |9.618875% ~1 m remaining  |=====                                               |9.800363% ~1 m remaining  |=====                                               |9.981851% ~1 m remaining  |=====                                               |10.16334% ~1 m remaining  |=====                                               |10.34483% ~1 m remaining  |=====                                               |10.52632% ~1 m remaining  |=====                                               |10.7078% ~1 m remaining   |=====                                               |10.88929% ~1 m remaining  |=====                                               |11.07078% ~1 m remaining  |=====                                               |11.25227% ~1 m remaining  |=====                                               |11.43376% ~1 m remaining  |======                                              |11.61525% ~1 m remaining  |======                                              |11.79673% ~1 m remaining  |======                                              |11.97822% ~1 m remaining  |======                                              |12.15971% ~1 m remaining  |======                                              |12.3412% ~1 m remaining   |======                                              |12.52269% ~1 m remaining  |======                                              |12.70417% ~1 m remaining  |======                                              |12.88566% ~1 m remaining  |======                                              |13.06715% ~1 m remaining  |======                                              |13.24864% ~1 m remaining  |======                                              |13.43013% ~1 m remaining  |=======                                             |13.61162% ~1 m remaining  |=======                                             |13.7931% ~1 m remaining   |=======                                             |13.97459% ~1 m remaining  |=======                                             |14.15608% ~1 m remaining  |=======                                             |14.33757% ~1 m remaining  |=======                                             |14.51906% ~1 m remaining  |=======                                             |14.70054% ~1 m remaining  |=======                                             |14.88203% ~1 m remaining  |=======                                             |15.06352% ~1 m remaining  |=======                                             |15.24501% ~1 m remaining  |========                                            |15.4265% ~1 m remaining   |========                                            |15.60799% ~1 m remaining  |========                                            |15.78947% ~1 m remaining  |========                                            |15.97096% ~1 m remaining  |========                                            |16.15245% ~1 m remaining  |========                                            |16.33394% ~1 m remaining  |========                                            |16.51543% ~60 s remaining |========                                            |16.69691% ~59 s remaining |========                                            |16.8784% ~59 s remaining  |========                                            |17.05989% ~59 s remaining |========                                            |17.24138% ~58 s remaining |=========                                           |17.42287% ~58 s remaining |=========                                           |17.60436% ~57 s remaining |=========                                           |17.78584% ~57 s remaining |=========                                           |17.96733% ~57 s remaining |=========                                           |18.14882% ~57 s remaining |=========                                           |18.33031% ~56 s remaining |=========                                           |18.5118% ~56 s remaining  |=========                                           |18.69328% ~56 s remaining |=========                                           |18.87477% ~55 s remaining |=========                                           |19.05626% ~55 s remaining |==========                                          |19.23775% ~55 s remaining |==========                                          |19.41924% ~55 s remaining |==========                                          |19.60073% ~1 m remaining  |==========                                          |19.78221% ~1 m remaining  |==========                                          |19.9637% ~60 s remaining  |==========                                          |20.14519% ~1 m remaining  |==========                                          |20.32668% ~1 m remaining  |==========                                          |20.50817% ~1 m remaining  |==========                                          |20.68966% ~60 s remaining |==========                                          |20.87114% ~59 s remaining |==========                                          |21.05263% ~59 s remaining |===========                                         |21.23412% ~1 m remaining  |===========                                         |21.41561% ~1 m remaining  |===========                                         |21.5971% ~1 m remaining   |===========                                         |21.77858% ~1 m remaining  |===========                                         |21.96007% ~1 m remaining  |===========                                         |22.14156% ~1 m remaining  |===========                                         |22.32305% ~1 m remaining  |===========                                         |22.50454% ~1 m remaining  |===========                                         |22.68603% ~1 m remaining  |===========                                         |22.86751% ~1 m remaining  |===========                                         |23.049% ~1 m remaining    |============                                        |23.23049% ~1 m remaining  |============                                        |23.41198% ~1 m remaining  |============                                        |23.59347% ~1 m remaining  |============                                        |23.77495% ~1 m remaining  |============                                        |23.95644% ~1 m remaining  |============                                        |24.13793% ~1 m remaining  |============                                        |24.31942% ~1 m remaining  |============                                        |24.50091% ~1 m remaining  |============                                        |24.6824% ~1 m remaining   |============                                        |24.86388% ~1 m remaining  |=============                                       |25.04537% ~1 m remaining  |=============                                       |25.22686% ~60 s remaining |=============                                       |25.40835% ~59 s remaining |=============                                       |25.58984% ~59 s remaining |=============                                       |25.77132% ~59 s remaining |=============                                       |25.95281% ~1 m remaining  |=============                                       |26.1343% ~1 m remaining   |=============                                       |26.31579% ~1 m remaining  |=============                                       |26.49728% ~1 m remaining  |=============                                       |26.67877% ~60 s remaining |=============                                       |26.86025% ~59 s remaining |==============                                      |27.04174% ~59 s remaining |==============                                      |27.22323% ~59 s remaining |==============                                      |27.40472% ~59 s remaining |==============                                      |27.58621% ~58 s remaining |==============                                      |27.7677% ~58 s remaining  |==============                                      |27.94918% ~58 s remaining |==============                                      |28.13067% ~57 s remaining |==============                                      |28.31216% ~57 s remaining |==============                                      |28.49365% ~57 s remaining |==============                                      |28.67514% ~57 s remaining |===============                                     |28.85662% ~56 s remaining |===============                                     |29.03811% ~56 s remaining |===============                                     |29.2196% ~56 s remaining  |===============                                     |29.40109% ~55 s remaining |===============                                     |29.58258% ~55 s remaining |===============                                     |29.76407% ~55 s remaining |===============                                     |29.94555% ~54 s remaining |===============                                     |30.12704% ~54 s remaining |===============                                     |30.30853% ~54 s remaining |===============                                     |30.49002% ~54 s remaining |===============                                     |30.67151% ~53 s remaining |================                                    |30.85299% ~53 s remaining |================                                    |31.03448% ~53 s remaining |================                                    |31.21597% ~53 s remaining |================                                    |31.39746% ~53 s remaining |================                                    |31.57895% ~52 s remaining |================                                    |31.76044% ~52 s remaining |================                                    |31.94192% ~52 s remaining |================                                    |32.12341% ~52 s remaining |================                                    |32.3049% ~52 s remaining  |================                                    |32.48639% ~52 s remaining |================                                    |32.66788% ~52 s remaining |=================                                   |32.84936% ~51 s remaining |=================                                   |33.03085% ~51 s remaining |=================                                   |33.21234% ~51 s remaining |=================                                   |33.39383% ~51 s remaining |=================                                   |33.57532% ~50 s remaining |=================                                   |33.75681% ~50 s remaining |=================                                   |33.93829% ~50 s remaining |=================                                   |34.11978% ~50 s remaining |=================                                   |34.30127% ~50 s remaining |=================                                   |34.48276% ~50 s remaining |==================                                  |34.66425% ~49 s remaining |==================                                  |34.84574% ~49 s remaining |==================                                  |35.02722% ~49 s remaining |==================                                  |35.20871% ~49 s remaining |==================                                  |35.3902% ~49 s remaining  |==================                                  |35.57169% ~49 s remaining |==================                                  |35.75318% ~48 s remaining |==================                                  |35.93466% ~48 s remaining |==================                                  |36.11615% ~48 s remaining |==================                                  |36.29764% ~48 s remaining |==================                                  |36.47913% ~50 s remaining |===================                                 |36.66062% ~50 s remaining |===================                                 |36.84211% ~50 s remaining |===================                                 |37.02359% ~49 s remaining |===================                                 |37.20508% ~49 s remaining |===================                                 |37.38657% ~49 s remaining |===================                                 |37.56806% ~49 s remaining |===================                                 |37.74955% ~48 s remaining |===================                                 |37.93103% ~48 s remaining |===================                                 |38.11252% ~48 s remaining |===================                                 |38.29401% ~48 s remaining |====================                                |38.4755% ~48 s remaining  |====================                                |38.65699% ~47 s remaining |====================                                |38.83848% ~47 s remaining |====================                                |39.01996% ~47 s remaining |====================                                |39.20145% ~47 s remaining |====================                                |39.38294% ~46 s remaining |====================                                |39.56443% ~46 s remaining |====================                                |39.74592% ~46 s remaining |====================                                |39.9274% ~46 s remaining  |====================                                |40.10889% ~46 s remaining |====================                                |40.29038% ~45 s remaining |=====================                               |40.47187% ~45 s remaining |=====================                               |40.65336% ~45 s remaining |=====================                               |40.83485% ~45 s remaining |=====================                               |41.01633% ~45 s remaining |=====================                               |41.19782% ~44 s remaining |=====================                               |41.37931% ~44 s remaining |=====================                               |41.5608% ~44 s remaining  |=====================                               |41.74229% ~44 s remaining |=====================                               |41.92377% ~44 s remaining |=====================                               |42.10526% ~43 s remaining |=====================                               |42.28675% ~43 s remaining |======================                              |42.46824% ~43 s remaining |======================                              |42.64973% ~43 s remaining |======================                              |42.83122% ~43 s remaining |======================                              |43.0127% ~42 s remaining  |======================                              |43.19419% ~42 s remaining |======================                              |43.37568% ~42 s remaining |======================                              |43.55717% ~42 s remaining |======================                              |43.73866% ~42 s remaining |======================                              |43.92015% ~41 s remaining |======================                              |44.10163% ~41 s remaining |=======================                             |44.28312% ~41 s remaining |=======================                             |44.46461% ~41 s remaining |=======================                             |44.6461% ~41 s remaining  |=======================                             |44.82759% ~41 s remaining |=======================                             |45.00907% ~40 s remaining |=======================                             |45.19056% ~40 s remaining |=======================                             |45.37205% ~40 s remaining |=======================                             |45.55354% ~40 s remaining |=======================                             |45.73503% ~40 s remaining |=======================                             |45.91652% ~39 s remaining |=======================                             |46.098% ~39 s remaining   |========================                            |46.27949% ~39 s remaining |========================                            |46.46098% ~39 s remaining |========================                            |46.64247% ~39 s remaining |========================                            |46.82396% ~39 s remaining |========================                            |47.00544% ~38 s remaining |========================                            |47.18693% ~38 s remaining |========================                            |47.36842% ~38 s remaining |========================                            |47.54991% ~38 s remaining |========================                            |47.7314% ~38 s remaining  |========================                            |47.91289% ~38 s remaining |=========================                           |48.09437% ~37 s remaining |=========================                           |48.27586% ~38 s remaining |=========================                           |48.45735% ~38 s remaining |=========================                           |48.63884% ~38 s remaining |=========================                           |48.82033% ~37 s remaining |=========================                           |49.00181% ~37 s remaining |=========================                           |49.1833% ~37 s remaining  |=========================                           |49.36479% ~37 s remaining |=========================                           |49.54628% ~37 s remaining |=========================                           |49.72777% ~37 s remaining |=========================                           |49.90926% ~36 s remaining |==========================                          |50.09074% ~36 s remaining |==========================                          |50.27223% ~36 s remaining |==========================                          |50.45372% ~36 s remaining |==========================                          |50.63521% ~36 s remaining |==========================                          |50.8167% ~35 s remaining  |==========================                          |50.99819% ~35 s remaining |==========================                          |51.17967% ~35 s remaining |==========================                          |51.36116% ~35 s remaining |==========================                          |51.54265% ~35 s remaining |==========================                          |51.72414% ~35 s remaining |==========================                          |51.90563% ~34 s remaining |===========================                         |52.08711% ~34 s remaining |===========================                         |52.2686% ~34 s remaining  |===========================                         |52.45009% ~34 s remaining |===========================                         |52.63158% ~34 s remaining |===========================                         |52.81307% ~34 s remaining |===========================                         |52.99456% ~33 s remaining |===========================                         |53.17604% ~33 s remaining |===========================                         |53.35753% ~33 s remaining |===========================                         |53.53902% ~33 s remaining |===========================                         |53.72051% ~33 s remaining |============================                        |53.902% ~32 s remaining   |============================                        |54.08348% ~32 s remaining |============================                        |54.26497% ~32 s remaining |============================                        |54.44646% ~32 s remaining |============================                        |54.62795% ~32 s remaining |============================                        |54.80944% ~32 s remaining |============================                        |54.99093% ~31 s remaining |============================                        |55.17241% ~31 s remaining |============================                        |55.3539% ~31 s remaining  |============================                        |55.53539% ~31 s remaining |============================                        |55.71688% ~31 s remaining |=============================                       |55.89837% ~31 s remaining |=============================                       |56.07985% ~30 s remaining |=============================                       |56.26134% ~30 s remaining |=============================                       |56.44283% ~30 s remaining |=============================                       |56.62432% ~30 s remaining |=============================                       |56.80581% ~30 s remaining |=============================                       |56.9873% ~30 s remaining  |=============================                       |57.16878% ~29 s remaining |=============================                       |57.35027% ~29 s remaining |=============================                       |57.53176% ~29 s remaining |==============================                      |57.71325% ~29 s remaining |==============================                      |57.89474% ~29 s remaining |==============================                      |58.07623% ~29 s remaining |==============================                      |58.25771% ~29 s remaining |==============================                      |58.4392% ~28 s remaining  |==============================                      |58.62069% ~28 s remaining |==============================                      |58.80218% ~28 s remaining |==============================                      |58.98367% ~28 s remaining |==============================                      |59.16515% ~28 s remaining |==============================                      |59.34664% ~28 s remaining |==============================                      |59.52813% ~27 s remaining |===============================                     |59.70962% ~27 s remaining |===============================                     |59.89111% ~27 s remaining |===============================                     |60.0726% ~27 s remaining  |===============================                     |60.25408% ~27 s remaining |===============================                     |60.43557% ~27 s remaining |===============================                     |60.61706% ~27 s remaining |===============================                     |60.79855% ~26 s remaining |===============================                     |60.98004% ~26 s remaining |===============================                     |61.16152% ~26 s remaining |===============================                     |61.34301% ~26 s remaining |===============================                     |61.5245% ~26 s remaining  |================================                    |61.70599% ~26 s remaining |================================                    |61.88748% ~25 s remaining |================================                    |62.06897% ~26 s remaining |================================                    |62.25045% ~26 s remaining |================================                    |62.43194% ~26 s remaining |================================                    |62.61343% ~25 s remaining |================================                    |62.79492% ~25 s remaining |================================                    |62.97641% ~25 s remaining |================================                    |63.15789% ~25 s remaining |================================                    |63.33938% ~25 s remaining |=================================                   |63.52087% ~25 s remaining |=================================                   |63.70236% ~25 s remaining |=================================                   |63.88385% ~24 s remaining |=================================                   |64.06534% ~24 s remaining |=================================                   |64.24682% ~24 s remaining |=================================                   |64.42831% ~24 s remaining |=================================                   |64.6098% ~24 s remaining  |=================================                   |64.79129% ~24 s remaining |=================================                   |64.97278% ~23 s remaining |=================================                   |65.15426% ~23 s remaining |=================================                   |65.33575% ~23 s remaining |==================================                  |65.51724% ~23 s remaining |==================================                  |65.69873% ~23 s remaining |==================================                  |65.88022% ~23 s remaining |==================================                  |66.06171% ~23 s remaining |==================================                  |66.24319% ~22 s remaining |==================================                  |66.42468% ~22 s remaining |==================================                  |66.60617% ~22 s remaining |==================================                  |66.78766% ~22 s remaining |==================================                  |66.96915% ~22 s remaining |==================================                  |67.15064% ~22 s remaining |===================================                 |67.33212% ~22 s remaining |===================================                 |67.51361% ~21 s remaining |===================================                 |67.6951% ~21 s remaining  |===================================                 |67.87659% ~21 s remaining |===================================                 |68.05808% ~21 s remaining |===================================                 |68.23956% ~21 s remaining |===================================                 |68.42105% ~21 s remaining |===================================                 |68.60254% ~21 s remaining |===================================                 |68.78403% ~20 s remaining |===================================                 |68.96552% ~20 s remaining |===================================                 |69.14701% ~20 s remaining |====================================                |69.32849% ~20 s remaining |====================================                |69.50998% ~20 s remaining |====================================                |69.69147% ~20 s remaining |====================================                |69.87296% ~20 s remaining |====================================                |70.05445% ~19 s remaining |====================================                |70.23593% ~19 s remaining |====================================                |70.41742% ~19 s remaining |====================================                |70.59891% ~19 s remaining |====================================                |70.7804% ~19 s remaining  |====================================                |70.96189% ~19 s remaining |====================================                |71.14338% ~19 s remaining |=====================================               |71.32486% ~19 s remaining |=====================================               |71.50635% ~18 s remaining |=====================================               |71.68784% ~18 s remaining |=====================================               |71.86933% ~18 s remaining |=====================================               |72.05082% ~18 s remaining |=====================================               |72.2323% ~18 s remaining  |=====================================               |72.41379% ~18 s remaining |=====================================               |72.59528% ~18 s remaining |=====================================               |72.77677% ~17 s remaining |=====================================               |72.95826% ~17 s remaining |======================================              |73.13975% ~17 s remaining |======================================              |73.32123% ~17 s remaining |======================================              |73.50272% ~17 s remaining |======================================              |73.68421% ~17 s remaining |======================================              |73.8657% ~17 s remaining  |======================================              |74.04719% ~17 s remaining |======================================              |74.22868% ~16 s remaining |======================================              |74.41016% ~16 s remaining |======================================              |74.59165% ~16 s remaining |======================================              |74.77314% ~16 s remaining |======================================              |74.95463% ~16 s remaining |=======================================             |75.13612% ~16 s remaining |=======================================             |75.3176% ~16 s remaining  |=======================================             |75.49909% ~16 s remaining |=======================================             |75.68058% ~15 s remaining |=======================================             |75.86207% ~15 s remaining |=======================================             |76.04356% ~15 s remaining |=======================================             |76.22505% ~15 s remaining |=======================================             |76.40653% ~15 s remaining |=======================================             |76.58802% ~15 s remaining |=======================================             |76.76951% ~15 s remaining |========================================            |76.951% ~14 s remaining   |========================================            |77.13249% ~14 s remaining |========================================            |77.31397% ~14 s remaining |========================================            |77.49546% ~14 s remaining |========================================            |77.67695% ~14 s remaining |========================================            |77.85844% ~14 s remaining |========================================            |78.03993% ~14 s remaining |========================================            |78.22142% ~14 s remaining |========================================            |78.4029% ~14 s remaining  |========================================            |78.58439% ~14 s remaining |========================================            |78.76588% ~13 s remaining |=========================================           |78.94737% ~13 s remaining |=========================================           |79.12886% ~13 s remaining |=========================================           |79.31034% ~13 s remaining |=========================================           |79.49183% ~13 s remaining |=========================================           |79.67332% ~13 s remaining |=========================================           |79.85481% ~13 s remaining |=========================================           |80.0363% ~13 s remaining  |=========================================           |80.21779% ~12 s remaining |=========================================           |80.39927% ~12 s remaining |=========================================           |80.58076% ~12 s remaining |=========================================           |80.76225% ~12 s remaining |==========================================          |80.94374% ~12 s remaining |==========================================          |81.12523% ~12 s remaining |==========================================          |81.30672% ~12 s remaining |==========================================          |81.4882% ~12 s remaining  |==========================================          |81.66969% ~11 s remaining |==========================================          |81.85118% ~11 s remaining |==========================================          |82.03267% ~11 s remaining |==========================================          |82.21416% ~11 s remaining |==========================================          |82.39564% ~11 s remaining |==========================================          |82.57713% ~11 s remaining |===========================================         |82.75862% ~11 s remaining |===========================================         |82.94011% ~11 s remaining |===========================================         |83.1216% ~10 s remaining  |===========================================         |83.30309% ~10 s remaining |===========================================         |83.48457% ~10 s remaining |===========================================         |83.66606% ~10 s remaining |===========================================         |83.84755% ~10 s remaining |===========================================         |84.02904% ~10 s remaining |===========================================         |84.21053% ~10 s remaining |===========================================         |84.39201% ~10 s remaining |===========================================         |84.5735% ~10 s remaining  |============================================        |84.75499% ~9 s remaining  |============================================        |84.93648% ~9 s remaining  |============================================        |85.11797% ~9 s remaining  |============================================        |85.29946% ~9 s remaining  |============================================        |85.48094% ~9 s remaining  |============================================        |85.66243% ~9 s remaining  |============================================        |85.84392% ~9 s remaining  |============================================        |86.02541% ~9 s remaining  |============================================        |86.2069% ~8 s remaining   |============================================        |86.38838% ~8 s remaining  |=============================================       |86.56987% ~8 s remaining  |=============================================       |86.75136% ~8 s remaining  |=============================================       |86.93285% ~8 s remaining  |=============================================       |87.11434% ~8 s remaining  |=============================================       |87.29583% ~8 s remaining  |=============================================       |87.47731% ~8 s remaining  |=============================================       |87.6588% ~8 s remaining   |=============================================       |87.84029% ~7 s remaining  |=============================================       |88.02178% ~7 s remaining  |=============================================       |88.20327% ~7 s remaining  |=============================================       |88.38475% ~7 s remaining  |==============================================      |88.56624% ~7 s remaining  |==============================================      |88.74773% ~7 s remaining  |==============================================      |88.92922% ~7 s remaining  |==============================================      |89.11071% ~7 s remaining  |==============================================      |89.2922% ~7 s remaining   |==============================================      |89.47368% ~6 s remaining  |==============================================      |89.65517% ~6 s remaining  |==============================================      |89.83666% ~6 s remaining  |==============================================      |90.01815% ~6 s remaining  |==============================================      |90.19964% ~6 s remaining  |==============================================      |90.38113% ~6 s remaining  |===============================================     |90.56261% ~6 s remaining  |===============================================     |90.7441% ~6 s remaining   |===============================================     |90.92559% ~5 s remaining  |===============================================     |91.10708% ~5 s remaining  |===============================================     |91.28857% ~5 s remaining  |===============================================     |91.47005% ~5 s remaining  |===============================================     |91.65154% ~5 s remaining  |===============================================     |91.83303% ~5 s remaining  |===============================================     |92.01452% ~5 s remaining  |===============================================     |92.19601% ~5 s remaining  |================================================    |92.3775% ~5 s remaining   |================================================    |92.55898% ~4 s remaining  |================================================    |92.74047% ~4 s remaining  |================================================    |92.92196% ~4 s remaining  |================================================    |93.10345% ~4 s remaining  |================================================    |93.28494% ~4 s remaining  |================================================    |93.46642% ~4 s remaining  |================================================    |93.64791% ~4 s remaining  |================================================    |93.8294% ~4 s remaining   |================================================    |94.01089% ~4 s remaining  |================================================    |94.19238% ~3 s remaining  |=================================================   |94.37387% ~3 s remaining  |=================================================   |94.55535% ~3 s remaining  |=================================================   |94.73684% ~3 s remaining  |=================================================   |94.91833% ~3 s remaining  |=================================================   |95.09982% ~3 s remaining  |=================================================   |95.28131% ~3 s remaining  |=================================================   |95.46279% ~3 s remaining  |=================================================   |95.64428% ~3 s remaining  |=================================================   |95.82577% ~2 s remaining  |=================================================   |96.00726% ~2 s remaining  |==================================================  |96.18875% ~2 s remaining  |==================================================  |96.37024% ~2 s remaining  |==================================================  |96.55172% ~2 s remaining  |==================================================  |96.73321% ~2 s remaining  |==================================================  |96.9147% ~2 s remaining   |==================================================  |97.09619% ~2 s remaining  |==================================================  |97.27768% ~2 s remaining  |==================================================  |97.45917% ~2 s remaining  |==================================================  |97.64065% ~1 s remaining  |==================================================  |97.82214% ~1 s remaining  |==================================================  |98.00363% ~1 s remaining  |=================================================== |98.18512% ~1 s remaining  |=================================================== |98.36661% ~1 s remaining  |=================================================== |98.54809% ~1 s remaining  |=================================================== |98.72958% ~1 s remaining  |=================================================== |98.91107% ~1 s remaining  |=================================================== |99.09256% ~1 s remaining  |=================================================== |99.27405% ~0 s remaining  |=================================================== |99.45554% ~0 s remaining  |=================================================== |99.63702% ~0 s remaining  |=================================================== |99.81851% ~0 s remaining  |====================================================|100% ~0 s remaining       |====================================================|100%                      Completed after 60 s 
```

Download TCGA Clinical Data
========================================================
class: small-code

```r
# Biospecimen Data
query.biospecimen <- GDCquery(project = "TCGA-PRAD",
                           data.category = "Biospecimen",
                           file.type = "xml")
GDCdownload(query.biospecimen)
prad.biospec <- GDCprepare_clinic(query.biospecimen,clinical.info = "sample")
```

```
  |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%
```

```r
#write.csv(prad.biospec, "Biospecimen_prad.csv")
```

Prepare a count matrix
========================================================
class: small-code

```r
counts <- assay(data.counts, "HTSeq - Counts")
dim(counts)
```

```
[1] 56512   551
```

```r
#write.csv(counts, file = "all.prad.htseq.counts.csv")
coldata<-(colData(data.counts))
```
- Genes are in rows and samples are in columns
- Data is raw counts generated using a program called [HTSeq](https://htseq.readthedocs.io/en/release_0.11.1/)

Sample Barcodes
========================================================
class: small-code

```r
colnames(counts)[1:5]
```

```
[1] "TCGA-YL-A8SI-01A-11R-A41O-07" "TCGA-V1-A8WW-01A-11R-A37L-07"
[3] "TCGA-KK-A6E5-01A-11R-A311-07" "TCGA-HC-A6HY-01A-11R-A31N-07"
[5] "TCGA-KK-A8IG-01A-11R-A36G-07"
```
- Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29.  See [codes](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes) for more details
<div class="midcenter" style="margin-left:-375px; margin-top:50px;">
<img src="./AnalysisOfGeneExpressionDataUsingR-figure/barcode.png"></img>
</div>

Separate the Samples
===========================================================
class: small-code
- Keep tumor samples only

```r
tumor.counts <- counts[ , grepl( "-01A" , colnames( counts ) ) ]
dim(tumor.counts)
```

```
[1] 56512   481
```

```r
colnames(tumor.counts) <- sapply(strsplit(colnames(tumor.counts), "-"), function(x) paste0(x[1:3], collapse = "-"))
#write.csv(prad.biospec, "Biospecimen_prad.csv")
col.gleason<-coldata[!is.na(coldata$subtype_Clinical_Gleason_sum),]
#write.csv(coldata, "prad.coldata.csv")
```

Separate the Samples
===========================================================
class: small-code
- Separate the samples into high and low gleason score

```r
#make matrix of patients with gleason>=8
hi.gleason <- col.gleason[col.gleason$subtype_Clinical_Gleason_sum >= 8,]
dim(hi.gleason)
```

```
[1] 106 152
```

Separate the Samples
===========================================================
class: small-code
- Separate the samples into high and low gleason score

```r
#make a matrix of patients with gleason<8
lo.gleason <- col.gleason[col.gleason$subtype_Clinical_Gleason_sum < 8,]
dim(lo.gleason)
```

```
[1] 175 152
```

Separate the Samples
===========================================================
class: small-code

```r
#make a count matrix of low and high gleason score patients
lo.counts<- tumor.counts[,which(colnames(tumor.counts) %in% lo.gleason$bcr_patient_barcode)]
hi.counts<- tumor.counts[,which(colnames(tumor.counts) %in% hi.gleason$bcr_patient_barcode)]
dim(lo.counts)
```

```
[1] 56512   170
```

```r
dim(hi.counts)
```

```
[1] 56512    97
```

```r
prad.counts<-cbind(lo.counts,hi.counts)
```

edgeR
===========================================================
class: small-code

```r
#BiocManager::install("edgeR")
library(edgeR)
#citation("edgeR")
library(ggplot2)
```
-Originally developed for differential expression analysis of RNA-seq data with biological replication
  - Can also be used for a variety of count based data genomic data, such as Bisulfite-seq and ChIP-seq
- Implements statistical methods based on the negative binomial distribution as a model for count variability
  - Includes empirical Bayes methods, exact test, and GLMs

Remove zero count genes
===========================================================
class: small-code

```r
countsNozero=prad.counts[rowSums(prad.counts)!=0,]
dim(prad.counts)
```

```
[1] 56512   267
```

```r
dim(countsNozero)
```

```
[1] 53891   267
```

Prepare sample information
==========================================================
class: small-code

```r
samples <- data.frame(row.names=colnames(countsNozero),
                      Type=as.factor(c(rep("low", dim(lo.counts)[2]),rep("high",dim(hi.counts)[2]))))
```

Set up experimental design specifications
==========================================================
class: small-code

```r
Type <- factor(samples$Type)
design<-model.matrix(~0+Type)
colnames(design)<-levels(Type)
```

Set up DGE list and filter genes
===========================================================
class: small-code

```r
d.filt <- DGEList(counts=countsNozero,group=Type)
# Filter out low expression genes, and readjust lib sizes 
keep <- rowSums(cpm(d.filt)>1) >= 1
d.filt <- d.filt[keep, , keep.lib.sizes=FALSE]
dim(d.filt)
```

```
[1] 23525   267
```

Look at library sizes
==========================================================
class: small-code
- The RNA that was sequenced is called the "RNA library"
- Library size is the total number of mapped reads

```r
barplot(d.filt$samples$lib.size*1e-6, ylab="Library size (millions)")
```

<img src="AnalysisOfGeneExpressionDataUsingR-figure/unnamed-chunk-16-1.png" title="plot of chunk unnamed-chunk-16" alt="plot of chunk unnamed-chunk-16" style="display: block; margin: auto;" />


Normalize the Data
===========================================================
class: small-code

```r
# Calculate normalization factors to scale the library 
d.filt <- calcNormFactors(d.filt)
# Estimate common, trended, and then tagwise dispersions in one run
d.filt <- estimateDisp(d.filt, design)
```

plot the biological coefficient of variation
==========================================================
class: small-code

```r
plotBCV(d.filt)
```

<img src="AnalysisOfGeneExpressionDataUsingR-figure/unnamed-chunk-18-1.png" title="plot of chunk unnamed-chunk-18" alt="plot of chunk unnamed-chunk-18" style="display: block; margin: auto;" />

Make an MDS plot
===========================================================
class: small-code

```r
plotMDS(d.filt, pch = 16, labels = as.numeric(samples$Type))
```

<img src="AnalysisOfGeneExpressionDataUsingR-figure/unnamed-chunk-19-1.png" title="plot of chunk unnamed-chunk-19" alt="plot of chunk unnamed-chunk-19" style="display: block; margin: auto;" />

```r
mds <- plotMDS(d.filt)
```

<img src="AnalysisOfGeneExpressionDataUsingR-figure/unnamed-chunk-19-2.png" title="plot of chunk unnamed-chunk-19" alt="plot of chunk unnamed-chunk-19" style="display: block; margin: auto;" />

```r
plot(mds,pch = 16, col= samples$Type, xlab = "dim1",ylab = "dim2")
```

<img src="AnalysisOfGeneExpressionDataUsingR-figure/unnamed-chunk-19-3.png" title="plot of chunk unnamed-chunk-19" alt="plot of chunk unnamed-chunk-19" style="display: block; margin: auto;" />
- x- and y-axes will be representative of Euclidean distances between your samples.
- These Euclidean distances will be produced/repeated across multiple dimensions in different ways

Differential Expression Analysis
===========================================================
class: small-code

```r
#the first group is baseline
et<-exactTest(d.filt,c("low","high"))
topTags(et)
```

```
Comparison of groups:  high-low 
                    logFC     logCPM       PValue          FDR
ENSG00000168505  5.505245 -1.2969014 1.700429e-40 4.000259e-36
ENSG00000261678  3.138630 -1.8824331 8.292590e-39 9.754159e-35
ENSG00000167244  3.132083  6.3595615 8.887182e-37 6.969032e-33
ENSG00000129990  2.918731 -1.0464215 1.270540e-29 7.472366e-26
ENSG00000175832 -4.866400  5.6435064 4.394285e-28 2.067511e-24
ENSG00000163285  4.586958 -1.9654280 9.505556e-28 3.726970e-24
ENSG00000053438  3.017353  2.7307715 2.719504e-27 9.139475e-24
ENSG00000175063  1.618642  2.4149746 7.298697e-27 2.146273e-23
ENSG00000256148 -4.273862  1.1647698 5.736873e-26 1.306994e-22
ENSG00000108309  1.823040  0.8682131 6.441899e-26 1.306994e-22
```

Summarize Results
===========================================================
class: small-code

```r
summary(summaryde <- decideTestsDGE(et, p=0.05))
```

```
       high-low
Down       2874
NotSig    15027
Up         5624
```

```r
summary(summaryde <- decideTestsDGE(et, p=0.05,lfc=1))
```

```
       high-low
Down        192
NotSig    22562
Up          771
```

Correct p-values
==========================================================
class: small-code

```r
FDR <- p.adjust(et$table$PValue, method="BH")
table_et <- cbind(et$table,FDR)
```

Save DEGs in a table
==========================================================
class: small-code

```r
table.de_et <- table_et[table_et$FDR<=0.05,]
#write.csv(table.de_et,file="de_PDAC_hilo_gleason.csv")
deg_FC <- table_et[table_et$FDR<=0.05 & abs(table_et$logFC)>=1,]
#write.csv(deg_FC,file="de_PDAC_hilo_gleason_FC2.csv")
```

Make an MA plot
==========================================================
class: small-code

```r
labelled_points <- table_et[rownames(table.de_et),]

g <- ggplot(table_et , aes(y = logFC, x = logCPM)) +
  geom_point(shape = ".") +
  geom_point(colour = "red", data = labelled_points,shape = ".") +
  ggtitle("MAplot", subtitle = "high-low")
```

Make an MA plot
==========================================================
class: small-code

```r
print(g)
```

<img src="AnalysisOfGeneExpressionDataUsingR-figure/unnamed-chunk-25-1.png" title="plot of chunk unnamed-chunk-25" alt="plot of chunk unnamed-chunk-25" style="display: block; margin: auto;" />

Downstream analysis
=========================================================
- Observe in genome browser
- Annotate DEG
- Find biological processes that are enriched amongst the DEGs
- Pathway analysis
- Clustering analyses
- Biological Validation
- Network analyses

Annotate genes using biomaRt
==========================================================
class: small-code

```r
library(biomaRt)
listMarts()
```

```
               biomart               version
1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 98
2   ENSEMBL_MART_MOUSE      Mouse strains 98
3     ENSEMBL_MART_SNP  Ensembl Variation 98
4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 98
```

```r
#choose ensembl dataset and filters to use
ensembl=useMart("ensembl")
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = 'www.ensembl.org')
filters = listFilters(ensembl)
filters[1:5,]
```

```
             name              description
1 chromosome_name Chromosome/scaffold name
2           start                    Start
3             end                      End
4      band_start               Band Start
5        band_end                 Band End
```
Annotate genes using biomaRt
==========================================================
class: small-code

```r
attributes = listAttributes(ensembl)
attributes[1:15,]
```

```
                            name                  description         page
1                ensembl_gene_id               Gene stable ID feature_page
2        ensembl_gene_id_version       Gene stable ID version feature_page
3          ensembl_transcript_id         Transcript stable ID feature_page
4  ensembl_transcript_id_version Transcript stable ID version feature_page
5             ensembl_peptide_id            Protein stable ID feature_page
6     ensembl_peptide_id_version    Protein stable ID version feature_page
7                ensembl_exon_id               Exon stable ID feature_page
8                    description             Gene description feature_page
9                chromosome_name     Chromosome/scaffold name feature_page
10                start_position              Gene start (bp) feature_page
11                  end_position                Gene end (bp) feature_page
12                        strand                       Strand feature_page
13                          band               Karyotype band feature_page
14              transcript_start        Transcript start (bp) feature_page
15                transcript_end          Transcript end (bp) feature_page
```

```r
filterOptions("ensembl_gene_id",ensembl)
```

```
[1] "[]"
```

Annotate genes using biomaRt
==========================================================
class: small-code

```r
background_genes <- rownames(table_et)
degs <- rownames(table.de_et)

bckgd.annot <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol',                                'entrezgene_id','description','external_gene_name',
                'chromosome_name','start_position','end_position','strand'), 
                filters='ensembl_gene_id',values=background_genes,
                mart=ensembl)
```

Annotate genes using biomaRt
==========================================================
class: small-code

```r
deg.annot <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','entrezgene_id',                  'description','external_gene_name','chromosome_name','start_position','end_position','strand'),
                filters='ensembl_gene_id', values=degs, mart=ensembl)
colnames(bckgd.annot)
```

```
[1] "ensembl_gene_id"    "hgnc_symbol"        "entrezgene_id"     
[4] "description"        "external_gene_name" "chromosome_name"   
[7] "start_position"     "end_position"       "strand"            
```

Annotate genes using biomaRt
==========================================================
class: small-code

```r
deg_ann.tb = merge(x = table_et, y = bckgd.annot, by.x="row.names",by.y='ensembl_gene_id',all.x=TRUE)
bckgd_ann.tb = merge(x = table.de_et, y = deg.annot, by.x="row.names",by.y='ensembl_gene_id',all.x=TRUE)
#write.csv(bckgd_ann.tb, file="allgenesAnnotation.csv")
#write.csv(deg_ann.tb, file="DEGAnnotation.csv")
```

Perform a pathway analysis and enrichment analysis
==========================================================
class: small-code

```r
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
```

Gene Ontology
==========================================================
- Ontology: a structured vocabulary
  - Describes concepts that exist in an area of knowledge
  - Describes relationships that exist
- Gene Ontology (GO):
  - Describes possible functions of genes
  - Describes relationships between genes
  - Describes where in the cell gene products are localized
- GO is independent of organism!

Gene Ontology
==========================================================
- structure: DAG 
  - directed acyclic graph
- relationships:
  - is a
  - has a
  - part of
<div class="midcenter" style="margin-left:-45px; margin-top:-200px;">
<img src="./AnalysisOfGeneExpressionDataUsingR-figure/BPgo.png"></img>
</div>

Perform a GO enrichment analysis 
==========================================================
class: small-code

```r
cc.go <- enrichGO(gene          = deg_ann.tb$entrezgene_id,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(cc.go)
```

```
                   ID                      Description GeneRatio   BgRatio
GO:0030055 GO:0030055          cell-substrate junction 410/15826 412/19717
GO:0005924 GO:0005924 cell-substrate adherens junction 406/15826 408/19717
GO:0005925 GO:0005925                   focal adhesion 403/15826 405/19717
GO:0010008 GO:0010008                endosome membrane 468/15826 479/19717
GO:0005759 GO:0005759             mitochondrial matrix 457/15826 469/19717
GO:0005774 GO:0005774                vacuolar membrane 403/15826 412/19717
                 pvalue     p.adjust       qvalue
GO:0030055 8.727243e-37 6.815977e-34 2.884583e-34
GO:0005924 2.104809e-36 8.219280e-34 3.478474e-34
GO:0005925 4.072408e-36 1.060183e-33 4.486793e-34
GO:0010008 8.895154e-32 1.736779e-29 7.350206e-30
GO:0005759 6.520861e-30 1.018559e-27 4.310632e-28
GO:0005774 6.174534e-28 8.037186e-26 3.401410e-26
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   geneID
GO:0030055                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           LAP3/CD99/LASP1/ABCB4/ITGA3/ITGA2B/RALA/CD9/MRC2/TSPAN9/PLAUR/EHD3/CAPN1/FHL1/VIM/CD44/ARHGAP31/VCL/TNC/CTNNA1/HSPA5/LIMA1/BCAR1/CYBA/SYNE2/GDI2/PPP1R12A/NCKAP1/RPL18/CNN2/SLC9A3R2/COL17A1/TLE2/RHOA/FGFR3/PABPC1/CDC42/MAP4K4/RPL31/ACTN1/LIMS2/PVR/FERMT2/CLASP1/HACD3/ACTB/REXO2/MCAM/USP33/APBB1IP/ACTN2/ITGA8/FAP/TNS1/SENP1/DNM2/EPB41L2/RAB21/PTPRC/ITGB5/RPS5/FAT1/RAB10/CD59/CPNE3/CTTN/NOX4/TRIP6/ADD1/CASS4/ASAP3/RPL6/RPLP0/PXN/SLC9A1/ICAM1/ITGA6/SMPX/SNAP23/EZR/SORBS1/JAK2/NRP1/MISP/MAPK1/TRIOBP/PACSIN2/RPL3/MYH9/ZFYVE21/PROCR/FERMT1/HCK/MAPRE1/CD99L2/ARHGEF7/FLT1/MAPK3/CORO2B/RPS16/RPS19/ITGB8/CAV2/CAV1/HSPB1/LIMK1/ENG/PDLIM1/GIT1/RPL19/PFN1/SLC6A4/YWHAE/LAMTOR3/HSPA8/LPXN/CBL/CD81/RPS13/PPFIBP1/CORO1C/TNS2/TRPV4/ARPC3/NEDD9/WASF1/PTK7/ERBIN/HSPA9/PDGFRB/PRKAR2A/ACTR3/EPB41L5/ITGB6/ITGA4/RPS15/FHL2/RND3/RPL22/RHOU/ARHGEF2/DOCK7/CD46/CNN3/GNA13/TEK/SORBS3/PTK2B/CAT/RPL5/PLAU/ADGRE5/LRP1/SDC4/STX16/RPS10/AHNAK/EFNB2/RPL23/VASP/FLRT3/RRAS/FLRT1/AIF1L/NUP214/MAP2K2/PTPN12/CDC42EP1/RAC2/FLNC/ARHGAP22/PALLD/AJUBA/STARD8/CNN1/NECTIN2/ACTN4/ARPC1B/PAK4/AKAP12/CAP1/RPL27/PPFIA1/TNS4/ITGB4/FLOT2/DCTN4/MPRIP/KRAS/RRAS2/NUMB/YWHAQ/ANXA1/TES/FLNB/LMO7/LCP1/TNS3/RAC1/CIB2/ARPC5L/TLN1/HMGA1/FLOT1/MDC1/ATAT1/SDCBP/RDX/KIF23/ITGA11/RPLP1/ADAM10/ACTR2/ITGAV/ARHGAP24/SCARB2/PARVG/GIT2/ITGB7/IQGAP1/TGFB1I1/ARMC5/CDH13/RPS2/PDPK1/CLTC/GRB7/RPS11/RPL13A/EPHA2/HSPG2/RPS8/MTF2/DCAF6/PRUNE1/PIP5K1A/S100A7/ARF1/RHOB/FBLN7/LIMD1/PHLDB2/LPP/RPS3A/ARHGAP26/G3BP1/GNA12/EGFR/SH3KBP1/CASK/MSN/ZNF185/RPL7/GSN/RPL7A/RSU1/CAPN5/PAK1/RPS3/HYOU1/ITGB1/DIXDC1/TWF1/ADAM17/DST/ARL14EP/TADA1/GJA1/DAB2/THY1/PGM5/ENAH/SORBS2/RPL30/MMP14/FZD1/CSRP1/ACTC1/ZYX/ITGB2/RPL8/NPHS1/ITGA5/JAK1/FBLIM1/NEXN/ARPC5/DDR2/NCSTN/CAPN2/XIRP2/ARPC2/NFASC/CLASP2/RPL9/ANXA5/ITGA2/RPS14/DLC1/SLC4A2/YWHAZ/HNRNPK/ARF6/ILK/HSP90B1/B2M/PPIB/YWHAB/MAPRE2/PDIA3/TPM4/SRP68/CTNNB1/FAM107A/IRF2/XIRP1/ADAM9/SNTB2/TM4SF20/MAP2K1/PTK2/LIMS1/ALCAM/YWHAG/PDCD6IP/CDH2/RPS9/RPS7/TLN2/KLF11/SNTB1/BSG/GNB2/SYNPO2/CORO1B/CFL1/RPL38/DAG1/PEAK1/CSPG4/JUP/RPL4/CSRP2/YES1/RHOG/RPLP2/CD151/FLII/PLEC/GAK/CALR/FZD2/NPM1/ADGRB1/FES/CAV3/RPS17/FHL3/ACTG1/UBOX5/FLRT2/LMLN/P4HB/ATP6V0C/PIP5K1C/PPP1CC/CHP1/SPRY4/NHS/FOCAD/PARVB/PPIA/EVL/AFAP1/MME/PDLIM7/FLNA/ANXA6/IGF2R/PCBP2/SVIL/DPP4/PARVA/RPL37A/NRAP/RPL12/MPZL1/RPS4X/ITGBL1/RPL10A/L1CAM/DMD/TGM2/LAYN/HSPA1B/HSPA1A/ARL2/PPP1CB/RPS29/ITGA1/TSPAN4/RPS18/ALKBH6/PI4KA/SCARF2/ITGB3/EPPK1/CYFIP1/PRAG1/MARCKS
GO:0005924                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   LAP3/CD99/LASP1/ABCB4/ITGA3/ITGA2B/RALA/CD9/MRC2/TSPAN9/PLAUR/EHD3/CAPN1/FHL1/VIM/CD44/ARHGAP31/VCL/TNC/CTNNA1/HSPA5/LIMA1/BCAR1/CYBA/SYNE2/GDI2/PPP1R12A/NCKAP1/RPL18/CNN2/SLC9A3R2/TLE2/RHOA/FGFR3/PABPC1/CDC42/MAP4K4/RPL31/ACTN1/LIMS2/PVR/FERMT2/CLASP1/HACD3/ACTB/REXO2/MCAM/USP33/APBB1IP/ACTN2/ITGA8/FAP/TNS1/SENP1/DNM2/EPB41L2/RAB21/PTPRC/ITGB5/RPS5/FAT1/RAB10/CD59/CPNE3/CTTN/NOX4/TRIP6/ADD1/CASS4/ASAP3/RPL6/RPLP0/PXN/SLC9A1/ICAM1/ITGA6/SMPX/SNAP23/EZR/SORBS1/JAK2/NRP1/MISP/MAPK1/TRIOBP/PACSIN2/RPL3/MYH9/ZFYVE21/PROCR/FERMT1/HCK/MAPRE1/CD99L2/ARHGEF7/FLT1/MAPK3/CORO2B/RPS16/RPS19/ITGB8/CAV2/CAV1/HSPB1/LIMK1/ENG/PDLIM1/GIT1/RPL19/PFN1/SLC6A4/YWHAE/LAMTOR3/HSPA8/LPXN/CBL/CD81/RPS13/PPFIBP1/CORO1C/TNS2/TRPV4/ARPC3/NEDD9/WASF1/PTK7/HSPA9/PDGFRB/PRKAR2A/ACTR3/EPB41L5/ITGB6/ITGA4/RPS15/FHL2/RND3/RPL22/RHOU/ARHGEF2/DOCK7/CD46/CNN3/GNA13/TEK/SORBS3/PTK2B/CAT/RPL5/PLAU/ADGRE5/LRP1/SDC4/STX16/RPS10/AHNAK/EFNB2/RPL23/VASP/FLRT3/RRAS/FLRT1/AIF1L/NUP214/MAP2K2/PTPN12/CDC42EP1/RAC2/FLNC/ARHGAP22/PALLD/AJUBA/STARD8/CNN1/NECTIN2/ACTN4/ARPC1B/PAK4/AKAP12/CAP1/RPL27/PPFIA1/TNS4/ITGB4/FLOT2/DCTN4/MPRIP/KRAS/RRAS2/NUMB/YWHAQ/ANXA1/TES/FLNB/LMO7/LCP1/TNS3/RAC1/CIB2/ARPC5L/TLN1/HMGA1/FLOT1/MDC1/ATAT1/SDCBP/RDX/KIF23/ITGA11/RPLP1/ADAM10/ACTR2/ITGAV/ARHGAP24/SCARB2/PARVG/GIT2/ITGB7/IQGAP1/TGFB1I1/ARMC5/CDH13/RPS2/PDPK1/CLTC/GRB7/RPS11/RPL13A/EPHA2/HSPG2/RPS8/MTF2/DCAF6/PRUNE1/PIP5K1A/S100A7/ARF1/RHOB/FBLN7/LIMD1/PHLDB2/LPP/RPS3A/ARHGAP26/G3BP1/GNA12/EGFR/SH3KBP1/CASK/MSN/ZNF185/RPL7/GSN/RPL7A/RSU1/CAPN5/PAK1/RPS3/HYOU1/ITGB1/DIXDC1/TWF1/ADAM17/DST/ARL14EP/TADA1/GJA1/DAB2/THY1/PGM5/ENAH/SORBS2/RPL30/MMP14/FZD1/CSRP1/ACTC1/ZYX/ITGB2/RPL8/NPHS1/ITGA5/JAK1/FBLIM1/NEXN/ARPC5/DDR2/NCSTN/CAPN2/XIRP2/ARPC2/NFASC/CLASP2/RPL9/ANXA5/ITGA2/RPS14/DLC1/SLC4A2/YWHAZ/HNRNPK/ARF6/ILK/HSP90B1/B2M/PPIB/YWHAB/MAPRE2/PDIA3/TPM4/SRP68/CTNNB1/FAM107A/IRF2/XIRP1/ADAM9/SNTB2/TM4SF20/MAP2K1/PTK2/LIMS1/ALCAM/YWHAG/PDCD6IP/CDH2/RPS9/RPS7/TLN2/KLF11/SNTB1/BSG/GNB2/SYNPO2/CORO1B/CFL1/RPL38/DAG1/PEAK1/CSPG4/JUP/RPL4/CSRP2/YES1/RHOG/RPLP2/CD151/FLII/PLEC/GAK/CALR/FZD2/NPM1/ADGRB1/FES/CAV3/RPS17/FHL3/ACTG1/UBOX5/FLRT2/LMLN/P4HB/ATP6V0C/PIP5K1C/PPP1CC/CHP1/SPRY4/NHS/FOCAD/PARVB/PPIA/EVL/AFAP1/MME/PDLIM7/FLNA/ANXA6/IGF2R/PCBP2/SVIL/DPP4/PARVA/RPL37A/NRAP/RPL12/MPZL1/RPS4X/ITGBL1/RPL10A/L1CAM/TGM2/LAYN/HSPA1B/HSPA1A/ARL2/PPP1CB/RPS29/ITGA1/TSPAN4/RPS18/ALKBH6/PI4KA/SCARF2/ITGB3/CYFIP1/PRAG1/MARCKS
GO:0005925                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  LAP3/CD99/LASP1/ABCB4/ITGA3/ITGA2B/RALA/CD9/MRC2/TSPAN9/PLAUR/EHD3/CAPN1/FHL1/VIM/CD44/ARHGAP31/VCL/TNC/CTNNA1/HSPA5/LIMA1/BCAR1/CYBA/SYNE2/GDI2/PPP1R12A/NCKAP1/RPL18/CNN2/SLC9A3R2/TLE2/RHOA/FGFR3/PABPC1/CDC42/MAP4K4/RPL31/ACTN1/LIMS2/PVR/FERMT2/CLASP1/HACD3/ACTB/REXO2/MCAM/USP33/APBB1IP/ACTN2/ITGA8/FAP/TNS1/SENP1/DNM2/EPB41L2/RAB21/PTPRC/ITGB5/RPS5/FAT1/RAB10/CD59/CPNE3/CTTN/NOX4/TRIP6/ADD1/CASS4/ASAP3/RPL6/RPLP0/PXN/SLC9A1/ICAM1/ITGA6/SNAP23/EZR/SORBS1/JAK2/NRP1/MISP/MAPK1/TRIOBP/PACSIN2/RPL3/MYH9/ZFYVE21/PROCR/FERMT1/HCK/MAPRE1/CD99L2/ARHGEF7/FLT1/MAPK3/CORO2B/RPS16/RPS19/ITGB8/CAV2/CAV1/HSPB1/LIMK1/ENG/PDLIM1/GIT1/RPL19/PFN1/SLC6A4/YWHAE/LAMTOR3/HSPA8/LPXN/CBL/CD81/RPS13/PPFIBP1/CORO1C/TNS2/TRPV4/ARPC3/NEDD9/WASF1/PTK7/HSPA9/PDGFRB/PRKAR2A/ACTR3/EPB41L5/ITGB6/ITGA4/RPS15/FHL2/RND3/RPL22/RHOU/ARHGEF2/DOCK7/CD46/CNN3/GNA13/TEK/SORBS3/PTK2B/CAT/RPL5/PLAU/ADGRE5/LRP1/SDC4/STX16/RPS10/AHNAK/EFNB2/RPL23/VASP/FLRT3/RRAS/FLRT1/AIF1L/NUP214/MAP2K2/PTPN12/CDC42EP1/RAC2/FLNC/ARHGAP22/PALLD/AJUBA/STARD8/CNN1/NECTIN2/ACTN4/ARPC1B/PAK4/AKAP12/CAP1/RPL27/PPFIA1/TNS4/ITGB4/FLOT2/DCTN4/MPRIP/KRAS/RRAS2/NUMB/YWHAQ/ANXA1/TES/FLNB/LMO7/LCP1/TNS3/RAC1/ARPC5L/TLN1/HMGA1/FLOT1/MDC1/ATAT1/SDCBP/RDX/KIF23/ITGA11/RPLP1/ADAM10/ACTR2/ITGAV/ARHGAP24/SCARB2/PARVG/GIT2/ITGB7/IQGAP1/TGFB1I1/ARMC5/CDH13/RPS2/PDPK1/CLTC/GRB7/RPS11/RPL13A/EPHA2/HSPG2/RPS8/MTF2/DCAF6/PRUNE1/PIP5K1A/S100A7/ARF1/RHOB/FBLN7/LIMD1/PHLDB2/LPP/RPS3A/ARHGAP26/G3BP1/GNA12/EGFR/SH3KBP1/CASK/MSN/ZNF185/RPL7/GSN/RPL7A/RSU1/CAPN5/PAK1/RPS3/HYOU1/ITGB1/DIXDC1/TWF1/ADAM17/DST/ARL14EP/TADA1/GJA1/DAB2/THY1/PGM5/ENAH/SORBS2/RPL30/MMP14/FZD1/CSRP1/ACTC1/ZYX/ITGB2/RPL8/NPHS1/ITGA5/JAK1/FBLIM1/NEXN/ARPC5/DDR2/NCSTN/CAPN2/XIRP2/ARPC2/NFASC/CLASP2/RPL9/ANXA5/ITGA2/RPS14/DLC1/SLC4A2/YWHAZ/HNRNPK/ARF6/ILK/HSP90B1/B2M/PPIB/YWHAB/MAPRE2/PDIA3/TPM4/SRP68/CTNNB1/FAM107A/IRF2/XIRP1/ADAM9/SNTB2/TM4SF20/MAP2K1/PTK2/LIMS1/ALCAM/YWHAG/PDCD6IP/CDH2/RPS9/RPS7/TLN2/KLF11/SNTB1/BSG/GNB2/SYNPO2/CORO1B/CFL1/RPL38/DAG1/PEAK1/CSPG4/JUP/RPL4/CSRP2/YES1/RHOG/RPLP2/CD151/FLII/PLEC/GAK/CALR/FZD2/NPM1/ADGRB1/FES/CAV3/RPS17/FHL3/ACTG1/UBOX5/FLRT2/LMLN/P4HB/ATP6V0C/PIP5K1C/PPP1CC/CHP1/SPRY4/NHS/FOCAD/PARVB/PPIA/EVL/AFAP1/MME/PDLIM7/FLNA/ANXA6/IGF2R/PCBP2/SVIL/DPP4/PARVA/RPL37A/RPL12/MPZL1/RPS4X/ITGBL1/RPL10A/L1CAM/TGM2/LAYN/HSPA1B/HSPA1A/ARL2/PPP1CB/RPS29/ITGA1/TSPAN4/RPS18/ALKBH6/PI4KA/SCARF2/ITGB3/CYFIP1/PRAG1/MARCKS
GO:0010008 CFTR/SPPL2B/LAMP2/AP2B1/VPS41/BAIAP3/VTA1/STARD3NL/CLCN6/EHD3/SLC11A1/EHD2/HSD17B6/SNX1/RABEP1/ATP6V0A1/CLEC16A/ZFYVE16/RAB27B/AP2S1/DTNBP1/VAMP3/AP5M1/ATP9A/TAB2/MCOLN3/ATP11B/YIPF1/PSD/ABCA7/AP3D1/SLC9A7/STOML1/GRIPAP1/LAPTM4A/VPS35/RAB27A/LRP6/WIPI1/EPHA8/SNX13/ATP6AP1/TFRC/ACAP1/CLCN4/TSG101/PLD1/RAB7A/LAMP3/ITCH/STX7/RAB21/PSEN1/RNF13/CHMP2B/APOB/RAB10/SCAMP1/EPS15/CHMP5/SNX10/MTMR2/SNX5/KIF16B/CDIP1/SPG21/RAB11FIP3/MCOLN1/CMTM6/TF/RFFL/TSPAN15/WASHC2A/MARCHF2/LZTR1/GGA1/MICALL1/SUN2/PACSIN2/TAB1/NCF4/GNPNAT1/VTI1B/SLA2/CHMP4B/SYNDIG1/TLR8/EEA1/NDFIP2/VAC14/GGA2/RAB11A/EHD4/OCA2/RHOV/VPS18/SLC30A4/ZDHHC2/LAPTM4B/NDRG1/SNX16/CLIP3/PLIN3/WDR91/ATP6V0A4/CAV1/SNX8/TMEM106B/SLC1A1/TYRP1/ABCA2/ABHD17B/RAB11FIP2/CUBN/DKK1/MTMR4/GOSR2/SLC6A4/RAB5C/MMD/LAMTOR3/UBE2D3/CLCN3/SNX25/EHD1/TCIRG1/CORO1C/SLC11A2/VPS29/RAB5B/RAB35/RAB23/SNX3/FIG4/ATP6V0E1/ACAP2/SNX4/PIKFYVE/STEAP3/STAM2/SLC30A3/SNX17/GRB14/CHMP3/ABCB6/PLEKHB2/CD207/EPHA4/TMEM59/KIAA1324/SCAMP3/LAMTOR2/WLS/PLEKHM2/TMEM9/ECE1/APH1A/ATP6V0B/STX12/OSBPL9/CTSD/MREG/VAMP8/RAB14/VPS4B/CD274/SNX19/OCRL/CLTA/VPS26A/WDR83/OPTN/RAP2C/SNX21/RAB22A/PMEPA1/VAMP7/RAB17/RAP2A/TM9SF2/GGA3/AP5S1/LAMP5/BECN1/TMEM175/TICAM1/STEAP4/MYO1B/CD68/SNX6/ABHD17A/SYT5/LDLR/CHMP2A/YIPF2/UBA1/CHMP1A/RAB11FIP4/TBC1D5/RBSN/PDLIM4/VPS25/NDFIP1/KREMEN2/WDR44/STARD3/LLGL1/VPS4A/SBF2/IRAK2/ARL8B/SORT1/LAMTOR5/TMEM165/ARHGAP32/ANXA1/SNX14/MAP3K7/CD63/LTV1/CD164/RAB11FIP5/TMBIM1/SCYL2/APPL2/VPS36/ITM2B/GPNMB/RAC1/ABHD17C/VPS45/STAM/TLR4/RABEPK/KIF13A/KIAA0319/RIPK1/SORL1/DNAJC13/USP8/SPPL2A/APH1B/SCARB2/VPS33A/VPS37B/RAB15/PML/SCAMP2/FURIN/SH3GL3/TOM1L1/VPS53/CLTC/RMC1/NPC1/ERBB2/MVB12A/SH3GL1/IFITM3/SNX27/KCNH1/RAB13/ARL8A/RHOB/RPS27A/RAB5A/TMEM108/OSBPL11/MARCHF1/FGD2/CLVS2/EGFR/CHMP7/SLC26A7/ATP6V0D2/SLC39A4/NTRK2/LAMTOR1/FCGR1A/UBC/VTI1A/RNF144A/TMEM163/MCOLN2/ABCA5/WNT3A/LY96/RABGEF1/EPHB1/PIP4P2/ZFYVE27/GRIA1/FZD7/VPS37A/ZFYVE9/STEAP2/APPL1/TAB3/MITD1/CD1D/CD1A/CD1C/CD1B/SNF8/ATP13A2/ATP6V0D1/ZFYVE28/VPS11/VPS28/NAPEPLD/AP2M1/CD300LG/CYB561A3/TPCN2/NCSTN/FZD5/ANTXR2/ABHD6/DTX3L/SLC9B2/TLR3/STEAP1/CHMP4C/TMEM184A/DCSTAMP/UBAP1/MARCHF8/ARF6/RET/PIP4P1/AMN/RAB8B/HPS6/PLEKHF1/B2M/VPS39/PDIA3/SNX20/OR51E2/RAB8A/UBXN6/RILP/WDR81/CD320/VPS37C/RAB4A/RAB31/LDLRAD4/NSG1/PARM1/PCSK9/ADRB2/MMGT1/ANTXR1/NSG2/STX8/UBB/CD14/INSR/ATP6V0E2/CLCN5/AQP4/CD8B/WASHC2C/MYD88/RHOD/COMMD1/GOLIM4/MARCHF3/REP15/TRAF6/ARHGAP1/TMEM9B/PLEKHF2/CHMP6/SPHK1/VPS37D/BOK/RUFY1/SLC38A9/CLVS1/SNX18/RCC2/HLA-DQB1/TMEM150B/ZNRF2/GPR62/WASHC1/RAP2B/GPR135/SLC9A9/ANXA2/AP2A2/TBK1/VPS33B/IRAK1/SORCS2/RAB11B/ATP6V0A2/HGS/IRF7/PMEL/ANKFY1/ATP6V0C/LAMP1/PIP5K1C/BACE1/TPCN1/LAMTOR4/KIR2DL4/LITAF/HLA-DRB1/GIMAP5/TLR7/HLA-DQA1/MVB12B/AP2A1/ANXA6/DIO3/IRAK4/FCGR1B/SLC29A3/NTRK1/HLA-DRB5/ARC/SLC9A6/SCAMP5/INPP5F/ATG9A/INPP5B/HLA-DOA/HLA-DMA/HLA-DRA/HLA-C/HLA-E/HLA-G/HLA-F/PSENEN/SNX2/RAB12/HLA-A/SLC48A1/LEPROT/AP1G2/VPS16/UBA52/VPS52/HLA-DPB1/CPTP/PLEKHM1/SCAMP4/HLA-DPA1/HLA-DQB2/HLA-B/SPAAR/HLA-DQA2/HLA-DOB/HLA-DMB/PRAF2/PLA2G4B/UBAP1L/CHMP4A/CHMP1B/MRC1/IKBKE/ANXA8/TBC1D3/RAB7B
GO:0005759                                                                                                                                                                                                  NDUFAF7/POLDIP2/NDUFAB1/PDK4/SLC25A5/ACSM3/PDK2/ELAC2/AASS/LARS2/CPS1/NDUFS1/ALAS1/GLRX2/MIPEP/TFB1M/TRIT1/MRPS10/RAD51/MRPL43/MRPS35/CS/MRPS24/HAGH/OAT/KARS1/MTHFD2/ATXN3/FECH/IARS2/IDH3G/PDK3/COASY/HSD17B10/TRNT1/ACADVL/MRPS34/ACAT1/REXO2/MCCC1/ARG2/ME2/MRPL22/BCKDHB/OXCT1/HADHA/DIMT1/MRPL28/SIRT4/MRPS33/PDPR/DLD/FH/MRPS18A/ATP5F1D/POLRMT/PRODH/MCAT/TXN2/ACO2/GSTZ1/DHRS2/PCK2/PRORP/MTG2/IDH3B/PIN4/MRPS31/ARL2BP/MLYCD/MPG/NME4/EARS2/DNAJA3/BCKDK/DECR1/GSR/SARS2/TIMM44/MRPL4/ETFB/BCAT2/GCDH/ETHE1/PMPCB/OGDH/SSBP1/HIBADH/GARS1/PTCD1/NUDT1/MRPL32/TWNK/PITRM1/TFAM/PPIF/C1QBP/MRPL27/LRRC59/GRPEL1/AADAT/PDHX/NDUFS8/COQ5/ATP5F1B/ACSS3/ACAD10/ALDH2/MRPL51/MTRF1L/SOD2/MRPL18/ALDH5A1/MRPL2/HARS2/MRPS30/HSPA9/MRPS27/NR3C1/PCCB/MRPL3/DGUOK/NDUFS7/ACADL/MRPL19/GLS/HSPE1/MRPL37/PARK7/MECR/ATP5PB/ADPRHL2/WARS2/MRPS15/ACADM/HMGCL/NSUN4/DARS2/FASTKD2/CREB1/FOXO3/ACOT2/DLST/ALDH6A1/MTHFD1L/MRPS14/TP53AIP1/DNAJC15/MTERF2/LIAS/MTERF4/MRPS2/MRM2/ACADS/ACOT9/TSFM/SARDH/CHPF/ATP5F1E/MCEE/SIRT5/GOT2/MRPS7/ALKBH7/CLPP/MRPS26/PRDX5/TRAP1/TRMT5/ECHS1/STYXL1/MTERF1/MPST/TST/MRPS12/IVD/MRPL34/NDUFA10/ATG4D/C12orf65/ACSS2/PDHA1/MCCC2/DHX30/MRPL35/COQ3/GRSF1/RIDA/ERAL1/DAP3/DMGDH/DGLUCY/MRPS36/CCNB1/HMGCS2/CARS2/ISCA1/GLS2/BLOC1S1/RNASEL/MRPL44/CYP27A1/MRPS9/ISCU/ALDH1L2/SUCLA2/TBRG4/MTHFS/MRPL47/FPGS/MRPL50/ALDH1B1/UQCC2/NARS2/MRPL15/FDX1/NDUFAF1/DBT/HADHB/LRPPRC/DNA2/PPA2/HADH/YARS2/ETFBKMT/NDUFA9/MMAB/SUOX/ETFA/CYP11A1/POLG/UQCRC2/GCSH/TP53/SIRT3/SOD1/MRPL24/CASQ1/TARS2/MRPL9/PYCR2/MRPS5/LIPT1/HSPD1/ABHD10/AMT/SNCA/CBR4/FARS2/MMUT/RARS2/MDH2/ADHFE1/MRPS28/LACTB2/AK3/AUH/MRRF/PDSS1/GLUD1/MTG1/RPS3/MRPL49/DLAT/ME3/NUBPL/ACAD8/MMAA/GUF1/ATP5F1A/PDK1/NADK2/SDHAF4/MRPL39/ACSS1/OXA1L/SUPV3L1/MALSU1/RPUSD3/MRPL17/ALAS2/NDUFS2/CCAR2/MRPL10/ALDH4A1/THEM4/FLAD1/SHC1/BDH1/FDXR/NAGS/ACOT11/PARS2/AK4/ACP6/TFB2M/MRPL55/MAIP1/MRPS18C/NAXE/SUCLG1/PPM1K/MTHFD2L/NMNAT3/ETNPPL/ABCE1/GRPEL2/PRIMPOL/GFM2/TERT/PDSS2/FASTK/ALDH7A1/PDP1/NUDT2/FXN/RPUSD4/ATP5F1C/VDAC2/PRDX3/PMPCA/ISCA2/GPT2/NDUFB8/NUDT13/IDH3A/TK2/ACSM1/CLPX/MRPL16/ACSF2/TRUB2/ACAA2/MRPL58/ECI1/SDHAF2/PDHB/DTYMK/GFM1/CA5B/BTD/CDK1/NUDT9/TRMT61B/MRPL36/ETFDH/BCL2L1/MRM3/QARS1/TEFM/MRPL13/SUCLG2/AGXT/CHCHD1/MRPL52/PDP2/ACSM6/MRPL57/PC/TRMT10C/MRPL11/PABPC5/PDE12/MRPS22/PCCA/PHYKPL/LIPT2/MRPL48/ACSF3/TYMS/PUS1/TOP3A/GLDC/ERBB4/NSUN3/DHFR2/TUFM/GADD45GIP1/FAHD1/D2HGDH/MRPL14/DHTKD1/IBA57/MRPS11/IDH2/MRPL41/MRPS16/SHMT2/GLRX5/DDX28/PYCR1/ABAT/SMDT1/ACSM5/MRPL54/ACSM2A/TOP1MT/TXNRD2/MRPL30/MRPL40/NAT8L/TMLHE/PDE2A/LYRM7/HYKK/LRRK2/TDRD7/ACADSB/LONP1/THEM5/SDHAF3/PPTC7/MRPL21/OGDHL/GSTK1/BCO2/ATAD3A/MRPL42/HIBCH/OXCT2/HSD17B8/MRPL38/HSPA1L/MRPS18B/MRPL53/SDHAF1/NT5M/DNAJC19/ARL2/VDAC1/NDUFS3/NAXD/MRPL23/LYRM4/ACSM4/FASTKD5/PAM16/PARG/AKR1B15/GPX1/MRPS17/HOGA1/MRPL20/MRPL33/MRPS6/NFS1/MARS2/BCKDHA/MPV17L2/POLG2/MRPL46/MRPL12/MRPS21/FDX2/NDUFA7/MRM1/MRPL45
GO:0005774                                                                                                                                                                                                                                                                                                                                                                        CFTR/M6PR/SPPL2B/LAMP2/AP2B1/VPS41/CYB561/SPAG9/MGST1/STARD3NL/CLCN6/SYT7/MAN2B2/SLC7A14/ACPP/CD74/RRAGD/GRN/ATP6V0A1/GABARAPL2/CLEC16A/SLC66A1/CTNS/AP2S1/ATP6V1H/CP/GOPC/AP5M1/MCOLN3/ATP11B/VMP1/SPHK2/CTSA/AP3D1/GNAI3/ATP11A/LAPTM4A/VPS35/WIPI1/CYBRD1/ZFYVE26/AP1M1/SCARB1/NSF/PLD1/RAB7A/GPR137B/LAMP3/GNB1/STX7/PSEN1/OSTM1/LRP2/RNF13/RRAGB/SEH1L/CPNE3/ATG16L1/B4GALT1/CEACAM6/MTMR2/GNA11/CDIP1/MCOLN1/ABCC6/CMTM6/SLC22A17/HSP90AB1/SH3GLB1/MARCHF2/SNAP29/DEPDC5/AP1B1/TOM1/SYNGR1/ATP6V1D/VTI1B/TM9SF1/DNAJC5/EEF1A2/MAP1LC3A/VAPA/GPR143/TLR8/ATP11C/ABCD1/MAGT1/RUBCNL/CLN5/VAC14/WDR59/NPRL3/ABCC1/CLCN7/OCA2/SPG11/VPS18/SLC30A4/LAPTM4B/RAB2A/MAN2B1/NAPA/RAB3D/KXD1/ATP6V0A4/AP1S1/TMEM106B/RHEB/ABCA2/CUBN/RAB5C/ABCC3/MMD/LAMTOR3/MANBA/HSPA8/SLC15A3/TCIRG1/SLC11A2/ITFG2/GPLD1/VNN1/TMEM30A/LNPEP/TRIM23/NPRL2/GNB4/ATP6V1A/SLC30A3/EVA1A/ABCB6/ATP6V1B1/TMEM59/KIAA1324/LAMTOR2/TMEM9/RRAGC/ECE1/ATP6V0B/TSPAN1/CTSD/KPTN/MREG/LRMP/VAMP8/RPN2/RAB14/ABCD4/SLC17A5/WDR11/ABCC11/CCZ1/CLTA/LRP1/SLC36A1/SLC12A4/VAMP7/ATP8A1/ABCC10/AHNAK/NAPB/AP5S1/LAMP5/ATG14/RAP1B/TMEM175/WDR24/GNAI1/ATP6V1F/CD68/SLC44A2/AP1M2/PGAP6/LDLR/BST2/PNPLA7/UBA1/ATP6V1E1/STARD3/VPS4A/AP3B1/SLC39A11/SBF2/ARL8B/SORT1/LAMTOR5/NAPG/ACP2/TMEM165/P2RX4/SNX14/CD63/BLOC1S1/LTV1/CD164/ITM2C/TMBIM1/CKAP4/DRAM1/TM6SF1/CALCOCO2/CCDC115/SLC2A8/STX17/ATP6V1G1/FLOT1/MYO7A/ARRB1/PRCP/SLC3A1/ATRAID/DNAJC13/SLC49A4/SPPL2A/SCARB2/ENPEP/EGF/GABARAPL1/MARCHF9/GLIPR1/SLC15A4/VPS33A/MAP1LC3B/MEAK7/CLTC/RMC1/NPC1/RPTOR/IFITM3/COL6A1/TRPM2/SNAPIN/ARL8A/ATP6V1C2/STK11IP/MARCHF1/CCZ1B/ATP6V1B2/ATP6V0D2/VLDLR/STOM/SURF4/LAMTOR1/TMEM138/SIDT2/ABCB9/NDUFC2/ATP6V1G3/AP1S3/DAB2/LPCAT1/ABCA5/ATP6V1C1/PIP4P2/PI4K2A/RRAGA/GNAQ/DRAM2/EEF1A1/SEC13/SLC30A2/CD1D/CD1B/ATP13A2/ATP6V0D1/SLC2A6/VPS11/MFSD12/AP2M1/CYB561A3/TPCN2/LAPTM5/NCSTN/PIGR/TMEM79/WDFY3/ABHD6/DTX3L/MFSD8/TEX264/WDR41/TLR3/DAGLB/MIOS/BRI3/TMEM74/HGSNAT/RNF183/MARCHF8/BORCS5/PIP4P1/HPS6/BORCS7/PLEKHF1/AP1G1/ANPEP/VPS39/YWHAB/UBXN6/RILP/WDR81/ATG16L2/VASN/LMBRD1/PLA2G4F/PCSK9/SPNS1/TMEM192/GABARAP/STX8/FPR1/GAA/CLCN5/C3AR1/GNB2/RAB37/HPSE/SLCO4C1/C12orf66/TMEM9B/RNF152/SLC38A9/ULK1/GBA/THBD/HLA-DQB1/ABCA13/TMEM150B/MYLPF/ZNRF2/SLC26A11/ATG9B/AP1S2/ANXA2/AP2A2/VPS33B/PRKD1/STING1/AP3M1/ATP6V0A2/TMEM179B/ANKFY1/ATP6V0C/LAMP1/TPCN1/PLA2G4E/LAMTOR4/CLN3/LITAF/BLOC1S2/HLA-DRB1/TMEM63A/GIMAP5/GDAP2/BORCS6/MYO6/TLR7/HLA-DQA1/AP2A1/ANXA6/IGF2R/ENTPD4/ENPP1/DPP4/PSAP/SZT2/SLC29A3/UVRAG/HLA-DRB5/GLMP/MTOR/HLA-DOA/HLA-DMA/HLA-DRA/NEU1/HLA-F/TECPR1/RAB12/SLC48A1/SLC35F6/ATP6V1G2/CPNE1/VPS16/UBA52/HLA-DPB1/PLEKHM1/HLA-DPA1/HLA-DQB2/SPAAR/NFAM1/HLA-DQA2/IRGM/RDH14/HLA-DOB/HLA-DMB/DDOST/TMEM199/TMEM150C/AP5B1/BORCS8/RAB44/MAP1LC3B2/RAB7B
           Count
GO:0030055   410
GO:0005924   406
GO:0005925   403
GO:0010008   468
GO:0005759   457
GO:0005774   403
```

```r
#write.csv(cc.go,file="cc.go.csv")
```

Perform a GO enrichment analysis
==========================================================
class: small-code

```r
bp.go <- enrichGO(gene          = deg_ann.tb$entrezgene_id,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(bp.go)
```

```
                   ID
GO:0022604 GO:0022604
GO:0007265 GO:0007265
GO:0043161 GO:0043161
GO:0010498 GO:0010498
GO:0006914 GO:0006914
GO:0061919 GO:0061919
                                                                 Description
GO:0022604                                  regulation of cell morphogenesis
GO:0007265                                   Ras protein signal transduction
GO:0043161 proteasome-mediated ubiquitin-dependent protein catabolic process
GO:0010498                             proteasomal protein catabolic process
GO:0006914                                                         autophagy
GO:0061919                            process utilizing autophagic mechanism
           GeneRatio   BgRatio       pvalue     p.adjust       qvalue
GO:0022604 476/15108 484/18670 6.545194e-34 4.230159e-30 2.215031e-30
GO:0007265 440/15108 448/18670 8.785327e-31 2.838978e-27 1.486570e-27
GO:0043161 413/15108 419/18670 1.522027e-30 3.278953e-27 1.716953e-27
GO:0010498 466/15108 477/18670 4.065504e-30 6.568837e-27 3.439630e-27
GO:0006914 482/15108 496/18670 7.697369e-29 8.291349e-26 4.341586e-26
GO:0061919 482/15108 496/18670 7.697369e-29 8.291349e-26 4.341586e-26
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              geneID
GO:0022604                                                                                                                                                                                  FGR/SEMA3F/HECW1/SARM1/PLXND1/ARHGAP33/DVL2/MYCBP2/TNFRSF12A/IFRD1/ARHGAP44/CDKL3/PAFAH1B1/CDKL5/SS18L2/SEMA3G/FYN/SEMA3B/RUFY3/ZRANB1/PLEKHO1/CD44/MYOC/RTN4R/CUL7/XK/NEDD4L/HEXB/FSTL4/CAMK2B/LZTS1/NGFR/NTN1/NGEF/MYO9A/RHOA/SYT1/NEDD4/CDC42/WDR1/LIMS2/TRPC5/MARK2/FERMT2/MAP3K13/SEMA3A/SEMA3C/DLG1/ARHGAP15/PLXNA2/PAK3/CAPZB/FBLN1/MAP2/CARMIL1/DNM2/RIMS1/RAB21/CHRNA3/PSEN1/EPB41L3/SEMA5B/GSK3B/ULK2/ZMPSTE24/CTTN/MT3/DNM1L/CASS4/LZTS3/ARHGAP4/ICAM1/MUL1/NRCAM/SEMA6A/EZR/UNC13D/DPYSL2/SEMA4G/BAMBI/ABL1/NRP1/PALMD/PALM/CRKL/TRIOBP/MYH9/COCH/NIN/CPNE6/NFATC4/ADNP/HCK/MCF2/FMR1/FGD1/ARHGEF7/OLFM4/CORO1A/METRN/FBXO31/EEF2K/SYT17/RAB11A/RHOV/ARHGEF18/ANKRD27/MYH14/PTPRS/MEGF8/MAG/HPN/PTN/CCL24/LIMK1/KANK1/DVL1/DNMBP/CXCL12/WNT3/C1QBP/CCL2/RND2/EFNB3/SLC9A3R1/SH3D19/RAPGEF2/MDK/CORO1C/CAPRIN2/CUX2/RASAL1/NEDD9/UST/BVES/NR2E1/SRF/VEGFA/SEMA5A/SPARC/DBN1/TBCCD1/UBE3A/WNT5A/PLXNA1/CSPG5/GORASP1/TTL/TLX2/RTN4/REG1A/FN1/STK25/RND3/TRAK2/EPHA4/TNR/RHOU/GBP1/MPL/STK11/APOA1/SPP1/ITGB1BP1/RHOQ/BCL11A/GNA13/PTK2B/ZMYM2/CXCR4/SEPTIN7/NEUROG3/LRP1/OBSL1/PREX1/VAMP7/CPNE5/RREB1/RAP2A/ID1/RHOJ/OMG/FGD3/MACF1/VIL1/RNF6/KDR/YWHAH/CDC42EP1/STRIP2/MKLN1/CHN1/FGF13/APOE/ACTN4/ZSWIM6/UNC13A/OLFM1/CAMSAP1/PLXNA3/DIAPH1/MAP1B/ZSWIM4/ENAM/DLG4/SYT4/ZMYM5/WASF3/MYH10/SPART/POSTN/EPHB2/NGF/LRP4/AGO4/WDR36/ANXA1/ADGRB3/EPHA7/CAPRIN1/ITGA7/SEMA4F/PLXNC1/RAC1/DBNL/ERMN/SKIL/ENPP2/PLAA/TUBB2B/TRPC6/RDX/PARP6/ITPKA/ADAM10/SEMA6D/ACTR2/ANXA7/HECW2/ABI2/SEMA7A/AP1AR/SHROOM3/PPP3CA/FGD4/RHOF/TPM1/NTRK3/ARHGDIA/RNF157/WTIP/BARHL2/STRIP1/CRABP2/SEMA6C/SYT2/RHOB/WDPCP/CPNE9/GOLGA4/LIMD1/SLIT2/MYO10/RASA1/FGD2/PHIP/ARHGAP18/ZMYM4/GNA12/SH3KBP1/MSN/ZMYM3/DOCK5/NTRK2/FAM171A1/GAS2/LRRC4C/PAK1/YTHDF1/CDC42EP2/ALDOA/TAOK2/DOCK1/PTPRO/EPS8/FEZF2/PTPRD/CFDP1/SEMA3D/THY1/TNIK/WNT3A/CHODL/WNT7A/FGD5/ZFYVE27/MOV10/RHOC/GRIP1/TIAM1/LRP8/KIT/PALM2AKAP2/FMNL2/ARMCX5-GPRASP2/PPP1R9A/DMTN/CDC42SE2/EPB41/ARHGAP35/KALRN/CFAP410/ITGB2/S100B/BRSK1/CHRNB2/DVL3/FMNL3/LARP4/ZSWIM5/FBLIM1/DRAXIN/PDPN/DISC1/PDLIM5/CDC42EP3/TRIM46/ARPC2/NKX6-1/RYK/ZMYM6/MELTF/PLXNB1/RHOBTB3/DLC1/CDK5/BRWD3/PDZD8/RET/NSMF/ILK/NDEL1/ISLR2/CRK/CDC42EP5/DAPK3/SEMA6B/RHOH/SEMA4C/SHOX2/MAP2K1/EFNA1/IL1RAPL1/SIN3A/PTK2/SDC2/ZEB2/RAC3/LIMS1/NLGN1/LINGO1/ROBO1/P2RY1/SEMA3E/CDH2/SLC26A5/HAS2/LPAR3/MAP6/FGG/FGA/FGB/DSCAM/KNDC1/PTEN/CCL11/RND1/WASHC2C/CORO1B/CFL1/RHOD/DAG1/DAB1/BRSK2/FZD4/DHX36/FBXW8/BAIAP2/ANAPC2/ZNF135/RIMS2/SYNE3/BDNF/CDK5R1/GRIN1/RHOG/SLITRK1/RCC2/CALR/CDH4/CDC42EP4/MYADM/F2/FGD6/CCL13/AMIGO1/IST1/FES/CSF1R/EPHB3/TRAK1/CNTN2/TACSTD2/EFNA5/SS18L1/POU3F2/FMNL1/ROBO2/SEMA4B/CIB1/PRKN/BRWD1/BCL9L/ARAP1/MAPT/SHTN1/DCC/TRPV2/SEMA4D/C15orf62/PARVB/OSTN/LRRK2/RELN/S100A13/SEMA4A/NLGN3/PRPF40A/PLXNB2/FLNA/ZMYM1/SRC/SYNGAP1/FITM2/SIPA1L1/CDC42SE1/PARVA/S100A10/KIF13B/DNM3/KEL/LPAR1/QRICH1/ITSN2/ARC/SMURF1/PLXNB3/OPA1/BHLHB9/L1CAM/GDI1/BMPR2/LST1/ATP10A/SYT3/ZSWIM8/PLXNA4/ARPIN/TWF2/SHANK3/CUX1/CYFIP1/PRAG1/SRCIN1/NEFL/CCL3
GO:0007265                                                                                                                                                                                                                                                                                                       ALS2/ITGA3/RALA/FARP2/ARHGAP44/CYTH3/PLEKHG6/RHOBTB2/NISCH/WAS/IGF1/MYOC/TIMP2/CUL3/TRIO/RIPOR1/RTN4R/RAB27B/PREX2/ARHGAP6/USP28/KITLG/ARHGEF5/MCF2L2/RASGRF1/PSD/NCKAP1/NGFR/RFXANK/NTN1/NGEF/ARFGEF1/RHOA/ROCK1/RASSF1/RASGRP2/RAB27A/FGF10/CDC42/MAP4K4/RHOBTB1/HACD3/ARHGEF10L/CELSR1/PLD1/RAB7A/ARHGEF1/GNB1/P2RY10/RAPGEF3/OPHN1/DNM2/RABL2B/RAB21/RAB10/HACE1/MAPKAPK5/BRAP/ARHGAP4/PLEKHG2/TGFB2/JAK2/ABL1/RAB18/NRP1/MYO9B/CRKL/LZTR1/CYTH4/RAB36/GRAP2/SGSM3/IFT27/SOS2/ARHGAP5/ARFGAP1/MCF2/FGD1/ARHGEF7/RAB11A/RHOV/RAB2A/PPP2CB/ARHGEF10/ARHGEF18/RASAL3/CYTH2/RAB3D/RASIP1/RAB3A/CADM4/RASA4/MET/LIMK1/KANK1/RAPGEF1/DNMBP/GBF1/SHOC2/KPNB1/CYTH1/RAB5C/RND2/RANGRF/TNFAIP1/RAB34/RAPGEF2/ARHGEF17/APOC3/CBL/MADD/KCTD10/FOXM1/RASAL1/ARHGDIB/RAB5B/RAB35/RIPOR2/MAPK14/RAB23/WASF1/ARFGEF3/RASGRF2/PDGFRB/BCL6/ECT2/ROPN1B/ARHGEF26/RTKN/DNAJC27/DOK1/SOS1/RND3/DHCR24/RALGPS2/PARK7/RAP1A/RHOU/ARHGEF2/MFN2/SSX2IP/RAB29/STMN1/APOA1/KIF14/RAB32/CTNNAL1/RAB14/MAPKAP1/DENND1A/RHOQ/GNA13/PLEKHG1/IQSEC3/ADRA1A/GPSM2/CDK2/RAB9B/RAB9A/RAP2C/RAB38/PREX1/ARFGEF2/RAB22A/IQSEC2/STAMBP/CDKN1A/RREB1/RAB17/GPR18/RAP2A/PSD4/MCF2L/RRAS/RHOJ/PLEKHG3/RHOT1/FGD3/RAP1B/F2RL3/CDC42EP1/RAC2/IFT22/RAB2B/ARHGEF6/SHC2/APOE/EPO/TRIM28/ARHGEF16/DNMT1/EPS8L1/ARHGEF9/ZNF304/RAF1/ARHGEF11/RAB25/EPHB2/KRAS/RRAS2/VAV3/NOTCH2/NGF/ROCK2/RERG/RAB33A/GPR55/ARHGEF4/SPRY2/RAC1/DBNL/MMD2/RALGPS1/DAB2IP/ARHGEF39/FLOT1/ARRB1/RAB30/RDX/ARHGAP29/RAB1A/PLCE1/ABI2/USP8/FGF2/G3BP2/FGD4/LPAR6/RB1/RHOF/RAB20/RAB15/GPR65/CDH13/RHOT2/KSR1/TP53/ARHGDIA/RAB40B/VAV1/ARHGEF19/CNKSR1/ABL2/SETDB1/RAB13/RIT1/ITPKB/RHOB/RALB/RABL2A/GPR17/RAB5A/IQSEC1/RABL3/AGTR1/CCNA2/PLK2/RASA1/G3BP1/PSD2/DOK3/FGD2/DYNLT1/TIAM2/GNA12/RAB19/RAB41/GPR174/OGT/DOK2/CDKN2A/SHC3/NOTCH1/CDC42EP2/ITGB1/ADRA2A/EPS8/RIT2/RASGRP3/FARP1/RAB3C/PLEKHG4B/GRAP/OBSCN/RABGEF1/FGD5/FLCN/RAB6B/RHOC/ELMO1/RASA2/RAB39B/PSD3/TIAM1/NRG1/DGKI/RAB28/MRAS/WASF2/AUTS2/CDC42SE2/RAPGEF6/ABR/ARHGAP35/KALRN/RALGDS/VAV2/SHC1/SH2B2/SQSTM1/PDPN/VANGL2/CDC42EP3/ARHGEF3/SPRY1/FBXO8/F2RL2/F2RL1/RHOBTB3/COL1A2/DLC1/ABCA1/FBP1/ARF6/ARHGEF40/ARHGAP42/RAB8B/SYNPO2L/MAPRE2/CRK/RAB8A/RAB4B/CDC42EP5/RAB26/RAB4A/RHOH/RAB31/COL3A1/SPRY3/RAB3B/RAB24/RAC3/ROBO1/CAVIN4/AKAP13/KSR2/PLEKHG5/RASGRP4/RAB33B/RAB40A/RASGRP1/RND1/CFL1/RAB43/RAB37/RHOD/SCAI/HEG1/NET1/TNK1/ABRA/HRAS/RAB1B/KCTD13/ARHGAP1/RAB6A/RHOG/EPS8L2/GPR4/JUN/GRB2/GPR35/RAB39A/CDC42EP4/FGD6/F2R/RAP2B/P2RY8/SPATA13/ARHGEF37/CCDC125/IQGAP3/PRKD1/CSF1/RAB11B/MAPK11/IRS2/RASA3/BCR/SHTN1/EPOR/PLEKHG7/SPRY4/ERAS/C15orf62/PLEKHG4/NF1/ARHGEF12/KANK2/SYNGAP1/STMN3/RAB40C/CDC42SE1/LPAR1/ITSN2/NTRK1/EPS8L3/DENND4B/ARHGEF15/RASGEF1A/ECT2L/GDI1/PHACTR4/STK19/AIF1/GPR20/ADGRG1/ITSN1/RAB12/NUP62/NRAS/CHUK/SYNJ2BP/LAT/TAX1BP3/ARHGEF33/ARHGEF28/RAB6C/RAB6D/ARHGEF38/RGL2/ARHGEF25/ARHGDIG/LYN/BRK1/RAB44/PECAM1/CYFIP1/PRAG1/RAB7B
GO:0043161                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          GCLC/ANKIB1/HECW1/CDC27/FBXL3/PSMB1/HFE/RNF216/PSMC4/NUB1/RNF14/UFL1/WWTR1/UBR2/PIAS1/RNF19A/CUL3/CLEC16A/PSMA4/HSPA5/ERCC8/NEDD4L/JKAMP/MAPK9/HERPUD1/ANAPC4/CUL1/RNF4/ATXN3/PSME4/SIRT2/ERLEC1/NEDD4/UFD1/RNF126/SEL1L/TRIB2/UBE2D1/CHFR/FBXW11/DERL2/TBX21/KLHL20/DNAJC10/UBE2A/ARAF/UBE2K/PPP2R5C/ITCH/SENP1/KEAP1/PSEN1/PCNP/GSK3B/PSMC5/KLHL42/AURKA/EDEM2/MMP24-AS1-EDEM2/NSFL1C/ANAPC5/MAEA/PSME1/RFFL/CDC23/PSMD5/WAC/HSP90AB1/SIRT1/PSMD8/FBXL19/CDC34/KLHL22/DERL3/FOXRED2/KCTD17/RBX1/TRIM9/PSMC6/PSMA3/PSMC1/PSMB5/DCAF11/PSMA6/PSME2/PSMA7/SEL1L2/TRIB3/USP14/SMAD7/PSMD10/TBL1X/N4BP1/PSMD7/FBXO31/STUB1/BFAR/UBE2W/UBXN8/PPP2CB/SGTA/FZR1/DMAC2/NOP53/GSK3A/CAV1/PSMA2/SEC61B/UBE2R2/DVL1/ERLIN1/FBXW4/FBXL15/CUL2/UBE2S/FBXL20/PSMD3/CCDC47/PSMD11/SMURF2/TNFAIP1/UBE2D3/WFS1/TRIM2/FBXW7/TRIM3/ANAPC15/BIRC2/UBE4A/PSMD9/KCTD10/USP5/SPSB2/FBXO5/FBXO9/FBXL4/TRIM38/FAF2/SKP1/CSNK1A1/CRBN/UBE3A/ARMC8/RNF7/PSMD14/AUP1/SUMO1/NFE2L2/PARK7/EDEM3/RNF19B/FBXO2/FBXO6/UCHL5/CDC20/KIF14/UBE3D/FBXL5/UBE2B/ECRG4/RAD23B/AREL1/CLU/TRIM25/SPOP/WWP1/GIPC1/CDK2/LRRC29/PSMF1/RBCK1/PSMB2/HECTD3/FBXL12/FBXL16/RFPL1/DNAJB9/HERC2/PSMA1/CBFA2T3/CDC16/DDA1/UBE4B/PSME3/SELENOS/UBE2G1/FBXO44/HSPBP1/AGAP3/ARNTL/RNF122/CCNB1/EDEM1/CLOCK/UBQLN1/OS9/MDM2/FBXL8/DNAJB2/USP44/ECPAS/PSMB7/DERL1/IL33/PLAA/FOXF2/RNF144B/ARRB1/RNF121/SDCBP/HECW2/BBS7/RNF185/PRICKLE1/CUL4A/PML/DET1/GID4/ARRB2/CSNK1D/ANAPC11/AKT1/PSMB6/RPL11/PSMA5/COP1/RPS27A/UBXN4/SPOPL/UBR3/RNF175/MARCHF6/NKD2/SKP2/PLK2/FBXL17/FBXO38/RMND5B/RNF217/GNA12/TLK2/TAF1/OGT/ERLIN2/DNAJB12/UBC/RNF144A/AMN1/FBXO4/DAB2/ANAPC1/FBXL2/RMND5A/UCHL1/BUB3/PSMA8/FBXL18/HSPA13/SH3RF2/BUB1B/KLHL40/UBXN11/CUL4B/ZFAND2B/CCAR2/FBXW5/PSMD4/PSMB4/UBR1/AMFR/UBE2J2/UBQLN4/FBXL13/PSMC2/FBXO27/SPSB3/CCNF/UBXN1/SYVN1/ANKZF1/STT3B/RYBP/PSMD6/RCHY1/DNAJB14/MAD2L1/ANAPC10/RNF180/PTTG1/CAMLG/SHH/TMUB1/TMEM67/VCP/FBXO33/NEMF/PSMC3/IFI27/BTRC/BAG5/ARIH1/CUL5/ANAPC16/HSP90B1/PLK1/MAP1A/FBXO22/KCTD5/DDB1/CTNNB1/PBK/RNF187/TMUB2/TMEM129/WNT10B/UBE2E1/CDK1/UBB/DNAJC18/RNF34/MTM1/SOCS5/SPSB1/FBXL14/USP19/KAT5/COMMD1/TRIB1/SMARCC1/PSMD1/PLK3/UBXN2A/KLHL15/FBXO45/GLMN/KCTD13/UBE2C/SPSB4/PSMD2/DDIT3/ANAPC2/CDC26/TRIM72/MAN1B1/FBXO39/ARIH2/TBL1XR1/GBA/ZFAND2A/RNF186/AURKB/RAD23A/SHARPIN/SOCS4/YOD1/KCTD2/SIAH2/FBXL6/NPLOC4/MTA1/FBXL7/UBE2G2/FAF1/PRKN/PSMD13/BCAP31/ZNRF1/UBE2H/USP7/NHLRC1/FAM122A/UBQLN2/NCCRP1/SUMO2/NHLRC3/LRRK2/FHIT/SIAH1/ANAPC7/PCBP2/PSMD12/PELI1/FBXL22/TOPORS/SGTB/CD2AP/SVIP/DDRGK1/STYX/WWP2/SMURF1/UBE2J1/LTN1/PSMB8/RNF5/HSPA1B/HSPA1A/BAG6/RACK1/FBXO48/TRIM13/PSMB10/CSNK1E/UBXN2B/SIAH3/UBA52/RNF103/PSMB9/CEBPA/FBXO17/TAF9/PANO1/PSMB3
GO:0010498                                                                                                                                                                                    GCLC/ANKIB1/HECW1/CDC27/FBXL3/OSBPL7/RHBDF1/PSMB1/PKD1/POMT2/HFE/RNF216/PSMC4/NUB1/RNF14/UFL1/WWTR1/UBR2/PIAS1/RNF19A/GABARAPL2/CUL3/CLEC16A/PSMA4/HSPA5/ERCC8/NEDD4L/JKAMP/MAPK9/HERPUD1/ANAPC4/CUL1/USP13/RNF4/ATXN3/PSME4/SIRT2/ERLEC1/NEDD4/UFD1/RNF126/SEL1L/TRIB2/PRKACA/UBE2D1/CHFR/FBXW11/DERL2/TBX21/KLHL20/DNAJC10/UBE2A/ARAF/LAMP3/UBE2K/PPP2R5C/ITCH/SENP1/OPHN1/KEAP1/PSEN1/PCNP/NFE2L1/GSK3B/PSMC5/KLHL42/AURKA/EDEM2/MMP24-AS1-EDEM2/NSFL1C/ANAPC5/MAEA/PSME1/RFFL/CDC23/PSMD5/WAC/TREM2/HSP90AB1/SIRT1/PSMD8/FBXL19/CDC34/KLHL22/DERL3/FOXRED2/KCTD17/RBX1/TRIM9/PSMC6/PSMA3/PSMC1/PSMB5/DCAF11/PSMA6/PSME2/PSMA7/SEL1L2/TRIB3/HM13/USP14/SMAD7/PSMD10/TBL1X/FMR1/UGGT2/N4BP1/PSMD7/FBXO31/STUB1/BFAR/RNF40/UBE2W/UBXN8/PPP2CB/SGTA/FZR1/DMAC2/NOP53/GSK3A/CAV1/PSMA2/SEC61B/UBE2R2/DVL1/ERLIN1/ATE1/FBXW4/FBXL15/CUL2/UBE2S/FBXL20/PSMD3/CCDC47/PSMD11/SMURF2/TNFAIP1/UBE2D3/WFS1/TRIM2/FBXW7/PRPF19/TRIM3/ANAPC15/BIRC2/UBE4A/PSMD9/KCTD10/USP5/SPSB2/FBXO5/FBXO9/BAG2/FBXL4/TRIM38/FAF2/SKP1/CSNK1A1/CRBN/UBE3A/ARMC8/RNF7/PSMD14/AUP1/SUMO1/NFE2L2/PARK7/EDEM3/RNF19B/FBXO2/FBXO6/UCHL5/CDC20/KIF14/UBE3D/FBXL5/UBE2B/ECRG4/RAD23B/AREL1/CLU/TRIM25/SPOP/WWP1/GIPC1/CDK2/LRRC29/PSMF1/RBCK1/PSMB2/HECTD3/NR1D1/PRKCG/FBXL12/FBXL16/SDF2L1/RFPL1/DNAJB9/HERC2/PSMA1/CBFA2T3/CDC16/APOE/DDA1/UBE4B/PSME3/SELENOS/UBE2G1/SDF2/FBXO44/SCO1/HSPBP1/AGAP3/ARNTL/RNF122/CCNB1/EDEM1/CLOCK/UBAC2/UBQLN1/OS9/MDM2/FBXL8/DNAJB2/USP44/NUDT15/MARCHF7/UGGT1/ECPAS/TOR1A/DAB2IP/PSMB7/DERL1/IL33/PLAA/FOXF2/RNF144B/ARRB1/RNF121/SDCBP/HECW2/BBS7/RNF185/PRICKLE1/TMTC3/CUL4A/PML/DET1/GID4/ARRB2/CSNK1D/ANAPC11/PMAIP1/AKT1/PSMB6/RCN3/RPL11/PSMA5/COP1/RPS27A/UBXN4/SPOPL/UBR3/RHBDD1/TMF1/RNF175/MARCHF6/NKD2/SKP2/PLK2/FBXL17/FBXO38/RMND5B/RNF217/GNA12/TLK2/TAF1/OGT/ERLIN2/ALAD/DNAJB12/UBC/RNF144A/AMN1/FBXO4/DAB2/ANAPC1/FBXL2/RMND5A/UCHL1/BUB3/PSMA8/FBXL18/HSPA13/USP25/SH3RF2/BUB1B/KLHL40/UBXN11/CUL4B/ZFAND2B/PINK1/CCAR2/FBXW5/PSMD4/PSMB4/UBR1/AMFR/UBE2J2/UBQLN4/FBXL13/PSMC2/FBXO27/SPSB3/CCNF/UBXN1/SYVN1/ANKZF1/STT3B/RYBP/PSMD6/RCHY1/DNAJB14/MAD2L1/ANAPC10/RNF180/PTTG1/CAMLG/SHH/TMUB1/TMEM67/VCP/FBXO33/NEMF/PSMC3/IFI27/BTRC/BAG5/ARIH1/CUL5/ANAPC16/HSP90B1/PLK1/MAP1A/FBXO22/UBXN6/KCTD5/DDB1/CTNNB1/PBK/RNF187/TMUB2/TMEM129/WNT10B/UBE2E1/CDK1/UBB/DNAJC18/RNF34/SOCS6/RNF139/MTM1/SOCS5/ENC1/SPSB1/FBXL14/USP19/PSME3IP1/KAT5/COMMD1/TRIB1/SMARCC1/PSMD1/PLK3/UBXN2A/KLHL15/FBXO45/BRSK2/GLMN/KCTD13/UBE2C/SPSB4/PSMD2/DDIT3/NUPR1/ANAPC2/CDC26/TRIM72/MAN1B1/FBXO39/ARIH2/TBL1XR1/GBA/ZFAND2A/RNF186/AURKB/RAD23A/SHARPIN/SOCS4/YOD1/KCTD2/SIAH2/RNF41/TMEM259/FBXL6/NPLOC4/MTA1/FBXL7/UBE2G2/FAF1/PRKN/PSMD13/BCAP31/ZNRF1/UBE2H/USP7/NHLRC1/FAM122A/UBQLN2/NCCRP1/SUMO2/NHLRC3/LRRK2/RNFT1/FHIT/SIAH1/ANAPC7/PCBP2/PSMD12/DDI2/PELI1/FBXL22/TOPORS/SGTB/CD2AP/SVIP/DDRGK1/STYX/WWP2/SMURF1/UBE2J1/LTN1/PSMB8/RNF5/HSPA1B/HSPA1A/BAG6/RACK1/FBXO48/TRIM13/PSMB10/CSNK1E/UBXN2B/SIAH3/UBA52/GPX1/RNF103/GET4/PSMB9/CEBPA/ECSCR/DNAAF4/FBXO17/TAF9/PANO1/PSMB3
GO:0006914 BAD/POLDIP2/LAMP2/OSBPL7/VPS41/VTA1/DCN/CAPN1/HGF/RB1CC1/RRAGD/IFT88/ATP6V0A1/GABARAPL2/DAPK2/CLEC16A/PHF23/ATP6V1H/VPS13D/PIK3CB/TAB2/USP36/ATG5/USP13/LZTS1/VMP1/CTSA/GNAI3/KDM4A/CD84/ATG2B/ROCK1/TBC1D25/SIRT2/VPS35/FUNDC1/NEDD4/WIPI1/CSNK2A2/PRKACA/SREBF1/MARK2/TSG101/RAB7A/USP33/LAMP3/PIK3C3/TIGAR/TP53INP2/TOLLIP/SESN1/MID2/PSEN1/HSP90AA1/GSK3B/ULK2/RRAGB/CHMP2B/ZMPSTE24/CTTN/ATG16L1/MT3/DNM1L/NSFL1C/DYNLL1/SNX5/SPTLC1/MUL1/MCOLN1/EXOC1/SCFD1/HDAC6/WAC/TREM2/SIRT1/ABL1/SH3GLB1/KLHL22/SNAP29/TOMM22/FBXO7/HMOX1/TSPO/MTMR3/EP300/HDAC10/ATP6V1D/SPTLC2/HIF1A/TM9SF1/EEF1A2/TRIB3/CSNK2A1/CHMP4B/MAP1LC3A/ATG4A/SRPX/MTMR8/PIM2/RUBCNL/SUPT20H/MAPK3/NPRL3/USP10/TSC2/STUB1/MEFV/PYCARD/KAT8/HERC1/BMF/VPS18/RIPK2/MTMR9/BNIP3L/SNRNP70/TBC1D17/CDC37/RASIP1/PIK3R2/GSK3A/NAMPT/MET/NOD1/HSPB1/RHEB/PRKAG2/MAPK8/ACBD5/IFT20/LAMTOR3/FBXW7/PPARGC1A/HSPA8/ATG2A/EIF4G2/IL10RA/AMBRA1/TCIRG1/CAMKK2/IFNG/GAPDH/PRKAB1/RAB23/EPM2A/TFEB/DAP/ARSB/ATP6V0E1/NPRL2/ATP6V1A/EIF4G1/PIKFYVE/STAM2/HTRA2/EVA1A/CHMP3/PRKAG3/RAB3GAP1/ATP6V1B1/TMEM59/QSOX1/PARK7/KIAA1324/LAMTOR2/LEPR/MFN2/EXOC8/RRAGC/LGALS8/ATP6V0B/STX12/STK11/VAMP8/FOXO3/STBD1/RAB3GAP2/VPS4B/CLU/ADRA1A/PIK3CA/VPS26A/OPTN/ATG101/KIF25/ATG4C/CAPNS1/BECN1/ATG14/WDR24/TICAM1/EMC6/KDR/EIF2AK4/VPS13C/SNX6/SH3BP4/TOMM40/MAP1S/CHMP2A/ATG4D/SESN2/ATP5IF1/GFAP/ATP6V1E1/TBC1D5/VPS25/EXOC4/PRKAB2/C19orf12/TRIM21/TRIM5/TRIM22/PRKAA1/TBC1D14/ITGB4/VPS4A/PIK3C2B/LARS1/SBF2/PTPN22/LAMTOR5/UBQLN1/USP30/SNX14/MAP3K7/HTR2B/DRAM1/VPS36/CALCOCO2/GATA4/IL10/STAM/STX17/ATP6V1G1/XPA/PLAA/FNBP1L/RAB1A/ANXA7/RNF185/GABARAPL1/TMBIM6/VPS33A/VPS37B/ULK3/MAP1LC3B/IRF8/CLTC/RMC1/NPC1/TP53/RPTOR/FOXK2/TRIM65/WDR45B/MVB12A/AKT1/ABL2/MCL1/CTSK/S100A8/SNAPIN/HAX1/ATP6V1C2/RALB/RAB5A/ATG3/RUBCN/SNCA/CISD2/PLK2/LIX1/ATG12/KLHL3/TLK2/ATP6V1B2/ATP6V0D2/MTDH/C9orf72/LRSAM1/NRBF2/SESN3/DRD2/ATM/LAMTOR1/EI24/FEZ1/SOGA1/VPS51/PIP4K2A/FOXO1/ITPR1/VPS26B/VTI1A/BAG3/HSPB8/EPG5/ATG10/BCL2L11/SCOC/TRAPPC8/FBXL2/RETREG1/TOMM70/UCHL1/FLCN/ATP6V1C1/LARP1/RRAGA/RAB39B/VPS37A/DRAM2/MTERF3/EEF1A1/TAB3/WIPI2/PINK1/SNF8/ATP13A2/HK2/ATP6V0D1/GPSM1/VPS11/UBQLN4/VPS28/SQSTM1/TPCN2/PRKAA2/TRIM17/SPATA18/S100A9/FZD5/DAPL1/IFI16/WDFY3/MTMR14/FYCO1/ZC3H12A/MFSD8/TEX264/WDR41/CASP3/CHMP4C/TMEM74/CDK5/FOXK1/TP53INP1/VCP/DEPP1/TSC1/FUNDC2/ZFYVE1/PLEKHF1/TMEM41B/VPS39/PIP4K2C/GOLGA2/RAB8A/DAPK3/UBXN6/WDR81/MLST8/VPS37C/ATG16L2/PAFAH1B2/ATG4B/MTCL1/STAT3/TMEM208/TMEM150A/DHRSX/RAB24/ADRB2/PTK2/GABARAP/FEZ2/MTM1/MFN1/ATP6V0E2/TRIM8/RGS19/BCL2/RAB33B/SYNPO2/SNX32/KAT5/TOMM20/HAP1/PLK3/UBXN2A/LEP/NLRP6/RAB1B/DDIT3/ATG13/ERCC4/TOMM5/NUPR1/CHMP6/BNIP3/VPS37D/RNF152/BOK/CDK5R1/SMCR8/SLC38A9/ACER2/ULK1/GBA/WDR6/ERN1/RAB39A/PACS2/TMEM150B/YOD1/MAPK15/WASHC1/ATG9B/RNF41/PRKAG1/EXOC7/TBK1/PRKD1/STING1/NRBP2/ATP6V0A2/PRKN/HGS/ATP6V0C/TPCN1/MAPT/NHLRC1/UBQLN2/LAMTOR4/RUFY4/NBR1/CLN3/LRRK2/ZKSCAN3/HMGB1/SUPT5H/PIK3R4/TECPR2/TOMM7/DAPK1/ILRUN/WDR45/SRC/HTT/ATG7/PSAP/VPS13A/SVIP/UVRAG/SMURF1/MTOR/SREBF2/ATG9A/RNF5/CSNK2B/TRIM13/TECPR1/RAB12/VDAC1/ATP6V1G2/FIS1/UBXN2B/VPS16/CPTP/PLEKHM1/IRGM/PGAM5/TMEM150C/ATP6V1E2/CHMP4A/MAP1LC3B2/DYNLL2/SEC22B/IKBKG/LIX1L/PIP4K2B
GO:0061919 BAD/POLDIP2/LAMP2/OSBPL7/VPS41/VTA1/DCN/CAPN1/HGF/RB1CC1/RRAGD/IFT88/ATP6V0A1/GABARAPL2/DAPK2/CLEC16A/PHF23/ATP6V1H/VPS13D/PIK3CB/TAB2/USP36/ATG5/USP13/LZTS1/VMP1/CTSA/GNAI3/KDM4A/CD84/ATG2B/ROCK1/TBC1D25/SIRT2/VPS35/FUNDC1/NEDD4/WIPI1/CSNK2A2/PRKACA/SREBF1/MARK2/TSG101/RAB7A/USP33/LAMP3/PIK3C3/TIGAR/TP53INP2/TOLLIP/SESN1/MID2/PSEN1/HSP90AA1/GSK3B/ULK2/RRAGB/CHMP2B/ZMPSTE24/CTTN/ATG16L1/MT3/DNM1L/NSFL1C/DYNLL1/SNX5/SPTLC1/MUL1/MCOLN1/EXOC1/SCFD1/HDAC6/WAC/TREM2/SIRT1/ABL1/SH3GLB1/KLHL22/SNAP29/TOMM22/FBXO7/HMOX1/TSPO/MTMR3/EP300/HDAC10/ATP6V1D/SPTLC2/HIF1A/TM9SF1/EEF1A2/TRIB3/CSNK2A1/CHMP4B/MAP1LC3A/ATG4A/SRPX/MTMR8/PIM2/RUBCNL/SUPT20H/MAPK3/NPRL3/USP10/TSC2/STUB1/MEFV/PYCARD/KAT8/HERC1/BMF/VPS18/RIPK2/MTMR9/BNIP3L/SNRNP70/TBC1D17/CDC37/RASIP1/PIK3R2/GSK3A/NAMPT/MET/NOD1/HSPB1/RHEB/PRKAG2/MAPK8/ACBD5/IFT20/LAMTOR3/FBXW7/PPARGC1A/HSPA8/ATG2A/EIF4G2/IL10RA/AMBRA1/TCIRG1/CAMKK2/IFNG/GAPDH/PRKAB1/RAB23/EPM2A/TFEB/DAP/ARSB/ATP6V0E1/NPRL2/ATP6V1A/EIF4G1/PIKFYVE/STAM2/HTRA2/EVA1A/CHMP3/PRKAG3/RAB3GAP1/ATP6V1B1/TMEM59/QSOX1/PARK7/KIAA1324/LAMTOR2/LEPR/MFN2/EXOC8/RRAGC/LGALS8/ATP6V0B/STX12/STK11/VAMP8/FOXO3/STBD1/RAB3GAP2/VPS4B/CLU/ADRA1A/PIK3CA/VPS26A/OPTN/ATG101/KIF25/ATG4C/CAPNS1/BECN1/ATG14/WDR24/TICAM1/EMC6/KDR/EIF2AK4/VPS13C/SNX6/SH3BP4/TOMM40/MAP1S/CHMP2A/ATG4D/SESN2/ATP5IF1/GFAP/ATP6V1E1/TBC1D5/VPS25/EXOC4/PRKAB2/C19orf12/TRIM21/TRIM5/TRIM22/PRKAA1/TBC1D14/ITGB4/VPS4A/PIK3C2B/LARS1/SBF2/PTPN22/LAMTOR5/UBQLN1/USP30/SNX14/MAP3K7/HTR2B/DRAM1/VPS36/CALCOCO2/GATA4/IL10/STAM/STX17/ATP6V1G1/XPA/PLAA/FNBP1L/RAB1A/ANXA7/RNF185/GABARAPL1/TMBIM6/VPS33A/VPS37B/ULK3/MAP1LC3B/IRF8/CLTC/RMC1/NPC1/TP53/RPTOR/FOXK2/TRIM65/WDR45B/MVB12A/AKT1/ABL2/MCL1/CTSK/S100A8/SNAPIN/HAX1/ATP6V1C2/RALB/RAB5A/ATG3/RUBCN/SNCA/CISD2/PLK2/LIX1/ATG12/KLHL3/TLK2/ATP6V1B2/ATP6V0D2/MTDH/C9orf72/LRSAM1/NRBF2/SESN3/DRD2/ATM/LAMTOR1/EI24/FEZ1/SOGA1/VPS51/PIP4K2A/FOXO1/ITPR1/VPS26B/VTI1A/BAG3/HSPB8/EPG5/ATG10/BCL2L11/SCOC/TRAPPC8/FBXL2/RETREG1/TOMM70/UCHL1/FLCN/ATP6V1C1/LARP1/RRAGA/RAB39B/VPS37A/DRAM2/MTERF3/EEF1A1/TAB3/WIPI2/PINK1/SNF8/ATP13A2/HK2/ATP6V0D1/GPSM1/VPS11/UBQLN4/VPS28/SQSTM1/TPCN2/PRKAA2/TRIM17/SPATA18/S100A9/FZD5/DAPL1/IFI16/WDFY3/MTMR14/FYCO1/ZC3H12A/MFSD8/TEX264/WDR41/CASP3/CHMP4C/TMEM74/CDK5/FOXK1/TP53INP1/VCP/DEPP1/TSC1/FUNDC2/ZFYVE1/PLEKHF1/TMEM41B/VPS39/PIP4K2C/GOLGA2/RAB8A/DAPK3/UBXN6/WDR81/MLST8/VPS37C/ATG16L2/PAFAH1B2/ATG4B/MTCL1/STAT3/TMEM208/TMEM150A/DHRSX/RAB24/ADRB2/PTK2/GABARAP/FEZ2/MTM1/MFN1/ATP6V0E2/TRIM8/RGS19/BCL2/RAB33B/SYNPO2/SNX32/KAT5/TOMM20/HAP1/PLK3/UBXN2A/LEP/NLRP6/RAB1B/DDIT3/ATG13/ERCC4/TOMM5/NUPR1/CHMP6/BNIP3/VPS37D/RNF152/BOK/CDK5R1/SMCR8/SLC38A9/ACER2/ULK1/GBA/WDR6/ERN1/RAB39A/PACS2/TMEM150B/YOD1/MAPK15/WASHC1/ATG9B/RNF41/PRKAG1/EXOC7/TBK1/PRKD1/STING1/NRBP2/ATP6V0A2/PRKN/HGS/ATP6V0C/TPCN1/MAPT/NHLRC1/UBQLN2/LAMTOR4/RUFY4/NBR1/CLN3/LRRK2/ZKSCAN3/HMGB1/SUPT5H/PIK3R4/TECPR2/TOMM7/DAPK1/ILRUN/WDR45/SRC/HTT/ATG7/PSAP/VPS13A/SVIP/UVRAG/SMURF1/MTOR/SREBF2/ATG9A/RNF5/CSNK2B/TRIM13/TECPR1/RAB12/VDAC1/ATP6V1G2/FIS1/UBXN2B/VPS16/CPTP/PLEKHM1/IRGM/PGAM5/TMEM150C/ATP6V1E2/CHMP4A/MAP1LC3B2/DYNLL2/SEC22B/IKBKG/LIX1L/PIP4K2B
           Count
GO:0022604   476
GO:0007265   440
GO:0043161   413
GO:0010498   466
GO:0006914   482
GO:0061919   482
```

```r
#write.csv(bp.go,file="bp.go.csv")
```

Perform a GO enrichment analysis
==========================================================
class: small-code

```r
barplot(bp.go, showCategory=20)
```

<img src="AnalysisOfGeneExpressionDataUsingR-figure/unnamed-chunk-34-1.png" title="plot of chunk unnamed-chunk-34" alt="plot of chunk unnamed-chunk-34" style="display: block; margin: auto;" />

Perform a GO enrichment analysis
==========================================================
class: small-code

```r
mf.go <- enrichGO(gene          = deg_ann.tb$entrezgene_id,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(mf.go)
```

```
                   ID                                 Description GeneRatio
GO:0050839 GO:0050839              cell adhesion molecule binding 493/15137
GO:0045296 GO:0045296                            cadherin binding 328/15137
GO:0004674 GO:0004674    protein serine/threonine kinase activity 427/15137
GO:0003779 GO:0003779                               actin binding 416/15137
GO:0031267 GO:0031267                        small GTPase binding 425/15137
GO:0019787 GO:0019787 ubiquitin-like protein transferase activity 392/15137
             BgRatio       pvalue     p.adjust       qvalue
GO:0050839 499/17697 2.594115e-26 3.107749e-23 2.348357e-23
GO:0045296 331/17697 6.634980e-19 3.974353e-16 3.003202e-16
GO:0004674 439/17697 5.250043e-17 2.096517e-14 1.584224e-14
GO:0003779 431/17697 2.206193e-14 6.607548e-12 4.992963e-12
GO:0031267 443/17697 4.257279e-13 8.551020e-11 6.461539e-11
GO:0019787 407/17697 4.282648e-13 8.551020e-11 6.461539e-11
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  geneID
GO:0050839 LASP1/ITGAL/ITGA3/CX3CL1/BAIAP2L1/PROM1/ADAM22/TENM1/CD9/NISCH/LYPLA2/ANLN/GPRC5A/IGF1/MRE11/NRXN3/SNX1/IBSP/VCL/TIMP2/CASR/CDH1/CDH10/CAPG/CTNNA1/HSPA5/DSG2/LIMA1/LAMA3/PKP2/CDH3/CCN5/CNN2/SLC9A3R2/PKN2/ZC3H15/SLK/CTNNA2/PFKP/PKM/DHX29/CDH19/ADGRL1/ACTN1/EPN2/MARK2/PVR/ADAM11/PICALM/NOTCH3/CDHR2/ENO1/FSCN1/DLG1/ICAM3/ACTN2/CAPZB/FAP/CDH17/PTPRH/PSEN1/CDH7/PKP1/BZW1/COBLL1/ITGB5/ERC1/CHMP2B/COL16A1/RAB10/CTTN/EPS15/CHMP5/ADD1/DOCK9/SNX5/RPL6/GCN1/PXN/FXYD5/LTBP4/NUDC/ICAM1/GOLGA3/LAMB1/ITGA6/EZR/TBC1D2/HSP90AB1/DSP/SH3GLB1/CDHR5/MADCAM1/ARVCF/RANBP1/CRKL/TBC1D10A/MICALL1/PACSIN2/MYH9/RANGAP1/AHSA1/EIF5/NOP56/MAPRE1/CHMP4B/CDH20/VAPA/EMD/STK24/OLFM4/UBFD1/EHD4/TJP1/EIF3E/CCN4/NDRG1/EEF1D/ADAM2/PPP1R13L/FCER2/PLIN3/ICAM5/COMP/ITGB8/PTN/ZC3HAV1/PTPRZ1/YKT6/EIF4H/ECM2/PDLIM1/CXCL12/PFN1/ICAM2/LRRC59/YWHAE/VTN/CPE/RPL34/HSPA8/EHD1/NRXN2/EIF4G2/DDX6/CBL/NECTIN1/CD81/VWF/PPFIBP1/KRT18/PTPN6/STK38/CCN6/CLINT1/THBS4/CDH6/FGF1/RARS1/DBN1/RPL24/ADAM23/ACVR1/ITGB6/ITGA4/SPTBN1/RTN4/FN1/STAT1/TXNDC9/HDLBP/PARK7/CAPZA1/LGALS8/IGSF21/CD46/PRDX1/CNN3/PRDX6/ESYT2/CCN2/SPP1/PPL/TJP2/ITGB1BP1/CTNNAL1/EPCAM/PCMT1/TNN/TGFBI/TMPO/TRIM25/SEPTIN7/CALD1/GIPC1/VAPB/CDH26/LRRFIP1/AHNAK/IL1B/VASP/MMP24/KTN1/TSPAN8/PTPRB/EPS15L1/MACF1/PAICS/KDR/CDC42EP1/CDHR3/MYO1B/CDH15/NECTIN2/SNX9/AFDN/ACTN4/COL5A1/PAK4/LAMA5/EIF2S3/ARHGEF16/EPS8L1/LILRB2/DNAJB1/CC2D1A/ARFIP2/RAN/ITGB4/ANGPTL3/SLC14A2/MPRIP/TRPC4/POSTN/SWAP70/NUMB/LDHA/ARGLU1/ADAMTS8/ANXA1/CEMIP2/P2RX4/TES/FLNB/LCP1/GPNMB/BZW2/DBNL/CIB2/ABI1/NIBAN2/DAB2IP/CCN3/UBAP2/TLN1/TRIM29/RDX/THBS1/PAK6/ADAM10/FNBP1L/RAB1A/EMILIN1/ANXA7/ATIC/IDH1/OLA1/ITGAV/USP8/TMOD3/SEMA7A/FGF2/USO1/ITGB7/DIAPH3/CDH24/FBLN5/MFGE8/UNC45A/IQGAP1/CDH11/CDH13/RPS2/CSNK1D/SH3GL1/SCYL1/PSMB6/EPHA2/EFHD2/SERBP1/CCN1/PTPRF/NECTIN4/HMCN1/CGN/ADAM15/PKP4/PHLDB2/EIF2A/SFRP2/CDH18/TENM2/ARHGAP18/EGFR/NLGN4X/MSN/ITGB1BP2/RPL7A/SH3GLB2/LRRC4C/TNKS1BP1/PTPRJ/TENM4/CDH22/FERMT3/PLCB3/ALDOA/MPP7/ITGB1/CDH8/CD226/ANK3/TWF1/PTPRO/ADAM8/ADAM17/DST/BAG3/UTRN/CAST/ASAP1/PTPRD/THY1/CDH12/PRKCA/CXADR/ADAMTS5/LARP1/CCT8/NPTN/NRG1/MMP14/CCNB2/FMNL2/NCK1/WASF2/CD1D/TAGLN2/F11R/LAD1/ITGB2/ADAMTS13/JAML/ABCF3/ITGA5/STX5/RPL29/NTNG1/VCAM1/PDLIM5/S100A11/CIP2A/S100P/ITGA2/EDIL3/ESM1/CAMLG/YWHAZ/SYK/HNRNPK/GAPVD1/NLGN4Y/JAM3/FBN1/ILK/YWHAB/GOLGA2/IGF2/EEF2/GLOD4/EVPL/SLC3A2/CTNNB1/VASN/SEPTIN2/TNXB/ATXN2L/COL3A1/ADAM9/NPNT/STXBP6/COL4A3/PCBP1/FASN/NLGN1/CTNND2/NLGN2/CDH2/KIF5B/RSL1D1/FGG/FGA/FGB/FRMD5/LAMB2/BSG/PPP1CA/HCFC1/CORO1B/PTPRM/JUP/SPTBN2/CCS/RPL15/KLC2/CKAP5/ARHGAP1/PTPN2/RUVBL1/SFN/BAIAP2/CDK5R1/EPS8L2/CD151/NECTIN3/PLEC/CALR/CDH4/PTPN11/CDH5/PDXDC1/NRXN1/PUF60/EXOC3/PAK2/MB21D2/SCRIB/IST1/ANXA2/CTNNA3/TACSTD2/PKP3/SEPTIN9/H1-10/RAB11B/P4HB/MRTFB/SHTN1/ISG15/RPL14/HMGB1/PTPRT/NLGN3/PTPN1/FLNA/SRC/H3C12/SND1/H3C4/DCHS2/SPTAN1/PARVA/RPS26/CD2AP/RPL23A/ITGBL1/CTNND1/EGFL6/GIGYF2/BMPR2/HSPA1A/RACK1/CD177/SNX2/S1PR3/CLIC1/EMP2/TSPAN4/PPME1/DDX3X/TENM3/PI4KA/TWF2/LYN/EEF1G/ITGB3/S1PR2/H3C8/H3C6/H3C11/H3C1/H3C7/H3C10
GO:0045296                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               LASP1/BAIAP2L1/PROM1/LYPLA2/ANLN/GPRC5A/MRE11/SNX1/VCL/CDH1/CDH10/CAPG/CTNNA1/HSPA5/LIMA1/PKP2/CDH3/CNN2/SLC9A3R2/PKN2/ZC3H15/SLK/CTNNA2/PFKP/PKM/DHX29/CDH19/EPN2/MARK2/PICALM/NOTCH3/ENO1/FSCN1/DLG1/CAPZB/CDH17/PTPRH/PSEN1/CDH7/PKP1/BZW1/COBLL1/ERC1/CHMP2B/RAB10/CTTN/EPS15/CHMP5/ADD1/DOCK9/SNX5/RPL6/GCN1/FXYD5/NUDC/GOLGA3/ITGA6/EZR/TBC1D2/HSP90AB1/SH3GLB1/ARVCF/RANBP1/CRKL/TBC1D10A/MICALL1/PACSIN2/MYH9/RANGAP1/AHSA1/EIF5/NOP56/MAPRE1/CHMP4B/CDH20/VAPA/EMD/STK24/OLFM4/UBFD1/EHD4/TJP1/EIF3E/NDRG1/EEF1D/PPP1R13L/PLIN3/ZC3HAV1/YKT6/EIF4H/PDLIM1/PFN1/LRRC59/YWHAE/RPL34/HSPA8/EHD1/EIF4G2/DDX6/CBL/PPFIBP1/KRT18/STK38/CLINT1/CDH6/RARS1/DBN1/RPL24/ACVR1/SPTBN1/RTN4/STAT1/TXNDC9/HDLBP/PARK7/CAPZA1/CD46/PRDX1/CNN3/PRDX6/ESYT2/PPL/TJP2/CTNNAL1/EPCAM/PCMT1/TMPO/TRIM25/SEPTIN7/CALD1/GIPC1/VAPB/CDH26/LRRFIP1/AHNAK/VASP/MMP24/KTN1/PTPRB/EPS15L1/MACF1/PAICS/KDR/CDC42EP1/CDHR3/MYO1B/CDH15/SNX9/AFDN/PAK4/EIF2S3/ARHGEF16/EPS8L1/DNAJB1/CC2D1A/ARFIP2/RAN/MPRIP/TRPC4/SWAP70/NUMB/LDHA/ARGLU1/ANXA1/CEMIP2/P2RX4/TES/FLNB/BZW2/DBNL/ABI1/NIBAN2/DAB2IP/UBAP2/TLN1/TRIM29/RDX/PAK6/FNBP1L/RAB1A/ATIC/IDH1/OLA1/USP8/TMOD3/USO1/DIAPH3/CDH24/UNC45A/IQGAP1/CDH11/CDH13/RPS2/CSNK1D/SH3GL1/SCYL1/PSMB6/EPHA2/EFHD2/SERBP1/CGN/PKP4/PHLDB2/EIF2A/CDH18/ARHGAP18/EGFR/RPL7A/SH3GLB2/TNKS1BP1/PTPRJ/CDH22/PLCB3/ALDOA/MPP7/ITGB1/CDH8/ANK3/TWF1/PTPRO/BAG3/CAST/ASAP1/CDH12/LARP1/CCT8/CCNB2/FMNL2/NCK1/WASF2/TAGLN2/F11R/LAD1/ABCF3/STX5/RPL29/PDLIM5/S100A11/CIP2A/S100P/YWHAZ/HNRNPK/GAPVD1/YWHAB/GOLGA2/EEF2/GLOD4/EVPL/SLC3A2/CTNNB1/VASN/SEPTIN2/ATXN2L/STXBP6/PCBP1/FASN/CTNND2/CDH2/KIF5B/RSL1D1/BSG/PPP1CA/HCFC1/CORO1B/PTPRM/JUP/SPTBN2/CCS/RPL15/KLC2/CKAP5/ARHGAP1/RUVBL1/SFN/BAIAP2/CDK5R1/EPS8L2/PLEC/CDH4/CDH5/PDXDC1/PUF60/EXOC3/PAK2/MB21D2/SCRIB/IST1/ANXA2/CTNNA3/TACSTD2/PKP3/SEPTIN9/H1-10/RAB11B/MRTFB/SHTN1/RPL14/PTPRT/PTPN1/FLNA/SRC/H3C12/SND1/H3C4/DCHS2/SPTAN1/PARVA/RPS26/CD2AP/RPL23A/CTNND1/GIGYF2/BMPR2/HSPA1A/RACK1/SNX2/CLIC1/PPME1/DDX3X/PI4KA/TWF2/EEF1G/H3C8/H3C6/H3C11/H3C1/H3C7/H3C10
GO:0004674                                                                                                                                                                                                                                                                                                                                                                 CAMKK1/MAP3K14/MAP3K9/CDKL3/MARK4/CDKL5/CAMK1G/CDK11A/DYRK4/MAP4K3/MAP4K5/CLK1/MNAT1/PRKCH/VRK2/MAP2K3/DAPK2/TRIO/PHKA2/LTBP1/MAPK9/EIF2AK2/CDK14/CAMK2B/RIOK2/CDK17/WNK1/TNK2/ADCK1/HIPK2/PKN2/MYLK/MAP2K4/SLK/PRKCQ/CDK13/PHKA1/PRKCZ/ROCK1/PDK3/MAST4/TGFBR3/MAPK6/TESK2/CSNK2A2/CAMK2A/MAP4K4/RPS6KA2/MYO3B/PRKACA/RPS6KA6/SPEG/MARK2/STK10/ALPK1/MAP3K13/NUAK1/MARK3/MAP2K7/PAK3/ARAF/MKNK1/MOK/STK17B/STRADB/GSK3B/ULK2/MAP3K4/CPNE3/MAST2/EIF2AK1/NLK/AURKA/MAPKAPK5/LTBP4/CCNK/IRAK3/MAP3K20/TRPM7/MAP3K1/MYO3A/SRPK1/CDC7/MAST3/MKNK2/MAPK1/GRK3/CDKL1/VRK1/RPS6KA5/SGK2/STK4/CSNK2A1/MYLK2/PAK5/RIOK3/PIM2/CDK16/CAB39L/STK24/MAPK3/PHKB/EEF2K/BCKDK/SGK3/C8orf44-SGK3/RIPK2/IKBKB/STK3/MAP4K1/ERCC2/DMPK/VRK3/AURKC/DYRK1B/AKT2/PRKD2/MAST1/GSK3A/CDK6/PIK3CG/PRKAG2/LIMK1/TGFBR1/TESK1/MAPK8/BMPR1A/MAP3K8/RPS6KB1/MAP2K6/MAPK10/HIPK3/GTF2H1/CAMKK2/GTF2H3/PRKAB1/MAK/MAPK14/STK38/ICK/CCND3/PRPF4B/TTK/CLK4/CSNK1A1/NEK11/MAPKAPK3/ACVR2B/NEK4/ACVR1/PRKAG3/STK16/PASK/STK25/PRKD3/AAK1/MARK1/AKT3/CDK18/NEK2/RPS6KA1/STK11/SGK1/TRPM6/NEK6/NEK9/MASTL/PTK2B/PIK3CA/ACVR2A/CIT/PKN1/CDK2/NRK/ACVR1C/RIOK1/GRK4/STK35/WNK4/PRKCG/MAP2K2/DYRK2/PKMYT1/EIF2AK4/TTBK2/RIPK3/STK33/PAK4/MAP3K10/PNCK/LATS1/PRKAB2/RAF1/PRKAA1/CDK8/DSTYK/DCLK1/CSNK1G2/ADCK2/CDK7/IRAK2/CAMK1/ROCK2/ERN2/CCNH/STK26/TAOK3/SRPK2/MAP3K7/AMHR2/CDK4/ACVR1B/CAB39/NEK3/ALPK3/RPS6KC1/CDK9/PIM1/RIPK1/NEK1/MAP2K5/ITPKA/PAK6/CDK15/PRKG2/BMPR1B/BMP2K/CDKL2/ACVRL1/MAP3K12/TSSK4/ULK3/MYLK3/PDPK1/MINK1/CSNK1D/MAPK4/BRD4/HUNK/SIK1/AKT1/LMTK3/PLK4/MAP3K6/PRKACB/DYRK3/MAP3K21/CDC42BPA/DCAF1/CAMK2D/PLK2/GTF2H2/MYLK4/TTBK1/EGFR/TLK2/CASK/TAF1/CAMK2G/PAK1/ATM/CHEK1/TAOK2/LATS2/CSNK1G3/NEK7/UHMK1/CAMK4/STK32B/PRKCA/LRRK1/TNIK/OBSCN/CDK19/TTN/CDK20/MAPK13/PHKG2/BUB1B/SMG1/DYRK1A/BRAF/PINK1/PSKH1/KALRN/HIPK4/PKN3/BRSK1/TAOK1/SIK3/NEK8/SQSTM1/RPS6KA4/PRKAA2/TSSK3/MAPKAPK2/ERCC3/HIPK1/STK36/NEK10/TGFBR2/NUAK2/PRKCI/DCLK3/SNRK/PRKCD/CAMKV/STK17A/LMTK2/PHKG1/CDK5/FASTK/SYK/WNK2/MELK/STK32C/ILK/MAPK7/PRKCB/PLK1/ATP23/CDK12/RSKR/DAPK3/ULK4/MAP4K2/PBK/MLKL/HJV/MAP2K1/CSNK1G1/STK32A/BUB1/MAP3K2/SIK2/ANKK1/CDK1/DCLK2/AKAP13/PPM1D/PRKCE/CDC42BPG/SOSTDC1/KSR2/EIF2AK3/TP53RK/OXSR1/GRK2/ADCK5/MAP3K11/PLK3/BRSK2/ATR/PDIK1L/RPS6KB2/CLK2/MAP3K19/ULK1/RPS6KA3/NIM1K/HASPIN/FAM20C/TSSK6/ERN1/GAK/AURKB/CLK3/PAK2/MAP3K15/MAPK15/AATK/PRKAG1/LIMK2/CAMK1D/RIPK4/TBK1/CHEK2/PRKX/IRAK1/PRKD1/SRPK3/NRBP2/CDK10/MAPK11/PRKG1/BCR/SBK2/MAPK12/SBK1/LRRK2/STK40/STK31/PIK3R4/WNK3/DAPK1/NEK5/MAP3K5/IRAK4/GRK6/RPS6KL1/PIM3/TLK1/STK39/CDC42BPB/MTOR/ALPK2/STKLD1/GRK5/TOP1/MAP3K3/BMPR2/STK19/CSNK2B/CDKL4/STK38L/CHUK/GTF2H4/CSNK1E/SBK3/NME2/CDK11B/CDK3/PRKDC/IKBKE/STRADA/SIK1B/CCL3
GO:0003779                                                                                                                                                                                                                                                                                                                                                                                                                                                          LASP1/BAIAP2L1/USH1C/SCIN/MYH13/GAS7/SYN1/ANLN/WAS/TMSB10/VCL/MYOM2/MYO16/CAPG/CTNNA1/LIMA1/SYNE2/CAMK2B/CROCC/LIMCH1/CNN2/MYLK/CTNNA2/MYO9A/PFN2/SPTB/WDR1/VASH1/MYO3B/ACTN1/TRPC5/FERMT2/ADD2/FSCN1/KLHL20/CAMSAP3/ACTN2/CAPZB/NEBL/MYH7B/TNS1/CEACAM1/OPHN1/EPB41L2/EPB41L3/COBLL1/SSH1/CTTN/MYBPC2/ADD1/PHACTR3/EPB41L1/FXYD5/MAEA/MYO15A/MYH7/TRPM7/EZR/HDAC6/SORBS1/MYO3A/ABL1/ABLIM1/MYO9B/MISP/TRIOBP/PICK1/MYH9/DAAM1/FERMT1/SNTA1/TNNC2/MYOM1/PLS3/EMD/KLHL4/CORO1A/COTL1/MEFV/SLC6A2/CORO2B/HOMER2/VPS18/BLOC1S6/NCALD/MYH14/COBL/WASL/ACTR3C/BCL7B/CORO2A/PDLIM1/PFN1/SLC6A4/PPP1R9B/MYH1/MYH3/KLHL2/FRG1/KLHL5/PANX1/CORO1C/TRPV4/ARPC3/MYL2/TULP1/PHACTR1/CAP2/WASF1/PHACTR2/PACRG/DBN1/TNNC1/ACTR3/SPTBN1/CCDC88A/MLPH/WIPF1/CAPZA1/GBP1/CNN3/KPTN/TNNT2/MED28/FKBP15/CTNNAL1/ABITRAM/MYOT/PLS1/PDLIM2/KIF18A/HPCA/CXCR4/LDB3/WIPF3/CALD1/GIPC1/OBSL1/POF1B/MYH2/VASP/DSTN/AIF1L/MACF1/VIL1/HIP1/YWHAH/FLNC/MYO1B/MYO5C/TMOD2/INO80/PALLD/CSRP3/AJUBA/TNNI3/CNN1/AFDN/ACTN4/ARPC1B/MAP1S/LSP1/TNNT3/TNNI2/GMFG/HIP1R/SYNE1/EPS8L1/CAP1/PDLIM4/DIAPH1/MAP1B/TNS4/MTSS2/WASF3/MYH8/MYH10/MPRIP/MYH11/MYO18B/ACTR3B/MICAL2/SPIRE1/MYBPC3/HOOK1/FHOD3/AVIL/MICAL1/FHOD1/VILL/FLNB/LCP1/DBNL/MYO1G/ERMN/BIN1/TMOD1/ARPC5L/TLN1/FHDC1/MYO7A/TRPC6/RDX/SPTBN5/ACTR2/MYPN/ITPRID2/TMOD3/SHROOM3/PARVG/FGD4/GAS2L3/DIAPH3/PSTPIP1/TLNRD1/TPM1/SSH2/MYO1F/MYOM3/ABL2/CGN/VASH2/TPM3/ESPNL/MYH15/GC/SNCA/MYO10/IQGAP2/KLHL3/PPP1R18/DAAM2/EGFR/SHROOM2/MSN/DIAPH2/FBXO25/ZNF185/GSN/MSRB2/ADD3/TAGLN/ALDOA/ITGB1/KLHL1/DIXDC1/TWF1/EPS8/DST/PSTPIP2/JMY/PLEKHH2/UTRN/HNRNPU/IMPACT/ENAH/PDLIM3/TMSB4Y/TTN/FMN2/ADCY8/KCNMA1/WHAMM/MYO1E/FMNL2/TMSB15A/WASF2/SHROOM4/PPP1R9A/DMTN/EPB41/TNNI1/PIP/SPTBN4/RUSC1/CCR5/MYL3/ALKBH4/FMNL3/NEXN/ARPC5/XIRP2/PDLIM5/LMOD3/LMOD1/ARPC2/ABLIM2/SHROOM1/MYOZ3/NOS3/MICALL2/CFL2/INPPL1/PKNOX2/CLMN/CACNB2/SYNPO2L/MYO1A/MAP1A/SAMD14/NOD2/MYO5B/TPM4/CORO6/EEF2/PXK/FAM107A/XIRP1/SNTB2/PTK2/ANTXR1/MYO7B/MYRIP/LMOD2/MTSS1/PRKCE/CDK5R2/WIPF2/ENC1/TLN2/SYNPO/SNTB1/MYOZ2/SYNPO2/SNTG2/CORO1B/CFL1/SSH3/ABLIM3/DAG1/SPTBN2/ABRA/GCSAM/MYO1H/MARCKSL1/SYNE3/MYO1D/PFN4/CDK5R1/EPS8L2/PAWR/FLII/MYOZ1/PPP1R42/PLEC/IGSF22/TIGD5/HCLS1/WASHC1/NEB/CTNNA3/FHL3/SETD3/SMTN/SPATA32/FMNL1/ADSS1/PRKN/NF2/FSCN2/MAPT/ESPN/SHTN1/KLHL17/PARVB/TMEM201/LRRK2/MYBPC1/S100A4/EVL/AFAP1/MYO18A/MYO6/MRTFA/PDLIM7/FLNA/ANXA6/GMFB/SVIL/IPP/MIB2/MYO5A/MYH6/SPTAN1/PARVA/MYO1C/NRAP/MYL4/TPM2/MSRB1/EPS8L3/RCSD1/CAPZA2/DMD/INF2/PHACTR4/MACO1/AIF1/SPIRE2/TMSB4X/STK38L/DNASE1/ANG/VPS16/ARPC4/ARPC1A/MICAL3/TWF2/FMN1/SHANK3/CORO7/ANXA8/CYFIP1/MARCKS/PRICKLE4/MYO19
GO:0031267                                                                                              ALS2/DVL2/MYCBP2/FARP2/ARHGAP44/NOX1/CDKL5/PLEKHG6/RANBP9/ANLN/RABGAP1/WAS/NPC1L1/RALBP1/RANBP3/TRIO/CLEC16A/ADRB1/PREX2/ARHGEF5/SIKE1/MCF2L2/TRAPPC3/TBC1D22A/GDI2/RASGRF1/YIPF1/YBX3/NCKAP1/SPHK2/NGFR/IPO5/PKN2/TBC1D22B/TBC1D1/NGEF/GOLGA5/EVI5/ROCK1/TBC1D25/EXOC5/ATP6AP1/LLGL2/PICALM/NSF/ARHGEF10L/RAP1GAP/ARHGEF1/USP33/PAK3/ADCYAP1R1/RAPGEF3/RIMS1/ERC1/XPO1/TNPO1/HACE1/IPO11/DNM1L/DOCK3/RPH3A/KIF16B/BIRC5/ARHGAP4/RAB11FIP3/PLEKHG2/EXOC1/RAPGEF4/UNC13D/NUP50/TBC1D2/MYO9B/RANBP1/LZTR1/TBC1D10A/GGA1/SH3BP1/HPS4/TRIOBP/MICALL1/SBF1/SGSM3/RANGAP1/SOS2/DAAM1/RIN3/RABGGTA/RIMS4/PAK5/KIF3B/MCF2/FGD1/SYTL4/ARHGEF7/NUTF2/PARD6A/TSC2/GGA2/DMXL2/NDRG1/ARHGEF10/ARHGEF18/TBC1D17/PPP6R1/ANKRD27/DENND3/TNPO2/CAV1/TBC1D13/RIC1/RGP1/RAPGEF1/DVL1/DNMBP/RAB11FIP2/TBC1D12/KPNB1/PFN1/SLC6A4/RANGRF/TNFAIP1/IFT20/RAB34/TBC1D9/RAPGEF2/EHD1/ARHGEF17/MADD/EXPH5/CORO1C/KCTD10/ARHGDIB/TBC1D30/WASF1/EXOC2/RASGRF2/ACAP2/ECT2/PEX5L/ARHGEF26/RTKN/MLPH/RAB3GAP1/SOS1/ERRFI1/RAP1A/ARHGEF2/DOCK7/NCF2/EXOC8/RIMS3/RAB29/IPO13/RAB3GAP2/MAPKAP1/DENND1A/PLEKHG1/TBC1D15/OCRL/PKN1/OPTN/NCKAP1L/PREX1/PARD6B/CSE1L/XPO5/GGA3/TBC1D20/MCF2L/PLEKHG3/FGD3/RAB3IP/CDC42EP1/DOCK4/SERGEF/USP6/ARHGEF6/SH3BP4/XPO7/AFDN/PAK4/YIPF2/ARHGEF16/EPS8L1/ARHGEF9/RAB11FIP4/SH3BP5/TBC1D5/RBSN/KIF3A/DIAPH1/EXOC4/WDR44/LLGL1/ARFIP2/TBC1D14/MTSS2/NXT1/RIN2/ARHGEF11/TBC1D8B/IPO8/SBF2/VAV3/ROCK2/DOCK2/BICDL1/MICAL1/RAB11FIP5/TMEM127/ARHGEF4/TBC1D4/RCBTB2/RAPGEF5/RAC1/RALGPS1/ARHGEF39/DENND4C/SYTL2/SORL1/ITPKA/PAK6/RABGGTB/PLCE1/ABI2/RASGEF1B/EGF/FGD4/PEX5/DIAPH3/GRTP1/IQGAP1/ARHGAP17/GAS8/SGSM2/ARHGDIA/VAV1/SOD1/EVI5L/ARHGEF19/SYTL1/STRIP1/RGL1/ADPRH/STXBP5L/TBCK/IQGAP2/TBC1D7/TBC1D7-LOC100130357/DAAM2/FGD2/TIAM2/DENND2A/SYTL5/DIAPH2/DOCK11/C9orf72/USP6NL/PAK1/CDC42EP2/DOCK1/PIH1D2/RILPL2/EPS8/BICD1/RABGAP1L/RASGRP3/FARP1/RGPD3/RANBP2/PLEKHG4B/OBSCN/RABGEF1/FGD5/WHAMM/TIAM1/RAB11FIP1/PDE6D/DGKI/FMNL2/CDC42SE2/RAPGEF6/ABR/KALRN/VAV2/PKN3/DVL3/FMNL3/GRASP/BICDL2/DENND2D/CDC42EP3/ARHGEF3/RANBP3L/WDR41/RHOBTB3/STXBP5/SYTL3/MICALL2/ABCA1/ATP7A/ARHGEF40/HPS6/DENND2B/AP1G1/SGSM1/TBC1D2B/TBC1D16/MYO5B/RAB8A/CDC42EP5/DAPK3/RILP/RAB3IL1/RHOH/STXBP6/XPO6/TBC1D10B/RGPD8/MYRIP/TRAPPC1/DENND5B/AKAP13/PLEKHG5/RASGRP4/KNDC1/SRGAP2C/RASGRP1/MAP3K11/NET1/PIFO/DENND4A/RIN1/DENND6A/KCTD13/SH3BP5L/ARHGAP1/TBC1D10C/DENND2C/RIMS2/RNF152/SMCR8/EPS8L2/ULK1/PARD6G/RCC2/CDC42EP4/RCC1/FGD6/PAK2/TRAPPC5/RPH3AL/RNF41/ANXA2/SPATA13/ARHGEF37/IQGAP3/DENND5A/XPOT/FMNL1/AP3M1/CIB1/RGPD2/ANKFY1/BICD2/BCR/PLEKHG7/SBK2/RILPL1/C15orf62/CHM/NOXA1/LRRK2/PLEKHG4/SRGAP3/SRGAP2B/IPO4/TRAPPC4/STRN3/ARHGEF12/FLNA/SRGAP1/TBC1D9B/MYO5A/CDC42SE1/MYO1C/ITSN2/IPO9/UNC13B/CDC42BPB/EPS8L3/DENND4B/ARHGEF15/RUSC2/RUNDC1/RASGEF1A/INF2/CHML/ECT2L/GDI1/STK19/LSM2/TBC1D8/RANBP17/IPO7/RGL3/DENND6B/ITSN1/DENND1C/DENND1B/ARHGEF33/ARHGEF28/VPS52/ARHGEF38/RGL2/ARHGEF25/ARHGDIG/MICAL3/BRK1/SRGAP2/RASSF5/CYFIP1/TBC1D3D/TBC1D3/TBC1D3E
GO:0019787                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ANKIB1/HECW1/KLHL13/MYCBP2/FBXL3/ASB4/MYLIP/UBE3C/RNF216/BRCA1/UBR7/RNF14/UFL1/RNF10/BIRC3/UBR2/PIAS1/RNF19A/CUL3/TSPAN17/ERCC8/NEDD4L/ANAPC4/RC3H2/RNF4/ASB1/PDZD4/HDAC4/NEDD4/RNF126/CDC42/RAD18/SEL1L/HLTF/LNX1/UBE2D1/CHFR/FBXW11/KLHL20/UBE2T/UBE2A/PIAS2/UBE2K/ITCH/UBE2D4/CNOT4/RNF13/HACE1/ATG16L1/NFX1/HUWE1/KLHL42/ANAPC5/BRAP/MAEA/MUL1/RNF31/G2E3/HECTD1/RFFL/CDC23/MARCHF2/CDC34/RNF215/PPIL2/FBXO7/RBX1/TRIM9/ASB2/BIRC7/RNF24/SEL1L2/RNF125/MIB1/XIAP/UBL4A/MGRN1/STUB1/UBE2I/BFAR/RNF40/HERC1/UBE2W/UBR5/PIAS4/DMAC2/CBLL1/FBXO24/UBE2R2/FBXL15/NEURL1/UBE2S/ZMIZ1/FBXL20/RNF43/TRIM37/RNF167/MED31/PEX12/SMURF2/TNFAIP1/UBE2D3/TRIM2/DTX4/PRPF19/TRIM3/RNF141/BIRC2/UBE4A/CBL/FBXO3/KCTD10/RNF8/FBXO9/ZNF451/FBXL4/RNF130/SKP1/TRIM23/UBE3A/RNF7/CBLB/FANCL/BIRC6/RNF19B/TRIM62/FBXO2/FBXO6/UBE3D/FBXO30/TNFAIP3/RNF146/FBXL5/UBE2B/TRIM32/FBXW2/IRF2BPL/AREL1/RNF170/TRIM25/TRIM6/PDZRN3/RNF2/RBBP6/ZMIZ2/TRIM24/EGR2/RNF11/WWP1/NEURL2/MED20/LRRC29/RNF113A/MED1/RBCK1/HECTD3/TRAF2/FBXL12/UBR4/FBXL16/RNF6/RNF112/HERC2/UBE2M/TRIM28/UBE4B/RLIM/TRAF3/UBE2D2/TRAF7/PIAS3/TRIM21/TRIM5/UBE2G1/TRIM47/FBXO44/RNF128/MED10/MKRN1/RNF122/VHL/DDB2/RNF138/FBXO21/MDM2/FBXL8/RC3H1/LMO7/RNF38/RNF144B/RNF121/FBXO11/BARD1/HECW2/HERC3/HERC6/HERC5/RNF185/PELI2/UBE2Q2/GID4/ANAPC11/RNF157/CBX4/RNF165/FEM1A/CBLC/NOSIP/COP1/UFC1/DTL/UBR3/MARCHF4/UBA3/ATG3/MARCHF1/RNF175/MARCHF6/SKP2/FBXL17/RNF145/RMND5B/RNF217/SHPRH/TAF1/FBXO25/UHRF2/FBXO10/LRSAM1/HERC4/UBE3B/RNF144A/AMN1/FBXO4/ATG10/TRIM36/MED21/ANAPC1/RANBP2/FBXL2/RMND5A/TRIP12/TRIM11/SH3RF1/ELOC/RABGEF1/FBXL18/RNF20/MED7/SH3RF2/UBE2L6/NSMCE2/RNF111/TRIM63/UBE2Z/UBR1/AMFR/UBE2J2/MED27/UBE2Q1/FBXL13/FBXO27/MED11/SYVN1/KLHL21/TRIM58/TRIM17/NEURL3/RNF149/DCST1/RNF25/RCHY1/FBXO40/DTX3L/RNF168/RNF123/RNF180/MED30/RNF183/HECTD2/MARCHF8/PDZRN4/BTRC/ARIH1/CUL5/RAG1/FBXO22/RNF214/TRIM68/ZNF598/KIAA1586/RNF187/RFWD3/RNF181/TMEM129/FEM1B/UBE2V2/NSMCE1/TRIM56/UBE2E3/UBE2E1/RNF150/IRF2BP1/RNF34/RNF139/FBXL14/MALT1/ZNF738/SH3RF3/HECTD4/RNF26/RNF213/MARCHF3/PELI3/MSL2/KCTD13/FBXW8/UBE2C/TRAF6/UBE2O/KCMF1/RNF152/TRIM72/ARIH2/UBE2N/RNF212/RNF186/SHARPIN/ZNRF2/RNF182/PJA1/RNF135/SIAH2/RNF41/UBA7/UBE2E2/FBXL6/TTC3/BCOR/ZNRF3/FBXL7/TRAIP/UBE2F/MED12/UBE2G2/UBOX5/PRKN/UBE2L3/TRIM69/ZNRF1/UBE2H/ZFP91/RNF220/NHLRC1/NCCRP1/NHLRC3/RNFT1/WDSUB1/SIAH1/TRIM33/PELI1/FBXL22/MIB2/TOPORS/MARCHF5/WWP2/KLHL9/SMURF1/UBE2J1/LTN1/ASB12/DZIP3/PJA2/RING1/RNF5/PPP1R11/TRIM27/TRIM13/RNF208/TRIM59/NEURL1B/ZBED1/NEURL4/UBE2QL1/SIAH3/RNF148/RNF103/UBE2V1/MAGEL2/RNF115/FBXO17/UHRF1
           Count
GO:0050839   493
GO:0045296   328
GO:0004674   427
GO:0003779   416
GO:0031267   425
GO:0019787   392
```

```r
#write.csv(bp.go,file="mf.go.csv")
```


Perform a KEGG pathway enrichment analysis
==========================================================
class: small-code

```r
kk <- enrichKEGG(gene         = deg_ann.tb$entrezgene_id,
                 organism     = 'hsa',
                 pAdjustMethod = "BH",
                 qvalueCutoff  = 0.05)
head(kk)
```

```
               ID                    Description GeneRatio  BgRatio
hsa04510 hsa04510                 Focal adhesion  198/6809 199/7946
hsa04360 hsa04360                  Axon guidance  180/6809 181/7946
hsa04144 hsa04144                    Endocytosis  243/6809 249/7946
hsa04010 hsa04010         MAPK signaling pathway  284/6809 295/7946
hsa04722 hsa04722 Neurotrophin signaling pathway  119/6809 119/7946
hsa04150 hsa04150         mTOR signaling pathway  151/6809 153/7946
               pvalue     p.adjust       qvalue
hsa04510 1.041735e-12 3.354385e-10 1.546153e-10
hsa04360 1.643594e-11 2.646186e-09 1.219720e-09
hsa04144 9.927879e-11 1.065592e-08 4.911687e-09
hsa04010 1.086962e-09 8.750041e-08 4.033200e-08
hsa04722 9.004526e-09 5.798915e-07 2.672922e-07
hsa04150 1.566294e-08 8.405780e-07 3.874518e-07
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    geneID
hsa04510                                                                                                                                                                                                                                                                                                                                                                                                                                                     572/3675/3674/2534/3479/3082/330/3381/7414/2324/3371/1298/10319/5601/9564/5291/3909/3918/4659/5923/4638/387/6093/998/87/60/4660/5063/8516/1286/2932/3693/5829/22798/3912/3655/1299/1399/5594/5155/6655/394/85366/10398/57144/10627/284217/331/2321/5595/208/5296/1311/3696/858/857/4233/93408/2889/5599/1277/7448/5602/595/329/7450/4633/1297/896/7422/3910/7060/5159/3694/3676/2335/6654/7143/5906/10000/8503/103910/6696/894/5228/63923/5290/54776/7408/5582/2002/5908/3791/5880/2318/25759/81/10298/3911/1729/5894/3691/10451/9475/5156/1284/3679/3915/2317/5879/1101/7094/7057/22801/56924/106821730/3685/1950/64098/1280/3695/3480/91807/5170/2064/7409/1291/1292/207/8515/3680/56034/5295/340156/1956/53358/5058/3688/7424/1793/5578/673/7791/2909/7410/6464/3678/824/1293/3673/1278/2277/3611/5579/1398/1499/7148/1285/5604/7059/5747/5881/80310/5293/596/5728/83660/3913/5499/256076/7423/3265/3725/2885/29895/5062/859/71/399694/23396/5501/7058/1282/1287/29780/5649/3908/3914/2316/6714/5154/1288/55742/131873/5500/3672/4636/3690
hsa04360                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      6405/56920/2534/7869/2042/5291/816/2773/9423/25791/4756/387/5590/6093/815/998/2046/8440/7224/4775/59277/10371/10512/5362/5063/285220/54437/2932/54434/1947/57556/1808/57715/25/3983/8829/1943/5594/4776/4773/655/10398/57144/10627/5595/50855/7976/5296/4233/2051/3984/6387/219699/5532/1949/81029/9037/90249/7474/2771/5361/2043/8503/103910/5533/5290/7852/84612/5335/1948/6237/2770/5880/64101/6608/10298/55558/5894/7223/2048/3845/9475/2045/10505/10154/5879/7225/56924/106821730/80031/8482/658/7222/5530/5998/85464/1969/10500/1944/151449/91653/7220/9353/2044/817/5295/5921/2041/56288/818/57689/5058/3688/5163/223117/64221/5578/55740/2047/137970/56896/4690/54361/22854/6259/84448/5364/6469/1020/1073/3611/10501/54910/1942/5747/5881/6091/9723/5293/27289/1072/54961/29984/22885/3265/84552/5781/5062/8633/56963/2242/3985/2049/6586/1946/6092/10509/5727/6585/1630/10507/5535/64218/9901/84628/2050/23654/23365/57522/6714/5336/5365/3897/659/4893/4636/5534/91584/1945/23380
hsa04144                                                                                                           381/163/9265/51534/30845/30846/6642/9135/9765/1175/116984/55040/23327/5662/29924/2263/10565/387/5590/2261/55737/4734/8218/998/22905/7037/9744/7251/5337/7879/832/273/83737/1785/25978/64744/10890/2060/51510/55616/27131/51324/9727/51100/387680/27128/157/3560/55738/128866/8411/50855/8766/30844/9266/858/857/8976/7046/1759/8395/6456/22841/8729/28964/9267/5878/64750/3312/10938/2348/867/10094/51699/5869/11021/60682/8724/23527/868/8723/10254/51652/7456/829/9829/9525/64411/440073/51019/7852/644150/1211/9559/11059/84612/10564/57403/23096/10617/2868/23550/3306/58513/5338/58533/3949/10095/27243/5119/84440/64145/84313/116983/27183/23111/116988/377/3559/5156/116986/26056/4193/23325/51028/11311/274/8027/81873/408/9101/9815/79720/3480/5371/6457/1213/409/93343/6455/23624/8394/375/5868/9922/84249/1956/30011/29934/3561/91782/56904/56288/84364/112936/8853/1601/50807/118813/137492/3798/23362/80223/9372/26119/116987/11267/1234/2264/51160/1173/10092/3577/10109/7048/5584/92421/9897/2350/382/4088/138429/4218/55048/5867/3800/11031/10015/7314/3799/147179/729092/253725/254122/156/3310/3265/7189/4087/1212/79643/155382/80230/84552/3579/100287171/10193/9798/859/161/390243/9230/9146/23396/116985/119016/89853/160/3482/6714/26052/2870/57154/2869/830/414189/642517/3304/3303/3305/3107/3133/3135/3134/6643/3105/3106/10093/10552/26286/29082/57132/119385
hsa04010 9020/8913/9254/23542/5536/8491/3479/3082/355/5606/2324/4254/5601/23118/5923/4804/2065/6416/2263/7132/782/2261/10235/2255/27006/998/9448/6196/5566/27330/4775/9175/2249/27092/5609/4791/2260/369/26281/8569/4208/779/2122/4216/51701/8550/2323/51776/285/7042/4214/1943/4616/2872/1399/5594/5155/10454/8911/6655/5494/9252/6789/994/51378/778/2254/2321/8717/5595/3551/6788/11184/5971/208/7040/8605/4233/3315/7046/5599/5532/1326/1845/786/5608/8649/4790/374/5602/9693/3312/84867/80824/1432/6722/7422/5924/2246/5159/7867/3552/3554/6654/1386/5906/5321/1647/10000/356/3925/6195/5228/7043/1843/7010/1846/5533/9479/2322/3164/2069/3553/4149/6237/5582/2002/3306/5605/7186/5908/3791/468/5880/2318/51295/10912/59285/3727/1852/4772/5894/3845/22800/4803/5156/51347/6885/2317/5879/4609/408/5607/5495/1847/2247/1950/5530/23162/1848/7786/2252/3480/409/7157/2064/773/207/59283/1969/9064/5567/11221/1944/5778/56034/5921/2768/1956/4915/774/5058/9344/7424/93589/775/25780/5801/781/284/5578/5922/8817/5603/776/3815/55799/673/1844/22808/2005/8822/283748/57551/2264/8986/9965/9261/7039/7048/1849/836/2277/783/5598/5579/1398/3481/784/5871/255189/8681/5604/1942/5881/10746/2353/929/80310/3643/115727/55970/10125/4615/5970/3310/4296/7423/3265/7189/1649/627/6197/3725/2885/2066/5062/785/1436/3654/1946/1435/1850/5600/4908/4137/2248/123745/6300/5535/3556/8912/4763/2316/4217/5154/51135/777/4914/4215/1616/3304/3303/3305/4893/1147/5534/4909/100506012/7124/1945/100137049/3630/8517
hsa04722                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     572/4145/5601/5291/816/4804/387/815/998/6196/27330/5609/7161/5663/2932/581/11213/4214/25/1399/5594/6655/9252/4792/5595/8767/3551/4793/208/5296/2889/5599/7531/4790/5602/2549/11108/10019/397/1432/6654/5906/10000/8503/356/6195/2309/5290/5335/5605/5908/468/91860/25759/5894/3845/3656/6272/4803/57498/5879/5607/4916/5170/7157/396/207/5664/805/817/5295/4794/4915/53358/818/814/5603/673/808/6464/10603/9261/5580/10818/5598/27018/1398/5604/3667/163688/10782/5293/596/5970/3265/7189/627/6197/3725/2885/25970/810/51806/9500/5781/3654/5600/399694/4908/6300/4217/5336/51135/4914/801/4215/4893/4909/398
hsa04150                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            51384/1856/3479/58528/51606/5291/57600/1975/7132/387/4040/6196/27330/7479/55437/2932/9706/10325/81929/7481/5594/9681/6655/51382/9894/9675/81617/5595/79726/8131/7249/8140/7976/3551/208/5296/7472/2887/6009/1855/7473/6198/8649/81029/11211/7474/10641/523/7475/6654/525/28956/64121/10000/8503/6195/6794/6446/79109/5290/7471/5582/5605/84219/9296/83667/6249/529/5894/5562/3845/7482/10542/23175/80326/9470/51719/9550/6194/3480/5170/57521/207/7483/245973/6502/5295/51256/107080638/526/55004/1977/127124/5578/89780/7476/201163/528/8324/64798/10670/6396/8321/673/7484/1857/4041/5563/54361/7855/253260/54468/8323/7248/5579/64223/6520/54541/5604/3667/7480/3643/5293/5728/3265/8322/6199/220441/153129/8408/6197/8325/2885/2535/55615/1978/7477/389541/8326/2475/84335/4893/1147/534/96459/7124/652968/90423/3630/92335/729438
         Count
hsa04510   198
hsa04360   180
hsa04144   243
hsa04010   284
hsa04722   119
hsa04150   151
```

```r
#write.csv(kk,file="enrichedKeggPathways.csv")
```

Perform a KEGG pathway enrichment analysis
==========================================================
class: small-code

```r
dotplot(kk, showCategory=30) + ggtitle("dotplot for KEGG")
```

<img src="AnalysisOfGeneExpressionDataUsingR-figure/unnamed-chunk-37-1.png" title="plot of chunk unnamed-chunk-37" alt="plot of chunk unnamed-chunk-37" style="display: block; margin: auto;" />
