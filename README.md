# NCBoost v2
[]()  

NCBoost is a pathogenicity score of non-coding variants to be used in the study of Mendelian diseases. It is based on supervised learning on a comprehensive set of ancient, recent and ongoing purifying selection signals in humans. NCBoost was trained on a collection of 2336 high-confidence pathogenic non-coding variants associated with monogenic Mendelian diseases. NCBoost performs consistently across diverse independent testing data sets and outperforms other existing reference methods. Further information can be found at the [NCBoost paper]().

Of note, the NCBoost software can score any type of genomic position, provided that the required puryfing selection features used by the model are available. However, it is important to realize that, among the set of high-confidence pathogenic non-coding variants that were used to train NCBoost, more than 98%  were found at proximal cis-regulatory regions, with only 27 variants overlapping more distal intergenic regions. Thus, for consistency with the training process, the most reliable genomic contexts for the use of the NCBoost score are the proximal cis-regulatory regions of protein-coding genes.

## Precomputed NCBoost scores in proximal cis-regulatory regions of protein-coding genes

We precomputed the NCBoost score for XXXXXX non-coding genomic positions overlapping intronic, 5'UTR, 3'UTR, upstream and downstream regions -i.e. closer than 1kb from the Transcription Start Site (TSS) and the Transcription End Site (TES), respectively- associated with a background set of [19433 protein-coding genes](https://github.com/RausellLab/NCBoost-2/tree/master/data#file-genedb_ncboost2.tsv) for which we could retrieve annotation features. Variant mapping and annotation of non-coding genomic positions was done through [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/user-guide/download/) software using the gene-based annotation option based on RefSeq (assembly version hg38). In the case of positions overlapping several types of regions associated with different genes and transcripts (either coding or non-coding), a number of criteria were adopted as described in the [NCBoost v2 paper]().

The precomputed hg38 NCBoost 2 scores in proximal cis-regulatory regions of protein-coding genes can be downloaded [here]() as a tabix indexed file (gz):
https://storage.googleapis.com/ncboost-cbl/ncboost_score_hg38_v2025XXXX.tsv.gz
and the corresponding index file is available [here]() (gz.tbi):
https://storage.googleapis.com/ncboost-cbl/ncboost_score_hg38_v2025XXXX.tsv.gz.tbi

The file contains the following columns:  
*chr*, chromosome name, as [1:22,X,Y]  
*pos*, 1-based genomic position (GrCh38 genome assembly)  
*region*, type of non-coding region overlapping the position, as provided by ANNOVAR (see above)  
*closest_gene_name*, name of the associated protein-coding gene  
*NCBoost_Score*, NCBoost v2 score. NCBoost score ranges from 0 to 1. The higher the score, the higher the pathogenicity potential of the position.  
*NCBoost_chr_rank_perc*, chromosome-wise rank percentile (ranging from 0 to 1) of the corresponding NCBoost v2 score. The higher the rank percentile, the higher the pathogenic potential of the position.  

## NCBoost gene database
The NCBoost gene database, integrating several database identifiers (Ensembl, HGNC, NCBI), OMIM disease-gene status and the gene-level conservation features used in this work are available [here](https://github.com/RausellLab/NCBoost-2/tree/master/data/geneDB_ncboost2.tsv).

## NCBoost software
The NCBoost software is also provided in this repository in case you are interested in training the NCBoost framework on your own variants, or assessing the NCBoost scores for genomic positions other than those included in the precomputed file.

The following sections will guide you through the steps needed for the annotation of variants, training and execution of NCBoost-2 pretrained models to obtain the pathogenicity score.


## Downloads, installation and processing of input files

### 1. Download NCBoost v2 software

NCBoost scripts and associated data may be cloned from the NCBoost github repository:
```
git clone https://github.com/RausellLab/NCBoost-2.git
cd NCBoost-2
```
The required python libraries are detailed in [libraries.txt](https://github.com/RausellLab/NCBoost-2/blob/master/libraries.txt).

Python3 libraries can be installed using:  
`pip install -r libraries.txt`

### 2. Download the feature file

NCBoost 2 features are available [here](). Compressed tabix-indexed files are provided for each chromosomes (total size = XXX Go), and can be downloaded using the following command:
```
wget ...
```
Move the downloaded data to data/WG_annotated/
```
mkdir data/WG_annotated
mv xxx data/WG_annotated/
```

Complete details about each features are available at [NCBoost v2 paper]().

## Variant input format
Variants have to be reported in 1-based, GrCh38 genomic coordinates. The variant file required is a tab-delimited textfile with column headers following the format:
```
chr start   ref  alt
1   12589   G   A
```

The *chr* column should not contain the 'chr' prefix.
Other columns can be added in addition to such first four columns.

## NCBoost training
NCBoost framework can be trained using the ncboost_train.ipynb script. It loads and annotate a set of pathogenic variants and the corresponding set of region-matched random common variants, train the 10 models and produce the corresponding feature importance plot, as well as ROC and PR curves. The trained models are saved and can be used for later scoring.

The annotation requires to download the full set of features used by NCBoost (XXXGo). For convenience, we also provide the set of pathogenic and common variants already annotated with NCBoost features.

## NCBoost scoring
NCBoost framework can be applied to score any variant using the ncboost_score.ipynb script. It will apply the trained framework used to generate the resutls in [NCBoost v2 paper].

The output file is a tab-delimited text file displaying by columns the following fields (in this order): The chromosome, position, reference and alternative allele of the variant, the name and Ensembl Gene ID of the nearest gene to which the variant was associated and the corresponding non-coding region (upstream, downstream, UTR5, UTR3, intronic and intergenic), the gene type and 11 gene-based features (slr_dnds, gene_age, pLI, zscore_mis, zscore_syn, loeuf, GDI, ncRVIS,
ncGERP, RVIS_percentile, pcGERP), using a reference of [19433 protein-coding genes](https://github.com/RausellLab/NCBoost-2/tree/master/data#file-genedb_ncboost2.tsv), 6 one-hot encoded non-coding region types, 18 features extracted from CADD annotation files [[5]](https://github.com/RausellLab/NCBoost#references) (GC, CpG, pri/mam/verPhCons, pri/mam/verPhyloP, GerpN, GerpS, GerpRS, GerpRSpval, ZooPriPhyloP, ZooVerPhyloP, bStatistic, ZooRoCC, ZooUCE, Roulette-AR), 9 positive-selection scores (TajimasD_YRI/CEU/CHB_pvalue, FuLisD_YRI/CEU/CHB_pvalue, FuLisF_YRI/CEU/CHB_pvalue), the mean DAF and mean Het, the 9 MAF from the 1000GP or GnomAD [[6]](https://github.com/RausellLab/NCBoost#references) (mean_MAF, mean_MAF_AFR/AMI/AMR/ASJ/EAS/FIN/MID/NFE/SAS), the CDTS score, the SpliceAI score and the NCBoost score and the extra columns provided by the user in the input file.  
NCBoost score range from 0 to 1 (the higher the value, the higher the predicted pathogenicity).  

More information about the can be found in [NCBoost 2 paper]().  


## References
1: Wang and Hakonarson; (2010). ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. Nucleic Acids Res. 38, e164-e164.

2: di Iulio et al. (2018). The human noncoding genome defined by genetic diversity. Nat. Genet. 50, 333-337.

3: Ritchie et al. (2014). Functional annotation of noncoding sequence variants. Nat. Methods 11, 294-296.

4: Pybus et al. (2014). 1000 Genomes Selection Browser 1.0: a genome browser dedicated to signatures of natural selection in modern humans. Nucleic Acids Res. 42, D903-D909.

5: Kircher et al. (2014). A general framework for estimating the relative pathogenicity of human genetic variants. Nat. Genet. 46, 310-315.

6: Lek et al. (2016). Analysis of protein-coding genetic variation in 60,706 humans. Nature 536, 285-291.

## Contact
Please address comments and questions about NCBoost to barthelemy.caron@institutimagine.org and antonio.rausell@inserm.fr

## License
NCBoost 2 scripts, framework and databases are available under the [Apache License 2.0](https://github.com/RausellLab/NCBoost-2/tree/master/LICENSE).

Copyright 2025 Clinical BioInformatics Laboratory - Institut Imagine

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
