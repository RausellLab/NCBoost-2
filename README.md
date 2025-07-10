# NCBoost v2
[]()  

NCBoost is a pathogenicity score of non-coding variants to be used in the study of Mendelian diseases. It is based on supervised learning on a comprehensive set of ancient, recent and ongoing purifying selection signals in humans. NCBoost was trained on a collection of 2336 high-confidence pathogenic non-coding variants associated with monogenic Mendelian diseases. NCBoost performs consistently across diverse independent testing data sets and outperforms other existing reference methods. Further information can be found at the [NCBoost paper]().

Of note, the NCBoost software can score any type of genomic position, provided that the required puryfing selection features used by the model are available. However, it is important to realize that, among the set of high-confidence pathogenic non-coding variants that were used to train NCBoost, more than 98%  were found at proximal cis-regulatory regions, with only 27 variants overlapping more distal intergenic regions. Thus, for consistency with the training process, the most reliable genomic contexts for the use of the NCBoost score are the proximal cis-regulatory regions of protein-coding genes.

## Precomputed NCBoost scores in proximal cis-regulatory regions of protein-coding genes

We precomputed the NCBoost score for XXXXXX non-coding genomic positions overlapping intronic, 5'UTR, 3'UTR, upstream and downstream regions -i.e. closer than 1kb from the Transcription Start Site (TSS) and the Transcription End Site (TES), respectively- associated with a background set of [19433 protein-coding genes](https://github.com/RausellLab/NCBoost/blob/master/data/geneDB_ncboost2.tsv) for which we could retrieve annotation features. Variant mapping and annotation of non-coding genomic positions was done through [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/user-guide/download/) software using the gene-based annotation option based on RefSeq (assembly version hg38). In the case of positions overlapping several types of regions associated with different genes and transcripts (either coding or non-coding), a number of criteria were adopted as described in the [NCBoost v2 paper]().

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

### 1. Software and database downloads

NCBoost scripts and associated data may be cloned from the NCBoost github repository:
```
git clone https://github.com/RausellLab/NCBoost-2.git
cd NCBoost-2
```
The required python libraries are detailed in [libraries.txt](https://github.com/RausellLab/NCBoost-2/blob/master/libraries.txt).

Python3 libraries can be installed using:  
`pip install -r libraries.txt`

### 2. Download and processing of the feature file

#### A - Download of feature files

NCBoost 2 features are available [here](). Compressed tabix-indexed files are provided for each chromosomes (total size = XXX Go), and can be downloaded using the following command:
```
wget ...
```
Once downloaded, 
The gene-level features used by NCBoost are provided as part of this repository and described [here](NCBoost_data).  

The script uses wget and gsutil to download the different data files from the original sources in the NCBoost_features/ folder.
It may be run using the following command line from the root NCBoost folder:

The total size of the feature data files is about 450Go, so it might be long...
Complete details about the associated source feature files are provided [here](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_features).

#### B - Indexing of feature files

To make the variant annotation more efficient, most of the feature mapping steps are done using tabix or pytabix. Some feature files need to be indexed at this step. 
To index the files using [tabix](http://www.htslib.org/doc/tabix.html), you may run the script `prepare_feature_databases.sh` from the main NCBoost folder.
```
./NCBoost_scripts/prepare_feature_databases.sh
```
The script will index the bed and vcf files containing CDTS scores [[2]](https://github.com/RausellLab/NCBoost#references), mean Het and mean DAF [[3]](https://github.com/RausellLab/NCBoost#references) and the features provided by the 1000 Genome Selection Browser [[4]](https://github.com/RausellLab/NCBoost#references). See [NCBoost paper](https://rdcu.be/bmlxX) Table 1 and Methods section for further details.


## Variant input format
Variants have to be reported in 1-based, GrCh37 genomic coordinates. The variant file required is a tab-delimited textfile with column headers following the format:
```
chr start   end   ref  alt
1   12589   12589   G   A
```

The *chr* column should not contain the 'chr' prefix.
For Single Nucleotide Variants (SNV), *start* and *end* should be equal (ANNOVAR's requirement), both corresponding to the position of the SNV.
Other columns can be added in addition to such first five columns. However, extra colums may slow down the annotation pipeline.

## Feature annotation and scoring of variants
The whole annotation process and the final scoring by NCBoost can be performed using the script `ncboost_annotate.sh`, with two arguments: (i) the path to the variants input file; and (ii) the path to the output file.

```
./NCBoost_scripts/ncboost_annotate.sh /folder1/subfolderB/inF.vcf /folder1/subfolderB/outF.tsv
```
The output file is a tab-delimited text file displaying by columns the following fields (in this order): The chromosome, position, reference and alternative allele of the variant, the nearest gene to which the variant was associated and the corresponding non-coding region (as determined through ANNOVAR, see above), the gene type and 8 gene-based features (pLI, familyMemberCount, ncRVIS, ncGERP, RVIS_percentile, dnds, GDI and gene_age), using a reference of [18404 protein-coding genes](https://github.com/RausellLab/NCBoost/blob/master/NCBoost_data/NCBoost_geneDB.tsv), 6 one-hot encoded non-coding region types, 11 features extracted from CADD annotation files [[5]](https://github.com/RausellLab/NCBoost#references) (GC, CpG, pri/mam/verPhCons, pri/mam/verPhyloP, GerpN, GerpS, bStatistic), 9 positive-selection scores (TajimasD_YRI/CEU/CHB_pvalue, FuLisD_YRI/CEU/CHB_pvalue, FuLisF_YRI/CEU/CHB_pvalue), the mean DAF and mean Het, the 9 MAF from the 1000GP or GnomAD [[6]](https://github.com/RausellLab/NCBoost#references) (meanMAF1000G, meanMAFGnomAD, meanMAF_AFR/AMR/ASJ/EAS/FIN/NFE/OTHGnomAD), the CDTS score, the NCBoost score and the extra columns provided by the user in the input file.  
NCBoost score range from 0 to 1 (the higher the value, the higher the predicted pathogenicity).  
More information about the gene-based features [here](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_data#file-ncboost_genedbtsv), about the other conservation features [here](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_features#feature-details), or in [NCBoost paper](https://rdcu.be/bmlxX).  

## Example
An example on how to run NCBoost software is available [here](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_example).

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
NCBoost scripts, framework and databases are available under the [Apache License 2.0](https://github.com/RausellLab/NCBoost/tree/master/LICENSE).

Copyright 2018 Clinical BioInformatics Laboratory - Institut Imagine

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
