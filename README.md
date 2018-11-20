![glass-logo](https://user-images.githubusercontent.com/9220167/48782377-13789600-ecac-11e8-8598-771b6e5d75c6.png)

## The GLASS consortium

### Overview
The Glioma Longitudinal AnalySiS (GLASS) consortium consists of clinical, bioinformaticians, and basic science researchers from leading institutions across the world striving to better understand glioma tumor evolution and to expose its therapeutic vulnerabilities. For more information, please read our position paper: **The Glioma Longitudinal Analysis Consortium. Glioma Through the Looking GLASS: Molecular Evolution of Diffuse Gliomas and the Glioma Longitudinal AnalySiS Consortium.** *Neuro Oncol. 2018* PMID: 29432615. 

#### Analytical contributors
*In no particular order*: Floris Barthel, Kevin Johnson, Samir Amin, Fred Varn, Anzhela Moskalik, Hoon Kim, Roel Verhaak

#### Data contributors
*In no particular order*: Mustafa Khasraw (U Sydney), Joe Costello (UCSF), Tom Mikkelsen (Henry Ford Hospital), Dohyun Nam (Samsung Medical Center), Pieter Wesseling (VUMC), Antonio Iavarone (Columbia), Gaetano Finochiaro (Instituto Besta), Lucy Stead (Leeds U), Adelheid Wöhrer (Vienna), Hiromichi Suzuki (Kyoto U), Priscilla Brastianos (MGH), Jason Huse (MD Anderson), John de Groot (MD Anderson), Kristin Alfaro-Munoz (MD Anderson), and NIH/NCI The Cancer Genome Atlas.

### The GLASS Data Model

The GLASS data model was built to carefully take into the account the temporal nature of the dataset. The layers of the dataset is best explained using the diagram below. Here we see a single hypothetical patient (left-most column) from which we have obtained tumor samples at three different time points (second column), including a second tumor sub-sample (eg. multisector sample) from the primary tumor (third column). Moreover, a simultaneous DNA-RNA extraction was performed on 1st tumor recurrence and both were sequenced (fourth column). Lastly, some of the DNA aliquots were split off and separately sent off for whole genome and whole exome sequencing (fifth column).

The gray colored blocks in the diagram represent the abstraction layers we are using to represent this complexity. Firstly, the `case` level represents patient in the cohort. Secondly, the `sample` level represent bulk samples (tumor and control) in the cohort and finally the `aliquot` level is the aggregate of the sample portion used for differentiating between multi-sector samples, the analyte used for differentiating between DNA and RNA extractions and finally the analysis to which the aliquot was subjected, such as whole genome, whole exome or RNA sequencing.

![data-model](https://user-images.githubusercontent.com/9220167/48782460-3dca5380-ecac-11e8-8ac5-c3c2d71bb94a.png)

Internally, we are using the PostgreSQL database management system to manage the complex relationships between data entities. More information on the database scheme can be found [here](https://www.synapse.org/#!Synapse:syn17038081/wiki/585706).

#### GLASS Barcode

The GLASS barcodes are inspired by the TCGA barcodes and are constructed in a recognizable fashion:

![barcode](https://user-images.githubusercontent.com/9220167/48782475-491d7f00-ecac-11e8-9bba-c4ba02e8aa07.png)

Key parts of the data model are represented in the barcode.

### Data Release version 2018-11-15
The data available here marks the initial release of the GLASS project dataset, which is managed internally using the PostgreSQL database management system. The data available here represents the `2018-11-15` snapshot of our dataset. These data are under active curation so future versions will include additional data as well as correct potential errors. Please report any inconsistencies you find on the discussion board here or on our Github issue tracker https://github.com/TheJacksonLaboratory/GLASS.

#### Current Data by the Numbers

The following infographic shows the some summary statistics of the `2018-11-15` release of the GLASS dataset.

![glass-numbers](https://user-images.githubusercontent.com/9220167/48782498-53d81400-ecac-11e8-807a-45f9b67e0e5f.png)

#### Data Download
The GLASS data can be downloaded from the `Tables` page [here](https://www.synapse.org/#!Synapse:syn17038081/tables/). It is also possible to query the data directly using the the API by using queries. You can read more about that [here](https://docs.synapse.org/articles/tables.html).

Royalty-free icon made by Darius Dan of <a href="https://www.flaticon.com/" title="Flaticon">www.flaticon.com</a> is licensed by <a href="http://creativecommons.org/licenses/by/3.0/" title="Creative Commons BY 3.0" target="_blank">CC 3.0 BY</a></div>
