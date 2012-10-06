# Welcome to HMM sample classification
=========
## Abstruct
Multiple myeloma is a cancer of plasma cells that normally produces antibodies. Aneuploidy, the alterations in the number of whole or partial chromosomes, is consistently observed in many cancer types including multiple myeloma. According to the status of chromosomes in multiple myeloma, this cancer can be subdivided into two main subtypes: hyperdiploid multiple myeloma (HMM) and non-hyperdiploid multiple myeloma (NHMM). The two subtypes have different survival prognosis, possibly due to their different but converging paths to myeloma. The existing method to identify the two subtypes is fluorescence in situ hybridization (FISH), and there is no effective method using gene expression profiles to classify the two subtypes. We have built a nearest-neighbor based classification method to separate multiple myeloma into HMM and NHMM with gene expression profiles, providing a useful tool for further studying of multiple myeloma. We also provide publicly the R packages and the processed training data sets at [there](http://hmmc.github.com).

```Classifying hyperdiploidy status of multiple myeloma samples using gene expression profiles```

## Download

Packages: window [SampleClassify_1.0.zip](data/sampleClassify.zip) Mac [SampleClassify_1.0.tar.gz](data/SampleClassify_1.0.tar.gz)

Training data: [TrainData.tar.bz2](data/TrainData.tar.bz2)

Test data(not whole): [GSE6401_RAW.tar.bz2](data/GSE6401_RAW.tar.bz2)

Source code: [mmsc.r](data/mmsc.r)

## Usage

### *install*

Package affy
	
	> source("http://bioconductor.org/biocLite.R")
	> biocLite("affy")

Package AnnotationDbi

	> source("http://bioconductor.org/biocLite.R")
	> biocLite("AnnotationDbi")
	
Data CDF file <http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp>

	run R GUI --> package/data --> package installer --> local source package.	
```Notes:Don't use command install.packages(), It's not work!```

Install SampleClassify

	run R GUI --> package/data --> package installer --> local source package.

	
### input

*Raw_Dir* | *Exp_File* : Raw_Dir for raw data directory and Exp_File for expression data file name. The format of raw data is compressed .cel from GEO database by affy matrix technology. You can input raw data directly, then the program will first calculate the expression by R package _affy_. Otherwise, you must provide common style expression data. You is not smart with provide both Raw_Dir and Exp_File.

training data : This baic data. 

cdf name : ref package affy.

[out put] : A directory contain all of result and temp file. By default, It's "./OutPut" 

[k] : .. By default, it's 16.

[run time] : .. By default, it's 1:140.

`Note:[ parameter ] parameter is optional.`

### output

ExprData_MM_sort.txt :

Expression_Temp.XLS :

MeanDistPlot :

MMSCresult_k16.txt :

### example


	#set workshop
	setwd("~/Documents/YingxiangL/");
	
	#overview the file and directory
	dir();
	
	#load package
	library(SampleClassify);
	
	#analysie data
	mmsc(Raw_Dir=“GSE6401_RAW/”, Tra_Dir="TrainData/", CDF_Name="HGU133a_hs_refseq");

### screenshot

![image](image/example.png)

## About

This a instest work. Mou, github, R, biology, statistic and so on… 
