mmsc <-
function(
				Raw_Dir  = "",
				Exp_File = "",
				CDF_Name,
				Tra_Data,
				Out_Dir  = "./OutPut",
				K_Const  = 16,
				Run_Time = 1:140)
{
	# switch
		switch_path <- function(string)
		{
			if(substr(string,nchar(string),nchar(string))!=.Platform$file.sep)
				string=paste(string,"",sep=.Platform$file.sep)
			return(string)
		}
	# raw(cel) --> expression(tab)
		rp <- function(rawdir,cdf,expfile){
			RawData         <- ReadAffy(celfile.path=rawdir)
			RawData@cdfName <- cdf
			RMAData         <- rma(RawData)
			ExprData_MM     <- exprs(RMAData)
			write.table(ExprData_MM, file=expfile)
		}
	# Classify for expression
		cf <- function(expfile,outdir,k,AllSampleid,td){
			print("multiple myeloma subtype classify")
			ExprData_MM         <- read.table(expfile)
			TrisomyChro_Gene    <- td$tri[which(td$tri[, 1] ==  'TrisomyChroGene'), 2]
			NontrisomyChro_Gene <- td$tri[which(td$tri[, 1] ==  'NontrisomyChroGene'), 2]
			ExprData_MM_scale   <- scale(ExprData_MM)
			##Get the tri-chro mean and adjusted nontri-chro mean of dataset.
			TriNontriChro_Mean <- function(Data_MM_scale, FeatureGene){
				TriChro_Gene <- intersect(FeatureGene, TrisomyChro_Gene)
				NontriChro_Gene <- intersect(FeatureGene, NontrisomyChro_Gene)
				Data_MM_scale_Tri <- Data_MM_scale[na.omit(match(TriChro_Gene, rownames(Data_MM_scale))), ]
				Data_MM_scale_Nontri <- Data_MM_scale[na.omit(match(NontriChro_Gene, rownames(Data_MM_scale))), ]
				#tri-chro mean = ∑ (value of tri-chro gene)/ amount of tri-chro gene
				TriChro_Mean <- apply(Data_MM_scale_Tri, 2, function(x) mean(x))
				#nontri-chro mean = ∑ (value of nontri-chro gene)/ amount of nontri-chro gene
				NontriChro_Mean <- apply(Data_MM_scale_Nontri, 2, function(x) mean(x))
				TriNontriChro_Mean_Output <- cbind(TriChro_Mean, NontriChro_Mean)
				TriNontriChro_Mean_Output
			}
			##Get the adjusted tri-chro mean and adjusted nontri-chro mean of dataset.
			Adjusted_TriNontriChro_Mean <- function(Data_MM_scale, FeatureGene){
				TriNontriChro_Mean_Input <- TriNontriChro_Mean(Data_MM_scale, FeatureGene)
				TriChro_Mean <- as.numeric(TriNontriChro_Mean_Input[, 1])
				NontriChro_Mean <- as.numeric(TriNontriChro_Mean_Input[, 2])
				#adjusted tri-chro mean = tri-chro mean value - ∑ (value of tri-chro mean)/ amount of tri-chro mean
				Adjusted_TriChro_Mean <- TriChro_Mean- mean(TriChro_Mean)
				#adjusted nontri-chro mean = nontri-chro mean value - ∑ (value of tri-chro mean)/ amount of tri-chro mean
				Adjusted_NontriChro_Mean <- NontriChro_Mean- mean(NontriChro_Mean)
				Adjusted_TriNontriChro_Mean_Output <- cbind(Adjusted_TriChro_Mean, Adjusted_NontriChro_Mean)
				Adjusted_TriNontriChro_Mean_Output
			}
			##Use training data to predict test data with KNN in different k value.
			KNN_k_Predict <- function(TrainingData, TestData, TrainingDataClass, k){
				PredictedSubtype <- knn(TrainingData, TestData, TrainingDataClass, as.numeric(k), prob = T)
				PredictedSubtype
			}
			##We predict each sample of input data for 140 times, because we assume each sample was one sample of leave-one-sample-out analysis.
			##We defined that in the 140 times of prediction, the subtype which is predicted more is the final subtype of this sample.
			PredictSubtype <- function(KNN_k_Classify){
				PredictedSummary <- c()
				for (i in 1:length(KNN_k_Classify[, 1])){
					OneSummary <- summary.factor(KNN_k_Classify[i, (1:140)* 2- 1])
					OnePredictedSubtype <- attributes(OneSummary)$names[which(OneSummary == max(OneSummary))[1]]
					OneAns <- c(OnePredictedSubtype, max(OneSummary)/140)
					PredictedSummary <- rbind(PredictedSummary, OneAns)
				}
				PredictedSummary
			}
			##Classify each sample for 140 times.
			KNN_k_Classification <- c()
			for (i in AllSampleid){
				FeatureGene <- unlist(strsplit(as.vector(td$dif[[1]][i]), ','))
				GSE6477_Adjusted_TriNontriChro_Mean_Data <- Adjusted_TriNontriChro_Mean(td$exp, FeatureGene)
				Adjusted_TriNontriChro_Mean_Data <- Adjusted_TriNontriChro_Mean(ExprData_MM_scale, FeatureGene)
				TrainingData <- GSE6477_Adjusted_TriNontriChro_Mean_Data[-i, ]
				TestData <- Adjusted_TriNontriChro_Mean_Data
				TrainingDataClass <- factor(c(rep('HMM', 70), rep('NHMM', 70))[-i])
				PredictedResult <- KNN_k_Predict(TrainingData, TestData, TrainingDataClass, k)
				KNN_k_Classification <- cbind(KNN_k_Classification, as.character(PredictedResult), attr(PredictedResult, 'prob'))
			}
			KNN_k_Summary <- cbind(colnames(ExprData_MM), PredictSubtype(KNN_k_Classification))
			HMMAmount <- c('PS:', 'HMM_amount', sum(KNN_k_Summary[, 2] == 'HMM'))
			NHMMAmount <- c('PS:', 'NHMM_amount', sum(KNN_k_Summary[, 2] == 'NHMM'))
			KNN_k_Summary <- rbind(KNN_k_Summary, HMMAmount, NHMMAmount)
			colnames(KNN_k_Summary) <- c('SampleName', 'PredictedSubtype','PredictedProb' )
			rownames(KNN_k_Summary) <- NULL
			##Output the classification result.
			output <- paste(outdir,'MMSCresult_k', k, '.txt', sep = '')
			write.table(KNN_k_Summary, file = output, sep = '\t', quote = F)
			###################################
			###################################
			print('adjusted trisomy&nontrisomy chro mean distribution')
			exp_sort <- paste(outdir,'ExprData_MM_sort.txt',sep='')
			Adjusted_TriNontriChroMean_Plot_Fold <- paste(outdir,'MeanDistrPlot/',sep='')
			dir.create(Adjusted_TriNontriChroMean_Plot_Fold)
			##Sort the data in the order: NHMM, HMM and output with limit '\t'.
			HMM_Name <- KNN_k_Summary[KNN_k_Summary[, 2] == 'HMM', 1]
			NHMM_Name <- KNN_k_Summary[KNN_k_Summary[, 2] == 'NHMM', 1]
			ExprData_MM_sort <- cbind(ExprData_MM[, match(NHMM_Name, colnames(ExprData_MM))], ExprData_MM[, match(HMM_Name, colnames(ExprData_MM))])
			write.table(ExprData_MM_sort, file = exp_sort, sep = '\t', quote = F)
			##Function of ploting ‘Adjusted Trisomy&Nontrisomy Chro Mean Distribution’
			Adjustd_TriNontriChroMean_Plot <- function(GSE6477_TriNontriChro_Mean_Input, GSE6477_HMMid, GSE6477_NHMMid, TriNontriChro_Mean_Input, HMMid, NHMMid){
				GSE6477_TriChro_Mean <- as.numeric(GSE6477_TriNontriChro_Mean_Input[c(GSE6477_HMMid, GSE6477_NHMMid), 1])
				GSE6477_NontriChro_Mean <- as.numeric(GSE6477_TriNontriChro_Mean_Input[c(GSE6477_HMMid, GSE6477_NHMMid), 2])
				TriChro_Mean <- as.numeric(TriNontriChro_Mean_Input[c(HMMid, NHMMid), 1])
				NontriChro_Mean <- as.numeric(TriNontriChro_Mean_Input[c(HMMid, NHMMid), 2])
				plot(c(GSE6477_TriChro_Mean, TriChro_Mean), c(GSE6477_NontriChro_Mean, NontriChro_Mean),
					col = c(rep(rgb(1, 0.5, 1), length(GSE6477_HMMid)), rep(rgb(0.5, 0, 1), length(GSE6477_NHMMid)), rep('red', length(HMMid)), rep('blue', length(NHMMid))),
					pch = c(rep(1, length(GSE6477_HMMid)), rep(1, length(GSE6477_NHMMid)), rep(2, length(HMMid)), rep(2, length(NHMMid))), 
					xlab = 'adjusted tri-chro mean', ylab = 'adjusted nontri-chro mean', main = 'Adjusted Trisomy&Nontrisomy Chro Mean Distribution')
				legend('topright', c('GSE6477 HMM', 'GSE6477 NHMM', 'UserData HMM', 'UserData NHMM'), pch = c(1,1,2,2), col = c(rgb(1, 0.5, 1), rgb(0.5, 0, 1), 'red', 'blue'))
				abline(h = 0, v = 0, col = 'gray', lty = 2)}
			##Plot the ‘Adjusted Trisomy&Nontrisomy Chro Mean Distribution’ of each sample.
			UserData_HMMid <- match(KNN_k_Summary[(KNN_k_Summary[, 2] == 'HMM'), 1], colnames(ExprData_MM))
			UserData_NHMMid <- match(KNN_k_Summary[(KNN_k_Summary[, 2] == 'NHMM'), 1], colnames(ExprData_MM))
			for (i in AllSampleid){
				GSE6477_HMMid <- intersect(1:70, AllSampleid[-i])
				GSE6477_NHMMid <- intersect(71:140, AllSampleid[-i])
				FeatureGene <- unlist(strsplit(as.vector(td$dif[[1]][i]), ','))
				GSE6477_Adjusted_TriNontriChro_Mean_Data <- Adjusted_TriNontriChro_Mean(td$exp, FeatureGene)
				Adjusted_TriNontriChro_Mean_Data <- Adjusted_TriNontriChro_Mean(ExprData_MM_scale, FeatureGene)
				PDF_Output <- paste(Adjusted_TriNontriChroMean_Plot_Fold, 'AdjustedMeanDistri_Plot', i, '.pdf', sep = '')
				pdf(PDF_Output)
				Adjustd_TriNontriChroMean_Plot(GSE6477_Adjusted_TriNontriChro_Mean_Data, GSE6477_HMMid, GSE6477_NHMMid, Adjusted_TriNontriChro_Mean_Data, UserData_HMMid, UserData_NHMMid)
				dev.off()}
		}
	# main function
		library(affy)
		library(AnnotationDbi)
		library(class)
		print("========================================")
		print(paste('Start time:', Sys.time(), sep = ''))
		#step 1
		print('Check the data ...')
		if (Raw_Dir == "" & Exp_File == ""){
			print('no data. EXIT!')
			Exp_File_Sign=F}
		if (Raw_Dir != ""){
			Raw_Dir=switch_path(Raw_Dir)
			Out_Dir=switch_path(Out_Dir)
			dir.create(Out_Dir)
			Exp_File="Expression_Temp.XLS"
			Exp_File=paste(Out_Dir,Exp_File,sep="")
			print('raw data --> expression start:')
			print('please wait about 1min~1h depend on the size of data')
			rp(Raw_Dir,CDF_Name,Exp_File)
			print('raw --> expression finished.')
			Exp_File_Sign=T}
		if (Raw_Dir == "" & Exp_File != ""){
			print("exist expression data")
			Exp_File_Sign=T}
		#step 2
		if (Exp_File_Sign)
			cf(Exp_File,Out_Dir,K_Const,Run_Time,Tra_Data)
		print(paste('End time:', Sys.time(), sep = ''))
		print("========================================")
}
