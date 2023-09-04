# Authors       : Myron G. Best & Sjors G.J.G. In 't Veld
# Email         : m.best@vumc.nl; g.intveld1@vumc.nl
# Summary       : Script to load thromboSeq processed HTSeq, perform data QCs, differential expression 
#                 of splice junction analysis, and PSO-enhanced development of thromboSeq classification
#                 algorithms
# Date          : 14th of November 2018
# Revision      : 19th of November 2018
# Revision v1.1 : 28th of November 2018
# Revision v1.2 : 29th of December 2018
# Revision v1.3 : 20th of March 2019
# Revision v1.4 : 29th of August 2019
# Revision v1.5 : 17th of April 2021
# req. packages : methods       [general functions]
#                 reshape       [visualize data in a heatmap]
#                 affy          [visualize data in a heatmap]
#                 ggplot2       [visualize data in a heatmap]
#                 dendroextras  [visualize data in a heatmap]
#                 e1071         [SVM-algorithm and classification]
#                 ROCR          [generate ROC curves and calculate AUC-values]
#                 caret         [calculate highly correlated transcripts during spliced RNA biomarker panel selection]
#                 pROC          [calculate ROC 95% confidence intervals]
#                 ppso          [run particle-swarm optimization algorithms, 
#                                included as ppsoThromboSeq.tar in the Supplemental Code (/bin/ppsoThromboSeq.tar); 
#                                install via: install.packages("ppsoThromboSeq.tar", type="source", repos=NULL)]
#                 edgeR         [store data in a data object called a DGEList, calculate normalization factor, 
#                                perform ANOVA Likelihood-ratio analysis]
#                 doSNOW/doMC   [parallel computing in order to enable efficient data analysis]
#                 foreach       [enable for repeated executions]
#                 RUVSeq        [perform RUVg data correction]
#
# note: If you do not have one of these packages installed, install by source("https://bioconductor.org/biocLite.R") and then biocLite("packageName")

# open R by running in terminal the following command: R
#rm(list = ls())
# continue in R
# load the required functions by running the following commands: 
source('bin/thromboSeqTools_PreProcessing_2.R')
source('bin/thromboSeqTools_PrePocessing_independent.R')
source('bin/thromboSeqTools_ANOVA.R')
source('bin/thromboSeqTools_PSO.R')



# sample.info.table <- read.csv("./sampleInfo1.csv", sep = ",")
# sampleInfo <- subset(dataCTEPH, dataCTEPH$discard == 0)
# write.csv(sampleInfo, file = "sampleInfo.csv")
# Collect the intron-spanning read counts outputted by the mapping pipeline 
# by running the following command:
counts <- collect.read.counts(inputDir = "HTSeq_CTEPH/")
# load("./allSamples_readCounts_ISreads.RData")
# Prepare a DGE object for subsequent data analysis.


# For this, run the following command:
dge <- prepare.dge.object(sample.info = "./SampleInfo1.csv",
                          gene.info = "bin/dgeGenesEnsembl75.RData",
                          read.counts = counts)
colnames(dge) <- paste(colnames(dge),"-",as.character(dge$samples$id),sep="")
summary(dge$samples$group)

#####################
# dge.Pre_Post_Capillary <- dge







# print('Rscript is is running...')
# 
# if(exists("./dataset/DGE_raw.Rdata")){
#   print('Load data from dge.Rdata...')
#   load("./dataset/DGE_raw.Rdata")
# }else{
#   print('loading data from countmatrix.mat and samplesInfo.csv...')
#   dge = read.mat.files(counts='./dataset/countMatrix.mat',sample.info='./dataset/sampleInfo.csv',genes='./dataset/genes.mat')
#   
#   dge <- dge[,colnames(dge)[dge$samples$discard==0]] # remove samplesnot udes in analyses before determining the plateletome
#   dge$samples <- droplevels(dge$samples)
#   summary(dge$samples$group)
#   
#   # check if columnNames for RUV are available
#   
#   # save DGE
#   save(dge, file = "./dataset/DGE_raw.Rdata")
# }
# 
# #####################
# dge <- dge[,colnames(dge)[which(dge$samples$discard == 0)]]
summary(dge$samples$group)

dge.CTEPH <- dge


predictedGroup = 'CTEPH'
controlGroup = 'HD'
predictedGroupColor = '#a8ddb5'
controlGroupColor = '#43a2ca'
trainingEvaluationColor = '#95A5A6'
validationColor = '#212F3C'
analyses.settings <- data.frame(predictedGroup,controlGroup,predictedGroupColor,controlGroupColor,trainingEvaluationColor,validationColor,stringsAsFactors=FALSE)
save(analyses.settings, file = 'bin/analysesSettings.RData')

rownames(dge.CTEPH$samples) <- colnames(dge.CTEPH)
# dge.CTEPH <- dge.CTEPH[,colnames(dge.CTEPH)[which(dge.CTEPH$samples$discard == 0)]]
dge.CTEPH$samples$isTrainingEvaluation <- (dge.CTEPH$samples$isTraining==1 | dge.CTEPH$samples$isEvaluation==1)
dge.CTEPH$samples <- droplevels(dge.CTEPH$samples)
summary(dge.CTEPH$samples$group)

# Filter the dataset for low abundant RNAs.
# For this, run the following command:
dgeFiltered <- filter.for.platelet.transcriptome.TrainEval.group(dge=dge.CTEPH, 
                                                         minimum.read.counts = 30,
                                                         minimum.prct.cohort = 90,
                                                         training.series.only = TRUE,
                                                         training.series.only.samples = colnames(dge.CTEPH)[dge.CTEPH$samples$isTrainingEvaluation==T],
                                                         groupPleteletomeSelection = TRUE,
                                                         verbose = TRUE)

# Perform thromboSeqQC pipeline.

dgeIncludedSamples <- thromboSeqQC.TrainEval(dge = dgeFiltered, 
                                            #min.number.reads.per.RNAs.detected = 0, 
                                             # min.number.total.RNAs.detected = 1500, 
                                             k.variables = 2,
                                             #variable.to.assess = c("lib.size"),
                                             #variable.threshold = c(0.8),
                                             ruvg.pvalue.threshold.group = 1e-2,
                                             ruvg.pvalue.threshold.strongest.variable = 1e-2,
                                             training.series.only = TRUE,
                                             training.series.only.samples = colnames(dgeFiltered)[dgeFiltered$samples$isTrainingEvaluation==T],
                                             leave.sample.out.threshold = 0.5, 
                                             figureDir = "figureOutputFolder/QC", 
                                             number.cores = 8, 
                                             TrainEvalBased = TRUE,
                                             analyses.settings = analyses.settings,
                                             verbose = TRUE)

# Perform thromboSeq ANOVA differential expression analysis of splice junctions
# For this, run the following command:
thromboSeq.anova <- thromboSeqANOVA(dge = dgeIncludedSamples[,colnames(dgeIncludedSamples)[dgeIncludedSamples$samples$isTrainingEvaluation==T]],
                                    k.variables = 2,
                                    #variable.to.assess = c("lib.size"),
                                    #variable.threshold = c(0.8),
                                    ruvg.pvalue.threshold.group = 1e-2,
                                    ruvg.pvalue.threshold.strongest.variable = 1e-2,
                                    training.series.only = FALSE,
                                    select.biomarker.FDR = TRUE,
                                    plot = TRUE,
                                    clinical.info.heatmap = c("group","group","group"),
                                    swarm.optimization = TRUE,
                                    #n.particles = 60,
                                    #n.iterations = 3,
                                    #iteration = NULL,
                                    figureDir = "figureOutputFolder/HeatMaps",
                                    number.cores = 8,
                                    verbose = TRUE)



#[1] "PSO-selected fdr threshold:  0.327879"
#[1] "Fisher exact test PSO-enhanced heatmap: 3e-06"

# Perform PSO-enhanced thromboSeq classifier development.
# For this, run the following command:
# 
# 
thromboPSO <- thromboSeqPSO(dge = dgeIncludedSamples, 
                            n.particles = 100,
                            n.iterations = 10,
                            number.cores = 8,
                            swarm.parameters = c("lib.size","fdr","correlatedTranscripts","rankedTranscripts"), #Added by Mo
                            swarm.boundaries = c(-0.1, 1.0, 0.00001, 1.0, 0.5, 1.0, 200, "all"),#Added by Mo
                            training.samples.provided = colnames(dgeIncludedSamples)[dgeIncludedSamples$samples$isTraining==1],#Added by Mo
                            evaluation.samples.provided = colnames(dgeIncludedSamples)[dgeIncludedSamples$samples$isEvaluation==1]#Added by Mo
)

# Validate the developed thromboSeq algorithm.
# For this, run the following command:
thromboPSOreadout <- thromboSeqPSO.readout(dge = dgeIncludedSamples,
                                           replace.counts.validation = 0, 
                                           filter.clinical.characteristics.validation = NA,
                                           filter.clinical.characteristics.group = NA,
                                           filter.clinical.characteristics.specified = NA,
                                           readout.training = T, 
                                           readout.evaluation = T, 
                                           readout.validation = T, 
                                           apply.rule.in.readout = F, 
                                           rule.in.setting = NA, 
                                           apply.rule.out.readout = F, 
                                           rule.out.setting = NA,
                                           additional.dge.to.predict = NA,
                                           number.cores = 8, 
                                           clinical.data.in.output = c(as.character(colnames(dgeIncludedSamples$samples)))
)

# Run LOOCV analysis on all samples
thromboSeqLOOCV <- thromboSeq.LOOCV(dge = dge)

# Perform control experiments for thromboSeq classifier development.
# For this, run the following command:
control <- thromboSeqPSO.controls(dge = dgeIncludedSamples,
                                  filter.clinical.characteristics.validation = NA,
                                  filter.clinical.characteristics.group = NA,
                                  filter.clinical.characteristics.specified = NA,
                                  thromboSeqPSO.shuffled = T,
                                  thromboSeqPSO.iterations = T,
                                  number.cores = 8
)

# Perform digitalSWARM analysis.
# For this, run the following command:

# 
# ####
# digitalSWARM.output <- digitalSWARM( dge = dgeIncludedSamples,
#                                      variable.to.assess = c("lib.size"),
#                                      variable.to.assess.thresholds = c(0.0, 1.0),
#                                      select.biomarker.FDR = FALSE,
#                                      FDR.value.threshold = 0.1, 
#                                      percentage.for.training = 30,
#                                      percentage.for.evaluation = 30,
#                                      percentage.for.verification = 20,
#                                      percentage.for.validation = 20,
#                                      training.samples.provided = NULL,
#                                      evaluation.samples.provided = NULL,
#                                      verification.samples.provided = NULL,
#                                      validation.samples.provided = NULL,
#                                      n.training.iterations = 200,
#                                      minimum.n.transcripts.biomarkerpanel = 2,
#                                      n.particles = 100,
#                                      n.iterations = 50,
#                                      number.cores = 8, 
#                                      verbose = TRUE
# )
# ####

# Perform control experiments for digitalSWARM analysis.
# For this, run the following command:

# ####
# digitalSWARMshuffled <- digitalSWARM.shuffled(dge = dgeIncludedSamples,
#                                               FDR.value.threshold = 0.1,
#                                               variable.to.assess = c("lib.size"),
#                                               variable.to.assess.thresholds = c(0.0, 1.0),
#                                               select.biomarker.FDR = FALSE,
#                                               percentage.for.training = 30,
#                                               percentage.for.evaluation = 30,
#                                               percentage.for.verification = 20,
#                                               percentage.for.validation = 20,
#                                               training.samples.provided = NULL,
#                                               evaluation.samples.provided = NULL,
#                                               verification.samples.provided = NULL,
#                                               validation.samples.provided = NULL,
#                                               n.training.iterations = 200,
#                                               minimum.n.transcripts.biomarkerpanel = 2,
#                                               n.particles = 100,
#                                               n.iterations = 50,
#                                               n.shuffled = 100,
#                                               number.cores = 8, 
#                                               verbose = TRUE,
#                                               matrix = digitalSWARM.output)
# #####



# Optionally save and store the R-session
sessionInfo <- sessionInfo()
save.image(paste(timestamp(),".RData",sep=""))
# End

# thromboPSO <- thromboSeqPSO(dge = dgeIncludedSamples, 
#                             training.samples.provided = colnames(dgeIncludedSamples)[dgeIncludedSamples$samples$isTraining==1],
#                             evaluation.samples.provided = colnames(dgeIncludedSamples)[dgeIncludedSamples$samples$isEvaluation==1],
#                             swarm.parameters = c("lib.size","fdr","correlatedTranscripts","rankedTranscripts"),
#                             swarm.boundaries = c(-0.1, 1.0, 0.00001, 1.0, 0.5, 1.0, 200, "all"),
#                             k.variables = 2,
#                             #variable.to.assess = c("lib.size"),
#                             #variable.threshold = c(0.8),
#                             #ruvg.pvalue.threshold.group = 1e-8,
#                             #ruvg.pvalue.threshold.strongest.variable = 1e-2,
#                             #select.biomarker.FDR = TRUE,
#                             #minimum.n.transcripts.biomarkerpanel = 2,
#                             #svm.gamma.range = 2^(-20:0),
#                             #svm.cost.range = 2^(0:20),
#                             #number.cross.splits = 2,
#                             #n.particles.gamma.cost.optimization = 50,
#                             #n.iterations.gamma.cost.optimization = 4,
#                             #n.particles = 60,
#                             #n.iterations = 8,
#                             figureDir = "figureOutputFolder/PSO",
#                             number.cores = 8, 
#                             verbose = TRUE
#                             )
# 
# # Validate the developed thromboSeq algorithm.
# # For this, run the following command:
# thromboPSOreadout <- thromboSeqPSO.readout(dge = dgeIncludedSamples,
#                                            replace.counts.validation = 0, 
#                                            filter.clinical.characteristics.validation = NA,
#                                            filter.clinical.characteristics.group = NA,
#                                            filter.clinical.characteristics.specified = NA,
#                                            readout.training = T, 
#                                            readout.evaluation = T, 
#                                            readout.validation = T, 
#                                            apply.rule.in.readout = F, 
#                                            rule.in.setting = NA, 
#                                            apply.rule.out.readout = F, 
#                                            rule.out.setting = NA,
#                                            additional.dge.to.predict = NA,
#                                            number.cores = 8, 
#                                            clinical.data.in.output = NA
#                                            )
# 
# # Run LOOCV analysis on all samples
# thromboSeqLOOCV <- thromboSeq.LOOCV(dge = dge)
# 
# # Perform control experiments for thromboSeq classifier development.
# # For this, run the following command:
# control <- thromboSeqPSO.controls(dge = dgeIncludedSamples,
#                                   filter.clinical.characteristics.validation = NA,
#                                   filter.clinical.characteristics.group = NA,
#                                   filter.clinical.characteristics.specified = NA,
#                                   thromboSeqPSO.shuffled = T,
#                                   thromboSeqPSO.iterations = T,
#                                   number.cores = 8
#                                   )
# 
# # Perform digitalSWARM analysis.
# # For this, run the following command:
# digitalSWARM.output <- digitalSWARM( dge = dgeIncludedSamples,
#                                      variable.to.assess = c("lib.size"),
#                                      variable.to.assess.thresholds = c(0.0, 1.0),
#                                      select.biomarker.FDR = FALSE,
#                                      FDR.value.threshold = 0.1, 
#                                      percentage.for.training = 30,
#                                      percentage.for.evaluation = 30,
#                                      percentage.for.verification = 20,
#                                      percentage.for.validation = 20,
#                                      training.samples.provided = NULL,
#                                      evaluation.samples.provided = NULL,
#                                      verification.samples.provided = NULL,
#                                      validation.samples.provided = NULL,
#                                      n.training.iterations = 200,
#                                      minimum.n.transcripts.biomarkerpanel = 2,
#                                      n.particles = 100,
#                                      n.iterations = 50,
#                                      number.cores = 8, 
#                                      verbose = TRUE
# )
# 
# # Perform control experiments for digitalSWARM analysis.
# # For this, run the following command:
# digitalSWARMshuffled <- digitalSWARM.shuffled(dge = dgeIncludedSamples,
#                                               FDR.value.threshold = 0.1,
#                                               variable.to.assess = c("lib.size"),
#                                               variable.to.assess.thresholds = c(0.0, 1.0),
#                                               select.biomarker.FDR = FALSE,
#                                               percentage.for.training = 30,
#                                               percentage.for.evaluation = 30,
#                                               percentage.for.verification = 20,
#                                               percentage.for.validation = 20,
#                                               training.samples.provided = NULL,
#                                               evaluation.samples.provided = NULL,
#                                               verification.samples.provided = NULL,
#                                               validation.samples.provided = NULL,
#                                               n.training.iterations = 200,
#                                               minimum.n.transcripts.biomarkerpanel = 2,
#                                               n.particles = 100,
#                                               n.iterations = 50,
#                                               n.shuffled = 100,
#                                               number.cores = 8, 
#                                               verbose = TRUE,
#                                               matrix = digitalSWARM.output)
# 
# # Optionally save and store the R-session
# sessionInfo <- sessionInfo()
# save.image(paste(timestamp(),".RData",sep=""))
# # End