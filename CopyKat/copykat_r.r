
library(devtools)
library(Seurat)
library(SeuratDisk)
library(copykat)
library(dplyr)


tnbc <- readRDS('/homefolder/royl/BRCA_TNBC_atlas_major_mid_minor.rds') 

# samples <- c('BIOKEY_10_Pre', 'BIOKEY_16_Pre', 'BIOKEY_15_Pre', 'BIOKEY_9_Pre', 'BIOKEY_31_Pre', 'BIOKEY_26_Pre', 'BIOKEY_11_Pre','BIOKEY_11_On', 'BIOKEY_8_Pre', 'BIOKEY_25_Pre'
# ,'CID44971','sc5rJUQ042','GSM4909282','scrJUQ070','CID44041','Pre_P022_t')
# samples <- c('BIOKEY_11_Pre','BIOKEY_25_Pre','GSM4909282','scrJUQ070','CID44041','Pre_P022_t','BIOKEY_2_Pre','BIOKEY_14_Pre','BIOKEY_1_On') first 2 already done
#04/04/24
# samples <- c('BIOKEY_8_Pre','BIOKEY_26_Pre','BIOKEY_9_Pre','BIOKEY_10_Pre', 'BIOKEY_16_Pre', 'BIOKEY_15_Pre',  'BIOKEY_31_Pre', 'BIOKEY_11_On', 'BIOKEY_2_Pre','BIOKEY_14_Pre','BIOKEY_1_Pre','CID44971','sc5rJUQ042','scrJUQ070','CID44041','Pre_P022_t','BIOKEY_1_On','GSM4909282')

# samples <- c('GSM4909282','CID44991','scrJUQ068','CID4515','CID4523','CID4513',
# 'GSM4909284','GSM4909283','sc5rJUQ058','sc5rJUQ039','sc5rJUQ045','CID4465','sc5rJUQ051',
# 'CID3946','CID4495','TBB165','sc5rJUQ026','sc5rJUQ053','sc5rJUQ033','GSM5457208','GSM4909281','GSM5457199','scrJUQ059') 5/4/24 - done

for (sample in samples) {
    
    copykat.results <- 0

    #Import of seurat object that contains raw data
    print(' ')
    print(paste('Begin run',sample,sep=' '))
    data_subset=subset(x = tnbc, SampleId == sample)

    data_subset[["Cells"]] <- Cells(data_subset)

    mainDir = "/homefolder/royl/CNV/CopyKat/TNBC/BIOKEYs"

    dir.create(file.path(mainDir, sample), showWarnings = FALSE) # subDir = sample
    setwd(file.path(mainDir, sample))

    data_subset[["reference"]] <- 'not_ref'
    data_subset$reference[data_subset$major_cell_type %in% c('B cells','Myeloid','Plasma','T cells/NK','pDC')] <- 'ref'
    #Ref by cell type

    ref_vector <- as.vector(data_subset$Cells[data_subset$reference == 'ref'])
    print(paste('ref check: ',ref_vector[1:2]))
    exp.rawdata <- as.matrix(data_subset@assays$RNA$counts)

    copykat.results <- copykat(rawmat=exp.rawdata,norm.cell.names =ref_vector, sam.name=sample, n.cores=30,ngene.chr = 3,LOW.DR = 0.01,KS.cut = 0.01, plot.genes="TRUE",)

    prediction = read.csv(paste(sample,'copykat_prediction.txt',sep='_'),sep = '\t', )
    rownames(prediction) <- prediction$cell.names
    pred_info <- prediction[match(rownames(data_subset[[]]), rownames(prediction)), "copykat.pred"]
    data_subset[["prediction"]] <- pred_info

    print(subset(x = data_subset, prediction == 'aneuploid')[[]]%>%  count(major_cell_type))
    write.table(subset(x = data_subset, prediction == 'aneuploid')[[]]%>%  count(major_cell_type), file = paste(sample, "cell_type_counts_aneuploid.txt",sep="_"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(subset(x = data_subset, prediction == 'diploid')[[]]%>%  count(major_cell_type), file = paste(sample,"cell_type_counts_diploid.txt",sep="_"), sep = "\t", row.names = FALSE, quote = FALSE)
}
# ___________________________________________________________________________________

# print('Begin run # 2')

# data_subset=subset(x = tnbc, SampleId == 'BIOKEY_2_Pre')

# data_subset[["Cells"]] <- Cells(data_subset)

# exp.rawdata <- as.matrix(data_subset@assays$RNA$counts)

# # **Run CopyKAT**

# mainDir = "/homefolder/royl/CNV/CopyKat/TNBC/BIOKEYs"
# subDir ="BIOKEY_2_Pre"

# dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
# setwd(file.path(mainDir, subDir))

# # **Define reference cells** c('CAF','B cells','Myeloid','Plasma','T cells/NK','pDC','Epithel'), infercnv score c('high CNV score','not Epithel','low CNV score','intermediate CNV score')
# data_subset[["reference"]] <- 'not_ref'
# data_subset$reference[data_subset$major_cell_type %in% c('B cells','Myeloid','Plasma','T cells/NK','pDC')] <= 'ref' #Ref by cell type
# # data_subset$reference[data_subset$cnv_status %in% c('not Epithel','low CNV score')] <- 'ref' #Ref by cell type and inferCNV score


# ref_vector <- as.vector(data_subset$Cells[data_subset$reference == 'ref'])

# samplename = subDir
# copykat.BIOKEY_2_Pre <- copykat(rawmat=exp.rawdata,norm.cell.names =ref_vector, sam.name=samplename, n.cores=16,ngene.chr = 3,LOW.DR = 0.01,KS.cut = 0.01)

# prediction = read.csv(paste(samplename,'copykat_prediction.txt',sep='_'),sep = '\t', )
# rownames(prediction) <- prediction$cell.names
# pred_info <- prediction[match(rownames(data_subset[[]]), rownames(prediction)), "copykat.pred"]
# data_subset[["prediction"]] <- pred_info

# print(subset(x = data_subset, prediction == 'aneuploid')[[]]%>%  count(major_cell_type))
# write.table(subset(x = data_subset, prediction == 'aneuploid')[[]]%>%  count(major_cell_type), file = paste(samplename, "cell_type_counts_aneuploid.txt",sep="_"), sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(subset(x = data_subset, prediction == 'diploid')[[]]%>%  count(major_cell_type), file = paste(samplename,"cell_type_counts_diploid.txt",sep="_"), sep = "\t", row.names = FALSE, quote = FALSE)
