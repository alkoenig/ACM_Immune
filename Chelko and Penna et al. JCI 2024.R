library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)
library(harmony)
library(ArchR)
library(ggsci)
library(ggpubr)
library(pheatmap)
library(Matrix)
library(RColorBrewer)
library(scales)
library(Nebulosa)
library(SeuratDisk)
library(rmarkdown)

#Following code is outline for processing raw filtered matrix files
#Input data
wt_data <- Read10X( data.dir = "filtered_gene_wt2")
wt <- CreateSeuratObject(counts = wt_data[["Gene Expression"]])
wt[["ADT"]] <- CreateAssayObject(counts= wt_data[["Antibody Capture"]])
wt$phenotype <- "wt"

dsg_data <- Read10X(data.dir = "filtered gene_dsg2")
dsg <- CreateSeuratObject(counts= dsg_data[["Gene Expression"]])
dsg[["ADT"]] <- CreateAssayObject(counts= dsg_data[["Antibody Capture"]])
dsg$phenotype <- "dsg"

ccr2_data <- Read10X(data.dir = "filtered gene_ccr2-2")
ccr2 <- CreateSeuratObject(counts= ccr2_data[["Gene Expression"]])
ccr2[["ADT"]] <- CreateAssayObject(counts=ccr2_data[["Antibody Capture"]])
ccr2$phenotype <- "ccr2"

#merge samples
sample <- merge(wt, y= c(dsg, ccr2))

#Following code is for inputing single nuclei seq files and scrublet filtering
wt6wk_data <- Read10X(data.dir = "filtered gene WT 6wk")
wt6wk <- CreateSeuratObject(counts = wt6wk_data)
wt6wk$genotype <- "wt6wk"

wt16wk_data <- Read10X(data.dir = "filtered gene WT 16wk")
wt16wk <- CreateSeuratObject(counts = wt16wk_data)
wt16wk$genotype <- "wt16wk"

dsg6wk_data <- Read10X(data.dir = "filtered gene dsg 6wk")
dsg6wk <- CreateSeuratObject(counts = dsg6wk_data)
dsg6wk$genotype <- "dsg6wk"

dsg16wk_data <- Read10X(data.dir = "filtered gene dsg 16wk")
dsg16wk <- CreateSeuratObject(counts = dsg16wk_data)
dsg16wk$genotype <- "dsg16wk"

ccr2_6wk_data <- Read10X(data.dir = "filtered gene dsgxccr2 6wk")
ccr2_6wk <- CreateSeuratObject(counts = ccr2_6wk_data)
ccr2_6wk$genotype <- "ccr2_6wk"

ccr2_16wk_data <- Read10X(data.dir = "filtered gene dsgxccr2 16wk")
ccr2_16wk <- CreateSeuratObject(counts = ccr2_16wk_data)
ccr2_16wk$genotype <- "ccr2_16wk"

csf2_6wk_data <- Read10X(data.dir = "filtered gene dsgxcsf2 6wk")
csf2_6wk <- CreateSeuratObject(counts = csf2_6wk_data)
csf2_6wk$genotype <- "csf2_6wk"

csf2_16wk_data <- Read10X(data.dir = "filtered gene dsgxcsf2 16wk")
csf2_16wk <- CreateSeuratObject(counts = csf2_16wk_data)
csf2_16wk$genotype <- "csf2_16wk"

IKB_16wk_data <- Read10X(data.dir = "filtered gene dsgxIKB 16wk")
IKB_16wk <- CreateSeuratObject(counts = IKB_16wk_data)
IKB_16wk$genotype <- "IKB_16wk"

sample_singlenuc <- merge(wt6wk, y= c(wt16wk, dsg6wk, dsg16wk, ccr2_6wk, ccr2_16wk, 
                                         csf2_6wk, csf2_16wk, IKB_16wk))

#scrublet
#make 5had file for scrublet
SaveH5Seurat(sample_singlenuc, filename= "sample_singlenuc.h5Seurat")
Convert("sample_singlenuc.h5Seurat", dest = "h5ad")

#Load scrublet data and remove doublets
scrub = read.csv("all.csv", header = T, row.names = 1)
sample_singlenuc@meta.data$scrublet_score = scrub$scrublet_score5
sample_singlenuc@meta.data$scrublet_cluster_score = scrub$scrublet_cluster_score5
sample_singlenuc@meta.data$bh_pval = scrub$bh_pval5

VlnPlot(sample_singlenuc, group.by = "genotype", features = "scrublet_score", 
        pt.size = 0) + NoLegend()

sample_singlenuc <- subset(
  x = sample_singlenuc,
  subset = scrublet_score < 0.25)

#Following code applies to both single nuc and CITE-seq unless otherwise stated

#QC filtering
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^mt-")
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                             "nCount_ADT"), ncol = 4, group.by = "phenotype")
sample <- subset(sample, subset = nFeature_RNA > 200 & 
                   nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 25)
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                             "nCount_ADT"), ncol = 4, group.by = "phenotype")

#normalization
DefaultAssay(sample) <- 'RNA'
sample <- SCTransform(sample, vars.to.regress = 
                        c("percent.mt", "nCount_RNA"))
sample <- RunPCA(sample, 
                 features = VariableFeatures(object = sample), 
                 npcs=100, verbose=TRUE)

#Make umap and gene list for RNA data
DefaultAssay(sample) <- "SCT"
sample <- FindNeighbors(sample, dims = 1:15, verbose = FALSE, 
                        reduction = "pca")
sample <- FindClusters(sample, resolution = 
                         c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), 
                       verbose = FALSE)
sample <- RunUMAP(sample, dims = 1:15, verbose = FALSE, reduction = "pca")

#Make umap and list for protein data (CITE-seq only)
DefaultAssay(sample) <- 'ADT'
VariableFeatures(sample) <- rownames(sample[["ADT"]])
sample <- NormalizeData(sample, normalization.method = 'CLR', 
                        margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca', npcs=100)

#Wknn analsyis, multimodal, combining RNA and protein data (CITE-seq only)
DefaultAssay(sample) <-"SCT"
sample <- FindMultiModalNeighbors(
  sample, reduction.list = list("pca", "apca"), 
  dims.list = list(1:15, 1:10), modality.weight.name = "RNA.weight"
)

sample <- RunUMAP(sample, nn.name = "weighted.nn", 
                  reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

sample <- FindClusters(sample, graph.name = "wsnn", algorithm = 3, 
                       resolution = c(0.1, 0.2, 0.3, 0.4), verbose = FALSE)

#Generating DEG list for clusters (same method used for generating all DEG lists, just change the ident)
Idents(sample) <- "wsnn_res.0.1"
sample.rnamarkers <- FindAllMarkers(sample, 
                                    only.pos = TRUE, 
                                    min.pct = 0.1, logfc.threshold = 0.25)
write.csv(sample.rnamarkers, file ="combined_data2.csv", quote = FALSE)

#Labeling groups
fun <- function(x) {
  if (x == "0") {"Fibroblast"} 
  else if (x == "1") {"Fibroblast"}
  else if (x == "2") {"B-cell"}
  else if (x == "3") {"Macrophage/monocyte"}
  else if (x == "4") {"Neutrophil"}
  else if (x == "5") {"T-cell"}
  else if (x == "6") {"Endothelial"}
  else if (x == "7") {"B-cell"}
  else if (x == "8") {"NK cell"}}

sample$annotations <- mapply(fun, sample$wsnn_res.0.1)

#subsetting cell types from global object (example of one cell type)
Idents(sample) <- "annotations"
fibroblast <- subset(sample, idents = c("Fibroblast"))

#Following code is used to generate various plots used for figures

#Violin plots
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                             "nCount_ADT"), ncol = 4, group.by = "phenotype")

#umaps
DimPlot(sample, group.by = "annotations")

#composition/stack plots
ggplot(sample@meta.data, aes(x=condition, fill=annotations)) + geom_bar(position = "fill") + theme_linedraw() + theme(axis.text.x = element_text(angle = 90)) +  scale_fill_manual(values=as.vector(paletteDiscrete(unique(sample$annotations), set = "stallion"))) + theme(axis.line = element_line(colour = "black"),
                                                                                                                                                                                                                                                                                        panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                        panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                                                                        panel.background = element_blank()) 



#Z-scores and umap plotting of z-scores
DefaultAssay(sample) <- "SCT"
expdata <- GetAssayData(sample)
Pop1 <- c(#enter genes of interest)
pops<-list(Pop1)
#Z-Scores
z_scores<-NULL

for (i in 1:length(pops)) {
  genes <- pops[[i]]
  zz <- which(tolower(rownames(expdata)) %in% tolower(genes))
  av <- numeric(ncol(expdata))
  
  geneExp <- as.matrix(expdata[zz, ])
  geneExp <- t(scale(t(geneExp)))
  geneExp[is.nan(geneExp)] <- 0
  z_scores <- rbind(z_scores,(av + colSums(geneExp) / length(zz)))
}
sample@meta.data$z_score<-z_scores[1,]
FeaturePlot(sample, features = "z_score", pt.size = 1, reduction = "umap") +
  theme(legend.position ='right')+
  theme(legend.box.margin = margin(6, 6, 6, 6))+
  theme(legend.key.size = unit(.8, "cm"))+
  theme(axis.ticks.y.right= element_blank())+
  coord_fixed() & scale_color_gradientn(colors=c("blue","turquoise2","yellow","red","red4"), oob=scales::squish, limits=c(0,1.5)) & NoAxes()

#heatmap generation
sample_avg <- AverageExpression(sample, assays = "SCT",
sample_avgmat <- GetAssayData(sample_avg, assay = "SCT")
pheatmap(sample_avgmat, scale = "row", cluster_rows = T, cluster_cols = F)

#constructing DEG dotplots

#Finding DEG in each cell type between conditions (one cell type used as example, rest would follow similarly)
cardiomyocyte <- readRDS("subcard_new.rds")
Idents(cardiomyocyte) <- "condition"
card_wt_dsg <- subset(cardiomyocyte, idents= c("wt", "dsg"))
card_wt_dsg.rnamarkers <- FindAllMarkers(card_wt_dsg, 
                                         only.pos = TRUE, 
                                         min.pct = 0.1, logfc.threshold = 0.25)
write.csv(card_wt_dsg.rnamarkers, file ="card_wt_dsg data.csv", quote = FALSE)

#dotplot for one cell type DEG
cardio <- read_delim("card_wt_dsg data.csv", ",", escape_double = FALSE, trim_ws = TRUE)
cardio$cell <- "Cardiomyocyte"
cardio$sigpvalue <- ifelse(cardio$p_val_adj < 0.05, "p < 0.05","p > 0.05")
cardio$sig <- ifelse(cardio$p_val_adj < 0.05 & abs(cardio$avg_log2FC) > 0, "Significant","Not Significant")

#combine all cell type dotplots into one graph
data <- data.frame(rbind(BCells,cardio,fibro,mac,TCells,endo,epi,peri,fat))
data$cell <- factor(data$cell, levels = c("BCells","Cardiomyocyte", "Fibroblast","Myeloid",
                                                "TCells","Endothelial","Epicardium",
                                                "Pericyte/SMC","Adipocyte"))



df_Count <- data %>% group_by(sig, cell) %>% dplyr::count()
df_Count <- data.frame(df_Count)

x <- df_Count[with(df_Count,order(n,decreasing = T)) ,][df_Count[with(df_Count,order(n, decreasing = T)) ,]$sig=="Significant",]$cell
df_Count$cell <- factor(df_Count$cell, levels = x)

data$cell <- factor(data$cell, levels = x)
data %>%
  ggplot(aes(x=cell, y=avg_log2FC, fill=cell, color=sig)) +
  geom_jitter(size=1, alpha=0.5, position=position_jitter(0.2)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="none") + scale_y_continuous(limits = c(0, 3)) +
  ggtitle("WT vs Dsg2 DEG") +
  xlab("Cell types") +
  scale_shape_manual(values=c(1,1))+
  scale_color_manual(values=c("blue", "red"))






