BiocManager::install("tximport")
library('tidyverse')
library('Seurat')
library('tximport')
BiocManager::install('org.Hs.eg.db')
BiocManager::install("EnsDb.Hsapiens.v79")
library(EnsDb.Hsapiens.v79)
library(org.Hs.eg.db)
library(biomaRt)
library(RColorBrewer)


#loading the counts matrix
path <- file.path("./alevin_output/alevin/quants_mat.gz")
file<- tximport(path, type="alevin")
counts<- file$counts

#changing ensembl gene ID to gene symbol
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl",host = 'useast.ensembl.org')

gene_ids <- counts@Dimnames[[1]]
names <- getBM(attributes=c('ensembl_gene_id',
                          'hgnc_symbol'),
            filters = 'ensembl_gene_id',
            values = sub("[.][0-9]*$", "", gene_ids),
             mart = ensembl)
names=names[!duplicated(names$ensembl_gene_id),]
final_map=c()
for (row in 1:nrow(names)) {
  if(names[row,]$hgnc_symbol!='') {
    final_map <- append(final_map, names[row,]$hgnc_symbol)
  }
  else {
    final_map <- append(final_map, names[row,]$ensembl_gene_id)
  }
}

counts@Dimnames[[1]] <- final_map

#original genes and cells
num_genes <- counts@Dimnames[[1]]
num_cells <- counts@Dimnames[[2]]

# keep genes that have more than 10 nonzero counts
genes <- rowSums(counts > 0) >= 10
filtered_counts <- counts[genes, ]

#creating the seurat object
scrna <- CreateSeuratObject(counts = filtered_counts, project = "scrna", min.cells = 3, min.features = 200)
scrna

#visualizing QC metrics
scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")
head(scrna@meta.data, 5)
VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(scrna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scrna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Filter cells
scrna <- subset(scrna,subset = (nFeature_RNA > 200) & (nFeature_RNA < 3500) & (nCount_RNA > 500) & (percent.mt < 0.9))

#genes and cells after filtering
num_genes_filtered <- scrna@assays$RNA@counts@Dim[2]
num_cells_filtered <- scrna@assays$RNA@counts@Dim[1]

#normalizing the data
scrna <- NormalizeData(scrna)

#feature selection
scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scrna), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scrna)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#scaling the data
all.genes <- rownames(scrna)
scrna <- ScaleData(scrna, features = all.genes)

#Perform linear dimensional reduction
scrna <- RunPCA(scrna, features = VariableFeatures(object = scrna))
VizDimLoadings(scrna, dims = 1:2, reduction = "pca")
DimPlot(scrna, reduction = "pca")

#heatmap
DimHeatmap(scrna, dims = 1, cells = 500, balanced = TRUE)

scrna <- JackStraw(scrna, num.replicate = 100)
scrna <- ScoreJackStraw(scrna, dims = 1:20)
JackStrawPlot(scrna, dims = 1:15)

ElbowPlot(scrna)

#clustering the cells
scrna <- FindNeighbors(scrna, dims = 1:10)
scrna <- FindClusters(scrna, resolution = 0.5)

#make bar chart or pie chart
clusters <- as.data.frame(table(Idents(scrna)))

clusters$proportion <- clusters$Freq / sum(clusters$Freq) * 100 

ggplot(data=clusters, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()+
  labs(title="Bar chart for number of cells in each cluster",x="Cluster",y='Frequency')

ggplot(data=clusters, aes(x=Var1, y=proportion)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()+
  labs(title="Bar chart for relative proportion of cell numbers in each cluster",x="Cluster",y='Percentage')
