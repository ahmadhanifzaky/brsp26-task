#R Script : Whole blood Transcriptional Profiles as a Prognostic and Diagnostic Tool in Complete and Incomplete Kawasaki Disease  
#Dataset: GSE68004 (Healthy vs Patient)
#Platform: Microarray (Illumina HumanHT-12 V4.0 - GPL10558)

#1. Install Packages
#Grouping all the library in BioManager and other package
biopackages <- c("GEOquery", "limma", "illuminaHumanv4.db", "AnnotationDbi")
packages2 <- c("pheatmap", "ggplot2", "dplyr", "umap")

#Installing Biocmanager and Packages
if (!require("BiocManager", quietly = TRUE))  {
  install.packages("BiocManager")
}
BiocManager::install(biopackages, ask = FALSE, update = FALSE)
install.packages(packages2)

#2. Import Library
sapply(biopackages, library, character.only=TRUE)
sapply(packages2, library, character.only=TRUE)

#3. Export GEO Data from GEO2R
dataset <- getGEO("GSE68004", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]] 

#4. Pre-Process Data
#Splitting the dataset to six equal quantile
exp <- exprs(dataset) 
qexp <- as.numeric(quantile(exp, c(0, 0.25, 0.5, 0.75, 0.99, 1),  na.rm = TRUE))

#Using log2 to detect and deleting NA
LogTransform <- (qexp[5] > 100) || (qexp[6] - qexp[1] > 50 && qexp[2] > 0) 
if (LogTransform) {
  exp[exp <= 0] <- NA
  exp <- log2(exp)
} 

#5 Define Sample Group
#Making group of data based on 'Source Name' Coloumn
groupinfo <- pData(dataset)[["source_name_ch1"]]
allgroups <- make.names(groupinfo)

#Classifying the data into two Category
#'whole.blood..healthy' : "Healthy"
#'whole blood, cKD', 'whole blood, HAdV', 'whole blood, inKD', 'whole blood, GAS', 'whole blood, GAS/SF'
allgroups <- ifelse (allgroups == "whole.blood..healthy", "Healthy", "Patient")
dataset$group <- factor(allgroups)
groupname <- levels(dataset$group)
print(groupname)

#6 Design Matrix
design <- model.matrix(~0 + dataset$group)
colnames(design) <- levels(dataset$group)

Group1 <- "Patient"
Group2 <- "Healthy"
contrast_formula <- paste(Group1, "-", Group2)
print(paste("Analyzed Group Contrast:", contrast_formula))

#7 Analyse Differentional Expression
fit <- lmFit(exp, design)
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design) 
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

topTableResults <- topTable(
  fit2,   adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.05
) 
head(topTableResults)

#8 Annotate Gene Name
probe_ids <- rownames(topTableResults)
gene_annotation <- AnnotationDbi::select(
  illuminaHumanv4.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")]) 

#9a Boxplot
group_colors <- as.numeric(dataset$group) 

boxplot(
  exp,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot",
  ylab = "Expression Value (log2)"
) 

legend(
  "topright",
  legend = levels(dataset$group),
  fill = unique(group_colors),
  cex = 0.8
) 

#9b Density Plot
expr_long <- data.frame (
  exp = as.vector(exp),
  Group = rep(dataset$group),
  each = nrow(exp)
 )

ggplot(expr_long, aes(x = exp, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs( 
    title = "Gen Expression  Value Distribution",
    x = "Expression Value (log2)",
    y = "Density"
  )

#9c UMAP
umap_input <- t(exp)
umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = dataset$group
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) + 
  geom_point(size = 3, alpha = 0.8) + 
  theme_minimal() + 
  labs( 
    title = "UMAP Plot",
    x = "UMAP 1",
    y = "UMAP 2"
  )

#10a Volcano Plot
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 0.1 & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -0.1 & volcano_data$adj.P.Val < 0.05] <- "DOWN" 

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "red", "NO" = "grey", "UP" = "green")) + 
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +   theme_minimal() +
  ggtitle("Volcano Plot")

#10b Heatmap
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]
top50 <- head(topTableResults, 50)

mat_heatmap <- exp[top50$PROBEID, ]

gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,      # jika SYMBOL kosong → probe ID
  top50$SYMBOL        # jika ada → gene symbol
)
rownames(mat_heatmap) <- gene_label

mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

annotation_col <- data.frame(
  Group = dataset$group
)
rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale = "row",         
  annotation_col = annotation_col,
  show_colnames = FALSE, 
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
) 

#11 Save File
write.csv(topTableResults, "testing.csv")
message("Analysis has been completed.") 
