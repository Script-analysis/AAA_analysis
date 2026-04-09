library(Seurat)
library(tidyr)
library(ggplot2)
library(scCustomize) 
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(tidyverse)
load("scRNA_annotation.RData")
features<-"FKBP11"
plot1<-FeaturePlot(scRNA, features =features, pt.size=0.5, reduction="umap",cols =rev(brewer.pal(n = 10, name = "RdBu")))
plot1
ggsave("umap_featureplot.pdf", plot = plot1, width = 5, height = 4.5)
plot2<-FeaturePlot(scRNA, features =features, pt.size=0.5, reduction="umap",cols =rev(brewer.pal(n = 10, name = "RdBu")),split.by = "group")
plot2
ggsave("umap_featureplot_split_by_group.pdf", plot = plot2, width = 10, height = 4.5)

