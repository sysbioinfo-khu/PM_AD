setwd("/Figure_4/Conserved")


immune.combined <- readRDS("/RDS/Conserved_all.RDS")

Idents(immune.combined)<-"Status"
immune.combined<-subset(immune.combined, idents = c("Placenta_early","Placenta_term","Fetal_skin", "Eczema"))

Idents(immune.combined)<-"Annotation"
DefaultAssay(immune.combined) <- "RNA"

immune.combined <- NormalizeData(immune.combined, normalization.method = "LogNormalize", scale.factor = 10000)

macro_markers <- FindConservedMarkers(immune.combined, ident.1 = "Macro_2", 
                                      grouping.var = "Status", verbose = FALSE, test.use="MAST", min.pct = 0.5, logfc.threshold = 0.25)



write.csv(macro_markers, "conserved_marker.csv")


