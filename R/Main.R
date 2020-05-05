# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl , Shift , B'
#   Check Package:             'Ctrl , Shift , E'
#   Test Package:              'Ctrl , Shift , T'

hello <- function() {
  print("Hello, world!")
}

OpenData <- function(dir="data/",
                     project.name="project",
                     min.cels=5,
                     min.feat=100,
                     outdir="output/"){

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  data <- Read10X(data.dir=dir)

  data.object<-CreateSeuratObject(counts=data,
                                  project=project.name,
                                  min.cells=min.cels,
                                  min.features=min.feat)
  data.object[["percent.mt"]]<-PercentageFeatureSet(data.object, pattern = "^MT-")

  Vplot<-VlnPlot(data.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

  plot1 <- FeatureScatter(data.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(data.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- plot1 + plot2

  png(filename=paste0(outdir,"FeatViolinPlot.png"))
  print(Vplot)
  dev.off()

  png(filename=paste0(outdir,"FeatScatterPlot.png"))
  print(plot3)
  dev.off()


  return(data.object)
}


NormalizeAndScale <- function(data.object,
                              nfeatures=500,
                              outdir="output/"){

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  data.object <- NormalizeData(data.object,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000)

  data.object <- FindVariableFeatures(data.object,
                                      selection.method = "vst",
                                      nfeatures = nfeatures)

  all.genes <- rownames(data.object)
  data.object <- ScaleData(data.object, features = all.genes)


  top20<-head(VariableFeatures(data.object),20)
  plot1 <- VariableFeaturePlot(data.object)
  plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)

  png(filename=paste0( outdir , "VarFeatPlot.png"))
  print(plot2)
  dev.off()

  return(data.object)
}



LinearAnalysis <- function (data.object,
                            dims=20,
                            cells=1000,
                            balance=TRUE,
                            JackStrawReplicates=100,
                            outdir="output/"){


  data.object <- RunPCA(data.object, features = VariableFeatures(object = data.object))
  data.object <- JackStraw(data.object,num.replicate = JackStrawReplicates)
  data.object <- ScoreJackStraw(data.object, dims = 1:dims)

  heatm <- DimHeatmap(data.object, dims = 1:dims, cells = 1000, balanced = TRUE)
  vizdimplot <- VizDimLoadings(data.object, dims = 1:4, reduction = "pca")
  dimplot <- DimPlot(data.object, reduction = "pca")
  JSPlot <- JackStrawPlot(data.object, dims = 1:dims)
  elbow <- ElbowPlot(data.object)

  png(filename=paste0(outdir,"VizdimPlot4dimsPCA.png"))
  print(vizdimplot)
  dev.off()

  png(filename=paste0(outdir,"DimplotPCA.png"))
  print(dimplot)
  dev.off()

  png(filename=paste0(outdir,"JackStrawPlot",dims,"dims.png"))
  print(JSPlot)
  dev.off()

  png(filename=paste0(outdir,"ElbowPlotPCA.png"))
  print(elbow)
  dev.off()


  return(data.object)
}


Cluster <-  function (data.object,
                      dims=5,
                      UMAPres=0.1,
                      tSNEREs=0.1,
                      outdir="output/"){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  data.object <- FindNeighbors(data.object,dims = 1:dims)

  if(UMAPres==tSNEREs){
    data.object <- FindClusters(data.object, resolution = UMAPres)

    data.object <- RunUMAP(data.object, dims = 1:dims)
    UMAPplot <- DimPlot(data.object, reduction = "umap")

    data.object <- RunTSNE(data.object, dims = 1:dims)
    tSNEPlot <- DimPlot(data.object, reduction = "tsne")
  }
  else{
    data.object <- FindClusters(data.object, resolution = UMAPres)
    data.object <- RunUMAP(data.object, dims = 1:dims)
    UMAPplot <- DimPlot(data.object, reduction = "umap")

    data.object <- FindClusters(data.object, resolution = tSNEREs)
    data.object <- RunTSNE(data.object, dims = 1:dims)
    tSNEPlot <- DimPlot(data.object, reduction = "tsne")

  }

  png(filename=paste0(outdir,"UMAPPlot_",dims,"dims_",UMAPres,"Res.png"))
  print(UMAPplot)
  dev.off()

  png(filename=paste0(outdir,"tSNEplot_",dims,"dims_",tSNEREs,"Res.png"))
  print(tSNEPlot)
  dev.off()



  return(data.object)

}





