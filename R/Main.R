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
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

OpenData <- function(dir="data/",
                     project.name="project",
                     min.cels=5,
                     min.feat=100,
                     outdir="output/"){

  dir.create("./output", showWarnings = FALSE)

  data <- Read10X(data.dir=dir)

  data.object<-CreateSeuratObject(counts=data,
                                  project=project.name,
                                  min.cells=min.cels,
                                  min.features=min.feat)
  data.object[["percent.mt"]]<-PercentageFeatureSet(day16, pattern = "^MT-")

  Vplot<-VlnPlot(day16, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

  plot1 <- FeatureScatter(data.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(data.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- plot1 + plot2

  png(filename = outdir+"FeatViolinPlot.png")
  Vplot
  dev.off()

  png(filename = outdir+"FeatScatterPlot.png")
  plot3
  dev.off()

  par(mfrow=c(2,1))
  Vplot
  plot3
  par(mfrow=c(1,1))

  return(data.object)
}


NormalizeAndScale <- function(data.object,
                              nfeatures=500,
                              outdir="output/"){

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

  png(filename = outdir + "VarFeatPlot.png")
  plot2
  dev.off()

  plot2


  return(data.object)
}


