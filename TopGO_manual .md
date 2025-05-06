## **GO enrichment using TopGO** 

Written by Gamaliel López-Leal

gamaliel.lopez@docentes.uaem.edu.com

```
 library ("topGO")
 library ("Rgraphviz")
```

Read your gene background or GO universe

```
geneID2GO <- readMappings(file = "MG_GO_universe.txt" , sep = "\t", IDsep = ",")

geneNames <- names(geneID2GO)

geneUniverse <- names(geneID2GO)
```

Read your gene list of interest

```
genesOfInterest <- read.table("Chicken_Chicken-Breeders.txt",header=FALSE)

genesOfInterest <- as.character(genesOfInterest$V1)

geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))

names(geneList) <- geneUniverse
```

Build the TopGO annotation set

```
myGOdata <- new("topGOdata", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
```

Run the GO enrichment

```
resultTopgo <- runTest(myGOdata, algorithm="classic", statistic="fisher")
```

Make gene table with the GO enrichment results

```
allResTopGO2 <- GenTable(myGOdata, topgoFisher = resultTopgo, ranksOf = "classicFisher", topNodes = 100)
```

Plot the results

```
showSigOfNodes(myGOdata, score(resultTopgo), firstSigNodes = 5, useInfo = 'all')
```















