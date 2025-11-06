## Quality controls

First you have to remove sequencing adapters and  those reads which have poor quality, before downstream analysis. To do this, you want to use TRIM_GALORE (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).

To display  the menu:

```
trim_galore --help
```

 Run TRIM_GALORE using some default parameters.

```
trim_galore --fastqc --paired --retain_unpaired R1.fastq R2.fastq 
```

If you want to look the results using strigment quality criteria you can use the `-q/—quality` option and or

 `--length_1` y `--length_2` creating a new output directory

```
trim_galore --paired --retain_unpaired --length_1 40 R1.fastq R2.fastq  --output_dir  OutDir_paired_40
```



## **Read Mapping**

Create a new work directory using `mkdir` named as **ReadMapping**  First, you need to format the reference genome with the command `bowtie2-build`. To display the menu use the option `--help` 

```
bowtie2-build --help
```

Formatting the reference genome.

```
bowtie2-build reference.fa REF
```

Mapping quality reads.

```
bowtie2 --no-unal -x REF -1 QualityReads_R1.fastq -2 QualityReads_R2.fastq -S outFile.sam
```

You can use the options  `--sensitive` and `--very sensitive` to use strigment criteria.

The file  outFile.sam contains only the mapped reads because we used the option `--no-unal`.

You can use the `grep -c` option to count the mapped reads or converter the sam file in a .bam  file (binary) to count the mapped using the -F 0x04 or -F 4 flag

```
samtools view -bSo MappedReads.bam MappedReads.sam
samtools sort MappedReads.bam -o MappedReads_SORTED.bam
samtools view -F 0x04 -c MappedReads_SORTED.bam
samtools index MappedReads_SORTED.bam
samtools view -h -o MappedReads_SORTED.sam MappedReads_SORTED.bam
```



# Read Count

An important aspect in RNA-Seq analyzes is to incorporate in a table the number of reads mapped for each gene. For this, we fisrt we have to perform a gene calling and annotation of the reference genome. This time we are going to use prokka (https://github.com/tseemann/prokka)

##### PROKKA

```
pharokka ReferenceGenome.fna
```

##### PROKKA Output Files

| Extension | Description                              |
| --------- | ---------------------------------------- |
| .gff      | This is the master annotation in GFF3 format, containing both sequences and annotations. It can be viewed directly in Artemis or IGV. |
| .gbk      | This is a standard Genbank file derived from the master .gff. If the input to prokka was a multi-FASTA, then this will be a multi-Genbank, with one record for each sequence. |
| .fna      | Nucleotide FASTA file of the input contig sequences. |
| .faa      | Protein FASTA file of the translated CDS sequences. |
| .ffn      | Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA) |
| .sqn      | An ASN1 format "Sequin" file for submission to Genbank. It needs to be edited to set the correct taxonomy, authors, related publication etc. |
| .fsa      | Nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file. It is mostly the same as the .fna file, but with extra Sequin tags in the sequence description lines. |
| .tbl      | Feature Table file, used by "tbl2asn" to create the .sqn file. |
| .err      | Unacceptable annotations - the NCBI discrepancy report. |
| .log      | Contains all the output that Prokka produced during its run. This is a record of what settings you used, even if the --quiet option was enabled. |
| .txt      | Statistics relating to the annotated features found. |
| .tsv      | Tab-separated file of all features: locus_tag,ftype,len_bp,gene,EC_number,COG,product |

##### Vizualize the reads using Artemis and analyzed the read count

Open in Artemis your reference annotated genome ( .gbk files from NCBI) and load mapped reads.

```
File -> Read BAM / VCF option
Select -> Select CDS Features ->  click on the pink region and select the option Analyzed Read Count 
```

##### Read Count using HTSeq

If you are working with a very large dataset (mouse genome) we recommend using HTSeq (https://github.com/htseq). Before use HTSeq, the prokka .gff file need to be formatted using the script `formatGFF.pl`

```
perl formatGFF.pl PROKKA.gff > new-formatted_prokka.gff
```

run `htseq-count`

```
htseq-count -m union -a 2 -s no --idattr=gene_id aln-28a-SORTED.sam new-formatted_prokka.gff -o COUNTS > Counts-28a.txt
```



# Differential gene expression analysis (DESeq2)

Open a terminal and accsess to  `R` 

use the commands `setwd` , `getwd` , `ls` and  `list.files`

`getwd` print the path of the working directory.
`setwd` move to an other directory.
`ls` list all the elements created in R.
`list.files` list all the files in the working directory.

```
getwd() 
setwd()
ls() 
list.files() 
```

Use DESeq2 from bioconductor package (http://bioconductor.org/packages/release/bioc/html/DESeq2.html). Load the libraries **DESeq2**, **gplots** and  **RColorBrewer**

> ```
> library("DESeq2")
> library("gplots")
> library("RColorBrewer")
> ```

With the function `read.csv` read your gene expression table

```
countData <- read.csv("A_B.table.txt", header=T, row.names=1, sep="\t") 
```

```
head (countData)
```

Assign the condition labels

```
colData <- DataFrame(condition=factor(c("C28","C28", "C37", "C37")))

(condition <- factor(c(rep("C28", 2), rep("C37", 2))))
```

Create the expression file

```
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
```

Run the analysis

```
dds <- DESeq(dds)
```

```
res <- results(dds)
```

You can sort the results using the adjusted *p*value 

```
res <- res[order(res$padj), ] 
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head (resdata)
```

Write and save the results as a .csv file

```
write.csv(resdata, file="diffexpr-results-deseq.csv")
```

An important aspect when analyzing the RNA-Seq data lies in knowing the dispersion of the data with respect to the experimental replicates.

```
rld <- rlogTransformation(dds)
plotPCA(rld, intgroup = "condition")

#save as a pdf file
pdf ("PCA_plot.pdf")
plotPCA(rld, intgroup = "condition")
dev.off()

```



```
plotDispEsts( dds, ylim = c(1e-6, 1e1))
```



```
plotMA(res, ylim=c(-2,2))
```



If you want to consult more aprroaches to analyzed your data you can go to: https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html



# Differential gene expression analysis (NOISeq)

Load the library

```
 library("NOISeq")

# or this way:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("NOISeq")
library("NOISeq")

```

Read the gene expression table

```
countData <- read.csv("B_A.table.txt", header=T, row.names=1, sep="\t")
head(countData)
```

Create the gene expression set

```
myfactors <- data.frame(Strain= c("C28", "C28", "C37", "C37"))
mydata <- readData(data = countData, factors = myfactors)
head(assayData(mydata)$exprs)
head(pData(mydata))
```

Estimate the  global saturation. **NOTE: This results could be note reliable if you are using the workshop dataset** 

```
mysaturation <- dat(mydata, k = 0, type = "saturation")
explo.plot(mysaturation, samples = 1:4, yleftlim = NULL, yrightlim = NULL)
```

Plot the gene expression values distribution for each sample

```
mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, samples = NULL, plottype = "boxplot")
```



Plot the proportion of low expresion genes 

```
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")
```

Normalyzed the data 

```
mycd = dat(mydata, type = "cd", norm = FALSE, refColumn = 1)
explo.plot(mycd)
```

Differential gene expression

```
mynoiseqbio = noiseqbio(mydata, k=0.5, norm="tmm", factor="Strain", lc=1, r=20, adj=1.5, plot=FALSE, a0per=0.9, random.seed=12345, filter=1)
```

Select the  DE genes

```
mynoiseq.deg= degenes(mynoiseqbio, q=0.95, M=NULL)
mynoiseq.deg= degenes(mynoiseqbio, q=0.90, M=NULL)
mynoiseq.deg= degenes(mynoiseqbio, q=0.85, M=NULL)
mynoiseq.deg= degenes(mynoiseqbio, q=0.95, M="up")
```

> "differentially expressed features (down in first condition)"

```
DE.plot(mynoiseqbio, q = 0.95, graphic = "MD")
```

Write and save the results as a .csv file

```
write.csv(mynoiseq.deg, file="diffexpr-results-noiseq.csv")
```

If you want to consult more aprroaches to analyzed your data you can go to: https://bioconductor.org/packages/devel/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
