# **Phylogenomics and pangenome analysis**



*<u>Written by  Gamaliel López-Leal</u>*

*<u>Mail: gamaliel.lopez@docentes.uaem.edu.mx</u>*



### Get homologous groups and usnig SGF approach

In this protocol we start from genome sequence (.fna files), first we run `prokka` to annotate the selected genomes, for this we used this command line:

```
prokka infile.fna --kingdom Bacteria --outdir $file.dir --locustag $file
```

If you want to run prokka in a single process. I have written a script that calls all the fna files and runs prokka for each of the files.

```
nohup perl Running_PROKKA_toMultipleGenomes.pl &
```



### Orthologous Single Gene Families (using roary)

In order to get the SGF we used Roary with the following parameters.

```
nohup roary  -b "blastp -qcov_hsp_perc 90" -f Roary-Results_80Identity-90Coverage  -i 80 -t 4  -s  -e -n -v ../*.gff -p 10 &
```

If you want to use the gene presence-absence matrix to creat multiple files based on the HG information. First, you need to create a files using the information in the presence-absence matrix file and adding a unique tag in the first column. The file should be like this:



Additionally, we perform several plots in R according to the Roary manual (using the Roary out files):

```R
library(ggplot2)


mydata = read.table("number_of_new_genes.Rtab")
png(filename = "number_new_genes.png", width = 6, height = 4, units = "in", res = 300)
boxplot(mydata, data=mydata, main="Number of new genes",
         xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(mydata)), outline=FALSE)
dev.off()

mydata = read.table("number_of_conserved_genes.Rtab")
png(filename = "number_of_conserved_genes.png", width = 6, height = 4, units = "in", res = 300)
boxplot(mydata, data=mydata, main="Number of conserved genes",
          xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(mydata)), outline=FALSE)
dev.off()


######################

unique_genes = colMeans(read.table("number_of_unique_genes.Rtab"))
new_genes = colMeans(read.table("number_of_new_genes.Rtab"))

genes = data.frame( genes_to_genomes = c(unique_genes,new_genes),
                    genomes = c(c(1:length(unique_genes)),c(1:length(unique_genes))),
                    Key = c(rep("Unique genes",length(unique_genes)), rep("New genes",length(new_genes))) )
                    
ggplot(data = genes, aes(x = genomes, y = genes_to_genomes, group = Key, linetype=Key)) +geom_line()+
theme_classic() +
ylim(c(1,max(unique_genes)))+
xlim(c(1,length(unique_genes)))+
xlab("No. of genomes") +
ylab("No. of genes")+ theme_bw(base_size = 16) +  theme(legend.justification=c(1,1),legend.position=c(1,1))
ggsave(filename="unique_vs_new_genes.png", scale=1)
```





































