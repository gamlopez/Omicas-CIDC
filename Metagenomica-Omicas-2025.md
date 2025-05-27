#### Written by Gamaliel López Leal

gamaliel.lopez@docentes.uaem.edu.mx	



### Clasificación Taxónomica con KRAKEN

Ahora usaremos Kraken para clasificar nuestros reads. Cuando usas Karken, debes de formatear la o las bases de datos de interes. Para esta páctica usaremos la  miniKraken de 4GB que se encuentra en: <u>/home/krakenDB/minikraken_20171013_4GB/</u>. Para mayor información consulta: https://ccb.jhu.edu/software/kraken/MANUAL.html

Corramos Kraken usando un archivo fastq

```
kraken --db /home/krakenDB/minikraken_20171013_4GB/ --fastq-input R1.fastq --unclassified-out UNCLASS_Reads-R1  --classified-out CLASS_Reads-R1 --output out_dirKraken_R1
```

También puedes usar el flag `—paired` para indicar que los archivos de entrada son el R1 y R2.

```
kraken --db /home/krakenDB/minikraken_20171013_4GB/ --fastq-input --paired R1.fastq R2.fastq --unclassified-out UNCLASS_Reads  --classified-out CLASS_Reads --output out_dirKraken
```

Para obtener la clasificación tenemos que usar `kraken-traslate`

```
kraken-translate --db /home/krakenDB/minikraken_20171013_4GB/ out_dirKraken > sequences.labels
```

Reporte de general de Kraken:

```
kraken-report --db /home/krakenDB/minikraken_20171013_4GB/ out_dirKraken
```

Contenido en las columnas:

1. Percentage of reads covered by the clade rooted at this taxon
2. Number of reads covered by the clade rooted at this taxon
3. Number of reads assigned directly to this taxon
4. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily,(G)enus, or (S)pecies. All other ranks are simply '-'.
5. NCBI taxonomy ID
6. indented scientific name

El módulo kraken-traslate también tiene la opción de generar un archivo de salida con formato tipo mpa, con el flag `--mpa-format`, que sólo informará de los niveles de la taxonomía con asignaciones de rango estándar (superreino, reino, filo, clase, orden, familia, género, especie).		

```
kraken-translate --mpa-format --db /home/krakenDB/minikraken_20171013_4GB/ out_dirKraken > sequences.labels.mpa		
```

Solo el reporte mpa sin información de la lectura.

```
kraken-mpa-report --db /home/krakenDB/minikraken_20171013_4GB/ out_dirKraken > Kraken_Abundance_mpa.txt
```

Seleccionar las abundancias de acuerdo al nivel taxonómico, es este caso será a nivel clase.

```
grep -E "(c__)|(^ID)" Kraken_Abundance_mpa.txt | grep -v "t__" | sed 's/^.*c__//g' > Kraken_Abundance_mpa_class.txt
```

​				
Usando LCA y valores de corte para Kraken

```
kraken-filter --db /home/krakenDB/minikraken_20171013_4GB/ [--threshold 0.1]  out_dirKraken > LCA_kraken
```



​	

### Clasificación Taxónomica con MetaPhlAn2

Usaremos Metaphla2 para clasificar por taxonomia nuestras lecturas. Para mayor información consulta: https://github.com/biobakery/MetaPhlAn2

```
metaphlan2.py R1.fastq  --input_type fastq > metaphlan2_profile1.txt
metaphlan2.py R2.fastq --input_type fastq > metaphlan2_profile2.txt
```

Pra concatenar varias muestras debes usar el script `merge_metaphlan_tables.py`

```
merge_metaphlan_tables.py metaphlan2_profile*.txt > merged_abundance_table.txt
```

Para obtener el reporte del número de lecturas clasificadas

```
metaphlan2.py SRR5689219.fastq  -t rel_ab_w_read_stats  --input_type fastq > metaphlan2_profileB_nreads.txt
```

Al igual que con Kraken, pudes seleccionar las abundancias de acuerdo al nivel taxonomico, es este caso será a nivel clase.

```
grep -E "(c__)|(^ID)" merged_abundance_table.txt | grep -v "t__" | sed 's/^.*c__//g' > merged_abundance_table_class.txt				
```

NOTA: Finalmente, la tabla de abundancias "merged_abundance_table.txt" se puede graficar en R o Excel.



```
  if(!requireNamespace("BiocManager")){
    install.packages("BiocManager")
  }
  BiocManager::install("phyloseq")
  
  library (phyloseq)
  
  
  > library (phyloseq)
  > Class_Abundance <- read.csv("AbundanceTable.txt", header=T, row.names=1, sep="\t")
  > OTU = otu_table(Class_Abundance, taxa_are_rows = TRUE)
  > taxmat = matrix(row.names(Class_Abundance))
  > TAX = taxa_names(taxmat)
  > physeq = phyloseq(OTU, TAX)
  > plot_richness(physeq, measures = c("Observed", "Chao1", "Shannon", "simpson"))
```

​			

























