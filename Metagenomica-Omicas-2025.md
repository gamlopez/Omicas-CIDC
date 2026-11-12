#### Written by Gamaliel López Leal

gamaliel.lopez@docentes.uaem.edu.mx	



### Clasificación Taxónomica con KRAKEN

Ahora usaremos Kraken para clasificar nuestros reads. Cuando usas Karken, debes de formatear la o las bases de datos de interes. Para esta páctica usaremos la  miniKraken de 4GB que se encuentra en: <u>/home/krakenDB/minikraken_20171013_4GB/</u>. Para mayor información consulta: https://ccb.jhu.edu/software/kraken/MANUAL.html

Las lecturas que usaremos estan en: /home/gama/Omicas/Metagenomica

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
merge_metaphlan_tables.py metaphlan2_profile1.txt metaphlan2_profile2.txt metaphlan2_profile1.1.txt > merged_abundance_table.txt

#o de esta manera:

paste metaphlan2_profile1.txt <(awk '{print $2}' metaphlan2_profile2.txt)
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



### Análisis de diversidad con `phyloseq` y `vegan`

phyloseq: https://joey711.github.io/phyloseq/index.html
vegan: https://peat-clark.github.io/BIO381/veganTutorial.html

```
  if(!requireNamespace("BiocManager")){
    install.packages("BiocManager")
  }
  BiocManager::install("phyloseq")
  
library (phyloseq)
library(phyloseq)
library(vegan)
library(ggplot2)
  
#version 1  
Class_Abundance <- read.csv("AbundanceTable.txt", header=T, row.names=1, sep="\t")
OTU = otu_table(Class_Abundance, taxa_are_rows = TRUE)
taxmat = matrix(row.names(Class_Abundance))
TAX = taxa_names(taxmat)
physeq = phyloseq(OTU, TAX)
plot_richness(physeq, measures = c("Observed", "Chao1", "Shannon", "simpson"))

#version 2
Order_Abundance <- read.delim(
  "AbundanceTable.txt",
  header = TRUE,
  row.names = 1
)

counts <- as.matrix(Order_Abundance)
OTU <- otu_table(counts, taxa_are_rows = TRUE)

taxmat <- matrix(
  rownames(counts),
  ncol = 1,
  dimnames = list(rownames(counts), "Order")
)
TAX <- tax_table(taxmat)

physeq <- phyloseq(OTU, TAX)

# Composición por muestra en número de lecturas
plot_bar(physeq, fill = "Order")

#Abundancia relativa
physeq_rel <- transform_sample_counts(
  physeq,
  function(x) 100 * x / sum(x)
)

plot_bar(physeq_rel, x = "Sample", fill = "Order") +
  ggplot2::labs(
    x = "Muestra",
    y = "Abundancia relativa (%)",
    fill = "Orden bacteriano"
  )

# Diversidad alfa a nivel de orden
plot_richness(
  physeq,
  measures = c("Observed", "Shannon", "Simpson")
)

# Diversidad beta entre muestras
conteos <- read.delim("AbundanceTable.txt", row.names = 1)
meta <- read.delim("Metadata.txt", row.names = 1)
comunidad <- t(as.matrix(conteos))
stopifnot(setequal(rownames(comunidad), rownames(meta)))
meta <- meta[rownames(comunidad), , drop = FALSE]
profundidad <- min(rowSums(comunidad))
paleta <- c(Sedimento = "#8C510A", Suelo = "#4D9221", Agua = "#2166AC")

rarecurve(
  comunidad, step = 10, sample = profundidad,
  col = unname(paleta[meta$Ambiente]), label = TRUE,
  xlab = "Número de lecturas",
  ylab = "Número esperado de órdenes"
)
legend("bottomright", legend = names(paleta), col = paleta, lty = 1, bty = "n")

bray <- vegdist(decostand(comunidad, "total"), method = "bray")
jaccard <- vegdist(comunidad, method = "jaccard", binary = TRUE)

round(as.matrix(bray), 3)
round(as.matrix(jaccard), 3)

# Riqueza gamma global y riqueza acumulada por ambiente
sum(colSums(comunidad) > 0)
specnumber(comunidad, groups = meta$Ambiente)

#plot para beta diversidad

beta <- rbind(
  transform(as.data.frame(as.table(as.matrix(bray))),
            Indice = "Bray–Curtis"),
  transform(as.data.frame(as.table(as.matrix(jaccard))),
            Indice = "Jaccard")
)

ggplot(beta, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "white") +
  facet_wrap(~ Indice) +
  scale_fill_gradient(
    low = "#F0F7FC", high = "#084594", limits = c(0, 1)
  ) +
  scale_y_discrete(limits = rev(rownames(comunidad))) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "Disimilitud") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
#plot distancias beta diversidad

distancias <- list(
  "Bray–Curtis" = bray,
  "Jaccard binario" = jaccard
)

coordenadas <- do.call(rbind, lapply(names(distancias), function(indice) {
  ajuste <- cmdscale(distancias[[indice]], k = 2, add = TRUE)
  puntos <- ajuste$points
  ids <- rownames(puntos)

  data.frame(
    Muestra = ids,
    Ambiente = meta[ids, "Ambiente"],
    PCoA1 = puntos[, 1],
    PCoA2 = puntos[, 2],
    Indice = indice,
    row.names = NULL
  )
}))

coordenadas$Indice <- factor(coordenadas$Indice, levels = names(distancias))

ggplot(coordenadas, aes(PCoA1, PCoA2, color = Ambiente)) +
  geom_point(size = 3) +
  geom_text(aes(label = Muestra), vjust = -0.9, size = 3,
            show.legend = FALSE) +
  facet_wrap(~ Indice) +
  coord_equal() +
  scale_color_manual(values = c(
    Sedimento = "#8C510A",
    Suelo = "#4D9221",
    Agua = "#2166AC"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  labs(x = "Eje PCoA 1", y = "Eje PCoA 2", color = "Ambiente") +
  theme_minimal()
  

#plot para gamma diversidad

gamma_ambiente <- specnumber(comunidad, groups = meta$Ambiente)

gamma <- data.frame(
  Grupo = c(names(gamma_ambiente), "Global"),
  Ordenes = c(unname(gamma_ambiente), sum(colSums(comunidad) > 0))
)

gamma$Grupo <- factor(
  gamma$Grupo, levels = c("Sedimento", "Suelo", "Agua", "Global")
)

ggplot(gamma, aes(x = Grupo, y = Ordenes, fill = Grupo)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = Ordenes), vjust = -0.4) +
  scale_fill_manual(values = c(
    Sedimento = "#8C510A", Suelo = "#4D9221",
    Agua = "#2166AC", Global = "#636363"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Número de órdenes observados") +
  theme_minimal()

#¿cuántos órdenes distintos se observaron al reunir las nueve muestras?
sum(colSums(comunidad) > 0)




```

​ 


		

























