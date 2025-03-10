# de Novo Assembly tutorial

#### Gamaliel López-Leal

###### Laboratorio de Biología Computacional y Virómica Integrativa

###### Centro de Investigación en Dinámica Celular-UAEM

gamaliel.lopez@docentes.uaem.edu.mx	



Crea en tu directorio `home` tres subdirectorios llamados `Reads`, `Genomes` y `Scripts`. 

 `home/Reads`: Contiene las lecturas crudas en formato fastq (archivos .fastq) que se usarán en eltaller para ensamblar genomas virales.

 `home/Genomes`: Contiene algunos genomas de fagos y hospederos en formato fasta (archivos .fna).

 `/home/Scripts`:Contiene algunos programas caseros escritos en perl (archivos .pl) que nosayudará a manejar los archivos durante el taller.



En nuestro directorio de Reads, vamos a encontrar las lecturas crudas producto de la secuenciación. Antes de proceder a realizar el ensamblado de los genomas, es necesario evaluar la calidad de las lecturas y eliminar las lecturas o nucleótidos de baja calidad. Para esto, usaremos las herramientas de FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) y TrimGalore (https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md).

Podemos verificar cuantas lecturas (contar cuantas lecturas hay en cada archivo) tenemos por archivo (R1 y R2) usando el comando `grep` con el flag `-c`. Para esto tenemos que buscar un patro único para cada lectura.

```
grep -c "@SRR" *.fastq
```

### Calidad de la Secuencias

Para revisar la calidad de las secuencias puedes usar la herramienta `fastqc`

```
fastqc Reads-R1.fastq Reads-R2.fastq 
```

### Filtrado de Lecturas por Calidad

```
trim_galore --fastqc -q 30 --paired Reads-R1-v3.fastq Reads-R2-v3.fastq -o Quality_Reads
```

Ahora las lecturas de calidad estará en el directorio `Quality_Reads`. Se recomienda verificar que el proceso se ejecuto bien contando cuantas lecturas de calidad tenemos.

```
grep -c "@SRR" *.fastq
```



### Ensambles *de Novo* de genomas virales

Para este paso usaremos un programa que se llama **metaviralsapdes**. Es de la colección de herramientas
para en ensambles de **spades**. Para mayor información pueden consultar en siguiente enlace: https://github.com/ablab/spades

Antes de ejecutar `metaviralspades.py`  debes de correr el siguiente comando:

```
PATH=/home/App/SPAdes-3.15.5:$PATH
```

Nota: Esta es una particularidad de nuestro servidor y no es necesario correr ese comando en
otros servidores o equipos de computo.

Como cualquier otra herramienta bioinformática podemos consultar el menú con los comandos: `--help`,
`-h`, `man`, por mencionar algunos.

```
metaviralspades.py --only-assembler  -k 21,33,55,77,99,127 -1 Reads-R1-v3_val_1.fq -2 Reads-R2-v3_val_2.fq -o MetaViral_results
```

Ahora en el subdirectorio `MetaViral_results` se encuentran los resultados del ensamble. Revisemos el archivo **contigs.fasta**, que es el contiene los genomas virales.

Es importante saber que existen varias herramientas para revisar la estadística general de los ensambles. Para esto, usaremos una herramienta llamada `quast`. Para mayor información consultar el siguiente enlace:  https://github.com/ablab/quast

```
quast.py contigs.fasta -o quast_results
```



















































































