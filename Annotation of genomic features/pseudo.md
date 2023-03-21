## How to run ideel to calculate the fraction of pseudogenised genes in the Arsenophonus genomes.

1. Clone ideel repo from [https://github.com/mw55309/ideel]

2. Make a directory called "proteins" and put all protein fasta files (.faa)

`ln -s $(pwd)/prokka_annotations/${NAME}.prokka/${NAME}.faa proteins/`

3. Download and format for [Diamond](https://github.com/bbuchfink/diamond/wiki) the uniprot Swiss-Prot database from [https://www.uniprot.org/help/downloads]

4. Modify the Snakefile as shown [here](./Snakefile_mod)

5. run snakemake within the ideel directory

`snakemake` 

6. Get the proportion of interupted genes from ideel .data output files

```
cd lengths
for i in *.data; do cat $i | awk -v VAR=$i  '{ c+=($1/$2 <0.8) } END { print VAR, c, NR, c/NR }' >> summary.txt ; done
```
