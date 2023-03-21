## How to run ideel to calculate the fraction of pseudogenised genes in the Arsenophonus genomes.

1. Clone ideel repo from [https://github.com/mw55309/ideel]
<p>2. Make a directory called "proteins" and put all protein fasta files (.faa)<br>
`ln -s $(pwd)/prokka_annotations/${NAME}.prokka/${NAME}.faa proteins/`
3. Download and format for [Diamond](https://github.com/bbuchfink/diamond/wiki) the uniprot Swiss-Prot database from [https://www.uniprot.org/help/downloads]
4. Modify the Snakefile as shown [here](./snakefile_mod)
5. run snakemake within the ideel directory

`snakemake` 
