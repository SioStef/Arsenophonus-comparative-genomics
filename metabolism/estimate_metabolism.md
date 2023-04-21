## Workflow for the estimation of metabolic potential across Arsenophonus genomes using [anvio](https://merenlab.org/tutorials/fmt-mag-metabolism/).

1. Starting with a set of genome files run anvi-script-reformat-fasta to prepare the data
>`for i in *.fasta ; do anvi-script-reformat-fasta ${i} -o ${i}-fixed.fasta --seq-type NT --simplify-names; done`


2. Generate a contig database for each genome
>```for i in *-fixed.fasta; do d=`basename ${i} .fasta-fixed.fasta`; anvi-gen-contigs-database -f ${i} -o ${d}.db -n ${d}_database; done```


3. Annotate the genomes with KOfam hits
>```for i in *.db; do anvi-run-kegg-kofams -c ${i} --num-threads 6; done```


4. Prepare an external-genome descriptions file  
>```echo -e "name\tcontigs_db_path" > genomes```
>
>```for i in *.db; do d=`basename ${i} .db`; echo -e "${d}\t${PWD}/${i}" >> genomes; done```


5. Create a new anvi’o genomes storage
>`anvi-gen-genomes-storage -e genomes -o GENOMES.db`


6. Estimate metabolism using the default module completion ratio
>`anvi-estimate-metabolism -e genomes --matrix-format --include-metadata`
