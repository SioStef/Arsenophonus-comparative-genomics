### This workflow describes the steps used to screen Arthropode SRA datasets for Arsenophonus reads and the subsequent assembly of draft Arsenophonus genomes.

The following tools are used:

* [mash](https://mash.readthedocs.io/en/latest/index.html)
* [sra_download.pl](https://github.com/gerthmicha/perlscripts/blob/master/sra_download.pl)
* [MEGAHIT](https://github.com/voutcn/megahit)
* [minimap2](https://github.com/lh3/minimap2)
* [MetaBat](https://bitbucket.org/berkeleylab/metabat/src/master/)
* [CheckM](https://github.com/Ecogenomics/CheckM/wiki)
* [anvio](https://anvio.org/install/)


---

1. We need to build a reference sketch from the available Arsenophonus genomes using [mash](https://mash.readthedocs.io/en/latest/index.html). We will use this reference sketch to screen the SRA datasets for containment of Arsenophonus reads.
```
mash sketch -o ArsRef Arsenophonus_reference_genomes/*.fna
```

2. Download SRA datasets from the European Nuleaotide archive (ENA) database using the sra_download.pl script written by [Michael Gerth](https://github.com/gerthmicha/perlscripts/blob/master/sra_download.pl)
```
for i in $(cat sra_list.txt)
do
mkdir ${i}
echo ${i} > tmp
perl sra_download.pl tmp
mv *.gz ${i}
rm tmp
done
```

3. Screen the read set for containment of Arsenophonus reference genomes using mash
```
for i in $(cat sra_list.txt)
do
cat ${i}/*.fastq.gz | mash screen ArsReference.msh - > ${i}.tab
done
```

4. Retreive SRA accessions where Arsenophonus containment is at least 80%
```
for i in *.tab; do s=`basename ${i} .tab`; sort -gr ${i} | head -n1 | awk -v SRA=$s -F'[\t/]' '$2/$3 > 0.8  {print SRA}' >> SRA_with_Ars.txt; done
```

5. The following commands will a) performe metagenomic assembly using the [MEGAHIT](https://github.com/voutcn/megahit) assembler, b) cluster contigs into bins using [MetaBat2](https://bitbucket.org/berkeleylab/metabat/src/master/) and c) Inspect the identified bins for Completeness and Contamination using [CkeckM](https://github.com/Ecogenomics/CheckM/wiki)
```
for i in $(cat SRA_with_Ars.txt)
do
cd ${i}

# Run the assembler and re-map raw reads back to the assembly 
megahit -1 ${i}_1.fastq.gz -2 ${i}_2.fastq.gz -t 32 -o megahitOut_${i} --out-prefix ${i}
minimap2 -t 48 -ax sr megahitOut_${i}/*.contigs.fa ./${i}_1.fastq.gz ./${i}_2.fastq.gz | samtools view -Sbh | samtools sort -o ${i}.sorted.bam

# cluster contigs into bins
runMetaBat.sh -m 1500 megahitOut_${i}/${i}.contigs.fa ${i}.sorted.bam 

# Inspect bins with CheckM
cd ${i}.contigs.fa.metabat-bins1500
checkm lineage_wf -f ${i}.checkm.txt -t 32 -x fa ../${i}.contigs.fa.metabat-bins1500 ./bin.checkm
cd ../..
done
```

6. Further manual refinment of The Arsenophonus bins can be done in [anvio v7](https://anvio.org/install/) by identifying and removing potential contaminant contigs based on atypical coverage and gene-level taxonomic classification following the [Metagenomic Workflow](https://merenlab.org/2016/06/22/anvio-tutorial-v2/) and the tutorial for [Importing taxonomy into contigs database](https://merenlab.org/2016/06/18/importing-taxonomy/).
