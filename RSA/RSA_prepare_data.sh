#!/bin/bash
# This pipeline prepares the data for relaxation of selection (RELAX) analysis for a set of highly conserved single-copy orthologs between a set of genomes of interest.
# It takes as an input a directory with genome sequences (e.g. .fna) and:
# 1) annotates each genome file using prokka
# 2) identifies Orthologous groups of proteins between the genomes
# 3) extracts single copy core proteins and
# 4) performs a phylogenetic analysis on the concatenated dataset using iqtree
#
# The following packages are needed to be installed and in your PATH.
# prokka (https://github.com/tseemann/prokka)
# orthofinder (https://github.com/davidemms/OrthoFinder)
# seqkit (https://github.com/shenwei356/seqkit)
# busco (http://busco.ezlab.org/)
# clipkit (https://github.com/jlsteenwyk/clipkit)
# iqtree (http://www.iqtree.org/)
# pal2nal (http://www.bork.embl.de/pal2nal/)
# All packages can be installed through conda using the command "conda install -c bioconda prokka orthofinder seqkit busco clipkit iqtree pal2nal"


echo Please provide the directory containing the genomes # note: genome file names should be simplified and without underscores.

read -p 'Genomes folder: ' datag

###########################################################
## Perform genome annotation with prokka for consistency ##
###########################################################
mkdir -p prokka_annotations
for i in ${datag}/*.*
do
NAME=`basename $i .fna`
echo "annotating " ${NAME}
prokka --outdir prokka_annotations/${NAME}.prokka --norrna --notrna --quiet --locustag ${NAME} --prefix ${NAME} ${i}
mkdir -p proteome
ln -s $(pwd)/prokka_annotations/${NAME}.prokka/${NAME}.faa proteome/
done

############################################################################################################################
## Search for orthologous groups of proteins between the genomes using Orthofinder and extract single_copy core orthologs ##
############################################################################################################################
echo ""
echo "Running Orthofider"
orthofinder -M msa -oa -f proteome -S blast -t 16

mkdir core_aln
for i in $(ls proteome/OrthoFinder/Results_*/Single_Copy_Orthologue_Sequences); do cp MultipleSequenceAlignments/${i} core_aln; done

# optional step: identify and extract only highly conserved BUSCO gammaproteobacteria_odb10 markers
mkdir core_aln_busco
busco -c 8 -e 0.00001 -l gammaproteobacteria_odb10 -m prot -i $(find prokka_annotations/*.prokka -type d | shuf -n 1)/*.faa -o busco_tmp #run busco on any genome at protein mode
awk 'NR>3' busco_tmp/run_gammaproteobacteria_odb10/full_table.tsv | cut -f3 > bids
cp $(grep -l -f bids core_aln/*.fa) core_aln_busco
rm bids
rm -r busco_tmp

#######################################################################################################################
## Get cds fasta files and convert protein sequence alignments into the corresponding codon alignments using pal2nal ##
#######################################################################################################################
mkdir core_cds_busco
for i in core_aln_busco/*.fa; do 
OG=`basename $i .fa`;
grep ">" ${i} | awk -F\_ '{print $(NF-1) "_" $NF}' > tmp;
find prokka_annotations/. -name "*.ffn" -exec cat {} \; | seqkit grep -f tmp > core_cds_busco/${OG}.fasta;
rm tmp
done

# rename headers in protein and cds files
for i in core_aln_busco/*.fa
do
sed -i -r "s/_.*$//g" $i
done
for i in core_cds_busco/*.fasta
do
sed -i -r "s/_.*$//g" $i
done

mkdir core_pal2nal_busco
for i in core_cds_busco/*.fasta; do
d=`basename ${i} .fasta`
pal2nal.pl  core_aln_busco/${d}.fa  ${i} -output fasta  > core_pal2nal_busco/${d}.pal2nal
done

########################################################################
## Perform ML core genome phylogeny in iqtree using the busco dataset ##
########################################################################
cd core_aln_busco
#concatenate core protein sequences
seqkit concat *.fa -o core_concatenated.fa --quiet
# trim alignment using clipkit
clipkit core_concatenated.fa -m smart-gap
#phylogenetic analysis on the concatenated dataset using iqtree 
mkdir phylo
iqtree -s core_concatenated.fa.clipkit -bb 1000 -nt AUTO -pre phylo/core_phylo_busco
cd ..
