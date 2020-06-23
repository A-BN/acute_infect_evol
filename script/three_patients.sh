#! /bin/bash

# This script will annotate the reference genomes of the 'three patients' study
# and run breseq with the query reads.



###############################################################################
###############   PATH and other preliminary var settings   ###################
###############################################################################
export PATH=${PATH}:~/Documents/tool/barrnap/bin/:~/Documents/tool/prokka/bin/
export PATH=${PATH}:~/Documents/tool/breseq-0.32.0/bin/
export PATH=${PATH}:~/Documents/tool/bowtie2-2.3.4.1-linux-x86_64/
n_threadus=12

# Message:
echo -e "§§§ Launched:\n\t${0} ${@} on the $(date "+%D at %H:%M:%S")"

if [ "$#" == 0 ]; then
	echo -e "--filter:\tcontig filtering w/ fastaUtils"
	echo -e "--annotate:\tannotation w/ prokka"
	echo -e "--variants:\tvariants calling w/ breseq"
fi

###############################################################################
#####################     Contigs filtering    ################################
###############################################################################
ref_dir="../data/reference"
min_length=100
if [[ $@ =~ "--filter" ]];then
	echo "§§§ Starting Filtering with fastaUtils"
	ref_fasta=$(find ${ref_dir} -type f -name "*.fasta")
	for fastus in ${ref_fasta}; do
		echo "Current fasta = ${fastus}"
		Rscript --slave -e "library(fastaUtils);fastaUtils::fastagrab(fasta = '${fastus}', min_size=${min_length})"
	done
fi

###############################################################################
####################   Annotation with prokka    ##############################
###############################################################################
ref_dir="../data/reference"
if [[ $@ =~ "--annotate" ]];then
	echo "§§§ Starting Annotation with prokka"
	ref_fasta=$(find ${ref_dir} -type f -name "*filtered.fasta")
	for fastus in ${ref_fasta}; do
		echo "Current fasta = ${fastus}"
		curr_patient=`echo ${fastus} | sed -r 's:.*(patient_[0-9]+).*:\1:'`
		echo ${curr_patient}
		annot_dir="$(dirname ${fastus})/annotation"
		if [ ! -d ${annot_dir} ]; then mkdir ${annot_dir};fi
		prokka \
			${fastus} \
			--outdir ${annot_dir} \
			--prefix ${curr_patient}\
			--force \
			--genus Escherichia \
			--species coli \
			--usegenus \
			--mincontiglen 2000 \
			--cpus ${n_threadus}
	done
fi

###############################################################################
####################   Variant calling  with breseq    ########################
###############################################################################
var_dir="../data/variants"
read_dir="../data/reads"
if [[ $@ =~ "--variants" ]]; then
	echo "Starting variant calling with breseq"
	ref_gff=$(find ${ref_dir} -type f -name "patient*.gff")
	echo -e "References found:\n 	${ref_gff}"
	for curr_ref in ${ref_gff}; do
		curr_patient=$(echo ${curr_ref} | sed -r 's:.*(patient_[0-9]+).*:\1:')
		all_reada=$(find "${read_dir}/${curr_patient}" -type f -name "*fastq*" ! -name "*R2*" )
		echo "§§§ Treating ${curr_patient}"
		echo -e "\t reference: ${curr_ref}"
		echo "${all_reada}"
		for curr_readus in ${all_reada}; do
			if [[ "${curr_readus}" =~ .*R1.* ]]; then
				echo "pair-end reads"
				curr_R1=${curr_readus}
				curr_R2=$(echo ${curr_R1} | sed 's/R1/R2/')
				curr_readus="${curr_R1} ${curr_R2}"
				curr_outdir=${var_dir}/${curr_patient}/$(basename ${curr_R1})
			echo "*** Treating ${curr_readus}"
			echo "OUT: ${curr_outdir}"
			elif [[ "$curr_readus" =~ .*R2.* ]]; then
				echo "*** Treating ${curr_readus}"
				echo "This is an R2, skipping"
				curr_outdir="skip"
			else
				curr_outdir=${var_dir}/${curr_patient}/$(basename ${curr_readus})				
				echo "*** Treating ${curr_readus}"
				echo "OUT: ${curr_outdir}"
			fi
			if [ ${curr_outdir} != "skip" ]; then
				breseq \
					--reference ${curr_ref}\
					--name ${curr_patient} \
					--output ${curr_outdir} \
					--num-processors ${n_threadus} \
					${curr_readus}
				mv ${curr_outdir}/output/output.gd ${curr_outdir}/output/$(basename ${curr_readus}).gd
				echo "*** Breseq done for ${curr_readus}"
			fi
		done
	done
fi
