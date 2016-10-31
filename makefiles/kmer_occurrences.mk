################################################################
## Prepare the data for the ASG course about motif discovery in
## ChIP-seq peaks.

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/kmer_occurrences.mk

QUICK_OPT=
#QUICK_OPT=-quick
OL=6
NOOV=-noov
STRANDS=-2str
OLIGOS=mm10_genome_${OL}nt${NOOV}${STRANDS}
OLIGO_CMD=oligo-analysis -v ${V} \
		-i mm10_genome.fasta \
		-return occ,freq \
		${QUICK_OPT} ${NOOV} ${STRANDS} -l ${OL} -grouprc \
		-o ${OLIGOS}.tab 2> ${OLIGOS}_log.txt
oligos_one_len:
	@echo 
	@echo "Computing genomic frequencies of ${OL}-mers"
	@echo "Log file"
	@echo "	${OLIGOS}_log.txt"
	${MAKE} my_command MY_COMMAND="${OLIGO_CMD}"
	@echo "Result file"
	@echo "	${OLIGOS}.tab"



################################################################
## Count k-mer occurrences in peak sequences
PEAK_DIR=data/kmer_occurrences/CEBPA_mm9_peaks_Ballester_2010
PEAKS=${PEAK_DIR}/CEBPA_mm9_SWEMBL_R0.12
PEAK_OLIGOS=${PEAKS}_${OL}nt${NOOV}${STRANDS}
peak_oligos_one_len:
	@echo
	@echo "Computing ${OL}-mer occurrences in peaks ${PEAKS}.fasta"
	oligo-analysis -v ${V} -l ${OL} -return occ,freq ${NOOV} ${STRANDS}  -i ${PEAKS}.fasta -o ${PEAK_OLIGOS}.tab
	@echo "	${PEAK_OLIGOS}.tab"

peak_oligos:
	${MAKE} oligos OLIGO_TASK=peak_oligos_one_len

## Compute k-mer frequencies for various lengths
OLIGO_LENGTHS=6 7 5 4 3 2 1 8
OLIGO_TASK=peak_oligos_one_len
oligos:
	@for ol in ${OLIGO_LENGTHS}; do \
		${MAKE} ${OLIGO_TASK} OL=$${ol} ; \
	done


################################################################
## Count k-mer  occurrences in random genomic regions
RAND_DIR=data/kmer_occurrences/random_fragments_mm10/
REPEAT=01
RAND_SEQ=${RAND_DIR}/random-genome-fragments_mm10_repeat${REPEAT}
rand_oligos_one_rep:
	${MAKE} oligos PEAKS=${RAND_SEQ}

REPETITIONS=01 02 03 04 05 06 07 08
rand_oligos:
	for rep in ${REPETITIONS}; do \
		${MAKE} rand_oligos_one_rep REPEAT=$${rep}; \
	done
