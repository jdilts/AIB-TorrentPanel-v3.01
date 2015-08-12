echo ////////////////////////////////////////////////////////////
echo Company:	 American International Biotechnology, LLC
echo Department: Bioinformatics
echo Author:	 Jeri Dilts and William Budd
echo Date:	 April 21st, 2014
echo ////////////////////////////////////////////////////////////

echo ---------------------------------------------
echo DESCRIPTION
echo ---------------------------------------------
echo This plugin filters mapped reads, executes Life Techs command-line variantCaller, and executes AIBiotechs genotyper.
echo Life Techs Torrent Suite software only allows for the automatic execution of a single plugin after a run analyses.
echo Essentially, this plugin is a combination of 3.
echo Clinical no longer has to manually click/execute any plugins whatsoever...
echo As long as this plugin is selected in their run template.
echo ---------------------------------------------

#!/bin/bash

VERSION="3.00"
echo version: ${VERSION}
#AUTORUNDISABLE

run ()
{
        echo "running: $*";
        eval $*;
        EXIT_CODE="$?";
}

#-----TRANSFER?no/yes-
TRANSFER=yes
#-----METRICS---------
SAMPLE_READ_FILTER=5500
SNP_COVERAGE=20
#-----PLUGIN NAME-----
PLUGIN_NAME=AIB-TorrentPanel-v${VERSION}
#-----REFERENCE-------
NREF=AIBTorrentPanelv300REFERENCE
#-----BED FILES-------
NTARGET=${PLUGIN_NAME}_TARGET
NHOT=${PLUGIN_NAME}_HOT
NNUMTARGET=51
NNUMHOT=52

#-----VARIABLES-------
BGZIP_VCF=/results/uploads/BED/${NNUMHOT}/${NREF}/unmerged/detail/${NHOT}.vcf.gz
HOTSPOT=/results/uploads/BED/${NNUMHOT}/${NREF}/unmerged/detail/${NHOT}.bed
JSON=/results/plugins/variantCaller/pluginMedia/configs/germline_high_stringency.json
LEFT_HOTSPOT=/results/uploads/BED/${NNUMHOT}/${NREF}/unmerged/detail/${NHOT}.left.bed
MP_PYTHON=/results/plugins/variantCaller/scripts/allele_count_mpileup_stdin.py
PA_PYTHON=/results/plugins/variantCaller/scripts/print_allele_counts.py
VCF=/results/uploads/BED/${NNUMHOT}/${NREF}/unmerged/detail/${NHOT}.vcf
REFERENCE=/results/referenceLibrary/tmap-f3/${NREF}/${NREF}.fasta
TARGET=/results/uploads/BED/${NNUMTARGET}/${NREF}/merged/plain/${NTARGET}.bed
VC_PATH=/results/plugins/variantCaller/
VC_PYTHON=/results/plugins/variantCaller/variant_caller_pipeline.py
MODULE_LOCATION=/results/plugins/${PLUGIN_NAME}/analyze/modules
IP_ADDR=`/sbin/ifconfig eth0 | grep "inet addr:*" | awk -F ':' '{ print $2 }' | awk -F ' ' '{ print $1 }'`
LOG=${TSP_FILEPATH_PLUGIN_DIR}/logs
INDICATOR=${TSP_FILEPATH_PLUGIN_DIR}/indicator.txt #don't re-run variantCaller if already ran
RUN=`echo ${TSP_FILEPATH_PLUGIN_DIR} | sed -e 's#.*/Home/\(.*\)/plugin_out/.*#\1#g'`;
AIB_RESULTS=${TSP_FILEPATH_PLUGIN_DIR}/AIB_RESULTS

#-----FLAT_FILES--------
DROPOUT_BARCODES=${AIB_RESULTS}/dropout_barcodes.txt
NON_DROPOUT_BARCODES=${AIB_RESULTS}/non-dropout_barcodes.txt
NT_FREQ_TABLE=${AIB_RESULTS}/nt_freq.xls
GENOTYPE=${AIB_RESULTS}/genotype.txt
SAMPLE_COVERAGE=${AIB_RESULTS}/sample_coverage.xls
SAMPLE_ID=${AIB_RESULTS}/sampleID.txt
POSITIVE_CONTROL_CHECK=${AIB_RESULTS}/positive_control_check.txt

#-----EXPORTS-----------
export IP_ADDR
export VERSION
export MODULE_LOCATION
export DROPOUT_BARCODES
export NON_DROPOUT_BARCODES
export GENOTYPE
export NT_FREQ_TABLE
export SAMPLE_COVERAGE
export SAMPLE_ID
export POSITIVE_CONTROL_CHECK

echo
echo ////////////////////////////////////////////////////////////
echo ONLY NEEDS TO BE EXECUTED WHEN NEW BED FILES ARE UPLOADED
echo ////////////////////////////////////////////////////////////
#-----Convert BED or VCF file into a valid hotspot file-----$tvcutils prepare_hotspots... to see options
run "tvcutils prepare_hotspots -b ${HOTSPOT} -d ${LEFT_HOTSPOT} -o ${VCF} -r ${REFERENCE} -a";
#-----Block compression/decompression utility
run "bgzip -c ${VCF} > ${BGZIP_VCF}";
#-----Generic indexer for TAB-delimited genome position files
run "tabix -p vcf ${BGZIP_VCF}";
echo ////////////////////////////////////////////////////////////


#-----CLEAR OLD ANALYSIS-----
if [ ! -d ${LOG} ]; then
	mkdir ${LOG}
	else
		rm ${LOG}/*
fi

if [ ! -d ${AIB_RESULTS} ]; then
	mkdir ${AIB_RESULTS}
	else
		rm ${AIB_RESULTS}/*
fi

echo
echo ////////////////////////////////////////////////////////////
echo BUILDING HEADERS
echo ////////////////////////////////////////////////////////////
perl -e 'use lib ($ENV{MODULE_LOCATION}); use printTorrent qw(printHeader); printHeader($ENV{IP_ADDR},$ENV{VERSION},$ENV{DROPOUT_BARCODES});'
perl -e 'use lib ($ENV{MODULE_LOCATION}); use printTorrent qw(printHeader); printHeader($ENV{IP_ADDR},$ENV{VERSION},$ENV{NON_DROPOUT_BARCODES});'
perl -e 'use lib ($ENV{MODULE_LOCATION}); use printTorrent qw(printHeader); printHeader($ENV{IP_ADDR},$ENV{VERSION},$ENV{GENOTYPE});'
perl -e 'use lib ($ENV{MODULE_LOCATION}); use printTorrent qw(printHeader); printHeader($ENV{IP_ADDR},$ENV{VERSION},$ENV{NT_FREQ_TABLE});'
perl -e 'use lib ($ENV{MODULE_LOCATION}); use printTorrent qw(printHeader); printHeader($ENV{IP_ADDR},$ENV{VERSION},$ENV{SAMPLE_COVERAGE});'
perl -e 'use lib ($ENV{MODULE_LOCATION}); use printTorrent qw(printHeader); printHeader($ENV{IP_ADDR},$ENV{VERSION},$ENV{SAMPLE_ID});'
perl -e 'use lib ($ENV{MODULE_LOCATION}); use printTorrent qw(printHeader); printHeader($ENV{IP_ADDR},$ENV{VERSION},$ENV{POSITIVE_CONTROL_CHECK});'


echo
echo ////////////////////////////////////////////////////////////
echo FIND DROPOUTS
echo ////////////////////////////////////////////////////////////
run "perl ${DIRNAME}/analyze/dropouts.pl ${RUN} ${DROPOUT_BARCODES} ${NON_DROPOUT_BARCODES} &> ${LOG}/dropouts_error.log";
echo ////////////////////////////////////////////////////////////

if [ ! -f ${INDICATOR} ]; then

	#----RUN CLVC PER SAMPLE-----
	if [ -f ${TSP_FILEPATH_BARCODE_TXT} ]; then

		grep -vE '^#|^$' < ${NON_DROPOUT_BARCODES} | while read SAMPLE

			do
				echo
				echo
				echo
				echo ${SAMPLE}

				#-----DYNAMIC PATHS-----
				ALLELE_COUNTS_TXT=${TSP_FILEPATH_PLUGIN_DIR}/${SAMPLE}/allele_counts.txt
				ALLELE_COUNTS_XLS=${TSP_FILEPATH_PLUGIN_DIR}/${SAMPLE}/allele_counts.xls
				SAMPLE_OUT_PATH=${TSP_FILEPATH_PLUGIN_DIR}/${SAMPLE}
				BAM=${SAMPLE_OUT_PATH}/${SAMPLE}_rawlib.bam

				if [ -d ${SAMPLE_OUT_PATH} ]; then

					rm -r ${SAMPLE_OUT_PATH}
					mkdir ${SAMPLE_OUT_PATH}

				else
					mkdir ${SAMPLE_OUT_PATH}
				fi

				echo
				echo ////////////////////////////////////////////////////////////
				echo EXECUTE MAPPING FILTER
				echo ////////////////////////////////////////////////////////////
				run "perl ${DIRNAME}/analyze/filter.pl ${TSP_FILEPATH_PLUGIN_DIR} ${SAMPLE} &>> ${LOG}/filter_error.log";

				echo ////////////////////////////////////////////////////////////
				echo VARIANT CALLER
				echo ////////////////////////////////////////////////////////////
				DOMAIN=`echo ${SAMPLE} | awk -F- '{print $2}'`;

				if [ "${DOMAIN}" != "blank" ]; then

					#-----Executes variant_caller_pipeline.py
					run "${VC_PYTHON} -b ${TARGET} -s ${VCF} -i ${BAM} -r ${REFERENCE} -o ${SAMPLE_OUT_PATH} -p ${JSON} -B ${VC_PATH} &>> ${LOG}/clvc_error.log";
					#-----Generate BCF or pileup for one or multiple BAM files-----$samtools mpileup... to see options
					run "samtools mpileup -BQ0 -d1000000 -f ${REFERENCE} -l ${VCF} ${BAM} | ${MP_PYTHON} > ${ALLELE_COUNTS_TXT}";
					#-----Executes print_allele_counts.py
					run "${PA_PYTHON} ${ALLELE_COUNTS_TXT} ${ALLELE_COUNTS_XLS} ${LEFT_HOTSPOT} ${HOTSPOT}";
				fi
			done

	fi

	touch "${INDICATOR}";
fi


echo
echo ////////////////////////////////////////////////////////////
echo CALCULATING SAMPLE COVERAGE
echo ////////////////////////////////////////////////////////////
run "perl ${DIRNAME}/analyze/coverage.pl ${NON_DROPOUT_BARCODES} ${SAMPLE_COVERAGE} ${TSP_FILEPATH_PLUGIN_DIR}&> ${LOG}/coverage_error.log";

echo
echo ////////////////////////////////////////////////////////////
echo CREATING NT FREQ FILE
echo ////////////////////////////////////////////////////////////
run "perl ${DIRNAME}/analyze/nt_freq.pl ${NON_DROPOUT_BARCODES} ${NT_FREQ_TABLE} ${TSP_FILEPATH_PLUGIN_DIR} &> ${LOG}/nt_freq.log";

echo
echo ////////////////////////////////////////////////////////////
echo GENOTYPER
echo ////////////////////////////////////////////////////////////
run "perl ${DIRNAME}/analyze/genotyper.pl ${PLUGIN_NAME} ${TSP_FILEPATH_PLUGIN_DIR} ${SAMPLE_READ_FILTER} ${SNP_COVERAGE} ${NT_FREQ_TABLE} ${GENOTYPE} ${SAMPLE_COVERAGE} &> ${LOG}/genotyper_error.log";

echo
echo ////////////////////////////////////////////////////////////
echo POSITIVE CONTROL CHECK
echo ////////////////////////////////////////////////////////////
run "perl ${DIRNAME}/analyze/positive_control_check.pl ${GENOTYPE} ${POSITIVE_CONTROL_CHECK} &> ${LOG}/positive_control_check.log";

echo
echo ////////////////////////////////////////////////////////////
echo INTERFACE REPORT
echo ////////////////////////////////////////////////////////////
run "perl ${DIRNAME}/report/report.pl ${TSP_FILEPATH_PLUGIN_DIR} ${IP_ADDR} ${RUN} ${PLUGIN_NAME} &> ${LOG}/report_error.log";

echo
echo ////////////////////////////////////////////////////////////
echo CONVERT BACK TO OLD/NEW FORMAT
echo ////////////////////////////////////////////////////////////
run "perl ${DIRNAME}/analyze/make_genotype_perfect.pl ${GENOTYPE} ${TSP_FILEPATH_PLUGIN_DIR} &> ${LOG}/make_genotype_perfect.log";

echo
echo ////////////////////////////////////////////////////////////
echo FILE TRANSFERS
echo ////////////////////////////////////////////////////////////

if [ "${TRANSFER}" == "yes" ]; then

	## actual names
	GENOTYPE_PERFECT="${AIB_RESULTS}/genotype.txt"
	GENOTYPE_ERRORS="${AIB_RESULTS}/genotype_errors.txt"

	## transfer names
	NEW_FORMAT="${TSP_RUN_NAME}_GENOTYPE"
	ERROR_FORMAT="${TSP_RUN_NAME}_ERROR"

	## transferring if not empty
	if [ -s ${GENOTYPE_PERFECT} ]; then
		echo cp ${GENOTYPE_PERFECT} /mnt/WindowsShared/Groups/DNA/AGI_pipeline/Torrent/PGX/${NEW_FORMAT}.txt
		cp ${GENOTYPE_PERFECT} /mnt/WindowsShared/Groups/DNA/AGI_pipeline/Torrent/PGX/${NEW_FORMAT}.txt
#		cp ${GENOTYPE_PERFECT} /mnt/WindowsShared/IONTorrent/Results/Processed/${NEW_FORMAT}.txt

	fi

	if [ -s ${GENOTYPE_ERRORS} ]; then
		echo cp ${GENOTYPE_ERRORS} /mnt/WindowsShared/IONTorrent/Results/Errored/${ERROR_FORMAT}.txt
		cp ${GENOTYPE_ERRORS} /mnt/WindowsShared/IONTorrent/Results/Errored/${ERROR_FORMAT}.txt
	fi

	## store data to aibitest
	run "expect ${DIRNAME}/transfer.exp ${RUN} ${TSP_FILEPATH_PLUGIN_DIR}";

	else
	echo no transfer
fi









