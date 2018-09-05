
if [ ! -d "cds_splice_work" ]; then
	echo Making new directory
	mkdir cds_splice_work
fi

raxmlHPC-PTHREADS -f a -x 12345 -p 12345 -# 1000 -m GTRGAMMA -s EBOV_CDS_splice_39_seqs.fa -o EBOV_AF272001_Mayinga_Yambuku-DRC_1976,EBOV_KC242801_deRoover_DRC_1976,EBOV_KC242791_Bonduni_Tandala_DRC_1977 -n cds_splice_work/raxml_cds_splice_tree -T 4

