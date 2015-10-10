This program is searching for tRNA genes inside genome. The program is divided into two scripts. The script get_tRNA.py is searching for a tRNA genes and compare_tRNA is comparing founded genes with records from Genomic tRNA Database.

Run:
	python get_tRNA.py Seq/sequence.fasta
	python compare_tRNA.py Seq/found_tRNA.multifasta Seq/known_tRNA.multifasta