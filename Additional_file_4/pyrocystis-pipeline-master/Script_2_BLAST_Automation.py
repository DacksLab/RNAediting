# Python 3.4
# Author:  Robert Newby
# Description: Pt. 2 of 3 script pipeline to process dinoflagellate RNAseq output from Trinity. Takes
#              ORFs extracted by the first script and BLASTp searches each against a custom database
#              of dinoflagellate protein sequences obtained from MMETSP project paper.
# License:  Author holds no IP on the content within. This software is provided freely for use and/or
#           modification, as long as appropriate citations to the parent article are included in any
#           projects/publications/patents.




import subprocess
import os

# List the names of the genes that are to be used, as they appear in the input file names.
name_list = ['ATPA', 'ATPB', 'PETB', 'PETD', 'PSAA', 'PSAB', 'PSBA', 'PSBB', 'PSBC', 'PSBD', 'PSBE']

# Use the name to access the filtered_ORFs.fasta file and send it to the BLAST subprocess.
# Capture output to XML file for each gene ran.
for name in name_list:
    print("BLASTing %s..." % name)
    in_file_name = os.getcwd() + '/pyrocystis_plastid_ORFs/%s_filtered_ORFs.fasta' % (name)
    out_file_name = os.getcwd() + '/pyrocystis_xml/%s_trinity_output.xml' % name

    if name == 'PSBE' or 'PETD':
        subprocess.call(['/Users/rob/NCBI/ncbi-blast-2.2.31+/bin/blastp', '-query',
                         in_file_name, '-out', out_file_name, '-outfmt', '5',
                         '-db', 'MMETSP_db.fasta', '-num_threads', '8', '-evalue', '1e-25',
                         '-max_target_seqs', '1'])
    else:
        subprocess.call(['/Users/rob/NCBI/ncbi-blast-2.2.31+/bin/blastp', '-query',
                         in_file_name, '-out', out_file_name, '-outfmt', '5',
                         '-db', 'MMETSP_db.fasta', '-num_threads', '8', '-evalue', '1e-150',
                         '-max_target_seqs', '1'])
    print("Complete")