# Python 3.4
# Author:  Robert Newby
# Description: Pt. 1 of 3 script pipeline to process dinoflagellate RNAseq output from Trinity. Takes
#              Trinity output file and extracts ORFs >40aa from all 6 reading frames.
# License:  Author holds no IP on the content within. This software is provided freely for use and/or
#           modification, as long as appropriate citations to the parent article are included in any
#           projects/publications/patents.




from Bio import SeqIO
import os

# List the names of the genes that are to be used as they appear in the Trinity output file.
name_list = ['ATPA', 'ATPB', 'PETB', 'PETD', 'PSAA', 'PSAB', 'PSBA', 'PSBB', 'PSBC', 'PSBD', 'PSBE']

# Set the minimum protein length, and the translation table to be used.
min_pro_len = 40
table = 1


# Takes each record from the queryfile (Trinity output file) and finds all ORFs >= the
# min_pro_len argument frome the call, which it then creates a fasta entry for in
# the out_file.

def ORF_finder(seq, min_pro_len):
    counter = 0
    print(str(protein_name), str(record.id))
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            for pro in nuc[frame:].translate(table).split("*"):
                pro = pro.strip('F')
                pro = pro.strip('K')
                if len(pro) >= min_pro_len:
                    counter += 1
                    out_file.write('>found_in_' + record.id + '_' + str(counter) + '_strand_' + \
                                   str(strand) + '_frame_' + str(frame + 1) + '\n')
                    out_file.write(str(pro) + '\n')
    out_file.close()


# FILE HANDLING AND SENDING SEQS TO ORF_FINDER
# Iterate through the name list, adding the names into the file name template
# and parse the sequences from each file into SeqIO records. Send the mRNA
# sequence from each record to the ORF_finder function.
for protein_name in name_list:
    queryfile = 'pyrocystis_assemblies/%s contig_trinityRNA_chloro_e001-%s.txt' % (protein_name, \
                                                                                   protein_name)
    # open output file in output directory
    file = "/pyrocystis_plastid_ORFs/%s_filtered_ORFs.fasta" % protein_name
    path = os.getcwd() + file
    out_file = open(path, 'w')
    for record in SeqIO.parse(queryfile, 'fasta'):
        orf_prot = ORF_finder(record.seq, min_pro_len)