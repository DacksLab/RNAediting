# Python 3.4
# Author:  Robert Newby
# Description: Pt. 3 of 3 script pipeline to process dinoflagellate RNAseq output from Trinity. Parses
#              BLASTp results from Step 2, and selects the most likely protein coding contig based on
#              BLASTp evalue and contig length (longer is ranked higher to get complete mRNA)
# License:  Author holds no IP on the content within. This software is provided freely for use and/or
#           modification, as long as appropriate citations to the parent article are included in any
#           projects/publications/patents.




import os
from Bio import SeqIO
from Bio.Seq import Seq

# List the names of the genes to be looked at. These names need to appear as they are found in the
# file names of the BLAST output XML files.
name_list = ['ATPA', 'ATPB', 'PETB', 'PETD', 'PSAA', 'PSAB', 'PSBA', 'PSBB', 'PSBC', 'PSBD', 'PSBE']

# Initialize variables.
match_length = {}
full_match_id = {}
output_file = open(os.getcwd() + '/pyrocystis_coding_sequences/pyrocystis_coding_seqs.fasta', 'w')
prot_out_file = open(os.getcwd() + '/pyrocystis_coding_sequences/pyrocystis_protein_seqs.fasta', 'w')


# Uses name to build a path to the ORF file for the protein, opens and checks line-byline for def-id
# if line.startswith('>'). When match is found, the header and following sequence are written to
# output file in fasta format. Parse the ORF file into SeqIO, and retrieve the record.seq with def_id
def get_ORF_entry(name, def_id):
    ORF_path = os.getcwd() + '/pyrocystis_plastid_ORFs/%s_filtered_ORFs.fasta' % (name)
    list_of_hits = []

    for record in SeqIO.parse(ORF_path, 'fasta'):
        if def_id in record.id:
            ORF = str(record.seq)
            ORF_length = len(ORF)
            match_length[name].append([def_id, ORF_length, ORF])


# Takes the current working gene name and the orf_id, and uses them to access the filtered_ORFs.fasta
# file, get the matched ORF sequence, and write it to a fasta entry in the protein output file.
def get_ORF(name, orf_id):
    ORF_path = os.getcwd() + '/pyrocystis_plastid_ORFs/%s_filtered_ORFs.fasta' % name
    for record in SeqIO.parse(ORF_path, 'fasta'):
        if orf_id in record.id:
            prot_out_file.write('>%s %s\n %s \n' % (name, record.id, record.seq))


# Takes the current working gene name, the contig_id, and coding strand, and uses them to access the
# Trinity output file, get the contig sequence, and write it to a fasta entry in the transcript output
# file. If the strand attribute is -1, then the reverse complement of the sequence is written instead.
def get_contig(name, contig_id, strand):
    print("Checking %s" % name)
    contig_path = os.getcwd() + '/pyrocystis_assemblies/%s contig_trinityRNA_chloro_e001-%s.txt' % (name,
                                                                                                    name)
    # print(contig_path)
    for record in SeqIO.parse(contig_path, 'fasta'):
        if contig_id in record.id:
            contig_seq = record.seq
            if strand == '1':
                print("found in forward")
                output_file.write('>%s %s\n %s \n' % (name, record.id, contig_seq))
            elif strand == '-1':
                print('found in reverse')
                output_file.write('>%s %s_RC\n %s \n' % (name, record.id, contig_seq.reverse_complement()))


# STAGE I
# Takes the gene names from the name_list, adds them to the match_length dict for later, and adds them into
# the file path string to open their BLAST XML output file. Determines matches by the presence of <Hsp_qseq>
# tag, and sends current gene name and the last recorded Iteration_query-def to the get_ORF_entry() function.
for name in name_list:
    print(name)
    match_length[name] = []
    query_list = []
    in_file_name = os.getcwd() + '/pyrocystis_xml/%s_trinity_output.xml' % name
    xml_out = open(in_file_name, 'r')

    for line in xml_out:
        if '<Iteration_query-def>' in line:
            query_list.append(line.strip().lstrip('<Iteration_query-def>').rstrip('</Iteration_query-def>'))
        elif '<Hsp_qseq>' in line:
            get_ORF_entry(name, query_list[-1])  # function here to get the full translation from the ORF file

# STAGE II
# Determines the length of the longest match(es) and adds their unique id string to the full_match_id
# dictionary, so they may be used to locate the fasta entries for the matched contig and ORF in their
# respective files.
for k, v in match_length.items():
    full_match_id[k] = []
    len_matches = []

    # Make a list of the lengths of the matches for each key
    for i in range(len(v)):
        len_matches.append(v[i][1])
    # print(k, max(len_matches))

    # Finds the longest matched seqs and appends the id (I think) to the full_match_id dictionary.
    for i in range(len(v)):
        if v[i][1] == max(len_matches):
            full_match_id[k].append(v[i][0])
            # print(k, len(v))

# STAGE III
# Gets the contig and protein sequence of the matches stored in full_match_id using functions for
# each. The functions write the data as fasta entries in transcript and protein sequence files.
for k, v in full_match_id.items():
    print(k, v)
    get_ORF(k, v[0])
    contig_id = v[0].split('_')
    record_id = "_".join(contig_id[:3])
    # print(contig_id)
    print(contig_id[5])
    get_contig(k, record_id, contig_id[5])

output_file.close()
prot_out_file.close()