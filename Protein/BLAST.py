from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

#Fazer os BLASTs!
#   TRAF3
TRAF3_blast = SeqIO.read(open('TRAF3.fa'), format ='fasta')
TRAF3_result = NCBIWWW.qblast('blastp', 'nr', TRAF3_blast.format('fasta'), hitlist_size=20, expect=1)
save_TRAF3_blast = open('TRAF3_BLAST.xml', 'w')
save_TRAF3_blast.writelines(TRAF3_result.read())
save_TRAF3_blast.close()
TRAF3_result.close()

#   MAVS
MAVS_blast = SeqIO.read(open('MAVS.fa'), format ='fasta')
MAVS_result = NCBIWWW.qblast('blastp', 'nr', MAVS_blast.format('fasta'), hitlist_size=20, expect=1)
save_MAVS_blast = open('MAVS_BLAST.xml', 'w')
save_MAVS_blast.writelines(MAVS_result.read())
save_MAVS_blast.close()
MAVS_result.close()

#   ORF9b
ORF9b_blast = SeqIO.read(open('ORF9b_COVID19.fa'), format ='fasta')
ORF9b_result = NCBIWWW.qblast('blastp', 'nr', ORF9b_blast.format('fasta'), hitlist_size=20, expect=1)
save_ORF9b_blast = open('ORF9b_BLAST.xml', 'w')
save_ORF9b_blast.writelines(ORF9b_result.read())
save_ORF9b_blast.close()
ORF9b_result.close()


# Ler os resultados do BLAST!
#   TRAF3
count=0
result_handle = open('TRAF3_BLAST.xml')
blast_records = NCBIXML.read(result_handle)

TRAF3_BLAST_file = open('TRAF3_BLAST_results.txt', 'w')

for alignment in blast_records.alignments:
    count += 1
    TRAF3_BLAST_file.writelines('\n')
    TRAF3_BLAST_file.writelines('***Alignment***\n')
    TRAF3_BLAST_file.writelines(f'sequence : {alignment.title}\n')
    TRAF3_BLAST_file.writelines(f'length: {alignment.length}\n')
    TRAF3_BLAST_file.writelines(f'id: {alignment.hit_id}\n')
    
    for hsp in alignment.hsps:
        TRAF3_BLAST_file.writelines(f'e-value: {hsp.expect}\n')
        TRAF3_BLAST_file.writelines(f'Query: {hsp.query[0:30]}' + '...' + f'{hsp.query[-31:-1]}\n')

TRAF3_BLAST_file.writelines(f'Number of alignments: {count}')
TRAF3_BLAST_file.close()
result_handle.close()

#   MAVS
count1=0
result_handle1 = open('MAVS_BLAST.xml')
blast_records1 = NCBIXML.read(result_handle1)

MAVS_BLAST_file = open('MAVS_BLAST_results.txt', 'w')

for alignment in blast_records1.alignments:
    count1 += 1
    MAVS_BLAST_file.writelines('\n')
    MAVS_BLAST_file.writelines('***Alignment***\n')
    MAVS_BLAST_file.writelines(f'sequence : {alignment.title}\n')
    MAVS_BLAST_file.writelines(f'length: {alignment.length}\n')
    MAVS_BLAST_file.writelines(f'id: {alignment.hit_id}\n')
    
    for hsp in alignment.hsps:
        MAVS_BLAST_file.writelines(f'e-value: {hsp.expect}\n')
        MAVS_BLAST_file.writelines(f'Query: {hsp.query[0:30]}' + '...' + f'{hsp.query[-31:-1]}\n')

MAVS_BLAST_file.writelines(f'Number of alignments: {count1}')

MAVS_BLAST_file.close()
result_handle1.close()


#   ORF9b
count2=0
result_handle2 = open('ORF9b_BLAST.xml')
blast_records2 = NCBIXML.read(result_handle2)

ORF9b_BLAST_file = open('ORF9b_BLAST_results.txt', 'w')

for alignment in blast_records2.alignments:
    count2 += 1
    ORF9b_BLAST_file.writelines('\n')
    ORF9b_BLAST_file.writelines('***Alignment***\n')
    ORF9b_BLAST_file.writelines(f'sequence : {alignment.title}\n')
    ORF9b_BLAST_file.writelines(f'length: {alignment.length}\n')
    ORF9b_BLAST_file.writelines(f'id: {alignment.hit_id}\n')
    
    for hsp in alignment.hsps:
        ORF9b_BLAST_file.writelines(f'e-value: {hsp.expect}\n')
        ORF9b_BLAST_file.writelines(f'Query: {hsp.query[0:30]}' + '...' + f'{hsp.query[-31:-1]}\n')

ORF9b_BLAST_file.writelines(f'Number of alignments: {count2}')
ORF9b_BLAST_file.close()
result_handle2.close()

