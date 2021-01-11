from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

#Fazer os BLASTs!

orfs = ['TRAF3', 'MAVS', 'ORF9b']

for orf in orfs:
    blast = SeqIO.read(open(orf + '.fa'), format='fasta')
    result = NCBIWWW.qblast('blastp', 'nr', blast.format('fasta'), hitlist_size=9, expect=1)
    save_blast = open(orf + '_blast.xml', 'w')
    save_blast.writelines(result.read())
    save_blast.close()
    result.close()


# Ler os resultados do BLAST!

for orf in orfs:
    count = 0
    result_handle = open(orf + '_blast.xml')
    blast_records = NCBIXML.read(result_handle)

    blast_file = open(orf + '_blast_results.txt', 'w')

    for alignment in blast_records.alignments:
        count += 1
        blast_file.writelines('\n')
        blast_file.writelines('***Alignment***\n')
        blast_file.writelines(f'sequence : {alignment.title}\n')
        blast_file.writelines(f'length: {alignment.length}\n')
        blast_file.writelines(f'id: {alignment.hit_id}\n')

        for hsp in alignment.hsps:
            blast_file.writelines(f'e-value: {hsp.expect}\n')
            blast_file.writelines(f'Query: {hsp.query[0:30]}' + '...' + f'{hsp.query[-31:-1]}\n')

    blast_file.writelines(f'Number of alignments: {count}')
    blast_file.close()
    result_handle.close()

# Homologia ;

for orf in orfs:
    orf_file = open(orf + '.fa', 'r')
    a = orf_file.readline().split()
    seq = ''
    for f in orf_file.readlines()[0::]:
        r = f.replace('\n', '')
        seq += r

    result_handle = open(orf + '_blast.xml')
    blast_records = NCBIXML.read(result_handle)
    file = open(orf + '_homology.fa', 'w')
    file.writelines(a[0] + '\n')
    file.writelines(seq + '\n')

    for alignment in blast_records.alignments:
        for hsp in range(len(alignment.hsps)):
            if hsp != 0:
                file.writelines('>' + alignment.hit_id + '\n')
                file.writelines(alignment.hsps[hsp].sbjct + '\n')
            else:
                file.writelines('>' + alignment.hit_id + '\n')
                file.writelines(alignment.hsps[hsp].sbjct + '\n')
    file.close()
    result_handle.close()
    orf_file.close()
