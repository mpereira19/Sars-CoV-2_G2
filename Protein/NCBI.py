from Bio import Entrez

Entrez.email = 'pg42875@alunos.uminho.pt'

# TRAF3 - aminoacid sequence       Uniprot id = Q13114
handle = Entrez.efetch(db='Protein', id='Q13114', rettype='fasta', retmode = 'text')
file = open('TRAF3.fa', 'w')
file.writelines(handle)
file.close()
handle.close()

handle = Entrez.efetch(db='Protein', id='Q13114', rettype='gb', retmode = 'text')
file = open('TRAF3.gb', 'w')
file.writelines(handle)
file.close()
handle.close()

# MAVS - aminoacid sequence       Uniprot id = Q7Z434
handle = Entrez.efetch(db='Protein', id='Q7Z434', rettype='fasta', retmode = 'text')
file = open('MAVS.fa', 'w')
file.writelines(handle)
file.close()
handle.close()

handle = Entrez.efetch(db='Protein', id='Q7Z434', rettype='gb', retmode = 'text')
file = open('MAVS.gb', 'w')
file.writelines(handle)
file.close()
handle.close()

# ORF9b - SARS-CoV-2 aminoacid sequence       Uniprot id = P0DTD2
handle = Entrez.efetch(db='Protein',id='P0DTD2', rettype='fasta', retmode = 'text')
file = open('ORF9b_COVID19.fa', 'w')
file.writelines(handle)
file.close()
handle.close()

handle = Entrez.efetch(db='Protein', id='P0DTD2', rettype='gb', retmode = 'text')
file = open('ORF9b_COVID19.gb', 'w')
file.writelines(handle)
file.close()
handle.close()


# Annotations e Features
#   TRAF3
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

record = SeqIO.read("TRAF3.gb", "genbank")

print('id:', record.id)
print('designação:', record.description)
print('nome:', record.name)
print('tamanho da sequência:', len(record.seq))
print('referências a bases de dados (externas):', record.dbxrefs)
print('Número de anotações:', len(record.annotations))
print('Fonte de anotação:', record.annotations["source"] )

print('Número de features', len(record.features))
print('Taxonomia:', *record.annotations['taxonomy'])
print('Organismo:', record.annotations['organism'])

featregion = [i for i in range(len(record.features)) if record.features[i].type=='Region']

for  k in featregion:
    print (record.features[k].location)
    print(record.features[k].extract(record.seq))


#   MAVS
record1 = SeqIO.read("MAVS.gb", "genbank")
print()
print('id:', record1.id)
print('designação:', record1.description)
print('nome:', record1.name)
print('tamanho da sequência:', len(record1.seq))
print('referências a bases de dados (externas):', record1.dbxrefs)
print('Número de anotações:', len(record1.annotations))
print('Fonte de anotação:', record1.annotations["source"] )

print('Número de features', len(record1.features))
print('Taxonomia:', *record1.annotations['taxonomy'])
print('Organismo:', record1.annotations['organism'])

featregion1 = [i for i in range(len(record1.features)) if record1.features[i].type=='Region']

for k in featregion1:
    print (record1.features[k].location)
    print(record1.features[k].extract(record1.seq))


#   ORF9b
record2 = SeqIO.read("ORF9b_COVID19.gb", "genbank")
print()
print('id:', record2.id)
print('designação:', record2.description)
print('nome:', record2.name)
print('tamanho da sequência:', len(record2.seq))
print('referências a bases de dados (externas):', record2.dbxrefs)
print('Número de anotações:', len(record2.annotations))
print('Fonte de anotação:', record2.annotations["source"] )

print('Número de features', len(record2.features))
print('Taxonomia:', *record2.annotations['taxonomy'])
print('Organismo:', record2.annotations['organism'])

featregion2 = [i for i in range(len(record2.features)) if record2.features[i].type=='Region']

for k in featregion2:
    print (record2.features[k].location)
    print(record2.features[k].extract(record2.seq))
    print(record2.features[k].sub_features)


#Fazer os BLASTs!
#   TRAF3
TRAF3_blast = SeqIO.read(open('TRAF3.fa'), format ='fasta')
TRAF3_result = NCBIWWW.qblast('blastp', 'nr', TRAF3_blast.format('fasta'))
save_TRAF3_blast = open('TRAF3_BLAST.xml', 'w')
save_TRAF3_blast.writelines(TRAF3_result.read())
save_TRAF3_blast.close()
TRAF3_result.close()

#   MAVS
MAVS_blast = SeqIO.read(open('MAVS.fa'), format ='fasta')
MAVS_result = NCBIWWW.qblast('blastp', 'nr', MAVS_blast.format('fasta'))
save_MAVS_blast = open('MAVS_BLAST.xml', 'w')
save_MAVS_blast.writelines(MAVS_result.read())
save_MAVS_blast.close()
MAVS_result.close()

#   ORF9b
ORF9b_blast = SeqIO.read(open('ORF9b_COVID19.fa'), format ='fasta')
ORF9b_result = NCBIWWW.qblast('blastp', 'nr', ORF9b_blast.format('fasta'))
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

