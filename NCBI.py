from Bio import Entrez

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

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SeqFeature
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
print('Record Annotation:' , record.letter_annotations)
print('Número de features', len(record.features))
print('Taxonomia:', *record.annotations['taxonomy'])
print('Organismo:', record.annotations['organism'])

featregion = [i for i in range(len(record.features)) if record.features[i].type=='Region']
seqs_location = [print (record.features[k].location) for k in featregion]
seqs_features = [print(record.features[k].extract(record.seq)) for k in featregion]


record1 = SeqIO.read("MAVS.gb", "genbank")
print()
print('id:', record1.id)
print('designação:', record1.description)
print('nome:', record1.name)
print('tamanho da sequência:', len(record1.seq))
print('referências a bases de dados (externas):', record1.dbxrefs)
print('Número de anotações:', len(record1.annotations))
print('Fonte de anotação:', record1.annotations["source"] )
print(record1.letter_annotations)
print('Número de features', len(record1.features))
print('Taxonomia:', *record1.annotations['taxonomy'])
print('Organismo:', record1.annotations['organism'])

featregion1 = [i for i in range(len(record1.features)) if record1.features[i].type=='Region']
seqs_location1 = [print (record1.features[k].location) for k in featregion1]
seqs_features1 = [print(record1.features[k].extract(record1.seq)) for k in featregion1]


record2 = SeqIO.read("ORF9b_COVID19.gb", "genbank")

print()
print('id:', record2.id)
print('designação:', record2.description)
print('nome:', record2.name)
print('tamanho da sequência:', len(record2.seq))
print('referências a bases de dados (externas):', record2.dbxrefs)
print('Número de anotações:', len(record2.annotations))
print('Fonte de anotação:', record2.annotations["source"] )
print(record2.letter_annotations)
print('Número de features', len(record2.features))
print('Taxonomia:', *record2.annotations['taxonomy'])
print('Organismo:', record2.annotations['organism'])

featregion2 = [i for i in range(len(record2.features)) if record2.features[i].type=='Region']
seqs_location2 = [print (record2.features[k].location) for k in featregion2]
seqs_features2 = [print(record2.features[k].extract(record2.seq)) for k in featregion2]


#Fazer os BLASTs!
TRAF3_blast = SeqIO.read(open('TRAF3.fa'), format ='fasta')
TRAF3_result = NCBIWWW.qblast('blastp', 'Protein', TRAF3_blast.format('fasta'), hitlist_size=25)
save_TRAF3_blast = open('TRAF3_BLAST.xml', 'w')
save_TRAF3_blast.write(TRAF3_result.read())
save_TRAF3_blast.close()
TRAF3_result.close()

MAVS_blast = SeqIO.read(open('MAVS.fa'), format ='fasta')
MAVS_result = NCBIWWW.qblast('blastp', 'Protein', MAVS_blast.format('fasta'), hitlist_size=25)
save_MAVS_blast = open('MAVS_BLAST.xml', 'w')
save_MAVS_blast.write(MAVS_result.read())
save_MAVS_blast.close()
MAVS_result.close()

ORF9b_blast = SeqIO.read(open('ORF9b_COVID19.fa'), format ='fasta')
ORF9b_result = NCBIWWW.qblast('blastp', 'Protein', ORF9b_blast.format('fasta'), hitlist_size=25)
save_ORF9b_blast = open('ORF9b_BLAST.xml', 'w')
save_ORF9b_blast.write(ORF9b_result.read())
save_ORF9b_blast.close()
ORF9b_result.close()




