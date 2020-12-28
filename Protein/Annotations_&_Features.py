# Annotations e Features
#   TRAF3
from Bio import SeqIO


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
    print (record.features[k])
    print(record.features[k].extract(record.seq))
    print('---------')

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
    print (record1.features[k])
    print(record1.features[k].extract(record1.seq))
    print('---------')

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
    print (record2.features[k])
    print(record2.features[k].extract(record2.seq))
    print('---------')
