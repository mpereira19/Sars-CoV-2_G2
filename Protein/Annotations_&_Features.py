# Annotations e Features
#   TRAF3
from Bio import SeqIO
orfs = ['TRAF3', 'MAVS', 'ORF9b']

for orf in orfs:
    record = SeqIO.read(orf + ".gb", "genbank")
    print('ID:', record.id)
    print('Nome:', record.name)
    print('Descrição:', record.description)
    print(f'Tamanho da sequência:{len(record.seq)} aa')
    print('Referências a bases de dados (externas):', record.dbxrefs)
    print('Número de anotações:', len(record.annotations))
    print('Fonte de anotação:', record.annotations["source"])

    print('Número de features', len(record.features))
    print('Taxonomia:', *record.annotations['taxonomy'])
    print('Organismo:', record.annotations['organism'])
    print('Comentário:')
    print(record.annotations['comment'])

    featregion = [i for i in range(len(record.features)) if record.features[i].type == 'Region']
    print()
    print('Features:')
    print()

    for k in featregion:
        print(*record.features[k].qualifiers['region_name'])
        print(record.features[k].location)
        print(record.features[k].extract(record.seq))
        print(record.features[k].qualifiers['note'])
        if k != featregion[-1]:
            print('---------')
    print()
    if orf != orfs[-1]:
        print('############################################')
    print()