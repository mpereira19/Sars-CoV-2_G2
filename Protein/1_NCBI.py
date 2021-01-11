from Bio import Entrez

Entrez.email = 'pg42875@alunos.uminho.pt'

orfs = ['TRAF3', 'MAVS', 'ORF9b']
id = ['Q13114', 'Q7Z434', 'P0DTD2']
# TRAF3 - aminoacid sequence                  Uniprot id = Q13114
# MAVS - aminoacid sequence                   Uniprot id = Q7Z434
# ORF9b - SARS-CoV-2 aminoacid sequence       Uniprot id = P0DTD2

for n in range(len(orfs)):
    handle = Entrez.efetch(db='Protein', id=id[n], rettype='fasta', retmode='text')
    file = open(orfs[n] + '.fa', 'w')
    file.writelines(handle)
    file.close()
    handle.close()

    handle = Entrez.efetch(db='Protein', id=id[n], rettype='gb', retmode='text')
    file = open(orfs[n] + '.gb', 'w')
    file.writelines(handle)
    file.close()
    handle.close()