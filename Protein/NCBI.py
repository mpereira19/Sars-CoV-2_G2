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
