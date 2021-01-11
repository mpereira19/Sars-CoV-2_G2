#from Bio.Align.Applications import ClustalOmegaCommandline
#from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Phylo

orfs = ['TRAF3', 'MAVS', 'ORF9b']

# ficheiros '..._alignment.clustal_num' e '..._tree.dnd' foram obtidos através de Clustal Omega online!

for orf in orfs:
    filename = orf + '_alignment.aln'
    format = 'clustal'
    align = AlignIO.read(orf + '_alignment.clustal_num', 'clustal')
    print(align)


for orf in orfs:
    tree = Phylo.read(orf + '_tree.dnd', 'newick')
    Phylo.draw_ascii(tree)
