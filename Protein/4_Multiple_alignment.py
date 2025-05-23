# Multiple Alignment e Filogenia

from Bio import AlignIO
from Bio import Phylo

orfs = ['TRAF3', 'MAVS', 'ORF9b']

# ficheiros '..._alignment.clustal_num' e '..._tree.dnd' foram obtidos através de Clustal Omega online!

for orf in orfs:
    align = AlignIO.read(orf + '_alignment.clustal_num', 'clustal')
    print(align)

    tree = Phylo.read(orf + '_tree.dnd', 'newick')
    Phylo.draw_ascii(tree)
