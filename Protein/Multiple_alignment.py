from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import Phylo

orfs = ['TRAF3', 'MAVS', 'ORF9b']

for orf in orfs:
    in_file = orf + '_homology.fa'
    out_file = orf + '_alignment.fa'
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
    print(clustalomega_cline)


#for orf in orfs:
#    tree = Phylo.read(orf + '_alignment.fa', 'newick')
#    Phylo.draw_ascii(tree)
