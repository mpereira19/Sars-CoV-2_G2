from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Phylo

orfs = ['TRAF3', 'MAVS', 'ORF9b']

for orf in orfs:
    in_file = orf + '_homology.fa'
    out_file = orf + '_alignment.aln'
    clustalw_cline = ClustalwCommandline(infile=in_file, outfile=out_file)
    print(clustalw_cline)
    clustalw_cline()


for orf in orfs:
    filename = orf + '_alignment.aln'
    format = 'clustal'
    align = AlignIO.read(orf + '_alignment.aln', 'clustal')
    print(align)


#for orf in orfs:
#    tree = Phylo.read(orf + '_homology.dnd', 'newick')
#    Phylo.draw_ascii(tree)
