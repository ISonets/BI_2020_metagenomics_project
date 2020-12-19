from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle
from Bio import SeqIO
import re

names = [rec.name for rec in SeqIO.parse('/home/ivan/Metagenom/Tree/set/all_fungi.fasta','fasta')]
t = Tree('/home/ivan/Metagenom/Tree/Tree_gblock/regions_tree.tre', format =0, quoted_node_names=True)

for node in t.traverse():
    nodestyle = NodeStyle()
    nodestyle["fgcolor"] = "green"
    nodestyle["size"] = 2
    node.set_style(nodestyle)

    style1 = NodeStyle()
    style1['bgcolor'] = "LightSteelBlue"
    style2 = NodeStyle()
    style2["bgcolor"] = "DarkSeaGreen"

    # run through the leaves
    for l in t.iter_leaves():
        if str(l)[3:] in set(names):
            l.img_style = style1
            l.scale = 1
        else:
            l.img_style = style2
# arrange the tree from bottom to top
t.ladderize(direction=1)
ts = TreeStyle()
ts.show_branch_support = True
t.show(tree_style=ts)
# save
t.render('/home/ivan/Metagenom/Tree/Tree_gblock/tree_18S_gblock.pdf', dpi=1200, tree_style=ts)