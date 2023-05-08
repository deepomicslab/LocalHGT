from ete3 import Tree



# Load a tree structure from a newick file.
t = Tree("/mnt/d/HGT/time_lines/distribution/bac120_iqtree.nwk")
t.prune(("GUT_GENOME275558","GUT_GENOME212961","GUT_GENOME210710"))
print (t)