import util.tools as tls

print tls.eval_until_licence_found("pdbconvert -brief -imae dock_sorted.mae -opdb dock_sorted.pdb 2> /dev/null")
