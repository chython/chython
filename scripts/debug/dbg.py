from chython import smiles


r = smiles('Nc1ccc(Br)cc1>>CC(C)(C)OC(=O)Nc1ccc(Br)cc1')
r.canonicalize()
r.reconstruct_mapping()

print(format(r, 'm'))
