from chython import MoleculeContainer, ReactionContainer, smiles
from rdkit.Chem import MolFromSmiles, MolToSmiles
import zipfile


# convert zipped PaChs to smiles and store in plain text file
with zipfile.ZipFile('SI.zip') as z, open('smiles.txt', 'w') as o:
    for x in z.filelist:
        if x.filename.endswith('pach'):
            m = MoleculeContainer.unpack(z.read(x))
            print(m)
            o.write(str(m))
            o.write('\n')


# reaction example
r = smiles('CC=O>>CCO')
p = r.pack()
u = ReactionContainer.unpack(p)
print(r == u)


# generate uncompressed PaCh
m = smiles('CCN')
p = m.pack(compressed=False)
u = MoleculeContainer.unpack(p, compressed=False)
print(m == u)


def MolFromPach(data, compressed=True):
    """
    parse pach to RDKit molecule
    """
    m = MoleculeContainer.unpack(data, compressed=compressed)
    return MolFromSmiles(str(m))


print(MolToSmiles(MolFromPach(smiles('CCO').pack())))
