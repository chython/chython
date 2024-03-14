from pandas import read_csv
from time import time
from chython import smiles, SDFWrite, mdl_mol, MoleculeContainer
from rdkit.Chem import MolFromSmiles, MolToSmiles, MolFromMolBlock, MolToMolBlock
from io import StringIO
from pickle import dumps, loads

data = read_csv('./lipophilicity.csv')

t = time()
chms = [smiles(x) for x in data.smiles]
print('chython smiles parsing', (time() - t) * 1000)

t = time()
rdms = [MolFromSmiles(x) for x in data.smiles]
print('rdkit smiles parsing', (time() - t) * 1000)

t = time()
for x in chms:
    str(x)
print('chython to smiles', (time() - t) * 1000)

t = time()
for x in rdms:
    MolToSmiles(x)
print('rdkit to smiles', (time() - t) * 1000)

t = time()
for x in chms:
    SDFWrite(StringIO()).write(x)
print('chython to mol', (time() - t) * 1000)

t = time()
mblk = [MolToMolBlock(x) for x in rdms]
print('rdkit to mol', (time() - t) * 1000)

t = time()
for x in mblk:
    mdl_mol(x)
print('chython mol parsing', (time() - t) * 1000)

t = time()
for x in mblk:
    MolFromMolBlock(x)
print('rdkit mol parsing', (time() - t) * 1000)

t = time()
chpk = [dumps(x) for x in chms]
print('chython to pickle', (time() - t) * 1000)

t = time()
rdpk = [dumps(x) for x in rdms]
print('rdkit to pickle', (time() - t) * 1000)

t = time()
for x in chpk:
    loads(x)
print('chython pickle parsing', (time() - t) * 1000)

t = time()
for x in rdpk:
    loads(x)
print('rdkit pickle parsing', (time() - t) * 1000)

t = time()
pach = [x.pack() for x in chms]
print('chython to PaCh zlib', (time() - t) * 1000)

t = time()
pachu = [x.pack(compressed=False) for x in chms]
print('chython to PaCh', (time() - t) * 1000)

t = time()
for x in pach:
    MoleculeContainer.unpack(x)
print('chython PaCh zlib parsing', (time() - t) * 1000)

t = time()
for x in pachu:
    MoleculeContainer.unpack(x, compressed=False)
print('chython PaCh parsing', (time() - t) * 1000)

print('smiles size', sum(len(x) for x in data.smiles) / 1024)
print('mol size', sum(len(x) for x in mblk) / 1024)
print('rdkit pickle size', sum(len(x) for x in rdpk) / 1024)
print('chython pickle size', sum(len(x) for x in chpk) / 1024)
print('pack compressed', sum(len(x) for x in pach) / 1024)
print('pack', sum(len(x) for x in pachu) / 1024)
