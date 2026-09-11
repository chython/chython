# Scratch generator for test/wedge_stereo.sdf -- public / textbook stereochemistry only.
# Run once with rdkit available; the SDF it writes is what the test suite reads.
from rdkit import Chem
from rdkit.Chem import AllChem

CORPUS = [
 ('menthol', 'C[C@@H]1CC[C@H](C(C)C)[C@@H](O)C1'),
 ('neomenthol', 'C[C@@H]1CC[C@H](C(C)C)[C@H](O)C1'),
 ('camphor', 'CC1(C)[C@@H]2CC[C@]1(C)C(=O)C2'),
 ('borneol', 'CC1(C)[C@@H]2CC[C@]1(C)[C@@H](O)C2'),
 ('isoborneol', 'CC1(C)[C@@H]2CC[C@]1(C)[C@H](O)C2'),
 ('camphene', 'CC1(C)[C@@H]2CC[C@H]1CC2=C'),
 ('alpha-pinene', 'CC1=CC[C@H]2C[C@@H]1C2(C)C'),
 ('beta-pinene', 'C=C1CC[C@H]2C[C@@H]1C2(C)C'),
 ('fenchone', 'CC1(C)[C@@H]2CC[C@](C)(C2)C1=O'),
 ('carvone', 'CC(=C)[C@@H]1CC=C(C)C(=O)C1'),
 ('limonene', 'CC(=C)[C@@H]1CCC(C)=CC1'),
 ('D-glucose', 'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O'),
 ('beta-D-glucose', 'OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O'),
 ('D-galactose', 'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O'),
 ('D-mannose', 'OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O'),
 ('D-fructose', 'OC[C@H]1O[C@](O)(CO)[C@@H](O)[C@@H]1O'),
 ('D-ribose', 'OC[C@H]1O[C@@H](O)[C@H](O)[C@H]1O'),
 ('sucrose', 'OC[C@H]1O[C@@](CO)(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O'),
 ('lactose', 'OC[C@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O'),
 ('glucosamine', 'OC[C@H]1O[C@@H](O)[C@H](N)[C@@H](O)[C@@H]1O'),
 ('L-tartaric acid', 'O[C@@H](C(=O)O)[C@H](O)C(=O)O'),
 ('meso-tartaric acid', 'O[C@@H](C(=O)O)[C@@H](O)C(=O)O'),
 ('L-threonine', 'C[C@@H](O)[C@H](N)C(=O)O'),
 ('L-isoleucine', 'CC[C@H](C)[C@H](N)C(=O)O'),
 ('L-alanine', 'C[C@H](N)C(=O)O'),
 ('L-serine', 'OC[C@H](N)C(=O)O'),
 ('L-cysteine', 'SC[C@H](N)C(=O)O'),
 ('L-proline', 'OC(=O)[C@@H]1CCCN1'),
 ('L-tryptophan', 'N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O'),
 ('L-histidine', 'N[C@@H](Cc1c[nH]cn1)C(=O)O'),
 ('cholesterol', 'C[C@]12CC[C@H]3[C@@H](CC=C1C[C@@H](O)CC2)[C@@]1(C)CC[C@@H]([C@H](C)CCCC(C)C)[C@@H]1CC3'),
 ('cholic acid', 'C[C@H](CCC(=O)O)[C@H]1CC[C@H]2[C@@H]1[C@H](O)C[C@H]1[C@@]2(C)[C@H](O)C[C@H]2[C@@]1(C)CC[C@H](O)C2'),
 ('testosterone', 'C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O'),
 ('progesterone', 'CC(=O)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C'),
 ('estradiol', 'C[C@]12CC[C@H]3c4ccc(O)cc4CC[C@H]3[C@@H]1CC[C@@H]2O'),
 ('cortisol', 'C[C@]12C[C@H](O)[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@]2(O)C(=O)CO'),
 ('androsterone', 'C[C@]12CC[C@H]3[C@@H](CC[C@H]4CC(=O)CC[C@]34C)[C@@H]1CCC2=O'),
 ('morphine', 'CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5'),
 ('quinine', 'C=C[C@H]1CN2CC[C@H]1C[C@H]2[C@@H](O)c1ccnc2ccc(OC)cc12'),
 ('quinidine', 'C=C[C@H]1CN2CC[C@H]1C[C@@H]2[C@H](O)c1ccnc2ccc(OC)cc12'),
 ('atropine', 'CN1[C@H]2CC[C@@H]1C[C@@H](C2)OC(=O)C(CO)c1ccccc1'),
 ('nicotine', 'CN1CCC[C@H]1c1cccnc1'),
 ('pseudoephedrine', 'C[C@H](NC)[C@@H](O)c1ccccc1'),
 ('adrenaline', 'CNC[C@H](O)c1ccc(O)c(O)c1'),
 ('penicillin-G', 'CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O'),
 ('cephalosporin-core', 'CC1=C(C(=O)O)N2C(=O)[C@@H](N)[C@H]2SC1'),
 ('ascorbic acid', 'OC[C@H](O)[C@H]1OC(=O)C(O)=C1O'),
 ('biotin', 'O=C1N[C@@H]2CS[C@@H](CCCCC(=O)O)[C@@H]2N1'),
 ('camptothecin', 'CC[C@@]1(O)C(=O)OCc2c1cc1n(c2=O)Cc2cc3ccccc3nc21'),
 ('artemisinin', 'C[C@@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@]4(C)CC[C@@H]1[C@]23OO4'),
 ('norbornane-diol', 'O[C@H]1C[C@@H]2CC[C@H]1C2'),
 ('norbornanol', 'O[C@H]1C[C@@H]2CC[C@H]1CC2'),
 ('bicyclo221-amine', 'N[C@H]1C[C@@H]2CC[C@H]1C2'),
 ('trans-decalol', 'O[C@H]1CC[C@H]2CCCC[C@H]2C1'),
 ('cis-decalol', 'O[C@H]1CC[C@@H]2CCCC[C@H]2C1'),
 ('hydrindane', 'C[C@H]1CC[C@H]2CCCC[C@H]12'),
 ('carvomenthol', 'CC(C)[C@H]1CC[C@@H](C)C[C@H]1O'),
 ('isopulegol', 'CC(=C)[C@@H]1CC[C@@H](C)C[C@H]1O'),
 ('trans-4-methylcyclohexanol', 'C[C@H]1CC[C@@H](O)CC1'),
 ('cis-2-methylcyclohexanol', 'C[C@H]1CCCC[C@@H]1O'),
 ('shikimic acid', 'O[C@@H]1C=C(C(=O)O)C[C@H](O)[C@H]1O'),
 ('quinic acid', 'O[C@H]1C[C@](O)(C(=O)O)C[C@H](O)[C@H]1O'),
 ('inositol', 'O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O'),
 ('sphingosine-core', 'CCCCCCCCCCCCCC=C[C@@H](O)[C@@H](N)CO'),
 ('ibuprofen-S', 'CC(C)Cc1ccc([C@@H](C)C(=O)O)cc1'),
 ('naproxen-S', 'COc1ccc2cc([C@@H](C)C(=O)O)ccc2c1'),
 ('thalidomide-S', 'O=C1CC[C@H](N2C(=O)c3ccccc3C2=O)C(=O)N1'),
 ('warfarin-S', 'C[C@H](CC(=O)c1ccccc1)c1c(O)c2ccccc2oc1=O'),
 ('captopril', 'C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O'),
 ('oseltamivir-core', 'CCOC(=O)C1=C[C@@H](N)[C@H](NC(C)=O)[C@@H](OCC)C1'),
 ('chloramphenicol', 'OC[C@H](NC(=O)C(Cl)Cl)[C@H](O)c1ccc([N+](=O)[O-])cc1'),
 ('pantothenic acid', 'OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)O'),
 ('lactic acid', 'C[C@H](O)C(=O)O'),
 ('glyceraldehyde', 'OC[C@@H](O)C=O'),
 ('malic acid', 'O[C@@H](CC(=O)O)C(=O)O'),
 ('citramalic acid', 'C[C@](O)(CC(=O)O)C(=O)O'),
 ('mandelic acid', 'O[C@H](C(=O)O)c1ccccc1'),
 ('trans-cyclohexanediol', 'O[C@H]1CCCC[C@@H]1O'),
 ('cis-cyclohexanediol', 'O[C@H]1CCCC[C@H]1O'),
 ('trans-stilbene-oxide', 'c1ccc([C@H]2O[C@@H]2c2ccccc2)cc1'),
 ('cyclopropane-acid', 'O=C(O)[C@H]1C[C@@H]1c1ccccc1'),
 ('aziridine-sub', 'C[C@H]1N[C@@H]1c1ccccc1'),
 ('longifolene-frag', 'CC1(C)CC[C@H]2CC[C@]1(C)C2'),
 ('caryophyllene-frag', 'CC1(C)C[C@@H]2CC[C@H]1CC2'),
 ('steviol-frag', 'C[C@]12CC[C@H](C1)[C@@H]1CC[C@H]2C1'),
 ('gibberellin-frag', 'O=C1O[C@@H]2CC[C@H]3CC[C@]1(C3)C2'),
 ('cedrol-frag', 'C[C@]12CC[C@@H](C)[C@H]1CC[C@H]2O'),
 ('brucine-frag', 'COc1ccc2c(c1)N1C(=O)C[C@H]3OCC=C4CN5CC[C@]2(C1)[C@H]5C[C@H]34'),
 ('cyclitol', 'O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O'),
 ('bridged-lactone', 'O=C1O[C@H]2C[C@@H]1CC2'),
 ('oxanorbornane', 'O1[C@@H]2CC[C@H]1[C@H]1CC[C@@H]21'),
 ('tropane-diol', 'CN1[C@H]2C[C@H](O)[C@@H](O)C[C@@H]1CC2'),
]


def main(path):
    w = Chem.SDWriter(path)
    kept = dropped = 0
    for name, smi in CORPUS:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            print('DROPPED (parse)', name)
            dropped += 1
            continue
        Chem.AssignStereochemistry(m, cleanIt=True, force=True)
        AllChem.Compute2DCoords(m)
        Chem.WedgeMolBonds(m, m.GetConformer())
        m.SetProp('_Name', name)
        w.write(m)
        kept += 1
    w.close()
    print('wrote', kept, 'dropped', dropped)


if __name__ == '__main__':
    main('test/wedge_stereo.sdf')
