from chython import smiles
from chython.algorithms.groups._roles import roles, transformers
from chython.reactor import Transformer


def test_roles_is_plain_dict_of_pairs():
    assert isinstance(roles, dict)
    assert roles, 'glossary must not be empty'
    for role_name, entries in roles.items():
        assert isinstance(role_name, str)
        for entry in entries:
            fg_name, cap_template = entry  # exactly a 2-tuple
            assert isinstance(fg_name, str)
            assert isinstance(cap_template, str)


def test_transformers_compiled_per_entry():
    # transformers Proxy mirrors the roles structure, values are Transformers
    for role_name, entries in roles.items():
        compiled = transformers[role_name]
        assert len(compiled) == len(entries)
        for (fg_name, _cap), (c_fg_name, t) in zip(entries, compiled):
            assert fg_name == c_fg_name
            assert isinstance(t, Transformer)


def test_halide_roles_exclude_fluoride():
    # F is its own role, never grouped with Cl/Br/I
    for role in ('aryl_halide', 'alkyl_halide', 'alkenyl_halide', 'alkynyl_halide'):
        fg_names = {fg for fg, _ in roles[role]}
        assert not any('fluoride' in fg for fg in fg_names)


def test_every_entry_fires_on_a_real_molecule():
    # functional spot-check: each compiled transformer must produce at least one
    # capped product on a representative substrate (guards against cap templates
    # referencing unmapped atoms, e.g. the alkynyl_bromide :2 bug).
    substrates = {
        'aryl_chloride': 'Clc1ccccc1', 'aryl_bromide': 'Brc1ccccc1', 'aryl_iodide': 'Ic1ccccc1',
        'alkyl_chloride': 'ClCCC', 'alkyl_bromide': 'BrCCC', 'alkyl_iodide': 'ICCC',
        'alkenyl_chloride': 'ClC=CC', 'alkenyl_bromide': 'BrC=CC', 'alkenyl_iodide': 'IC=CC',
        'alkynyl_chloride': 'ClC#CC', 'alkynyl_bromide': 'BrC#CC', 'alkynyl_iodide': 'IC#CC',
        'aryl_fluoride': 'Fc1ccccc1', 'alkyl_fluoride': 'FCCC',
        'alkenyl_fluoride': 'FC=CC', 'alkynyl_fluoride': 'FC#CC',
        'aryl_triflate': 'O=S(=O)(Oc1ccccc1)C(F)(F)F', 'aryl_mesylate': 'O=S(=O)(Oc1ccccc1)C',
        'aryl_tosylate': 'O=S(=O)(Oc1ccccc1)c1ccc(C)cc1', 'alkyl_triflate': 'O=S(=O)(OCCC)C(F)(F)F',
        'alkyl_mesylate': 'O=S(=O)(OCCC)C', 'alkyl_tosylate': 'O=S(=O)(OCCC)c1ccc(C)cc1',
        'aryl_boronic_acid': 'OB(O)c1ccccc1', 'aryl_boronic_ester': 'CC1(C)OB(c2ccccc2)OC1(C)C',
        'aryl_molander_salt': '[B-](F)(F)(F)c1ccccc1',
        'alkyl_boronic_acid': 'OB(O)CCC', 'alkyl_boronic_ester': 'CC1(C)OB(CCC)OC1(C)C',
        'alkyl_molander_salt': '[B-](F)(F)(F)CCC',
        'alkenyl_boronic_acid': 'OB(O)C=CC', 'alkenyl_boronic_ester': 'CC1(C)OB(C=CC)OC1(C)C',
        'alkenyl_molander_salt': '[B-](F)(F)(F)C=CC',
        'alkynyl_boronic_acid': 'OB(O)C#CC', 'alkynyl_boronic_ester': 'CC1(C)OB(C#CC)OC1(C)C',
        'alkynyl_molander_salt': '[B-](F)(F)(F)C#CC',
        'aryl_grignard': '[Mg](Br)c1ccccc1', 'alkyl_grignard': '[Mg](Br)CCC', 'alkenyl_grignard': '[Mg](Br)C=CC',
        'aryl_zinc': '[Zn](Br)c1ccccc1', 'alkyl_zinc': '[Zn](Br)CCC', 'alkenyl_zinc': '[Zn](Br)C=CC',
        'aryl_stannane': '[Sn](C)(C)(C)c1ccccc1', 'alkyl_stannane': '[Sn](C)(C)(C)CCC',
        'alkenyl_stannane': '[Sn](C)(C)(C)C=CC',
        'aryl_silane': '[Si](C)(C)(C)c1ccccc1', 'alkenyl_silane': '[Si](C)(C)(C)C=CC',
        'alkynyl_silane': '[Si](C)(C)(C)C#Cc1ccccc1',
        'alkyl_carboxylic_acid': 'OC(=O)CCC', 'aryl_carboxylic_acid': 'OC(=O)c1ccccc1',
        # cinnamic acid: beta-vinyl carbon (:4) carries a phenyl, so a cap that
        # drops :4 would visibly truncate the ring (regression guard).
        'alkenyl_carboxylic_acid': 'OC(=O)C=Cc1ccccc1', 'alkynyl_carboxylic_acid': 'OC(=O)C#CC',
        'acyl_chloride': 'ClC(=O)CCC', 'acyl_bromide': 'BrC(=O)CCC', 'acyl_fluoride': 'FC(=O)CCC',
        'primary_amine': 'NCCC', 'secondary_amine': 'CNCC',
        'primary_aniline': 'Nc1ccccc1', 'secondary_aniline': 'CNc1ccccc1',
        'primary_alcohol': 'OCCC', 'secondary_alcohol': 'OC(C)CC', 'tertiary_alcohol': 'OC(C)(C)CC',
        'phenol': 'Oc1ccccc1', 'carboxylic_acid': 'OC(=O)CCC', 'thiol': 'SCCC',
    }
    for role_name, entries in transformers.items():
        for fg_name, t in entries:
            mol = smiles(substrates[fg_name])
            products = list(t(mol))
            assert products, f'{role_name}/{fg_name} produced no capped product'
            assert any(a.atomic_symbol == 'At' for p in products for _, a in p.atoms()), \
                f'{role_name}/{fg_name} produced no astatine cap'


def test_alkenyl_acid_roles_keep_the_vinyl_arm():
    # regression: alkenyl_carboxylic_acid maps the beta-vinyl carbon as :4; caps
    # omitting :4 silently delete it and everything beyond (e.g. cinnamic acid's
    # phenyl). The capped fragment must still contain the aromatic ring.
    cinnamic = smiles('OC(=O)C=Cc1ccccc1')
    for role in ('alkenyl_acyl', 'alkenyl_decarboxy'):
        products = [p for _fg, t in transformers[role] for p in t(cinnamic)]
        assert products, f'{role} produced no product'
        for p in products:
            aromatic = sum(1 for _, a in p.atoms() if a.atomic_symbol == 'C' and a.hybridization == 4)
            assert aromatic == 6, f'{role} lost the phenyl ring: {p}'
