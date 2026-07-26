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


def test_leaving_halogen_excludes_fluoride():
    # F is its own role, never grouped with Cl/Br/I
    fg_names = {fg for fg, _ in roles.get('leaving_halogen', [])}
    assert 'aryl_fluoride' not in fg_names
    assert 'alkyl_fluoride' not in fg_names
