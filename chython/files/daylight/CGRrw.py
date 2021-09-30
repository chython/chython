from collections import defaultdict
from itertools import permutations
from re import split, compile, fullmatch, findall, search
from ._parser import DaylightParser
from ...exceptions import IncorrectSmiles
from ...containers import CGRContainer
from ...containers.bonds import DynamicBond
from ...periodictable import DynamicElement


dyn_charge_dict = {'-4': -4, '-3': -3, '-2': -2, '--': -2, '-': -1, '0': 0, '+': 1, '++': 2, '+2': 2, '+3': 3, '+4': 4}
tmp = {f'{i}>{j}': (x, y) for (i, x), (j, y) in permutations(dyn_charge_dict.items(), 2) if x != y}
dyn_charge_dict = {k: (v, v) for k, v in dyn_charge_dict.items()}
dyn_charge_dict.update(tmp)
dyn_radical_dict = {'*': (True, True), '*>^': (True, False), '^>*': (False, True)}

dyn_atom_re = compile(r'([1-9][0-9]{0,2})?([A-IK-PR-Zacnopsbt][a-ik-pr-vy]?)([+-0][1-4+-]?(>[+-0][1-4+-]?)?)?'
                      r'([*^](>[*^])?)?')


class CGRRead(DaylightParser):
    def parse(self, string: str):
        ...

    @staticmethod
    def __dynatom_parse(token):
        # [isotope]Element[element][+-charge[>+-charge]][*^[>*^]]
        match = fullmatch(dyn_atom_re, token)
        if match is None:
            raise IncorrectSmiles(f'atom token invalid {{{token}}}')
        isotope, element, charge, _, is_radical, _ = match.groups()

        if isotope:
            isotope = int(isotope)

        if charge:
            try:
                charge, p_charge = dyn_charge_dict[charge]
            except KeyError:
                raise IncorrectSmiles('charge token invalid')
        else:
            charge = p_charge = 0

        if is_radical:
            try:
                is_radical, p_is_radical = dyn_radical_dict[is_radical]
            except KeyError:
                raise IncorrectSmiles('invalid dynamic radical token')
        else:
            is_radical = p_is_radical = False

        if element in ('c', 'n', 'o', 'p', 's', 'as', 'se', 'b', 'te'):
            _type = 12
            element = element.capitalize()
        else:
            _type = 11
        return _type, {'element': element, 'charge': charge, 'isotope': isotope, 'is_radical': is_radical,
                       'p_charge': p_charge, 'p_is_radical': p_is_radical}

    def _convert_cgr(self, data):
        atoms = data['atoms']
        bonds = defaultdict(dict)

        for n, m, value in data['cgr']:
            bonds[n][m] = bonds[m][n] = DynamicBond(*value)

        g = object.__new__(CGRContainer)
        g_atoms = {}
        g_bonds = {}
        plane = {}
        charges = {}
        radicals = {}
        p_charges = {}
        p_radicals = {}
        for n, atom in enumerate(atoms, 1):
            g_atoms[n] = DynamicElement.from_symbol(atom['element'])(atom['isotope'])
            g_bonds[n] = {}
            charges[n] = atom['charge']
            radicals[n] = atom['is_radical']
            if 'p_charge' in atom:
                p_charges[n] = atom['p_charge']
                p_radicals[n] = atom['p_is_radical']
            else:
                p_charges[n] = atom['charge']
                p_radicals[n] = atom['is_radical']
            plane[n] = (0., 0.)
        for n, m, b in data['bonds']:
            if m in bonds[n]:
                if b != 8:
                    raise ValueError('CGR spec invalid')
                b = bonds[n][m]
            else:
                b = DynamicBond(b, b)
            n += 1
            m += 1
            if n in g_bonds[m]:
                raise ValueError('atoms already bonded')
            g_bonds[n][m] = g_bonds[m][n] = b
        g.__setstate__({'atoms': g_atoms, 'bonds': g_bonds, 'plane': plane, 'charges': charges, 'radicals': radicals,
                        'p_charges': p_charges, 'p_radicals': p_radicals, 'conformers': []})
        return g


__all__ = ['CGRRead']
