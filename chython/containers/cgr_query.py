from typing import Union
from . import cgr  # cyclic imports resolve
from .bonds import Bond, DynamicBond
from .graph import Graph
from ..algorithms.calculate2d import Calculate2DCGR
from ..algorithms.depict import DepictQueryCGR
from ..algorithms.smiles import QueryCGRSmiles
from ..exceptions import IsConnectedAtom
from ..periodictable import (Element, DynamicElement, QueryElement, DynamicQueryElement, AnyElement, DynamicAnyElement,
                             AnyAtom)


class QueryCGRContainer(Graph, QueryCGRSmiles, DepictQueryCGR, Calculate2DCGR):
    __slots__ = ('_p_charges', '_p_radicals', '_neighbors', '_hybridizations', '_p_neighbors', '_p_hybridizations', '_plane', '_charges', '_radicals')

    def __init__(self):
        self._p_charges: dict[int, int] = {}
        self._p_radicals: dict[int, bool] = {}
        self._neighbors: dict[int, tuple[int, ...]] = {}
        self._hybridizations: dict[int, tuple[int, ...]] = {}
        self._p_neighbors: dict[int, tuple[int, ...]] = {}
        self._p_hybridizations: dict[int, tuple[int, ...]] = {}
        self._plane: dict[int, tuple[float, float]] = {}
        self._charges: dict[int, int] = {}
        self._radicals: dict[int, bool] = {}
        super().__init__()

    def _validate_charge(self, charge: int) -> int:
        if not isinstance(charge, int):
            raise TypeError(f"Charge must be an integer, got {type(charge)}")
        return charge

    def _validate_radical(self, is_radical: bool) -> bool:
        if not isinstance(is_radical, bool):
            raise TypeError(f"is_radical must be a boolean, got {type(is_radical)}")
        return is_radical

    def add_atom(self, atom: Union[AnyAtom, int, str], *args,
                 p_charge: int = 0, p_is_radical: bool = False,
                 neighbors: Union[int, list[int], tuple[int, ...], None] = None,
                 hybridization: Union[int, list[int], tuple[int, ...], None] = None,
                 p_neighbors: Union[int, list[int], tuple[int, ...], None] = None,
                 p_hybridization: Union[int, list[int], tuple[int, ...], None] = None, **kwargs):
        xy_coord = kwargs.pop('xy', None)
        neighbors = self._validate_neighbors(neighbors)
        p_neighbors = self._validate_neighbors(p_neighbors)
        hybridization = self._validate_hybridization(hybridization)
        p_hybridization = self._validate_hybridization(p_hybridization)
        neighbors, p_neighbors = self._validate_neighbors_pairing(neighbors, p_neighbors)
        hybridization, p_hybridization = self._validate_hybridization_pairing(hybridization, p_hybridization)

        p_charge = self._validate_charge(p_charge)
        p_is_radical = self._validate_radical(p_is_radical)

        if not isinstance(atom, (DynamicQueryElement, DynamicAnyElement)):
            if isinstance(atom, AnyElement):
                atom = DynamicAnyElement()
            elif isinstance(atom, (Element, QueryElement, DynamicElement)):
                atom = DynamicQueryElement.from_atomic_number(atom.atomic_number, isotope=atom.isotope)
            elif isinstance(atom, str):
                atom = DynamicQueryElement.from_symbol(atom, isotope=None)
            elif isinstance(atom, int):
                atom = DynamicQueryElement.from_atomic_number(atom, isotope=None)
            else:
                raise TypeError('QueryElement object expected')

        _map = super().add_atom(atom, *args, **kwargs)
        atom._attach_to_graph(self, _map)
        if xy_coord is not None:
            self._plane[_map] = xy_coord
            if hasattr(atom, 'xy'):
                atom.xy = xy_coord
        self._charges[_map] = getattr(atom, '_charge', 0)
        self._radicals[_map] = getattr(atom, '_is_radical', False)
        self._p_charges[_map] = p_charge
        self._p_radicals[_map] = p_is_radical
        self._neighbors[_map] = neighbors
        self._hybridizations[_map] = hybridization
        self._p_neighbors[_map] = p_neighbors
        self._p_hybridizations[_map] = p_hybridization
        return _map

    def add_bond(self, n, m, bond: Union[DynamicBond, Bond, int]):
        if isinstance(bond, Bond):
            bond = DynamicBond.from_bond(bond)
        elif not isinstance(bond, DynamicBond):
            bond = DynamicBond(bond, bond)
        super().add_bond(n, m, bond)

    def delete_atom(self, n):
        super().delete_atom(n)
        self._charges.pop(n, None)
        self._radicals.pop(n, None)
        del self._p_charges[n]
        del self._p_radicals[n]
        del self._neighbors[n]
        del self._hybridizations[n]
        del self._p_neighbors[n]
        del self._p_hybridizations[n]
        if n in self._plane:
            del self._plane[n]

    def remap(self, mapping, *, copy=False) -> 'QueryCGRContainer':
        h = super().remap(mapping, copy=copy)
        mg = mapping.get
        spr = self._p_radicals
        sn = self._neighbors
        sh = self._hybridizations
        spn = self._p_neighbors
        sph = self._p_hybridizations
        if copy:
            hch = getattr(h, '_charges', {})
            hrd = getattr(h, '_radicals', {})
            hpc = h._p_charges
            hpr = h._p_radicals
            hn = h._neighbors
            hh = h._hybridizations
            hpn = h._p_neighbors
            hph = h._p_hybridizations
        else:
            hch = {}
            hrd = {}
            hpc = {}
            hpr = {}
            hn = {}
            hh = {}
            hpn = {}
            hph = {}

        for n, c in self._p_charges.items():
            m = mg(n, n)
            hch[m] = self._charges.get(n, 0)
            hrd[m] = self._radicals.get(n, False)
            hpc[m] = c
            hpr[m] = spr[n]
            hn[m] = sn[n]
            hpn[m] = spn[n]
            hh[m] = sh[n]
            hph[m] = sph[n]

        if copy:
            return h

        self._charges = hch
        self._radicals = hrd
        self._p_charges = hpc
        self._p_radicals = hpr
        self._neighbors = hn
        self._hybridizations = hh
        self._p_neighbors = hpn
        self._p_hybridizations = hph
        return self

    def copy(self, **kwargs) -> 'QueryCGRContainer':
        copy = super().copy(**kwargs)
        copy._neighbors = self._neighbors.copy()
        copy._hybridizations = self._hybridizations.copy()
        copy._p_neighbors = self._p_neighbors.copy()
        copy._p_hybridizations = self._p_hybridizations.copy()
        copy._p_radicals = self._p_radicals.copy()
        copy._p_charges = self._p_charges.copy()
        copy._plane = self._plane.copy()
        copy._charges = self._charges.copy()
        copy._radicals = self._radicals.copy()
        for n, atom in copy._atoms.items():
            try:
                atom._attach_to_graph(copy, n)
            except (AttributeError, IsConnectedAtom):
                pass
            if n in copy._plane and hasattr(atom, 'xy'):
                atom.xy = copy._plane[n]
        return copy

    def substructure(self, atoms, **kwargs) -> 'QueryCGRContainer':
        """
        create substructure containing atoms from atoms list

        :param atoms: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        """
        atoms_set = {n for n in atoms if n in self._atoms}
        sub = self.__class__()

        for n in atoms_set:
            a = DynamicQueryElement.from_atom(self._atoms[n])
            sub.add_atom(a, n, xy=self._plane.get(n),
                         p_charge=self._p_charges.get(n, 0), p_is_radical=self._p_radicals.get(n, False),
                         neighbors=self._neighbors.get(n, ()), hybridization=self._hybridizations.get(n, ()),
                         p_neighbors=self._p_neighbors.get(n, ()), p_hybridization=self._p_hybridizations.get(n, ()))

        for n in atoms_set:
            for m, bond in self._bonds[n].items():
                if m in atoms_set and n < m:
                    sub.add_bond(n, m, bond.copy())
        return sub

    def union(self, other, *, remap: bool = False, copy: bool = True) -> 'QueryCGRContainer':
        if not isinstance(other, QueryCGRContainer):
            raise TypeError('QueryCGRContainer expected')

        u = super().union(other, remap=remap, copy=copy)
        if copy:
            # when copy, super().union already returned new instance (u)
            u._charges = {**self._charges, **other._charges}
            u._radicals = {**self._radicals, **other._radicals}
            u._p_charges = {**self._p_charges, **other._p_charges}
            u._p_radicals = {**self._p_radicals, **other._p_radicals}
            u._neighbors = {**self._neighbors, **other._neighbors}
            u._hybridizations = {**self._hybridizations, **other._hybridizations}
            u._p_neighbors = {**self._p_neighbors, **other._p_neighbors}
            u._p_hybridizations = {**self._p_hybridizations, **other._p_hybridizations}
            u._plane = {**self._plane, **other._plane}
        else:
            self._charges.update(other._charges)
            self._radicals.update(other._radicals)
            self._p_charges.update(other._p_charges)
            self._p_radicals.update(other._p_radicals)
            self._neighbors.update(other._neighbors)
            self._hybridizations.update(other._hybridizations)
            self._p_neighbors.update(other._p_neighbors)
            self._p_hybridizations.update(other._p_hybridizations)
            self._plane.update(other._plane)
            u = self
        return u

    def get_mapping(self, other: Union['QueryCGRContainer', 'cgr.CGRContainer'], **kwargs):
        if isinstance(other, (QueryCGRContainer, cgr.CGRContainer)):
            return super().get_mapping(other, **kwargs)
        raise TypeError('CGRContainer or QueryCGRContainer expected')

    def get_mcs_mapping(self, other: Union['QueryCGRContainer', 'cgr.CGRContainer'], **kwargs):
        if isinstance(other, (QueryCGRContainer, cgr.CGRContainer)):
            return super().get_mcs_mapping(other, **kwargs)
        raise TypeError('CGRContainer or QueryCGRContainer expected')

    @staticmethod
    def _validate_neighbors(neighbors):
        if neighbors is None:
            neighbors = ()
        elif isinstance(neighbors, int):
            if neighbors < 0 or neighbors > 14:
                raise ValueError('neighbors should be in range [0, 14]')
            neighbors = (neighbors,)
        elif isinstance(neighbors, (tuple, list)):
            if not all(isinstance(n, int) for n in neighbors):
                raise TypeError('neighbors should be list or tuple of ints')
            if any(n < 0 or n > 14 for n in neighbors):
                raise ValueError('neighbors should be in range [0, 14]')
            if len(set(neighbors)) != len(neighbors):
                raise ValueError('neighbors should be unique')
            neighbors = tuple(sorted(neighbors))
        else:
            raise TypeError('neighbors should be int or list or tuple of ints')
        return neighbors

    @staticmethod
    def _validate_hybridization(hybridization):
        if hybridization is None:
            hybridization = ()
        elif isinstance(hybridization, int):
            if hybridization < 1 or hybridization > 4:
                raise ValueError('hybridization should be in range [1, 4]')
            hybridization = (hybridization,)
        elif isinstance(hybridization, (tuple, list)):
            if not all(isinstance(h, int) for h in hybridization):
                raise TypeError('hybridizations should be list or tuple of ints')
            if any(h < 1 or h > 4 for h in hybridization):
                raise ValueError('hybridizations should be in range [1, 4]')
            if len(set(hybridization)) != len(hybridization):
                raise ValueError('hybridizations should be unique')
            hybridization = tuple(sorted(hybridization))
        else:
            raise TypeError('hybridization should be int or list or tuple of ints')
        return hybridization

    @staticmethod
    def _validate_neighbors_pairing(neighbors, p_neighbors):
        if len(neighbors) != len(p_neighbors):
            raise ValueError('neighbors and p_neighbors should be same length')
        if neighbors:
            if len(set(zip(neighbors, p_neighbors))) != len(neighbors):
                raise ValueError('paired neighbors and p_neighbors should be unique')
            neighbors, p_neighbors = zip(*sorted(zip(neighbors, p_neighbors)))
        return neighbors, p_neighbors

    @staticmethod
    def _validate_hybridization_pairing(hybridization, p_hybridization):
        if len(hybridization) != len(p_hybridization):
            raise ValueError('hybridization and p_hybridization should be same length')
        if hybridization:
            if len(set(zip(hybridization, p_hybridization))) != len(hybridization):
                raise ValueError('paired hybridization and p_hybridization should be unique')
            hybridization, p_hybridization = zip(*sorted(zip(hybridization, p_hybridization)))
        return hybridization, p_hybridization

    def __getstate__(self):
        atomic_numbers = {n: getattr(a, '_atomic_number', getattr(a, 'atomic_number', 0)) for n, a in self._atoms.items()}
        isotopes = {n: getattr(a, '_isotope', getattr(a, 'isotope', None)) for n, a in self._atoms.items()}
        return {'_atoms': self._atoms, '_bonds': self._bonds, '_charges': self._charges, '_radicals': self._radicals,
                'p_charges': self._p_charges, 'p_radicals': self._p_radicals, 'neighbors': self._neighbors,
                'hybridizations': self._hybridizations, 'p_neighbors': self._p_neighbors,
                'p_hybridizations': self._p_hybridizations, 'plane': self._plane,
                'atomic_numbers': atomic_numbers, 'isotopes': isotopes, **super().__getstate__()}

    def __setstate__(self, state):
        super().__setstate__(state)
        self._atoms = state['_atoms']
        self._bonds = state['_bonds']
        self._charges = state.get('_charges', {})
        self._radicals = state.get('_radicals', {})
        self._neighbors = state['neighbors']
        self._hybridizations = state['hybridizations']
        self._p_neighbors = state['p_neighbors']
        self._p_hybridizations = state['p_hybridizations']
        self._p_charges = state['p_charges']
        self._p_radicals = state['p_radicals']
        self._plane = state.get('plane', {})
        atomic_numbers = state.get('atomic_numbers', {})
        isotopes = state.get('isotopes', {})
        for n, atom in self._atoms.items():
            try:
                atom._atomic_number = atomic_numbers.get(n, getattr(atom, '_atomic_number', 0))
                atom._isotope = isotopes.get(n, getattr(atom, '_isotope', None))
                atom._attach_to_graph(self, n)
            except AttributeError:
                pass
            if n in self._plane and hasattr(atom, 'xy'):
                atom.xy = self._plane[n]

    @property
    def atoms_order(self) -> dict[int, int]:
        """
        Lightweight order mapping for SMILES generation when Morgan is absent.
        """
        return {n: i for i, n in enumerate(self._atoms, 1)}


__all__ = ['QueryCGRContainer']
