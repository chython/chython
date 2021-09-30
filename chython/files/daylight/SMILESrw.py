from itertools import chain
from functools import reduce
from operator import or_
from re import split, compile, fullmatch, findall, search
from ._parser import DaylightParser
from ...containers import ReactionContainer

cx_fragments = compile(r'f:(?:[0-9]+(?:\.[0-9]+)+)(?:,(?:[0-9]+(?:\.[0-9]+)+))*')
cx_radicals = compile(r'\^[1-7]:[0-9]+(?:,[0-9]+)*')
delimiter = compile(r'[=:]')


class SMILESRead(DaylightParser):
    def parse(self, smiles: str):
        """SMILES string parser."""
        if not smiles:
            raise ValueError('Empty string')

        smi, *data = smiles.split()
        if data and data[0].startswith('|') and data[0].endswith('|'):
            fr = search(cx_fragments, data[0])
            if fr is not None:
                contract = [sorted(int(x) for x in x.split('.')) for x in fr.group()[2:].split(',')]
                if len({x for x in contract for x in x}) < len([x for x in contract for x in x]):
                    self._info(f'collisions in cxsmiles fragments description: {data[0]}')
                    contract = None
                elif any(x[0] < 0 for x in contract):
                    self._info(f'invalid cxsmiles fragments description: {data[0]}')
                    contract = None

                radicals = [int(x) for x in findall(cx_radicals, data[0]) for x in x[3:].split(',')]
                if any(x < 0 for x in radicals):
                    self._info(f'invalid cxsmiles radicals description: {data[0]}')
                    radicals = []
                if len(set(radicals)) != len(radicals):
                    self._info(f'collisions in cxsmiles radicals description: {data[0]}')
                    radicals = []
                data = data[1:]
            else:
                radicals = [int(x) for x in findall(cx_radicals, data[0]) for x in x[3:].split(',')]
                if radicals:
                    if any(x < 0 for x in radicals):
                        self._info(f'invalid cxsmiles radicals description: {data[0]}')
                        radicals = []
                    if len(set(radicals)) != len(radicals):
                        self._info(f'collisions in cxsmiles radicals description: {data[0]}')
                        radicals = []
                    data = data[1:]
                contract = None
        else:
            radicals = []
            contract = None

        if self.__header is None:
            meta = {}
            for x in data:
                try:
                    k, v = split(delimiter, x, 1)
                    meta[k] = v
                except ValueError:
                    self._info(f'invalid metadata entry: {x}')
        else:
            meta = dict(zip(self.__header, data))

        if '>' in smi and (smi[smi.index('>') + 1] in '>([' or smi[smi.index('>') + 1].isalpha()):
            record = {'reactants': [], 'reagents': [], 'products': [], 'meta': meta, 'title': ''}
            try:
                reactants, reagents, products = smi.split('>')
            except ValueError as e:
                raise ValueError('invalid reaction smiles') from e

            for k, d in zip(('reactants', 'products', 'reagents'), (reactants, products, reagents)):
                if d:
                    for x in d.split('.'):
                        if not x:
                            if self._ignore:
                                self._info('two dots in line ignored')
                            else:
                                raise ValueError('invalid reaction smiles. two dots in line')
                        else:
                            record[k].append(self.__parse_tokens(x))

            if radicals:
                atom_map = dict(enumerate(a for m in chain(record['reactants'], record['reagents'], record['products'])
                                          for a in m['atoms']))
                for x in radicals:
                    atom_map[x]['is_radical'] = True

            container = self._convert_reaction(record)
            if contract:
                if max(x for x in contract for x in x) >= len(container):
                    self._info(f'skipped invalid contract data: {contract}')
                lr = len(container.reactants)
                reactants = set(range(lr))
                reagents = set(range(lr, lr + len(container.reagents)))
                products = set(range(lr + len(container.reagents),
                                     lr + len(container.reagents) + len(container.products)))
                new_molecules = [None] * len(container)
                for c in contract:
                    if reactants.issuperset(c):
                        new_molecules[c[0]] = reduce(or_, (container.reactants[x] for x in c))
                        reactants.difference_update(c)
                    elif products.issuperset(c):
                        new_molecules[c[0]] = reduce(or_, (container.products[x - len(container)] for x in c))
                        products.difference_update(c)
                    elif reagents.issuperset(c):
                        new_molecules[c[0]] = reduce(or_, (container.reagents[x - lr] for x in c))
                        reagents.difference_update(c)
                    else:
                        self._info(f'impossible to contract different parts of reaction: {contract}')
                for x in reactants:
                    new_molecules[x] = container.reactants[x]
                for x in products:
                    new_molecules[x] = container.products[x - len(container)]
                for x in reagents:
                    new_molecules[x] = container.reagents[x - lr]

                meta = container.meta
                if self._store_log:
                    if log := self._format_log():
                        if 'ParserLog' in meta:
                            meta['ParserLog'] += '\n' + log
                        else:
                            meta['ParserLog'] = log
                else:
                    self._flush_log()
                return ReactionContainer([x for x in new_molecules[:lr] if x is not None],
                                         [x for x in new_molecules[-len(container.products):] if x is not None],
                                         [x for x in new_molecules[lr: -len(container.products)] if x is not None],
                                         meta=meta)
            return container
        else:
            record = self.__parse_tokens(smi)
            record['meta'].update(meta)
            if 'cgr' in record:  # CGR smiles parser
                return self._convert_cgr(record)
            for x in radicals:
                record['atoms'][x]['is_radical'] = True
            return self._convert_molecule(record)


__all__ = ['SMILESRead']
