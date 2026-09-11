"""Scratch measurement harness for task 10 fix round 1 (not part of the library)."""
from chython.core import MoleculeContainer, QueryContainer


def target(parity, fourth=None, group=None, explicit_h=False):
    m = MoleculeContainer()
    with m.edit():
        els = ['C', 'F', 'Cl', 'Br']
        if fourth is not None:
            els.append(fourth)
        if explicit_h:
            els.append('H')
        sids = [m.add_atom(e, implicit_h=1 if (k == 0 and fourth is None and not explicit_h) else 0)
                for k, e in enumerate(els)]
        for s in sids[1:]:
            m.add_bond(sids[0], s, 1)
        if parity:
            m.set_parity(sids[0], parity)
        if group is not None:
            m.set_stereo_group(sids[0], group[0], group[1])
    return m, sids


def q_tokens(numbers, tokens, orders=None):
    """numbers: element per query atom (atom 0 is the stereo atom). tokens: extra ops on atom 0."""
    q = QueryContainer()
    qids = [q.add_atom() for _ in numbers]
    for sid, number in zip(qids, numbers):
        q.atom_primitive(sid, 'element', number)
    if orders is None:
        orders = [1] * (len(numbers) - 1)
    for s, order in zip(qids[1:], orders):
        q.add_bond(qids[0], s)
        q.bond_primitive(qids[0], s, 'bond_order', order)
    for kind, arg in tokens:
        if kind == 'op':
            q.atom_operator(qids[0], arg)
        else:
            q.atom_primitive(qids[0], kind, arg)
    return q


def show(label, fn):
    try:
        print('%-58s %s' % (label, fn()))
    except Exception as exc:
        print('%-58s %s: %s' % (label, type(exc).__name__, exc))


print('--- F87 case 1: [C@,N] against N(F)(Cl)Br ---')
n_target = MoleculeContainer()
with n_target.edit():
    s = [n_target.add_atom('N'), n_target.add_atom('F'), n_target.add_atom('Cl'),
         n_target.add_atom('Br')]
    for x in s[1:]:
        n_target.add_bond(s[0], x, 1)
# [C@,N]: element 6 AND_HIGH stereo 1 OR element 7
cn = q_tokens((6, 9, 17, 35), [('op', 'and_high'), ('stereo', 1), ('op', 'or'), ('element', 7)])
plain = q_tokens((6, 9, 17, 35), [('op', 'or'), ('element', 7)])
show('[C@,N] vs N(F)(Cl)Br  (want True)', lambda: cn.is_substructure(n_target))
show('[C,N]  vs N(F)(Cl)Br  (control True)', lambda: plain.is_substructure(n_target))
show('[C@,N] box_counts', lambda: cn.box_counts())

print('--- F87 case 2: [C;@,D3] against a parity-less C(F)(Cl)Br ---')
bare, _ = target(0)
cd = q_tokens((6, 9, 17, 35), [('op', 'and_low'), ('stereo', 1), ('op', 'or'), ('degree', 3)])
cc = q_tokens((6, 9, 17, 35), [('op', 'and_low'), ('charge', 0), ('op', 'or'), ('degree', 3)])
show('[C;@,D3]  vs C(F)(Cl)Br no parity (want True)', lambda: cd.is_substructure(bare))
show('[C;+0,D3] vs same (control True)', lambda: cc.is_substructure(bare))
show('[C;@,D3] box_counts (merge must be blocked)', lambda: cd.box_counts())

print('--- F87 case 3: [C;@,N;@@] against a parity-2 carbon ---')
p2, _ = target(2)
show('[C;@,N;@@] vs parity-2 C (want False/raise)',
     lambda: q_tokens((6, 9, 17, 35),
                      [('op', 'and_low'), ('stereo', 1), ('op', 'or'), ('element', 7),
                       ('op', 'and_low'), ('stereo', 2)]).is_substructure(p2))
show('[C;@;@@] vs parity-1 C (want False/raise)',
     lambda: q_tokens((6, 9, 17, 35),
                      [('op', 'and_low'), ('stereo', 1), ('op', 'and_low'),
                       ('stereo', 2)]).is_substructure(target(1)[0]))

print('--- F88: explicit vs implicit hydrogen ---')
eh, ehs = target(1, explicit_h=True)
print('explicit-H unit  ', eh.unit_of(ehs[0]))
ih, ihs = target(1)
print('implicit-H unit  ', ih.unit_of(ihs[0]))
for sign in (1, 2):
    show('explicit-H target vs 3-named query sign %d' % sign,
         lambda sign=sign: q_tokens((6, 9, 17, 35),
                                    [('op', 'and_low'), ('stereo', sign)]).is_substructure(eh))
    show('implicit-H target vs 3-named query sign %d' % sign,
         lambda sign=sign: q_tokens((6, 9, 17, 35),
                                    [('op', 'and_low'), ('stereo', sign)]).is_substructure(ih))
print('named heavy unaccounted (F67, want False both signs):')
for sign in (1, 2):
    for parity in (1, 2):
        show('  I-target parity %d vs sign %d' % (parity, sign),
             lambda sign=sign, parity=parity: q_tokens(
                 (6, 9, 17, 35), [('op', 'and_low'), ('stereo', sign)]
             ).is_substructure(target(parity, fourth='I')[0]))

print('--- verification 3: two explicit hydrogens ---')
two = MoleculeContainer()
with two.edit():
    s = [two.add_atom('C'), two.add_atom('F'), two.add_atom('Cl'), two.add_atom('H'),
         two.add_atom('H')]
    for x in s[1:]:
        two.add_bond(s[0], x, 1)
    two.set_parity(s[0], 1)
print('units          ', two.stereo_units())
print('unit_of(centre)', two.unit_of(s[0]))
show('CH2(F)Cl 2xH explicit vs 3-named query @',
     lambda: q_tokens((6, 9, 17, 1), [('op', 'and_low'), ('stereo', 1)]).is_substructure(two))

print('--- h / H primitives are independent of F88 ---')
for name, val in (('implicit_h', 1), ('total_h', 1)):
    for label, mol in (('implicit-H twin', ih), ('explicit-H twin', eh)):
        show('[C;%s%d;@] vs %s' % (name, val, label),
             lambda name=name, val=val, mol=mol: q_tokens(
                 (6, 9, 17, 35),
                 [('op', 'and_low'), (name, val), ('op', 'and_low'), ('stereo', 1)]
             ).is_substructure(mol))
