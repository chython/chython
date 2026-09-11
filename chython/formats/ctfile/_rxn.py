# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of chython.
#
#  chython is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""MDL reaction records: `$RXN` V2000, `$RXN V3000`, and the record object an RDfile yields.

Framing only, no column offsets: a molfile component goes to `parse_record`, a bare CTAB to
`parse_ctab_block`.  The two versions frame differently -- V2000 delimits components with `$MOL` and
gives each a three-line header, V3000 nests bare CTABs in `BEGIN REACTANT`/`PRODUCT`/`AGENT` with no
per-component header -- so writing a V3000 reaction as V2000 invents three blank header lines.
"""

from ._errors import MalformedCtfile
from ._sdf import V2000_STAMP, V3000_STAMP, parse_record
from ._sgroup import merge_log
from ._tokens import join_continuations, tokenize
from ._v2000 import emit_v2000
from ._v3000 import emit_v3000, parse_ctab_block
from ...core import LogRecord, LOST, REPAIRED, ReactionContainer
from ...core._reaction_passes import _located


RXN_HEADER_LINES = 4          # `$RXN` line plus name, program, comment -- counts is line index 4
#: `emit_v3000` writes a molfile: title, program, comment, blank, then the CTAB.  A reaction's
#: component has no header of its own, so the header is sliced off.  Not `RXN_HEADER_LINES`: both
#: being 4 is a coincidence of two unrelated formats.
_MOLFILE_HEADER_LINES = 4
_MOL_TAG = '$MOL'


def sniff_rxn_version(lines, log=None):
    """`V3000_STAMP` if the `$RXN` line names V3000 or the body holds `M  V30`, else `V2000_STAMP`.

    The body outranks the tag: a `$RXN` with no V3000 token whose body is all `M  V30` is a file some
    tool really produced, and reading it as V2000 finds no `$MOL`.

    The body scan stops at the first `$MOL`, because framing version and component-CTAB version are
    different questions: a V2000-framed reaction may hold V3000 component molfiles, and a `M  V30`
    inside one of them would otherwise mis-identify the framing.  V3000 framing has no `$MOL`.
    """
    log = [] if log is None else log
    # Tolerant by design: two spaces, a trailing field or a different case all sniff correctly.
    tagged = 'V3000' in lines[0].upper() if lines else False

    body = ''
    for line in lines[RXN_HEADER_LINES:]:
        if line.startswith(_MOL_TAG):
            # A $MOL line is proof of V2000 framing; stop here.
            break
        if line.startswith('M  V30'):
            body = V3000_STAMP
            break

    if tagged and not body:
        log.append(LogRecord('rxn:v3000-tag-no-body', (),
                             '$RXN V3000 tag present but no M  V30 line in the framing body; read as V3000',
                             REPAIRED))
    elif body and not tagged:
        log.append(LogRecord('rxn:v3000-body-no-tag', (),
                             '$RXN without the V3000 token but with M  V30 lines; read as V3000',
                             REPAIRED))
    return V3000_STAMP if (tagged or body) else V2000_STAMP


def _counts_v2000(line, log):
    """`(reactants, products, agents)` from the V2000 counts line.

    Three-character fields.  The third is not in the specification -- the counts line officially
    carries reactants and products -- but writing agents there is what the ecosystem does, so it is
    read and its presence is logged.
    """
    try:
        reactants = int(line[0:3] or 0)
        products = int(line[3:6] or 0)
    except ValueError:
        raise MalformedCtfile(f'unreadable $RXN counts line {line!r:.40}')
    third = line[6:9].strip()
    if not third:
        return reactants, products, 0
    try:
        agents = int(third)
    except ValueError:
        log.append(LogRecord('rxn:counts-third-unreadable', (),
                             f'$RXN counts line third field {third!r} unreadable, no agents read', LOST))
        return reactants, products, 0
    log.append(LogRecord('rxn:counts-third-as-agents', (),
                         f'$RXN counts line has a third count ({agents}); read as agents, which is a '
                         f'widespread convention and not in the specification', REPAIRED))
    return reactants, products, agents


def split_rxn_v2000(lines, log=None):
    """`((reactants, products, agents), [[molfile lines], ...])` for a V2000 reaction record."""
    log = [] if log is None else log
    if len(lines) <= RXN_HEADER_LINES:
        raise MalformedCtfile(f'$RXN record has {len(lines)} lines, too few for a counts line')
    counts = _counts_v2000(lines[RXN_HEADER_LINES], log)

    blocks = []
    current = None
    for line in lines[RXN_HEADER_LINES + 1:]:
        if line.startswith(_MOL_TAG):
            current = []
            blocks.append(current)
        elif current is not None:
            current.append(line)
    return counts, blocks


_V3000_ROLES = {'REACTANT': 0, 'PRODUCT': 1, 'AGENT': 2}
_ROLE_NAMES = ('reactant', 'product', 'agent')


def _read_v3000_sides(lines, log, spans, ignore_stereo):
    """`(counts, (reactants, products, agents))` for a `$RXN V3000` record.

    The body is joined and prefix-stripped once, then walked for role blocks.  Each `BEGIN CTAB` is
    handed to `parse_ctab_block` as a slice starting at that line; it stops at the matching `END CTAB`
    itself.  Component logs are prefixed by role and ordinal (`reactant 1: `), since V3000 states the
    role on the block, so two components logging the same string stay distinguishable.

    `spans` collects `(start, end)` slices of `log` that are one component's prefixed lines; see
    :func:`_mirror_components` for what reads them.
    """
    body = join_continuations(lines[RXN_HEADER_LINES:], log)
    counts = (0, 0, 0)
    sides = ([], [], [])
    role_counts = [0, 0, 0]   # per-role component counter for prefix numbering
    role = None
    i = 0
    while i < len(body):
        line = body[i].strip()
        upper = line.upper()
        if upper.startswith('COUNTS'):
            fields = tokenize(line)[1:]
            values = []
            for field in fields[:3]:
                try:
                    values.append(int(field))
                except ValueError:
                    log.append(LogRecord('rxn:counts-field-unreadable', (),
                                         f'M  V30 COUNTS field {field!r} unreadable, read as 0', REPAIRED))
                    values.append(0)
            counts = tuple(values + [0] * (3 - len(values)))
            i += 1
        elif upper.startswith('BEGIN ') and upper[6:].strip() in _V3000_ROLES:
            new_role = _V3000_ROLES[upper[6:].strip()]
            if role is not None:
                # A second role block opened before the first was closed.  V3000 requires explicit
                # END REACTANT/PRODUCT/AGENT framing; a file that omits the close is broken.
                log.append(LogRecord('rxn:role-framing-broken', (),
                                     f'M  V30 BEGIN {upper[6:].strip()} while {_ROLE_NAMES[role]} role is '
                                     f'already open; role framing is broken'))
            role = new_role
            i += 1
        elif upper.startswith('END ') and upper[4:].strip() in _V3000_ROLES:
            role = None
            i += 1
        elif upper.startswith('BEGIN CTAB'):
            if role is None:
                # A CTAB outside any role block.  V3000 always states the role, so this is a broken
                # file; reading it as a reactant is a guess, and saying so beats dropping it.
                log.append(LogRecord('rxn:ctab-outside-role', (),
                                     'M  V30 BEGIN CTAB outside a REACTANT/PRODUCT/AGENT block; read as a '
                                     'reactant', REPAIRED))
                role = 0
            role_counts[role] += 1
            prefix = f'{_ROLE_NAMES[role]} {role_counts[role]}: '
            component_log = []
            ctab = parse_ctab_block(body[i:], component_log)
            ctab.log = component_log
            molecule, _, build_log = ctab.build(ignore_stereo=ignore_stereo)
            # build_log is list(component_log) + build-time entries; merge all with role prefix.
            start = len(log)
            merge_log(log, build_log, prefix)
            spans.append((start, len(log)))
            sides[role].append(molecule)
            # Skip past this CTAB's END so the outer walk does not re-enter it.
            depth = 0
            while i < len(body):
                probe = body[i].strip().upper()
                if probe.startswith('BEGIN CTAB'):
                    depth += 1
                elif probe.startswith('END CTAB'):
                    depth -= 1
                    if depth <= 0:
                        i += 1
                        break
                i += 1
            else:
                # END CTAB was never found.  parse_ctab_block already logged the cause; here we
                # report the consequence: everything after this CTAB in the record is unreachable.
                log.append(LogRecord('rxn:end-ctab-missing', (),
                                     f'{prefix}END CTAB missing; components after this point in the record '
                                     f'are lost', LOST))
        else:
            i += 1
    return counts, tuple(sides)


def parse_rxn(lines, log=None, *, ignore_stereo=False):
    """One reaction record's lines -- `$RXN` or `$RXN V3000` -- into a `ReactionContainer`.

    An empty reaction is not an error: an RDfile with a placeholder `$RFMT` is a real file, and the
    container permits one so no record is lost to say that it declared nothing.

    Everything read lands on ``reaction.log`` as well as on `log`: the record's own framing lines
    directly, each component's under the ``subject`` that names it.  See :func:`_mirror_components`.
    """
    log = [] if log is None else log
    # THIS RECORD'S OWN LIST.  Written in file order, handed to the caller whole at the end, and split
    # by `spans` into what is about the reaction record and what is about one of its components.
    own = []
    spans = []
    version = sniff_rxn_version(lines, own)
    if version == V3000_STAMP:
        counts, sides = _read_v3000_sides(lines, own, spans, ignore_stereo)
        found = tuple(len(x) for x in sides)
        if counts != found:
            # V3000 names each component's role, so the role blocks are evidence and the COUNTS
            # line is a claim.  Log the disagreement and keep the roles -- unlike V2000, where
            # the counts line is all there is to split on.
            own.append(LogRecord('rxn:counts-role-mismatch', (),
                                 f'M  V30 COUNTS states {counts} and the role blocks hold {found}; the role '
                                 f'blocks decide the sides'))
    else:
        counts, blocks = split_rxn_v2000(lines, own)
        molecules = []
        for i, block in enumerate(blocks):
            component_log = []
            try:
                molecule = parse_record(block, component_log, ignore_stereo=ignore_stereo)
            except MalformedCtfile as e:
                # One unreadable component does not cost the record its other three.
                own.append(LogRecord('rxn:component-read-failed', (),
                                     f'component {i + 1} could not be read and is dropped: {e}', LOST))
                continue
            # Prefix by file position, not by role: in V2000 the sides are not known until
            # `_apportion`.  Two components logging the same thing stay distinguishable.
            start = len(own)
            merge_log(own, component_log, f'component {i + 1}: ')
            spans.append((start, len(own)))
            molecules.append(molecule)
        sides = _apportion(molecules, counts, own)

    reaction = ReactionContainer(*sides, title=_title_of(lines))
    log.extend(own)
    _mirror_components(reaction, own, spans)
    return reaction


def _mirror_components(reaction, own, spans):
    """`reaction.log` gets the record's own lines; each component's are mirrored under its `subject`.

    The mirrored copy is read off the component's own log rather than out of `own`, and the slices
    `spans` names are left out of the direct absorb, so no component event reaches `reaction.log`
    twice.  The two spellings differ on purpose: the caller's flat list carries the role prefix, which
    is the only way to tell two components apart in one sequence, while a record on `reaction.log`
    says which molecule it is about in `subject` -- `LogRecord.atoms` are stable ids in ONE container,
    so a reaction log pooling three sides' records unstamped hands back numbers naming a different
    atom depending on which component you read them against.  Same arrangement
    `core/_reaction_passes.py` gives a pass, and the reason `ReactionContainer.log` documents.
    """
    covered = {i for start, end in spans for i in range(start, end)}
    reaction.log.absorb('read', [r for i, r in enumerate(own) if i not in covered])
    for where, molecule in _located(reaction):
        with reaction.log.stage('read', subject=where) as out:
            out.extend(molecule.log)


def _apportion(molecules, counts, log):
    """Split a flat component list into `(reactants, products, agents)` by the counts line.

    File order is reactants, then products, then agents.  When the counts disagree with what was
    found the counts line decides the split and the shortfall is logged: a truncated file is commoner
    than a mis-counted one, and re-deriving the split turns a missing reactant into a product.
    """
    reactants, products, agents = counts
    promised = reactants + products + agents
    if promised != len(molecules):
        log.append(LogRecord('rxn:counts-molecule-mismatch', (),
                             f'$RXN counts line promises {promised} component(s) and {len(molecules)} were '
                             f'read; the counts line decides the sides', REPAIRED))
    return (molecules[:reactants],
            molecules[reactants:reactants + products],
            molecules[reactants + products:reactants + products + agents])


def _title_of(lines):
    """The reaction name line.

    The name line is stored as the ``str`` the file decoded to, which is what
    :attr:`ReactionContainer.title` is.
    """
    return lines[1].rstrip() if len(lines) > 1 else ''


def parse_rxn_record(body, fields, log=None, *, ignore_stereo=False, header=None):
    """One RDfile ``$RFMT`` record as a `ReactionContainer`, from the body lines and its data fields.

    ``body`` is the part of the record before the first ``$DTYPE`` -- the ``$RXN`` header and the
    component molfiles; ``fields`` is the ``{name: value}`` tail from ``parse_rdf_fields``, which lands
    on ``reaction.meta``.  The name line is ``reaction.title``.

    `header`, when given a dict, takes ``version``, ``program`` and ``comment`` -- the two header lines
    no container holds, and the version actually read.  It is sniffed into a throwaway log, because
    ``parse_rxn`` sniffs the same body into the real one and would double the message.
    """
    log = [] if log is None else log
    version = sniff_rxn_version(body, [])   # throwaway: parse_rxn sniffs authoritatively below
    reaction = parse_rxn(body, log, ignore_stereo=ignore_stereo)
    if fields:
        reaction.meta.update(fields)
    if header is not None:
        header['version'] = version
        header['program'] = body[2].rstrip() if len(body) > 2 else ''
        header['comment'] = body[3].rstrip() if len(body) > 3 else ''
    return reaction


def emit_rxn(reaction, *, version=V2000_STAMP, title=None, program='', comment='', log=None):
    """Render one reaction record.  `(lines, log)`, no trailing newlines, no `$RFMT`.

    The RDfile's `$RFMT` line is the RDfile's business and is added by `_rdf.py`, so this function
    writes something that is a valid standalone `.rxn` file.
    """
    log = [] if log is None else log
    if version not in (V2000_STAMP, V3000_STAMP):
        raise MalformedCtfile(f'unknown reaction file version {version!r}; expected {V2000_STAMP} '
                              f'or {V3000_STAMP}')
    name = reaction.title if title is None else title
    reactants, products, agents = reaction.reactants, reaction.products, reaction.agents

    if version == V2000_STAMP:
        n_r, n_p, n_a = len(reactants), len(products), len(agents)
        for label, count in (('reactants', n_r), ('products', n_p), ('agents', n_a)):
            if count > 999:
                raise MalformedCtfile(f'{count} {label} will not fit the V2000 3-character count '
                                      f'field; write this reaction as V3000')
        if agents:
            # Not an escalation: a V2000 RXN carrying a third count is what the ecosystem reads, and
            # the unofficial field is logged rather than switching the caller's chosen version.
            log.append(LogRecord('rxn:agents-unofficial-count', (),
                                 f'{n_a} agent(s) written in the counts line\'s third count, which is '
                                 f'a widespread convention and not in the specification'))
            counts = f'{n_r:3d}{n_p:3d}{n_a:3d}'
        else:
            counts = f'{n_r:3d}{n_p:3d}'
        lines = ['$RXN', name, program, comment, counts]
        for molecule in (*reactants, *products, *agents):
            block, log = emit_v2000(molecule, log=log)
            lines.append('$MOL')
            lines.extend(block)
        return lines, log

    counts = f'M  V30 COUNTS {len(reactants)} {len(products)}'
    if agents:
        counts += f' {len(agents)}'
    lines = ['$RXN V3000', name, program, comment, counts]
    for role, side in (('REACTANT', reactants), ('PRODUCT', products), ('AGENT', agents)):
        if not side:
            continue
        lines.append(f'M  V30 BEGIN {role}')
        for molecule in side:
            block, log = emit_v3000(molecule, log=log)
            # `emit_v3000` writes a molfile: four header lines, the CTAB, then `M  END`.  A reaction's
            # component has no header and the record carries one `M  END` at the end, so both are
            # sliced off rather than adding a second V3000 CTAB writer.
            body = block[_MOLFILE_HEADER_LINES:]
            while body and body[-1].startswith('M  END'):
                body.pop()
            lines.extend(body)
        lines.append(f'M  V30 END {role}')
    lines.append('M  END')
    return lines, log
