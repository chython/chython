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
"""The reaction container on the chython 3 core: three sides, metadata, and the ML view.

WHY THIS IS PYTHON AND NOT A `.pxi`.  A reaction holds three tuples, a dict and a title; there is no
loop in it for C to make faster, and RULES.md's chapter 3 is about not relying on the optimiser
INSIDE the arena rather than an instruction to compile glue.  The dict assembly in `modeling_view`
unpacks the kernel's arrays -- the per-atom work is already in C.  If a measurement ever says
otherwise the class moves; until then this file is where the reaction's SEMANTICS are, in the
language a reader can follow, and `chython/core/RULES.md` chapter 10 says a claim like "this would
be faster in C" is not actionable until someone has measured it.

WHAT IS DELIBERATELY NOT HERE, WITH REASONS -- because the absences are the design.

`CGRContainer` and `compose()` / `~rxn`.  A CGR exists in this library for reaction machine learning
and nothing else, and the actual downstream consumer (a tokenizer) never builds one: it reads a
reaction and wants per-atom, per-side numbers.  So the ML path is served directly by
:meth:`modeling_view`, and the general-purpose overlay container -- dynamic bonds, dynamic atoms,
`decompose`, CGR SMILES, CGR isomorphism, CGR depiction, `centers_list` -- does not exist on this
release.  A caller who genuinely wants the union of two graphs has `MoleculeContainer.union`.

WHAT `__eq__` AND `__hash__` REST ON
------------------------------------

Not a canonical reaction SMILES: a canonical writer that oscillates on symmetric stereocentres makes a
container compare unequal to itself across a round trip.  Molecule `__eq__`/`__hash__` rest on
`mol_identity_bytes`, and that oscillation -- the mirror automorphism -- is closed in the canonical
search, so `MoleculeContainer` is sound to use as a dict key and a multiset comparison over the three
sides is a few lines of work.

Three SEMANTIC questions are decided deliberately, because this is public API and every wrong answer
here is silent.  The answers live on `__eq__` where a reader meets them:

  * do AGENTS participate?  YES.  A reaction is a record of what was done, so the same transformation
    run with and without a named catalyst is two reactions.  The alternative reading, in which a
    reaction is a transformation and agents are circumstance, is not lost: the sides are compared
    separately, so `rxn._identity()[0::2]` IS that comparison and a caller who wants it can spell it
    in one expression;
  * does atom-to-atom MAPPING participate?  NO -- two mappings of one reaction are the same chemistry
    and different data.  The mapping epic's own metric, "did we reproduce this mapping", therefore
    cannot be spelled `a == b` and compares map numbers explicitly.  That is the correct home for it
    -- it is a question about an annotation, not about a reaction;
  * does `title` participate?  NO.  It is not chemistry, and neither is `meta`.  Written down because
    it is the kind of field that gets folded in by accident, and an absence with no reason beside it
    is what invites the accident.

Reaction SMILES.  `write_smiles` writes ONE molecule including its CXSMILES tail, and a tail's
indices count atoms from the start of the whole string -- so concatenating three sides' output with
`>` produces a string with tails stranded in the middle of it, which every reader rejects.  Writing
a correct one means aggregating radical indices and fragment groups across the whole reaction, and
that is the writer's job, not the container's: `write_reaction_smiles` lives next to `smw_tail_text`
and `DetachedSmiles.tail`, where that machinery already is.  `__str__`, `__format__` and `.smiles`
here are three doors onto that one function and hold no formatting of their own, which is why a
nearly-right reaction SMILES cannot appear on this side of the boundary.  `__repr__` stays a
structural summary: it answers what the container holds, where `str` answers what the reaction is.
"""
import zlib
from collections.abc import Iterator, Sequence
from itertools import chain
from warnings import warn

from ._core import (H_UNKNOWN, MoleculeContainer, TensorEncoding, pach_dump, pach_load,
                    pach_record_length, reaction_transition_view, write_reaction_smiles,
                    _reaction_depict_fns, _reaction_draw_fns,
                    _reaction_attention_fn, _reaction_interop_fn, _reaction_reconstruct_fn,
                    _set_reaction_factory)
from ._log import Log
from ._reaction_passes import (reaction_canonicalize, reaction_clean_isotopes, reaction_clean_stereo,
                               reaction_contract_ions, reaction_explicify_hydrogens,
                               reaction_implicify_hydrogens, reaction_kekule, reaction_neutralize,
                               reaction_remove_reagents, reaction_reset_mapping, reaction_standardize,
                               reaction_thiele)


# --------------------------------------------------------------------------------------------------
# THE REACTION-LEVEL pach FORMAT, which is four bytes in front of a concatenation of molecule records:
#
#   byte 0     0x01, the version 1 format.  Version 5 is the same four bytes; see below.
#   byte 1     reactant count, uint8
#   byte 2     AGENT count, uint8       <- the middle field.  chython 2 spelled this side `reagents`.
#   byte 3     product count, uint8
#   the rest   the molecules' pach records, UNCOMPRESSED, concatenated in `molecules()` order:
#              reactants, then agents, then products.
#
# zlib on top of the whole thing by default, exactly as on the molecule side.
#
# THERE IS NO LENGTH FIELD, anywhere.  A molecule record's length is a function of its own header --
# for versions 0 and 2: atom count, neighbour counts and cis/trans count; for versions 3 and 4: atom
# count, bond count, stereo count, sgroup count and the map flag -- so a reader walks the stream with
# `pach_record_length` and a wrong count and a truncated buffer are the same failure seen twice.
#
# WHY THIS IS RESTORED RATHER THAN REPLACED.  Nothing new should be stored in it; it exists so that
# buffers chython 2 wrote can still be read, and there are corpora of them that cannot be regenerated.
# The bytes are pinned by `core/test/reaction_pach_v2_corpus.bin.gz`, written by an installed 2.24.
#
# VERSION 5 is the same four bytes over third-generation molecule records:
#
#   byte 0     0x05, and bytes 1-3 the three counts, exactly as version 1
#   the rest   version 3 or version 4 molecule records, uncompressed, in `molecules()` order
#
# The difference is entirely in the mapping: the molecule record has a map-number field, so nothing is
# relabelled on the way out and nothing is un-relabelled on the way in.  That removes version 1's one
# ambiguity -- a record it wrote cannot distinguish "mapped 1..N" from "unmapped and auto-numbered
# 1..N" -- and makes a PARTIALLY mapped molecule writable, which one number field per atom cannot
# express.  A version 5 reader must not reverse a relabelling that did not happen, which is what the
# version byte is for.
#
# ==================================================================================================
# THE ONE HARD PART: WHERE THE ATOM-TO-ATOM MAPPING LIVES.
# ==================================================================================================
#
# chython 2 had NO map-number field on an atom.  An atom's NUMBER was its mapping -- its SMILES reader
# put `[C:7]`'s 7 straight into the graph's key and numbered unmapped atoms consecutively across the
# whole record -- and the reaction packer wrote those numbers into pach's 12-bit atom-number field.
# The arena splits the two: `n` is a container's private label and `map_number` is chemistry,
# and versions 0 and 2 of `pach_dump` refuse a molecule carrying a map number because those molecule
# versions have no field for one.  Versions 3 and 4 have a map block.
#
# So the reaction writer puts the map numbers into that field, and the reader takes them out again.
# VERSION 1 ONLY:
#
#   WRITING   a molecule whose every atom carries a map number is relabelled onto those numbers and
#             written with `drop=['map_number']` -- dropped from the molecule layer because the
#             reaction layer is carrying it.  A molecule with NO mapping is written on its stable ids,
#             which is what chython 2 would have written for it.  A PARTIALLY mapped molecule has no
#             honest spelling in a format with one number field per atom and is refused by name.
#
#   READING   every atom's map number is set to the number the record gave it.  UNCONDITIONALLY, and
#             this is the one place the format is genuinely ambiguous: a chython 2 record cannot
#             distinguish "mapped 1..N" from "unmapped and auto-numbered 1..N", because chython 2
#             could not distinguish them either.  Reporting the number is the choice that cannot
#             silently corrupt a real corpus -- dropping it would lose the mapping of every mapped
#             reaction ever stored, which is the whole reason reactions get packed, while setting it
#             on an unmapped record hands back a number the record actually contained.  A caller who
#             knows its corpus is unmapped can clear it; a caller whose mapping was dropped has no way
#             to get it back.  Pinned by `test_reaction_pach.py`.
#
#   A version 5 record's molecule carries the map block itself, so re-deriving a mapping from atom
#   numbers would overwrite a real one with a renumbering.
# ==================================================================================================

_PACH_REACTION_VERSION_1 = 1
_PACH_REACTION_VERSION_5 = 5
_PACH_REACTION_VERSIONS = frozenset({_PACH_REACTION_VERSION_1, _PACH_REACTION_VERSION_5})
# uint8 per side.  chython 2 built the header with `bytearray((1, len(reactants), ...))`, which raises
# `ValueError: byte must be in range(0, 256)` -- so it never silently truncated, and NO STORED BUFFER
# CAN HOLD such a record.  V3 refuses too, saying which side and how many; widening the field would
# invent a version byte for data that does not exist.
_PACH_SIDE_MAX = 255
# The 12-bit atom-number field. Named here because the writer checks map numbers against it and owes a
# message about map numbers rather than about stable ids.
_PACH_MAX_NUMBER = 4095

# The molecule writer's waivers. `title` and `meta` each cover both layers -- the reaction's name line
# and metadata and each component's -- because pach has no field for either on either layer, so one
# `drop=['meta']` waives the reaction's own dict and its molecules' through `molecule_drop` below. The
# `| {'meta'}` is therefore a no-op in value and kept as the statement of intent.
_PACH_MOLECULE_DROP_NAMES = frozenset({'map_number', 'title', 'sgroups', 'cip', 'wedges',
                                       'stereo_groups', 'stereo', 'coordinates', 'conformers',
                                       'meta'})
_PACH_DROP_NAMES = _PACH_MOLECULE_DROP_NAMES | {'meta'}


def _as_title(value) -> str:
    """Coerce a title to `str`, and refuse the rest.

    `bytes` is DECODED with `surrogateescape` rather than refused, because a caller holding a raw name
    line is honoured: every surrogate the handler makes re-encodes to the byte it came from, so honouring
    it costs no second spelling of the field.  The molecule's `_as_title_bytes` is the same rule read the
    other way.
    """
    if isinstance(value, str):
        return value
    elif isinstance(value, bytes):
        return value.decode('utf8', 'surrogateescape')
    elif isinstance(value, (bytearray, memoryview)):
        return bytes(value).decode('utf8', 'surrogateescape')
    raise TypeError(f'title must be str, bytes, bytearray or memoryview, '
                    f'got {type(value).__name__}')


class ReactionModelingView:
    """Per-atom, per-side numbers for a mapped reaction -- the whole of the ML surface.

    This is what replaces the CGR for modelling.  It is keyed by ATOM MAP NUMBER, because a map
    number is the only thing that identifies "the same atom" across the arrow -- keying by atom number
    works only where a reader makes the two coincide.

    IT REPORTS WHAT THE RECORD DID NOT STATE instead of leaving it undefined.  Calling unmapped atoms
    and colliding map numbers "undefined results" means, for a training set, silent corruption at
    whatever rate the corpus happens to contain.  Here they land in :attr:`unmapped` and
    :attr:`collisions`, both empty for a well-mapped reaction, so a pipeline can decide -- and a
    refusal, if one is wanted, happens at the pipeline's boundary rather than inside a container that
    was handed a record it did not choose.
    """
    __slots__ = ('_states', '_union_bonds', '_unmapped', '_collisions')

    def __init__(self, states, union_bonds, unmapped, collisions):
        self._states = states
        self._union_bonds = union_bonds
        self._unmapped = unmapped
        self._collisions = collisions

    @property
    def states(self) -> dict[int, tuple]:
        """`{map_number: (element, h_before, n_before, h_after, n_after)}`, in union order.

        A NEGATIVE KEY IS AN ATOM THE RECORD DOES NOT NUMBER, allocated `-1` onwards in that order.  It
        is placed like any other atom, as reactant-only or product-only; the negative key says the
        pairing is this view's reading of an absence rather than something the record stated.

        `element` is the atomic number.  `h_*` are IMPLICIT hydrogen counts per side and `n_*` are
        heavy-atom neighbour counts per side -- the per-side hydrogen count being the thing the
        tokenizer actually consumes, and the reason this view exists at all.

        H_UNKNOWN (15) WHERE THE RECORD DOES NOT STATE A COUNT, not 0.  An `h or 0` here trains a
        molecule with an unstated hydrogen count as though it had none -- a plausible wrong number,
        which is the worst kind.  15 is the core's own sentinel and is
        outside the 0..14 domain, so a consumer that ignores it gets an obviously broken value and
        one that checks for it gets the truth.

        A REACTANT-ONLY ATOM keeps its hydrogen count on the after side and counts only its bonds to
        other reactant-only atoms: the fragment leaves intact, the bonds joining it to the retained
        part do not survive.  A product-only atom is the mirror image.  That is a modelling
        convention rather than chemistry, it is the one the reference consumer was written against,
        and it is stated here so that changing it is a decision.
        """
        return self._states

    @property
    def union_bonds(self) -> dict[tuple[int, int], tuple[int, int]]:
        """`{(n, m): (order_before, order_after)}` with `n < m`, `0` meaning the bond is absent.

        A bond appears when it exists on at least one side, so `(1, 0)` is a broken bond, `(0, 1)` a
        formed one and `(1, 2)` an order change.  This is the whole of what a CGR's dynamic bond
        carried, without a container to hold it.
        """
        return self._union_bonds

    @property
    def unmapped(self) -> dict[str, int]:
        """How many atoms carried no map number, per side: `{'reactants': n, 'products': n}`.

        Not an error and not a refusal.  Each is placed on the side it came from and keyed negatively in
        :attr:`states`, so no neighbour of one loses a degree; what the count is for is the atom the
        record leaves bare on BOTH sides, which becomes two rows because nothing pairs them.  A consumer
        that requires a fully mapped reaction checks this rather than discovering the hole in its loss
        curve.
        """
        return self._unmapped

    @property
    def collisions(self) -> dict[str, tuple[int, ...]]:
        """Map numbers used by more than one atom on the same side, per side.

        A collision means the union silently merged two atoms.  Reported for the same reason as
        :attr:`unmapped`: it is the caller's record and the caller's decision.
        """
        return self._collisions

    def __repr__(self):
        return (f'{type(self).__name__}({len(self._states)} atoms, {len(self._union_bonds)} bonds, '
                f'unmapped={self._unmapped}, collisions={self._collisions})')


class MappingResult:
    """What :meth:`ReactionContainer.attention_mapping` did, beyond the map numbers it wrote.

    FOUR FIELDS AND NOT A BOOL.  "Ran and changed nothing", "refused to run" and "ran and could not
    place two atoms" are three different answers, and a single `bool` -- or a `bool | float` chosen by
    a keyword -- reports the first for all three.  A caller filtering a corpus needs to tell them
    apart, and the score is what it filters on.

    Lives here rather than in `chython.reactions` because the core owns the method name, so it owns the
    type that method's signature names; the package one layer up owns the numbers in it.
    """
    __slots__ = ('_changed', '_score', '_unplaced', '_skipped')

    def __init__(self, changed, score, unplaced, skipped):
        self._changed = changed
        self._score = score
        self._unplaced = unplaced
        self._skipped = skipped

    @property
    def changed(self) -> bool:
        """The mapping written differs from the one the record carried."""
        return self._changed

    @property
    def score(self) -> float:
        """Mean attention over the placements that were accepted; `0.0` when none were.

        A MODEL CONFIDENCE AND NOT AN ACCURACY.  It is the mean of the raw attention at each accepted
        cell, read before the neighbourhood bonus scales that cell's surroundings, so it says how
        strongly the weights preferred the correspondences taken -- not how many of them are right.
        """
        return self._score

    @property
    def unplaced(self) -> tuple[tuple[int, int], ...]:
        """`(product index, atom id)` per product atom left at map number 0.

        A PAIR BECAUSE AN ATOM ID IS UNIQUE WITHIN A CONTAINER AND NOT WITHIN A RECORD: two product
        molecules both hold an atom 1.  The same key :func:`chython.reactions.mapping_agrees` compares
        on.
        """
        return self._unplaced

    @property
    def skipped(self) -> str | None:
        """Why the model was not run, or None when it was.  No mapping is written at all when it is set.

        =============== =================================================================================
        `'empty'`       one of the two sides holds no molecule, so there is no correspondence to find
        `'hypervalent'` an atom carries more than 14 heavy neighbours, outside the domain the weights
                        were trained on
        =============== =================================================================================
        """
        return self._skipped

    def __bool__(self):
        """The mapping changed -- so `if rxn.attention_mapping():` reads as the pass convention does."""
        return self._changed

    def __repr__(self):
        return (f'{type(self).__name__}(changed={self._changed}, score={self._score:.3f}, '
                f'unplaced={self._unplaced}, skipped={self._skipped!r})')


class ReactionContainer:
    """Reactants, products and agents, plus the record's metadata.

    :param reactants: the left side
    :param products: the right side
    :param agents: the middle side -- solvents, catalysts, anything present but not consumed
    :param meta: record metadata, e.g. an RDF's DTYPE/DATUM pairs
    :param title: the record's name line, as bytes

    AN EMPTY REACTION IS A REACTION.  Refusing one -- `ValueError('At least one graph object
    required')` -- leaves a reader meeting an empty `$RFMT` record in a real RDF unable to store it, so
    it has to either drop the record or crash, and dropping a record is the one thing a reader may never
    do.  Input is garbage by default; a container's job is to hold it and let a caller ask for it to
    be repaired.  Refusals belong at the answer boundary, not here.

    THE TYPE CHECK STAYS, and the distinction is deliberate: an empty side is a RECORD, while a
    string in the reactants list is a PROGRAMMING ERROR.  Storing one would defer the failure to
    whichever accessor touched it first, several layers away from the line that caused it.
    """
    __slots__ = ('_reactants', '_products', '_agents', '_meta', '_title', '_log')

    def __init__(self, reactants: Sequence[MoleculeContainer] = (),
                 products: Sequence[MoleculeContainer] = (),
                 agents: Sequence[MoleculeContainer] = (), *,
                 meta: dict | None = None, title='',
                 reagents: Sequence[MoleculeContainer] | None = None):
        if reagents is not None:
            # the chython 2 keyword.  The POSITIONAL third slot is unchanged, so only a caller who
            # spelled it out needs this, and they get told what to spell instead.
            if agents:
                raise TypeError('pass either `agents` or its former spelling `reagents`, not both')
            _renamed('reagents', 'agents')
            agents = reagents

        reactants = tuple(reactants)
        products = tuple(products)
        agents = tuple(agents)
        for side, name in ((reactants, 'reactants'), (products, 'products'), (agents, 'agents')):
            for m in side:
                if not isinstance(m, MoleculeContainer):
                    raise TypeError(f'{name} must hold MoleculeContainer, got {type(m).__name__}')

        self._reactants = reactants
        self._products = products
        self._agents = agents
        self._meta = dict(meta) if meta else None  # lazy: most records have none
        self._log = None  # lazy
        self._title = _as_title(title)

    # --- the three sides -------------------------------------------------------------------------

    @property
    def reactants(self) -> tuple[MoleculeContainer, ...]:
        return self._reactants

    @property
    def products(self) -> tuple[MoleculeContainer, ...]:
        return self._products

    @property
    def agents(self) -> tuple[MoleculeContainer, ...]:
        """The middle side: solvents, catalysts, anything present but not consumed.

        `agents` and not `reagents` for two reasons.  It is the word the formats already use -- MRV
        spells this side `agentList` -- and `reagents` was doing a second job in this library as the
        name of a lookup set of common small molecules, so one word named two things.
        """
        return self._agents

    @property
    def reagents(self) -> tuple[MoleculeContainer, ...]:
        """DEPRECATED -- the former spelling of :attr:`agents`.

        A PURE rename, which is what earns an alias: the same tuple object, the same type, the same
        empty-not-absent behaviour.  Contrast `name`/`title`, which is also a pure rename and still has
        no alias, because chython 2's `name` was SETTABLE and this has no setter -- see
        `chython/core/test/test_alternative_spellings.py`.
        """
        _renamed('reagents', 'agents')
        return self._agents

    def molecules(self) -> Iterator[MoleculeContainer]:
        """Every molecule, reactants then agents then products.

        THE ORDER IS LOAD-BEARING: it is the order the sides' counts are written in, so a serialiser
        that walks this and a reader that trusts the counts have to agree.
        """
        return chain(self._reactants, self._agents, self._products)

    # --- record metadata -------------------------------------------------------------------------

    @property
    def meta(self) -> dict:
        """Record metadata, e.g. an RDF's DTYPE/DATUM pairs.  Created on first access."""
        if self._meta is None:
            self._meta = {}
        return self._meta

    @property
    def log(self) -> Log:
        """What was read, repaired or lost about this reaction.  Created on first access.

        Same property `MoleculeContainer.log` is, unconditional in the same way.  A component's own
        records stay on the component AS WELL: `LogRecord.subject` is how the copy here names which
        molecule it is about, so `rxn.log.by_subject('products[0]')` and `rxn.products[0].log` answer
        the same question from the two ends.
        """
        if self._log is None:
            self._log = Log()
        return self._log

    @property
    def title(self) -> str:
        """The record's name line, as `str`.

        `str` for the same reason the molecule's title is, and by the same handler: see
        :attr:`MoleculeContainer.title`.  `''` for a record never given one.

        Carried faithfully in both directions.  A reader does not tidy it and a writer does not
        repair it, however ugly it is, because a title is what the file said and not what the file
        should have said.
        """
        return self._title

    def set_title(self, title) -> None:
        """Replace the name line.  `set_title` and not a property setter, to match the molecule."""
        self._title = _as_title(title)

    # --- 2D layout, registered by `chython.depict` ------------------------------------------------

    def layout2d(self, *, engine=None, force: bool = False):
        """Arrange this reaction left to right and return `(planes, arrow, signs)`, storing nothing.

        `planes` is one `{n: (x, y)}` per molecule in `molecules()` order, `arrow` is
        `(x1, x2, y)`, and `signs` is one `(x, y)` per `+` between two members of a side.  A member
        that already carries a layout keeps it and is only shifted; `force=True` relays every one.

        Registered by `chython.depict` -- the body is a JavaScript layout engine or a third-party
        toolkit, neither of which this layer knows about.
        """
        _, layout = _reaction_depict_fns()
        return layout(self, engine=engine, force=force)

    def clean2d(self, *, engine=None, force: bool = False):
        """Arrange this reaction, STORE each molecule's plane, and return `(arrow, signs)`.

        THE ARROW AND THE SIGNS ARE RETURNED AND NOT STORED, and that is a decision.  An arrow is a
        property of one drawing at one style, so two figures of one reaction at different bond lengths
        have different arrows and a stored `_arrow`/`_signs` pair would make the second wrong.
        `__slots__` has no field for either, deliberately.  The molecules' coordinates ARE stored,
        because coordinates belong to a molecule.

        Registered by `chython.depict`.
        """
        clean, _ = _reaction_depict_fns()
        return clean(self, engine=engine, force=force)

    # --- drawing, registered by `chython.depict` ---------------------------------------------------

    def scene(self, *, style=None, overlays=None, log=None):
        """This reaction as ONE backend-independent `Scene`: the molecules, the arrow, the `+` signs.

        No `plane` argument, unlike the molecule's: a reaction's picture is an ARRANGEMENT, and a caller
        handing in one plane per member would also be handing in the arrow's span -- which is derived
        from the members' extents and would then have to be an argument too.  `layout2d()` is the door
        for a caller who wants the arrangement without a drawing.

        `overlays` is a `{index: [overlay, ...]}` dict keyed by POSITION in `molecules()`, which yields
        reactants->agents->products in a load-bearing order.

        DRAWING STORES NOTHING: a member with no layout is drawn against a temporary, and the arrow and
        the signs are recomputed per figure because they belong to the drawing and not to the reaction.

        Registered by `chython.depict`.
        """
        _, scene = _reaction_draw_fns()
        return scene(self, style=style, overlays=overlays, log=log)

    def depict(self, *, style=None, overlays=None, log=None):
        """This reaction as an SVG document.  `scene()` plus serialization, and not cached.

        Registered by `chython.depict`.
        """
        depict, _ = _reaction_draw_fns()
        return depict(self, style=style, overlays=overlays, log=log)

    def _repr_svg_(self):
        """Jupyter's hook.  The PROCESS DEFAULT style, because a notebook cell states none."""
        depict, _ = _reaction_draw_fns()
        return depict(self)

    # --- toolkit conversion ----------------------------------------------------------------------

    def to_rdkit(self, **kwargs):
        """This reaction as an RDKit `rdChemReactions.ChemicalReaction`.

        THE ONLY TOOLKIT METHOD ON A REACTION, because it is the only one of the five with a reaction
        form at all: `to_indigo` and the rest are molecule methods, and a reaction has no Indigo,
        OpenBabel, CDK or CDPKit shape to be converted into.  `interop.indigo(rxn)` says so with the
        message the converter itself raises; there is no method here to make the gap look smaller.

        The three sides go to the three template lists in `molecules()` order, and `keep_mapping=True`
        -- the default -- carries `map_number` into RDKit's atom map field, so a mapped record survives
        the round trip back through `interop.rdkit`.  Keywords are `to_rdkit`'s, forwarded per molecule.

        Registered by `chython.interop`, not implemented here -- see `_set_interop_fns`.
        """
        return _reaction_interop_fn('rdkit')(self, **kwargs)

    # --- copies ----------------------------------------------------------------------------------

    def copy(self) -> 'ReactionContainer':
        """A deep copy: every molecule is copied, metadata and title come along."""
        copy = object.__new__(type(self))
        copy._reactants = tuple(m.copy() for m in self._reactants)
        copy._products = tuple(m.copy() for m in self._products)
        copy._agents = tuple(m.copy() for m in self._agents)
        copy._meta = None if self._meta is None else self._meta.copy()
        copy._log = None  # per handle, exactly as the molecule's is
        copy._title = self._title
        return copy

    # --- the ML view -----------------------------------------------------------------------------

    def modeling_view(self) -> ReactionModelingView:
        """Per-atom transition state keyed by map number, and the union's bond orders per side.

        A dict assembly over `reaction_transition_view`'s arrays, so there is ONE union derivation.
        The encoding passes `unknown_h=H_UNKNOWN`: this view reports an underivable hydrogen count as
        the sentinel, and the ML default of 0 is for a trained vocabulary rather than a chemical
        answer.

        AGENTS ARE NOT IN IT.  An agent is by definition not consumed, so it contributes the same
        state to both sides and no bond change; including it would add tokens carrying no signal and
        would make the atom count depend on how the record's author chose to split the left side.
        Folding agents into the reactant side is a different decision for a different purpose -- a CGR
        that is also a drawing.  Stated rather than assumed, because it is the kind of choice a model
        silently trains around.

        AN UNMAPPED ATOM IS KEYED NEGATIVELY, `-1` onwards in union order, and still counted in
        `unmapped`.  It holds a union row like any other atom -- on the side it came from -- so the key
        is only the identity handle a dict needs: 0 would put every such atom on one entry, and a
        positive one would be indistinguishable from a number the record stated.

        A COLLIDING MAP NUMBER KEEPS ITS FIRST CLAIM.  The second atom to claim a number is counted in
        `collisions` and contributes no state and no bond; a union cannot hold two atoms at one key.
        """
        view = reaction_transition_view(self._reactants, self._products,
                                        TensorEncoding(unknown_h=H_UNKNOWN))
        keys, bare = [], 0
        for mn in view.map_numbers.tolist():
            if not mn:
                bare -= 1
                mn = bare
            keys.append(mn)
        states = {mn: (int(view.elements[i]), int(view.h_before[i]), int(view.n_before[i]),
                       int(view.h_after[i]), int(view.n_after[i]))
                  for i, mn in enumerate(keys)}
        union = {}
        for (i, j), before, after in zip(view.bonds.tolist(), view.bond_before, view.bond_after):
            a, b = keys[i], keys[j]
            union[(a, b) if a < b else (b, a)] = (int(before), int(after))
        return ReactionModelingView(states, union, view.unmapped, view.collisions)

    def transition_view(self, encoding=None):
        """Per-atom state on each side of the reaction, over the mapped union, as int32 arrays.

        Agents contribute nothing: only the two sides being unioned are handed to the kernel. A badly
        mapped record is reported through `unmapped` and `collisions`, never refused. See `docs/ml.rst`.
        """
        return reaction_transition_view(self._reactants, self._products, encoding)

    # --- the standardization passes ---------------------------------------------------------------
    #
    # Ten one-line forwards to `chython/core/_reaction_passes.py`, which is where their reasoning is.
    # The two things every one of them shares:
    #
    # THEY MUTATE IN PLACE AND ANSWER "DID ANYTHING CHANGE?", never a new reaction.  A reaction record
    # arrives from a file, is repaired, and is written or stored; copying three sides on every pass in
    # a loop over a million records would cost more than every pass put together.  A caller who wants
    # the original keeps `rxn.copy()`.
    #
    # NONE OF THEM TAKES A `log=`.  Every record goes to `self.log`, unconditionally, and to the
    # component's own `.log` as well: `rxn.log.by_subject('products[0]')` and `rxn.products[0].log`
    # answer the same question from the two ends.  `subject` is why the reaction-level copy means
    # anything -- a record's `atoms` are stable ids in ONE container and say nothing without it.

    def standardize(self, *, fix_hydrogens: bool = True, fix_tautomers: bool = True) -> bool:
        """Run the functional-group repair table over every molecule on every side."""
        return reaction_standardize(self, fix_hydrogens=fix_hydrogens, fix_tautomers=fix_tautomers)

    def canonicalize(self, *, fix_tautomers: bool = True, keep_kekule: bool = False) -> bool:
        """The full repair sequence on every molecule: standardize, hydrogens, aromatic form.

        NO MAPPING REPAIR: this does not call `fix_mapping`, the table of rules that repairs a
        mis-drawn atom-to-atom mapping.  Mapping repair is the mapping epic's, it is not a chemistry
        pass, and a caller asking for a canonical representation has not asked for their mapping to be
        rewritten.  The two are one call apart when both are wanted.
        """
        return reaction_canonicalize(self, fix_tautomers=fix_tautomers, keep_kekule=keep_kekule)

    def reconstruct_mapping(self, *, max_size_ratio: float = 5.,
                            min_filter_size: int = 42) -> tuple[str, ...]:
        """Assign an atom-atom mapping by reconstructing the recorded product from the recorded inputs.

        Canonicalizes both sides IN PLACE and then numbers them, so the reaction comes back normalized
        as well as mapped -- the two are one operation because a mapping written over an un-normalized
        record is a mapping of a structure the caller is about to change.  Returns a 1-tuple with the
        label of the explanation that was applied (`'react:amidation'`, `'deprotect:amine_boc'`,
        `'purification'`), or empty when nothing explained the record.

        When no rung explains the record, a ``LOST`` log line with rule ``'reconstruct:unexplained'``
        is emitted and the product's map numbers are cleared and not replaced -- the product arrives
        with its incoming numbers gone and none derived in their place.

        REFUSES RATHER THAN GUESSES, at the answer boundary where a refusal belongs: a multi-product
        record, a record with no inputs or no products, and a record whose product is grossly larger
        than everything that went in all come back empty with a `REFUSED` record on `self.log`.

        Registered by `chython.reactions`, not implemented here -- see `_reaction_reconstruct_fn`.
        """
        return _reaction_reconstruct_fn()(self, max_size_ratio=max_size_ratio,
                                          min_filter_size=min_filter_size)

    def attention_mapping(self, *, multiplier: float = 1.75, keep_reactant_mapping: bool = False,
                          threads: int | None = None) -> MappingResult:
        """Assign an atom-atom mapping from a transformer's attention over the two sides.

        Writes map numbers IN PLACE and touches nothing else: atom ids, bonds, charges and the
        structures themselves come back as they went in.  Reactant atoms take 1..N in container order,
        each product atom takes the number of the reactant atom the model matched it to, a product atom
        the model matched to nothing keeps 0, and agents are numbered last.

        `keep_reactant_mapping=True` leaves the reactant side's existing numbers alone and numbers the
        products against them -- for a record whose inputs are already mapped by something else.

        AGENTS ARE NUMBERED AND NEVER MODELLED.  Only the reactants and the products are encoded; the
        weights were trained without a catalyst on the input side.

        NO RULE-BASED REPAIR RUNS AFTERWARDS.  This method is the model and nothing else, so that its
        own accuracy is a number a caller can obtain -- composing it with a fixer is the caller's next
        line.

        `multiplier` scales the attention around an accepted correspondence, biasing the next choice
        towards a neighbour of it.  `threads` is the ONNX Runtime intra-op thread count, defaulting to
        `min(cpu_count(), 8)`; a second value costs a second loaded model.

        Two records are declined rather than mapped, and :attr:`MappingResult.skipped` says which: an
        empty side and an atom past 14 heavy neighbours.  Both leave every map number as it was.

        Needs `chython[mapping]` -- the runtime and the weights are an extra, and the weights are their
        own 80 MiB distribution.  Registered by `chython.reactions`, not implemented here -- see
        `_reaction_attention_fn`.
        """
        return _reaction_attention_fn()(self, multiplier=multiplier,
                                        keep_reactant_mapping=keep_reactant_mapping, threads=threads)

    def kekule(self) -> bool:
        """Give every aromatic ring on every side an alternating-bond form."""
        return reaction_kekule(self)

    def thiele(self) -> bool:
        """Find the aromatic rings on every side and mark them aromatic."""
        return reaction_thiele(self)

    def neutralize(self, *, keep_charge: bool = True) -> bool:
        """Move every proton the acid/base table can from a cation onto an anion, on every side.

        PER MOLECULE, NOT PER SIDE and never across the arrow: a reactant written `C[NH3+].[Cl-]` is one
        container and neutralizes, while a chloride written as a separate reactant does not pair with
        anything -- pairing two ions listed side by side is `contract_ions`' question.
        """
        return reaction_neutralize(self, keep_charge=keep_charge)

    def explicify_hydrogens(self) -> int:
        """Make every implicit hydrogen an atom.  How many were added?

        THE MAP NUMBERS ARE PAIRED ACROSS THE ARROW: a hydrogen added to a mapped atom on the left and
        one added to the atom with the same number on the right are given the SAME new number, because
        they are the same hydrogen and a mapping that numbered them differently would claim a C-H bond
        broke and an identical one formed.  Hydrogens in an unmapped molecule are left unmapped.
        """
        return reaction_explicify_hydrogens(self)

    def implicify_hydrogens(self) -> int:
        """Fold every hydrogen atom that can be a count back into a count.  How many were removed?"""
        return reaction_implicify_hydrogens(self)

    def clean_isotopes(self) -> bool:
        """Drop every isotope label on every side, in place.  Did the reaction carry one?

        PER MOLECULE, and the parity a dropped label was the only justification for goes with it --
        `MoleculeContainer.clean_isotopes` borrows `validate_stereo` for that, once per molecule, so no
        atom is judged against a constitution from another side of the arrow.
        """
        return reaction_clean_isotopes(self)

    def clean_stereo(self) -> dict[str, dict]:
        """Wipe every kind of stereo state on every side, unconditionally.  Returns what was wiped.

        KEYED BY LOCATION -- `{'reactants[0]': {'parities': [2]}, ...}`, the same string a log record's
        `subject` carries, so a key is also how the molecule is addressed again.  Each value is
        `MoleculeContainer.clean_stereo`'s own report unchanged, and a molecule that carried no stereo
        is absent rather than present and empty, so `{}` means the reaction had none anywhere.
        """
        return reaction_clean_stereo(self)

    def remove_reagents(self, *, keep_reagents: bool = False, mapping: bool = True,
                        common: Sequence[MoleculeContainer] | None = None) -> bool:
        """Move the molecules that are not part of the transformation out of reactants and products.

        With `keep_reagents` they become agents; without it they are dropped.  `mapping=True` reads the
        atom-to-atom mapping and moves whatever the reaction centre does not touch, and raises
        `ValueError` when the record has no mapping to read.  `mapping=False` takes the rule-based door:
        a molecule appearing on both sides, plus anything in `common`, which is a LIST THE CALLER PASSES
        because a table of solvents is chemistry knowledge and does not belong in `core`.

        Neither door will empty a side.  A reaction whose every reactant looks like a reagent is a
        record this pass cannot improve, and it is returned unchanged with `False`.
        """
        return reaction_remove_reagents(self, keep_reagents=keep_reagents, mapping=mapping,
                                        common=common)

    def contract_ions(self) -> bool:
        """Join the free ions on each side into salts, when which pairs with which is determined.

        `[Na+].[OH-]` on one side becomes one molecule of two components.  Two different cations and
        one anion do not, because nothing in the record says which the anion belongs to -- a refusal
        to guess, reported as `False`, not an error.
        """
        return reaction_contract_ions(self)

    def reset_mapping(self) -> bool:
        """Number every atom in the reaction 1..N, from a single counter, in `molecules()` order.

        Only when the numbering is not already a unique numbering of every atom -- a record that
        arrived correctly mapped is left exactly as it is, since renumbering it would destroy a real
        atom-to-atom mapping to fix nothing.  This is a reaction method and not a molecule one
        precisely because the counter has to span the sides.
        """
        return reaction_reset_mapping(self)

    # --- the wire format -------------------------------------------------------------------------

    def pack(self, *, compressed=True, drop=None, version=None) -> bytes:
        """One reaction pach record.  See :func:`reaction_pach_dump`, which this forwards to.

        `version` is None for the current record, 5, or 1 for the legacy one.  Version 1 moves each
        molecule's map numbers into pach's atom-number field and refuses a partially mapped molecule;
        version 5's molecule records carry the field themselves.

        `pack`/`unpack`/`drop` and no `check=`: the naming owes its consistency to the molecule side
        of this release, where `MoleculeContainer.pack` takes `drop=` and refuses by field name.
        """
        return reaction_pach_dump(self, compressed=compressed, drop=drop, version=version)

    @staticmethod
    def unpack(data, *, compressed=None) -> 'ReactionContainer':
        """A reaction from a reaction pach record.  `compressed` defaults to sniffing.

        THIS IS AN ANSWER BOUNDARY AND IT RAISES, exactly as `MoleculeContainer.unpack` does:
        `ValueError` names what was wrong with the record, with every problem the decoder found
        appended -- including the ones it recovered from, since a caller who cannot have a reaction is
        owed the whole story.  A caller walking a store who wants the complaints instead of an
        exception wants :func:`reaction_pach_load`.

        THERE IS NO `to_bytes` COUNTERPART AND SO NO `__bytes__`.  chython 2 spelled `bytes(rxn)` as
        `rxn.pack()`, but on this release `bytes(mol)` is the ARENA and not pach -- so a `__bytes__`
        here meaning pach would make the same expression mean two different formats one layer apart.
        A reaction has no arena of its own to return instead, so the spelling is simply absent.
        """
        rxn, problems = reaction_pach_load(data, compressed=compressed)
        if rxn is None:
            raise ValueError('this is not a readable reaction pach record: %s' % '; '.join(problems))
        if problems:
            raise ValueError('this reaction pach record is damaged: %s. reaction_pach_load() returns '
                             'the reaction that could be recovered from it along with these problems'
                             % '; '.join(problems))
        return rxn

    def pach(self, *, compressed=True, drop=None, version=None) -> bytes:
        """chython 2's name for `pack`, and the same record byte for byte.

        chython 2's `check=` and `order=` are a `TypeError` rather than accepted and ignored, exactly
        as on `MoleculeContainer.pach`.
        """
        return reaction_pach_dump(self, compressed=compressed, drop=drop, version=version)

    @staticmethod
    def unpach(data, *, compressed=None) -> 'ReactionContainer':
        """chython 2's name for `unpack`, with its behaviour: an answer boundary that raises."""
        return ReactionContainer.unpack(data, compressed=compressed)

    @staticmethod
    def pack_len(data, *, compressed=None) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
        """Each molecule's ATOM COUNT, per side: `(reactants, agents, products)`, without decoding.

        The counts sit in each molecule record's own header, so this walks the stream and reads them --
        which is what a caller sizing a batch, or picking records by size out of a column of stored
        buffers, wants instead of unpacking several thousand reactions.

        An answer boundary: it returns numbers and has no way to say "unknown", so it raises
        `ValueError` where `reaction_pach_load` would report.  An empty side is read from the counts and
        never by slicing -- `molecules[-products:]` with no products takes the whole list.
        """
        problems: list[str] = []
        raw = _pach_reaction_raw(data, compressed, problems)
        if raw is None:
            raise ValueError(problems[0])
        counts = (raw[1], raw[2], raw[3])
        atoms: list[int] = []
        shift = 4
        for index in range(sum(counts)):
            if shift + 4 > len(raw):
                raise ValueError('the header declares %d molecules and the buffer ends after %d'
                                 % (sum(counts), index))
            atoms.append(_pach_molecule_atom_count(raw, shift))
            shift += pach_record_length(raw[shift:], compressed=False)
        if shift > len(raw):
            raise ValueError('the declared molecule records overrun the buffer by %d byte(s)'
                             % (shift - len(raw)))
        first, second = counts[0], counts[0] + counts[1]
        return (tuple(atoms[:first]), tuple(atoms[first:second]), tuple(atoms[second:]))

    # --- protocol --------------------------------------------------------------------------------

    def __len__(self):
        return len(self._reactants) + len(self._agents) + len(self._products)

    def __bool__(self):
        """True when there is both a left and a right side -- i.e. something actually happens.

        NOT `__len__ != 0`, and the difference is the point: a record with reactants and no products
        is storable (see the class docstring) and is not a reaction.  `len()` answers how many
        molecules are held; `bool()` answers whether the record describes a transformation.
        """
        return bool(self._reactants and self._products)

    def __repr__(self):
        """A STRUCTURAL SUMMARY, and deliberately still one now that `__str__` writes a real SMILES.

        `repr` answers "what am I holding" -- three counts and a title -- where `str` answers "what
        reaction is this".  Both questions get asked, and a debugger printing a hundred reactions
        wants the counts, not a hundred canonical searches.
        """
        return (f'{type(self).__name__}({len(self._reactants)} reactants, {len(self._agents)} agents,'
                f' {len(self._products)} products, title={self._title!r})')

    def __str__(self):
        """The reaction SMILES, which for a reaction IS the chemical identifier.

        One writer behind all three doors -- this, `__format__` and `.smiles` -- so a string cannot
        depend on which one a caller reached for.
        """
        return write_reaction_smiles(self)

    def __format__(self, format_spec):
        """`format(rxn, spec)`, with the writer's spec keys.

        `!c` keeps the container's order of molecules within each side; the default sorts each side by
        the molecules' own strings, which is what makes the result an identifier.  Every other key --
        `a`, `!s`, `A`, `m`, `h`, `!b`, `!x`, `!z` -- goes to each molecule's `write_smiles` unchanged,
        so a spec that means something per molecule means the same thing here.
        """
        return write_reaction_smiles(self, format_spec)

    @property
    def smiles(self) -> str:
        """The reaction SMILES with default options, beside `MoleculeContainer.smiles`."""
        return write_reaction_smiles(self)

    def _identity(self):
        """The reaction's identity: one sorted tuple of molecule identities PER SIDE.

        Three tuples and not one, because a side is part of the chemistry: moving a molecule from
        the reactant side to the product side is a different reaction, and a single pooled multiset
        would call the two equal.

        SORTED, so it is a multiset and not a sequence.  Reactant ORDER is not chemistry -- `A + B`
        and `B + A` are one reaction -- but multiplicity is: `2 A -> B` is not `A -> B`, which is
        what rules out a `frozenset` here.  The sort key is `canonical_bytes`, the molecule's own
        identity, so it is a function of the molecules and not of the order they were added in.

        Sorting on the BYTES and not on `hash(mol)`: `hash` of a `bytes` is salted per interpreter
        under PYTHONHASHSEED, so a hash-ordered tuple would differ between processes.  Nothing here
        is persisted, so that would not be a wrong answer today, but it would make this value
        untrustworthy the moment anyone wrote it down, and the bytes cost nothing extra -- they are
        already cached on each molecule.

        Raises `AutomorphismBudgetExceeded` from any molecule whose canonical search truncates, for
        the reason `canonical_bytes` gives: there is no degraded identity.
        """
        return (tuple(sorted(m.canonical_bytes for m in self._reactants)),
                tuple(sorted(m.canonical_bytes for m in self._agents)),
                tuple(sorted(m.canonical_bytes for m in self._products)))

    def __eq__(self, other):
        """Equal when the three sides hold the same molecules with the same multiplicities.

        WHAT PARTICIPATES, and each of these was decided rather than fallen into -- the module
        docstring carries the reasoning and this is the summary:

          * the three SIDES, separately.  AGENTS INCLUDED (ruled 2026-09-03): a reaction is a record
            of what was done, so the same transformation run with and without a named catalyst is two
            reactions.  A caller who wants the transformation alone compares
            `(reactants, products)` -- which is what `_identity`'s first and third element are, so it
            costs nothing to spell;
          * each molecule's own identity, which carries constitution, isotopes, charges, radicals,
            implicit hydrogen counts, bond orders and configured parities.

        WHAT DOES NOT, and both absences are load-bearing:

          * ATOM-TO-ATOM MAPPING (ruled 2026-09-03: "mapping is not needed for reaction
            comparison").  This is free rather than arranged: a map number is not in the atom
            invariant word `mol_identity_bytes` reads, so two mappings of one reaction already
            compare equal molecule by molecule.  The consequence the mapping epic must design
            around is that "did we reproduce this mapping" CANNOT be spelled `a == b`; it is a
            question about an annotation and needs its own orbit-aware comparison;
          * `title` and `meta`.  Not chemistry.  A record read from a file and the same record read
            from a different file with a different name are one reaction.

        Returns `NotImplemented` for a non-reaction rather than `False`, so Python can try the
        other operand's `__eq__` and `!=` stays consistent with it.
        """
        if not isinstance(other, ReactionContainer):
            return NotImplemented
        if self is other:
            return True
        # Cheap discriminators first: the identity of even one side is a canonical search per
        # molecule, and unequal side sizes are the common case among unequal reactions.
        if (len(self._reactants) != len(other._reactants)
                or len(self._agents) != len(other._agents)
                or len(self._products) != len(other._products)):
            return False
        return self._identity() == other._identity()

    def __hash__(self):
        """Hashes exactly what `__eq__` compares, so the two cannot drift apart.

        Defined explicitly because a class that defines `__eq__` and not `__hash__` is unhashable in
        Python 3 -- and an unhashable reaction is the thing this whole decision was for.
        """
        return hash(self._identity())


def _resolve_drop(drop):
    """`drop=` to a frozenset of names, refusing an unrecognised one.

    Refused rather than ignored for the reason `pach_dump` gives: a misspelt waiver that silently did
    nothing would turn back into a raise on some later record.
    """
    if drop is None:
        return frozenset()
    if drop == '*':
        return _PACH_DROP_NAMES
    names = frozenset(drop)
    unknown = names - _PACH_DROP_NAMES
    if unknown:
        raise ValueError('%s is not a droppable field; the drop names are %s'
                         % (', '.join(repr(n) for n in sorted(unknown)),
                            ', '.join(sorted(_PACH_DROP_NAMES))))
    return names


def _pach_molecule_record(mol, index, drop):
    """One molecule's pach bytes with its map numbers moved into the atom-number field."""
    numbers = [(a.n, a.map_number) for a in mol.atoms()]
    mapped = [n for _, n in numbers if n]
    if mapped and 'map_number' not in drop:
        if len(mapped) != len(numbers):
            raise ValueError('molecule %d is PARTIALLY mapped -- %d of its %d atoms carry a '
                             'map_number -- and the reaction pach format has one number field per '
                             'atom, so there is no spelling for that; map the rest, or pass '
                             'drop=[\'map_number\'] to write the record without the mapping'
                             % (index, len(mapped), len(numbers)))
        if len(set(mapped)) != len(mapped):
            raise ValueError('molecule %d has two atoms sharing a map_number, and the reaction pach '
                             'format keys an atom by that number; pass drop=[\'map_number\'] to '
                             'write the record without the mapping' % index)
        high = max(mapped)
        if high > _PACH_MAX_NUMBER:
            raise ValueError('molecule %d carries map_number %d and the pach atom-number field is 12 '
                             'bits, so it holds 1..%d; pass drop=[\'map_number\'] to write the record '
                             'without the mapping' % (index, high, _PACH_MAX_NUMBER))
        mol = mol.copy()
        mol.remap(dict(numbers))
    # 'map_number' unconditionally: the relabelled copy still carries the field, and an unmapped
    # molecule has nothing there for the waiver to waive.
    return pach_dump(mol, compressed=False, drop=sorted(drop | {'map_number'}), version=2)


def reaction_pach_dump(rxn: 'ReactionContainer', *, compressed=True, drop=None, version=None) -> bytes:
    """Write one reaction pach record.  `ReactionContainer.pack` forwards here.

    `version` is `None` for the current record, 5, and `1` for the legacy one.  Version 5's molecules
    carry their own map numbers; version 1 moves them into the atom-number field, which has no
    spelling for a partially mapped molecule and refuses one by name.

    Raises `ValueError` naming any field the container holds and the format cannot -- the reaction's
    `meta` and `title`, and everything `pach_dump` refuses on each molecule.  `drop` waives those:
    an iterable of names, or `'*'` for all of them.  The names are `pach_dump`'s nine plus `meta`.

    Also raises when a side holds more than 255 molecules, because the count field is a uint8.
    """
    names = _resolve_drop(drop)
    if 'meta' not in names and rxn._meta:
        raise ValueError('the reaction carries %d metadata key(s) and the pach format has no field '
                         'for any of them; pass drop=[\'meta\'] to write the record without them'
                         % len(rxn._meta))
    if 'title' not in names and rxn._title:
        raise ValueError('the reaction carries the title %r and the pach format has no text of any '
                         'kind; pass drop=[\'title\'] to write the record without it' % rxn._title)
    molecule_drop = names & _PACH_MOLECULE_DROP_NAMES

    counts = []
    for side, name in ((rxn._reactants, 'reactants'), (rxn._agents, 'agents'),
                       (rxn._products, 'products')):
        if len(side) > _PACH_SIDE_MAX:
            raise ValueError('the %s side holds %d molecules and the reaction pach count field is a '
                             'uint8, so it holds at most %d; the format cannot store this reaction'
                             % (name, len(side), _PACH_SIDE_MAX))
        counts.append(len(side))

    if version is None:
        version = _PACH_REACTION_VERSION_5
    elif version.__class__ is not int or version not in _PACH_REACTION_VERSIONS:
        raise ValueError('%r is not a writable reaction pach version; they are 1, 5 and None for 5'
                         % (version,))

    out = bytearray((version, counts[0], counts[1], counts[2]))
    for index, mol in enumerate(rxn.molecules()):
        if version == _PACH_REACTION_VERSION_1:
            out += _pach_molecule_record(mol, index, molecule_drop)
        elif version == _PACH_REACTION_VERSION_5:
            out += pach_dump(mol, compressed=False, drop=sorted(molecule_drop))
        else:
            raise ValueError('%r is in the writable version set but has no writer; the set and the '
                             'dispatch must move together' % (version,))
    if compressed:
        return zlib.compress(bytes(out), 9)
    return bytes(out)


def _pach_molecule_atom_count(raw, shift):
    """The atom count out of a molecule record's header: 12 bits split across two bytes in versions 0
    and 2, a little-endian uint16 in versions 3 and 4."""
    if raw[shift] == 3 or raw[shift] == 4:
        return raw[shift + 2] | (raw[shift + 3] << 8)
    return (raw[shift + 1] << 4) | (raw[shift + 2] >> 4)


def _pach_reaction_raw(data, compressed, problems):
    """The buffer as raw reaction pach bytes, or None with a sentence saying why not.

    The sniff is exact rather than heuristic, the same way the molecule side's is: a raw record's first
    byte is 1 or 5, and a zlib header's low nibble is its compression method, always 8, so neither is a
    byte a zlib header can begin with.
    """
    raw = bytes(data)
    if not len(raw):
        problems.append('the buffer is empty; a reaction pach record is at least a 4 byte header')
        return None
    looks_raw = raw[0] in _PACH_REACTION_VERSIONS
    if compressed is True and looks_raw:
        problems.append('compressed=True was stated and the buffer begins with %d, which is a '
                        'reaction pach version, so it is a raw record' % raw[0])
        return None
    if compressed is False and not looks_raw:
        problems.append('compressed=False was stated and the buffer begins with %d, which is not a '
                        'reaction pach version' % raw[0])
        return None
    if not looks_raw:
        try:
            raw = zlib.decompress(raw)
        except Exception as err:
            problems.append('the buffer begins with %d, so it is neither a raw reaction pach record '
                            'nor a readable zlib stream: %s' % (raw[0], err))
            return None
    if len(raw) < 4:
        problems.append('a reaction pach record is at least a 4 byte header and this buffer is %d '
                        'byte(s)' % len(raw))
        return None
    if raw[0] not in _PACH_REACTION_VERSIONS:
        problems.append('byte 0 is %d, which is not a reaction pach version; they are 1 and 5'
                        % raw[0])
        return None
    return raw


def reaction_pach_load(data, *, compressed=None):
    """Read one reaction pach record, version 1 or 5.  `(ReactionContainer or None, problems)`; never raises.

    The loop-safe door, mirroring `pach_load`: a store of forty thousand reaction records must not be
    stopped by one of them, so this reports what was wrong instead of raising.  `ReactionContainer.
    unpack` is the answer boundary and raises.

    A PARTIAL REACTION IS NEVER RETURNED.  Where `pach_load` hands back the molecule it could recover,
    this hands back `None` as soon as one declared molecule is missing or unreadable, because a
    reaction is a relation between its sides: a side with a hole in it is not a smaller reaction, it is
    a wrong one, and a caller comparing mappings would have no way to notice.  The molecules' own
    recoverable damage IS carried through -- a bit-flipped bond in a readable record gives a reaction
    and a non-empty `problems`.

    `compressed` defaults to sniffing; `True` and `False` state it instead.
    """
    problems: list[str] = []
    raw = _pach_reaction_raw(data, compressed, problems)
    if raw is None:
        return None, problems

    reactants, agents, products = raw[1], raw[2], raw[3]
    total = reactants + agents + products
    molecules = []
    shift = 4
    for index in range(total):
        if shift >= len(raw):
            problems.append('the header declares %d molecules and the buffer ends after %d'
                            % (total, index))
            return None, problems
        try:
            length = pach_record_length(raw[shift:], compressed=False)
        except ValueError as err:
            problems.append('molecule %d, at byte %d, is not a measurable pach record: %s'
                            % (index, shift, err))
            return None, problems
        if shift + length > len(raw):
            problems.append('molecule %d, at byte %d, declares a %d byte record and only %d byte(s) '
                            'are left' % (index, shift, length, len(raw) - shift))
            return None, problems
        mol, mol_problems = pach_load(raw[shift:shift + length], compressed=False)
        problems.extend('molecule %d: %s' % (index, p) for p in mol_problems)
        if mol is None:
            problems.append('molecule %d, at byte %d, could not be read at all' % (index, shift))
            return None, problems
        # VERSION 1 ONLY: its writer put the map numbers into the atom-number field, so its reader takes
        # them out again.  A version 5 record carries the field, and re-deriving it from atom numbers
        # would overwrite a real mapping with a renumbering.  A new version decides for itself whether
        # its reader restores, so it adds its own branch rather than inheriting this fall-through.
        if raw[0] == _PACH_REACTION_VERSION_1:
            _restore_map_numbers(mol)
        molecules.append(mol)
        shift += length
    if shift != len(raw):
        problems.append('the %d declared molecule(s) end at byte %d and the buffer is %d bytes, so '
                        '%d trailing byte(s) were ignored' % (total, shift, len(raw), len(raw) - shift))
    return (ReactionContainer(molecules[:reactants], molecules[reactants + agents:],
                              molecules[reactants:reactants + agents]),
            problems)


def _restore_map_numbers(mol):
    """Set every atom's map number to the number its pach record gave it.

    VERSION 1 ONLY.  See the module's format notes for why that version's reader must do this and why a
    version 5 reader must not.  One `edit()` block for the whole molecule, so the arena is rebuilt once
    rather than once per atom.
    """
    ids = list(mol.atom_numbers)
    with mol.edit():
        for n in ids:
            mol.set_map_number(n, n)


def _renamed(old, new):
    """Announce a superseded spelling, naming what replaced it.

    `stacklevel=3` -- `warn` -> here -> the property -> the caller.  MEASURED AND ASSERTED, not
    reasoned: the same intent needs 1 inside the compiled core, because neither a `cdef` helper nor a
    compiled `def` pushes a Python frame.  The right number is a property of the call chain, so a
    test pins the blamed line.
    """
    warn(f'`{old}` was renamed to `{new}` and will be removed in a later release; use `{new}`',
         DeprecationWarning, stacklevel=3)


# INJECTION, not inheritance -- the same shape `_ich_set_kekule_fn` uses, and for the same reason: the
# reaction SMILES reader lives in the extension, this class does not, and the extension cannot import
# upwards.  So the reader is told what to build rather than knowing it.
_set_reaction_factory(ReactionContainer)


__all__ = ['MappingResult', 'ReactionContainer', 'ReactionModelingView', 'reaction_pach_dump',
           'reaction_pach_load']
