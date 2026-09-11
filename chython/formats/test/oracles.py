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
"""Four reference toolkits as format oracles: RDKit, Indigo, CDK and ChemAxon Marvin.

WHY THIS MODULE EXISTS.  A file format is only useful if two programs agree about it, so every
conformance claim in :mod:`.test_conformance` is measured against a toolkit that is not this one.  The
matrix is (format x feature x toolkit x direction) and BOTH directions are asserted: a reference writer
emits the feature and chython must recover it, and chython emits the feature and the reference reader
must recover it.  The second is the release-critical half.

THE ASYMMETRY, taken from :mod:`chython.core.test.oracle`: an **absent** oracle SKIPS -- a machine
without Marvin still runs the whole suite green -- while a **present but unusable** one FAILS loudly,
because a provisioned oracle that answers nothing turns a differential into a tautology.

VERSIONS ARE PINNED IN THE REPORT.  Every failure message and the ``matrix`` verb carry the version
string each oracle reports for itself, so a difference is attributable to a named release rather than to
"some toolkit".  Nothing here judges another implementation: a cell records what a toolkit can and
cannot represent, and where the specification is silent it says so and picks no winner.

PROVISIONING.  ``rdkit``, ``indigo`` (with ``indigo.inchi``) and ``jpype`` come from pip.  The CDK jar
and Marvin are not shipped and not in the tree:

    export CHYTHON_CDK_JAR=/path/to/cdk-2.12.jar
    export MARVIN_BIN=/Applications/MarvinSuite/bin/molconvert

Both have defaults that are allowed to be absent.  ``python -m chython.formats.test.oracles`` reports
which of the four this machine has.
"""

from functools import cache
from importlib.util import find_spec
from io import StringIO
from os import environ
from pathlib import Path
from subprocess import run
from sys import exit as sys_exit
from typing import NamedTuple

from pytest import skip

from ..ctfile import RDFRead, SDFRead, SDFWrite, mol, rxn
from ..mol2 import read_mol2
from ..pdb import read_mmcif, read_pdb
from ..xml import read_cml, read_mrv
from ..xyz import xyz
from ..ctfile import SGroup, SGroupStore
from ..xml import write_cml, write_mrv
from ...core import (inchi_to_molecule, molecule_to_inchi, read_reaction_smiles, read_smiles,
                     write_reaction_smiles, write_smiles)


__all__ = ['PINS', 'CDK_JAR', 'MOLCONVERT', 'TOOLKITS',
           'rdkit_version', 'indigo_version', 'cdk_version', 'marvin_version',
           'cdk', 'require', 'versions', 'main',
           'Parsed', 'Unrepresentable', 'READERS', 'READS', 'skeleton_of', 'read_chython',
           'FEATURE_NAMES', 'FEATURE_FORMATS', 'PROBES', 'CHYTHON_PROBES', 'RDKIT_PROBES',
           'INDIGO_PROBES', 'CDK_PROBES', 'MARVIN_PROBES', 'Seed', 'SEEDS', 'SEED_KEYS',
           'seed_of', 'is_reaction', 'WRITERS', 'write_chython', 'write_once',
           'attach_dat_sgroup', 'build_chython', 'parse_once', 'prewarm_marvin',
           'seed_formats', 'warm', 'prewarm_marvin_writes', 'SD_FIELDS', 'ALIAS_TEXT', 'container_of',
           'molconvert']


#: The four oracles, in report order.
TOOLKITS = ('rdkit', 'indigo', 'cdk', 'marvin')

#: Env var per out-of-tree oracle.  `rdkit` and `indigo` are importable or not, so they have none.
PINS = {'cdk': 'CHYTHON_CDK_JAR', 'marvin': 'MARVIN_BIN'}


def _root():
    """The repository root of this checkout: `chython/formats/test/oracles.py` is four levels down."""
    return Path(__file__).resolve().parents[3]


def _find_cdk_jar():
    """The CDK jar, from the env var or from a `java/` directory at or above this checkout.

    `java/` is untracked, so a git worktree of this repository does not have one while the main
    checkout beside it does -- hence the walk up the ancestors rather than a single fixed path.
    """
    named = environ.get('CHYTHON_CDK_JAR') or environ.get('CDK_PATH')
    if named:
        return Path(named)
    root = _root()
    for base in (root, *root.parents[:4]):
        found = sorted((base / 'java').glob('cdk-*.jar'))
        if found:
            return found[-1]
    return root / 'java' / 'cdk-2.12.jar'


#: The committed fixture directory: bytes for the formats no oracle here writes.
DATA = _root() / 'test'

#: Where each out-of-tree oracle is found.  Neither default is shipped and both may be absent.
CDK_JAR = _find_cdk_jar()
MOLCONVERT = Path(environ.get('MARVIN_BIN') or '/Applications/MarvinSuite/bin/molconvert')


@cache
def rdkit_version():
    """`rdkit.__version__`, or None.  `find_spec` first so collection never pays for the import."""
    if find_spec('rdkit') is None:
        return None
    import rdkit
    return rdkit.__version__


@cache
def indigo_version():
    """Indigo AND its InChI plugin, because the skeleton probe needs both.

    The distribution is imported as `indigo`, not `epam.indigo`.
    """
    if find_spec('indigo') is None:
        return None
    from indigo import Indigo
    from indigo.inchi import IndigoInchi
    i = Indigo()
    return f'{i.version()} + {IndigoInchi(i).version()}'


@cache
def cdk():
    """`jpype.JClass`, with the JVM up and CDK on its classpath -- or None.

    One JVM per process and no restart, so a JVM someone else started without the jar is a MISSING
    oracle, not a failure.  `chython/interop/_java.py` starts one from `config.class_paths`; importing
    it here would execute the facade, which `test_isolation.py` forbids, so this is a second launcher
    on purpose.

    `jpype.JClass` and not `jpype.imports`: `from java.io import StringReader` raises
    `ModuleNotFoundError` in this venv even after `import jpype.imports`.
    """
    if find_spec('jpype') is None or not CDK_JAR.is_file():
        return None
    import jpype
    if not jpype.isJVMStarted():
        jpype.startJVM('--enable-native-access=ALL-UNNAMED', classpath=[str(CDK_JAR)])
    try:
        jpype.JClass('org.openscience.cdk.silent.SilentChemObjectBuilder')
    except Exception:
        return None
    return jpype.JClass


@cache
def cdk_version():
    """`org.openscience.cdk.CDK.getVersion()`, or None."""
    JClass = cdk()
    return None if JClass is None else str(JClass('org.openscience.cdk.CDK').getVersion())


@cache
def marvin_version():
    """The banner line `molconvert` prints for `-h`, verbatim, or None."""
    if not MOLCONVERT.is_file():
        return None
    done = run([str(MOLCONVERT), '-h'], capture_output=True, text=True)
    lines = [x.strip() for x in done.stdout.splitlines() if x.strip()]
    return lines[0] if lines else None


_PROBES = {'rdkit': rdkit_version, 'indigo': indigo_version, 'cdk': cdk_version,
           'marvin': marvin_version}


def versions():
    """All four version strings, `None` for an absent oracle -- the report header."""
    return {name: probe() for name, probe in _PROBES.items()}


def require(name):
    """`pytest.skip` unless *name* is present.  A present oracle with an empty version is a failure."""
    try:
        probe = _PROBES[name]
    except KeyError:
        raise ValueError(f'unknown oracle {name!r}; have {sorted(_PROBES)}')
    version = probe()
    if version is None:
        pin = PINS.get(name)
        where = f'; set ${pin}' if pin else f'; pip install {name}'
        skip(f'{name} oracle absent{where}')
    if not version.strip():
        raise AssertionError(f'the {name} oracle is present and reports no version -- a provisioned '
                             f'oracle that cannot name itself makes every cell unattributable')


# ---------------------------------------------------------------------------------------------------
# The readback layer.  A file goes in, a parsed handle comes out -- one handle per (toolkit, format,
# text), because the feature probes below run against the handle rather than re-parsing per cell.  For
# Marvin, whose reader is a subprocess, that is the difference between one `molconvert` pair per fixture
# and one per matrix cell.
# ---------------------------------------------------------------------------------------------------

class Unrepresentable(Exception):
    """A toolkit cannot express this feature at all -- a declared capability gap, not a silent loss.

    Raised by a probe, caught by the classifier, reported as its own outcome.  The message is a
    statement of capability and nothing else.
    """
    def __init__(self, feature, toolkit, reason):
        super().__init__(f'{toolkit} cannot represent {feature}: {reason}')
        self.feature, self.toolkit, self.reason = feature, toolkit, reason


class Parsed(NamedTuple):
    """What a reference reader made of a file.

    `obj` is the toolkit's own handle, shaped for the probes: an `_RD` for RDKit, an `IndigoObject`, a
    `_CDK` carrying `JClass` alongside the container, and for Marvin the MRV text it re-emitted.
    `error` set means the reader refused the file, and that IS the cell's answer.  `evidence` is the
    literal call or argv that produced it, so a failure message can be re-run by hand.
    """
    toolkit: str
    obj: object
    skeleton: str
    error: str
    evidence: str


class _RD(NamedTuple):
    """RDKit hands out two molecules per file: what the file said, and a sanitized copy.

    `sanitize=False` first because RDKit's full sanitization refuses a molfile stating bond order 4 --
    which chython, Indigo 1.45 and `molconvert` all write for an unkekulized aromatic ring -- and the
    harness measures what the file said.  `clean` is `None` when the partial sanitize refuses; probes
    that need a perceived hydrogen count use `clean`, probes that read the file's own statement use
    `raw`.
    """
    raw: object
    clean: object


class _CDK(NamedTuple):
    """CDK's handle plus the `JClass` factory, because its constants are read off Java interfaces."""
    JClass: object
    obj: object


def skeleton_of(inchi):
    """Formula + `/c` + `/h` from an InChI, whatever prefix or extra layers it carries.

    Constitution only, on purpose: this is the "same molecule" half of a cell and the FEATURE is
    measured separately through each toolkit's own API.  The `/t`, `/m` and `/s` layers are not
    comparable across these four -- Indigo's `getInchi` is non-standard by default, emitting
    `InChI=1/...` with `/s3` where RDKit and CDK emit `InChI=1S/...` with `/m0/s1` -- so comparing them
    would report a difference about InChI rather than about the file.
    """
    if not inchi:
        return None
    body = inchi.split('/', 1)[1] if inchi.startswith('InChI=') and '/' in inchi else inchi
    parts = body.split('/')
    kept = [parts[0]] + [p for p in parts[1:] if p[:1] in 'ch']
    return '/'.join(kept)


# ------------------------------------------------------------------ chython, the library under test

def read_chython(text, fmt, log=None):
    """chython's own reader for *fmt*.  A list-returning reader is reduced to its single record."""
    log = [] if log is None else log
    if fmt in ('v2000', 'v3000'):
        return mol(text, log=log)
    if fmt == 'sdf':
        with SDFRead(StringIO(text)) as f:
            record, = f
            # the data fields belong to the frame and the core stores none, so the reader's view rides
            # alongside the container -- when a container grows a `meta` slot, `_ch_meta` reads that
            return _CH(record, dict(f.meta))
    if fmt == 'mrv':
        record, = read_mrv(text, log=log)
        return record
    if fmt == 'cml':
        record, = read_cml(text, log=log)
        return record
    if fmt == 'cxsmiles':
        return read_smiles(text)
    if fmt == 'inchi':
        return inchi_to_molecule(text)
    if fmt == 'rxn':
        return rxn(text, log=log)
    if fmt == 'rdf':
        with RDFRead(StringIO(text)) as f:
            record, = f
        return record
    if fmt == 'xyz':
        frame, = xyz(text, log=log)
        return frame
    if fmt == 'pdb':
        record, = read_pdb(text.splitlines(), log=log)
        return record
    if fmt == 'mmcif':
        record, = read_mmcif(text.splitlines(), log=log)
        return record
    if fmt == 'mol2':
        (record, _), = read_mol2(StringIO(text))
        return record
    raise KeyError(f'chython has no reader for {fmt!r}')


def _parse_chython(text, fmt):
    log = []
    call = f'chython.formats.test.oracles.read_chython(text, {fmt!r})'
    try:
        obj = read_chython(text, fmt, log=log)
    except Exception as e:
        return Parsed('chython', None, None, f'{type(e).__name__}: {e}', call)
    try:
        skeleton = skeleton_of(molecule_to_inchi(container_of(obj)))
    except Exception:
        skeleton = None
    return Parsed('chython', obj, skeleton, None, f'{call}  # log: {log}')


# ------------------------------------------------------------------------------------------- RDKit

#: The literal RDKit call per format, quoted verbatim in a failure message.
_RDKIT_READ = {
    'v2000': 'Chem.MolFromMolBlock(text, sanitize=False, removeHs=False)',
    'v3000': 'Chem.MolFromMolBlock(text, sanitize=False, removeHs=False)',
    'sdf': 'next(Chem.ForwardSDMolSupplier(BytesIO(text.encode()), sanitize=False, removeHs=False))',
    'cxsmiles': 'Chem.MolFromSmiles(text, sanitize=False)',
    'inchi': 'Chem.MolFromInchi(text, sanitize=False, removeHs=False)',
    'xyz': 'Chem.MolFromXYZBlock(text)',
    'pdb': 'Chem.MolFromPDBBlock(text, sanitize=False, removeHs=False)',
    'mol2': 'Chem.MolFromMol2Block(text, sanitize=False, removeHs=False)',
    'rxn': 'AllChem.ReactionFromRxnBlock(text, sanitize=False)',
}


def _parse_rdkit(text, fmt):
    from io import BytesIO
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog('rdApp.*')          # a parser diagnostic is data, not console noise
    call = _RDKIT_READ[fmt]
    try:
        obj = eval(call, {'Chem': Chem, 'AllChem': AllChem, 'BytesIO': BytesIO, 'text': text})
    except Exception as e:
        return Parsed('rdkit', None, None, f'{type(e).__name__}: {e}', call)
    if obj is None:
        return Parsed('rdkit', None, None, 'the reader returned None', call)
    if fmt == 'rxn':
        return Parsed('rdkit', obj, None, None, call)
    clean = Chem.Mol(obj)
    try:
        Chem.SanitizeMol(clean, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES ^ Chem.SANITIZE_KEKULIZE)
    except Exception:
        clean = None
    skeleton = None
    if clean is not None:
        try:
            skeleton = skeleton_of(Chem.MolToInchi(clean))
        except Exception:
            skeleton = None
    return Parsed('rdkit', _RD(obj, clean), skeleton, None, call)


# ------------------------------------------------------------------------------------------ Indigo

#: Indigo loads a molfile, a CML document and a SMILES through the same door.
_INDIGO_READ = {'v2000': 'loadMolecule', 'v3000': 'loadMolecule', 'sdf': 'loadMolecule',
                'cml': 'loadMolecule', 'cxsmiles': 'loadMolecule', 'inchi': 'inchi',
                'rxn': 'loadReaction', 'mol2': 'loadMolecule'}


def _parse_indigo(text, fmt):
    from indigo import Indigo
    from indigo.inchi import IndigoInchi
    session = Indigo()
    session.setOption('ignore-stereochemistry-errors', True)
    door = _INDIGO_READ[fmt]
    if door == 'inchi':
        call = 'IndigoInchi(Indigo()).loadMolecule(text)'
    else:
        call = f'Indigo().{door}(text)'
    try:
        if door == 'inchi':
            obj = IndigoInchi(session).loadMolecule(text)
        else:
            obj = getattr(session, door)(text)
    except Exception as e:
        return Parsed('indigo', None, None, f'{type(e).__name__}: {e}', call)
    try:
        skeleton = skeleton_of(IndigoInchi(session).getInchi(obj))
    except Exception:
        skeleton = None
    return Parsed('indigo', obj, skeleton, None, call)


# --------------------------------------------------------------------------------------------- CDK

_CDK_READER = {'v2000': 'org.openscience.cdk.io.MDLV2000Reader',
               'v3000': 'org.openscience.cdk.io.MDLV3000Reader',
               'sdf': 'org.openscience.cdk.io.MDLV2000Reader',
               'cml': 'org.openscience.cdk.io.CMLReader',
               'rxn': 'org.openscience.cdk.io.MDLRXNV3000Reader',
               'mol2': 'org.openscience.cdk.io.Mol2Reader',
               'pdb': 'org.openscience.cdk.io.PDBReader',
               'xyz': 'org.openscience.cdk.io.XYZReader'}

#: CDK 2.12 readers that read into a `ChemFile` and state `Only supported is reading of ChemFile
#: objects` for a container.  The container is taken back out with `ChemFileManipulator`.
_CDK_CHEMFILE = frozenset({'cml', 'xyz', 'pdb'})


def _parse_cdk(text, fmt):
    JClass = cdk()
    builder = JClass('org.openscience.cdk.silent.SilentChemObjectBuilder').getInstance()
    if fmt == 'inchi':
        call = 'InChIGeneratorFactory.getInstance().getInChIToStructure(text, SilentChemObjectBuilder)'
        try:
            factory = JClass('org.openscience.cdk.inchi.InChIGeneratorFactory').getInstance()
            obj = factory.getInChIToStructure(text, builder).getAtomContainer()
        except Exception as e:
            return Parsed('cdk', None, None, f'{type(e).__name__}: {e}', call)
    else:
        name = _CDK_READER[fmt]
        if fmt == 'rxn' and 'V3000' not in text.split('\n', 1)[0]:
            # CDK states one reader per RXN version and no sniffing wrapper, and its own MDLRXNWriter
            # writes V2000, so the version on the `$RXN` line picks the door
            name = 'org.openscience.cdk.io.MDLRXNReader'
        chemfile = fmt in _CDK_CHEMFILE
        # CDK 2.12's CMLReader states constructors over an InputStream and over a String filename and
        # none over a Reader, so that one door is opened with bytes
        source = 'ByteArrayInputStream(text.encode())' if fmt == 'cml' else 'StringReader(text)'
        seed = 'new ChemFile()' if chemfile else 'newReaction()' if fmt == 'rxn' \
            else 'newAtomContainer()'
        call = (f'{name.rsplit(".", 1)[1]}({source})'
                f'.read(SilentChemObjectBuilder.getInstance().{seed})')
        try:
            if fmt == 'cml':
                reader = JClass(name)(JClass('java.io.ByteArrayInputStream')(text.encode()))
            else:
                reader = JClass(name)(JClass('java.io.StringReader')(text))
            if chemfile:
                obj = reader.read(JClass('org.openscience.cdk.silent.ChemFile')())
                containers = JClass('org.openscience.cdk.tools.manipulator.ChemFileManipulator'
                                    ).getAllAtomContainers(obj)
                obj = containers.get(0) if containers.size() else None
                call = f'{call} -> ChemFileManipulator.getAllAtomContainers(...).get(0)'
            else:
                target = builder.newReaction() if fmt == 'rxn' else builder.newAtomContainer()
                obj = reader.read(target)
        except Exception as e:
            return Parsed('cdk', None, None, f'{type(e).__name__}: {e}', call)
        if obj is None:
            return Parsed('cdk', None, None, 'the reader returned no container', call)
    if fmt == 'rxn':
        return Parsed('cdk', _CDK(JClass, obj), None, None, call)
    try:
        gen = JClass('org.openscience.cdk.inchi.InChIGeneratorFactory').getInstance()
        skeleton = skeleton_of(str(gen.getInChIGenerator(obj).getInchi()))
    except Exception:
        skeleton = None
    return Parsed('cdk', _CDK(JClass, obj), skeleton, None, call)


# ------------------------------------------------------------------------------------------- Marvin

def molconvert(text, target, extra=()):
    """`molconvert` on stdin.  Returns `(returncode, stdout, stderr, argv)` -- argv IS the evidence."""
    argv = [str(MOLCONVERT), '-g', *extra, target]
    done = run(argv, input=text, capture_output=True, text=True)
    return done.returncode, done.stdout, done.stderr, ' '.join(argv) + '  < (stdin)'


def _parse_marvin(text, fmt):
    """Marvin's readback is its own MRV, and the skeleton comes from the molecule that MRV describes.

    `molconvert -g inchi` reports `InChI native library is not available for Mac.` and yields
    `InChI=1S//` on this platform, so the skeleton cannot come from Marvin's own InChI.  chython
    supplies the InChI *for the molecule Marvin reported* -- the cell still measures Marvin's reading,
    not chython's, because the constitution being compared is the one Marvin handed back.  When the MRV
    is not something chython reads (a reaction, an ill-formed document) the fallback is Marvin's
    `cxsmiles`, which it also produces from the same input.

    The `smiles` target refuses an atom carrying a stereo group, naming `cxsmiles`, `smarts`, `cxsmarts`
    and `mrv` as the targets that carry it, so `cxsmiles` and not `smiles` is the fallback here.
    """
    rc, mrv, err, argv = molconvert(text, 'mrv')
    if rc or not mrv.strip():
        return Parsed('marvin', None, None, (err or mrv).strip()[:400], argv)
    return _marvin_from_mrv(text, mrv, argv)


def _marvin_from_mrv(text, mrv, argv):
    """One Marvin answer, given the MRV it re-emitted for `text`."""
    skeleton = None
    try:
        record, = read_mrv(mrv)
        skeleton = skeleton_of(molecule_to_inchi(record))
    except Exception:
        rc, cx, _, argv2 = molconvert(text, 'cxsmiles')
        if not rc and cx.strip():
            argv = f'{argv} ; {argv2}'
            try:
                skeleton = skeleton_of(molecule_to_inchi(read_smiles(cx.strip().split('\t')[0])))
            except Exception:
                skeleton = None
    return Parsed('marvin', mrv, skeleton, None, argv)


#: How several records of one format are handed to `molconvert` in a single call.  A `molconvert` run
#: costs about 0.65 s of JVM start, and the matrix asks for about a hundred Marvin readbacks; batching
#: the formats that concatenate cleanly turns those into about eight calls.  A format not listed here is
#: read one record per call.
_MARVIN_BATCHABLE = {
    'v2000': lambda texts: ''.join(f'{t.rstrip()}\n$$$$\n' for t in texts),
    'v3000': lambda texts: ''.join(f'{t.rstrip()}\n$$$$\n' for t in texts),
    'sdf': lambda texts: ''.join(f'{t.rstrip().removesuffix("$$$$").rstrip()}\n$$$$\n' for t in texts),
    'cxsmiles': lambda texts: ''.join(f'{t.strip()}\n' for t in texts),
}


def _split_mrv(document):
    """One single-record MRV document per `<MDocument>` in a multi-record one, header kept."""
    from re import findall
    head = document.find('<cml')
    if head < 0:
        return []
    prefix = document[:document.index('>', head) + 1]
    return [f'{prefix}{found}</cml>' for found in findall(r'<MDocument>.*?</MDocument>', document, 16)]


def prewarm_marvin(parses, written):
    """Fill the Marvin readback cache in batches, so the matrix pays one JVM start per format.

    A batch whose record count does not come back intact is dropped rather than aligned by guesswork --
    `parse_once` then reads those records one call at a time.  The evidence string states the
    single-record command that reproduces the cell by hand.
    """
    groups = {}
    for (_, fmt, _), text in written.items():
        if isinstance(text, BaseException) or fmt not in _MARVIN_BATCHABLE:
            continue
        if ('marvin', fmt, hash(text)) not in parses:
            groups.setdefault(fmt, []).append(text)
    for fmt, texts in groups.items():
        texts = list(dict.fromkeys(texts))
        rc, out, _, argv = molconvert(_MARVIN_BATCHABLE[fmt](texts), 'mrv')
        if rc:
            continue
        documents = _split_mrv(out)
        if len(documents) != len(texts):
            continue
        for text, mrv in zip(texts, documents):
            evidence = f'{argv} (one of {len(texts)} {fmt} records in one batch)'
            parses['marvin', fmt, hash(text)] = _marvin_from_mrv(text, mrv, evidence)


#: How a multi-record `molconvert` output of one target is cut back into one record per seed.  A target
#: not listed here is written one seed per call: Marvin's CML target puts several `<molecule>` elements in
#: one document, and cutting that into single-molecule documents is a rewrite rather than a split.
_MARVIN_SPLIT = {
    'v2000': lambda out: _split_after(out, 'M  END'),
    'v3000': lambda out: _split_after(out, 'M  END'),
    'sdf': lambda out: _split_after(out, '$$$$'),
    'mrv': lambda out: _split_mrv(out),
    'cxsmiles': lambda out: [line for line in out.splitlines() if line.strip()],
}


def _split_after(document, terminator):
    """One record per `terminator` line, the line kept with the record it ends."""
    out, current = [], []
    for line in document.splitlines(keepends=True):
        current.append(line)
        if line.rstrip() == terminator:
            out.append(''.join(current))
            current = []
    return out


def prewarm_marvin_writes(written, seeds):
    """Write the Marvin fixtures for several seeds per call, one call per target.

    Same trade as the readback batch: `molconvert` costs a JVM start, and its targets that concatenate
    cleanly let one start serve every seed.  A short or long batch is dropped, and `write_once` then
    writes those seeds one at a time.
    """
    for fmt, split in _MARVIN_SPLIT.items():
        batch = [s for s in seeds if not is_reaction(s) and fmt in seed_formats(s)
                 and ('marvin', fmt, s.key) not in written
                 # a fixture with its own framed input is written one at a time
                 and _marvin_source(s, fmt) is None]
        if len(batch) < 2:
            continue
        code, out, _, _ = molconvert(''.join(f'{s.cxsmiles}\n' for s in batch),
                                     _MARVIN_TARGET[fmt], ('-2',))
        if code:
            continue
        records = split(out)
        if len(records) != len(batch):
            continue
        for seed, record in zip(batch, records):
            written['marvin', fmt, seed.key] = record


def warm(written, parses):
    """Write every fixture the matrix can ask for, then batch the Marvin readbacks over them.

    Eager on purpose: batching needs the whole set of files up front, and a full run writes them all
    anyway.  A writer that refuses is remembered as its exception, so this pass never raises.

    The Marvin batches are skipped when `molconvert` is absent -- this runs from a session fixture, and
    an exception here would ERROR every cell on a machine that simply does not have the oracle.
    """
    batched = marvin_version() is not None
    if batched:
        prewarm_marvin_writes(written, SEEDS)
    for seed in SEEDS:
        formats = seed_formats(seed) & READS['chython']
        for toolkit, fmt in WRITERS:
            if fmt in formats:
                try:
                    write_once(written, toolkit, fmt, seed)
                except Exception:
                    pass
    if batched:
        prewarm_marvin(parses, written)


#: One reader per toolkit.  `chython` is here too, so the classifier treats the five uniformly.
READERS = {'chython': _parse_chython, 'rdkit': _parse_rdkit, 'indigo': _parse_indigo,
           'cdk': _parse_cdk, 'marvin': _parse_marvin}

#: Which formats each reader has a door for.  A cell for a missing pair does not exist.
READS = {
    'chython': frozenset({'v2000', 'v3000', 'sdf', 'mrv', 'cml', 'cxsmiles', 'inchi', 'rxn', 'rdf',
                          'xyz', 'pdb', 'mmcif', 'mol2'}),
    'rdkit': frozenset(_RDKIT_READ),
    'indigo': frozenset(_INDIGO_READ),
    'cdk': frozenset(_CDK_READER) | {'inchi'},
    # `molconvert` sniffs its input format, so every format the harness produces goes in the same door.
    'marvin': frozenset({'v2000', 'v3000', 'sdf', 'mrv', 'cml', 'cxsmiles', 'rxn', 'rdf', 'xyz',
                         'pdb', 'mol2'}),
}


# ---------------------------------------------------------------------------------------------------
# The feature probes.  Five dicts with the same keys, one entry per (toolkit, feature), each returning a
# value comparable ACROSS toolkits.  Atom indices are forbidden in a probe value: on the same file RDKit
# reports an AND group as `&1:1,5` and Indigo as `&1:1,4` because their atom orders differ, so a probe
# reduces to shape -- sorted member counts, sorted (kind, count) pairs, a frozenset of (name, value)
# tuples, a multiset of charges.
# ---------------------------------------------------------------------------------------------------

#: Every feature the matrix knows.  All five probe dicts carry exactly these keys; `test_oracles.py`
#: asserts it, because a probe missing from one dict silently drops a whole toolkit's column.
FEATURE_NAMES = ('stereo_groups', 'dat_sgroup', 'other_sgroups', 'wedge', 'either_bond', 'charge',
                 'radical', 'isotope', 'dative', 'aromatic_input', 'implicit_h', 'alias',
                 'fragment_group', 'sd_fields', 'reaction_sides', 'mapping', 'element_multiset',
                 'coordinate_count')

#: Where a feature can be stated at all.  A CTfile wedge has no CXSMILES spelling and an SD field has no
#: molfile spelling, so those cells do not exist and are not findings.  This is what keeps the matrix
#: from asking `molconvert` a question no format can carry.
FEATURE_FORMATS = {
    'stereo_groups': frozenset({'v2000', 'v3000', 'sdf', 'mrv', 'cxsmiles'}),
    'dat_sgroup': frozenset({'v2000', 'v3000', 'sdf', 'mrv', 'cxsmiles'}),
    'other_sgroups': frozenset({'v2000', 'v3000', 'sdf', 'mrv'}),
    'wedge': frozenset({'v2000', 'v3000', 'sdf', 'mrv'}),
    'either_bond': frozenset({'v2000', 'v3000', 'sdf', 'mrv'}),
    'charge': frozenset({'v2000', 'v3000', 'sdf', 'mrv', 'cml', 'cxsmiles'}),
    'radical': frozenset({'v2000', 'v3000', 'sdf', 'mrv', 'cxsmiles'}),
    'isotope': frozenset({'v2000', 'v3000', 'sdf', 'mrv', 'cml', 'cxsmiles'}),
    'dative': frozenset({'v2000', 'v3000', 'sdf', 'mrv', 'cxsmiles'}),
    'aromatic_input': frozenset({'v2000', 'v3000', 'sdf'}),
    'implicit_h': frozenset({'v2000', 'v3000', 'sdf', 'mrv', 'cml'}),
    'alias': frozenset({'v2000', 'v3000', 'sdf', 'mrv'}),
    'fragment_group': frozenset({'v2000', 'v3000', 'sdf', 'cxsmiles'}),
    # MRV and CML carry the same fields as a `<propertyList>` inside the molecule, which is where
    # `molconvert` puts an SD field converting an SDF.
    'sd_fields': frozenset({'sdf', 'mrv', 'cml'}),
    'reaction_sides': frozenset({'rxn', 'rdf'}),
    'mapping': frozenset({'rxn', 'rdf'}),
    'element_multiset': frozenset({'xyz', 'pdb', 'mmcif', 'mol2'}),
    'coordinate_count': frozenset({'xyz', 'pdb', 'mmcif', 'mol2'}),
}

#: chython's wedge codes, and the names every probe reduces to.  Each toolkit spells the three states
#: with its own constants -- V2000 column `1/4/6`, V3000 `CFG=1/2/3`, `IBond.Stereo`, `<bondStereo>` --
#: so the comparable value is the state and not the number.
_WEDGE_NAMES = {1: 'up', 2: 'down', 3: 'either'}


def _mapping_value(left, right):
    """`(numbers, the two sides agree)`, or `()` for a record that states no mapping at all.

    The empty answer has to be falsy: the classifier reads a falsy probe as "the writer stated no such
    feature", and `((), True)` would make an unmapped record look like a recovered cell.
    """
    numbers = left | right
    return (tuple(sorted(numbers)), left == right) if numbers else ()


def _by_text(pair):
    """Sort key for a probe whose value mixes a count with a word -- `unknown` beside a 3."""
    return pair[0], str(pair[1])


#: CTfile Sgroup type per toolkit spelling, for the non-DAT types.  A type absent from a table keeps the
#: toolkit's own spelling, so an unmapped type shows up in the cell rather than being silently equated.
_CDK_SGROUP_TYPES = {'CtabAbbreviation': 'SUP', 'CtabMultipleGroup': 'MUL',
                     'CtabStructureRepeatUnit': 'SRU', 'CtabMonomer': 'MON', 'CtabCopolymer': 'COP',
                     'CtabCrossLink': 'CRO', 'CtabModified': 'MOD', 'CtabGraft': 'GRA',
                     'CtabComponent': 'COM', 'CtabMer': 'MER', 'CtabFormulation': 'FOR',
                     'CtabMixture': 'MIX', 'CtabAnyPolymer': 'ANY', 'CtabGeneric': 'GEN',
                     'CtabData': 'DAT'}
_MRV_SGROUP_ROLES = {'DataSgroup': 'DAT', 'SuperatomSgroup': 'SUP', 'MultipleSgroup': 'MUL',
                     'SruSgroup': 'SRU', 'MonomerSgroup': 'MON', 'CopolymerSgroup': 'COP',
                     'GenericSgroup': 'GEN', 'ComponentSgroup': 'COM', 'MixtureSgroup': 'MIX',
                     'FormulationSgroup': 'FOR', 'AnyPolymerSgroup': 'ANY'}


# ------------------------------------------------------------------------- Marvin's MRV, once, for all
# Marvin's readback is an MRV document and every Marvin probe reads it, so the document is normalised
# once here: namespaces stripped, and the two forms Marvin writes -- array attributes on `<atomArray>`
# and `<atom>` child elements -- reduced to one list of dicts.

def _mrv_tree(mrv):
    """The MRV document with namespaces stripped, so `find('atomArray')` works."""
    from xml.etree.ElementTree import fromstring
    root = fromstring(mrv)
    for el in root.iter():
        if '}' in el.tag:
            el.tag = el.tag.rsplit('}', 1)[1]
    return root


def _mrv_molecules(mrv):
    """Top-level `<molecule>` elements in document order; an S-group's nested molecule is not one.

    For a `<reaction>` the sides are flattened, since `_mrv_sides` reads the counts separately.
    """
    out = []
    for struct in _mrv_tree(mrv).iter('MChemicalStruct'):
        for child in struct:
            if child.tag == 'molecule':
                out.append(child)
            elif child.tag == 'reaction':
                for side in child:
                    out.extend(m for m in side if m.tag == 'molecule')
    return out


def _mrv_sides(mrv):
    """`{'reactantList': [<molecule>, ...], ...}` for a `<reaction>`, empty for a molecule document."""
    out = {}
    for reaction in _mrv_tree(mrv).iter('reaction'):
        for side in reaction:
            out.setdefault(side.tag, []).extend(m for m in side if m.tag == 'molecule')
    return out


def _mrv_atoms(el):
    """One dict per atom, whichever of MRV's two spellings the document used.

    Array form puts space-separated values on `<atomArray>` aligned with `atomID`; element form puts
    them on `<atom>` children.  `atomID` becomes `id` so a probe reads one key either way.
    """
    array = el.find('atomArray')
    if array is None:
        return []
    kids = [a for a in array if a.tag == 'atom']
    if kids:
        return [dict(a.attrib) for a in kids]
    columns = {('id' if k == 'atomID' else k): v.split() for k, v in array.attrib.items()}
    count = len(columns.get('id', ()))
    return [{k: v[i] for k, v in columns.items() if i < len(v)} for i in range(count)]


def _mrv_bonds(el):
    """One dict per bond, with `<bondStereo>`'s text and attributes folded in under `bondStereo`."""
    out = []
    array = el.find('bondArray')
    if array is None:
        return out
    for bond in array:
        if bond.tag != 'bond':
            continue
        row = dict(bond.attrib)
        stereo = bond.find('bondStereo')
        if stereo is not None:
            row['bondStereo'] = (stereo.text or '').strip()
            row['bondStereoValue'] = stereo.get('conventionValue', '')
        out.append(row)
    return out


def _mrv_sgroups(mrv):
    """`(role, atom count)` per nested S-group molecule, plus its attributes."""
    out = []
    for molecule in _mrv_molecules(mrv):
        for nested in molecule.iter('molecule'):
            role = nested.get('role')
            if role:
                out.append((role, dict(nested.attrib)))
    return out


def _mrv_order(text):
    """An MRV bond order as a number, `A` being CTfile's 4."""
    return 4 if text == 'A' else int(text or 1)


def _mrv_wedge(row):
    """`'up' | 'down' | 'either' | None` from one `_mrv_bonds` row.

    `W` and `H` are Marvin's wedge and hash; an either bond comes back as
    `<bondStereo convention="MDL" conventionValue="4">` with no text.
    """
    text = row.get('bondStereo', '')
    if text in ('W', 'w'):
        return 'up'
    if text in ('H', 'h'):
        return 'down'
    if row.get('bondStereoValue') == '4' or text in ('E', 'e'):
        return 'either'
    return None


# ------------------------------------------------------------------------------------ chython's probes

class _CH(NamedTuple):
    """chython's record for a framed format: the container plus the data fields the frame carried.

    `MoleculeContainer` has no `meta` slot on this base, so a reader's field view rides alongside.  When
    the container grows one, `_ch_meta` reads it and this wrapper collapses to the container.
    """
    obj: object
    meta: dict


def container_of(obj):
    """The container, whether it arrived bare or inside a `_CH`."""
    return obj.obj if isinstance(obj, _CH) else obj


#: The chython reader's own field view rides in `_CH`; `_ch` is the short spelling used by the probes.
_ch = container_of


def _ch_meta(obj):
    if isinstance(obj, _CH):
        return dict(obj.meta)
    return dict(getattr(obj, 'meta', None) or {})


def _ch_stereo_groups(obj):
    return tuple(sorted((kind, len(members))
                        for (kind, _), members in _ch(obj).canonical_stereo_groups().items()
                        if kind in (2, 3)))


def _ch_dat(obj):
    return frozenset((r['name'].decode(), b'\n'.join(r['data']).decode(), len(r['atoms']))
                     for r in _ch(obj).sgroups if r['type'] == b'DAT')


def _ch_other_sgroups(obj):
    return frozenset((r['type'].decode(), len(r['atoms']))
                     for r in _ch(obj).sgroups if r['type'] != b'DAT')


def _ch_wedge(obj):
    return tuple(sorted((b.order, _WEDGE_NAMES[b.wedge[1]])
                        for b in _ch(obj).bonds() if b.wedge and b.wedge[1] in _WEDGE_NAMES))


def _ch_either(obj):
    return sum(1 for b in _ch(obj).bonds() if b.wedge and b.wedge[1] == 3)


def _ch_charge(obj):
    return tuple(sorted(a.charge for a in _ch(obj).atoms() if a.charge))


def _ch_radical(obj):
    return sum(1 for a in _ch(obj).atoms() if a.is_radical)


def _ch_isotope(obj):
    return tuple(sorted((a.atomic_symbol, a.isotope) for a in _ch(obj).atoms() if a.isotope))


def _ch_dative(obj):
    return sum(1 for b in _ch(obj).bonds() if b.order == 8)


def _ch_aromatic(obj):
    return tuple(sorted(b.order for b in _ch(obj).bonds()))


def _ch_implicit_h(obj):
    """`(symbol, total H)` per heavy atom.  `H_UNKNOWN` becomes the string, never a number.

    TOTAL and not implicit: RDKit and CDK distribute a stated hydrogen count between an implicit and an
    explicit slot by their own rules, so the count comparable across five toolkits is the sum.
    """
    return tuple(sorted(((a.atomic_symbol, 'unknown' if a.total_h == 15 else a.total_h)
                         for a in _ch(obj).atoms()), key=_by_text))


def _ch_alias(obj):
    return frozenset(v.decode() for v in _ch(obj).aliases.values())


def _ch_fragment_group(obj):
    mol_ = _ch(obj)
    return tuple(sorted(len(c) for c in mol_.connected_components))


def _ch_sd_fields(obj):
    return frozenset(_ch_meta(obj).items())


def _ch_sides(obj):
    r = _ch(obj)
    return (len(r.reactants), len(r.agents), len(r.products))


def _ch_mapping(obj):
    r = _ch(obj)
    left = {a.map_number for m in r.reactants for a in m.atoms() if a.map_number}
    right = {a.map_number for m in r.products for a in m.atoms() if a.map_number}
    return _mapping_value(left, right)


def _ch_elements(obj):
    record = _ch(obj)
    if hasattr(record, 'atoms') and not callable(record.atoms):
        return tuple(sorted(a.element for a in record.atoms))
    return tuple(sorted(a.atomic_symbol for a in record.atoms()))


def _ch_coordinates(obj):
    record = _ch(obj)
    if hasattr(record, 'atoms') and not callable(record.atoms):
        return sum(1 for a in record.atoms if a.x is not None and a.y is not None)
    return sum(1 for a in record.atoms() if a.x is not None and a.y is not None)


CHYTHON_PROBES = {
    'stereo_groups': _ch_stereo_groups, 'dat_sgroup': _ch_dat, 'other_sgroups': _ch_other_sgroups,
    'wedge': _ch_wedge, 'either_bond': _ch_either, 'charge': _ch_charge, 'radical': _ch_radical,
    'isotope': _ch_isotope, 'dative': _ch_dative, 'aromatic_input': _ch_aromatic,
    'implicit_h': _ch_implicit_h, 'alias': _ch_alias, 'fragment_group': _ch_fragment_group,
    'sd_fields': _ch_sd_fields, 'reaction_sides': _ch_sides, 'mapping': _ch_mapping,
    'element_multiset': _ch_elements, 'coordinate_count': _ch_coordinates,
}


# -------------------------------------------------------------------------------------- RDKit's probes

def _rd_stereo_groups(h):
    from rdkit.Chem import StereoGroupType
    kinds = {StereoGroupType.STEREO_OR: 2, StereoGroupType.STEREO_AND: 3}
    return tuple(sorted((kinds[g.GetGroupType()], len(g.GetAtoms()))
                        for g in h.raw.GetStereoGroups() if g.GetGroupType() in kinds))


def _rd_substance_groups(h, want_dat):
    from rdkit import Chem
    out = set()
    for group in Chem.GetMolSubstanceGroups(h.raw):
        props = group.GetPropsAsDict()
        kind = props.get('TYPE', '')
        if (kind == 'DAT') != want_dat:
            continue
        atoms = len(list(group.GetAtoms()))
        if want_dat:
            try:
                data = '\n'.join(group.GetStringVectProp('DATAFIELDS'))
            except Exception:
                data = ''
            out.add((props.get('FIELDNAME', ''), data, atoms))
        else:
            out.add((kind, atoms))
    return frozenset(out)


#: V2000 writes the wedge in the bond block's stereo column, V3000 in `CFG=`, and the two numberings
#: differ; RDKit keeps whichever the file used under its own property name.
_RD_V2000_WEDGE = {1: 'up', 4: 'either', 6: 'down'}
_RD_V3000_WEDGE = {1: 'up', 2: 'either', 3: 'down'}


def _rd_wedges(h):
    out = []
    for bond in h.raw.GetBonds():
        props = bond.GetPropsAsDict()
        order = props.get('_MolFileBondType', int(bond.GetBondTypeAsDouble()))
        if '_MolFileBondStereo' in props:
            name = _RD_V2000_WEDGE.get(props['_MolFileBondStereo'])
        elif '_MolFileBondCfg' in props:
            name = _RD_V3000_WEDGE.get(props['_MolFileBondCfg'])
        else:
            name = None
        if name is not None:
            out.append((order, name))
    return tuple(sorted(out))


def _rd_isotope(h):
    return tuple(sorted((a.GetSymbol(), a.GetIsotope()) for a in h.raw.GetAtoms() if a.GetIsotope()))


def _rd_dative(h):
    from rdkit.Chem import BondType
    return sum(1 for b in h.raw.GetBonds() if b.GetBondType() == BondType.DATIVE)


def _rd_aromatic(h):
    from rdkit.Chem import BondType
    out = []
    for bond in h.raw.GetBonds():
        props = bond.GetPropsAsDict()
        if '_MolFileBondType' in props:
            out.append(props['_MolFileBondType'])
        elif bond.GetBondType() == BondType.AROMATIC:
            out.append(4)
        else:
            out.append(int(bond.GetBondTypeAsDouble()))
    return tuple(sorted(out))


def _rd_implicit_h(h):
    """RDKit perceives a hydrogen count only after a sanitize, so this reads the sanitized copy."""
    if h.clean is None:
        raise Unrepresentable('implicit_h', 'rdkit',
                              'a hydrogen count is perceived during sanitization, and the partial '
                              'sanitize this file admits did not complete')
    return tuple(sorted((a.GetSymbol(), a.GetTotalNumHs()) for a in h.clean.GetAtoms()))


def _rd_alias(h):
    return frozenset(a.GetPropsAsDict()['molFileAlias']
                     for a in h.raw.GetAtoms() if 'molFileAlias' in a.GetPropsAsDict())


def _rd_fragment_group(h):
    from rdkit import Chem
    return tuple(sorted(len(f) for f in Chem.GetMolFrags(h.raw)))


def _rd_sd_fields(h):
    return frozenset((k, v) for k, v in h.raw.GetPropsAsDict(includePrivate=False).items()
                     if isinstance(v, str) and not k.startswith('_'))


def _rd_sides(r):
    return (r.GetNumReactantTemplates(), r.GetNumAgentTemplates(), r.GetNumProductTemplates())


def _rd_mapping(r):
    left = {a.GetAtomMapNum() for m in r.GetReactants() for a in m.GetAtoms() if a.GetAtomMapNum()}
    right = {a.GetAtomMapNum() for m in r.GetProducts() for a in m.GetAtoms() if a.GetAtomMapNum()}
    return _mapping_value(left, right)


def _rd_coordinates(h):
    conf = h.raw.GetNumConformers() and h.raw.GetConformer()
    return h.raw.GetNumAtoms() if conf else 0


RDKIT_PROBES = {
    'stereo_groups': _rd_stereo_groups,
    'dat_sgroup': lambda h: _rd_substance_groups(h, True),
    'other_sgroups': lambda h: _rd_substance_groups(h, False),
    'wedge': _rd_wedges,
    'either_bond': lambda h: sum(1 for _, n in _rd_wedges(h) if n == 'either'),
    'charge': lambda h: tuple(sorted(a.GetFormalCharge() for a in h.raw.GetAtoms()
                                     if a.GetFormalCharge())),
    'radical': lambda h: sum(1 for a in h.raw.GetAtoms() if a.GetNumRadicalElectrons()),
    'isotope': _rd_isotope, 'dative': _rd_dative, 'aromatic_input': _rd_aromatic,
    'implicit_h': _rd_implicit_h, 'alias': _rd_alias, 'fragment_group': _rd_fragment_group,
    'sd_fields': _rd_sd_fields, 'reaction_sides': _rd_sides, 'mapping': _rd_mapping,
    'element_multiset': lambda h: tuple(sorted(a.GetSymbol() for a in h.raw.GetAtoms())),
    'coordinate_count': _rd_coordinates,
}


# ------------------------------------------------------------------------------------- Indigo's probes

def _in_stereo_groups(m):
    """Indigo states the collection per atom: `stereocenterType()` 1 ABS / 2 OR / 3 AND, plus
    `stereocenterGroup()`.  A non-stereocentre raises, which is the "not in a collection" answer."""
    seen = {}
    for atom in m.iterateAtoms():
        try:
            kind, group = atom.stereocenterType(), atom.stereocenterGroup()
        except Exception:
            continue
        if kind in (2, 3):
            seen[(kind, group)] = seen.get((kind, group), 0) + 1
    return tuple(sorted((kind, n) for (kind, _), n in seen.items()))


def _in_dat(m):
    return frozenset((d.description(), d.data(), len(list(d.iterateAtoms())))
                     for d in m.iterateDataSGroups())


def _in_other_sgroups(m):
    out = set()
    for door, kind in (('iterateSuperatoms', 'SUP'), ('iterateRepeatingUnits', 'SRU'),
                       ('iterateMultipleGroups', 'MUL'), ('iterateGenericSGroups', 'GEN')):
        try:
            groups = list(getattr(m, door)())
        except Exception:
            continue
        for group in groups:
            out.add((kind, len(list(group.iterateAtoms()))))
    return frozenset(out)


#: Indigo's bond stereo constants: 4 either, 5 up, 6 down.
_IN_WEDGE = {4: 'either', 5: 'up', 6: 'down'}


def _in_wedges(m):
    return tuple(sorted((b.bondOrder(), _IN_WEDGE[b.bondStereo()])
                        for b in m.iterateBonds() if b.bondStereo() in _IN_WEDGE))


def _in_alias(m):
    raise Unrepresentable('alias', 'indigo',
                          'the Indigo 1.45 Python API states no atom-alias accessor, so an `A` line '
                          'cannot be read back through it')


def _in_sd_fields(m):
    return frozenset((p.name(), p.rawData()) for p in m.iterateProperties())


def _in_fragment_group(m):
    return tuple(sorted(c.clone().countAtoms() for c in m.iterateComponents()))


def _in_sides(r):
    return (r.countReactants(), r.countCatalysts(), r.countProducts())


def _in_mapping(r):
    left = {r.atomMappingNumber(a) for m in r.iterateReactants() for a in m.iterateAtoms()
            if r.atomMappingNumber(a)}
    right = {r.atomMappingNumber(a) for m in r.iterateProducts() for a in m.iterateAtoms()
             if r.atomMappingNumber(a)}
    return _mapping_value(left, right)


INDIGO_PROBES = {
    'stereo_groups': _in_stereo_groups, 'dat_sgroup': _in_dat, 'other_sgroups': _in_other_sgroups,
    'wedge': _in_wedges,
    'either_bond': lambda m: sum(1 for b in m.iterateBonds() if b.bondStereo() == 4),
    'charge': lambda m: tuple(sorted(a.charge() for a in m.iterateAtoms() if a.charge())),
    'radical': lambda m: sum(1 for a in m.iterateAtoms() if a.radical()),
    'isotope': lambda m: tuple(sorted((a.symbol(), a.isotope())
                                      for a in m.iterateAtoms() if a.isotope())),
    'dative': lambda m: sum(1 for b in m.iterateBonds() if b.bondOrder() == 8),
    'aromatic_input': lambda m: tuple(sorted(b.bondOrder() for b in m.iterateBonds())),
    'implicit_h': lambda m: tuple(sorted((a.symbol(), a.countHydrogens())
                                         for a in m.iterateAtoms())),
    'alias': _in_alias, 'fragment_group': _in_fragment_group, 'sd_fields': _in_sd_fields,
    'reaction_sides': _in_sides, 'mapping': _in_mapping,
    'element_multiset': lambda m: tuple(sorted(a.symbol() for a in m.iterateAtoms())),
    'coordinate_count': lambda m: sum(1 for a in m.iterateAtoms() if any(a.xyz())),
}


# ---------------------------------------------------------------------------------------- CDK's probes

def _cdk_stereo_groups(h):
    """CDK packs the collection into `getGroupInfo()` on the stereo element.

    `GRP_ABS` is zero, so an unmarked stereocentre and an explicitly absolute one are the same value --
    which is why the comparable set is the OR and AND collections and not the absolute one, in all five
    probes.  The masks come from `IStereoElement`, not spelled as integers here.
    """
    SE = h.JClass('org.openscience.cdk.interfaces.IStereoElement')
    kinds = {int(SE.GRP_REL): 2, int(SE.GRP_RAC): 3}
    type_mask, num_mask = int(SE.GRP_TYPE_MASK), int(SE.GRP_NUM_MASK)
    seen = {}
    for element in h.obj.stereoElements():
        info = int(element.getGroupInfo())
        kind = kinds.get(info & type_mask)
        if kind is None:
            continue
        key = (kind, info & num_mask)
        seen[key] = seen.get(key, 0) + 1
    return tuple(sorted((kind, n) for (kind, _), n in seen.items()))


def _cdk_sgroups(h, want_dat):
    groups = h.obj.getProperty('cdk:CtabSgroups')
    if groups is None:
        return frozenset()
    key = h.JClass('org.openscience.cdk.sgroup.SgroupKey')
    out = set()
    for group in groups:
        kind = _CDK_SGROUP_TYPES.get(str(group.getType()), str(group.getType()))
        if (kind == 'DAT') != want_dat:
            continue
        atoms = group.getAtoms().size()
        if want_dat:
            name = group.getValue(key.DataFieldName)
            data = group.getValue(key.Data)
            out.add(('' if name is None else str(name), '' if data is None else str(data), atoms))
        else:
            out.add((kind, atoms))
    return frozenset(out)


#: `IBond.Stereo` names, reduced to the three states.  `E_Z_BY_COORDINATES` is not a wedge.
_CDK_WEDGE = {'UP': 'up', 'DOWN': 'down', 'UP_OR_DOWN': 'either',
              'UP_INVERTED': 'up', 'DOWN_INVERTED': 'down', 'UP_OR_DOWN_INVERTED': 'either'}


def _cdk_order(bond):
    """A CDK bond order as a CTfile number.  `UNSET` on an aromatic-flagged bond is CTfile's 4."""
    name = str(bond.getOrder())
    if name == 'UNSET' or bond.isAromatic():
        return 4
    return {'SINGLE': 1, 'DOUBLE': 2, 'TRIPLE': 3, 'QUADRUPLE': 4, 'QUINTUPLE': 5,
            'SEXTUPLE': 6}.get(name, name)


def _cdk_wedges(h):
    return tuple(sorted((_cdk_order(b), _CDK_WEDGE[str(b.getStereo())])
                        for b in h.obj.bonds() if str(b.getStereo()) in _CDK_WEDGE))


def _cdk_implicit_h(h):
    """CDK 2.12's readers leave `getImplicitHydrogenCount()` null and fill it in a separate
    configuration step, so this asks for the count the way CDK is driven to produce one.

    Two calls, in this order: the adder alone raises `CDKException: IAtom is not typed!`.  Null after
    both is CDK stating no count for that atom, which stays `'unstated'` -- reading it as 0 would turn
    "no answer" into an answer and hide exactly the cells this matrix is for.
    """
    try:
        h.JClass('org.openscience.cdk.tools.manipulator.AtomContainerManipulator') \
            .percieveAtomTypesAndConfigureAtoms(h.obj)
        h.JClass('org.openscience.cdk.tools.CDKHydrogenAdder') \
            .getInstance(h.obj.getBuilder()).addImplicitHydrogens(h.obj)
    except Exception:
        pass  # an atom CDK does not type: whatever the read left stands, null included
    return tuple(sorted(((str(a.getSymbol()),
                          'unstated' if a.getImplicitHydrogenCount() is None
                          else int(a.getImplicitHydrogenCount()))
                         for a in h.obj.atoms()), key=_by_text))


def _cdk_alias(h):
    out = set()
    for atom in h.obj.atoms():
        label = getattr(atom, 'getLabel', None)
        if label is None:
            continue
        try:
            text = label()
        except Exception:
            continue
        if text:
            out.add(str(text))
    return frozenset(out)


def _cdk_fragment_group(h):
    partition = h.JClass('org.openscience.cdk.graph.ConnectivityChecker')
    return tuple(sorted(int(c.getAtomCount())
                        for c in partition.partitionIntoMolecules(h.obj).atomContainers()))


def _cdk_dative(h):
    """CDK 2.12's `IBond.Order` enumerates SINGLE..SEXTUPLE and no dative member."""
    raise Unrepresentable('dative', 'cdk',
                          'IBond.Order in CDK 2.12 enumerates SINGLE through SEXTUPLE, so a '
                          'coordination bond has no order to be read into')


def _cdk_sd_fields(h):
    out = set()
    for entry in h.obj.getProperties().entrySet():
        key, value = str(entry.getKey()), entry.getValue()
        if key.startswith('cdk:') or value is None:
            continue
        out.add((key, str(value)))
    return frozenset(out)


def _cdk_sides(h):
    return (int(h.obj.getReactantCount()), int(h.obj.getAgents().getAtomContainerCount()),
            int(h.obj.getProductCount()))


#: CDK 2.12 states the map number under one key from its SMILES parser and another from its CTfile
#: readers, so a probe that reads one key reports no mapping for a file written by the other.
_CDK_MAP_KEYS = ('molAtomMapNumber', 'cdk:AtomAtomMapping')


def _cdk_map_numbers(h, side):
    out = set()
    for container in side.atomContainers():
        for atom in container.atoms():
            for key in _CDK_MAP_KEYS:
                value = atom.getProperty(key)
                if value:
                    out.add(int(str(value)))
                    break
    return out


def _cdk_mapping(h):
    left = _cdk_map_numbers(h, h.obj.getReactants())
    right = _cdk_map_numbers(h, h.obj.getProducts())
    return _mapping_value(left, right)


def _cdk_coordinates(h):
    return sum(1 for a in h.obj.atoms() if a.getPoint2d() is not None or a.getPoint3d() is not None)


CDK_PROBES = {
    'stereo_groups': _cdk_stereo_groups,
    'dat_sgroup': lambda h: _cdk_sgroups(h, True),
    'other_sgroups': lambda h: _cdk_sgroups(h, False),
    'wedge': _cdk_wedges,
    'either_bond': lambda h: sum(1 for _, n in _cdk_wedges(h) if n == 'either'),
    'charge': lambda h: tuple(sorted(int(a.getFormalCharge()) for a in h.obj.atoms()
                                     if a.getFormalCharge())),
    'radical': lambda h: int(h.obj.getSingleElectronCount()),
    'isotope': lambda h: tuple(sorted((str(a.getSymbol()), int(a.getMassNumber()))
                                      for a in h.obj.atoms() if a.getMassNumber() is not None)),
    'dative': _cdk_dative,
    'aromatic_input': lambda h: tuple(sorted(_cdk_order(b) for b in h.obj.bonds())),
    'implicit_h': _cdk_implicit_h, 'alias': _cdk_alias, 'fragment_group': _cdk_fragment_group,
    'sd_fields': _cdk_sd_fields, 'reaction_sides': _cdk_sides, 'mapping': _cdk_mapping,
    'element_multiset': lambda h: tuple(sorted(str(a.getSymbol()) for a in h.obj.atoms())),
    'coordinate_count': _cdk_coordinates,
}


# ------------------------------------------------------------------------------------- Marvin's probes
# Every Marvin probe reads the MRV `molconvert` re-emitted, because that is the one target on this
# platform that carries all of these features back out.

def _mv_stereo_groups(mrv):
    """`mrvStereoGroup="0 and1 0 0 0 and1 0 0"` -- column or element form, `0` meaning no collection."""
    from re import fullmatch
    kinds = {'or': 2, 'and': 3}
    seen = {}
    for molecule in _mrv_molecules(mrv):
        for atom in _mrv_atoms(molecule):
            token = atom.get('mrvStereoGroup', '0')
            matched = fullmatch(r'(abs|or|and)(\d*)', token)
            if matched is None or matched.group(1) not in kinds:
                continue
            key = (kinds[matched.group(1)], matched.group(2))
            seen[key] = seen.get(key, 0) + 1
    return tuple(sorted((kind, n) for (kind, _), n in seen.items()))


def _mv_dat(mrv):
    out = set()
    for role, attrs in _mrv_sgroups(mrv):
        if _MRV_SGROUP_ROLES.get(role) != 'DAT':
            continue
        out.add((attrs.get('fieldName', ''), attrs.get('fieldData', ''),
                 len(attrs.get('atomRefs', '').split())))
    return frozenset(out)


def _mv_other_sgroups(mrv):
    out = set()
    for role, attrs in _mrv_sgroups(mrv):
        kind = _MRV_SGROUP_ROLES.get(role, role)
        if kind == 'DAT':
            continue
        out.add((kind, len(attrs.get('atomRefs', '').split())))
    return frozenset(out)


def _mv_wedges(mrv):
    out = []
    for molecule in _mrv_molecules(mrv):
        for bond in _mrv_bonds(molecule):
            name = _mrv_wedge(bond)
            if name is not None:
                out.append((_mrv_order(bond.get('order', '1')), name))
    return tuple(sorted(out))


def _mv_atom_values(mrv, key, cast=int):
    """Every non-default value of one atom attribute, cast, across every top-level molecule."""
    out = []
    for molecule in _mrv_molecules(mrv):
        for atom in _mrv_atoms(molecule):
            raw = atom.get(key, '0')
            if raw in ('0', '', 'none'):
                continue
            out.append((atom.get('elementType', ''), cast(raw)))
    return out


def _mv_implicit_h(mrv):
    """MRV's `hydrogenCount` is the implicit count and Marvin omits it where it derives one.

    An omitted attribute is therefore "derive it", not "zero", and the value comparable with the other
    four is the total -- which for a molecule Marvin re-emits is what chython computes from the same
    MRV.  Marvin's own statement is kept where it made one.

    An MRV holding no molecule states nothing about any hydrogen, and `()` is that: the classifier then
    records the cell against whatever the other side wrote.
    """
    records = read_mrv(mrv)
    if len(records) != 1:
        return ()
    return _ch_implicit_h(records[0])


def _mv_sides(mrv):
    sides = _mrv_sides(mrv)
    return (len(sides.get('reactantList', ())), len(sides.get('agentList', ())),
            len(sides.get('productList', ())))


def _mv_mapping(mrv):
    sides = _mrv_sides(mrv)

    def numbers(name):
        out = set()
        for molecule in sides.get(name, ()):
            for atom in _mrv_atoms(molecule):
                value = int(atom.get('mrvMap', '0') or 0)
                if value:
                    out.add(value)
        return out

    left, right = numbers('reactantList'), numbers('productList')
    return _mapping_value(left, right)


def _mv_sd_fields(mrv):
    """Marvin carries SD fields into MRV as a `<propertyList>` of `<property>`/`<scalar>` pairs."""
    out = set()
    for prop in _mrv_tree(mrv).iter('property'):
        name = prop.get('title', '')
        for scalar in prop:
            out.add((name, (scalar.text or '').strip()))
    return frozenset(out)


def _mv_elements(mrv):
    return tuple(sorted(a.get('elementType', '') for m in _mrv_molecules(mrv) for a in _mrv_atoms(m)))


def _mv_coordinates(mrv):
    return sum(1 for m in _mrv_molecules(mrv) for a in _mrv_atoms(m)
               if 'x2' in a or 'x3' in a)


def _mv_fragment_group(mrv):
    records = read_mrv(mrv)
    if len(records) != 1:
        return ()
    return _ch_fragment_group(records[0])


MARVIN_PROBES = {
    'stereo_groups': _mv_stereo_groups, 'dat_sgroup': _mv_dat, 'other_sgroups': _mv_other_sgroups,
    'wedge': _mv_wedges,
    'either_bond': lambda mrv: sum(1 for _, n in _mv_wedges(mrv) if n == 'either'),
    'charge': lambda mrv: tuple(sorted(v for _, v in _mv_atom_values(mrv, 'formalCharge'))),
    'radical': lambda mrv: len(_mv_atom_values(mrv, 'radical', str)),
    'isotope': lambda mrv: tuple(sorted(_mv_atom_values(mrv, 'isotope'))),
    'dative': lambda mrv: sum(1 for m in _mrv_molecules(mrv) for b in _mrv_bonds(m)
                              if b.get('convention') == 'cxn:coord'),
    'aromatic_input': lambda mrv: tuple(sorted(_mrv_order(b.get('order', '1'))
                                               for m in _mrv_molecules(mrv) for b in _mrv_bonds(m))),
    'implicit_h': _mv_implicit_h,
    'alias': lambda mrv: frozenset(v for _, v in _mv_atom_values(mrv, 'mrvAlias', str)),
    'fragment_group': _mv_fragment_group, 'sd_fields': _mv_sd_fields,
    'reaction_sides': _mv_sides, 'mapping': _mv_mapping,
    'element_multiset': _mv_elements, 'coordinate_count': _mv_coordinates,
}


#: One probe table per toolkit, same keys in all five.
PROBES = {'chython': CHYTHON_PROBES, 'rdkit': RDKIT_PROBES, 'indigo': INDIGO_PROBES,
          'cdk': CDK_PROBES, 'marvin': MARVIN_PROBES}


# ---------------------------------------------------------------------------------------------------
# The emit layer.  One seed molecule, written by every toolkit that has a writer for the format, so the
# INBOUND half of the matrix reads a file a reference implementation produced rather than one written by
# hand -- a hand-written block is how a harness manufactures a finding about its own typing.
# ---------------------------------------------------------------------------------------------------

class Seed(NamedTuple):
    """One fixture structure.

    `cxsmiles` is the single source: every toolkit parses it with its own parser, so a seed no toolkit
    can build simply has no cells for that toolkit.  `citation` is required -- a structure without one
    is not a public compound as far as this module is concerned.  `post` is the preparation the writers
    need: `coords` is a 2D layout, which every CTfile wedge and MRV `x2` needs, and `kekule` is a layout
    plus an explicit Kekule form.
    """
    key: str
    cxsmiles: str
    citation: str
    features: frozenset
    post: str


SEEDS = (
    Seed('glycerol', 'OCC(O)CO', 'glycerol, PubChem CID 753',
         frozenset({'sd_fields'}), 'coords'),
    Seed('aspirin', 'CC(=O)Oc1ccccc1C(=O)O', 'acetylsalicylic acid, PubChem CID 2244',
         frozenset({'aromatic_input', 'dat_sgroup', 'implicit_h'}), 'coords'),
    Seed('hexanediol', 'C[C@H](O)CC[C@@H](O)C |&1:1,5|', 'hexane-2,5-diol, PubChem CID 12278',
         frozenset({'stereo_groups', 'wedge'}), 'coords'),
    Seed('alanine', 'C[C@@H](N)C(=O)O', 'L-alanine, PubChem CID 5950',
         frozenset({'wedge', 'stereo_groups'}), 'coords'),
    Seed('glycine_zwitterion', '[NH3+]CC(=O)[O-]', 'glycine zwitterion, PubChem CID 750',
         frozenset({'charge', 'implicit_h'}), 'coords'),
    Seed('nitric_oxide', '[N]=O |^1:0|', 'nitric oxide, PubChem CID 145068',
         frozenset({'radical'}), 'coords'),
    Seed('methanol_13c', '[13CH3]O', 'methanol-13C, PubChem CID 12220 labelled',
         frozenset({'isotope', 'implicit_h'}), 'coords'),
    Seed('ammonia_borane', '[NH3]->[BH3]', 'ammonia borane, PubChem CID 24863',
         frozenset({'dative'}), 'coords'),
    Seed('sodium_acetate', 'CC(=O)[O-].[Na+] |f:0.1|', 'sodium acetate, PubChem CID 517045',
         frozenset({'fragment_group', 'charge'}), 'coords'),
    Seed('pyrrole', 'c1cc[nH]c1', 'pyrrole, PubChem CID 8027',
         frozenset({'implicit_h', 'aromatic_input'}), 'coords'),
    Seed('caffeine', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'caffeine, PubChem CID 2519',
         frozenset({'aromatic_input'}), 'kekule'),
    Seed('ethanol', 'CCO', 'ethanol, PubChem CID 702', frozenset({'alias'}), 'coords'),
    # the coordinate formats: no bond block to carry a feature, so the seed declares the two things a
    # file of atoms and coordinates can state
    Seed('cysteine', 'N[C@@H](CS)C(=O)O', 'L-cysteine, PubChem CID 5862',
         frozenset({'element_multiset', 'coordinate_count'}), 'coords'),
    Seed('esterification',
         '[CH3:1][C:2](=[O:3])[OH:4].[OH:5][CH2:6][CH3:7]>[H+]>'
         '[CH3:1][C:2](=[O:3])[O:5][CH2:6][CH3:7].[OH2:4]',
         'Fischer esterification of acetic acid with ethanol; all five species public compounds '
         '(CIDs 176, 702, 1038, 8857, 962)',
         frozenset({'reaction_sides', 'mapping'}), 'coords'),
)

SEED_KEYS = tuple(s.key for s in SEEDS)


def seed_formats(seed):
    """Every format a cell for this seed can live in -- the union over the features it declares."""
    out = set()
    for feature in seed.features:
        out |= FEATURE_FORMATS[feature]
    return out


def seed_of(key):
    for seed in SEEDS:
        if seed.key == key:
            return seed
    raise KeyError(key)


def is_reaction(seed):
    return '>' in seed.cxsmiles


# ----------------------------------------------------------------------------------- chython's writers

#: The atom alias every writer that can state one attaches, on the first atom.
ALIAS_TEXT = 'Me'

#: The data fields every SD writer attaches for a seed that declares `sd_fields`: one single-line
#: value, and one whose second line is the thing a reader can lose.
SD_FIELDS = (('BATCH_ID', 'lot-42'), ('NOTES', 'first line\nsecond line'))


def attach_dat_sgroup(m, name='BATCH_ID', data=b'lot-42'):
    """A DAT record on the first two atoms.

    `disp` is left unset: chython supplies no anchor and the writer omits FIELDDISP, which is a file
    every one of the four readers accepts.  A hand-written anchor tail is what made Indigo report
    `Expected 'A' or 'D' but got ' '` -- the layout is a reference writer's job, not a fixture's.
    """
    group = SGroup(type='DAT', index=1)
    group.atoms = list(m.atom_numbers)[:2]
    group.name = name
    group.data = [data]
    SGroupStore([group]).to_molecule(m)
    return m


def build_chython(seed):
    """The seed as a chython container, prepared per `post` and carrying its DAT group when it has one."""
    if is_reaction(seed):
        record = read_reaction_smiles(seed.cxsmiles)
        for molecule in record.molecules():
            molecule.clean2d()
        return record
    record = read_smiles(seed.cxsmiles)
    if seed.post == 'kekule':
        record.kekule()
    record.clean2d()
    if 'dat_sgroup' in seed.features:
        attach_dat_sgroup(record)
    if 'alias' in seed.features:
        record.set_aliases({next(iter(record.atom_numbers)): ALIAS_TEXT.encode()})
    return record


#: `clean2d()` runs before every CTfile and MRV write.  MEASURED: with no coordinates chython's V3000
#: carries the `MDLV30/STERAC1` collection but no `CFG=` on any bond, so the parities are not in the file
#: and every outbound stereo cell would fail for a reason about the fixture rather than the writer.
_CHYTHON_WRITE = {
    'v2000': lambda r: mol(r, version=2000),
    'v3000': lambda r: mol(r, version=3000),
    'sdf': lambda r: mol(r, version=2000) + '\n$$$$\n',
    'mrv': write_mrv,
    'cml': write_cml,
    'cxsmiles': lambda r: write_reaction_smiles(r) if hasattr(r, 'reactants') else write_smiles(r),
    'inchi': molecule_to_inchi,
    'rxn': lambda r: rxn(r, version=3000),
}


def write_chython(seed, fmt):
    record = build_chython(seed)
    if fmt == 'sdf' and 'sd_fields' in seed.features:
        # the field layout under test is `SDFWrite`'s own.  The other SD cells are framed by `mol()`
        # above, because `emit_record` refuses an aromatic bond by design and the aromatic seeds are
        # there to measure exactly that bond reaching a reference reader
        buffer = StringIO()
        with SDFWrite(buffer) as handle:
            handle.write(record, meta=dict(SD_FIELDS))
        return buffer.getvalue()
    if 'sd_fields' in seed.features and fmt in ('mrv', 'cml'):
        # the XML writers take the fields off `meta`, which is the container's own channel for them
        record.meta.update(SD_FIELDS)
    return _CHYTHON_WRITE[fmt](record)


# ------------------------------------------------------------------------------------ RDKit's writers

def _rd_build(seed):
    from rdkit.Chem import AllChem, MolFromSmiles, SanitizeMol
    if is_reaction(seed):
        reaction = AllChem.ReactionFromSmarts(seed.cxsmiles, useSmiles=True)
        for template in (*reaction.GetReactants(), *reaction.GetAgents(), *reaction.GetProducts()):
            SanitizeMol(template)
            AllChem.Compute2DCoords(template)
        return reaction
    molecule = MolFromSmiles(seed.cxsmiles)
    if molecule is None:
        raise Unrepresentable(seed.key, 'rdkit', 'MolFromSmiles returned None for the seed')
    if 'sd_fields' in seed.features:
        for name, value in SD_FIELDS:
            molecule.SetProp(name, value)
    if 'alias' in seed.features:
        molecule.GetAtomWithIdx(0).SetProp('molFileAlias', ALIAS_TEXT)
    if seed.post == 'kekule':
        from rdkit.Chem import Kekulize
        Kekulize(molecule, clearAromaticFlags=True)
    AllChem.Compute2DCoords(molecule)
    return molecule


def _rd_sdf(molecule):
    from rdkit.Chem import SDWriter
    buffer = StringIO()
    writer = SDWriter(buffer)
    writer.write(molecule)
    writer.close()
    return buffer.getvalue()


#: The `rdkit.Chem` name of each writer.  RDKit states no CML and no RDF writer, so those cells do not
#: exist rather than failing.
_RDKIT_WRITE = {'v2000': 'MolToMolBlock', 'v3000': 'MolToV3KMolBlock', 'sdf': _rd_sdf,
                'cxsmiles': 'MolToCXSmiles', 'inchi': 'MolToInchi',
                'rxn': 'rdkit.Chem.AllChem:ReactionToV3KRxnBlock', 'xyz': 'MolToXYZBlock',
                'pdb': 'MolToPDBBlock'}


def _write_rdkit(seed, fmt):
    from importlib import import_module
    writer = _RDKIT_WRITE[fmt]
    if not isinstance(writer, str):
        return writer(_rd_build(seed))
    module, _, name = writer.rpartition(':')
    return getattr(import_module(module or 'rdkit.Chem'), name)(_rd_build(seed))


# ----------------------------------------------------------------------------------- Indigo's writers

def _in_build(seed, session):
    if is_reaction(seed):
        record = session.loadReaction(seed.cxsmiles)
    else:
        record = session.loadMolecule(seed.cxsmiles)
        if seed.post == 'kekule':
            record.dearomatize()
        if 'sd_fields' in seed.features:
            for name, value in SD_FIELDS:
                record.setProperty(name, value)
    record.layout()
    return record


def _in_saved(session, record, target):
    """Indigo's saver route, the one that carries a property list into an SD or RD frame."""
    buffer = session.writeBuffer()
    saver = session.createSaver(buffer, target)
    saver.append(record)
    saver.close()
    return buffer.toString()


def _write_indigo(seed, fmt):
    from indigo import Indigo
    from indigo.inchi import IndigoInchi
    session = Indigo()
    session.setOption('ignore-stereochemistry-errors', True)
    record = _in_build(seed, session)
    if fmt == 'v2000':
        session.setOption('molfile-saving-mode', '2000')
        return record.molfile()
    elif fmt == 'v3000':
        # a session-global option, so it is put back before returning
        session.setOption('molfile-saving-mode', '3000')
        try:
            return record.molfile()
        finally:
            session.setOption('molfile-saving-mode', '2000')
    elif fmt == 'sdf':
        return _in_saved(session, record, 'sdf')
    elif fmt == 'cml':
        return record.cml()
    elif fmt == 'cxsmiles':
        return record.smiles()
    elif fmt == 'inchi':
        return IndigoInchi(session).getInchi(record)
    elif fmt == 'rxn':
        return record.rxnfile()
    raise KeyError(fmt)


_INDIGO_WRITE_FORMATS = ('v2000', 'v3000', 'sdf', 'cml', 'cxsmiles', 'inchi', 'rxn')


# -------------------------------------------------------------------------------------- CDK's writers

def _cdk_build(seed, JClass, fmt=''):
    builder = JClass('org.openscience.cdk.silent.SilentChemObjectBuilder').getInstance()
    parser = JClass('org.openscience.cdk.smiles.SmilesParser')(builder)
    if is_reaction(seed):
        record = parser.parseReactionSmiles(seed.cxsmiles)
        layout = JClass('org.openscience.cdk.layout.StructureDiagramGenerator')()
        for side in (record.getReactants(), record.getAgents(), record.getProducts()):
            for container in side.atomContainers():
                layout.generateCoordinates(container)
        return record
    record = parser.parseSmiles(seed.cxsmiles)
    JClass('org.openscience.cdk.layout.StructureDiagramGenerator')().generateCoordinates(record)
    if 'sd_fields' in seed.features:
        for name, value in SD_FIELDS:
            record.setProperty(name, value)
    if 'alias' in seed.features:
        _cdk_alias_atom(JClass, record)
    if fmt in ('xyz', 'pdb', 'mol2'):
        # a coordinate format has no bond block, and CDK's writers for the three read `getPoint3d`;
        # the 2D layout is promoted with z = 0 rather than a conformer being generated, so the file
        # states exactly the geometry the layout produced
        point3d = JClass('javax.vecmath.Point3d')
        for atom in record.atoms():
            flat = atom.getPoint2d()
            atom.setPoint3d(point3d(flat.x, flat.y, 0.))
    return record


def _cdk_alias_atom(JClass, record):
    """Label the first atom.  CDK 2.12 states an atom label on an `IPseudoAtom`, and `setAtomicNumber`
    takes a boxed `Integer`, so the number is boxed explicitly rather than passed as a Python int."""
    from jpype import JInt, JObject

    old = record.getAtom(0)
    pseudo = JClass('org.openscience.cdk.silent.PseudoAtom')(ALIAS_TEXT)
    pseudo.setSymbol(old.getSymbol())
    pseudo.setAtomicNumber(JObject(int(old.getAtomicNumber()), JInt))
    pseudo.setPoint2d(old.getPoint2d())
    pseudo.setImplicitHydrogenCount(old.getImplicitHydrogenCount())
    JClass('org.openscience.cdk.tools.manipulator.AtomContainerManipulator'
           ).replaceAtomByAtom(record, old, pseudo)
    return record


def _cdk_written(JClass, writer_name, record):
    writer_out = JClass('java.io.StringWriter')()
    writer = JClass(f'org.openscience.cdk.io.{writer_name}')(writer_out)
    writer.write(record)
    writer.close()
    return str(writer_out.toString())


def _cdk_flavor(JClass):
    """`SmiFlavor` bits for a CXSMILES carrying the collections, the fragment groups and the DAT data.

    MEASURED: `SmiFlavor` on 2.12 states no `CxSmilesWithAtomLabels`, so the flavour is built from the
    members that are there.
    """
    flavor = JClass('org.openscience.cdk.smiles.SmiFlavor')
    return int(flavor.Absolute) | int(flavor.CxSmiles) | int(flavor.CxEnhancedStereo) \
        | int(flavor.CxFragmentGroup) | int(flavor.CxDataSgroups) | int(flavor.CxRadical) \
        | int(flavor.Cx2dCoordinates)


def _write_cdk(seed, fmt):
    JClass = cdk()
    if JClass is None:
        raise Unrepresentable(seed.key, 'cdk', 'the CDK jar is absent')
    record = _cdk_build(seed, JClass, fmt)
    if fmt == 'v2000':
        return _cdk_written(JClass, 'MDLV2000Writer', record)
    elif fmt == 'v3000':
        return _cdk_written(JClass, 'MDLV3000Writer', record)
    elif fmt == 'sdf':
        return _cdk_written(JClass, 'SDFWriter', record)
    elif fmt == 'cml':
        return _cdk_written(JClass, 'CMLWriter', record)
    elif fmt == 'rxn':
        # CDK 2.12 states `MDLRXNWriter` and no V3000 RXN writer, so this cell is a V2000 RXN
        return _cdk_written(JClass, 'MDLRXNWriter', record)
    elif fmt == 'xyz':
        return _cdk_written(JClass, 'XYZWriter', record)
    elif fmt == 'pdb':
        return _cdk_written(JClass, 'PDBWriter', record)
    elif fmt == 'mol2':
        return _cdk_written(JClass, 'Mol2Writer', record)
    elif fmt == 'cxsmiles':
        generator = JClass('org.openscience.cdk.smiles.SmilesGenerator')(_cdk_flavor(JClass))
        return str(generator.create(record))
    elif fmt == 'inchi':
        factory = JClass('org.openscience.cdk.inchi.InChIGeneratorFactory').getInstance()
        return str(factory.getInChIGenerator(record).getInchi())
    raise KeyError(fmt)


_CDK_WRITE_FORMATS = ('v2000', 'v3000', 'sdf', 'cml', 'cxsmiles', 'inchi', 'rxn', 'xyz', 'pdb', 'mol2')


# ----------------------------------------------------------------------------------- Marvin's writers

#: `molconvert` target names, all verified from stdin.
_MARVIN_TARGET = {'v2000': 'mol', 'v3000': 'mol:V3', 'sdf': 'sdf', 'mrv': 'mrv', 'cml': 'cml',
                  'cxsmiles': 'cxsmiles', 'rxn': 'rxn', 'rdf': 'rdf', 'xyz': 'xyz', 'pdb': 'pdb',
                  'mol2': 'mol2'}


def _marvin_source(seed, fmt):
    """The framed input for a feature a SMILES cannot state, or None for the seed's SMILES.

    A SMILES carries neither a data field nor an atom alias, so the only Marvin-written file with one in
    it is one `molconvert` re-emits from a framed input.  The layout is then Marvin's and the feature is
    chython's, which is stated wherever such a fixture appears.
    """
    if 'sd_fields' in seed.features and fmt in ('sdf', 'mrv', 'cml'):
        return write_chython(seed, 'sdf')
    if 'alias' in seed.features and fmt in ('v2000', 'v3000', 'sdf', 'mrv'):
        return write_chython(seed, 'mrv')
    return None


def _write_marvin(seed, fmt):
    """One `molconvert` call.  `-2` computes the 2D layout the CTfile wedge and MRV `x2` need; `-3` is
    used for the coordinate formats, which carry no bond block to hold a parity."""
    extra = ('-3',) if fmt in ('xyz', 'pdb', 'mol2') else ('-2',)
    source = _marvin_source(seed, fmt)
    if source is None:
        source = seed.cxsmiles
    else:
        extra = ()
    code, out, err, argv = molconvert(source, _MARVIN_TARGET[fmt], extra)
    if code or not out.strip():
        raise Unrepresentable(seed.key, 'marvin',
                              f'{" ".join(argv)} exited {code}: {err.strip()[:200]}')
    return out


#: Every writer the matrix has, keyed `(toolkit, fmt)`.  An absent key means that toolkit states no
#: writer for that format, so the cell does not exist -- MRV among these four is Marvin's and chython's
#: only, and RDF is Marvin's only.
WRITERS = {('chython', fmt): (lambda f: lambda s: write_chython(s, f))(fmt) for fmt in _CHYTHON_WRITE}
WRITERS.update({('rdkit', fmt): (lambda f: lambda s: _write_rdkit(s, f))(fmt)
                for fmt in _RDKIT_WRITE})
WRITERS.update({('indigo', fmt): (lambda f: lambda s: _write_indigo(s, f))(fmt)
                for fmt in _INDIGO_WRITE_FORMATS})
WRITERS.update({('cdk', fmt): (lambda f: lambda s: _write_cdk(s, f))(fmt)
                for fmt in _CDK_WRITE_FORMATS})
WRITERS.update({('marvin', fmt): (lambda f: lambda s: _write_marvin(s, f))(fmt)
                for fmt in _MARVIN_TARGET})


def write_once(cache, toolkit, fmt, seed):
    """The fixture text for one cell, written at most once per session.

    A refusal is cached as the exception object and re-raised, so a `molconvert` that already failed is
    not spent again by the four readers that would each ask for the same bytes.
    """
    key = (toolkit, fmt, seed.key)
    if key not in cache:
        try:
            cache[key] = WRITERS[toolkit, fmt](seed)
        except Exception as error:
            cache[key] = error
    got = cache[key]
    if isinstance(got, BaseException):
        raise got
    return got


def parse_once(cache, toolkit, fmt, text):
    """One reader's answer for one exact file, at most once per session.

    Keyed on the text itself, because the same bytes are asked about by every feature of a seed and a
    `molconvert` readback costs about 0.6 s -- the cache is what keeps the matrix inside the suite's
    time budget rather than beside it.
    """
    key = (toolkit, fmt, hash(text))
    if key not in cache:
        cache[key] = READERS[toolkit](text, fmt)
    return cache[key]


def main(argv=()):
    """`python -m chython.formats.test.oracles` -- which oracles this machine has, and their versions.

    This command IS the harness behind any conformance claim about formats: a coverage claim ships with
    the harness that produced it or it is not made.
    """
    if 'matrix' in argv:
        return _matrix()
    present = []
    for name in TOOLKITS:
        version = _PROBES[name]()
        if version is None:
            pin = PINS.get(name)
            hint = f'set ${pin}' if pin else f'pip install {name}'
            print(f'{name:8} absent ({hint})')
        else:
            print(f'{name:8} {version}')
            present.append(name)
    if 'cdk' in present:
        print(f'{"":8} jar {CDK_JAR}')
    if 'marvin' in present:
        print(f'{"":8} bin {MOLCONVERT}')
    if not present:
        print('no oracle present; every conformance cell would skip')
        return 1
    return 0


def _matrix():
    """The full matrix as markdown.  Defined in the report section below."""
    from .test_conformance import print_matrix
    return print_matrix()


if __name__ == '__main__':
    from sys import argv as _argv
    sys_exit(main(_argv[1:]))
