#!/usr/bin/env python3
"""QM property contour plots over a 2D depiction -- the gallery behind docs/depiction.rst's field section.

Every number here is COMPUTED, by the smallest quantum model that gives the picture its meaning: Hückel
π theory on the molecule's own adjacency matrix.  Four scalars come out of one diagonalization --

    π charge q_i      = 1 − 2 Σ_occ c_ir²        signed, so a diverging map with zero at its midpoint
    HOMO c_i          = coefficients of ψ_HOMO   signed, and the only field here with a NODAL LINE
    Fukui f⁻_i        = c_i,HOMO²                one-sided: where an electrophile attacks
    Fukui f⁺_i        = c_i,LUMO²                one-sided: where a nucleophile attacks
    self-polarizability π_ii (Coulson–Longuet-Higgins), the near-uniform case

-- and a DFT or semiempirical calculation drops in at exactly the same seam: `AtomField` takes a
`{atom_id: float}` dict and knows nothing about where the floats came from.

Run: python docs/demo_qm/render_fields.py   ->  docs/demo_qm/index.html
"""
from math import isclose
from pathlib import Path

import numpy as np

from chython import DepictStyle, smiles
from chython.depict.overlay import AtomField, AtomHalo, Highlight, ValueLabels


OUT = Path(__file__).parent            # the gallery is written beside this script, in `docs/demo_qm/`

#: Rendered at 9 mm per molecule unit: a contour band is a large smooth shape and reads better big.
STYLE = DepictStyle.preset('screen').tuned(**{'page.scale_mm': 9., 'atom.carbon': False})


class Huckel:
    """One diagonalization of a conjugated skeleton, and every scalar this demo draws.

    The π system is taken as every carbon of the molecule -- true for the two hydrocarbons here and the
    reason they were chosen; a heteroatom needs Streitwieser h/k parameters on the diagonal.
    """

    def __init__(self, smi: str, name: str):
        mol = smiles(smi)
        mol.clean2d()
        self.mol = mol
        self.name = name
        self.ids = [a.n for a in mol.atoms()]
        assert all(a.atomic_symbol == 'C' for a in mol.atoms()), 'π system is the carbons here'

        index = {n: i for i, n in enumerate(self.ids)}
        size = len(self.ids)
        adjacency = np.zeros((size, size))
        for bond in mol.bonds():
            adjacency[index[bond.n], index[bond.m]] = adjacency[index[bond.m], index[bond.n]] = 1.

        # E = α + xβ with β < 0, so the MOST stable orbital has the LARGEST x: sort descending and fill
        # from the front.  One π electron per carbon.
        x, c = np.linalg.eigh(adjacency)
        order = np.argsort(x)[::-1]
        self.x, self.c = x[order], c[:, order]
        self.nocc = size // 2
        assert size % 2 == 0, 'an odd π system has a half-filled orbital and no closed-shell HOMO'
        assert not isclose(self.x[self.nocc - 1], self.x[self.nocc], abs_tol=1e-6), \
            'HOMO and LUMO are degenerate: no single frontier orbital to plot'

    def _dict(self, values) -> dict:
        return dict(zip(self.ids, (float(v) for v in values)))

    @property
    def gap(self) -> float:
        """HOMO–LUMO separation in |β| -- the number azulene's colour is famous for being small."""
        return float(self.x[self.nocc - 1] - self.x[self.nocc])

    def charges(self) -> dict:
        """π charge: positive is electron-POOR.  Sums to zero, and is identically zero for an alternant."""
        return self._dict(1. - 2. * (self.c[:, :self.nocc] ** 2).sum(1))

    def homo(self) -> dict:
        """Signed HOMO coefficients.  Sign is arbitrary up to a global phase; the nodal line is not."""
        return self._dict(self.c[:, self.nocc - 1])

    def fukui_minus(self) -> dict:
        return self._dict(self.c[:, self.nocc - 1] ** 2)

    def fukui_plus(self) -> dict:
        return self._dict(self.c[:, self.nocc] ** 2)

    def self_polarizability(self) -> dict:
        """π_ii = 4 Σ_occ Σ_unocc c_ir² c_is² / (x_r − x_s), in |β|⁻¹.  Positive by construction."""
        size = len(self.ids)
        out = np.zeros(size)
        for r in range(self.nocc):
            for s in range(self.nocc, size):
                out += 4. * self.c[:, r] ** 2 * self.c[:, s] ** 2 / (self.x[r] - self.x[s])
        return self._dict(out)


def write(stem: str, svg: str, caption: str, panels: list):
    (OUT / f'{stem}.svg').write_text(svg)
    panels.append((stem, caption))
    print(f'  {stem}.svg  {len(svg) // 1024:3d} KiB')


def main():
    OUT.mkdir(exist_ok=True)
    panels: list = []

    azulene = Huckel('c1ccc2cccccc12', 'azulene')
    naphthalene = Huckel('c1ccc2ccccc2c1', 'naphthalene')
    print(f'azulene gap {azulene.gap:.3f} |β|, naphthalene gap {naphthalene.gap:.3f} |β|')

    # --- diverging: a signed field, zero at the colormap's neutral midpoint -----------------------------
    # `coolwarm` is diverging and the data spans zero, so `fitted()` symmetrises the domain itself: the
    # grey midpoint lands on q = 0 without the caller stating a domain.
    charges = azulene.charges()
    write('01_charge_azulene',
          azulene.mol.depict(style=STYLE, overlays=[AtomField(charges, colormap='coolwarm')]),
          'π charge on azulene — red electron-poor (seven-ring), blue electron-rich (five-ring). '
          'The dipole azulene is known for, as a field.', panels)

    # Same field, isolines only, and the levels stated: a contour AT zero is the charge-neutral curve.
    write('02_charge_azulene_lines',
          azulene.mol.depict(style=STYLE, overlays=[
              AtomField(charges, colormap='coolwarm', fill=False,
                        levels=[-.12, -.06, 0., .06, .12])]),
          'The same field as line contours (fill=False) at stated levels. The level at 0 is the '
          'charge-neutral curve, and the legend draws rules rather than blocks because no interval '
          'of values was filled.', panels)

    # --- a signed field with a NODAL LINE ---------------------------------------------------------------
    # Naphthalene's HOMO is zero at both fusion carbons; the zero contour is the node, not an artefact.
    homo = naphthalene.homo()
    write('03_homo_naphthalene',
          naphthalene.mol.depict(style=STYLE, overlays=[AtomField(homo, colormap='RdBu')]),
          'HOMO coefficients of naphthalene (RdBu). The two lobes meet along the nodal line through '
          'the fusion carbons, where the coefficient is exactly 0.', panels)

    # A zero level on a DIVERGING map is painted the map's neutral midpoint, which is the page colour: the
    # node is drawn and invisible.  A second overlay carrying only that level on `mono` -- grey at the
    # midpoint of a symmetric domain -- is the way to see it.  Two overlays with two different scales
    # cannot share one colorbar and the figure says so, hence `page.legend='none'`: panel 03 labelled
    # this field already.
    extent = max(abs(v) for v in homo.values())
    write('04_homo_naphthalene_node',
          naphthalene.mol.depict(style=STYLE.tuned(**{'page.legend': 'none'}), overlays=[
              AtomField(homo, colormap='RdBu', clip='box'),
              AtomField(homo, colormap='mono', levels=[0.], fill=False,
                        domain=(-extent, extent), clip='box')]),
          'The node drawn as a contour: level 0.0 on mono is mid-grey, so it is visible where the '
          'diverging map paints it the page colour. Four lobes of alternating sign, nodes along the '
          'fusion axis and across it.', panels)

    # --- one-sided: sequential single hue, and ONE domain for a pair meant to be compared -------------
    f_minus, f_plus = azulene.fukui_minus(), azulene.fukui_plus()
    shared = (0., max(max(f_minus.values()), max(f_plus.values())))
    for stem, values, what in [('05_fukui_minus_azulene', f_minus, 'f⁻ (electrophilic attack, HOMO²)'),
                               ('06_fukui_plus_azulene', f_plus, 'f⁺ (nucleophilic attack, LUMO²)')]:
        write(stem,
              azulene.mol.depict(style=STYLE,
                                 overlays=[AtomField(values, colormap='viridis', domain=shared)]),
              f'Azulene {what}, viridis, one-sided. Both panels carry the SAME domain 0–{shared[1]:.3f}, '
              'so the two are read against one scale — f⁻ peaks on the five-ring, f⁺ on the '
              'seven-ring.', panels)

    # WHEN NOT TO CONTOUR.  An alternant's self-polarizability is nearly uniform (0.330–0.443 |β|⁻¹), and
    # a contour of it is a flat blob with the whole ramp crammed into the decay rim: the field is a
    # Gaussian INTERPOLATION, so it falls to zero away from the atoms whatever the atom values are.  A
    # property with no spatial story gets the per-atom form instead.
    polarizability = naphthalene.self_polarizability()
    write('07_polarizability_contoured',
          naphthalene.mol.depict(style=STYLE, overlays=[
              AtomField(polarizability, colormap='cividis')]),
          'Self-polarizability π_ii contoured — and this is the panel that should not be a contour. '
          'Range 0.330–0.443 |β|⁻¹: the interior is one flat colour and every band sits in the rim '
          'where the interpolation decays to zero.', panels)

    write('08_polarizability_halo',
          naphthalene.mol.depict(style=STYLE, overlays=[
              AtomHalo(polarizability, colormap='cividis', encode='both')]),
          'The same numbers as AtomHalo(encode=\'both\') — radius AND colour per atom. A scalar with no '
          'spatial structure is a per-atom quantity, and this form does not invent one. The α positions '
          'are the polarizable ones.', panels)

    # --- clip modes, on one field so the difference is the clip and nothing else -----------------------
    for stem, clip, note in [('09_clip_none', None, "clip=None — the field's own cutoff closes it"),
                             ('10_clip_hull', 'hull', 'hull, drawn through the atoms: bands are sliced'),
                             ('11_clip_box', 'box', 'box over the label boxes')]:
        write(stem,
              azulene.mol.depict(style=STYLE, overlays=[AtomField(charges, clip=clip)]),
              f'clip={clip!r} — {note}.', panels)

    # --- composition: three overlays, one figure -------------------------------------------------------
    # Labelling every atom collides at the fusion bond and says nothing a reader could not see in the
    # bands: only the four extremes are numbered, which is what the numbers are for.
    five_ring = next(r for r in azulene.mol.sssr if len(r) == 5)
    ranked = sorted(charges, key=lambda n: abs(charges[n]), reverse=True)[:4]
    write('12_combined',
          azulene.mol.depict(style=STYLE, overlays=[
              AtomField(charges, colormap='coolwarm', opacity=.55),
              Highlight(atoms=five_ring, style='outline', label='five-ring'),
              ValueLabels({n: charges[n] for n in ranked}, fmt='{:+.2f}')]),
          'AtomField + Highlight + ValueLabels, with only the four extreme atoms numbered. The field '
          'is dropped to opacity 0.55 so the numbers stay legible on it; the numbers wear the ink '
          'colour, never the band colour.', panels)

    # --- legend off: the structure must not move -------------------------------------------------------
    no_bar = STYLE.tuned(**{'page.legend': 'none'})
    write('13_no_legend',
          azulene.mol.depict(style=no_bar, overlays=[AtomField(charges)]),
          'page.legend=\'none\'. The bar is placed outside the content box, so withholding it moves '
          'no atom.', panels)

    cards = '\n'.join(
        f'    <figure>\n      <img src="{stem}.svg" alt="{stem}">\n'
        f'      <figcaption><b>{stem}</b><br>{caption}</figcaption>\n    </figure>'
        for stem, caption in panels)
    (OUT / 'index.html').write_text(f"""<!doctype html>
<meta charset="utf-8"><title>chython — QM field depiction</title>
<style>
  body {{ font: 15px/1.55 -apple-system, system-ui, sans-serif; margin: 0; padding: 2.5rem;
         background: #fbfbfa; color: #1c1c1a; }}
  h1 {{ font-size: 1.5rem; margin: 0 0 .3rem; }}
  p.lead {{ color: #555; max-width: 62ch; margin: 0 0 2rem; }}
  .grid {{ display: grid; gap: 1.75rem; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); }}
  figure {{ margin: 0; background: #fff; border: 1px solid #e6e4e0; border-radius: 10px;
            padding: 1rem; }}
  figure img {{ width: 100%; height: auto; display: block; }}
  figcaption {{ font-size: 13px; color: #555; margin-top: .75rem; }}
  figcaption b {{ color: #1c1c1a; font-weight: 600; }}
</style>
<h1>QM property fields over a 2D depiction</h1>
<p class="lead">Hückel π theory on each molecule's own adjacency matrix, contoured by
<code>chython.depict.overlay.AtomField</code>. Azulene HOMO–LUMO gap
{azulene.gap:.3f}&nbsp;|β|, naphthalene {naphthalene.gap:.3f}&nbsp;|β|. Regenerate with
<code>python docs/demo_qm/render_fields.py</code>.</p>
<div class="grid">
{cards}
</div>
""")
    print(f'\n{len(panels)} panels -> {OUT / "index.html"}')


if __name__ == '__main__':
    main()
