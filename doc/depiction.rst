Depiction
=========

SVG rendering, 2D coordinate generation, and display settings.


SVG Output
----------

.. code-block:: python

    from chython import smiles

    mol = smiles('c1ccccc1O')

    # Generate SVG (auto-calculates 2D coords if needed)
    svg = mol.depict()

    # With explicit size
    svg = mol.depict(width='10cm', height='10cm')

    # PNG (requires pyppeteer)
    png = mol.depict(format='png', png_width=1000, png_heigh=1000)

    # Compressed SVG
    svgz = mol.depict(format='svgz')

    # Jupyter notebook: molecules render automatically via _repr_svg_

    # Reactions too
    rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    svg = rxn.depict()


2D Coordinate Generation
-------------------------

.. code-block:: python

    import chython

    # Set engine globally (before calling clean2d)
    chython.clean2d_engine = 'smilesdrawer'  # default, built-in
    # Other options: 'rdkit', 'cdk', 'obabel', 'indigo'

    mol = smiles('c1ccccc1')
    mol.clean2d()  # generates 2D coordinates


Depiction Settings
------------------

.. code-block:: python

    from chython import depict_settings

    depict_settings(
        carbon=False,        # hide C labels (default)
        aam=True,            # show atom-atom mapping
        monochrome=False,    # use CPK colors
        bond_color='black',
        font_size=0.5,
        bond_width=0.04,
    )

    # Restore defaults
    depict_settings()

After changing settings, flush cached depictions:

.. code-block:: python

    mol.flush_cache()
    svg = mol.depict()


Grid Depiction
--------------

.. code-block:: python

    from chython import grid_depict, smiles

    mols = [smiles('CCO'), smiles('c1ccccc1'), smiles('CC(=O)O')]
    svg = grid_depict(mols, cols=3)

    # In Jupyter, wrap with ipywidgets for display
    from ipywidgets import HTML
    HTML(grid_depict(mols, cols=3))
