.. include:: ../README.rst

Chython package API
===================

chython. **pickle_cache** = True - Store cached attributes in pickle. Effective for multiprocessing.

chython. **torch_device** = 'cpu' - Atom-to-Atom mapping model device in torch notation. Change before first `reset_mapping` call!

.. automodule:: chython
    :members: smiles, inchi, xyz, mdl_mol, smarts, depict_settings
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::
    :maxdepth: 4

    containers
    files
    reactor
    utils
    periodictable

Notebooks
=========

.. toctree::
    :caption: Tutorial
    :maxdepth: 1

    tutorial/notebook.ipynb
