Group Glossary
==============

Every row of the two group corpora, by name: what each pattern matches, and for a protecting group,
what removing it reveals. This page is generated from
``chython/reactions/tables/functional.tsv`` and ``chython/reactions/tables/protective.tsv``, which are
the authority — a name that appears here is a name the shipped table has.

The API is on :doc:`reactions`: :meth:`chython.MoleculeContainer.functional_groups` and
:meth:`chython.MoleculeContainer.protective_groups` report what a molecule carries,
:meth:`chython.MoleculeContainer.deprotect` applies the removals, and both corpora are also reachable
as tables through :func:`chython.functional_rules` and :func:`chython.protective_rules` — the whole
inventory, keyed by the name below, no molecule involved. The SMARTS column is the chython dialect,
documented in :doc:`substructure`; ``z`` and ``x`` there do not mean what the same letters mean in
another toolkit's dialect.

The reaction corpus and the role table are not listed here: a reaction row is a template rather than a
name a molecule carries, and :doc:`reactions` documents it where it is applied.

.. BEGIN GENERATED GLOSSARY: python -m chython.reactions.test.gen_corpus_glossary
.. Generated from chython/reactions/tables/functional.tsv and protective.tsv, which are
   the authority.  Do not edit below this line -- run the command above.

Functional groups
-----------------

420 functional groups, alphabetically.  The name is the key
:meth:`chython.MoleculeContainer.functional_groups` returns and the key
:func:`chython.functional_rules` is keyed on; the id is what a consumer stores.

.. list-table::
   :header-rows: 1
   :widths: 22 12 30 36

   * - Name
     - Id
     - SMARTS
     - What it matches
   * - ``1_2_dihaloethane``
     - ``functional:311``
     - ``[Cl,Br,I;D1]-[C;z1;x1;D2:1]-[C;z1;x1;D2:2]-[Cl,Br,I;D1]``
     - X-CH2-CH2-X, the cyclopropanating bis-electrophile
   * - ``1_2_diketone``
     - ``functional:75``
     - ``[O;z2;x0:1]=[C;D3;x1:2]-[C;z2;x1;D3:3]=[O:4]``
     - 1,2-diketone
   * - ``1_3_diketone``
     - ``functional:76``
     - ``[O;z2;x0:1]=[C;D3;x1:2]-[C;z1;D2,D3:3]-[C;z2;x1;D3:4]=[O:5]``
     - 1,3-diketone
   * - ``1_4_diketone``
     - ``functional:77``
     - ``[O;z2;x0:1]=[C;D3;x1:2]-[C;z1;D2,D3:3]-[C;z1;D2,D3:4]-[C;z2;x1;D3:5]=[O:6]``
     - 1,4-diketone, the Paal-Knorr substrate
   * - ``NO_dialkylhydroxylamine``
     - ``functional:164``
     - ``[N;D2;z1;x1;h1:1](-[C;z1:2])-[O;D2;z1;x1:3]-[C;z1:4]``
     - R'-O-NH-R'', the Weinreb amine
   * - ``O_alkylhydroxylamine``
     - ``functional:163``
     - ``[N;D1;z1;x1:1]-[O;D2;z1;x1:2]-[C:3]``
     - R-O-NH2
   * - ``acetal``
     - ``functional:186``
     - ``[O;D2;z1:1]-[C;z1;x2;h0,h1:2]-[O;D2;z1:3]``
     - R-O-CH(R)-O-R, the masked carbonyl.  ``h0,h1`` admits the ketal and excludes an orthoester
   * - ``acrylamide``
     - ``functional:195``
     - ``[C;z2;D1:1]=[C;z2:2]-[C;D3;z2;x2:3](=[O;D1;x0:4])-[N:5]``
     - H2C=CH-CO-N, the covalent-warhead acceptor.  ``D1`` on the terminal carbon: a beta-substituted acrylamide is a far weaker acceptor and is not this
   * - ``acrylate_ester``
     - ``functional:196``
     - ``[C;z2;D1:1]=[C;z2:2]-[C;D3;z2;x2:3](=[O;D1;x0:4])-[O;D2:5]``
     - H2C=CH-CO-O-R, the polymerizable monomer.  ``D1`` as in ``acrylamide``
   * - ``activated_alkyl_halide``
     - ``functional:390``
     - ``[Cl,Br,I;D1]-[C;z1;x1;h1,h2:1]-[C;z2,z4;M]``
     - X-CH or X-CH2 on an aromatic or sp2 carbon: benzylic, allylic, beside C=N or C=O.  The halide is the leaving group, unlike ``benzyl_halide`` and ``allyl_halide``, which number it
   * - ``activated_aza_ortho_halide``
     - ``functional:297``
     - ``[Cl,Br,I;D1]-[C;a;r6:1]:[N;a;D2;h0;r6;M]``
     - Ar-X (Cl, Br, I) beside a six-ring aza N: the SNAr electrophile ``aryl_halide`` also matches
   * - ``activated_aza_ortho_nitro``
     - ``functional:376``
     - ``[O;-;D1]-[N;+;D3](=[O;D1])-[C;a;r6:1]:[N;a;D2;h0;r6;M]``
     - Ar-NO2 beside a six-ring aza N: nitrite is the SNAr leaving group
   * - ``activated_aza_ortho_sulfonate``
     - ``functional:301``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;a;r6:1]:[N;a;D2;h0;r6;M])-[C]``
     - ArOSO2R beside a six-ring aza N
   * - ``activated_aza_ortho_sulfone``
     - ``functional:359``
     - ``[C;z1]-[S;D4](=[O;D1])(=[O;D1])-[C;a;r6:1]:[N;a;D2;h0;r6;M]``
     - Ar-SO2R beside a six-ring aza N: the alkylsulfonyl is the SNAr leaving group
   * - ``activated_aza_para_halide``
     - ``functional:298``
     - ``[Cl,Br,I;D1]-[C;a:1]1:[C,N;a;M]:[C,N;a;M]:[N;a;D2;h0;M]:[C,N;a;M]:[C,N;a;M]:1``
     - Ar-X (Cl, Br, I) para to a six-ring aza N
   * - ``activated_aza_para_nitro``
     - ``functional:377``
     - ``[O;-;D1]-[N;+;D3](=[O;D1])-[C;a:1]1:[C,N;a;M]:[C,N;a;M]:[N;a;D2;h0;M]:[C,N;a;M]:[C,N;a;M]:1``
     - Ar-NO2 para to a six-ring aza N: nitrite is the SNAr leaving group
   * - ``activated_aza_para_sulfonate``
     - ``functional:302``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;a:1]1:[C,N;a;M]:[C,N;a;M]:[N;a;D2;h0;M]:[C,N;a;M]:[C,N;a;M]:1)-[C]``
     - ArOSO2R para to a six-ring aza N
   * - ``activated_aza_para_sulfone``
     - ``functional:360``
     - ``[C;z1]-[S;D4](=[O;D1])(=[O;D1])-[C;a:1]1:[C,N;a;M]:[C,N;a;M]:[N;a;D2;h0;M]:[C,N;a;M]:[C,N;a;M]:1``
     - Ar-SO2R para to a six-ring aza N
   * - ``activated_azole_halide``
     - ``functional:305``
     - ``[Cl,Br,I;D1]-[C;a;r5:1](:[N;a;D2;h0;r5;M]):[N,O,S;a;r5;M]``
     - Ar-X (Cl, Br, I) on a five-ring carbon between two ring heteroatoms, one an aza N: 2-haloazoles and benzazoles, C8 of purine
   * - ``activated_azole_sulfonate``
     - ``functional:306``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;a;r5:1](:[N;a;D2;h0;r5;M]):[N,O,S;a;r5;M])-[C]``
     - ArOSO2R on a five-ring carbon between two ring heteroatoms, one an aza N
   * - ``activated_azole_sulfone``
     - ``functional:361``
     - ``[C;z1]-[S;D4](=[O;D1])(=[O;D1])-[C;a;r5:1](:[N;a;D2;h0;r5;M]):[N,O,S;a;r5;M]``
     - Ar-SO2R on a five-ring carbon between two ring heteroatoms, one an aza N: 2-sulfonyl oxadiazoles, thiazoles, benzazoles
   * - ``activated_ester``
     - ``functional:264``
     - ``[O;z2;x0:1]=[C;D3;x2;z2:2]-[O;D2;x1]-[N]``
     - an N-O ester -- NHS, HOBt, HOAt, oxime -- the acylating agent ``ester`` (``x0`` oxygen) leaves out
   * - ``activated_isocyanide``
     - ``functional:166``
     - ``[C;-;D1:1]#[N;+;D2:2]-[C;z1;h1,h2:3]``
     - isocyanide with an acidic alpha C-H; what activates it is left unstated
   * - ``activated_nitro_ortho_halide``
     - ``functional:299``
     - ``[Cl,Br,I;D1]-[C;a:1]:[C;a;M]-[N;+;D3;M](=[O;M])-[O;-;M]``
     - Ar-X (Cl, Br, I) ortho to a nitro group
   * - ``activated_nitro_ortho_sulfonate``
     - ``functional:303``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;a:1]:[C;a;M]-[N;+;D3;M](=[O;M])-[O;-;M])-[C]``
     - ArOSO2R ortho to a nitro group
   * - ``activated_nitro_para_halide``
     - ``functional:300``
     - ``[Cl,Br,I;D1]-[C;a:1]1:[C;a;M]:[C;a;M]:[C;a;M](-[N;+;D3;M](=[O;M])-[O;-;M]):[C;a;M]:[C;a;M]:1``
     - Ar-X (Cl, Br, I) para to a nitro group
   * - ``activated_nitro_para_sulfonate``
     - ``functional:304``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;a:1]1:[C;a;M]:[C;a;M]:[C;a;M](-[N;+;D3;M](=[O;M])-[O;-;M]):[C;a;M]:[C;a;M]:1)-[C]``
     - ArOSO2R para to a nitro group
   * - ``active_ch``
     - ``functional:309``
     - ``[C;z1;h1,h2:1](-[C,S,N;*;M]=,#[O,N;M])-[C,S,N;*;M]=,#[O,N;M]``
     - a C-H between two C=O, C#N, S=O or N=O groups, alkylated under carbonate or hydride base
   * - ``active_methylene``
     - ``functional:167``
     - ``[C;z1;D2,D3;x0:1](-[C;z2,z3;x1,x2:2])-[C;z2,z3;x1,x2:3]``
     - CH flanked by two EWGs, the Knoevenagel nucleophile
   * - ``acyclic_ester``
     - ``functional:259``
     - ``[O;z2;x0:1]=[C;D3;x2;z2:2]-;!@[O;D2;x0]``
     - ``ester`` outside a ring, so hydrolysis leaves one product; a lactone opens to a hydroxy acid instead
   * - ``acyl_azide``
     - ``functional:240``
     - ``[N;D1;z2;-:1]=[N;D2;+:2]=[N;D2;z2:3]-[C;D3;z2:4]=[O;D1;x0:5]``
     - R-CO-N3, the Curtius substrate, as distinct from the alkyl azide
   * - ``acyl_bromide``
     - ``functional:84``
     - ``[Br;D1]-[C;z2;x2;D3:1]=[O:2]``
     - R-COBr
   * - ``acyl_chloride``
     - ``functional:23``
     - ``[Cl;D1][C;z2;x2;D3:1]=[O:2]``
     - R-COCl
   * - ``acyl_fluoride``
     - ``functional:85``
     - ``[F;D1]-[C;z2;x2;D3:1]=[O:2]``
     - R-COF
   * - ``acyl_hemiaminal``
     - ``functional:398``
     - ``[O;D1;z1;x0:1]-[C;z1;x2:2]-[N;z1;D3;M]-[C;z2;M]=[O;M]``
     - HO-C-N(acyl), the N-acyliminium precursor
   * - ``aldehyde``
     - ``functional:20``
     - ``[O;z2;x0:1]=[C;D2;x1;z2:2]``
     - R-CHO
   * - ``aldoxime_ch``
     - ``functional:400``
     - ``[O;D1,D2;z1;x1:1]-[N;D2;z2;x1:2]=[C;D2;h1:3]``
     - an aldoxime or its O-alkyl ether with the C-H numbered, the hydroximoyl halide precursor
   * - ``alkene``
     - ``functional:2``
     - ``[C;z2;x0;D2,D3:1]=[C;z2;x0;D2,D3:2]``
     - internal C=C, no heteroatom on either carbon
   * - ``alkenyl_boronic_acid``
     - ``functional:62``
     - ``[B;D3;z1;x2](-[O;D1])(-[O;D1])-;!@[C;z2;x1,x2:1]=[C:2]``
     - vinylboronic acid
   * - ``alkenyl_boronic_ester``
     - ``functional:63``
     - ``[B;D3;z1;x2](-[O;D2;x1])(-[O;D2;x1])-;!@[C;z2;x1,x2:1]=[C:2]``
     - vinylboronate ester
   * - ``alkenyl_boroxine``
     - ``functional:288``
     - ``[B;D3;z1;x2;r6](-[O;D2;r6]-[B;r6])(-[O;D2;r6]-[B;r6])-;!@[C;z2:1]=[C:2]``
     - (CH2=CHBO)3, trivinylboroxine the common one: a vinyl donor for the Suzuki
   * - ``alkenyl_bromide``
     - ``functional:49``
     - ``[Br;D1][C;z2;x1,x2:1]=[C:2]``
     - vinylic C-Br
   * - ``alkenyl_carboxylic_acid``
     - ``functional:81``
     - ``[O;D1;z1;x0]-[C;z2;x2;D3:1](=[O:2])-[C;z2:3]=[C:4]``
     - alpha,beta-unsaturated acid
   * - ``alkenyl_chloride``
     - ``functional:48``
     - ``[Cl;D1][C;z2;x1,x2:1]=[C:2]``
     - vinylic C-Cl, the Heck/Negishi electrophile
   * - ``alkenyl_fluoride``
     - ``functional:47``
     - ``[F;D1][C;z2;x1:1]=[C:2]``
     - vinylic C-F
   * - ``alkenyl_grignard``
     - ``functional:112``
     - ``[Mg;D2](-[F,Cl,Br,I])-[C;z2:1]=[C:2]``
     - vinyl Grignard
   * - ``alkenyl_halide``
     - ``functional:255``
     - ``[Cl,Br,I;D1]-[C;z2;x1,x2:1]=[C:2]``
     - vinylic C-X, the three ``alkenyl_*`` halides as one class; F excluded as in ``aryl_halide``
   * - ``alkenyl_iodide``
     - ``functional:50``
     - ``[I;D1][C;z2;x1,x2:1]=[C:2]``
     - vinylic C-I
   * - ``alkenyl_molander_salt``
     - ``functional:68``
     - ``[B;D4;z1;x3;-](-[F])(-[F])(-[F])-;!@[C;z2:1]=[C:2]``
     - vinyl trifluoroborate
   * - ``alkenyl_nitrile``
     - ``functional:408``
     - ``[C;z2;x0:1]=[C;z2:2]-[C;D2;z3:3]#[N;D1:4]``
     - C=C-C#N, the unsaturated nitrile Michael acceptor
   * - ``alkenyl_silane``
     - ``functional:123``
     - ``[Si;D4]-;!@[C;z2:1]=[C:2]``
     - vinylsilane
   * - ``alkenyl_stannane``
     - ``functional:120``
     - ``[Sn;D4;z1]-;!@[C;z2:1]=[C:2]``
     - vinylstannane
   * - ``alkenyl_sulfonate``
     - ``functional:256``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;z2;x1,x2:1]=[C:2])-[C]``
     - enol triflate, nonaflate or tosylate, the vinylic pseudohalide
   * - ``alkenyl_sulfonium``
     - ``functional:413``
     - ``[C;z2;x0:1]=[C;z2;D2:2]-[S;+;D3:3]``
     - C=C-S+R2, the alkenyl sulfonium (thianthrenium) acceptor
   * - ``alkenyl_thiol``
     - ``functional:353``
     - ``[S;x0;D1;z1:1]-[C;z2:2]=[C;M]``
     - ene-thiol, the drawn tautomer of a non-aromatic ring thioamide such as a 2-mercaptopyridinone
   * - ``alkenyl_zinc``
     - ``functional:115``
     - ``[Zn;*]-[C;z2:1]=[C:2]``
     - vinylzinc
   * - ``alkoxide``
     - ``functional:232``
     - ``[O;D1;z1;x0;-:1]-[C;z1:2]``
     - RO-
   * - ``alkoxy_ester``
     - ``functional:294``
     - ``[O;z2;x0:1]=[C;D3;x2;z2:2]-[O;D2;x0:3]``
     - ``ester`` with the alkoxy O numbered, for the rows that keep it: partial reduction to the hemiacetal or lactol
   * - ``alkyl_boronic_acid``
     - ``functional:14``
     - ``[B;D3;z1;x2](-[O;D1])(-[O;D1])-;!@[C;z1;x1,x2,x3:1]``
     - RB(OH)2 on sp3 carbon
   * - ``alkyl_boronic_ester``
     - ``functional:61``
     - ``[B;D3;z1;x2](-[O;D2;x1])(-[O;D2;x1])-;!@[C;z1;x1,x2,x3:1]``
     - RB(OR')2 on sp3 carbon
   * - ``alkyl_boroxine``
     - ``functional:260``
     - ``[B;D3;z1;x2;r6](-[O;D2;r6]-[B;r6])(-[O;D2;r6]-[B;r6])-;!@[C;z1:1]``
     - (RBO)3, trimethylboroxine the common one: a methyl donor for the Suzuki
   * - ``alkyl_bromide``
     - ``functional:10``
     - ``[Br;D1][C;z1;x1:1]``
     - R-Br on sp3 carbon
   * - ``alkyl_carbamate``
     - ``functional:358``
     - ``[N;z1,a:1]-[C;D3;z2;x3:2](=[O;D1:3])-[O;D2;x0]-[C;z1]``
     - N-CO-OR on sp3 R, Boc, Cbz and methyl/ethyl carbamates included; the alkoxide leaves only in a ring closure
   * - ``alkyl_carboxylic_acid``
     - ``functional:79``
     - ``[O;D1;z1;x0]-[C;z2;x2;D3:1](=[O:2])-[C;z1:3]``
     - R-COOH on sp3 carbon
   * - ``alkyl_chloride``
     - ``functional:9``
     - ``[Cl;D1][C;z1;x1:1]``
     - R-Cl on sp3 carbon
   * - ``alkyl_fluoride``
     - ``functional:46``
     - ``[F;D1][C;z1;x1:1]``
     - R-F on sp3 carbon; not a leaving group, a metabolic-stability flag
   * - ``alkyl_germatrane``
     - ``functional:283``
     - ``[Ge;D4;R]-;!@[C;z1:1]``
     - an alkyl carbagermatrane, a bench-stable alkyl donor
   * - ``alkyl_grignard``
     - ``functional:110``
     - ``[Mg;D2](-[F,Cl,Br,I])-[C;z1:1]``
     - RMgX on sp3 carbon
   * - ``alkyl_halide``
     - ``functional:45``
     - ``[Cl,Br,I;D1][C;z1;x1:1]``
     - R-X on sp3 carbon, the SN2 electrophile
   * - ``alkyl_hydrazine``
     - ``functional:139``
     - ``[N;D1;z1;x1:1]-[N;D2;z1;x1:2]-[C;z1:3]``
     - R-NH-NH2 on sp3 carbon
   * - ``alkyl_iodide``
     - ``functional:11``
     - ``[I;D1][C;z1;x1:1]``
     - R-I on sp3 carbon
   * - ``alkyl_mesylate``
     - ``functional:59``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;z1;x1:1])-[C;D1]``
     - ROMs
   * - ``alkyl_molander_salt``
     - ``functional:67``
     - ``[B;D4;z1;x3;-](-[F])(-[F])(-[F])-;!@[C;z1:1]``
     - RBF3- on sp3 carbon
   * - ``alkyl_stannane``
     - ``functional:121``
     - ``[Sn;D4;z1]-;!@[C;z1;D2,D3,D4:1]``
     - RSnR3 on sp3 carbon, terminal methyls excluded
   * - ``alkyl_sulfonate``
     - ``functional:243``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;z1;x1:1])-[C]``
     - ROSO2R' on sp3 carbon, any R'.  The SN2 leaving group as one class
   * - ``alkyl_sulfonylhydrazide``
     - ``functional:284``
     - ``[S;D4;x3](=[O])(=[O])(-[C])-[N;D2;z1;x2]-[N;D2;z1;x1]-[C;z1:1]``
     - R-NH-NH-SO2Ar, an alkyl radical source; ``sulfonylhydrazide`` requires the terminal NH2
   * - ``alkyl_tosylate``
     - ``functional:60``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;z1;x1:1])-[C;a]:1:[C;a;D2]:[C;a;D2]:[C;a](-[C;D1]):[C;a;D2]:[C;a;D2]:1``
     - ROTs
   * - ``alkyl_triflate``
     - ``functional:58``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;z1;x1:1])-[C;D4](-[F])(-[F])-[F]``
     - ROTf on sp3 carbon, a strong SN2 leaving group
   * - ``alkyl_zinc``
     - ``functional:113``
     - ``[Zn;*]-[C;z1:1]``
     - RZnX, R2Zn or RZn+, the Negishi nucleophile; ``*`` takes the zinc at any charge
   * - ``alkyne``
     - ``functional:4``
     - ``[C;z3;x0;D2:1]#[C;x0;D2:2]``
     - internal C#C
   * - ``alkynyl_boronic_acid``
     - ``functional:64``
     - ``[B;D3;z1;x2](-[O;D1])(-[O;D1])-;!@[C;z3;x1:1]#[C:2]``
     - alkynylboronic acid
   * - ``alkynyl_boronic_ester``
     - ``functional:65``
     - ``[B;D3;z1;x2](-[O;D2;x1])(-[O;D2;x1])-;!@[C;z3;x1:1]#[C:2]``
     - alkynylboronate ester
   * - ``alkynyl_bromide``
     - ``functional:53``
     - ``[Br;D1][C;z3;x1:1]#[C:2]``
     - C(sp)-Br
   * - ``alkynyl_carboxylic_acid``
     - ``functional:82``
     - ``[O;D1;z1;x0]-[C;z2;x2;D3:1](=[O:2])-[C;z3;x0:3]#[C:4]``
     - propiolic-type acid
   * - ``alkynyl_chloride``
     - ``functional:52``
     - ``[Cl;D1][C;z3;x1:1]#[C:2]``
     - C(sp)-Cl
   * - ``alkynyl_fluoride``
     - ``functional:51``
     - ``[F;D1][C;z3;x1:1]#[C:2]``
     - C(sp)-F; ``z3`` is genuine sp here
   * - ``alkynyl_grignard``
     - ``functional:417``
     - ``[Mg;D2](-[F,Cl,Br,I])-[C;z3:1]#[C:2]``
     - alkynyl Grignard
   * - ``alkynyl_iodide``
     - ``functional:54``
     - ``[I;D1][C;z3;x1:1]#[C:2]``
     - C(sp)-I
   * - ``alkynyl_molander_salt``
     - ``functional:69``
     - ``[B;D4;z1;x3;-](-[F])(-[F])(-[F])-;!@[C;z3;x1:1]#[C:2]``
     - alkynyl trifluoroborate
   * - ``alkynyl_silane``
     - ``functional:124``
     - ``[Si;D4]-;!@[C;z3:1]#[C:2]``
     - TMS-protected alkyne
   * - ``allyl_halide``
     - ``functional:223``
     - ``[Cl,Br,I;D1:1]-[C;z1;x1:2]-[C;z2:3]=[C;z2:4]``
     - X-C-C=C, the allylic SN2/SN2' substrate
   * - ``allylic_methyl``
     - ``functional:345``
     - ``[C;z1;D1;x0:1]-[C;z2;M]=[C;z2;M]``
     - CH3-C=C, the allylic Wohl-Ziegler site
   * - ``alpha_amino_acid``
     - ``functional:292``
     - ``[N;D1;z1;x0:1]-[C;z1;D2,D3:2]-[C;z2;x2:3](=[O;D1:4])-[O;D1:5]``
     - H2N-CHR-COOH, the diazotization substrate that keeps its configuration
   * - ``alpha_ch_acyl``
     - ``functional:342``
     - ``[C;z1;x0;h1,h2:1]-[C;z2;x2;D3;M](=[O;z2;x0;M])-[O,N;z1;D2,D3;M]``
     - a C-H beside an ester, lactone or amide carbonyl, the site of a base-mediated enolate halogenation
   * - ``alpha_ch_carbonyl``
     - ``functional:341``
     - ``[C;z1;x0;h1,h2,h3:1]-[C;z2;x1;D2,D3;M]=[O;z2;x0;M]``
     - an enolizable C-H beside a ketone or aldehyde carbonyl, the site an alpha-halogenation takes
   * - ``alpha_dihalo_ketone``
     - ``functional:336``
     - ``[Cl,Br;D1]-[C;z1;x2;D4:1](-[Cl,Br;D1])-[C;z2;x1:2]=[O;z2;x0:3]``
     - alpha,alpha-dihaloketone, the dichloroketene [2+2] adduct the Zn reduction strips
   * - ``alpha_haloester``
     - ``functional:74``
     - ``[O;D2;x0:1]-[C;D3;x2;z2:2](=[O:3])-[C;z1;D2,D3;x1:4]-[Cl,Br,I;D1]``
     - alpha-halo ester
   * - ``alpha_haloketone``
     - ``functional:73``
     - ``[O;z2;x0:1]=[C;D3;x1:2]-[C;z1;D2,D3;x1:3]-[Cl,Br,I;D1]``
     - alpha-halo ketone, the Hantzsch thiazole electrophile
   * - ``alpha_hetero_alkyl_halide``
     - ``functional:270``
     - ``[Cl,Br,I;D1]-[C;z1;x2:1]-[O;M]-[C,P;M]=[O;M]``
     - X-CH-O-C(=O) or X-CH-O-P(=O): acyloxymethyl, carbonate and phosphate halides; MOM-, SEM- and BOM-Cl stay protecting groups
   * - ``alpha_hetero_alkyl_sulfonate``
     - ``functional:273``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;z1;x2:1]-[N,O,S,P;M])-[C]``
     - RSO2-O-C-Y with Y an N, O, S or P: TsOCH2P(O)(OEt)2; ``alkyl_sulfonate`` (``x1``) leaves it out
   * - ``alpha_heteroatom_halide``
     - ``functional:225``
     - ``[Cl,Br,I;D1:1]-[C;z1;h1,h2:2]-[N,O,S;z1:3]``
     - X-C-Het, the MOM-chloride class: an alkylating agent far more reactive than ``alkyl_halide`` says
   * - ``alpha_ketone``
     - ``functional:72``
     - ``[O;z2;x0:1]=[C;D3;x1:2]-[C;z1;D1,D2;x0:3]``
     - ketone with an enolizable alpha CH2/CH3
   * - ``amide_carbonyl``
     - ``functional:295``
     - ``[O;z2;x0:1]=[C;D2,D3;z2;x2:2]-[N:3]``
     - any amide, lactam or imide C=O, the N kept: aromatic N included, carbamate and urea (``x3``) excluded
   * - ``amidine``
     - ``functional:152``
     - ``[N;D1;z1;x0:1]-[C;D3;z2;x2:2]=[N;D1:3]``
     - RC(=NH)NH2, the pyrimidine partner
   * - ``amidine_nh``
     - ``functional:416``
     - ``[N;D1;z1;x0:1]-[C;D3;z2:2]=[N;D1;x0:3]``
     - H2N-C(=NH)-R for R = C, N or S: the free amidine, guanidine or isothiourea
   * - ``amidoxime``
     - ``functional:151``
     - ``[N;D1;z2;x0:1]=[C;D3;x2:2]-[N;D2;z1;x1:3]-[O;D1:4]``
     - RC(=NH)NHOH, the 1,2,4-oxadiazole precursor
   * - ``amino_alcohol``
     - ``functional:159``
     - ``[N;D1;z1;x0:1]-[C;z1:2]-[C;z1:3]-[O;D1:4]``
     - 1,2-amino alcohol
   * - ``aminopyrazole_nh``
     - ``functional:419``
     - ``[N;D1;z1;x0:1]-[C;a;r5:2]:[N;a;D2;r5:3]:[N;a;D2;r5:4]``
     - NH2 on an azole carbon beside an N-N pair with no N-substituent: 3(5)-aminopyrazole and aminotriazole in either drawing
   * - ``aminopyridine``
     - ``functional:158``
     - ``[N;D1;z1;x0:1]-[C;a:2]:[N;a;h0;D2:3]``
     - 2-aminoazine, the GBB partner
   * - ``ammonia``
     - ``functional:350``
     - ``[N;D0;h3:1]``
     - NH3, neutral; an ammonium salt drawn as NH3 plus an acid matches, a drawn [NH4+] does not
   * - ``ammonium``
     - ``functional:235``
     - ``[N;z1;x0;+;h1,h2,h3:1]``
     - The protonated amine.  ``h1,h2,h3`` is what distinguishes it from a quaternary salt
   * - ``anhydride``
     - ``functional:96``
     - ``[C;z2;x2;D3:1](=[O:2])-[O;D2;x0]-[C;z2;x2;D3:3]=[O:4]``
     - R-CO-O-CO-R; also matches ``ester`` twice
   * - ``aniline_ortho_ch``
     - ``functional:168``
     - ``[N;D1;z1;x0:1]-[C;a:2]:[C;a;D2:3]``
     - Ar-NH2 with a free ortho CH
   * - ``anthranilic_acid``
     - ``functional:150``
     - ``[N;D1;z1;x0:1]-[C;a:2]:[C;a:3]-[C;z2;x2;D3:4](=[O:5])-[O;D1]``
     - 2-aminobenzoic acid
   * - ``arene_ch``
     - ``functional:41``
     - ``[C;a;D2:1]``
     - an aromatic CH, the electrophilic-substitution site
   * - ``aryl_boronic_acid``
     - ``functional:12``
     - ``[B;D3;z1;x2](-[O;D1])(-[O;D1])-;!@[C;a:1]``
     - ArB(OH)2
   * - ``aryl_boronic_ester``
     - ``functional:13``
     - ``[B;D3;z1;x2](-[O;D2;x1])(-[O;D2;x1])-;!@[C;a:1]``
     - ArB(OR)2, pinacol boronate and its kin
   * - ``aryl_boronic_hemiester``
     - ``functional:285``
     - ``[B;D3;z1;x2](-[O;D1])(-[O;D2;x1])-;!@[C;a:1]``
     - ArB(OH)(OR), a boronic hemiester or ring-opened pinacol boronate
   * - ``aryl_bromide``
     - ``functional:7``
     - ``[Br;D1]-[C;a:1]``
     - Ar-Br
   * - ``aryl_bromide_iodide``
     - ``functional:44``
     - ``[Br,I;D1]-[C;a:1]``
     - Ar-Br or Ar-I, the halides that couple without a specialist ligand
   * - ``aryl_carbamate``
     - ``functional:356``
     - ``[N;z1:1]-[C;D3;z2;x3:2](=[O;D1:3])-[O;D2;x0]-[C;a]``
     - N-CO-OAr -- phenyl, 4-nitrophenyl, 4-fluorophenyl: the aryloxide leaves and an amine takes the carbonyl
   * - ``aryl_carbonate``
     - ``functional:89``
     - ``[O;D2:1]-[C;z2;x3;D3:2](=[O:3])-[O;D2;x0]-[C;a]``
     - aryl carbonate; the aryloxide is the leaving group
   * - ``aryl_carboxylic_acid``
     - ``functional:80``
     - ``[O;D1;z1;x0]-[C;z2;x2;D3:1](=[O:2])-[C;a:3]``
     - Ar-COOH
   * - ``aryl_chloride``
     - ``functional:6``
     - ``[Cl;D1]-[C;a:1]``
     - Ar-Cl
   * - ``aryl_dialkylborane``
     - ``functional:286``
     - ``[B;D3;z1;x0](-[C;z1])(-[C;z1])-;!@[C;a:1]``
     - ArBR2, diethyl(3-pyridyl)borane and its kin
   * - ``aryl_difluoroborane``
     - ``functional:287``
     - ``[B;D3](-[F;D1])(-[F;D1])-;!@[C;a:1]``
     - ArBF2, a trifluoroborate drawn without its fluoride and cation
   * - ``aryl_ester``
     - ``functional:307``
     - ``[O;z2;x0:1]=[C;D3;x2;z2:2]-[O;D2;x0]-[C;a]``
     - R-CO-OAr -- pentafluorophenyl, 4-nitrophenyl, phenyl -- the acylating agent; ``ester`` also matches
   * - ``aryl_ether``
     - ``functional:185``
     - ``[O;D2;z1;x0:1](-[C;a:2])[C;z1,a:3]``
     - Ar-O-R and Ar-O-Ar.  Separate from ``dialkyl_ether`` because only this one is an SNAr and Ullmann product and a demethylation substrate
   * - ``aryl_fluoride``
     - ``functional:5``
     - ``[F;D1]-[C;a:1]``
     - Ar-F, an SNAr electrophile rather than a cross-coupling one
   * - ``aryl_grignard``
     - ``functional:111``
     - ``[Mg;D2](-[F,Cl,Br,I])-[C;a:1]``
     - ArMgX
   * - ``aryl_halide``
     - ``functional:43``
     - ``[Cl,Br,I;D1]-[C;a:1]``
     - Ar-X, the cross-coupling electrophile; F excluded, it does SNAr instead
   * - ``aryl_hydrazine``
     - ``functional:140``
     - ``[N;D1;z1;x1:1]-[N;D2;z1;x1:2]-[C;a:3]``
     - Ar-NH-NH2, no constraint on the ring beyond the attachment carbon
   * - ``aryl_hydrazine_ortho_ch``
     - ``functional:141``
     - ``[N;D1;z1;x1:1]-[N;D2;z1;x1:2]-[C;a:3]:[C;a;D2:4]``
     - Ar-NH-NH2 with a free ortho CH, which an indolization consumes
   * - ``aryl_iodide``
     - ``functional:8``
     - ``[I;D1]-[C;a:1]``
     - Ar-I
   * - ``aryl_magnesium_cation``
     - ``functional:289``
     - ``[Mg;D1;+]-[C;a:1]``
     - Ar-Mg+, an aryl Grignard drawn without its halide
   * - ``aryl_mesylate``
     - ``functional:56``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;a:1])-[C;D1]``
     - ArOMs
   * - ``aryl_methyl``
     - ``functional:344``
     - ``[C;z1;D1;x0:1]-[C;a;M]``
     - Ar-CH3, the methyl ``benzylic_ch`` (``D2,D3``) leaves out
   * - ``aryl_molander_salt``
     - ``functional:66``
     - ``[B;D4;z1;x3;-](-[F])(-[F])(-[F])-;!@[C;a:1]``
     - ArBF3-, the air-stable trifluoroborate
   * - ``aryl_silane``
     - ``functional:122``
     - ``[Si;D4]-;!@[C;a:1]``
     - ArSiR3, the Hiyama nucleophile
   * - ``aryl_stannane``
     - ``functional:119``
     - ``[Sn;D4;z1]-;!@[C;a:1]``
     - ArSnR3, the Stille nucleophile
   * - ``aryl_sulfonate``
     - ``functional:242``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;a:1])-[C]``
     - ArOSO2R for any R -- triflate, mesylate, tosylate, besylate.  The three named rows answer "which pseudohalide"; this one answers "does this couple", which is the question a coupling row asks
   * - ``aryl_sulfoxide``
     - ``functional:269``
     - ``[S;D3;z2;!R](=[O;D1])(-[C;z1])-[C;a:1]``
     - Ar-S(O)R, a sulfinyl leaving group on an arene
   * - ``aryl_thioether``
     - ``functional:267``
     - ``[S;D2;z1;x0;!R](-[C;z1])-[C;a:1]``
     - Ar-SR, the C-S electrophile of the Liebeskind-Srogl and thioether Negishi couplings
   * - ``aryl_thiol``
     - ``functional:129``
     - ``[S;x0;D1;z1:1]-[C;a:2]``
     - Ar-SH
   * - ``aryl_tosylate``
     - ``functional:57``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;a:1])-[C;a]:1:[C;a;D2]:[C;a;D2]:[C;a](-[C;D1]):[C;a;D2]:[C;a;D2]:1``
     - ArOTs
   * - ``aryl_triflate``
     - ``functional:55``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;a:1])-[C;D4](-[F])(-[F])-[F]``
     - ArOTf, the pseudohalide that couples like an aryl bromide
   * - ``aryl_zinc``
     - ``functional:114``
     - ``[Zn;*]-[C;a:1]``
     - ArZnX, Ar2Zn or ArZn+
   * - ``aryloxymethyl_halide``
     - ``functional:380``
     - ``[Cl,Br,I;D1]-[C;z1;x2;D2:1]-[O;D2;M]-[C;a;M]``
     - ArO-CH2-X, the aryloxymethylating agent
   * - ``aza_ortho_chloride``
     - ``functional:338``
     - ``[Cl;D1]-[C;a;r6:1]:[N;a;D2;h0;r6:2]``
     - 2-chloroazine with the ring N numbered, for the row that puts the hydrogen on it; ``activated_aza_ortho_halide`` masks that N
   * - ``aza_ortho_thioether``
     - ``functional:354``
     - ``[S;D2;z1;x0;!R](-[C;z1])-[C;a:1]:[N;a;D2;h0;M]``
     - alkylthio on an aromatic carbon beside an aza N, the leaving group of a chlorinolysis
   * - ``aza_para_chloride``
     - ``functional:339``
     - ``[Cl;D1]-[C;a;r6:1]1:[C,N;a:2]:[C,N;a:3]:[N;a;D2;h0:4]:[C,N;a:5]:[C,N;a:6]:1``
     - 4-chloroazine with the whole ring numbered, for the row that rewrites its bonds to the 4-azinone
   * - ``azetidine``
     - ``functional:238``
     - ``[N;D2;z1;x0;r4;h1:1]([C;z1;r4:2])[C;z1;r4:3]``
     - The four-membered NH amine
   * - ``azide``
     - ``functional:102``
     - ``[N;x1;D2:1]=[N;+:2]=[N;-:3]``
     - R-N3
   * - ``azine_ch_ortho``
     - ``functional:320``
     - ``[C;a;D2;r6:1]:[N;a;D2;h0;r6;M]``
     - a six-ring aromatic C-H beside a pyridine-type N
   * - ``azine_ch_para``
     - ``functional:321``
     - ``[C;a;D2;r6:1]:[C;a;r6;M]:[C;a;r6;M]:[N;a;D2;h0;r6;M]``
     - a six-ring aromatic C-H para to a pyridine-type N
   * - ``azine_dione``
     - ``functional:334``
     - ``[O;z2;x0]=[C;z2;x2;D3;r6:1]1-[N;z1;D2;r6:2]-[C;z2;x3;D3;r6:3](=[O;z2;x0])-[N;z1;D2;r6:4]-[C;z2,z4;M]=,:[C;z2,z4;M]1``
     - uracil class, 2,4-dione with both N-H: the double deoxychlorination substrate
   * - ``azine_n_oxide``
     - ``functional:206``
     - ``[N;a;D3;+:1]-[O;D1;z1;-:2]``
     - The pyridine N-oxide, an activated ring for C-H functionalization.  The oxygen is NUMBERED: a deoxygenation omits it from its own product side, and leaving it unnumbered would delete it from every other template too
   * - ``azine_n_oxide_ch``
     - ``functional:319``
     - ``[C;a;D2:1]:[N;a;D3;+;M]-[O;D1;-;M]``
     - the C-H beside the nitrogen of an azine N-oxide
   * - ``azinone``
     - ``functional:132``
     - ``[O;z2;x0:1]=[C;z2;x2;D3;r6:2](-[C,N;z2,z4;r6:3])-[N;z1;D2;r6:4]``
     - aromatizable cyclic amide; the sp2/aromatic flank excludes a saturated lactam
   * - ``azinone_aza``
     - ``functional:332``
     - ``[O;z2;x0:1]=[C;z2;x3;D3;r6:2](-[N;z2,z4;r6:3])-[N;z1;D2;r6:4]``
     - 2-pyrimidinone class: a ring urea carbonyl flanked by an imine N and an N-H; ``azinone`` (``x2``) leaves it out
   * - ``aziridine_nh``
     - ``functional:94``
     - ``[N;D2;z1;x0;h1;r3:1]``
     - aziridine N-H, matched alone so a ring opening leaves the ring intact
   * - ``azo``
     - ``functional:105``
     - ``[N;D2;x1:1]=[N;D2;x1:2]``
     - R-N=N-R
   * - ``azole_ch``
     - ``functional:318``
     - ``[C;a;D2;r5:1]:[N,O,S;a;M]``
     - an aromatic five-ring C-H beside a ring heteroatom
   * - ``azolo_pyridine``
     - ``functional:347``
     - ``[N;a;D3;r5;r6:1]1:[C;a;R1;r6:2]:[C;a;R1;r6:3]:[C;a;R1;r6:4]:[C;a;R1;r6:5]:[C;a;D3;r5;r6:6]:1``
     - the azine ring of an azole fused through its bridgehead N: imidazo-, pyrazolo- and triazolo[1,5-a]pyridine
   * - ``azolone``
     - ``functional:374``
     - ``[O;z2;x0:1]=[C;z2;x2;D3;r5:2]1-[C;z2;r5:3]=[C;z2;r5;M]-[N,O,S;r5;M]-[N;z1;D2;h1;r5:4]1``
     - five-ring lactam whose O-substituted form is aromatic: isoxazol-3-one, pyrazol-3-one; ``azinone`` is the six-ring class
   * - ``azolone_aza``
     - ``functional:375``
     - ``[O;z2;x0:1]=[C;z2;x3;D3;r5:2](-[N;z2;r5:3])-[N;z1;D2;h1;r5:4]``
     - five-ring cyclic urea flanked by an imine N and an N-H, the 1,2,4-triazol-3-one; ``azinone_aza`` is the six-ring class
   * - ``azolone_vinylogous``
     - ``functional:388``
     - ``[O;z2;x0:1]=[C;z2;x2;D3;r5:2]1-[C;z2;r5:3]=[C;z2;r5:4]-[N;z1;D2;h1;r5:5]-[N,O,S;r5;M]1``
     - five-ring lactam with C=O and N-H across a ring C=C, the 1-substituted pyrazol-5-one; its O-substituted form is the aromatic 5-alkoxypyrazole
   * - ``benzo_pyridine_ring``
     - ``functional:406``
     - ``[N;a;D2;h0;R1:1]1:[C;a;R1:2]:[C;a;R1:3]:[C;a;R1:4]:[C;a;R2:5]:[C;a;R2:6]:1``
     - the pyridine ring of a quinoline, 7-azaindole or naphthyridine, nitrogen beside the fusion carbon
   * - ``benzyl_halide``
     - ``functional:222``
     - ``[Cl,Br,I;D1:1]-[C;z1;x1;h1,h2:2]-[C;a:3]``
     - Ar-CH2-X.  ``alkyl_halide`` already matches it; this says the SN2 is a benzylic one, which is a different rate and a different selectivity
   * - ``benzylic_alcohol``
     - ``functional:397``
     - ``[O;D1;z1;x0:1]-[C;z1;x1:2]-[C;a,z3;M]``
     - an alcohol on a carbon beside an arene or an alkyne, the stabilized carbocation silane reduction needs
   * - ``benzylic_ch``
     - ``functional:42``
     - ``[C;z1;h1,h2;D2,D3:1]-[C;a:2]``
     - an sp3 carbon bearing hydrogen next to a ring.  ``D2,D3`` excludes a TOLUENE METHYL, which is a D1 -- so ``benzylic_bromination`` does not reach one
   * - ``beta_arylethylamine``
     - ``functional:157``
     - ``[N;D1;z1;x0:1]-[C;z1:2]-[C;z1:3]-[C;a:4]:[C;a;D2:5]``
     - the Pictet-Spengler substrate
   * - ``beta_carbonyl_acid``
     - ``functional:403``
     - ``[O;D1;z1;x0]-[C;z2;x2;D3](=[O;z2;x0])-[C;z1:1]-[C,S,N;z2,z3,z5;M]=,#[O,N;M]``
     - a carboxylic acid beside a second carbonyl, nitrile or sulfonyl: malonic, beta-keto and cyanoacetic acids
   * - ``beta_carbonyl_ester``
     - ``functional:402``
     - ``[O;z2;x0]=[C;z2;x2;D3](-[O;D2;x0]-[C;z1;D1,D2,D3])-[C;z1:1]-[C,S,N;z2,z3,z5;M]=,#[O,N;M]``
     - an alkyl ester on a carbon beside a second carbonyl, nitrile or sulfonyl, the ester unnumbered so it leaves whole; a tert-butyl ester is left to ``carboxyl_tbu``
   * - ``beta_dicarbonyl_ketone``
     - ``functional:392``
     - ``[O;z2;x0:1]=[C;D3;x1;z2:2]-[C;z1;h1,h2:3]-[C;z2;x1,x2;M]=[O;z2;x0;M]``
     - a ketone whose alpha C-H is flanked by a second carbonyl (1,3-dione, beta-keto ester or amide), enolized by an amine base
   * - ``beta_ketoester``
     - ``functional:78``
     - ``[O;z2;x0:1]=[C;D3;x1:2]-[C;z1;D2;x0:3]-[C;z2;x2;D3:4](=[O:5])-[O;D2]``
     - beta-keto ester
   * - ``beta_lactam``
     - ``functional:239``
     - ``[N;D2,D3;z1;r4:1]-[C;D3;z2;r4:2]=[O;D1;x0:3]``
     - The four-membered lactam.  ``D2,D3`` because a fused penicillin or cephalosporin nitrogen is D3 and ``D2`` alone reaches only the monocyclic model compound
   * - ``biaryl_aniline``
     - ``functional:95``
     - ``[N;D2;z1;x0:1](-[C;a:2])-[C;a:3]``
     - Ar-NH-Ar
   * - ``boronate_alkyl_bromide``
     - ``functional:117``
     - ``[Br;D1]-[C;z1;D2;x2:1]-[B:2]``
     - B-CH2-Br
   * - ``boronate_alkyl_chloride``
     - ``functional:116``
     - ``[Cl;D1]-[C;z1;D2;x2:1]-[B:2]``
     - B-CH2-Cl, the one-carbon SN2 electrophile
   * - ``boronate_alkyl_iodide``
     - ``functional:118``
     - ``[I;D1]-[C;z1;D2;x2:1]-[B:2]``
     - B-CH2-I
   * - ``bridged_diaryl``
     - ``functional:181``
     - ``[C,N;a:1]-;!@[C,O,N,S:2]-;!@[C,N;a:3]``
     - two rings joined through one bridging atom
   * - ``bt_sulfone``
     - ``functional:291``
     - ``[C;z1:1]-[S;D4](=[O])(=[O])-[C;a]1:[N;a]:[C;a]:[C;a]:[S;a]:1``
     - an alkyl benzothiazolyl sulfone, the Julia-Kocienski reagent
   * - ``carbamate``
     - ``functional:191``
     - ``[N:1]-[C;D3;z2;x3:2](=[O;D1;z2;x0:3])-[O;D2;z1:4]-[C:5]``
     - N-CO-O-R, the Boc/Cbz/Fmoc backbone.  Reported wherever it occurs, unlike ``protective.tsv``, which only names the ones it can cleave
   * - ``carbamate_alpha_ch``
     - ``functional:322``
     - ``[C;z1;h1,h2,h3:1]-[N;D3;z1;M]-[C;z2;M](=[O;M])-[O;D2;M]``
     - a C-H on the nitrogen of an N,N-disubstituted carbamate, the HAT radical source
   * - ``carbamoyl_chloride``
     - ``functional:91``
     - ``[Cl;D1]-[C;z2;x3;D3:1](=[O:2])-[N;D2,D3:3]``
     - R2N-CO-Cl
   * - ``carbamoyl_fluoride``
     - ``functional:92``
     - ``[F;D1]-[C;z2;x3;D3:1](=[O:2])-[N;D2,D3:3]``
     - R2N-CO-F
   * - ``carbazate``
     - ``functional:371``
     - ``[N;D1;z1;x1:1]-[N;D2;z1;x1:2]-[C;D3;z2;x3:3](=[O;D1;M])-[O;D2;M]``
     - H2N-NH-CO2R, Boc- or ethoxycarbonyl hydrazine; the carbamate stays on the N it was on
   * - ``carbodiimide``
     - ``functional:217``
     - ``[N;D2;z2:1]=[C;D2;z5:2]=[N;D2;z2:3]``
     - N=C=N, the DCC/EDC coupling reagent.  ``z5`` is the cumulated-double centre; chython 2 wrote ``z3`` here and V3's ``z3`` is sp and nothing else
   * - ``carbonate``
     - ``functional:193``
     - ``[C:1]-[O;D2;z1:2]-[C;D3;z2;x3:3](=[O;D1;z2;x0:4])-[O;D2;z1:5]-[C:6]``
     - R-O-CO-O-R, both oxygens esterified
   * - ``carboxylate``
     - ``functional:230``
     - ``[O;D1;z1;x0;-:1]-[C;D3;z2:2]=[O;D1;z2;x0:3]``
     - The deprotonated acid.  ``carboxylic_acid`` cannot match it, its own oxygen being neutral by omission, so a registration salt reports no acid at all without this row
   * - ``carboxylic_acid``
     - ``functional:22``
     - ``[O;D1;z1;x0:3][C;z2;x2;D3:1]=[O:2]``
     - R-COOH; the hydroxyl is the leaving group
   * - ``catechol``
     - ``functional:174``
     - ``[O;D1;z1;x0:1]-[C;a:2]:[C;a:3]-[O;D1;z1;x0:4]``
     - 1,2-dihydroxyarene
   * - ``cdi``
     - ``functional:383``
     - ``[O;D1;z2:1]=[C;D3;x3:2](-[N;a;r5:3]1:[C;a;M]:[N;a;M]:[C;a;M]:[C;a;M]:1)-[N;a;r5]1:[C;a]:[N;a]:[C;a]:[C;a]:1``
     - carbonyldiimidazole as the reagent that keeps one imidazole: the carbamoyl imidazole
   * - ``chloroazine``
     - ``functional:133``
     - ``[Cl;D1]-[C;a:1]:[N;a;D2:2]``
     - Cl on an aromatic carbon next to a ring N; hydrolyses to the azinone
   * - ``chloroformate``
     - ``functional:86``
     - ``[Cl;D1]-[C;z2;x3;D3:1](=[O:2])-[O;D2:3]``
     - RO-CO-Cl, the carbamoylation reagent
   * - ``chlorosilane``
     - ``functional:229``
     - ``[F,Cl,Br,I;D1:1]-[Si;D4:2]``
     - Si-X, the silylating agent
   * - ``chlorothionoformate``
     - ``functional:276``
     - ``[Cl;D1]-[C;D3;z2:1](=[S:2])-[O;D2:3]``
     - Cl-C(=S)-OR, the O-thiocarbamate and thionocarbonate former
   * - ``curtius_acid``
     - ``functional:346``
     - ``[C:1]-[C;z2;x2;D3:2](=[O;D1;z2:3])-[O;D1;z1]``
     - R-COOH with R numbered, for the Curtius row that moves R from the carbonyl carbon to the new nitrogen
   * - ``cyanogen_halide``
     - ``functional:365``
     - ``[Cl,Br;D1]-[C;D2:1]#[N;D1:2]``
     - BrCN or ClCN, the 2-aminoazole carbon
   * - ``cyclic_carbamate``
     - ``functional:192``
     - ``[N;D2,D3;z1;r5,r6:1]-[C;D3;z2;x3;r5,r6:2](=[O;D1;z2:3])-[O;D2;z1;r5,r6:4]``
     - The oxazolidinone/cyclic-carbamate ring, which is a scaffold rather than a protected amine
   * - ``cyclic_carboxylic_acid``
     - ``functional:83``
     - ``[O;D1;z1;x0]-[C;z2;x2;D3:1](=[O:2])-[C;z1;r5,r6:3]``
     - COOH on a saturated 5- or 6-ring
   * - ``dialkyl_ether``
     - ``functional:184``
     - ``[O;D2;z1;x0:1]([C;z1:2])[C;z1:3]``
     - R-O-R, both sp3.  ``x0`` and both carbons stated: the bare ``[O;D2;z1]`` is also every ester O, every anisole and every silyl ether
   * - ``diaryl``
     - ``functional:180``
     - ``[C,N;a:1]-;!@[C,N;a:2]``
     - two rings joined by an acyclic single bond
   * - ``diazo``
     - ``functional:103``
     - ``[C;z2:1]=[N;+;D2:2]=[N;-;D1:3]``
     - R2C=N2
   * - ``diazonium``
     - ``functional:104``
     - ``[C;a:1]-[N;+;D2:2]#[N;D1:3]``
     - ArN2+
   * - ``diboron_ester``
     - ``functional:258``
     - ``[B;D3:1](-[O;D2:2])(-[O;D2:3])-[B;D3](-[O;D2])-[O;D2]``
     - B2pin2 and its kin, the Miyaura boron source; one boron and its two oxygens are kept
   * - ``difluoromethoxy``
     - ``functional:208``
     - ``[F:1][C;D3;z1;x3:2]([F:3])-[O;D2:4]``
     - OCHF2
   * - ``difluoromethyl``
     - ``functional:173``
     - ``[C;D3;z1;x2:1](-[F:2])-[F:3]``
     - -CHF2
   * - ``difluoromethyl_sulfonate``
     - ``functional:274``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;z1;x3;D3:1](-[F;M])-[F;M])-[C]``
     - TfO-CHF2, the difluoromethylating sulfonate
   * - ``dimethyl_sulfate``
     - ``functional:330``
     - ``[C;D1:1]-[O;D2]-[S;D4](=[O])(=[O])-[O;D2]-[C;D1]``
     - (MeO)2SO2, the methylating agent; one methyl transfers
   * - ``disulfide``
     - ``functional:130``
     - ``[S;D2;z1:1]-[S;D2;z1:2]``
     - R-S-S-R
   * - ``dmf_dma``
     - ``functional:372``
     - ``[C;D3;x3:1](-[O;D2]-[C;D1])(-[O;D2]-[C;D1])-[N;D3](-[C;D1])-[C;D1]``
     - Me2N-CH(OMe)2 (DMF-DMA): the one-carbon unit that turns a methyl or methylene ketone into its enaminone; NMe2 and both methoxides leave
   * - ``enal``
     - ``functional:71``
     - ``[O;z2;x0:1]=[C;D2;x1;z2:2]-[C;z2;x0:3]=[C;x0:4]``
     - alpha,beta-unsaturated aldehyde
   * - ``enamine``
     - ``functional:182``
     - ``[C;z2;D2,D3:1]=[C;z2;D2,D3:2]-[N;z1;D2,D3:3]``
     - C=C-N, the hydrogenation and Stork substrate
   * - ``enol``
     - ``functional:188``
     - ``[O;D1;z1;x0;h1:1]-[C;z2:2]=[C;z2:3]``
     - HO-C=C.  ``z2`` on both carbons, so a phenol -- aromatic, not ``z2`` -- is not one
   * - ``enol_ether``
     - ``functional:183``
     - ``[C;z2;D2,D3:1]=[C;z2;D2,D3:2]-[O;D2:3]``
     - C=C-O-R
   * - ``enolizable_ch``
     - ``functional:308``
     - ``[C;z1;h1,h2,h3:1]-[C;z2,z3;M]=,#[O,N;M]``
     - an sp3 C-H beside a C=O, C=N or C#N carbon, the enolate nucleophile; the activating group is context
   * - ``enolizable_ch2``
     - ``functional:310``
     - ``[C;z1;h2,h3:1]-[C;z2,z3;M]=,#[O,N;M]``
     - ``enolizable_ch`` with two hydrogens to give: the gem-dialkylation, condensation and olefination substrate
   * - ``enolizable_ketone``
     - ``functional:373``
     - ``[O;D1;z2]=[C;z2;x1;D3:1]-[C;z1;D1,D2:2]``
     - a ketone with an alpha CH2 or CH3, the site DMF-DMA formylates
   * - ``enone``
     - ``functional:194``
     - ``[C;z2:1]=[C;z2:2]-[C;D3;z2;x1:3]=[O;D1;x0:4]``
     - C=C-C=O, the Michael acceptor.  The conjugation is the group: neither ``ketone`` nor ``alkene`` alone says a molecule is electrophilic at the beta carbon
   * - ``enoyl``
     - ``functional:407``
     - ``[C;z2;x0:1]=[C;z2:2]-[C;D3;z2;x1,x2:3]=[O;D1;x0:4]``
     - C=C-C(=O)R for R = C, O or N -- enone, unsaturated ester, acid, amide, lactone -- numbered from a beta carbon with no heteroatom
   * - ``epoxide``
     - ``functional:39``
     - ``[O;D2;r3:1]1-[C;r3:2]-[C;r3:3]-1``
     - oxirane
   * - ``ester``
     - ``functional:24``
     - ``[O;z2;x0:1]=[C;D3;x2;z2:2]-[O;D2;x0]``
     - R-CO-OR'; the alkoxy is the leaving group
   * - ``fluoroformate``
     - ``functional:87``
     - ``[F;D1]-[C;z2;x3;D3:1](=[O:2])-[O;D2:3]``
     - RO-CO-F
   * - ``formaldehyde``
     - ``functional:315``
     - ``[O;z2;x0:1]=[C;D1:2]``
     - CH2=O
   * - ``formamide``
     - ``functional:323``
     - ``[N;D1;z1;x0:1]-[C;D2;x2;z2;h1:2]=[O:3]``
     - H2N-CHO, the carbamoyl radical source
   * - ``formate``
     - ``functional:378``
     - ``[O;D1;z1;x0;-:1]-[C;D2;z2:2]=[O;D1;z2;x0:3]``
     - HCOO-, which ``carboxylate`` (``D3``) leaves out
   * - ``furan_o``
     - ``functional:204``
     - ``[O;a;D2:1]``
     - The aromatic divalent oxygen of a furan, oxazole or isoxazole
   * - ``fused_aromatic``
     - ``functional:178``
     - ``[A;a;D3:1](:[A;a:2])(:[A;a:3]):[A;a:4]``
     - a ring-fusion atom: three aromatic neighbours
   * - ``gem_dihalide``
     - ``functional:224``
     - ``[Cl,Br,I;D1:1]-[C;z1:2]-[Cl,Br,I;D1:3]``
     - Two halides on one carbon.  Fluorine is excluded because a gem-difluoride is ``difluoromethyl`` and not an alkylating agent
   * - ``glycosyl_halide``
     - ``functional:296``
     - ``[Cl,Br,I;D1]-[C;z1;x2;R:1]-[O;D2;R;M]``
     - X on the ring carbon beside the ring O: the glycosyl halide, the Koenigs-Knorr donor; an acyclic X-C-O is a protecting-group reagent
   * - ``guanidine``
     - ``functional:107``
     - ``[N;z1;x0:1]-[C;!R:2](-[N;z1;x0:3])=[N;x0:4]``
     - acyclic guanidine
   * - ``halodifluoromethyl``
     - ``functional:272``
     - ``[Br,I;D1]-[C;z1;x3,x4:1](-[F;M])-[F;M]``
     - Br-CF2-R and I-CF2-R, the difluoroalkylating electrophiles
   * - ``hemiketal``
     - ``functional:187``
     - ``[O;D1;z1;x0;h1:1]-[C;z1;x2:2]-[O;D2;z1:3]``
     - HO-C(-OR), the ring-opened sugar form
   * - ``hetero_methyl_halide``
     - ``functional:379``
     - ``[Cl,Br,I;D1]-[C;z1;x2;D2:1]-[S,P,F;M]``
     - X-CH2 on S, P or F: a halomethyl sulfide, sulfone, sulfonamide or phosphonate, or FCH2X.  ``alkyl_halide`` (``x1``) leaves it out; the O-bonded halomethyls are ``alpha_hetero_alkyl_halide`` or protecting groups
   * - ``heteroaryl_ol``
     - ``functional:293``
     - ``[O;D1;z1;x0;h1:1]-[C;a:2]:[N,O,S;a;M]``
     - an OH on an aromatic carbon beside a ring heteroatom, the deoxychlorination substrate ``phenol`` does not distinguish
   * - ``heteroatom_alkene``
     - ``functional:278``
     - ``[C;z2;x1,x2:1]=[C;z2:2]``
     - C=C with a heteroatom on at least one carbon: enol ether, enamine, nitroalkene, vinyl fluoride; ``alkene`` (``x0``) leaves it out
   * - ``heteroatom_terminal_alkene``
     - ``functional:277``
     - ``[C;z2;x0;D1:1]=[C;z2;x1;D2,D3:2]``
     - CH2=C with one heteroatom on the inner carbon: vinyl ether, enamide, vinyl sulfonamide, 2-fluoroacrylate, vinyl boronate
   * - ``hexamethyldisilazane``
     - ``functional:351``
     - ``[N;D2;x2;*:1](-[Si;D4](-[C;D1])(-[C;D1])-[C;D1])-[Si;D4](-[C;D1])(-[C;D1])-[C;D1]``
     - (Me3Si)2NH or its anion (LiHMDS): an ammonia equivalent, both silyls leave
   * - ``hindered_aryloxyl``
     - ``functional:251``
     - ``[O;D1;h0;z1:1]-[C;a:2](:[C;a:3]-[C;D4;z1:4]):[C;a:5]-[C;D4;z1:6] |^1:0|``
     - Ar-O• with a quaternary carbon at BOTH ortho positions, which is what makes it isolable; the unhindered phenoxyl radical is a decoy
   * - ``hydrazide``
     - ``functional:161``
     - ``[N;D1;z1;x1:1]-[N;D2;z1;x1:2]-[C;z2;x2:3]=[O:4]``
     - R-CO-NH-NH2
   * - ``hydrazine``
     - ``functional:215``
     - ``[N;z1;x1:1]-[N;z1;x1:2]``
     - N-N single bond, both sp3.  A hydrazide matches this too and ``hydrazide`` names it more precisely
   * - ``hydrazine_nh2nh2``
     - ``functional:367``
     - ``[N;D1;z1:1]-[N;D1;z1;x1:2]``
     - H2N-NH2 itself.  ``x1`` on :2 alone holds for either N and breaks the pattern symmetry, so both N-H tautomers of a ring it closes are enumerated
   * - ``hydrazone``
     - ``functional:142``
     - ``[C;z2:1]=[N;D2;z2;x1:2]-[N;D2;z1;x1:3]``
     - C=N-NH-R
   * - ``hydrazyl``
     - ``functional:252``
     - ``[N;D2;h0;z1:1](-[C;a:2])-[N;D3;z1:3](-[C;a:4])-[C;a:5] |^1:0|``
     - Ar2N-N•-Ar, the DPPH class: three aryl rings to delocalise over, and the trialkyl hydrazyl with none of them is a decoy
   * - ``hydrolysable_amide_n``
     - ``functional:349``
     - ``[N;D3;z1,z4;x0:1]-[C;z2;x2;D3:2](=[O;z2;x0:3])-[C;D2,D3,D4;x0,x1,x2;M]``
     - R-CO-NR'2 and the N-acyl azole, R not methyl or CF3; the N-substituted lactam included
   * - ``hydrolysable_amide_nh``
     - ``functional:348``
     - ``[N;D2;z1;x0:1]-[C;z2;x2;D3:2](=[O;z2;x0:3])-[C;D2,D3,D4;x0,x1,x2;M]``
     - R-CO-NH-R' with R not methyl or CF3: the hydrolysis substrate the named ``amine_acyl`` and ``amine_tfa`` deprotections leave to ``protective.tsv``
   * - ``hydroxamic_acid``
     - ``functional:160``
     - ``[O;D1;z1;x1:1]-[N;D2;z1;x1:2]-[C;z2;x2:3]=[O:4]``
     - R-CO-NH-OH
   * - ``hydroxy_azinone``
     - ``functional:337``
     - ``[O;D1;z1;x0]-[C;z2;r6:1]=[C;z2;r6:2]-[C;z2;x2;D3;r6:3](=[O;z2;x0])-[N;z1;D2;r6:4]``
     - 4-hydroxy-2-pyridone class: an enol conjugated to the lactam carbonyl, the double deoxychlorination substrate
   * - ``hydroxylamine``
     - ``functional:316``
     - ``[N;D1;z1;x1:1]-[O;D1;z1;x1:2]``
     - H2N-OH
   * - ``imidate``
     - ``functional:364``
     - ``[N;D1;z2;x0]=[C;D3;z2;x2:1]-[O;D2]-[C;z1]``
     - R-C(=NH)-OR', the Pinner salt; ammonia and the alkoxide leave
   * - ``imidazole``
     - ``functional:136``
     - ``[N;h1;D2;a;r5:1]:[A:2]:[N;h0;D2;r5:3]``
     - N-H and N separated by one ring atom
   * - ``imidazolyl_carbonate``
     - ``functional:90``
     - ``[O;D2:1]-[C;z2;x3;D3:2](=[O:3])-[N;a;D3;r5]1:[C;a]:[N;a]:[C;a]:[C;a]:1``
     - CDI adduct of an alcohol; here the azole nitrogen leaves
   * - ``imidoyl_halide``
     - ``functional:266``
     - ``[Cl,Br,I;D1]-[C;z2:1]=[N;M]``
     - X-C=N outside an aromatic ring: chloropyridazinone, chloroisoxazoline, imidoyl chloride; the ring N is context, not consumed
   * - ``imidoyl_thioether``
     - ``functional:268``
     - ``[S;D2;z1;x0;!R](-[C;z1])-[C;z2;!a:1]=[N:2]``
     - RS-C=N outside an aromatic ring: an isothiourea or a 2-methylthio-pyrimidinone
   * - ``imine``
     - ``functional:213``
     - ``[N;D2;z2;x0:1]=[C;z2:2]``
     - R-N=C, the reductive-amination intermediate and the aza-Michael acceptor.  ``x0`` excludes the oxime and the hydrazone, which have their own rows
   * - ``imine_nh``
     - ``functional:214``
     - ``[N;D1;z2;x0;h1:1]=[C;z2:2]``
     - HN=C, the free imine
   * - ``isocyanate``
     - ``functional:34``
     - ``[N;z2;x0;D2:1]=[C:2]=[O:3]``
     - R-N=C=O
   * - ``isocyano``
     - ``functional:106``
     - ``[N;+;D2;x0:1]#[C;-;D1:2]``
     - R-NC, charge-separated
   * - ``isopropylamine``
     - ``functional:386``
     - ``[N;D1;z1;x0:1]-[C;D3;x1]([C;D1;x0])[C;D1;x0]``
     - iPrNH2, the transaminase amine donor; acetone leaves
   * - ``isothiocyanate``
     - ``functional:109``
     - ``[N;z2;x0;D2:1]=[C:2]=[S;D1:3]``
     - R-N=C=S
   * - ``isoxazole``
     - ``functional:137``
     - ``[O;a;D2;r5:1]:[N;a;D2;r5:2]``
     - adjacent aromatic O/N pair
   * - ``katritzky_pyridinium``
     - ``functional:282``
     - ``[C;z1:1]-;!@[N;a;+;D3](:[C;a]-;!@[C;a]):[C;a]-;!@[C;a]``
     - an N-alkyl-2,6-diarylpyridinium, the deaminative alkyl radical source
   * - ``ketimine``
     - ``functional:396``
     - ``[N;D2;z2;x0:1]=[C;z2;D3;x1:2]``
     - R2C=NR, a ketone-derived imine, cyclic 1-pyrroline included; ``x1`` on the carbon excludes the amidine
   * - ``ketone``
     - ``functional:21``
     - ``[O;z2;x0:1]=[C;D3;x1;z2:2]``
     - R2C=O
   * - ``ketone_alpha_ch``
     - ``functional:391``
     - ``[O;z2;x0:1]=[C;D3;x1;z2:2]-[C;z1;h1,h2,h3:3]``
     - a ketone with its carbonyl oxygen and an alpha C-H numbered, the enol triflate and silyl enol ether site
   * - ``knorr_dicarbonyl``
     - ``functional:368``
     - ``[O;D1;z2]=[C;z2;x1:1]-[C;z1;D2,D3:2]-[C;z2;x1;D2,D3:3]=[O;D1;z2]``
     - a 1,3-diketone or 2-acyl aldehyde whose central carbon can aromatize; both oxygens leave.  ``D2,D3`` on :3 alone restates what every carbonyl carbon is and breaks the pattern symmetry, so both condensation orientations are enumerated
   * - ``knorr_ketoester``
     - ``functional:370``
     - ``[O;D1;z2]=[C;z2;x1:1]-[C;z1;D2,D3:2]-[C;z2;x2:3](=[O;D1:4])-[O;D2]-[C]``
     - a beta-ketoester or 2-formyl ester: the carbonyl condenses, the ester acylates, giving the pyrazolone
   * - ``lactam_1_halide``
     - ``functional:244``
     - ``[F,Cl,Br,I;D1]-[C;z2;r6:1]1=[C,N;z2;M]-[C;D3;z2;M]-[N;D3;M]-[C,N;z2,z4;M]=,:[C,N;z2,z4;M]1``
     - Activated C(sp2)-X on an N-substituted six-membered lactam, halide at C4 of a 2-pyridinone.  Every ring atom but the reacting carbon is masked, so a reaction row bonds to :1 and restates nothing
   * - ``lactam_2_halide``
     - ``functional:245``
     - ``[F,Cl,Br,I;D1]-[C;z2;r6:1]=1-[N;D3;M]-[C;D3;z2;M]-[C,N;z2,z4;M]=,:[C,N;z2,z4;M]-[C,N;z2;M]=1``
     - Activated C(sp2)-X on an N-substituted six-membered lactam, halide at C6 of a 2-pyridinone, adjacent to the ring nitrogen
   * - ``lactam_3_halide``
     - ``functional:246``
     - ``[F,Cl,Br,I;D1]-[C;z2;r6:1]=1-[N;D3;M]-[C,N;z2,z4;M]=,:[C,N;z2,z4;M]-[C;D3;z2;M]-[C,N;z2;M]=1``
     - Activated C(sp2)-X on an N-substituted six-membered lactam, halide at C2 of a 4-pyridinone, adjacent to the ring nitrogen
   * - ``lactam_4_halide``
     - ``functional:247``
     - ``[F,Cl,Br,I;D1]-[C;z2;r6:1]=1-[C;D3;z2;M]-[N;D3;M]-[C,N;z2,z4;M]=,:[C,N;z2,z4;M]-[C,N;z2;M]=1``
     - Activated C(sp2)-X on an N-substituted six-membered lactam, halide at C3 of a 2-pyridinone, adjacent to the carbonyl
   * - ``lactam_5_halide``
     - ``functional:248``
     - ``[F,Cl,Br,I;D1]-[C;z2;r6:1]1=[C,N;z2;M]-[N;D3;M]-[C,N;z2,z4;M]=,:[C,N;z2,z4;M]-[C;D3;z2;M]1``
     - Activated C(sp2)-X on an N-substituted six-membered lactam, halide at C3 of a 4-pyridinone, adjacent to the carbonyl
   * - ``lactam_6_halide``
     - ``functional:249``
     - ``[F,Cl,Br,I;D1]-[C;z2;r6:1]1=[C,N;z2;M]-[N;D3;M]-[C;D3;z2;M]-[C,N;z2,z4;M]=,:[C,N;z2,z4;M]1``
     - Activated C(sp2)-X on an N-substituted six-membered lactam, halide at C5 of a 2-pyridinone
   * - ``lactim_ether``
     - ``functional:418``
     - ``[N;D2;z2;x0:1]=[C;D3;z2;x2:2]-[O;D2]-[C;z1]``
     - R-C(=NR)-OR', the N-substituted imidate or lactim ether; the alkoxide leaves
   * - ``lactone``
     - ``functional:389``
     - ``[O;z2;x0:1]=[C;D3;x2;z2:2]-[O;D2;x0;R:3]``
     - a ring ester whose ring oxygen stays on the product: opening it by a nucleophile gives the hydroxy acyl compound, which ``ester`` (alkoxy unnumbered) deletes
   * - ``maleimide``
     - ``functional:177``
     - ``[N;r5:1]1-[C;z2;r5:2](=[O:3])-[C;z2;r5:4]=[C;z2;r5:5]-[C;z2;r5:6]1=[O:7]``
     - maleimide, the thiol-conjugation warhead
   * - ``malonic_acid``
     - ``functional:382``
     - ``[O;D1;z1;x0]-[C;z2;x2;D3](=[O;D1])-[C;z1;x0;h1,h2:1]-[C;z2;x2;D3:2](=[O;D1:3])-[O;D1;z1:4]``
     - HOOC-CH2-COOH; one carboxyl is lost as CO2
   * - ``malonic_half_ester``
     - ``functional:324``
     - ``[O;D1;*]-[C;z2;x2](=[O])-[C;z1;h1,h2:1]-[C;z2;x2;M](=[O;M])-[O;D2;M]``
     - HOOC-CH2-COOR or its salt, the decarboxylative Claisen nucleophile
   * - ``metalate_carbanion``
     - ``functional:241``
     - ``[Mg,Zn,Cu,Li;D1;+:1]-[C:2]``
     - The ionic spelling of an organometallic, which the neutral ``alkyl_grignard`` row cannot match; the ``*_zinc`` rows take it too
   * - ``metalloid_methyl_halide``
     - ``functional:271``
     - ``[Cl,Br,I;D1]-[C;z1;x2;D2:1]-[B,Si,Sn;*;M]``
     - X-CH2-M with M boron, silicon or tin: ICH2SnBu3, BrCH2BF3K
   * - ``methanol``
     - ``functional:254``
     - ``[O;D1;z1;x0:1]-[C;D1:2]``
     - CH3-OH, which ``primary_alcohol`` (``D2``) leaves out
   * - ``methyl_aryl_thioether``
     - ``functional:401``
     - ``[C;z1;D1:1]-[S;D2;z1;x0;M]-[C;a;M]``
     - the S-methyl of a thioanisole, chlorinated by NCS through the sulfonium
   * - ``methyl_ester``
     - ``functional:170``
     - ``[O;z2;x0:1]=[C;D3;x2;z2:2]-[O;D2;x0:3]-[C;D1]``
     - R-CO-OMe; distinct saponification and transesterification behaviour
   * - ``n_alkoxy_amide_nh``
     - ``functional:281``
     - ``[N;D2;z1;x1:1](-[C;z2;x2:2]=[O:3])-[O;D2:4]``
     - a hydroxamate N-H, isoxazolidin-3-one included
   * - ``n_halo_imide``
     - ``functional:227``
     - ``[Cl,Br,I;D1:1]-[N;D3:2]``
     - N-X, the NBS/NCS halogen source
   * - ``n_hydroxylamine``
     - ``functional:216``
     - ``[O;D1;z1;x1;h1:1]-[N;z1;x1:2]-[C;z1,a:3]``
     - HO-N with carbon on the nitrogen.  The carbon is stated to exclude the hydroxamic acid, whose nitrogen carries a carbonyl
   * - ``n_methyl_amine``
     - ``functional:399``
     - ``[C;D1;z1;x1]-[N;D3;z1;x0:1](-[C;z1:2])-[C;z1:3]``
     - R2N-CH3 on an sp3 amine; the methyl is the leaving carbon
   * - ``n_substituted_azole``
     - ``functional:205``
     - ``[N;a;D3:1]-[C;z1,a:2]``
     - The substituted pyrrole-type nitrogen -- an N-alkylated or N-arylated azole, as distinct from the ``nh_azole`` that could still be alkylated
   * - ``nh2_hydrazine``
     - ``functional:317``
     - ``[N;D1;z1;x1:1]-[N;D2,D3;z1:2]``
     - a substituted hydrazine, hydrazide, carbazate or sulfonylhydrazide NH2
   * - ``nh2_hydrazone``
     - ``functional:384``
     - ``[N;D1;z1;x1:1]-[N;D2;z2;x1:2]=[C;z2:3]``
     - R2C=N-NH2; ``hydrazone`` (``D2``) is the N-substituted one
   * - ``nh_carbamate``
     - ``functional:261``
     - ``[N;D1,D2;z1;x0:1]-[C;D3;z2;x3:2](=[O;D1:3])-[O;D2:4]``
     - a carbamate N-H, BocNH2 and the oxazolidinone included; ``carbamate`` does not bound the N
   * - ``nh_hydrazine``
     - ``functional:279``
     - ``[N;D2;z1;x1;h1:1]-[N:2]``
     - a substituted hydrazine, hydrazide or hydrazone N-H: the N-N nitrogen ``secondary_amine`` (``x0``) leaves out
   * - ``nh_sulfonamide``
     - ``functional:263``
     - ``[N;D1,D2;z1;x1:1]-[S;D4:2](=[O:3])=[O:4]``
     - a sulfonamide N-H; ``sulfonamide`` matches the tertiary one too
   * - ``nh_sulfoximine``
     - ``functional:387``
     - ``[N;D1;h1:1]=[S;D4:2](=[O;D1:3])(-[C:4])-[C:5]``
     - an N-H sulfoximine
   * - ``nh_thiocarbonyl``
     - ``functional:352``
     - ``[S;D1;z2;x0:1]=[C;D3;z2;x2,x3:2]-[N;D2;h1:3]``
     - a thioamide or thiourea N-H, drawn as the thione: S-alkylation gives the thioimidate or isothiourea
   * - ``nh_thiourea``
     - ``functional:156``
     - ``[N;z1;h1,h2;!R:1]-[C;D3;z2;x3:2](=[S;D1:3])-[N;z1;h1,h2;!R:4]``
     - the thiourea half of the same
   * - ``nh_urea``
     - ``functional:155``
     - ``[N;z1;h1,h2;!R:1]-[C;D3;z2;x3:2](=[O:3])-[N;z1;h1,h2;!R:4]``
     - acyclic urea with an N-H on both nitrogens, the Biginelli subset
   * - ``nh_urea_any``
     - ``functional:262``
     - ``[N;D1,D2;z1;x0:1]-[C;D3;z2;x3:2](=[O;D1:3])-[N;z1:4]``
     - a urea N-H, ring or chain; ``nh_urea`` is the acyclic both-NH subset
   * - ``nitrile``
     - ``functional:32``
     - ``[N;D1;z3;x0:1]#[C;D2;x1:2]``
     - R-CN; ``z3`` here is a genuine sp nitrogen
   * - ``nitro``
     - ``functional:33``
     - ``[N;D3;x2;+:1]([O;-:2])=[O:3]``
     - R-NO2, charge-separated as chython stores it
   * - ``nitroalkane_ch``
     - ``functional:381``
     - ``[C;z1;x1;h2,h3:1]-[N;+;D3;M](=[O;D1;M])-[O;-;D1;M]``
     - MeNO2 or RCH2NO2: the carbon that condenses to the nitroalkene
   * - ``nitroalkene``
     - ``functional:412``
     - ``[C;z2;x0:1]=[C;z2:2]-[N;+;D3:3](=[O;D1:4])-[O;-;D1:5]``
     - C=C-NO2, numbered from the beta carbon
   * - ``nitroalkene_methyl``
     - ``functional:331``
     - ``[C;D1;z1:1]-[C;z2;M]=[C;z2;M]-[N;+;M](=[O;M])-[O;-;M]``
     - a methyl on the far carbon of a nitroalkene, vinylogously acidic
   * - ``nitroarene_ortho_methyl``
     - ``functional:312``
     - ``[C;D1;z1:1]-[C;a;M]:[C;a;M]-[N;+;M](=[O;M])-[O;-;M]``
     - an aryl methyl ortho to NO2, the Leimgruber-Batcho substrate
   * - ``nitroso``
     - ``functional:108``
     - ``[N;D2;z2:1]=[O;D1:2]``
     - R-N=O
   * - ``nitroxide``
     - ``functional:250``
     - ``[O;D1;h0;z1:1]-[N;D3;z1:2](-[C;z1,z4;D3,D4:3])-[C;z1,z4;D3,D4:4] |^1:0|``
     - R2N-O•, the radical on oxygen -- TEMPO and its kin.  ``D3,D4`` on both alpha carbons is the isolability: a nitroxide with an alpha hydrogen disproportionates, so dimethylnitroxide is a decoy and not a match.  The hydroxylamine and the oxoammonium cation are the same skeleton at the other two oxidation states
   * - ``o_amidoaniline``
     - ``functional:325``
     - ``[C;z2;x2;D2,D3:1](=[O])-[N;D2;z1;x0:2]-[C;a:3]:[C;a:4]-[N;D1,D2;z1;x0:5]``
     - an anilide with an ortho NH or NH2, the cyclodehydration substrate
   * - ``o_aminoanilide``
     - ``functional:362``
     - ``[N;D1;z1;x0:1]-[C;a:2]:[C;a:3]-[N;D2,D3;z1:4]-[C;D3;z2;x2:5]=[O;D1]``
     - an anilide with a free ortho amine N-H: cyclodehydration closes the benzimidazole, the amide O leaving as water
   * - ``o_aminobenzaldehyde``
     - ``functional:149``
     - ``[N;D1;z1;x0:1]-[C;a:2]:[C;a:3]-[C;D2;z2;x1:4]=[O:5]``
     - 2-aminobenzaldehyde, the Friedlander partner
   * - ``o_aminophenol``
     - ``functional:147``
     - ``[N;D1;z1;x0:1]-[C;a:2]:[C;a:3]-[O;D1:4]``
     - 2-aminophenol
   * - ``o_aminothiophenol``
     - ``functional:148``
     - ``[N;D1;z1;x0:1]-[C;a:2]:[C;a:3]-[S;D1:4]``
     - 2-aminothiophenol
   * - ``o_diaminoarene``
     - ``functional:146``
     - ``[N;D1;z1;x0:1]-[C;a:2]:[C;a:3]-[N;D1,D2;z1;x0:4]``
     - ortho-phenylenediamine
   * - ``o_haloaniline``
     - ``functional:169``
     - ``[N;D1;z1;x0:1]-[C;a:2]:[C;a:3]-[Cl,Br,I;D1]``
     - 2-haloaniline, the Larock indole substrate
   * - ``o_nitroaniline``
     - ``functional:366``
     - ``[N;D1,D2;z1;x0:1]-[C;a:2]:[C;a:3]-[N;+;D3:4](=[O;D1])-[O;D1;-]``
     - an o-nitroaniline N-H: reduced in situ, its nitro N becomes a benzimidazole ring N
   * - ``oriented_epoxide``
     - ``functional:355``
     - ``[O;D2;r3:1]1-[C;r3:2]-[C;r3;z1:3]-1``
     - an oxirane whose :3 is the carbon a nucleophile takes.  ``z1`` on :3 alone holds for every oxirane carbon and breaks the pattern symmetry, so a template yields both regioisomers rather than one orientation
   * - ``orthoester``
     - ``functional:363``
     - ``[C;D3,D4;x3;z1:1](-[O;D2]-[C;z1])(-[O;D2]-[C;z1])-[O;D2]-[C;z1]``
     - RC(OR')3, trimethyl orthoformate and kin; all three alkoxides leave
   * - ``oxetane``
     - ``functional:237``
     - ``[O;D2;z1;r4:1]([C;z1;r4:2])[C;z1;r4:3]``
     - The four-membered ether, a carbonyl bioisostere and a slow electrophile
   * - ``oxime``
     - ``functional:162``
     - ``[O;D1;z1;x1:1]-[N;D2;z2;x1:2]=[C:3]``
     - C=N-OH
   * - ``oxindole_ch2``
     - ``functional:314``
     - ``[C;z1;h2;r5:1](-[C;a;M])-[C;r5;M](=[O;M])-[N;M]``
     - the C3 methylene of an oxindole
   * - ``pentafluorosulfanyl``
     - ``functional:209``
     - ``[S;D6:1]([F:2])([F:3])([F:4])([F:5])-[C:6]``
     - SF5, the CF3 replacement
   * - ``peroxide``
     - ``functional:176``
     - ``[O;D2;z1:1]-[O;D2;z1:2]``
     - R-O-O-R
   * - ``phenol``
     - ``functional:18``
     - ``[O;D1;z1;x0:1]-[C;a:2]``
     - Ar-OH
   * - ``phenoxide``
     - ``functional:233``
     - ``[O;D1;z1;x0;-:1]-[C;a:2]``
     - ArO-
   * - ``phosphate_ester``
     - ``functional:210``
     - ``[C:1]-[O;D2:2]-[P;D4:3](=[O;D1:4])(-[O;D2:5])-[O;D2:6]``
     - R-O-PO(OR)(OR), three oxygens on phosphorus
   * - ``phosphine_oxide``
     - ``functional:211``
     - ``[P;D4:1](=[O;D1;x1:2])([C:3])([C:4])[C:5]``
     - R3P=O, the Wittig and Mitsunobu byproduct.  ``x1`` on the oxygen: its only neighbour is phosphorus, which is a heteroatom
   * - ``phosphonate``
     - ``functional:126``
     - ``[P;D4;x3](=[O])(-[O;D2;x1])(-[O;D2;x1])-[C:1]``
     - (RO)2P(=O)R, the HWE reagent
   * - ``phosphonic_acid``
     - ``functional:212``
     - ``[O;D1;z1;h1:1]-[P;D4:2](=[O;D1:3])-[O;D1;z1;h1:4]``
     - R-PO(OH)2
   * - ``phosphonium_salt``
     - ``functional:290``
     - ``[P;+;D4;x0](-[C;a])(-[C;a])(-[C;a])-[C;z1:1]``
     - R-CH2-PPh3+, the Wittig reagent recorded as the salt the ylide is made from
   * - ``phosphonium_ylide``
     - ``functional:125``
     - ``[P;D4;z2;x0]=[C:1]``
     - R3P=CR2, the Wittig reagent
   * - ``phosphorus_halide``
     - ``functional:228``
     - ``[F,Cl,Br,I;D1:1]-[P:2]``
     - P-X
   * - ``polyarene``
     - ``functional:179``
     - ``[A;a:1]:[A;a;D3:2](:[A;a:3]):[A;a:4]:[A;a;D3:5](:[A;a:6]):[A;a:7]``
     - two fusion atoms one atom apart, as in a linear triarene
   * - ``primary_alcohol``
     - ``functional:15``
     - ``[O;D1;z1;x0:1][C;D2;x1;z1:2]``
     - R-CH2-OH
   * - ``primary_amide``
     - ``functional:30``
     - ``[N;D1;z1;x0:1][C;z2;x2;D2,D3:2]=[O:3]``
     - R-CO-NH2.  ``D2,D3`` because a FORMAMIDE carbonyl is D2, hydrogens not counting toward ``D``, and ``D3`` alone reported formamide as carrying no group at all
   * - ``primary_amidine_amine``
     - ``functional:93``
     - ``[N;D1;z1;x0:1]-[C;z2:2]=[N:3]``
     - NH2 on an sp2 carbon that also carries C=N
   * - ``primary_amine``
     - ``functional:25``
     - ``[N;D1;z1;x0:1][C;z1:2]``
     - R-NH2 on sp3 carbon
   * - ``primary_aniline``
     - ``functional:26``
     - ``[N;D1;z1;x0:1][C;a:2]``
     - Ar-NH2
   * - ``primary_carbinol_alpha_ch``
     - ``functional:328``
     - ``[O;D1;z1;x0:1]-[C;D2;x1;z1:2]-[C;z1;h1,h2:3]``
     - RCH2-OH with an sp3 CH or CH2 beside the carbinol
   * - ``primary_sulfinamide``
     - ``functional:414``
     - ``[N;D1;z1:1]-[S;D3;z2:2](=[O;D1:3])-[C:4]``
     - R-S(=O)-NH2, the Ellman or Davis auxiliary
   * - ``propenyl``
     - ``functional:329``
     - ``[C;D1;z1]-[C;z2;x0;D2:1]=[C;z2;x0;D2:2]``
     - CH3-CH=CH-R, the internal alkene that metathesis cleaves to propene
   * - ``pyrazol_3_one``
     - ``functional:394``
     - ``[O;z2;x0:1]=[C;z2;x2;D3;r5:2]1-[N;z1;D2;h1;r5:3]-[N;z1;D3;r5;M]-[C;z2;r5;M]=[C;z2;r5;M]-1``
     - an N1-substituted 1,2-dihydropyrazol-3-one drawn as its NH tautomer, the 3-hydroxypyrazole in disguise
   * - ``pyrazol_5_one``
     - ``functional:393``
     - ``[O;z2;x0:1]=[C;z2;x2;D3;r5:2]1-[N;z1;D3;r5;M]-[N;z1;D2;h1;r5:3]-[C;z2;r5:4]=[C;z2;r5:5]-1``
     - an N1-substituted 1,2-dihydropyrazol-3-one drawn as its NH tautomer, the 5-hydroxypyrazole in disguise
   * - ``pyrazole``
     - ``functional:135``
     - ``[N;h1;D2;a;r5:1]:[N;h0;D2;r5:2]``
     - adjacent N-H/N pair
   * - ``pyridazine``
     - ``functional:138``
     - ``[N;a;D2;r6:1]:[N;a;D2;r6:2]``
     - adjacent aromatic N pair in a 6-ring
   * - ``pyridazine_dione``
     - ``functional:335``
     - ``[O;z2;x0]=[C;z2;x2;D3;r6:1]1-[N;z1;D2;r6:2]-[N;z1;D2;r6:3]-[C;z2;x2;D3;r6:4](=[O;z2;x0])-[C;z2,z4;M]=,:[C;z2,z4;M]1``
     - maleic hydrazide class, 3,6-dione with both N-H: the double deoxychlorination substrate
   * - ``pyridine_n``
     - ``functional:40``
     - ``[N;a;D2;h0:1]``
     - an aromatic nitrogen with no hydrogen, the N-oxidation substrate
   * - ``pyridine_ring``
     - ``functional:405``
     - ``[N;a;D2;h0;R1:1]1:[C;a;R1:2]:[C;a;R1:3]:[C;a;R1:4]:[C;a;R1:5]:[C;a;R1:6]:1``
     - an unfused pyridine ring, all six atoms numbered for the full hydrogenation
   * - ``pyrrole``
     - ``functional:134``
     - ``[N;h1;D2;a;r5:1]``
     - azole N-H; the H stays out of a template for tautomerism
   * - ``quaternary_ammonium``
     - ``functional:236``
     - ``[N;D4;z1;x0;+;h0:1]``
     - R4N+, a phase-transfer catalyst rather than a protonation state
   * - ``redox_active_ester``
     - ``functional:128``
     - ``[C;z1:1]-[C;z2;x2;D3:2](=[O:3])-[O;D2;x1]-[N;D3;x1;r5](-[C;z2;r5]=[O])-[C;z2;r5]=[O]``
     - NHPI/NHS ester; decarboxylative coupling transfers the alkyl
   * - ``ring_ch_ch``
     - ``functional:404``
     - ``[C;z1;h1,h2;R:1](-;@[C;z1;h1,h2;R:2])-[A;z2,z4;M]``
     - two adjacent ring sp3 C-H carbons, one of them beside an sp2 or aromatic atom: the desaturation site that conjugates or aromatizes
   * - ``ring_enamine_nh2``
     - ``functional:265``
     - ``[N;D1;z1;x0:1]-[C;z2;x1,x2;r5,r6:2]=[C;z2:3]``
     - an NH2 on a ring C=C that ``primary_aniline`` (``[C;a]``) does not see: aminopyridone, 3-aminocoumarin, 6-aminouracil
   * - ``ring_vinyl_ch``
     - ``functional:340``
     - ``[C;z2;D2;h1:1]=;@[C,N;z2;M]``
     - a C-H on a non-aromatic ring double bond: the azinone, enaminone and uracil C-H ``arene_ch`` (aromatic only) leaves out
   * - ``secondary_alcohol``
     - ``functional:16``
     - ``[O;D1;z1;x0:1][C;D3;x1;z1:2]``
     - R2CH-OH
   * - ``secondary_amide``
     - ``functional:31``
     - ``[N;D2;z1;x0:1][C;z2;x2;D2,D3:2]=[O:3]``
     - R-CO-NH-R', N-methylformamide included -- see ``primary_amide`` on ``D2,D3``
   * - ``secondary_amine``
     - ``functional:27``
     - ``[N;D2;z1;x0:1]([C;z1:2])[C;z1:3]``
     - R2NH, both substituents sp3 carbon
   * - ``secondary_aniline``
     - ``functional:28``
     - ``[N;D2;z1;x0:1]([C;a:2])[C;z1:3]``
     - Ar-NH-R
   * - ``secondary_carbinol_alpha_ch``
     - ``functional:327``
     - ``[O;D1;z1;x0:1]-[C;D3;x1;z1:2]-[C;z1;h1,h2,h3:3]``
     - R2CH-OH with an sp3 C-H beside the carbinol
   * - ``secondary_enamine``
     - ``functional:280``
     - ``[N;D2;z1;x0:1](-[C;z1:2])-[C;z2:3]=[C:4]``
     - R-NH-C=C, the vinylogous-amide N-H ``secondary_amine`` (both carbons ``z1``) leaves out
   * - ``silyl_enol_ether``
     - ``functional:343``
     - ``[C;z2;x0;D1,D2:1]=[C;z2;x1;D2,D3:2]-[O;D2:3]-[Si;D4]``
     - C=C-O-SiR3, the enol equivalent an electrophilic halogen source takes to the alpha-halo ketone
   * - ``silyl_ether``
     - ``functional:189``
     - ``[Si;D4:1]-[O;D2;z1:2]-[C:3]``
     - R3Si-O-R, a protected alcohol.  ``protective.tsv`` removes the named ones; this reports any
   * - ``succinimidyl_carbamate``
     - ``functional:357``
     - ``[N;z1:1]-[C;D3;z2;x3:2](=[O;D1:3])-[O;D2;x1]-[N;D3;x1;r5](-[C;z2;r5]=[O])-[C;z2;r5]=[O]``
     - N-CO-O-NSu, the bench-stable isocyanate equivalent; the N-hydroxysuccinimide leaves
   * - ``succinimidyl_carbonate``
     - ``functional:88``
     - ``[O;D2:1]-[C;z2;x3;D3:2](=[O:3])-[O;D2;x1]-[N;D3;x1;r5](-[C;z2;r5]=[O])-[C;z2;r5]=[O]``
     - NHS/DSC carbonate; numbered like chloroformate so it drops into the same rows
   * - ``sulfamate``
     - ``functional:199``
     - ``[N:1]-[S;D4;x4:2](=[O;D1:3])(=[O;D1:4])-[O;D2:5]-[C:6]``
     - N-SO2-O-R
   * - ``sulfamide``
     - ``functional:198``
     - ``[N:1]-[S;D4;x4:2](=[O;D1:3])(=[O;D1:4])-[N:5]``
     - N-SO2-N, nitrogen on both sides
   * - ``sulfinamide``
     - ``functional:201``
     - ``[S;D3;z2:1](=[O;D1:2])(-[N:3])[C:4]``
     - R-S(=O)-N, the Ellman auxiliary and the sulfoximine precursor
   * - ``sulfinate_ester``
     - ``functional:202``
     - ``[S;D3;z2:1](=[O;D1:2])(-[O;D2:3]-[C:4])[C:5]``
     - R-S(=O)-O-R
   * - ``sulfinimine``
     - ``functional:415``
     - ``[C;z2:1]=[N;D2;z2:2]-[S;D3:3](=[O;D1:4])-[C:5]``
     - R-S(=O)-N=CR2, the N-sulfinyl imine
   * - ``sulfinyl_ketimine``
     - ``functional:395``
     - ``[N;D2;z2;x1:1](=[C;z2;D3;x1:2])-[S;D3:3]=[O;D1:4]``
     - a ketone-derived N-sulfinyl imine (Ellman), the auxiliary S kept and numbered
   * - ``sulfo``
     - ``functional:97``
     - ``[S;D4;x3:1](=[O:2])(=[O:3])-[O;D1:4]``
     - R-SO3H, distinct from a sulfonyl halide or a sulfonamide
   * - ``sulfonamide``
     - ``functional:100``
     - ``[S;D4;x3:1](=[O:2])(=[O:3])-[N;z1:4]``
     - R-SO2-NR2
   * - ``sulfonate``
     - ``functional:231``
     - ``[O;D1;z1;-:1]-[S;D4:2](=[O;D1:3])=[O;D1:4]``
     - The deprotonated sulfonic acid
   * - ``sulfonate_ester``
     - ``functional:197``
     - ``[C:1]-[O;D2;z1:2]-[S;D4:3](=[O;D1:4])=[O;D1:5]``
     - R-O-SO2-R, any sulfonate ester.  The named ones each have a row; this catches the rest
   * - ``sulfone``
     - ``functional:38``
     - ``[S;D4:1](=[O:2])(=[O:3])([C:4])[C:5]``
     - R-SO2-R'
   * - ``sulfonyl_anhydride``
     - ``functional:101``
     - ``[S;D4;x3:1](=[O:2])(=[O:3])-[O;D2]-[S;D4;x3](=[O])=[O]``
     - R-SO2-O-SO2-R
   * - ``sulfonyl_chloride``
     - ``functional:98``
     - ``[Cl;D1]-[S;D4:1](=[O:2])=[O:3]``
     - R-SO2Cl
   * - ``sulfonyl_fluoride``
     - ``functional:99``
     - ``[F;D1]-[S;D4:1](=[O:2])=[O:3]``
     - R-SO2F
   * - ``sulfonylhydrazide``
     - ``functional:144``
     - ``[N;D1;z1;x1:1]-[N;D2;z1;x2:2]-[S;D4;x3:3](=[O:4])=[O:5]``
     - R-SO2-NH-NH2, the reagent that makes the above
   * - ``sulfonylhydrazone``
     - ``functional:143``
     - ``[C;z2:1]=[N;D2;z2;x1:2]-[N;D2;z1;x2:3]-[S;D4;x3](=[O])=[O]``
     - tosylhydrazone and kin, the Bamford-Stevens diazo precursor
   * - ``sulfoxide``
     - ``functional:37``
     - ``[S;D3;z2:1](=[O:2])([C:3])[C:4]``
     - R-S(=O)-R'
   * - ``sulfoximine``
     - ``functional:200``
     - ``[S;D4;z5:1](=[O;D1:2])(=[N:3])([C:4])[C:5]``
     - R2S(=O)=N, the sulfone bioisostere.  A stereocentre at sulfur, which is why the substituents are numbered rather than left to a template
   * - ``terminal_alkene``
     - ``functional:1``
     - ``[C;z2;x0;D1:1]=[C;z2;x0;D2,D3:2]``
     - monosubstituted C=C
   * - ``terminal_alkyne``
     - ``functional:3``
     - ``[C;z3;x0;D1:1]#[C;D2:2]``
     - monosubstituted C#C, the Sonogashira nucleophile
   * - ``terminal_epoxide``
     - ``functional:171``
     - ``[O;D2;r3:1]1-[C;r3:2]-[C;r3;D2:3]-1``
     - epoxide oriented for opening: ``:3`` is the CH2 a nucleophile attacks
   * - ``tertiary_alcohol``
     - ``functional:17``
     - ``[O;D1;z1;x0:1][C;D4;x1;z1:2]``
     - R3C-OH
   * - ``tertiary_alcohol_with_alpha_h``
     - ``functional:70``
     - ``[O;D1;z1;x0:1]-[C;D4;x1;z1:2]-[C;z1;h1,h2:3]``
     - R3C-OH with an eliminable alpha C-H, the dehydration substrate
   * - ``tertiary_amide``
     - ``functional:190``
     - ``[N;D3;z1;x0:1]-[C;D2,D3;z2;x2:2]=[O;D1;z2;x0:3]``
     - R-CO-NR2, the level ``primary_amide`` and ``secondary_amide`` left out.  ``x0`` on the nitrogen because an amide N's three neighbours are all CARBON; ``x2`` on the carbonyl is what excludes a carbamate and a urea, whose carbonyl is ``x3``
   * - ``tertiary_amine``
     - ``functional:29``
     - ``[N;D3;z1;x0:1]([C;z1,a:2])([C;z1,a:3])[C;z1,a:4]``
     - R3N, the N-oxidation substrate.  The three substituents are stated because the bare ``[N;D3;z1;x0]`` also matched every TERTIARY AMIDE and every carbamate -- an N-Boc amine reported as an amine, and ``nitrogen_oxidation`` offered to oxidize DMF.  ``primary_amine`` and ``secondary_amine`` always stated theirs
   * - ``tertiary_carbinol_alpha_ch``
     - ``functional:326``
     - ``[O;D1;z1;x0:1]-[C;D4;x1;z1:2]-[C;z1;h1,h2,h3:3]``
     - R3C-OH with an sp3 C-H beside the carbinol
   * - ``tetramethylstannane``
     - ``functional:257``
     - ``[Sn;D4](-[C;D1])(-[C;D1])(-[C;D1])-[C;D1:1]``
     - Me4Sn, the Stille methylating agent; ``alkyl_stannane`` excludes the methyls of a Me3Sn dummy
   * - ``thiazolidinone_ch2``
     - ``functional:313``
     - ``[C;z1;h2;r5:1](-[S;r5;M])-[C;r5;M]=[O;M]``
     - the ring CH2 between S and C=O of a thiazolidinone, rhodanine or thiazolone
   * - ``thioamide``
     - ``functional:145``
     - ``[S;z2;x0;D1:1]=[C;D2,D3;x2:2]-[N;D1:3]``
     - R-CS-NH2, the Hantzsch thiazole partner.  ``D2,D3`` as in ``primary_amide``
   * - ``thiocarbamate``
     - ``functional:219``
     - ``[N:1]-[C;D3;z2:2](=[S;D1;x0:3])-[O;D2:4]``
     - N-C(=S)-O-R, the Newman-Kwart substrate
   * - ``thiocyanate``
     - ``functional:221``
     - ``[N;D1;z3:1]#[C;D2;z3:2]-[S;D2:3]``
     - R-S-C#N, which is not the isothiocyanate ``S=C=N-R``
   * - ``thioester``
     - ``functional:131``
     - ``[O;z2;x0:1]=[C;D3;x2;z2:2]-[S;D2;z1;x0]``
     - R-CO-S-R', the Liebeskind-Srogl electrophile
   * - ``thioether``
     - ``functional:36``
     - ``[S;D2;z1;x0:1]([C:2])[C:3]``
     - R-S-R'
   * - ``thiol``
     - ``functional:35``
     - ``[S;x0;D1;z1:1][C;z1:2]``
     - R-SH
   * - ``thiolate``
     - ``functional:234``
     - ``[S;D1;z1;x0;-:1]-[C:2]``
     - RS-
   * - ``thione``
     - ``functional:218``
     - ``[S;D1;z2;x0:1]=[C;D3;z2;x1:2]``
     - C=S with one heteroatom on the carbon, the thioketone
   * - ``thiophene_s``
     - ``functional:203``
     - ``[S;a;D2:1]``
     - The aromatic divalent sulfur of a thiophene, thiazole or isothiazole
   * - ``thiourea``
     - ``functional:154``
     - ``[N;z1:1]-[C;D3;z2;x3:2](=[S;D1:3])-[N;z1:4]``
     - R2N-CS-NR2
   * - ``togni_cf3``
     - ``functional:275``
     - ``[I;D2,D3]-[C;z1;x4:1](-[F;M])(-[F;M])-[F;M]``
     - a hypervalent-iodine CF3 donor, the Togni reagents
   * - ``tosyl_isocyanide``
     - ``functional:165``
     - ``[C;-;D1:1]#[N;+;D2:2]-[C;D2,D3;z1;x2:3]-[S;D4;x2](=[O])=[O]``
     - TosMIC, the Van Leusen oxazole reagent
   * - ``trialkylborane``
     - ``functional:385``
     - ``[C;z1:1]-;!@[B;D3;x0](-[C;z1])-[C;z1]``
     - R3B: Et3B, or a B-alkyl-9-BBN whose exocyclic alkyl transfers
   * - ``triarylmethyl_radical``
     - ``functional:253``
     - ``[C;D3;h0;z1:1](-[C;a:2])(-[C;a:3])-[C;a:4] |^1:0|``
     - Ar3C•, Gomberg's radical and the trityl spin labels.  ``_radical`` is in the name because ``protective.tsv``'s trityl rows are this skeleton without one
   * - ``trifluoromethoxy``
     - ``functional:207``
     - ``[F:1][C;D4;z1;x4:2]([F:3])([F:4])-[O;D2:5]``
     - OCF3.  ``x4`` on the carbon: three fluorines and the oxygen all count
   * - ``trifluoromethyl``
     - ``functional:172``
     - ``[C;D4;z1;x3:1](-[F:2])(-[F:3])-[F:4]``
     - -CF3
   * - ``trihalomethyl``
     - ``functional:226``
     - ``[Cl,Br,I;D1:1]-[C;D4;z1:2](-[Cl,Br,I;D1:3])-[Cl,Br,I;D1:4]``
     - CX3 for X other than fluorine
   * - ``urea``
     - ``functional:153``
     - ``[N;z1:1]-[C;D3;z2;x3:2](=[O:3])-[N;z1:4]``
     - R2N-CO-NR2
   * - ``vicinal_diol``
     - ``functional:19``
     - ``[O;D1;z1;x0:1]-[C;z1;x1:2]-[C;z1;x1:3]-[O;D1;z1;x0:4]``
     - 1,2-diol, the dihydroxylation product
   * - ``vinyl_azine_ortho``
     - ``functional:410``
     - ``[C;z2;x0;D1:1]=[C;z2;D2:2]-[C;a;r6:3]:[N;a;D2;h0;r6;M]``
     - CH2=CH- beside a six-ring aza N, the 2-vinylazine Michael acceptor
   * - ``vinyl_azine_para``
     - ``functional:411``
     - ``[C;z2;x0;D1:1]=[C;z2;D2:2]-[C;a:3]1:[C,N;a;M]:[C,N;a;M]:[N;a;D2;h0;M]:[C,N;a;M]:[C,N;a;M]:1``
     - CH2=CH- para to a six-ring aza N, the 4-vinylazine Michael acceptor
   * - ``vinyl_sulfone``
     - ``functional:175``
     - ``[S;D4:1](=[O:2])(=[O:3])(-[O,N,C:4])-[C;z2:5]=[C:6]``
     - vinyl sulfone, sulfonamide or sulfonate -- a covalent warhead.  chython 2 wrote the third substituent ``[O,N]`` and so missed every true sulfone; widened here, which is why the ``x3`` it carried had to go
   * - ``vinyl_sulfonyl``
     - ``functional:409``
     - ``[C;z2;x0:1]=[C;z2:2]-[S;D4:3](=[O;D1:4])=[O;D1:5]``
     - C=C-SO2X for any X -- sulfone, sulfonamide, sulfonate, sulfonyl fluoride -- numbered from the beta carbon
   * - ``vinylogous_azinone``
     - ``functional:333``
     - ``[O;z2;x0:1]=[C;z2;x1,x2;D3;r6:2]1-[C;z2;r6:3]=[C;z2;r6:4]-[N;z1;D2;r6:5]-[A;M]~[A;M]1``
     - 4-pyridone class: C=O and N-H across a ring C=C, the six-ring lactam ``azinone`` does not name
   * - ``vinylogous_dielectrophile``
     - ``functional:369``
     - ``[O;D1;z2]=[C;z2;x1:1]-[C;z2:2]=[C;z2;x1:3]-[ND3,OD1,OD2;z1]``
     - an enaminone, alkoxymethylene ketone or acyl enol: a 1,3-dielectrophile whose C=C end carries the leaving NR2, OR or OH
   * - ``weinreb_amide``
     - ``functional:127``
     - ``[O:1]=[C;D3;x2:2]-[N;D3;x1](-[C;z1])-[O;D2;x1]``
     - N-methoxy-N-methyl amide; the extra C on N excludes an N-acyloxy imide
   * - ``xanthate``
     - ``functional:220``
     - ``[S;D1,D2:1]-[C;D3;z2:2](=[S;D1;x0:3])-[O;D2:4]``
     - S-C(=S)-O-R, the RAFT agent.  ``D1,D2`` because the second sulfur is a free thiol as often as it is substituted
   * - ``xec_aryl_halide``
     - ``functional:420``
     - ``[Cl,Br,I;D1]-[C;a:1]``
     - ``aryl_halide`` under a second name, for the second electrophile of a biaryl XEC: slots naming two groups never read one aryl halide as both

Protecting groups
-----------------

115 protecting groups, alphabetically.  *Protects* names the functional
group the row reveals, which is a row of the table above.
:meth:`chython.MoleculeContainer.protective_groups` reports a match and
:meth:`chython.MoleculeContainer.deprotect` applies the patch.

Two rows can match one substructure -- a Boc is also a tert-butyl -- and the larger
pattern is served first, so a name here is the most specific group that fits, not
every group that could.

.. list-table::
   :header-rows: 1
   :widths: 22 11 15 26 26

   * - Name
     - Id
     - Protects
     - SMARTS
     - What it removes
   * - ``amine_acyl``
     - ``protective:84``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x2](=O)-[C;D1]``
     - acetyl, revealing an amine
   * - ``amine_alloc``
     - ``protective:56``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;z1;x1][C;D2;x0;z2]=[C;D1]``
     - allyloxycarbonyl, revealing an amine
   * - ``amine_alpha_methyl_methoxybenzyl``
     - ``protective:109``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D3;z1;x1]([C;D1])-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1``
     - 1-(4-methoxyphenyl)ethyl, revealing an amine
   * - ``amine_alpha_methylbenzyl``
     - ``protective:108``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D3;z1;x1]([C;D1])-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - alpha-methylbenzyl, revealing an amine
   * - ``amine_benzhydryl``
     - ``protective:107``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D3;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2``
     - benzhydryl, revealing an amine
   * - ``amine_benzhydrylidene``
     - ``protective:85``
     - ``amine``
     - ``[N;D2:1]=;!@[C;D3;z2;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2``
     - benzhydrylidene (benzophenone imine), revealing an amine
   * - ``amine_benzoate``
     - ``protective:87``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x2](=O)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - benzoyl, revealing an amine
   * - ``amine_benzyl``
     - ``protective:75``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - benzyl, revealing an amine
   * - ``amine_bhoc``
     - ``protective:79``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D3;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - benzhydryloxycarbonyl, revealing an amine
   * - ``amine_boc``
     - ``protective:65``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x3](=O)-[O;D2;x0]-[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]``
     - tert-butoxycarbonyl, revealing an amine
   * - ``amine_bom``
     - ``protective:113``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - benzyloxymethyl, revealing an amine, amide or azole N-H
   * - ``amine_cbz``
     - ``protective:60``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x3](=O)-[O;D2;x0][C;D2;x1;z1][C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - benzyloxycarbonyl (Cbz/Z), revealing an amine
   * - ``amine_chloro_cbz``
     - ``protective:61``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x3](=O)-[O;D2;x0][C;D2;x1;z1][C;a;r6]:1:[C;D3;x1]([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - 4-chlorobenzyloxycarbonyl, revealing an amine
   * - ``amine_chloro_tritil``
     - ``protective:81``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D3;x1]([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - 2-chlorotrityl, revealing an amine
   * - ``amine_dde_enamine``
     - ``protective:70``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x1;D3]([C;D1])=[C;D3;r6;x0;z2]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O``
     - Dde, as the enamine tautomer, revealing an amine
   * - ``amine_dde_imine``
     - ``protective:71``
     - ``amine``
     - ``[N;D2:1]=;!@[C;x1;D3]([C;D1])-[C;D3;r6;x0;z1]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O``
     - Dde, as the imine tautomer, revealing an amine
   * - ``amine_dimethoxybenzyl``
     - ``protective:77``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1``
     - 2,4-dimethoxybenzyl (DMB), revealing an amine
   * - ``amine_ethylcarbamate``
     - ``protective:55``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0][C;D2;x1;z1][C;D1]``
     - ethyl carbamate, revealing an amine
   * - ``amine_fmoc``
     - ``protective:67``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;x1;z1][C;D3;z1;x0;r5]1[C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D3]:2-[C;a;r6]:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C1:3``
     - 9-fluorenylmethyloxycarbonyl, revealing an amine
   * - ``amine_ivdde_enamine``
     - ``protective:72``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x1;D3]([C;D2;x0;z1][C;D3;x0;z1]([C;D1])[C;D1])=[C;D3;r6;x0;z2]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O``
     - ivDde, as the enamine tautomer, revealing an amine
   * - ``amine_ivdde_imine``
     - ``protective:73``
     - ``amine``
     - ``[N;D2:1]=;!@[C;x1;D3]([C;D2;x0;z1][C;D3;x0;z1]([C;D1])[C;D1])-[C;D3;r6;x0;z1]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O``
     - ivDde, as the imine tautomer, revealing an amine
   * - ``amine_mem``
     - ``protective:115``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D2;z1;x1][C;D2;z1;x1][O;D2;x0][C;D1]``
     - 2-methoxyethoxymethyl (MEM), revealing an amine or azole N-H
   * - ``amine_methoxy_benzyl``
     - ``protective:76``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1``
     - 4-methoxybenzyl (PMB), revealing an amine
   * - ``amine_methylcarbamate``
     - ``protective:54``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0][C;D1]``
     - methyl carbamate, revealing an amine
   * - ``amine_mtr``
     - ``protective:69``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[S;D4;x3](=O)(=O)-[C;a;r6]:1:[C;D3;x0]([C;D1]):[C;D3;x0]([C;D1]):[C;D3;x1](-;!@[O;D2;x0][C;D1]):[C;D2]:[C;D3;x0]([C;D1]):1``
     - 4-methoxy-2,3,6-trimethylbenzenesulfonyl, revealing an amine
   * - ``amine_mtt``
     - ``protective:78``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x0]([C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - 4-methyltrityl (Mtt), revealing an amine
   * - ``amine_nosyl``
     - ``protective:64``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[S;D4;x3](=O)(=O)-[C;a;r6]:1:[C;D3;x1]([N+](=O)[O-]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - nitrobenzenesulfonyl (Ns), revealing an amine
   * - ``amine_o_nitrobenzyl``
     - ``protective:111``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D3;x1]([N+](=O)[O-]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - ortho-nitrobenzyl, revealing an amine
   * - ``amine_pbf``
     - ``protective:68``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[S;D4;x3](=O)(=O)-[C;a;r6]:1:[C;D3;x0]([C;D1]):[C;D3;x0]([C;D1]):[C;D3;x1]:2-[O;D2;x0;r5][C;D4;x1]([C;D1])([C;D1])[C;D2;x0;z1][C;D3]:2:[C;D3;x0]([C;D1]):1``
     - 2,2,4,6,7-pentamethyldihydrobenzofuran-5-sulfonyl, revealing an amine
   * - ``amine_phenylsulfonyl``
     - ``protective:62``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[S;D4;x3](=O)(=O)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - benzenesulfonyl, revealing an amine
   * - ``amine_phth``
     - ``protective:74``
     - ``amine``
     - ``[N;D3:1]1[C;z2;x2;D3](=O)[C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D3]:2[C;z2;x2;D3]1=O``
     - phthalimide, revealing an amine
   * - ``amine_pom``
     - ``protective:112``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D2;x2;z1][O;D2;x0]-[C;z2;x2](=[O;D1])-[C;D4]([C;D1])([C;D1])[C;D1]``
     - pivaloyloxymethyl (POM), revealing an amine or azole N-H
   * - ``amine_sem``
     - ``protective:58``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D2;x2;z1][O;D2;x0]-[C;D2;z1;x1][C;D2;x1;z1][Si;D4;z1;x0]([C;D1])([C;D1])[C;D1]``
     - 2-(trimethylsilyl)ethoxymethyl, revealing an amine
   * - ``amine_sulfinyl``
     - ``protective:83``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[S;D3;x2;z2](=O)[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]``
     - tert-butanesulfinyl, revealing an amine
   * - ``amine_tbu``
     - ``protective:86``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]``
     - tert-butyl, revealing an amine
   * - ``amine_teoc``
     - ``protective:57``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;z1;x1][C;D2;x1;z1][Si;D4;z1;x0]([C;D1])([C;D1])[C;D1]``
     - 2-(trimethylsilyl)ethoxycarbonyl, revealing an amine
   * - ``amine_tfa``
     - ``protective:66``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x2](=O)-[C;D4;z1;x3](F)(F)F``
     - trifluoroacetyl, revealing an amine
   * - ``amine_thp``
     - ``protective:82``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D3;x2;z1;r6]1[O;D2][C;D2][C;D2][C;D2][C;D2]1``
     - tetrahydropyranyl, revealing an amine
   * - ``amine_tosyl``
     - ``protective:63``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[S;D4;x3](=O)(=O)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x0]([C;D1]):[C;D2]:[C;D2]:1``
     - 4-toluenesulfonyl (Ts), revealing an amine
   * - ``amine_trimethoxybenzyl``
     - ``protective:110``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]:1[O;D2;x0][C;D1]``
     - 2,4,6-trimethoxybenzyl, revealing an amine
   * - ``amine_tritil``
     - ``protective:80``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - triphenylmethyl (trityl), revealing an amine
   * - ``amine_troc``
     - ``protective:59``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2][C;D4;x3]([Cl;D1])([Cl;D1])[Cl;D1]``
     - 2,2,2-trichloroethoxycarbonyl, revealing an amine
   * - ``amine_veratryl``
     - ``protective:114``
     - ``amine``
     - ``[N;D2,D3:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:1``
     - 3,4-dimethoxybenzyl, revealing an amine
   * - ``carbonyl_dimethoxy``
     - ``protective:52``
     - ``carbonyl``
     - ``[C;D3,D4;z1;x2:1](-;!@[O;D2;x0][C;D1])-;!@[O;D2;x0][C;D1]``
     - dimethyl acetal, revealing a ketone or aldehyde
   * - ``carbonyl_dimethylsulfide``
     - ``protective:49``
     - ``carbonyl``
     - ``[C;D3,D4;z1;x2:1](-;!@[S;D2][C;D1])-;!@[S;D2][C;D1]``
     - bis(methylthio) acetal, revealing a ketone or aldehyde
   * - ``carbonyl_dioxane``
     - ``protective:51``
     - ``carbonyl``
     - ``[C;D3,D4;z1;x2;r6:1]1[O;D2;x0][C;D2;x1;z1][C;D2;x0;z1][C;D2;x1;z1][O;D2;x0]1``
     - 1,3-dioxane, revealing a ketone or aldehyde
   * - ``carbonyl_dioxolane``
     - ``protective:50``
     - ``carbonyl``
     - ``[C;D3,D4;z1;x2;r5:1]1[O;D2;x0][C;D2;x1;z1][C;D2;x1;z1][O;D2;x0]1``
     - 1,3-dioxolane, revealing a ketone or aldehyde
   * - ``carbonyl_dithiane``
     - ``protective:48``
     - ``carbonyl``
     - ``[C;D3,D4;z1;x2;r6:1]1[S;D2;x0;z1][C;D2;x1;z1][C;D2;x0;z1][C;D2;x1;z1][S;D2;x0;z1]1``
     - 1,3-dithiane, revealing a ketone or aldehyde
   * - ``carbonyl_dithiolane``
     - ``protective:47``
     - ``carbonyl``
     - ``[C;D3,D4;z1;x2;r5:1]1[S;D2;x0;z1][C;D2;x1;z1][C;D2;x1;z1][S;D2;x0;z1]1``
     - 1,3-dithiolane, revealing a ketone or aldehyde
   * - ``carboxyl_allyl``
     - ``protective:106``
     - ``carboxyl``
     - ``[C;x2:1](=[O:2])-[O;D2:3]-;!@[C;D2;z1;x1][C;D2;x0;z2]=[C;D1]``
     - allyl ester, revealing a carboxylic acid
   * - ``carboxyl_benzyl``
     - ``protective:105``
     - ``carboxyl``
     - ``[C;x2:1](=[O:2])-[O;D2:3]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - benzyl ester, revealing a carboxylic acid
   * - ``carboxyl_tbu``
     - ``protective:104``
     - ``carboxyl``
     - ``[C;x2:1](=[O:2])-[O;D2:3]-;!@[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]``
     - tert-butyl ester, revealing a carboxylic acid
   * - ``carboxyl_trioxabicyclooctane``
     - ``protective:53``
     - ``carboxyl``
     - ``[C;D4;x3;r6:1]12-;@[O;D2][C;D2;x1;z1][C;D4;x0;z1]([C;D1])([C;D2;x1;z1][O;D2]1)[C;D2;x1;z1][O;D2]2``
     - OBO ester (2,6,7-trioxabicyclo[2.2.2]octane), revealing a carboxylic acid
   * - ``diol_12_acetone``
     - ``protective:34``
     - ``diol``
     - ``[O;D2;x0;r5:1]1-;@[C;D4;x2;z1]([C;D1])([C;D1])-[O;D2;x0:2][C:3]-[C:4]1``
     - acetonide (isopropylidene acetal), revealing a 1,2-diol
   * - ``diol_12_benzylidene``
     - ``protective:45``
     - ``diol``
     - ``[O;D2;x0;r5:1]1-;@[C;D3;x2;z1]([C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)-[O;D2;x0:2][C:3]-[C:4]1``
     - benzylidene acetal, revealing a 1,2-diol
   * - ``diol_12_cyclohexanone``
     - ``protective:41``
     - ``diol``
     - ``[O;D2;x0;r5:1]1-;@[C;D4;x2;z1]2([C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1]2)-[O;D2;x0:2][C:3]-[C:4]1``
     - cyclohexylidene acetal, revealing a 1,2-diol
   * - ``diol_12_cyclopentanone``
     - ``protective:39``
     - ``diol``
     - ``[O;D2;x0;r5:1]1-;@[C;D4;x2;z1]2([C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1]2)-[O;D2;x0:2][C:3]-[C:4]1``
     - cyclopentylidene acetal, revealing a 1,2-diol
   * - ``diol_12_diacetal``
     - ``protective:43``
     - ``diol``
     - ``[O;D2;x0;r6:1]1-;@[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])-[O;D2;x0:2][C:3]-[C:4]1``
     - butane-2,3-diacetal (BDA), revealing a 1,2-diol
   * - ``diol_12_formalin``
     - ``protective:37``
     - ``diol``
     - ``[O;D2;x0;r5:1]1-;@[C;D2;x2;z1]-[O;D2;x0:2][C:3]-[C:4]1``
     - methylene acetal, revealing a 1,2-diol
   * - ``diol_13_acetone``
     - ``protective:35``
     - ``diol``
     - ``[O;D2;x0;r6:1]1-;@[C;D4;x2;z1]([C;D1])([C;D1])-[O;D2;x0:2][C:3][C:4]-[C:5]1``
     - acetonide (isopropylidene acetal), revealing a 1,3-diol
   * - ``diol_13_benzylidene``
     - ``protective:46``
     - ``diol``
     - ``[O;D2;x0;r6:1]1-;@[C;D3;x2;z1]([C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)-[O;D2;x0:2][C:3][C:4]-[C:5]1``
     - benzylidene acetal, revealing a 1,3-diol
   * - ``diol_13_cyclohexanone``
     - ``protective:42``
     - ``diol``
     - ``[O;D2;x0;r6:1]1-;@[C;D4;x2;z1]2([C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1]2)-[O;D2;x0:2][C:3][C:4]-[C:5]1``
     - cyclohexylidene acetal, revealing a 1,3-diol
   * - ``diol_13_cyclopentanone``
     - ``protective:40``
     - ``diol``
     - ``[O;D2;x0;r6:1]1-;@[C;D4;x2;z1]2([C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1]2)-[O;D2;x0:2][C:3][C:4]-[C:5]1``
     - cyclopentylidene acetal, revealing a 1,3-diol
   * - ``diol_13_diacetal``
     - ``protective:44``
     - ``diol``
     - ``[O;D2;x0;r7:1]1-;@[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])-[O;D2;x0:2][C:3][C:4]-[C:5]1``
     - butane-2,3-diacetal (BDA), revealing a 1,3-diol
   * - ``diol_13_formalin``
     - ``protective:38``
     - ``diol``
     - ``[O;D2;x0;r6:1]1-;@[C;D2;x2;z1]-[O;D2;x0:2][C:3][C:4]-[C:5]1``
     - methylene acetal, revealing a 1,3-diol
   * - ``hydroxyl_acyl``
     - ``protective:101``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;z2;x2](=O)-[C;D1]``
     - acetyl, revealing an alcohol
   * - ``hydroxyl_alloc``
     - ``protective:5``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;z1;x1][C;D2;x0;z2]=[C;D1]``
     - allyloxycarbonyl, revealing an alcohol
   * - ``hydroxyl_allyl``
     - ``protective:99``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;z1;x1][C;D2;x0;z2]=[C;D1]``
     - allyl, revealing an alcohol
   * - ``hydroxyl_amine_acetone``
     - ``protective:36``
     - ``hydroxyl``, ``amine``
     - ``[O;D2;x0;r5:1]1-;@[C;D4;x2;z1]([C;D1])([C;D1])-[N;z1:2][C:3]-[C:4]1``
     - acetonide across the amino alcohol, revealing an amino alcohol
   * - ``hydroxyl_benzoate``
     - ``protective:18``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;z2;x2](=O)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - benzoyl, revealing an alcohol
   * - ``hydroxyl_benzyl``
     - ``protective:100``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - benzyl, revealing an alcohol
   * - ``hydroxyl_boc``
     - ``protective:97``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]``
     - tert-butoxycarbonyl, revealing an alcohol
   * - ``hydroxyl_bom``
     - ``protective:15``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - benzyloxymethyl, revealing an alcohol
   * - ``hydroxyl_chloro_tritil``
     - ``protective:28``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D3;x1]([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - 2-chlorotrityl, revealing an alcohol
   * - ``hydroxyl_dimethoxybenzyl``
     - ``protective:13``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1``
     - 2,4-dimethoxybenzyl (DMB), revealing an alcohol
   * - ``hydroxyl_dimetoxy_tritil``
     - ``protective:27``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - 4,4′-dimethoxytrityl (DMT), revealing an alcohol
   * - ``hydroxyl_dmab_enamine``
     - ``protective:32``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;x1;z1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1](:[C;D2]:[C;D2]:1)-[N;D2;x0;z1]-[C;z2;x1;D3]([C;D2;x0;z1][C;D3;x0;z1]([C;D1])[C;D1])=[C;D3;r6;x0;z2]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O``
     - Dmab, as the enamine tautomer, revealing an alcohol
   * - ``hydroxyl_dmab_imine``
     - ``protective:33``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;x1;z1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1](:[C;D2]:[C;D2]:1)-[N;D2;x0;z2]=[C;x1;D3]([C;D2;x0;z1][C;D3;x0;z1]([C;D1])[C;D1])-[C;D3;r6;x0;z1]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O``
     - Dmab, as the imine tautomer, revealing an alcohol
   * - ``hydroxyl_ee``
     - ``protective:23``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D3;x2;z1]([O;D2;x0][C;D2;x1;z1][C;D1])[C;D1]``
     - 1-ethoxyethyl, revealing an alcohol
   * - ``hydroxyl_ethyl``
     - ``protective:103``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;z1;x1]-[C;D1]``
     - ethyl, revealing an alcohol
   * - ``hydroxyl_fmoc``
     - ``protective:2``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;x1;z1][C;D3;z1;x0;r5]1[C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D3]:2-[C;a;r6]:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C1:3``
     - 9-fluorenylmethyloxycarbonyl, revealing an alcohol
   * - ``hydroxyl_mem``
     - ``protective:21``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D2;z1;x1][C;D2;z1;x1][O;D2;x0][C;D1]``
     - 2-methoxyethoxymethyl, revealing an alcohol
   * - ``hydroxyl_methoxy_benzoate``
     - ``protective:17``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;z2;x2](=O)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1``
     - 4-methoxybenzoyl (anisoyl), revealing an alcohol
   * - ``hydroxyl_methoxy_benzyl``
     - ``protective:12``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1``
     - 4-methoxybenzyl (PMB), revealing an alcohol
   * - ``hydroxyl_methyl``
     - ``protective:102``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D1]``
     - methyl, revealing an alcohol
   * - ``hydroxyl_mmt``
     - ``protective:29``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - 4-methoxytrityl (MMT), revealing an alcohol
   * - ``hydroxyl_mom``
     - ``protective:20``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D1]``
     - methoxymethyl, revealing an alcohol
   * - ``hydroxyl_mop``
     - ``protective:24``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])[C;D1]``
     - 2-methoxyprop-2-yl, revealing an alcohol
   * - ``hydroxyl_mpe``
     - ``protective:30``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D4;x1;z1]([C;D1])([C;D2;x0;z1][C;D1])[C;D2;x0;z1][C;D1]``
     - 4-methoxyphenylethyl, revealing an alcohol
   * - ``hydroxyl_naphthyl``
     - ``protective:14``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D3]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D3]:2:[C;D2]:[C;D2]:1``
     - 2-naphthylmethyl, revealing an alcohol
   * - ``hydroxyl_o_nitrobenzyl``
     - ``protective:11``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D3;x1]([N+](=O)[O-]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - ortho-nitrobenzyl, revealing an alcohol
   * - ``hydroxyl_piv``
     - ``protective:16``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;z2;x2](=O)-[C;D4;x0;z1]([C;D1])([C;D1])[C;D1]``
     - pivaloyl, revealing an alcohol
   * - ``hydroxyl_sem``
     - ``protective:25``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D2;z1;x1][C;D2;z1;x1][Si;D4;z1;x0]([C;D1])([C;D1])[C;D1]``
     - 2-(trimethylsilyl)ethoxymethyl, revealing an alcohol
   * - ``hydroxyl_tbdps``
     - ``protective:10``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[Si;D4;z1;x1]([C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)([C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]``
     - tert-butyldiphenylsilyl, revealing an alcohol
   * - ``hydroxyl_tbs``
     - ``protective:8``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[Si;D4;z1;x1]([C;D1])([C;D1])[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]``
     - tert-butyldimethylsilyl, revealing an alcohol
   * - ``hydroxyl_tbu``
     - ``protective:98``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]``
     - tert-butyl, revealing an alcohol
   * - ``hydroxyl_teoc``
     - ``protective:4``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;z1;x1][C;D2;x1;z1][Si;D4;z1;x0]([C;D1])([C;D1])[C;D1]``
     - 2-(trimethylsilyl)ethoxycarbonyl, revealing an alcohol
   * - ``hydroxyl_tes``
     - ``protective:7``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[Si;D4;z1;x1]([C;D2;x1;z1][C;D1])([C;D2;x1;z1][C;D1])[C;D2;x1;z1][C;D1]``
     - triethylsilyl, revealing an alcohol
   * - ``hydroxyl_tfa``
     - ``protective:19``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;z2;x2](=O)-[C;D4;z1;x3](F)(F)F``
     - trifluoroacetyl, revealing an alcohol
   * - ``hydroxyl_thiocarbamate``
     - ``protective:1``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;x3;z2](=[S;D1])[N;D3;x0]([C;D1])[C;D1]``
     - N,N-dimethylthiocarbamate, revealing an alcohol
   * - ``hydroxyl_thp``
     - ``protective:22``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D3;x2;z1;r6]1[O;D2][C;D2][C;D2][C;D2][C;D2]1``
     - tetrahydropyranyl, revealing an alcohol
   * - ``hydroxyl_tips``
     - ``protective:9``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[Si;D4;z1;x1]([C;D3;z1;x1]([C;D1])[C;D1])([C;D3;z1;x1]([C;D1])[C;D1])[C;D3;z1;x1]([C;D1])[C;D1]``
     - triisopropylsilyl, revealing an alcohol
   * - ``hydroxyl_tms``
     - ``protective:6``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[Si;D4;z1;x1]([C;D1])([C;D1])[C;D1]``
     - trimethylsilyl, revealing an alcohol
   * - ``hydroxyl_trifluoroethyl``
     - ``protective:31``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D2;x1;z1][C;D4;x3;z1](F)(F)F``
     - 2,2,2-trifluoroethyl, revealing an alcohol
   * - ``hydroxyl_tritil``
     - ``protective:26``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - triphenylmethyl (trityl), revealing an alcohol
   * - ``hydroxyl_troc``
     - ``protective:3``
     - ``hydroxyl``
     - ``[O;D2:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2][C;D4;x3]([Cl;D1])([Cl;D1])[Cl;D1]``
     - 2,2,2-trichloroethoxycarbonyl, revealing an alcohol
   * - ``thiol_amine_dimethoxybenzyl``
     - ``protective:96``
     - ``thiol``, ``amine``
     - ``[S;D2;r5;x0;z1:1]1[C:3][C:4][N;z1:2]-[C;D3;x2;z1]1-[C;a;r6]:1:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1``
     - 2,4-dimethoxybenzylidene across the aminothiol, revealing an aminothiol
   * - ``thiol_benzyl``
     - ``protective:92``
     - ``thiol``
     - ``[S;D2;x0;z1:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - benzyl, revealing a thiol
   * - ``thiol_chloro_tritil``
     - ``protective:91``
     - ``thiol``
     - ``[S;D2;x0;z1:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D3;x1]([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - 2-chlorotrityl, revealing a thiol
   * - ``thiol_dimetoxy_tritil``
     - ``protective:90``
     - ``thiol``
     - ``[S;D2;x0;z1:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - 4,4′-dimethoxytrityl (DMT), revealing a thiol
   * - ``thiol_mmt``
     - ``protective:89``
     - ``thiol``
     - ``[S;D2;x0;z1:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - 4-methoxytrityl (MMT), revealing a thiol
   * - ``thiol_stbu``
     - ``protective:94``
     - ``thiol``
     - ``[S;D2;x1;z1:1]-;!@[S;D2;z1;x1]-[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]``
     - tert-butyl disulfide, revealing a thiol
   * - ``thiol_strimethoxyphenyl``
     - ``protective:95``
     - ``thiol``
     - ``[S;D2;x1;z1:1]-;!@[S;D2;z1;x1]-[C;a;r6]:1:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]:1[O;D2;x0][C;D1]``
     - 2,4,6-trimethoxyphenyl disulfide, revealing a thiol
   * - ``thiol_tbu``
     - ``protective:93``
     - ``thiol``
     - ``[S;D2;x0;z1:1]-;!@[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]``
     - tert-butyl, revealing a thiol
   * - ``thiol_tritil``
     - ``protective:88``
     - ``thiol``
     - ``[S;D2;x0;z1:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1``
     - triphenylmethyl (trityl), revealing a thiol

.. END GENERATED GLOSSARY
