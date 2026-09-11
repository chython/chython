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

249 functional groups, alphabetically.  The name is the key
:meth:`chython.MoleculeContainer.functional_groups` returns and the key
:func:`chython.functional_rules` is keyed on; the id is what a consumer stores.

.. list-table::
   :header-rows: 1
   :widths: 22 12 30 36

   * - Name
     - Id
     - SMARTS
     - What it matches
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
   * - ``activated_isocyanide``
     - ``functional:166``
     - ``[C;-;D1:1]#[N;+;D2:2]-[C;z1;h1,h2:3]``
     - isocyanide with an acidic alpha C-H; what activates it is left unstated
   * - ``active_methylene``
     - ``functional:167``
     - ``[C;z1;D2,D3;x0:1](-[C;z2,z3;x1,x2:2])-[C;z2,z3;x1,x2:3]``
     - CH flanked by two EWGs, the Knoevenagel nucleophile
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
   * - ``aldehyde``
     - ``functional:20``
     - ``[O;z2;x0:1]=[C;D2;x1;z2:2]``
     - R-CHO
   * - ``alkene``
     - ``functional:2``
     - ``[C;z2;x0;D2,D3:1]=[C;z2;x0;D2,D3:2]``
     - internal C=C, no heteroatom on either carbon
   * - ``alkenyl_boronic_acid``
     - ``functional:62``
     - ``[B;D3;z1;x2](-[O;D1])(-[O;D1])-;!@[C;z2;x1:1]=[C:2]``
     - vinylboronic acid
   * - ``alkenyl_boronic_ester``
     - ``functional:63``
     - ``[B;D3;z1;x2](-[O;D2;x1])(-[O;D2;x1])-;!@[C;z2;x1:1]=[C:2]``
     - vinylboronate ester
   * - ``alkenyl_bromide``
     - ``functional:49``
     - ``[Br;D1][C;z2;x1:1]=[C:2]``
     - vinylic C-Br
   * - ``alkenyl_carboxylic_acid``
     - ``functional:81``
     - ``[O;D1;z1;x0]-[C;z2;x2;D3:1](=[O:2])-[C;z2:3]=[C:4]``
     - alpha,beta-unsaturated acid
   * - ``alkenyl_chloride``
     - ``functional:48``
     - ``[Cl;D1][C;z2;x1:1]=[C:2]``
     - vinylic C-Cl, the Heck/Negishi electrophile
   * - ``alkenyl_fluoride``
     - ``functional:47``
     - ``[F;D1][C;z2;x1:1]=[C:2]``
     - vinylic C-F
   * - ``alkenyl_grignard``
     - ``functional:112``
     - ``[Mg;D2](-[F,Cl,Br,I])-[C;z2:1]=[C:2]``
     - vinyl Grignard
   * - ``alkenyl_iodide``
     - ``functional:50``
     - ``[I;D1][C;z2;x1:1]=[C:2]``
     - vinylic C-I
   * - ``alkenyl_molander_salt``
     - ``functional:68``
     - ``[B;D4;z1;x3;-](-[F])(-[F])(-[F])-;!@[C;z2:1]=[C:2]``
     - vinyl trifluoroborate
   * - ``alkenyl_silane``
     - ``functional:123``
     - ``[Si;D4]-;!@[C;z2:1]=[C:2]``
     - vinylsilane
   * - ``alkenyl_stannane``
     - ``functional:120``
     - ``[Sn;D4;z1]-;!@[C;z2:1]=[C:2]``
     - vinylstannane
   * - ``alkenyl_zinc``
     - ``functional:115``
     - ``[Zn;D2](-[F,Cl,Br,I])-[C;z2:1]=[C:2]``
     - vinylzinc
   * - ``alkoxide``
     - ``functional:232``
     - ``[O;D1;z1;x0;-:1]-[C;z1:2]``
     - RO-
   * - ``alkyl_boronic_acid``
     - ``functional:14``
     - ``[B;D3;z1;x2](-[O;D1])(-[O;D1])-;!@[C;z1;x1:1]``
     - RB(OH)2 on sp3 carbon
   * - ``alkyl_boronic_ester``
     - ``functional:61``
     - ``[B;D3;z1;x2](-[O;D2;x1])(-[O;D2;x1])-;!@[C;z1;x1:1]``
     - RB(OR')2 on sp3 carbon
   * - ``alkyl_bromide``
     - ``functional:10``
     - ``[Br;D1][C;z1;x1:1]``
     - R-Br on sp3 carbon
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
     - ``[Zn;D2](-[F,Cl,Br,I])-[C;z1:1]``
     - RZnX, the Negishi nucleophile
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
   * - ``alpha_haloester``
     - ``functional:74``
     - ``[O;D2;x0:1]-[C;D3;x2;z2:2](=[O:3])-[C;z1;D2,D3;x1:4]-[Cl,Br,I;D1]``
     - alpha-halo ester
   * - ``alpha_haloketone``
     - ``functional:73``
     - ``[O;z2;x0:1]=[C;D3;x1:2]-[C;z1;D2,D3;x1:3]-[Cl,Br,I;D1]``
     - alpha-halo ketone, the Hantzsch thiazole electrophile
   * - ``alpha_heteroatom_halide``
     - ``functional:225``
     - ``[Cl,Br,I;D1:1]-[C;z1;h1,h2:2]-[N,O,S;z1:3]``
     - X-C-Het, the MOM-chloride class: an alkylating agent far more reactive than ``alkyl_halide`` says
   * - ``alpha_ketone``
     - ``functional:72``
     - ``[O;z2;x0:1]=[C;D3;x1:2]-[C;z1;D1,D2;x0:3]``
     - ketone with an enolizable alpha CH2/CH3
   * - ``amidine``
     - ``functional:152``
     - ``[N;D1;z1;x0:1]-[C;D3;z2;x2:2]=[N;D1:3]``
     - RC(=NH)NH2, the pyrimidine partner
   * - ``amidoxime``
     - ``functional:151``
     - ``[N;D1;z2;x0:1]=[C;D3;x2:2]-[N;D2;z1;x1:3]-[O;D1:4]``
     - RC(=NH)NHOH, the 1,2,4-oxadiazole precursor
   * - ``amino_alcohol``
     - ``functional:159``
     - ``[N;D1;z1;x0:1]-[C;z1:2]-[C;z1:3]-[O;D1:4]``
     - 1,2-amino alcohol
   * - ``aminopyridine``
     - ``functional:158``
     - ``[N;D1;z1;x0:1]-[C;a:2]:[N;a;h0;D2:3]``
     - 2-aminoazine, the GBB partner
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
   * - ``aryl_bromide``
     - ``functional:7``
     - ``[Br;D1]-[C;a:1]``
     - Ar-Br
   * - ``aryl_bromide_iodide``
     - ``functional:44``
     - ``[Br,I;D1]-[C;a:1]``
     - Ar-Br or Ar-I, the halides that couple without a specialist ligand
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
   * - ``aryl_mesylate``
     - ``functional:56``
     - ``[S;D4](=[O])(=[O])(-[O;D2]-;!@[C;a:1])-[C;D1]``
     - ArOMs
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
     - ``[Zn;D2](-[F,Cl,Br,I])-[C;a:1]``
     - ArZnX
   * - ``azetidine``
     - ``functional:238``
     - ``[N;D2;z1;x0;r4;h1:1]([C;z1;r4:2])[C;z1;r4:3]``
     - The four-membered NH amine
   * - ``azide``
     - ``functional:102``
     - ``[N;x1;D2:1]=[N;+:2]=[N;-:3]``
     - R-N3
   * - ``azine_n_oxide``
     - ``functional:206``
     - ``[N;a;D3;+:1]-[O;D1;z1;-:2]``
     - The pyridine N-oxide, an activated ring for C-H functionalization.  The oxygen is NUMBERED: a deoxygenation omits it from its own product side, and leaving it unnumbered would delete it from every other template too
   * - ``azinone``
     - ``functional:132``
     - ``[O;z2;x0:1]=[C;z2;x2;D3;r6:2](-[C,N;z2,z4;r6:3])-[N;z1;D2;r6:4]``
     - aromatizable cyclic amide; the sp2/aromatic flank excludes a saturated lactam
   * - ``aziridine_nh``
     - ``functional:94``
     - ``[N;D2;z1;x0;h1;r3:1]``
     - aziridine N-H, matched alone so a ring opening leaves the ring intact
   * - ``azo``
     - ``functional:105``
     - ``[N;D2;x1:1]=[N;D2;x1:2]``
     - R-N=N-R
   * - ``benzyl_halide``
     - ``functional:222``
     - ``[Cl,Br,I;D1:1]-[C;z1;x1;h1,h2:2]-[C;a:3]``
     - Ar-CH2-X.  ``alkyl_halide`` already matches it; this says the SN2 is a benzylic one, which is a different rate and a different selectivity
   * - ``benzylic_ch``
     - ``functional:42``
     - ``[C;z1;h1,h2;D2,D3:1]-[C;a:2]``
     - an sp3 carbon bearing hydrogen next to a ring.  ``D2,D3`` excludes a TOLUENE METHYL, which is a D1 -- so ``benzylic_bromination`` does not reach one
   * - ``beta_arylethylamine``
     - ``functional:157``
     - ``[N;D1;z1;x0:1]-[C;z1:2]-[C;z1:3]-[C;a:4]:[C;a;D2:5]``
     - the Pictet-Spengler substrate
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
   * - ``carbamate``
     - ``functional:191``
     - ``[N:1]-[C;D3;z2;x3:2](=[O;D1;z2;x0:3])-[O;D2;z1:4]-[C:5]``
     - N-CO-O-R, the Boc/Cbz/Fmoc backbone.  Reported wherever it occurs, unlike ``protective.tsv``, which only names the ones it can cleave
   * - ``carbamoyl_chloride``
     - ``functional:91``
     - ``[Cl;D1]-[C;z2;x3;D3:1](=[O:2])-[N;D2,D3:3]``
     - R2N-CO-Cl
   * - ``carbamoyl_fluoride``
     - ``functional:92``
     - ``[F;D1]-[C;z2;x3;D3:1](=[O:2])-[N;D2,D3:3]``
     - R2N-CO-F
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
   * - ``difluoromethoxy``
     - ``functional:208``
     - ``[F:1][C;D3;z1;x3:2]([F:3])-[O;D2:4]``
     - OCHF2
   * - ``difluoromethyl``
     - ``functional:173``
     - ``[C;D3;z1;x2:1](-[F:2])-[F:3]``
     - -CHF2
   * - ``disulfide``
     - ``functional:130``
     - ``[S;D2;z1:1]-[S;D2;z1:2]``
     - R-S-S-R
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
   * - ``enone``
     - ``functional:194``
     - ``[C;z2:1]=[C;z2:2]-[C;D3;z2;x1:3]=[O;D1;x0:4]``
     - C=C-C=O, the Michael acceptor.  The conjugation is the group: neither ``ketone`` nor ``alkene`` alone says a molecule is electrophilic at the beta carbon
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
   * - ``guanidine``
     - ``functional:107``
     - ``[N;z1;x0:1]-[C;!R:2](-[N;z1;x0:3])=[N;x0:4]``
     - acyclic guanidine
   * - ``hemiketal``
     - ``functional:187``
     - ``[O;D1;z1;x0;h1:1]-[C;z1;x2:2]-[O;D2;z1:3]``
     - HO-C(-OR), the ring-opened sugar form
   * - ``hydrazide``
     - ``functional:161``
     - ``[N;D1;z1;x1:1]-[N;D2;z1;x1:2]-[C;z2;x2:3]=[O:4]``
     - R-CO-NH-NH2
   * - ``hydrazine``
     - ``functional:215``
     - ``[N;z1;x1:1]-[N;z1;x1:2]``
     - N-N single bond, both sp3.  A hydrazide matches this too and ``hydrazide`` names it more precisely
   * - ``hydrazone``
     - ``functional:142``
     - ``[C;z2:1]=[N;D2;z2;x1:2]-[N;D2;z1;x1:3]``
     - C=N-NH-R
   * - ``hydroxamic_acid``
     - ``functional:160``
     - ``[O;D1;z1;x1:1]-[N;D2;z1;x1:2]-[C;z2;x2:3]=[O:4]``
     - R-CO-NH-OH
   * - ``imidazole``
     - ``functional:136``
     - ``[N;h1;D2;a;r5:1]:[A:2]:[N;h0;D2;r5:3]``
     - N-H and N separated by one ring atom
   * - ``imidazolyl_carbonate``
     - ``functional:90``
     - ``[O;D2:1]-[C;z2;x3;D3:2](=[O:3])-[N;a;D3;r5]1:[C;a]:[N;a]:[C;a]:[C;a]:1``
     - CDI adduct of an alcohol; here the azole nitrogen leaves
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
   * - ``isothiocyanate``
     - ``functional:109``
     - ``[N;z2;x0;D2:1]=[C:2]=[S;D1:3]``
     - R-N=C=S
   * - ``isoxazole``
     - ``functional:137``
     - ``[O;a;D2;r5:1]:[N;a;D2;r5:2]``
     - adjacent aromatic O/N pair
   * - ``ketone``
     - ``functional:21``
     - ``[O;z2;x0:1]=[C;D3;x1;z2:2]``
     - R2C=O
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
   * - ``maleimide``
     - ``functional:177``
     - ``[N;r5:1]1-[C;z2;r5:2](=[O:3])-[C;z2;r5:4]=[C;z2;r5:5]-[C;z2;r5:6]1=[O:7]``
     - maleimide, the thiol-conjugation warhead
   * - ``metalate_carbanion``
     - ``functional:241``
     - ``[Mg,Zn,Cu,Li;D1;+:1]-[C:2]``
     - The ionic spelling of an organometallic, which the neutral ``alkyl_grignard`` and ``alkyl_zinc`` rows cannot match
   * - ``methyl_ester``
     - ``functional:170``
     - ``[O;z2;x0:1]=[C;D3;x2;z2:2]-[O;D2;x0:3]-[C;D1]``
     - R-CO-OMe; distinct saponification and transesterification behaviour
   * - ``n_halo_imide``
     - ``functional:227``
     - ``[Cl,Br,I;D1:1]-[N;D3:2]``
     - N-X, the NBS/NCS halogen source
   * - ``n_hydroxylamine``
     - ``functional:216``
     - ``[O;D1;z1;x1;h1:1]-[N;z1;x1:2]-[C;z1,a:3]``
     - HO-N with carbon on the nitrogen.  The carbon is stated to exclude the hydroxamic acid, whose nitrogen carries a carbonyl
   * - ``n_substituted_azole``
     - ``functional:205``
     - ``[N;a;D3:1]-[C;z1,a:2]``
     - The substituted pyrrole-type nitrogen -- an N-alkylated or N-arylated azole, as distinct from the ``nh_azole`` that could still be alkylated
   * - ``nh_thiourea``
     - ``functional:156``
     - ``[N;z1;h1,h2;!R:1]-[C;D3;z2;x3:2](=[S;D1:3])-[N;z1;h1,h2;!R:4]``
     - the thiourea half of the same
   * - ``nh_urea``
     - ``functional:155``
     - ``[N;z1;h1,h2;!R:1]-[C;D3;z2;x3:2](=[O:3])-[N;z1;h1,h2;!R:4]``
     - acyclic urea with an N-H on both nitrogens, the Biginelli subset
   * - ``nitrile``
     - ``functional:32``
     - ``[N;D1;z3;x0:1]#[C;D2;x1:2]``
     - R-CN; ``z3`` here is a genuine sp nitrogen
   * - ``nitro``
     - ``functional:33``
     - ``[N;D3;x2;+:1]([O;-:2])=[O:3]``
     - R-NO2, charge-separated as chython stores it
   * - ``nitroso``
     - ``functional:108``
     - ``[N;D2;z2:1]=[O;D1:2]``
     - R-N=O
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
   * - ``oxetane``
     - ``functional:237``
     - ``[O;D2;z1;r4:1]([C;z1;r4:2])[C;z1;r4:3]``
     - The four-membered ether, a carbonyl bioisostere and a slow electrophile
   * - ``oxime``
     - ``functional:162``
     - ``[O;D1;z1;x1:1]-[N;D2;z2;x1:2]=[C:3]``
     - C=N-OH
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
   * - ``pyrazole``
     - ``functional:135``
     - ``[N;h1;D2;a;r5:1]:[N;h0;D2;r5:2]``
     - adjacent N-H/N pair
   * - ``pyridazine``
     - ``functional:138``
     - ``[N;a;D2;r6:1]:[N;a;D2;r6:2]``
     - adjacent aromatic N pair in a 6-ring
   * - ``pyridine_n``
     - ``functional:40``
     - ``[N;a;D2;h0:1]``
     - an aromatic nitrogen with no hydrogen, the N-oxidation substrate
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
   * - ``silyl_ether``
     - ``functional:189``
     - ``[Si;D4:1]-[O;D2;z1:2]-[C:3]``
     - R3Si-O-R, a protected alcohol.  ``protective.tsv`` removes the named ones; this reports any
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
     - ``[C;z3;x0;D1:1]#[C;x0;D2:2]``
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
   * - ``tosyl_isocyanide``
     - ``functional:165``
     - ``[C;-;D1:1]#[N;+;D2:2]-[C;D2,D3;z1;x2:3]-[S;D4;x2](=[O])=[O]``
     - TosMIC, the Van Leusen oxazole reagent
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
   * - ``vinyl_sulfone``
     - ``functional:175``
     - ``[S;D4:1](=[O:2])(=[O:3])(-[O,N,C:4])-[C;z2:5]=[C:6]``
     - vinyl sulfone, sulfonamide or sulfonate -- a covalent warhead.  chython 2 wrote the third substituent ``[O,N]`` and so missed every true sulfone; widened here, which is why the ``x3`` it carried had to go
   * - ``weinreb_amide``
     - ``functional:127``
     - ``[O:1]=[C;D3;x2:2]-[N;D3;x1](-[C;z1])-[O;D2;x1]``
     - N-methoxy-N-methyl amide; the extra C on N excludes an N-acyloxy imide
   * - ``xanthate``
     - ``functional:220``
     - ``[S;D1,D2:1]-[C;D3;z2:2](=[S;D1;x0:3])-[O;D2:4]``
     - S-C(=S)-O-R, the RAFT agent.  ``D1,D2`` because the second sulfur is a free thiol as often as it is substituted

Protecting groups
-----------------

103 protecting groups, alphabetically.  *Protects* names the functional
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
