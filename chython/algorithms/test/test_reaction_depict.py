from chython import smarts


REACTION_SMARTS_WITH_SINGLE_ATOM_FRAGMENTS = (
    '[C:1]-[C:2](=[O:6])-[O:7]>>[C:1]-[C:2](-[O:3]-[C:4]-[C:5])=[O:6].[O-:7]',
    '[C:3](=[O:4])(-[C:5])-[O:6]>>[O-:6].[C:1]-[O:2]-[C:3](=[O:4])-[C:5]',
    '[c:1]-[C:2](=[O:6])-[O:7]>>[O-:7].[c:1]-[C:2](-[O:3]-[C:4]-[C:5])=[O:6]',
    '[c:1]-[C:2](=[O:4])-[N:5]>>[N+:5].[c:1]-[C:2](-[O:3])=[O:4]',
    '[c:1]-[C:2](=[O:3])-[N:4]>>[c:1]-[C:2](=[O:3])-[O:5].[N:4]',
    '[c:2](:[c:3])(:[n:4])-[N:5]>>[N:5].[Cl:1]-[c:2](:[c:3]):[n:4]',
    '[C:1]-[C:2](=[O:4])-[N:5]>>[C:1]-[C:2](-[O:3])=[O:4].[N+:5]',
    '[C:1]-[C:2](=[O:4])-[N:5]>>[N:5].[C:1]-[C:2](-[O:3])=[O:4]',
    '[c:1]-[C:2](=[O:5])-[O:6]>>[c:1]-[C:2](-[O:3]-[C:4])=[O:5].[O-:6]',
    '[C:1]-[O:5].[c:2]-[C:3](=[O:4])-[N:6]>>[N:6].[C:1]-[O:5]-[C:3](-[c:2])=[O:4]',
    '[C:1]-[C:2](=[O:3])-[N:4]>>[C:1]-[C:2](=[O:3])-[O:5]-[C:6].[N:4]',
    '[c:1](:[c:2])(:[c:3])-[I:4]>>[I-:4].[c:1](:[c:2])(:[c:3])-[N:5]',
    '[c:1]-[C:2](=[O:6])-[N:7]>>[N:7].[c:1]-[C:2](-[O:3]-[C:4]-[C:5])=[O:6]',
    '[c:1]-[C:2](-[N:3])=[O:4]>>[c:1]-[C:2]#[N:3].[O-:4]',
    '[C:3](-[C:4])(=[O:5])-[N:7]>>[N:7].[C:1](-[O:2]-[C:3](-[C:4])=[O:5])-[C:6]',
)


def _seed_multatom_query_planes(reaction):
    for graph in reaction.molecules():
        if len(graph) <= 1 or not hasattr(graph, '_plane'):
            continue
        for i, n in enumerate(graph):
            xy = (float(i), 0.)
            graph._plane[n] = xy
            if hasattr(graph._atoms[n], 'xy'):
                graph._atoms[n].xy = xy


def test_reaction_depict_handles_single_atom_query_fragments_without_plane():
    for pattern in REACTION_SMARTS_WITH_SINGLE_ATOM_FRAGMENTS:
        reaction = smarts(pattern)
        _seed_multatom_query_planes(reaction)

        svg = reaction.depict(clean2d=False)

        assert svg.startswith('<svg')
