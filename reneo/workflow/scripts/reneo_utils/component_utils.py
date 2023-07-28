def get_components(
    assembly_graph,
    unitig_names,
    smg_unitigs,
    unitig_vogs,
    circular,
    edges_lengths,
    cicular_len,
):
    """
    Get connected components with VOGs and no SMGs.
    """

    pruned_vs = {}

    i = 0

    comp_vogs = {}

    for component in assembly_graph.components():
        vogs_found = set()
        vogs_present = False

        if len(component) > 1:
            for unitig in component:
                if unitig_names[unitig] in smg_unitigs:
                    break
                elif unitig_names[unitig] in unitig_vogs:
                    for vog in unitig_vogs[unitig_names[unitig]]:
                        vogs_found.add(vog)
                        vogs_present = True

            if vogs_present:
                pruned_vs[i] = component
                comp_vogs[i] = vogs_found
                i += 1

        if len(component) == 1:
            unitig = component[0]
            vogs_present = False

            if unitig_names[unitig] in unitig_vogs:
                for vog in unitig_vogs[unitig_names[unitig]]:
                    vogs_found.add(vog)
                    vogs_present = True

            if vogs_present:
                vogs_present = True

            if (
                vogs_present
                # and unitig_names[unitig] in circular
                and edges_lengths[unitig_names[unitig]] > cicular_len
            ):
                pruned_vs[i] = component
                comp_vogs[i] = vogs_found
                i += 1

    return pruned_vs, comp_vogs
