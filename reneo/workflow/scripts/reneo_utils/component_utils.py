def get_components(**kwargs):
    """
    Get connected components with VOGs and no SMGs.
    """

    pruned_vs = {}

    i = 0

    comp_vogs = {}

    for component in kwargs["assembly_graph"].components():
        vogs_found = set()
        vogs_present = False

        if len(component) > 1:
            for unitig in component:
                if kwargs["unitig_names"][unitig] in kwargs["smg_unitigs"]:
                    break
                elif kwargs["unitig_names"][unitig] in kwargs["unitig_vogs"]:
                    for vog in kwargs["unitig_vogs"][kwargs["unitig_names"][unitig]]:
                        vogs_found.add(vog)
                        vogs_present = True

            if vogs_present:
                pruned_vs[i] = component
                comp_vogs[i] = vogs_found
                i += 1

        if len(component) == 1:
            unitig = component[0]
            vogs_present = False

            if kwargs["unitig_names"][unitig] in kwargs["unitig_vogs"]:
                for vog in kwargs["unitig_vogs"][kwargs["unitig_names"][unitig]]:
                    vogs_found.add(vog)
                    vogs_present = True

            if vogs_present:
                vogs_present = True

            if (
                vogs_present
                and len(vogs_found) >= kwargs['nvogs']
                and kwargs["edges_lengths"][kwargs["unitig_names"][unitig]] > kwargs["minlength"]
            ):
                pruned_vs[i] = component
                comp_vogs[i] = vogs_found
                i += 1

    return pruned_vs, comp_vogs
