#!/usr/bin/env python3

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2023, Reneo Project"
__license__ = "MIT"
__version__ = "0.6.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Development"


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

        vog_functions = {
            "rdrp": 0,
            "reverse_transcriptase": 0,
            "integrase": 0,
            "dna_polymerase": 0,
            "major_capsid_protein": 0,
            "rna_polymerase": 0,
            "terl": 0,
            "portal_protein": 0,
            "tail": 0,
            "replication_associated_protein": 0,
            "helicase": 0,
            "protease": 0,
            "hypothetical": 0,
            "other_func": 0,
        }

        if len(component) > 1:
            for unitig in component:
                if kwargs["unitig_names"][unitig] in kwargs["smg_unitigs"]:
                    break
                elif kwargs["unitig_names"][unitig] in kwargs["unitig_vogs"]:
                    for vog in kwargs["unitig_vogs"][kwargs["unitig_names"][unitig]]:
                        
                        if "rna-dependent rna polymerase" in kwargs["vog_dict"][vog].lower():
                            vog_functions["rdrp"] += 1
                        elif "reverse transcriptase" in kwargs["vog_dict"][vog].lower():
                            vog_functions["reverse_transcriptase"] += 1
                        elif "integrase" in kwargs["vog_dict"][vog].lower():
                            vog_functions["integrase"] += 1
                        elif "dna polymerase" in kwargs["vog_dict"][vog].lower():
                            vog_functions["dna_polymerase"] += 1
                        elif "major capsid protein" in kwargs["vog_dict"][vog].lower():
                            vog_functions["major_capsid_protein"] += 1
                        elif "rna polymerase" in kwargs["vog_dict"][vog].lower():
                            vog_functions["rna_polymerase"] += 1
                        elif "terminase large subunit" in kwargs["vog_dict"][vog].lower():
                            vog_functions["terl"] += 1
                        elif "portal protein" in kwargs["vog_dict"][vog].lower():
                            vog_functions["portal_protein"] += 1
                        elif "tail" in kwargs["vog_dict"][vog].lower():
                            vog_functions["tail"] += 1
                        elif "replication associated protein" in kwargs["vog_dict"][vog].lower():
                            vog_functions["replication_associated_protein"] += 1
                        elif "helicase" in kwargs["vog_dict"][vog].lower():
                            vog_functions["helicase"] += 1
                        elif "protease" in kwargs["vog_dict"][vog].lower():
                            vog_functions["protease"] += 1
                        elif "hypothetical protein" in kwargs["vog_dict"][vog].lower():
                            vog_functions["hypothetical"] += 1
                        else:
                            vog_functions["other_func"] += 1
                        
                        vogs_found.add(vog)

            if vog_functions["rdrp"] > 0:
                vogs_present = True
            elif vog_functions["reverse_transcriptase"] > 0 and vog_functions["integrase"] > 0:
                vogs_present = True
            elif vog_functions["dna_polymerase"] > 0 and vog_functions["major_capsid_protein"] > 0 and vog_functions["rna_polymerase"] > 0:
                vogs_present = True
            elif vog_functions["terl"] > 0 and vog_functions["portal_protein"] > 0 and vog_functions["tail"] > 0 and vog_functions["major_capsid_protein"] > 0:
                vogs_present = True
            elif vog_functions["replication_associated_protein"] > 0:
                vogs_present = True
            elif vog_functions["rdrp"] > 0 and vog_functions["helicase"] > 0 and vog_functions["protease"] > 0:
                vogs_present = True
            elif vog_functions["terl"] > 0 and vog_functions["major_capsid_protein"] > 0 and vog_functions["tail"] > 0:
                vogs_present = True
            elif vog_functions["other_func"] > 0:
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

                    if "rna-dependent rna polymerase" in kwargs["vog_dict"][vog].lower():
                        vog_functions["rdrp"] += 1
                    elif "reverse transcriptase" in kwargs["vog_dict"][vog].lower():
                        vog_functions["reverse_transcriptase"] += 1
                    elif "integrase" in kwargs["vog_dict"][vog].lower():
                        vog_functions["integrase"] += 1
                    elif "dna polymerase" in kwargs["vog_dict"][vog].lower():
                        vog_functions["dna_polymerase"] += 1
                    elif "major capsid protein" in kwargs["vog_dict"][vog].lower():
                        vog_functions["major_capsid_protein"] += 1
                    elif "rna polymerase" in kwargs["vog_dict"][vog].lower():
                        vog_functions["rna_polymerase"] += 1
                    elif "terminase large subunit" in kwargs["vog_dict"][vog].lower():
                        vog_functions["terl"] += 1
                    elif "portal protein" in kwargs["vog_dict"][vog].lower():
                        vog_functions["portal_protein"] += 1
                    elif "tail" in kwargs["vog_dict"][vog].lower():
                        vog_functions["tail"] += 1
                    elif "replication associated protein" in kwargs["vog_dict"][vog].lower():
                        vog_functions["replication_associated_protein"] += 1
                    elif "helicase" in kwargs["vog_dict"][vog].lower():
                        vog_functions["helicase"] += 1
                    elif "protease" in kwargs["vog_dict"][vog].lower():
                        vog_functions["protease"] += 1
                    elif "hypothetical protein" in kwargs["vog_dict"][vog].lower():
                        vog_functions["hypothetical"] += 1
                    else:
                        vog_functions["other_func"] += 1

                    vogs_found.add(vog)
                    
            if vog_functions["rdrp"] > 0:
                vogs_present = True
            elif vog_functions["reverse_transcriptase"] > 0 and vog_functions["integrase"] > 0:
                vogs_present = True
            elif vog_functions["dna_polymerase"] > 0 and vog_functions["major_capsid_protein"] > 0 and vog_functions["rna_polymerase"] > 0:
                vogs_present = True
            elif vog_functions["terl"] > 0 and vog_functions["portal_protein"] > 0 and vog_functions["tail"] > 0 and vog_functions["major_capsid_protein"] > 0:
                vogs_present = True
            elif vog_functions["replication_associated_protein"] > 0:
                vogs_present = True
            elif vog_functions["rdrp"] > 0 and vog_functions["helicase"] > 0 and vog_functions["protease"] > 0:
                vogs_present = True
            elif vog_functions["terl"] > 0 and vog_functions["major_capsid_protein"] > 0 and vog_functions["tail"] > 0:
                vogs_present = True
            elif vog_functions["other_func"] > 0:
                vogs_present = True

            # if vogs_present:
            #     vogs_present = True

            if (
                vogs_present
                and kwargs["unitig_names"][unitig] in kwargs["circular"]
                and len(vogs_found) >= kwargs["nvogs"]
                and kwargs["edges_lengths"][kwargs["unitig_names"][unitig]]
                > kwargs["minlength"]
            ):
                pruned_vs[i] = component
                comp_vogs[i] = vogs_found
                i += 1

    return pruned_vs, comp_vogs
