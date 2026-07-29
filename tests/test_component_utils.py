from reneo_utils.component_utils import get_components


class AssemblyGraph:
    def __init__(self, components):
        self._components = components

    def components(self):
        return self._components


def base_kwargs(components):
    return {
        "assembly_graph": AssemblyGraph(components),
        "unitig_names": {0: "edge_1", 1: "edge_2", 2: "edge_3", 3: "edge_4"},
        "unitig_vogs": {
            "edge_1": {"VOG_RDRP"},
            "edge_2": {"VOG_CAPSID"},
            "edge_3": {"VOG_OTHER"},
        },
        "vog_dict": {
            "VOG_RDRP": "RNA-dependent RNA polymerase",
            "VOG_CAPSID": "major capsid protein",
            "VOG_OTHER": "uncharacterized enzyme",
        },
        "smg_unitigs": set(),
        "circular": {"edge_3": 7000},
        "edges_lengths": {"edge_1": 4000, "edge_2": 6000, "edge_3": 7000},
        "nvogs": 1,
        "minlength": 5000,
    }


def test_get_components_keeps_multinode_component_with_viral_vog():
    pruned_vs, comp_vogs = get_components(**base_kwargs([[0, 1], [3]]))

    assert pruned_vs == {0: [0, 1]}
    assert comp_vogs == {0: {"VOG_RDRP", "VOG_CAPSID"}}


def test_get_components_skips_multinode_component_when_smg_is_present():
    kwargs = base_kwargs([[0, 1]])
    kwargs["smg_unitigs"] = {"edge_1"}

    assert get_components(**kwargs) == ({}, {})


def test_get_components_keeps_singleton_only_when_circular_vog_rich_and_long():
    pruned_vs, comp_vogs = get_components(**base_kwargs([[2], [0]]))

    assert pruned_vs == {0: [2]}
    assert comp_vogs == {0: {"VOG_OTHER"}}
