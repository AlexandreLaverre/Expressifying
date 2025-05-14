import argparse
from ete3 import Tree
from neutrality_index import open_tree, prune_tree


def rename_tree(tree: Tree) -> Tree:
    for node in tree.traverse("levelorder"):
        if node.is_leaf():
            if node.name.endswith("_E"):
                node.name = "_".join(node.name.replace("_E", "").split("_")[:-1])
            elif "PD_" in node.name[:3]:
                node.name = "_".join(node.name.split("_")[2:])
                assert node.name.count("_") == 1
            else:
                node.name = node.name.replace(" ", "_")
    return tree


def main(primate_tree, mammal_tree, output_tree):
    """
    Prune the mammal tree to keep only the primate species
    :param primate_tree:
    :param mammal_tree:
    :param output_tree:
    :return:
    """
    t_primate = rename_tree(open_tree(primate_tree, format_ete3=1))
    discard_list = ["Galeopterus_variegatus", "Tupaia_belangeri", "Oryctolagus_cuniculus", "Mus_musculus"]
    keep_species = list(set(t_primate.get_leaf_names()).difference(set(discard_list)))
    t_primate = prune_tree(t_primate, keep_species)
    t_mammal = rename_tree(open_tree(mammal_tree, format_ete3=1))
    difference = list(set(t_mammal.get_leaf_names()) - set(t_primate.get_leaf_names()))
    print(f"Difference between mammals and primates has {len(difference)} leaves")
    assert len(difference) > 0, "No difference between primates and mammals"
    assert len(difference) < len(t_mammal.get_leaf_names()), "No primates in mammals"
    t_mammal = prune_tree(t_mammal, difference)
    assert len(set(t_mammal.get_leaf_names())) == len(difference)
    t_mammal.write(outfile=output_tree, format=1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--primate_tree", help="Input primate tree file", required=True)
    parser.add_argument("--mammal_tree", help="Input mammal neutral tree file", required=True)
    parser.add_argument("--output_tree", help="Output non primate tree file", required=True)
    args = parser.parse_args()
    main(args.primate_tree, args.mammal_tree, args.output_tree)
