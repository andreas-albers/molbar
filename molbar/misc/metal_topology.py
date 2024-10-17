import numpy as np


class MetalTopo:

    def __init__(self) -> None:

        pass

    def detach_ligands(self):

        is_metal = self.return_node_data(attribute="is_metal")

        metal_indices = np.where(is_metal == 1)[0]

        visible_nodes = self.return_visible_nodes()

        metal_nodes = np.array(visible_nodes)[metal_indices]

        for metal_node in metal_nodes:

            self._detach_formate(metal_node=metal_node, visible_nodes=visible_nodes)

    def _detach_formate(self, metal_node, visible_nodes):

        self.set_nodes_visible(visible_nodes=visible_nodes)

        adjacent_metal_nodes = self.return_adjacent_nodes(core_nodes=metal_node)

        self.set_nodes_visible(visible_nodes=[metal_node] + adjacent_metal_nodes)

        complex_nodes = self.return_visible_nodes()

        elements = self.return_node_data(attribute="elements")

        carbon_indices = np.where(elements == "C")[0]

        carbon_nodes = np.array(complex_nodes)[carbon_indices]

        for carbon_node in carbon_nodes:

            self._remove_edges_in_formate(
                metal_node=metal_node, carbon_node=carbon_node
            )

    def _remove_edges_in_formate(self, metal_node, carbon_node):

        adjacent_carbon_nodes = self.return_adjacent_nodes(core_nodes=carbon_node)

        self.set_nodes_visible(visible_nodes=adjacent_carbon_nodes)

        ligand_nodes = self.return_visible_nodes()

        is_metal = self.return_node_data(attribute="is_metal")

        ligand_indices = np.where(is_metal == 0)[0]

        ligand_nodes = np.array(ligand_nodes)[ligand_indices]

        if len(ligand_nodes) != 2:

            return

        for ligand_node in ligand_nodes:

            if ligand_node > metal_node:

                self.remove_edge(metal_node, ligand_node)

            else:

                self.remove_edge(ligand_node, metal_node)

        if carbon_node > metal_node:

            self.remove_edge(metal_node, carbon_node)

        else:

            self.remove_edge(carbon_node, metal_node)
