import os
import logging
from typing import List, Dict, Optional, Set
import pandas as pd
import hashlib


class PathwayDatabase:
    """
    Centralized manager for working with KEGG pathway database.

    Loads subsystem_matrix.csv (reactions × pathways) and provides convenient methods
    for retrieving information about reaction membership in metabolic pathways.

    Attributes:
        matrix (pd.DataFrame): Membership matrix (reactions × pathways)
        all_pathways (List[str]): List of all pathways in the database
        all_reactions (List[str]): List of all reactions in the database
    """

    def __init__(self, matrix_path: str = "data/models/subsystem_matrix.csv"):

        if not os.path.exists(matrix_path):
            raise FileNotFoundError(f"KEGG matrix not found: {matrix_path}")

        logging.info(f"Loading KEGG pathway database from {matrix_path}")

        self.matrix = pd.read_csv(matrix_path, index_col=0)

        if self.matrix.empty:
            raise ValueError(f"KEGG matrix is empty: {matrix_path}")

        self.all_reactions = list(self.matrix.index)
        self.all_pathways = list(self.matrix.columns)

        logging.info(
            f"Loaded KEGG database: {len(self.all_reactions)} reactions, "
            f"{len(self.all_pathways)} pathways"
        )

        reactions_with_pathways = (self.matrix.sum(axis=1) > 0).sum()
        avg_pathways_per_reaction = self.matrix.sum(axis=1).mean()

        logging.info(
            f"Coverage: {reactions_with_pathways}/{len(self.all_reactions)} reactions have pathways "
            f"(avg {avg_pathways_per_reaction:.2f} pathways/reaction)"
        )

    def get_pathways_for_reaction(self, rxn_id: str) -> List[str]:
        """
        Return list of metabolic pathways for a reaction.

        """
        if rxn_id not in self.matrix.index:
            return []

        pathway_mask = self.matrix.loc[rxn_id] == 1
        return self.matrix.columns[pathway_mask].tolist()

    def get_reactions_for_pathway(self, pathway: str) -> List[str]:
        """
        Return list of reactions in a metabolic pathway.

        """
        if pathway not in self.matrix.columns:
            pathway_lower = pathway.lower()
            for col in self.matrix.columns:
                if col.lower() == pathway_lower:
                    pathway = col
                    break
            else:
                return []

        rxn_mask = self.matrix[pathway] == 1
        return self.matrix.index[rxn_mask].tolist()

    def get_all_pathways(self, exclude_empty: bool = False) -> List[str]:
        """
        Return list of all metabolic pathways in the database.

        """
        if not exclude_empty:
            return self.all_pathways.copy()

        return [
            pw for pw in self.all_pathways
            if self.matrix[pw].sum() > 0
        ]

    def get_all_reactions(self, with_pathways_only: bool = False) -> List[str]:
        """
        Return list of all reactions in the database.

        """
        if not with_pathways_only:
            return self.all_reactions.copy()

        return [
            rxn for rxn in self.all_reactions
            if self.matrix.loc[rxn].sum() > 0
        ]

    def map_to_unified(
        self,
        kegg_pathway: str,
        unified_mapping: Optional[Dict[str, str]] = None
    ) -> str:
        """
        Map KEGG pathway names to UNIFIED_PATHWAYS (for backward compatibility).

        """
        if unified_mapping is None:
            unified_mapping = DEFAULT_UNIFIED_MAPPING

        return unified_mapping.get(kegg_pathway, kegg_pathway)

    def assign_global_colors(
        self,
        pathways: Optional[List[str]] = None,
        palette: str = "tab20"
    ) -> Dict[str, str]:
        """
        Assign unique colors to metabolic pathways.

        """
        try:
            import matplotlib.pyplot as plt
            import matplotlib.colors as mcolors
        except ImportError:
            logging.warning("matplotlib not available, using fallback colors")
            return self._assign_colors_fallback(pathways)

        if pathways is None:
            pathways = self.get_all_pathways(exclude_empty=True)

        try:
            cmap = plt.get_cmap(palette)
            if palette.startswith('tab'):
                base_colors = [mcolors.rgb2hex(cmap(i)) for i in range(cmap.N)]
            elif palette.startswith('Set'):
                base_colors = [mcolors.rgb2hex(cmap(i)) for i in range(cmap.N)]
            else:
                n_colors = min(len(pathways), 50)
                base_colors = [
                    mcolors.rgb2hex(cmap(i / n_colors))
                    for i in range(n_colors)
                ]
        except Exception as e:
            logging.warning(f"Palette {palette} not available, using default: {e}")
            base_colors = list(mcolors.TABLEAU_COLORS.values())

        while len(base_colors) < len(pathways):
            base_colors.extend(base_colors)

        colors = {}
        used_colors = set()

        for pathway in sorted(pathways):
            hash_val = int(hashlib.md5(pathway.encode()).hexdigest(), 16)
            color_idx = hash_val % len(base_colors)
            color = base_colors[color_idx]

            attempt = 0
            while color in used_colors and attempt < len(base_colors):
                color_idx = (color_idx + 1) % len(base_colors)
                color = base_colors[color_idx]
                attempt += 1

            if color in used_colors:
                import random
                random.seed(hash_val)
                color = f"#{random.randint(0, 0xFFFFFF):06x}"

            colors[pathway] = color
            used_colors.add(color)

        logging.info(f"Assigned {len(colors)} unique colors to pathways")

        return colors

    def _assign_colors_fallback(self, pathways: Optional[List[str]] = None) -> Dict[str, str]:
        """
        Fallback method for assigning colors without matplotlib.
        Uses basic CSS/HTML color palette.
        """
        if pathways is None:
            pathways = self.get_all_pathways(exclude_empty=True)

        base_colors = [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
            '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
            '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
            '#393b79', '#637939', '#8c6d31', '#843c39', '#7b4173',
            '#5254a3', '#8ca252', '#bd9e39', '#ad494a', '#a55194'
        ]

        colors = {}
        used_colors = set()

        for pathway in sorted(pathways):
            hash_val = int(hashlib.md5(pathway.encode()).hexdigest(), 16)
            color_idx = hash_val % len(base_colors)
            color = base_colors[color_idx]

            attempt = 0
            while color in used_colors and attempt < len(base_colors):
                color_idx = (color_idx + 1) % len(base_colors)
                color = base_colors[color_idx]
                attempt += 1

            if color in used_colors:
                import random
                random.seed(hash_val)
                color = f"#{random.randint(0, 0xFFFFFF):06x}"

            colors[pathway] = color
            used_colors.add(color)

        return colors

    def get_coverage_stats(self) -> Dict[str, any]:
        """
        Return database coverage statistics.

        Returns:
            Dictionary with statistics:
            - total_reactions: Total reactions
            - reactions_with_pathways: Reactions with pathways
            - total_pathways: Total pathways
            - pathways_with_reactions: Pathways with reactions
            - avg_pathways_per_reaction: Average pathways per reaction
            - avg_reactions_per_pathway: Average reactions per pathway
            - max_pathways_per_reaction: Maximum pathways for one reaction
            - max_reactions_per_pathway: Maximum reactions in one pathway
        """
        rxn_counts = self.matrix.sum(axis=1)
        pw_counts = self.matrix.sum(axis=0)

        return {
            "total_reactions": len(self.all_reactions),
            "reactions_with_pathways": (rxn_counts > 0).sum(),
            "total_pathways": len(self.all_pathways),
            "pathways_with_reactions": (pw_counts > 0).sum(),
            "avg_pathways_per_reaction": rxn_counts.mean(),
            "avg_reactions_per_pathway": pw_counts.mean(),
            "max_pathways_per_reaction": rxn_counts.max(),
            "max_reactions_per_pathway": pw_counts.max(),
        }

    def find_multi_pathway_reactions(self, min_pathways: int = 3) -> List[tuple]:
        """
        Find reactions belonging to multiple pathways.
        """
        results = []

        for rxn_id in self.all_reactions:
            pathways = self.get_pathways_for_reaction(rxn_id)
            if len(pathways) >= min_pathways:
                results.append((rxn_id, len(pathways), pathways))

        results.sort(key=lambda x: x[1], reverse=True)

        return results


def load_pathway_database(matrix_path: str = "data/models/subsystem_matrix.csv") -> PathwayDatabase:
    """
    Convenience function to load pathway database.
    """
    return PathwayDatabase(matrix_path)
