#!/usr/bin/env python3
from collections import Counter
from pathlib import Path

import cobra
import pandas as pd

MODEL_PATH = Path("data/models/curated_model51.xml")
SUBSYSTEM_MATRIX_PATH = Path("data/models/subsystem_matrix.csv")
OUTPUT_PATH = Path("data/models/subsystem_matrix_updated.csv")
NON_TRANSPORT_OUTPUT = Path("data/models/non_transport_reactions.txt")

TRANSPORT_CATEGORIES = {
    "Transport, Mitochondrial":       ["c", "m"],
    "Transport, Peroxisomal":         ["c", "x"],
    "Transport, Endoplasmic Reticulum": ["c", "r"],
    "Transport, Golgi":               ["c", "g"],
    "Transport, Lysosomal":           ["c", "l"],
    "Transport, Nuclear":             ["c", "n"],
    "Transport, Endosomal":           ["c", "e"],
    "Extracellular Transport":        ["c", "s"],
}


def get_compartments(reaction) -> list:
    return sorted({met.compartment for met in reaction.metabolites})


def classify_transport(reaction) -> str | None:
    compartments = set(get_compartments(reaction))
    if len(compartments) < 2:
        return None
    for category, pattern in TRANSPORT_CATEGORIES.items():
        if set(pattern) == compartments:
            return category
    for category, pattern in TRANSPORT_CATEGORIES.items():
        if set(pattern).issubset(compartments):
            return category
    return "Transport, Other"


def find_new_reactions(model, matrix: pd.DataFrame) -> set:
    return {rxn.id for rxn in model.reactions} - set(matrix.index)


def classify_new_reactions(model, new_ids: set) -> tuple[dict, list]:
    transport, non_transport = {}, []
    for rxn_id in new_ids:
        rxn = model.reactions.get_by_id(rxn_id)
        category = classify_transport(rxn)
        if category:
            transport[rxn_id] = category
        else:
            non_transport.append(rxn_id)
    return transport, non_transport


def update_subsystem_matrix(matrix: pd.DataFrame, transport_reactions: dict) -> pd.DataFrame:
    updated = matrix.copy()
    for category in set(transport_reactions.values()):
        if category not in updated.columns:
            updated[category] = 0
    for rxn_id, category in transport_reactions.items():
        row = pd.Series(0, index=updated.columns, name=rxn_id)
        row[category] = 1
        updated = pd.concat([updated, row.to_frame().T])
    return updated


def _show_transport_examples(model, transport_reactions: dict, n: int = 5):
    for rxn_id, category in list(transport_reactions.items())[:n]:
        rxn = model.reactions.get_by_id(rxn_id)
        print(f"\n  {rxn_id} [{category}]")
        for comp in get_compartments(rxn):
            mets = [f"{m.id} ({coef:+.1f})"
                    for m, coef in rxn.metabolites.items() if m.compartment == comp]
            print(f"    [{comp}]: {', '.join(mets)}")


def _save_non_transport(model, non_transport: list):
    with open(NON_TRANSPORT_OUTPUT, "w") as f:
        f.write("# reaction_id | name | reaction | compartments\n\n")
        for rxn_id in sorted(non_transport):
            rxn = model.reactions.get_by_id(rxn_id)
            comps = get_compartments(rxn)
            f.write(f"{rxn_id} | {rxn.name} | {rxn.reaction} | {comps}\n")
    print(f"Non-transport reactions saved: {NON_TRANSPORT_OUTPUT}")


def main():
    model = cobra.io.read_sbml_model(str(MODEL_PATH))
    matrix = pd.read_csv(SUBSYSTEM_MATRIX_PATH, index_col=0)

    new_ids = find_new_reactions(model, matrix)
    if not new_ids:
        print("No new reactions found.")
        return

    print(f"New reactions: {len(new_ids)}")
    transport, non_transport = classify_new_reactions(model, new_ids)

    if transport:
        counts = Counter(transport.values())
        print("\nTransport type distribution:")
        for t, c in sorted(counts.items()):
            print(f"  {t}: {c}")
        _show_transport_examples(model, transport, n=10)
        updated = update_subsystem_matrix(matrix, transport)
        updated.to_csv(OUTPUT_PATH)
        print(f"Updated matrix saved: {OUTPUT_PATH}")

    if non_transport:
        print(f"\nNon-transport reactions: {len(non_transport)}")
        for rxn_id in sorted(non_transport)[:20]:
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"  {rxn_id}: {rxn.name} [{', '.join(get_compartments(rxn))}]")
        _save_non_transport(model, non_transport)


if __name__ == "__main__":
    main()
