import pandas as pd
import cobra


def create_subsystem_dataframe(model) -> pd.DataFrame:
    all_subsystems = set()
    for rxn in model.reactions:
        if rxn.subsystem:
            all_subsystems.update(s.strip() for s in rxn.subsystem.split(';'))

    subsystem_columns = sorted(all_subsystems)
    reaction_ids = [rxn.id for rxn in model.reactions]

    data = {
        sub: [
            1 if rxn.subsystem and sub in {s.strip() for s in rxn.subsystem.split(';')} else 0
            for rxn in model.reactions
        ]
        for sub in subsystem_columns
    }

    return pd.DataFrame(data, index=reaction_ids)


if __name__ == "__main__":
    model = cobra.io.load_json_model("data/models/model_curated.json")
    df = create_subsystem_dataframe(model)
    df.to_csv('data/models/subsystem_matrix.csv', index=True)
