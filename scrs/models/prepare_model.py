#!/usr/bin/env python3
import os
import pickle
import logging
import time

import cobra
import pandas as pd
from bioservices import KEGG

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

kegg = KEGG()


def _to_title_case(text: str) -> str | None:
    if not text or pd.isna(text):
        return None
    return ' '.join(word.capitalize() for word in text.split())


def get_kegg_ids(reaction) -> list:
    ids = []
    for key in ('kegg.reaction', 'keggR', 'kegg'):
        val = reaction.annotation.get(key)
        if val is None:
            continue
        ids.extend(val if isinstance(val, list) else [val])
    return sorted(set(ids))


def subsystem_to_list(subsystem) -> list:
    if not subsystem:
        return []
    if isinstance(subsystem, list):
        return [s.strip() for s in subsystem if s and s.strip()]
    for sep in ('; ', ', '):
        if sep in subsystem:
            return [s.strip() for s in subsystem.split(sep[0]) if s.strip()]
    return [subsystem.strip()] if subsystem.strip() else []


def list_to_subsystem(pathways: list) -> str | None:
    unique = sorted(set(_to_title_case(p) for p in pathways if p))
    return '; '.join(unique) if unique else None


def get_pathways_from_kegg(kegg_id: str, cache_file: str = 'pathways_cache.pkl') -> list | None:
    cache = {}
    if os.path.exists(cache_file):
        with open(cache_file, 'rb') as f:
            cache = pickle.load(f)
        if kegg_id in cache:
            return cache[kegg_id]

    pathways = None
    try:
        info = kegg.get(kegg_id)
        if info and 'PATHWAY' in info:
            pathways = []
            in_pathway = False
            for line in info.split('\n'):
                if line.startswith('PATHWAY'):
                    in_pathway = True
                    parts = line.replace('PATHWAY', '').strip().split(maxsplit=1)
                    if len(parts) == 2:
                        pathways.append(parts[1])
                elif in_pathway:
                    if line.startswith(' ' * 12):
                        parts = line.strip().split(maxsplit=1)
                        if len(parts) == 2:
                            pathways.append(parts[1])
                    else:
                        break
        time.sleep(0.5)
    except Exception as e:
        logging.warning(f"KEGG error for {kegg_id}: {e}")

    cache[kegg_id] = pathways
    with open(cache_file, 'wb') as f:
        pickle.dump(cache, f)
    return pathways


def annotate_pathways_from_kegg(model, output_model_path: str = "prepared_model.xml"):
    records = []
    updated = 0

    for rxn in model.reactions:
        if rxn.id.startswith(('EX_', 'DM_', 'SK_')):
            continue

        kegg_ids = get_kegg_ids(rxn)
        old_pathways = subsystem_to_list(rxn.subsystem)

        new_from_kegg = []
        if kegg_ids:
            for kid in kegg_ids:
                try:
                    paths = get_pathways_from_kegg(kid)
                    if paths:
                        new_from_kegg.extend(paths)
                except Exception as e:
                    logging.error(f"{rxn.id} (KEGG: {kid}): {e}")

        all_pathways = set(old_pathways) | set(new_from_kegg)
        final = list_to_subsystem(list(all_pathways))

        if final != rxn.subsystem:
            rxn.subsystem = final
            updated += 1
            logging.info(f"{rxn.id}: updated ({len(all_pathways)} pathways)")
            status = 'updated'
        else:
            status = 'unchanged'

        records.append({
            'Reactions':           rxn.id,
            'Old_Pathways':        '; '.join(old_pathways) if old_pathways else None,
            'KEGG_IDs':            ', '.join(kegg_ids) if kegg_ids else None,
            'New_Pathways_From_KEGG': '; '.join(sorted(set(new_from_kegg))) if new_from_kegg else None,
            'Final_Pathways':      final,
            'Status':              status,
        })

    cobra.io.write_sbml_model(model, output_model_path)
    return model, pd.DataFrame(records)


def prepare_model(model_path: str, output_model_path: str, output_pathways_path: str):
    model = cobra.io.read_sbml_model(model_path)
    model, pathways_df = annotate_pathways_from_kegg(model, output_model_path)
    pathways_df.to_csv(output_pathways_path, index=False)
    return model, pathways_df


if __name__ == "__main__":
    model, pathways_df = prepare_model(
        model_path='data/models/genenametransport2.xml',
        output_model_path='data/models/iES1300_prepared3.xml',
        output_pathways_path='pathways_annotation3.csv',
    )
    updated = pathways_df[pathways_df['Status'] == 'updated']
    print(f"Updated reactions: {len(updated)}/{len(pathways_df)}")
