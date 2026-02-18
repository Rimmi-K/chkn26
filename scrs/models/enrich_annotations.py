#!/usr/bin/env python3
import re
import sys
import sqlite3
from pathlib import Path

import cobra


def normalize_model_id(model_id: str) -> str:
    clean = re.sub(r'^[MR]_', '', model_id)
    clean = re.sub(r'_[a-z]$', '', clean)
    return clean


def find_mnx_id(cursor, model_id: str, obj_type: str):
    variants = [
        model_id,
        normalize_model_id(model_id),
        model_id.replace('__', '_'),
        model_id.upper(),
        model_id.lower(),
    ]
    for variant in variants:
        cursor.execute(
            "SELECT DISTINCT mnx_id FROM xref WHERE LOWER(external_id) = LOWER(?) AND type = ? LIMIT 1",
            (variant, obj_type),
        )
        row = cursor.fetchone()
        if row:
            return row[0]
    return None


def get_all_xrefs(cursor, mnx_id: str) -> list:
    cursor.execute("SELECT source, external_id FROM xref WHERE mnx_id = ?", (mnx_id,))
    return cursor.fetchall()


def _annotation_key(source: str, obj_type: str) -> str:
    mapping = {
        'chebi':    'chebi',
        'kegg':     'kegg.reaction' if obj_type == 'reaction' else 'kegg.compound',
        'rhea':     'rhea',
        'reactome': 'reactome',
        'biocyc':   'biocyc',
        'metanetx': 'metanetx.reaction' if obj_type == 'reaction' else 'metanetx.chemical',
    }
    return mapping.get(source.lower(), source)


def enrich_object(obj, cursor, obj_type: str = 'reaction', verbose: bool = False) -> tuple[bool, int]:
    mnx_id = find_mnx_id(cursor, obj.id, obj_type)
    if not mnx_id:
        if verbose:
            print(f"Not found: {obj.id}")
        return False, 0

    xrefs = get_all_xrefs(cursor, mnx_id)
    if not xrefs:
        if verbose:
            print(f"No xrefs for {obj.id} (MNX: {mnx_id})")
        return False, 0

    added = 0
    for source, external_id in xrefs:
        key = _annotation_key(source, obj_type)
        current = obj.annotation.setdefault(key, [])
        if not isinstance(current, list):
            obj.annotation[key] = [current]
            current = obj.annotation[key]
        if external_id not in current:
            current.append(external_id)
            added += 1

    if verbose and added > 0:
        print(f"  {obj.id} -> MNX:{mnx_id}, +{added} annotations")
    return added > 0, added


def enrich_model(model_path: str, db_path: str, output_path: str = None, verbose_count: int = 5):
    model = cobra.io.read_sbml_model(model_path)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    rxn_enriched = rxn_missing = rxn_added = 0
    for i, rxn in enumerate(model.reactions):
        verbose = rxn_enriched < verbose_count
        if (i + 1) % 100 == 0:
            print(f"  Reactions {i+1}/{len(model.reactions)} "
                  f"(enriched: {rxn_enriched}, missing: {rxn_missing})")
        ok, count = enrich_object(rxn, cursor, 'reaction', verbose)
        if ok:
            rxn_enriched += 1
            rxn_added += count
        else:
            rxn_missing += 1

    met_enriched = met_missing = met_added = 0
    for met in model.metabolites:
        verbose = met_enriched < verbose_count
        ok, count = enrich_object(met, cursor, 'metabolite', verbose)
        if ok:
            met_enriched += 1
            met_added += count
        else:
            met_missing += 1

    conn.close()

    if output_path is None:
        output_path = model_path.replace('.xml', '_enriched.xml')
    cobra.io.write_sbml_model(model, output_path)

    print(f"\nReactions: {rxn_enriched} enriched, {rxn_missing} not found, {rxn_added} annotations added")
    print(f"Metabolites: {met_enriched} enriched, {met_missing} not found, {met_added} annotations added")
    print(f"Saved: {output_path}")

    for rxn in list(filter(lambda r: r.annotation, model.reactions))[:3]:
        print(f"\n  {rxn.id}: {rxn.name}")
        for key, vals in sorted(rxn.annotation.items()):
            vals_list = vals if isinstance(vals, list) else [vals]
            preview = ', '.join(map(str, vals_list[:3]))
            suffix = f" (+{len(vals_list)-3} more)" if len(vals_list) > 3 else ""
            print(f"    {key}: {preview}{suffix}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Enrich COBRA model annotations from xref database')
    parser.add_argument('--model',   default='iES1300_metpath2.xml')
    parser.add_argument('--db',      default='xref.db')
    parser.add_argument('--output',  default=None)
    parser.add_argument('--verbose', type=int, default=5)
    args = parser.parse_args()

    if not Path(args.db).exists():
        print(f"Database not found: {args.db}. Run tsv_to_db.py first.")
        sys.exit(1)

    try:
        enrich_model(args.model, args.db, args.output, args.verbose)
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
