#!/usr/bin/env python3
import sqlite3

SOURCE_NORMALIZATION = {
    'seedR':   'seed.reaction',
    'seedM':   'seed.compound',
    'biggR':   'bigg.reaction',
    'biggM':   'bigg.metabolite',
    'vmhR':    'vmh.reaction',
    'vmhM':    'vmh.metabolite',
    'rh':      'rhea',
    'sabiorkR':'sabiork.reaction',
    'slm':     'slm',
    'SLM':     'slm',
}


def _parse_source_id(source_id_str: str) -> tuple[str | None, str | None]:
    if not source_id_str or ':' not in source_id_str:
        return None, None
    source, external_id = source_id_str.split(':', 1)
    return source.strip(), external_id.strip()


def read_xref_file(filename: str, record_type: str) -> list:
    records = []
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            source, external_id = _parse_source_id(parts[0])
            mnx_id = parts[1]
            description = parts[2] if len(parts) > 2 else ''

            if not source or not external_id or not mnx_id:
                continue

            records.append({
                'mnx_id':      mnx_id,
                'source':      SOURCE_NORMALIZATION.get(source, source),
                'external_id': external_id,
                'description': description,
                'type':        record_type,
            })
    return records


def create_database(db_path: str, records: list):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS xref (
            mnx_id      TEXT NOT NULL,
            source      TEXT NOT NULL,
            external_id TEXT NOT NULL,
            description TEXT,
            type        TEXT
        )
    ''')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_mnx_id ON xref(mnx_id)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_external_id ON xref(external_id COLLATE NOCASE)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_source ON xref(source)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_source_external ON xref(source, external_id COLLATE NOCASE)')

    batch_size = 10000
    for i in range(0, len(records), batch_size):
        cursor.executemany(
            'INSERT INTO xref (mnx_id, source, external_id, description, type) '
            'VALUES (:mnx_id, :source, :external_id, :description, :type)',
            records[i:i + batch_size],
        )

    conn.commit()
    conn.close()


def _print_db_stats(db_path: str):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute('SELECT COUNT(*) FROM xref')
    print(f"Total records: {cursor.fetchone()[0]:,}")

    cursor.execute('SELECT COUNT(DISTINCT mnx_id) FROM xref WHERE type="reaction"')
    print(f"Unique reactions: {cursor.fetchone()[0]:,}")

    cursor.execute('SELECT COUNT(DISTINCT mnx_id) FROM xref WHERE type="metabolite"')
    print(f"Unique metabolites: {cursor.fetchone()[0]:,}")

    cursor.execute('SELECT DISTINCT source FROM xref ORDER BY source')
    for (source,) in cursor.fetchall()[:20]:
        cursor.execute('SELECT COUNT(*) FROM xref WHERE source=?', (source,))
        print(f"  {source:30s}: {cursor.fetchone()[0]:,}")

    for record_type in ('reaction', 'metabolite'):
        print(f"\nSample {record_type}s:")
        cursor.execute(f'SELECT mnx_id, source, external_id, description FROM xref WHERE type=? LIMIT 3', (record_type,))
        for mnx_id, source, external_id, desc in cursor.fetchall():
            print(f"  MNX:{mnx_id}  {source}:{external_id}")
            print(f"    {desc[:50]}{'...' if len(desc) > 50 else ''}")

    conn.close()


def main():
    records = read_xref_file('reac_xref.tsv', 'reaction') + read_xref_file('chem_xref.tsv', 'metabolite')
    db_path = 'xref.db'
    create_database(db_path, records)
    print(f"Database created: {db_path}")
    _print_db_stats(db_path)


if __name__ == '__main__':
    main()
