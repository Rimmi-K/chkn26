#!/usr/bin/env python3
"""
Конвертирует TSV файлы cross-references (reac_xref.tsv и chem_xref.tsv) в базу данных SQLite.

Структура xref файлов:
- source:external_id  MNX_ID  description

Создает таблицу xref с колонками:
- mnx_id: MetaNetX ID (MNXR... для реакций, MNXM... для метаболитов)
- source: источник данных (seed.reaction, CHEBI, kegg, etc.)
- external_id: внешний идентификатор
- description: описание и альтернативные имена
- type: тип записи (reaction или metabolite)
"""

import sqlite3
import sys
import re


def parse_source_id(source_id_str):
    """
    Разделяет source:id на source и external_id.
    Примеры:
      seed.reaction:rxn31178 -> ('seed.reaction', 'rxn31178')
      seedR:rxn08730 -> ('seedR', 'rxn08730')
      CHEBI:165749 -> ('CHEBI', '165749')
    """
    if not source_id_str or ':' not in source_id_str:
        return None, None

    parts = source_id_str.split(':', 1)
    return parts[0].strip(), parts[1].strip()


def normalize_source(source):
    """
    Нормализует название источника к стандартному формату.
    """
    # Словарь для нормализации
    normalization = {
        'seedR': 'seed.reaction',
        'seedM': 'seed.compound',
        'biggR': 'bigg.reaction',
        'biggM': 'bigg.metabolite',
        'vmhR': 'vmh.reaction',
        'vmhM': 'vmh.metabolite',
        'rh': 'rhea',
        'sabiorkR': 'sabiork.reaction',
        'slm': 'slm',
        'SLM': 'slm',
    }

    return normalization.get(source, source)


def read_xref_file(filename, record_type):
    """
    Читает xref TSV файл, пропуская строки преамбулы (начинающиеся с #).

    Args:
        filename: путь к файлу
        record_type: 'reaction' или 'metabolite'

    Returns:
        список словарей с данными
    """
    records = []

    with open(filename, 'r', encoding='utf-8') as f:
        # Пропускаем все комментарии
        for line in f:
            if line.startswith('#'):
                continue

            # Обрабатываем строку с данными
            parts = line.strip().split('\t')

            if len(parts) < 2:
                continue

            source_id_str = parts[0]
            mnx_id = parts[1]
            description = parts[2] if len(parts) > 2 else ''

            # Парсим source:external_id
            source, external_id = parse_source_id(source_id_str)

            if not source or not external_id or not mnx_id:
                continue

            # Нормализуем источник
            source = normalize_source(source)

            records.append({
                'mnx_id': mnx_id,
                'source': source,
                'external_id': external_id,
                'description': description,
                'type': record_type
            })

    return records


def create_database(db_path, records):
    """
    Создает SQLite базу данных и заполняет её данными.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Создаем таблицу
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS xref (
            mnx_id TEXT NOT NULL,
            source TEXT NOT NULL,
            external_id TEXT NOT NULL,
            description TEXT,
            type TEXT
        )
    ''')

    # Создаем индексы для быстрого поиска
    cursor.execute('''
        CREATE INDEX IF NOT EXISTS idx_mnx_id
        ON xref(mnx_id)
    ''')

    cursor.execute('''
        CREATE INDEX IF NOT EXISTS idx_external_id
        ON xref(external_id COLLATE NOCASE)
    ''')

    cursor.execute('''
        CREATE INDEX IF NOT EXISTS idx_source
        ON xref(source)
    ''')

    cursor.execute('''
        CREATE INDEX IF NOT EXISTS idx_source_external
        ON xref(source, external_id COLLATE NOCASE)
    ''')

    # Вставляем данные пакетами для ускорения
    batch_size = 10000
    for i in range(0, len(records), batch_size):
        batch = records[i:i+batch_size]
        cursor.executemany('''
            INSERT INTO xref (mnx_id, source, external_id, description, type)
            VALUES (:mnx_id, :source, :external_id, :description, :type)
        ''', batch)

        if (i + batch_size) % 100000 == 0:
            print(f"  Вставлено {i + batch_size} записей...")

    conn.commit()
    conn.close()


def main():
    print("=" * 60)
    print("Конвертация xref.tsv файлов в базу данных SQLite")
    print("=" * 60)

    print("\nЧитаем reac_xref.tsv...")
    reac_records = read_xref_file('reac_xref.tsv', 'reaction')
    print(f"  Найдено {len(reac_records):,} записей реакций")

    print("\nЧитаем chem_xref.tsv...")
    chem_records = read_xref_file('chem_xref.tsv', 'metabolite')
    print(f"  Найдено {len(chem_records):,} записей метаболитов")

    # Объединяем все записи
    all_records = reac_records + chem_records
    print(f"\nВсего записей: {len(all_records):,}")

    # Создаем базу данных
    db_path = 'xref.db'
    print(f"\nСоздаем базу данных: {db_path}...")
    create_database(db_path, all_records)

    print("\nГотово!")

    # Проверяем созданную базу данных
    print("\n" + "=" * 60)
    print("Статистика базы данных")
    print("=" * 60)

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute('SELECT COUNT(*) FROM xref')
    count = cursor.fetchone()[0]
    print(f"\nВсего записей в БД: {count:,}")

    cursor.execute('SELECT COUNT(DISTINCT mnx_id) FROM xref WHERE type="reaction"')
    unique_reactions = cursor.fetchone()[0]
    print(f"Уникальных реакций (MNX ID): {unique_reactions:,}")

    cursor.execute('SELECT COUNT(DISTINCT mnx_id) FROM xref WHERE type="metabolite"')
    unique_metabolites = cursor.fetchone()[0]
    print(f"Уникальных метаболитов (MNX ID): {unique_metabolites:,}")

    cursor.execute('SELECT DISTINCT source FROM xref ORDER BY source')
    sources = [row[0] for row in cursor.fetchall()]
    print(f"\nИсточники данных ({len(sources)}):")
    for source in sources[:20]:  # показываем первые 20
        cursor.execute('SELECT COUNT(*) FROM xref WHERE source=?', (source,))
        count = cursor.fetchone()[0]
        print(f"  {source:30s} : {count:,}")

    if len(sources) > 20:
        print(f"  ... и еще {len(sources) - 20} источников")

    # Показываем несколько примеров
    print("\n" + "=" * 60)
    print("Примеры записей")
    print("=" * 60)

    print("\nПример реакций:")
    cursor.execute('SELECT * FROM xref WHERE type="reaction" LIMIT 3')
    for row in cursor.fetchall():
        mnx_id, source, external_id, desc, rtype = row
        desc_short = desc[:50] + '...' if len(desc) > 50 else desc
        print(f"  MNX: {mnx_id}, {source}:{external_id}")
        print(f"    {desc_short}")

    print("\nПример метаболитов:")
    cursor.execute('SELECT * FROM xref WHERE type="metabolite" LIMIT 3')
    for row in cursor.fetchall():
        mnx_id, source, external_id, desc, rtype = row
        desc_short = desc[:50] + '...' if len(desc) > 50 else desc
        print(f"  MNX: {mnx_id}, {source}:{external_id}")
        print(f"    {desc_short}")

    conn.close()


if __name__ == '__main__':
    main()
