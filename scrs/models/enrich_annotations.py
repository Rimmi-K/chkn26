#!/usr/bin/env python3
"""
Обогащает аннотации COBRA модели из базы данных xref.

Для каждой реакции и метаболита:
1. Ищет MNX ID по внешнему идентификатору модели
2. Получает все cross-references для найденного MNX ID
3. Добавляет их как аннотации в модель
"""

import cobra
import sqlite3
import sys
import re
from pathlib import Path


def normalize_model_id(model_id):
    """
    Нормализует ID из модели для поиска в базе данных.
    Убирает префиксы типа M_, R_, суффиксы компартментов.

    Примеры:
      M_glc__D_c -> glc__D
      R_ATPS4r -> ATPS4r
      glc_D_c -> glc_D
    """
    # Убираем префиксы M_, R_
    clean_id = re.sub(r'^[MR]_', '', model_id)

    # Убираем суффикс компартмента (последний _X где X - одна буква)
    clean_id = re.sub(r'_[a-z]$', '', clean_id)

    return clean_id


def find_mnx_id(cursor, model_id, obj_type):
    """
    Ищет MNX ID для данного ID модели.

    Args:
        cursor: database cursor
        model_id: ID объекта в модели
        obj_type: 'reaction' или 'metabolite'

    Returns:
        MNX ID или None
    """
    # Пробуем разные варианты ID
    id_variants = [
        model_id,  # оригинальный ID
        normalize_model_id(model_id),  # нормализованный
        model_id.replace('__', '_'),  # с одинарным подчеркиванием
        model_id.upper(),  # в верхнем регистре
        model_id.lower(),  # в нижнем регистре
    ]

    for variant in id_variants:
        # Ищем в базе данных
        cursor.execute("""
            SELECT DISTINCT mnx_id
            FROM xref
            WHERE LOWER(external_id) = LOWER(?)
            AND type = ?
            LIMIT 1
        """, (variant, obj_type))

        result = cursor.fetchone()
        if result:
            return result[0]

    return None


def get_all_xrefs_for_mnx(cursor, mnx_id):
    """
    Получает все cross-references для данного MNX ID.

    Returns:
        список кортежей (source, external_id)
    """
    cursor.execute("""
        SELECT source, external_id
        FROM xref
        WHERE mnx_id = ?
    """, (mnx_id,))

    return cursor.fetchall()


def map_source_to_annotation_key(source, obj_type):
    """
    Преобразует название источника в ключ аннотации для COBRA модели.

    Args:
        source: источник из базы данных (e.g., 'seed.reaction', 'CHEBI')
        obj_type: 'reaction' или 'metabolite'

    Returns:
        ключ аннотации (e.g., 'seed.reaction', 'chebi')
    """
    source_lower = source.lower()

    # Специальные случаи
    mapping = {
        'chebi': 'chebi',
        'kegg': 'kegg.reaction' if obj_type == 'reaction' else 'kegg.compound',
        'rhea': 'rhea',
        'reactome': 'reactome',
        'biocyc': 'biocyc',
        'metanetx': 'metanetx.reaction' if obj_type == 'reaction' else 'metanetx.chemical',
    }

    # Проверяем точные совпадения
    if source_lower in mapping:
        return mapping[source_lower]

    # Для источников с уже правильным форматом (seed.reaction, bigg.metabolite)
    if '.' in source:
        return source

    # Для остальных - используем как есть
    return source


def enrich_object_annotation(obj, cursor, obj_type='reaction', verbose=False):
    """
    Обогащает аннотацию объекта (реакции или метаболита).

    Args:
        obj: объект COBRA (reaction или metabolite)
        cursor: database cursor
        obj_type: 'reaction' или 'metabolite'
        verbose: показывать детальное логирование

    Returns:
        (enriched, added_count): обогащён ли объект и кол-во добавленных аннотаций
    """
    # Сохраняем аннотации до обогащения для логирования
    if verbose:
        old_annotation = dict(obj.annotation)

    # Ищем MNX ID для данного объекта
    mnx_id = find_mnx_id(cursor, obj.id, obj_type)

    if not mnx_id:
        if verbose:
            print(f"    ❌ Не найдено в БД: {obj.id}")
        return False, 0

    # Получаем все cross-references для найденного MNX ID
    xrefs = get_all_xrefs_for_mnx(cursor, mnx_id)

    if not xrefs:
        if verbose:
            print(f"    ⚠️  MNX ID найден ({mnx_id}), но нет cross-references: {obj.id}")
        return False, 0

    # Добавляем аннотации
    added_count = 0
    for source, external_id in xrefs:
        # Преобразуем source в ключ аннотации
        annotation_key = map_source_to_annotation_key(source, obj_type)

        # Проверяем, существует ли уже этот ключ
        if annotation_key not in obj.annotation:
            obj.annotation[annotation_key] = []

        # Добавляем external_id, если его еще нет
        if isinstance(obj.annotation[annotation_key], list):
            if external_id not in obj.annotation[annotation_key]:
                obj.annotation[annotation_key].append(external_id)
                added_count += 1
        else:
            # Если это не список, преобразуем в список
            existing = obj.annotation[annotation_key]
            obj.annotation[annotation_key] = [existing]
            if external_id not in obj.annotation[annotation_key]:
                obj.annotation[annotation_key].append(external_id)
                added_count += 1

    if verbose and added_count > 0:
        print(f"    ✓ {obj.id} -> MNX: {mnx_id}: добавлено {added_count} аннотаций из {len(xrefs)} источников")
        # Показываем только новые аннотации
        new_keys = set(obj.annotation.keys()) - set(old_annotation.keys())
        if new_keys:
            print(f"      Новые источники: {', '.join(sorted(new_keys))}")

    return added_count > 0, added_count


def enrich_model(model_path, db_path, output_path=None, verbose_count=5):
    """
    Обогащает аннотации COBRA модели из базы данных xref.

    Args:
        model_path: путь к COBRA модели в формате SBML
        db_path: путь к базе данных xref.db
        output_path: путь для сохранения обогащённой модели
        verbose_count: количество объектов для детального логирования
    """
    print("=" * 70)
    print("Обогащение аннотаций COBRA модели")
    print("=" * 70)

    print(f"\nЗагружаем модель: {model_path}")
    model = cobra.io.read_sbml_model(model_path)

    print(f"  Реакций: {len(model.reactions)}")
    print(f"  Метаболитов: {len(model.metabolites)}")

    # Подключаемся к базе данных
    print(f"\nПодключаемся к базе данных: {db_path}")
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Проверяем структуру БД
    cursor.execute("SELECT COUNT(*) FROM xref")
    total_xrefs = cursor.fetchone()[0]
    print(f"  Записей в БД: {total_xrefs:,}")

    # Обогащаем реакции
    print("\n" + "=" * 70)
    print("Обогащаем аннотации реакций...")
    print("=" * 70)

    reactions_enriched = 0
    reactions_not_found = 0
    total_annotations_added = 0

    for i, reaction in enumerate(model.reactions):
        # Показываем детальное логирование для первых N обогащенных реакций
        verbose = (reactions_enriched < verbose_count)

        if (i + 1) % 100 == 0:
            print(f"  Обработано {i + 1}/{len(model.reactions)} реакций... "
                  f"(обогащено: {reactions_enriched}, не найдено: {reactions_not_found})")

        enriched, count = enrich_object_annotation(reaction, cursor, obj_type='reaction', verbose=verbose)
        if enriched:
            reactions_enriched += 1
            total_annotations_added += count
        else:
            reactions_not_found += 1

    print(f"\nРезультаты обогащения реакций:")
    print(f"  ✓ Обогащено: {reactions_enriched} ({reactions_enriched/len(model.reactions)*100:.1f}%)")
    print(f"  ✗ Не найдено: {reactions_not_found} ({reactions_not_found/len(model.reactions)*100:.1f}%)")
    print(f"  + Всего добавлено аннотаций: {total_annotations_added:,}")
    if reactions_enriched > 0:
        print(f"  ø Среднее аннотаций на реакцию: {total_annotations_added/reactions_enriched:.1f}")

    # Обогащаем метаболиты
    print("\n" + "=" * 70)
    print("Обогащаем аннотации метаболитов...")
    print("=" * 70)

    metabolites_enriched = 0
    metabolites_not_found = 0
    total_met_annotations_added = 0

    for i, metabolite in enumerate(model.metabolites):
        # Показываем детальное логирование для первых N обогащенных метаболитов
        verbose = (metabolites_enriched < verbose_count)

        if (i + 1) % 100 == 0:
            print(f"  Обработано {i + 1}/{len(model.metabolites)} метаболитов... "
                  f"(обогащено: {metabolites_enriched}, не найдено: {metabolites_not_found})")

        enriched, count = enrich_object_annotation(metabolite, cursor, obj_type='metabolite', verbose=verbose)
        if enriched:
            metabolites_enriched += 1
            total_met_annotations_added += count
        else:
            metabolites_not_found += 1

    print(f"\nРезультаты обогащения метаболитов:")
    print(f"  ✓ Обогащено: {metabolites_enriched} ({metabolites_enriched/len(model.metabolites)*100:.1f}%)")
    print(f"  ✗ Не найдено: {metabolites_not_found} ({metabolites_not_found/len(model.metabolites)*100:.1f}%)")
    print(f"  + Всего добавлено аннотаций: {total_met_annotations_added:,}")
    if metabolites_enriched > 0:
        print(f"  ø Среднее аннотаций на метаболит: {total_met_annotations_added/metabolites_enriched:.1f}")

    conn.close()

    # Сохраняем обогащенную модель
    if output_path is None:
        output_path = model_path.replace('.xml', '_enriched.xml')

    print("\n" + "=" * 70)
    print(f"Сохраняем обогащенную модель: {output_path}")
    cobra.io.write_sbml_model(model, output_path)

    print("\n" + "=" * 70)
    print("Готово! Статистика обогащения:")
    print("=" * 70)
    print(f"  Реакции: {reactions_enriched}/{len(model.reactions)} "
          f"({reactions_enriched/len(model.reactions)*100:.1f}%), "
          f"{total_annotations_added:,} аннотаций")
    print(f"  Метаболиты: {metabolites_enriched}/{len(model.metabolites)} "
          f"({metabolites_enriched/len(model.metabolites)*100:.1f}%), "
          f"{total_met_annotations_added:,} аннотаций")

    # Показываем примеры обогащенных аннотаций
    print("\n" + "=" * 70)
    print("Примеры обогащенных аннотаций:")
    print("=" * 70)

    # Находим реакции с аннотациями
    reactions_with_annotations = [r for r in model.reactions if r.annotation]
    if reactions_with_annotations:
        print(f"\nПример реакций (показано {min(3, len(reactions_with_annotations))} из {len(reactions_with_annotations)}):")
        for reaction in reactions_with_annotations[:3]:
            print(f"\n  {reaction.id}: {reaction.name}")
            for key, values in sorted(reaction.annotation.items()):
                if isinstance(values, list):
                    print(f"    {key}: {', '.join(map(str, values[:3]))}" +
                          (f" (+{len(values)-3} ещё)" if len(values) > 3 else ""))
                else:
                    print(f"    {key}: {values}")

    # Находим метаболиты с аннотациями
    metabolites_with_annotations = [m for m in model.metabolites if m.annotation]
    if metabolites_with_annotations:
        print(f"\nПример метаболитов (показано {min(3, len(metabolites_with_annotations))} из {len(metabolites_with_annotations)}):")
        for metabolite in metabolites_with_annotations[:3]:
            print(f"\n  {metabolite.id}: {metabolite.name}")
            for key, values in sorted(metabolite.annotation.items()):
                if isinstance(values, list):
                    print(f"    {key}: {', '.join(map(str, values[:3]))}" +
                          (f" (+{len(values)-3} ещё)" if len(values) > 3 else ""))
                else:
                    print(f"    {key}: {values}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Обогащает аннотации COBRA модели из базы данных xref',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Примеры использования:

  # Базовое использование (используются пути по умолчанию)
  python enrich_annotations.py

  # Указание пути к модели
  python enrich_annotations.py --model path/to/model.xml

  # Указание всех параметров
  python enrich_annotations.py --model model.xml --db xref.db --output enriched.xml

  # Показать больше примеров в логе
  python enrich_annotations.py --verbose 10
        """
    )
    parser.add_argument(
        '--model',
        default='iES1300_metpath2.xml',
        help='Путь к COBRA модели (по умолчанию: iES1300_metpath2.xml)'
    )
    parser.add_argument(
        '--db',
        default='xref.db',
        help='Путь к базе данных xref (по умолчанию: xref.db)'
    )
    parser.add_argument(
        '--output',
        help='Путь к выходному файлу (по умолчанию: <model>_enriched.xml)'
    )
    parser.add_argument(
        '--verbose',
        type=int,
        default=5,
        help='Количество объектов для детального логирования (по умолчанию: 5)'
    )

    args = parser.parse_args()

    # Проверяем существование файлов
    if not Path(args.model).exists():
        print(f"Ошибка: файл модели не найден: {args.model}")
        sys.exit(1)

    if not Path(args.db).exists():
        print(f"Ошибка: база данных не найдена: {args.db}")
        print(f"\nСначала запустите скрипт tsv_to_db.py для создания базы данных:")
        print(f"  python tsv_to_db.py")
        sys.exit(1)

    try:
        enrich_model(args.model, args.db, args.output, args.verbose)
    except Exception as e:
        print(f"\nОшибка при обогащении модели: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
