#!/usr/bin/env python3
"""
Исправленная версия prepare_model.py
Улучшения:
- Извлекает KEGG ID из аннотаций модели (kegg.reaction)
- Преобразует subsystem в список
- Добавляет новые пути к существующим
"""

import cobra
import pandas as pd
import requests
from bs4 import BeautifulSoup
import os
import logging
import pickle
from glob import glob
import time

# ===== НАСТРОЙКА ПРОКСИ ДЛЯ VPN =====
# Прокси отключен - KEGG доступен напрямую через HTTPS
# os.environ['HTTP_PROXY'] = 'socks5h://127.0.0.1:2080'
# os.environ['HTTPS_PROXY'] = 'socks5h://127.0.0.1:2080'
# os.environ['ALL_PROXY'] = 'socks5h://127.0.0.1:2080'

# Настройка прокси для requests
# PROXIES = {
#     'http': 'socks5h://127.0.0.1:2080',
#     'https': 'socks5h://127.0.0.1:2080',
# }
# ====================================

# Настройка логирования с явной кодировкой UTF-8
logging.basicConfig(
    level=logging.INFO,
    filename='prepare_model.log',
    filemode='w',
    format='%(asctime)s - %(levelname)s - %(message)s',
    encoding='utf-8'
)

# Настройка вывода в консоль с UTF-8
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

from bioservices import KEGG

# Инициализация KEGG API (использует HTTPS напрямую)
kegg = KEGG()
# Обновляем URL на HTTPS (bioservices использует старый HTTP URL)
if hasattr(kegg, 'services') and hasattr(kegg.services, 'url'):
    old_url = kegg.services.url
    kegg.services.url = old_url.replace('http://', 'https://')
    logging.info(f"KEGG API URL обновлен: {old_url} -> {kegg.services.url}")
else:
    logging.warning("Не удалось обновить KEGG URL")


def to_title_case(text):
    """Приводит строку к формату Title Case (каждое слово с большой буквы, остальные строчные)."""
    if not text or pd.isna(text):
        return None
    return ' '.join(word.capitalize() for word in text.split())


def get_kegg_ids_from_annotation(reaction):
    """
    Извлекает KEGG reaction IDs из аннотаций реакции.

    Returns:
        list: список KEGG IDs или пустой список
    """
    kegg_ids = []

    # Проверяем различные варианты ключей для KEGG
    kegg_keys = ['kegg.reaction', 'keggR', 'kegg']

    for key in kegg_keys:
        if key in reaction.annotation:
            value = reaction.annotation[key]

            # Обрабатываем как список, так и одиночное значение
            if isinstance(value, list):
                kegg_ids.extend(value)
            elif isinstance(value, str):
                kegg_ids.append(value)

    # Убираем дубликаты и сортируем
    return sorted(set(kegg_ids))


def subsystem_to_list(subsystem):
    """
    Преобразует subsystem в список.

    Args:
        subsystem: строка или список

    Returns:
        list: список путей или пустой список
    """
    if not subsystem:
        return []

    # Если уже список
    if isinstance(subsystem, list):
        return [s.strip() for s in subsystem if s and s.strip()]

    # Если строка
    if isinstance(subsystem, str):
        # Разделяем по "; " или ", "
        if '; ' in subsystem:
            return [s.strip() for s in subsystem.split(';') if s and s.strip()]
        elif ', ' in subsystem:
            return [s.strip() for s in subsystem.split(',') if s and s.strip()]
        else:
            return [subsystem.strip()] if subsystem.strip() else []

    return []


def list_to_subsystem(pathways_list):
    """
    Преобразует список путей в строку для subsystem.

    Args:
        pathways_list: список путей

    Returns:
        str: строка с путями через "; " или None
    """
    if not pathways_list:
        return None

    # Убираем дубликаты, сортируем и преобразуем в Title Case
    unique_pathways = sorted(set(to_title_case(p) for p in pathways_list if p))

    return '; '.join(unique_pathways) if unique_pathways else None


def get_pathways_from_kegg(kegg_id, cache_file='pathways_cache.pkl'):
    """
    Получает метаболические пути из KEGG по ID реакции.

    Args:
        kegg_id: KEGG reaction ID (например, R00303)
        cache_file: путь к файлу кэша

    Returns:
        list: список названий путей или None
    """
    # Кэш
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
            lines = info.split('\n')
            pathways = []
            in_pathway = False

            for line in lines:
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

            time.sleep(0.5)  # Rate limiting

            if pathways:
                cache[kegg_id] = pathways
                with open(cache_file, 'wb') as f:
                    pickle.dump(cache, f)
                return pathways
    except Exception as e:
        logging.warning(f"KEGG ошибка для {kegg_id}: {e}")

    # Сохраняем в кэш даже если ничего не найдено
    cache[kegg_id] = pathways
    with open(cache_file, 'wb') as f:
        pickle.dump(cache, f)

    return pathways


def load_model(model_path):
    """Загружает метаболическую модель из SBML файла."""
    try:
        model = cobra.io.read_sbml_model(model_path)
        logging.info(f"Модель загружена из {model_path}")
        return model
    except Exception as e:
        logging.error(f"Ошибка загрузки модели: {e}")
        raise


def annotate_pathways_from_kegg(model, output_model_path="prepared_model.xml"):
    """
    Аннотирует реакции модели метаболическими путями из KEGG.

    Args:
        model: COBRA модель
        output_model_path: путь для сохранения модели

    Returns:
        DataFrame с результатами
    """
    data = []
    reactions_processed = 0
    reactions_updated = 0
    reactions_kegg_found = 0

    # Работаем НАПРЯМУЮ с model.reactions, без промежуточных списков
    for reaction in model.reactions:
        # Пропускаем префиксы
        if reaction.id.startswith(('EX_', 'DM_', 'SK_')):
            continue

        reactions_processed += 1

        # Получаем KEGG IDs из аннотаций
        kegg_ids = get_kegg_ids_from_annotation(reaction)

        # Получаем существующие пути
        old_pathways = subsystem_to_list(reaction.subsystem)
        old_pathways_str = '; '.join(old_pathways) if old_pathways else None

        kegg_ids_str = ', '.join(kegg_ids) if kegg_ids else None

        # Собираем новые пути из KEGG
        new_pathways_from_kegg = []
        if kegg_ids:
            reactions_kegg_found += 1
            for kegg_id in kegg_ids:
                try:
                    pathways = get_pathways_from_kegg(kegg_id)
                    if pathways:
                        new_pathways_from_kegg.extend(pathways)
                        logging.info(f"{reaction.id}: найдено {len(pathways)} путей из KEGG ID {kegg_id}")
                except Exception as e:
                    logging.error(f"Ошибка для {reaction.id} (KEGG ID: {kegg_id}): {e}")

        new_pathways_str = '; '.join(sorted(set(new_pathways_from_kegg))) if new_pathways_from_kegg else None

        # Объединяем старые и новые пути
        all_pathways = set(old_pathways)

        if new_pathways_from_kegg:
            all_pathways.update(new_pathways_from_kegg)

        # Преобразуем в строку для subsystem
        final_pathways_str = list_to_subsystem(list(all_pathways))

        # Обновляем subsystem напрямую через get_by_id
        if final_pathways_str != reaction.subsystem:
            model.reactions.get_by_id(reaction.id).subsystem = final_pathways_str
            reactions_updated += 1
            status = 'updated'
            logging.info(f"{reaction.id}: subsystem обновлен ({len(all_pathways)} путей)")
        else:
            status = 'unchanged'

        data.append({
            'Reactions': reaction.id,
            'Old_Pathways': old_pathways_str,
            'KEGG_IDs': kegg_ids_str,
            'New_Pathways_From_KEGG': new_pathways_str,
            'Final_Pathways': final_pathways_str,
            'Status': status
        })
    save_model(model, output_model_path)

    df_pathways = pd.DataFrame(data)

    logging.info("=" * 70)
    logging.info("СТАТИСТИКА АННОТАЦИИ")
    logging.info("=" * 70)
    logging.info(f"Обработано реакций: {reactions_processed}")
    logging.info(f"Реакций с KEGG ID: {reactions_kegg_found}")
    logging.info(f"Реакций обновлено: {reactions_updated}")
    logging.info(f"Реакций без изменений: {reactions_processed - reactions_updated}")

    return model, df_pathways


def save_model(model, output_path):
    """
    Сохраняет модель в JSON формате (subsystems сохраняются).

    SBML не сохраняет subsystems корректно, поэтому используем JSON.
    """
    try:
        cobra.io.write_sbml_model(model, output_path)
    except Exception as e:
        logging.error(f"Ошибка сохранения модели: {e}")
        raise


def prepare_model(model_path, output_model_path, output_pathways_path):
    """
    Основная функция подготовки модели.

    Args:
        model_path: путь к входной модели
        output_model_path: путь для сохранения модели (будет .json)
        output_pathways_path: путь для сохранения CSV с путями

    Returns:
        модель с обогащенными путями и DataFrame с результатами
    """
    model = load_model(model_path)
    model, pathways_df = annotate_pathways_from_kegg(
        model,
        output_model_path=output_model_path
    )

    pathways_df.to_csv(output_pathways_path, index=False)
    logging.info(f"Пути сохранены в {output_pathways_path}")

    return model, pathways_df


if __name__ == "__main__":
    model_path = 'data/models/genenametransport2.xml'
    output_model_path = 'data/models/iES1300_prepared3.xml'
    output_pathways_path = 'pathways_annotation3.csv'

    logging.info("=" * 70)
    logging.info("ЗАПУСК ОБОГАЩЕНИЯ МОДЕЛИ")
    logging.info("=" * 70)

    model, pathways_df = prepare_model(
        model_path,
        output_model_path,
        output_pathways_path
    )

    # Показываем результаты
    print("\n" + "=" * 70)
    print("РЕЗУЛЬТАТЫ")
    print("=" * 70)

    updated = pathways_df[pathways_df['Status'] == 'updated']
    print(f"Обновлено реакций: {len(updated)}/{len(pathways_df)}")

    if len(updated) > 0:
        print(f"\nПримеры обновленных реакций (первые 5):")
        for idx, row in updated.head(5).iterrows():
            print(f"\n{row['Reactions']}:")
            print(f"  KEGG ID: {row['KEGG_IDs']}")
            print(f"  Стало: {row['Final_Pathways']}")
