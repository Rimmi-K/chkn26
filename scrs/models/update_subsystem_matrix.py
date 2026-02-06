#!/usr/bin/env python3
"""
Скрипт для обновления subsystem_matrix.csv с новыми транспортными реакциями из curated_model51.xml
Классифицирует транспортные реакции по компартментам метаболитов.
"""

import cobra
import pandas as pd
from pathlib import Path

# Пути к файлам
MODEL_PATH = Path("data/models/curated_model51.xml")
SUBSYSTEM_MATRIX_PATH = Path("data/models/subsystem_matrix.csv")
OUTPUT_PATH = Path("data/models/subsystem_matrix_updated.csv")

# Определение транспортных категорий по компартментам
TRANSPORT_CATEGORIES = {
    "Transport, Mitochondrial": ["c", "m"],        # цитозоль <-> митохондрия
    "Transport, Peroxisomal": ["c", "x"],          # цитозоль <-> пероксисома
    "Transport, Endoplasmic Reticulum": ["c", "r"], # цитозоль <-> ЭПР
    "Transport, Golgi": ["c", "g"],                # цитозоль <-> аппарат Гольджи
    "Transport, Lysosomal": ["c", "l"],            # цитозоль <-> лизосома
    "Transport, Nuclear": ["c", "n"],              # цитозоль <-> ядро
    "Transport, Endosomal": ["c", "e"],            # цитозоль <-> эндосома
    "Extracellular Transport": ["c", "s"],         # цитозоль <-> внеклеточное пространство
}


def load_model(model_path):
    """Загрузка метаболической модели"""
    print(f"Загрузка модели из {model_path}...")
    model = cobra.io.read_sbml_model(str(model_path))
    print(f"  Загружено {len(model.reactions)} реакций")
    return model


def load_subsystem_matrix(matrix_path):
    """Загрузка матрицы subsystem"""
    print(f"Загрузка матрицы subsystem из {matrix_path}...")
    df = pd.read_csv(matrix_path, index_col=0)
    print(f"  Загружено {len(df)} реакций в {len(df.columns)} путях")
    return df


def get_compartments(reaction):
    """Получить список уникальных компартментов метаболитов реакции"""
    compartments = set()
    for met in reaction.metabolites:
        compartments.add(met.compartment)
    return sorted(compartments)


def classify_transport_reaction(reaction):
    """
    Классифицировать транспортную реакцию по типу транспорта.

    Returns:
        str или None: название транспортной категории или None
    """
    compartments = get_compartments(reaction)

    # Реакция должна затрагивать как минимум 2 компартмента, чтобы быть транспортной
    if len(compartments) < 2:
        return None

    # Проверяем каждую транспортную категорию
    for category, comp_pattern in TRANSPORT_CATEGORIES.items():
        # Проверяем, соответствует ли набор компартментов паттерну
        if set(comp_pattern) == set(compartments):
            return category

    # Если компартментов больше 2, пробуем найти основной транспорт
    # (например, c-m-x классифицируем как митохондриальный, если есть c и m)
    for category, comp_pattern in TRANSPORT_CATEGORIES.items():
        if set(comp_pattern).issubset(set(compartments)):
            return category

    # Если не подошло ни одно правило, считаем просто "Transport"
    return "Transport, Other"


def find_new_reactions(model, subsystem_matrix):
    """Найти реакции, которых нет в subsystem_matrix"""
    existing_reactions = set(subsystem_matrix.index)
    all_reactions = set(rxn.id for rxn in model.reactions)
    new_reactions = all_reactions - existing_reactions

    print(f"\nНайдено {len(new_reactions)} новых реакций не в subsystem_matrix")
    return new_reactions


def analyze_new_reactions(model, new_reaction_ids):
    """
    Проанализировать новые реакции и классифицировать транспортные.

    Returns:
        dict: {reaction_id: transport_category}
    """
    transport_reactions = {}
    non_transport_reactions = []

    for rxn_id in new_reaction_ids:
        rxn = model.reactions.get_by_id(rxn_id)
        transport_type = classify_transport_reaction(rxn)

        if transport_type:
            transport_reactions[rxn_id] = transport_type
        else:
            non_transport_reactions.append(rxn_id)

    print(f"\nКлассификация новых реакций:")
    print(f"  Транспортных реакций: {len(transport_reactions)}")
    print(f"  Не-транспортных реакций: {len(non_transport_reactions)}")

    # Статистика по типам транспорта
    if transport_reactions:
        print("\nРаспределение по типам транспорта:")
        from collections import Counter
        transport_counts = Counter(transport_reactions.values())
        for transport_type, count in sorted(transport_counts.items()):
            print(f"  {transport_type}: {count}")

    return transport_reactions, non_transport_reactions


def update_subsystem_matrix(subsystem_matrix, transport_reactions):
    """
    Обновить subsystem_matrix с новыми транспортными реакциями.

    Args:
        subsystem_matrix: pandas DataFrame с текущей матрицей
        transport_reactions: dict {reaction_id: transport_category}

    Returns:
        pandas DataFrame: обновленная матрица
    """
    # Создаем копию
    updated_matrix = subsystem_matrix.copy()

    # Убедимся, что все транспортные колонки существуют
    for transport_type in set(transport_reactions.values()):
        if transport_type not in updated_matrix.columns:
            print(f"⚠ Предупреждение: колонка '{transport_type}' не найдена в матрице!")
            # Добавляем новую колонку
            updated_matrix[transport_type] = 0

    # Добавляем новые реакции
    for rxn_id, transport_type in transport_reactions.items():
        # Создаем новую строку с нулями
        new_row = pd.Series(0, index=updated_matrix.columns, name=rxn_id)
        # Устанавливаем 1 для соответствующего транспортного пути
        new_row[transport_type] = 1
        # Добавляем строку к матрице
        updated_matrix = pd.concat([updated_matrix, new_row.to_frame().T])

    print(f"\nДобавлено {len(transport_reactions)} реакций в subsystem_matrix")
    print(f"Новый размер матрицы: {len(updated_matrix)} реакций × {len(updated_matrix.columns)} путей")

    return updated_matrix


def show_examples(model, transport_reactions, n=5):
    """Показать примеры классифицированных реакций"""
    print(f"\n{'='*80}")
    print(f"Примеры классифицированных реакций (первые {n}):")
    print(f"{'='*80}")

    for i, (rxn_id, transport_type) in enumerate(list(transport_reactions.items())[:n]):
        rxn = model.reactions.get_by_id(rxn_id)
        compartments = get_compartments(rxn)

        print(f"\n{i+1}. {rxn_id}")
        print(f"   Название: {rxn.name}")
        print(f"   Реакция: {rxn.reaction}")
        print(f"   Компартменты: {', '.join(compartments)}")
        print(f"   Классификация: {transport_type}")

        # Показать метаболиты по компартментам
        print(f"   Метаболиты:")
        for comp in sorted(compartments):
            mets_in_comp = [f"{m.id} ({coef:+.1f})"
                           for m, coef in rxn.metabolites.items()
                           if m.compartment == comp]
            print(f"     [{comp}]: {', '.join(mets_in_comp)}")


def main():
    """Основная функция"""
    print("="*80)
    print("Обновление subsystem_matrix с новыми транспортными реакциями")
    print("="*80)

    # 1. Загрузка данных
    model = load_model(MODEL_PATH)
    subsystem_matrix = load_subsystem_matrix(SUBSYSTEM_MATRIX_PATH)

    # 2. Поиск новых реакций
    new_reaction_ids = find_new_reactions(model, subsystem_matrix)

    if not new_reaction_ids:
        print("\n✓ Все реакции из модели уже есть в subsystem_matrix!")
        return

    # 3. Анализ и классификация
    transport_reactions, non_transport = analyze_new_reactions(model, new_reaction_ids)

    # 4. Показать примеры
    if transport_reactions:
        show_examples(model, transport_reactions, n=10)

    # 5. Сохранение информации о не-транспортных реакциях
    if non_transport:
        print(f"\n{'='*80}")
        print(f"⚠ ВНИМАНИЕ: Найдено {len(non_transport)} не-транспортных реакций")
        print(f"{'='*80}")
        print("Эти реакции НЕ будут добавлены в матрицу (требуют ручной классификации):")

        non_transport_file = Path("data/models/non_transport_reactions.txt")
        with open(non_transport_file, "w") as f:
            f.write("# Не-транспортные реакции, требующие ручной классификации\n")
            f.write("# Формат: reaction_id | name | reaction | compartments\n\n")

            for rxn_id in sorted(non_transport)[:20]:  # Показываем первые 20
                rxn = model.reactions.get_by_id(rxn_id)
                compartments = get_compartments(rxn)
                print(f"  {rxn_id}: {rxn.name} [{', '.join(compartments)}]")
                f.write(f"{rxn_id} | {rxn.name} | {rxn.reaction} | {compartments}\n")

            if len(non_transport) > 20:
                print(f"  ... и еще {len(non_transport)-20} реакций")
                for rxn_id in sorted(non_transport)[20:]:
                    rxn = model.reactions.get_by_id(rxn_id)
                    compartments = get_compartments(rxn)
                    f.write(f"{rxn_id} | {rxn.name} | {rxn.reaction} | {compartments}\n")

        print(f"\nПолный список сохранен в: {non_transport_file}")

    # 6. Обновление матрицы
    if transport_reactions:
        updated_matrix = update_subsystem_matrix(subsystem_matrix, transport_reactions)

        # 7. Сохранение
        print(f"\nСохранение обновленной матрицы в {OUTPUT_PATH}...")
        updated_matrix.to_csv(OUTPUT_PATH)
        print(f"✓ Готово!")

        print(f"\n{'='*80}")
        print("Резюме:")
        print(f"  Исходная матрица: {len(subsystem_matrix)} реакций")
        print(f"  Добавлено транспортных: {len(transport_reactions)}")
        print(f"  Пропущено не-транспортных: {len(non_transport)}")
        print(f"  Итоговая матрица: {len(updated_matrix)} реакций")
        print(f"{'='*80}")
    else:
        print("\n⚠ Нет транспортных реакций для добавления")


if __name__ == "__main__":
    main()
