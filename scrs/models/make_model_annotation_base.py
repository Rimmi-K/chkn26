import pandas as pd
import cobra

# Загрузка модели (замените на путь к вашей модели)
# model = read_sbml_model('path_to_your_model.xml')

def create_subsystem_dataframe(model):
    """
    Создает датафрейм с реакциями и их принадлежностью к метаболическим путям.
    
    Parameters:
    -----------
    model : cobra.Model
        Модель COBRApy
    
    Returns:
    --------
    pd.DataFrame
        Датафрейм с реакциями (строки) и метаболическими путями (колонки)
    """
    
    # Шаг 1: Собираем все уникальные метаболические пути
    all_subsystems = set()
    print(str(model.reactions[0].subsystem))
    for reaction in model.reactions:
        
        if reaction.subsystem:  # Проверяем, что subsystem не пустой
            # Разделяем по точке с запятой и убираем пробелы
            subsystems = [s.strip() for s in reaction.subsystem.split(';')]
            print(subsystems)
            all_subsystems.update(subsystems)
    print(all_subsystems)
    # Преобразуем в отсортированный список для удобства
    subsystem_columns = sorted(list(all_subsystems))
    
    # Шаг 2: Создаем словарь для датафрейма
    data_dict = {subsystem: [] for subsystem in subsystem_columns}
    reaction_ids = []
    
    # Шаг 3: Заполняем данные для каждой реакции
    for reaction in model.reactions:
        reaction_ids.append(reaction.id)
        
        # Получаем список путей для текущей реакции
        if reaction.subsystem:
            current_subsystems = set(s.strip() for s in reaction.subsystem.split(';'))
        else:
            current_subsystems = set()
        
        # Для каждого метаболического пути проверяем принадлежность
        for subsystem in subsystem_columns:
            if subsystem in current_subsystems:
                data_dict[subsystem].append(1)
            else:
                data_dict[subsystem].append(0)
    
    # Шаг 4: Создаем датафрейм
    df = pd.DataFrame(data_dict, index=reaction_ids)
    
    return df

if __name__ == "__main__":
# Использование функции
    model = cobra.io.load_json_model("data/models/model_curated.json")
    df = create_subsystem_dataframe(model)

# Сохранение в CSV
    df.to_csv('data/models/subsystem_matrix.csv', index=True)

# Вывод информации о датафрейме
    print(f"Размер датафрейма: {df.shape}")
    print(f"Количество реакций: {df.shape[0]}")
    print(f"Количество уникальных метаболических путей: {df.shape[1]}")
    print("\nПервые 5 строк:")
    print(df.head())
