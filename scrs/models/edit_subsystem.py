import cobra

# Маппинг для конкретных реакций (шаг 1)
REACTION_SUBSYSTEM_FIXES = {
    'r0340g': 'Glycolysis / Gluconeogenesis',
    'r0384': 'Citric Acid Cycle',
    'r1391': 'Starch And Sucrose Metabolism',
    'r1392': 'Starch And Sucrose Metabolism',
    'r0510': 'Biosynthesis Of Unsaturated Fatty Acids',
    'r1174': 'Steroid Metabolism',
    'r1177': 'Steroid Metabolism',
    'r1167': 'Glycerophospholipid Metabolism',
    'r1479': 'Fatty Acid Oxidation',
    'R_group_phosphotase_1': 'Nucleotide Metabolism'
}

# Словарь группировки (шаг 2)
SUBSYSTEM_NORMALIZATION = {
    # Энергетический метаболизм
    'Citrate Cycle (tca Cycle)': 'Citric Acid Cycle',
    'Glycolysis/gluconeogenesis': 'Glycolysis / Gluconeogenesis',
    
    # Аминокислоты - Глицин/Серин/Треонин/Аланин
    'Glycine': 'Glycine, Serine And Threonine Metabolism',
    'Serine': 'Glycine, Serine And Threonine Metabolism',
    'Alanine And Threonine Metabolism': 'Glycine, Serine And Threonine Metabolism',
    
    # Аминокислоты - Аланин/Аспартат/Глутамат
    'Alanine And Aspartate Metabolism': 'Alanine, Aspartate And Glutamate Metabolism',
    'Glutamate Metabolism': 'Alanine, Aspartate And Glutamate Metabolism',
    
    # BCAA (Разветвленные аминокислоты)
    'Valine, Leucine And Isoleucine Degradation': 'Valine, Leucine And Isoleucine Metabolism',
    'Valine, Leucine And Isoleucine Biosynthesis': 'Valine, Leucine And Isoleucine Metabolism',
    'Valine': 'Valine, Leucine And Isoleucine Metabolism',
    'Leucine': 'Valine, Leucine And Isoleucine Metabolism',
    'And Isoleucine Metabolism': 'Valine, Leucine And Isoleucine Metabolism',
    
    # Серосодержащие аминокислоты
    'Methionine And Cysteine Metabolism': 'Cysteine And Methionine Metabolism',
    'Cysteine Metabolism': 'Cysteine And Methionine Metabolism',
    
    # Нуклеотиды (объединяем анаболизм/катаболизм)
    'Purine Synthesis': 'Purine Metabolism',
    'Purine Catabolism': 'Purine Metabolism',
    'Pyrimidine Synthesis': 'Pyrimidine Metabolism',
    'Pyrimidine Catabolism': 'Pyrimidine Metabolism',
    
    # Сахара
    'Amino Sugar And Nucleotide Sugar Metabolism': 'Aminosugar Metabolism',
    'Biosynthesis Of Nucleotide Sugars': 'Biosynthesis Of Various Nucleotide Sugars',
    
    # CoA метаболизм
    'Coa Synthesis': 'Pantothenate And Coa Biosynthesis',
    'Coa Catabolism': 'Pantothenate And Coa Biosynthesis',
    
    # Витамины
    'Vitamin A Metabolism': 'Retinol Metabolism',
    
    # Heme метаболизм
    'Heme Synthesis': 'Heme Metabolism',
    'Heme Degradation': 'Heme Metabolism',
    'Heme': 'Heme Metabolism',
    
    # Жирные кислоты - объединить дубликаты
    'Fatty Acid Synthesis': 'Fatty Acid Biosynthesis',
    'Fatty Acid Degradation': 'Fatty Acid Oxidation',
    
    # Желчные кислоты
    'Primary Bile Acid Biosynthesis': 'Bile Acid Synthesis',
    
    # Линолевая кислота
    'Linoleic Acid Metabolism': 'Linoleate Metabolism',

    # Стероиды
    'Steroid Biosynthesis': 'Steroid Metabolism',
    'Steroid Hormone Biosynthesis': 'Steroid Metabolism',
    
    # Лизин
    'Lysine Biosynthesis': 'Lysine Metabolism',
    'Lysine Degradation': 'Lysine Metabolism',
    
    # Убихинон
    'Ubiquinone Synthesis': 'Ubiquinone And Other Terpenoid-quinone Biosynthesis',
    
    # Сера
    'Sulfur Cycle': 'Sulfur Metabolism',
    
    # Компартменты → Транспорт
    'Mitochondrial': 'Mitochondrial Transport',
    'Peroxisomal': 'Peroxisomal Transport',
    'Endoplasmic Reticular': 'Endoplasmic Reticular Transport',
    'Lysosomal': 'Lysosomal Transport',
    'Nuclear': 'Nuclear Transport',
    'Golgi Apparatus': 'Golgi Transport',
    'Endosomal': 'Endosomal Transport',

        # Selenium
    'Selenoamino Acid Metabolism': 'Selenium Metabolism',
    'Selenocompound Metabolism': 'Selenium Metabolism',
}

SUBSYSTEMS_TO_REMOVE = {
    # Общие/мусорные категории
    'Metabolic Pathways',
    'Biosynthesis Of Secondary Metabolites',
    'Microbial Metabolism In Diverse Environments',
    'Miscellaneous',
    'Unassigned',
    'Exchange/demand Reaction',
    'Extracellular exchange',
    'Transport',
    
    # Бактериальные пути
    'Streptomycin Biosynthesis',
    'Neomycin, Kanamycin And Gentamicin Biosynthesis',
    'Carbapenem Biosynthesis',
    'Monobactam Biosynthesis',
    'Mycolic Acid Biosynthesis',
    'Novobiocin Biosynthesis',
    'Teichoic Acid Biosynthesis',
    
    # Растительные/грибные пути
    'Aflatoxin Biosynthesis',
    'Betalain Biosynthesis',
    'Glucosinolate Biosynthesis',
    'Tropane, Piperidine And Pyridine Alkaloid Biosynthesis',
    'Isoquinoline Alkaloid Biosynthesis',
    'Biosynthesis Of Various Plant Secondary Metabolites',
    
    # Неспецифичные для животных
    'Styrene Degradation',


}

def normalize_subsystems(model):
    """
    Нормализует подсистемы в метаболической модели:
    1. Заменяет пустые subsystem на значения из REACTION_SUBSYSTEM_FIXES
    2. Нормализует названия через SUBSYSTEM_NORMALIZATION
    3. Удаляет мусорные категории из SUBSYSTEMS_TO_REMOVE
    4. Удаляет дубликаты в пределах одной реакции
    """
    
    stats = {
        'fixed_empty': 0,
        'normalized_subsystems': 0,
        'removed_garbage': 0,
        'reactions_changed': 0,
        'unchanged': 0
    }
    
    for reaction in model.reactions:
        original_subsystem = reaction.subsystem
        changed = False
        
        # Шаг 1: Исправить пустые/Unassigned реакции
        if (not reaction.subsystem or reaction.subsystem == 'Unassigned') and reaction.id in REACTION_SUBSYSTEM_FIXES:
            reaction.subsystem = REACTION_SUBSYSTEM_FIXES[reaction.id]
            stats['fixed_empty'] += 1
            changed = True
            print(f"✓ Исправлено пустое: {reaction.id} → {reaction.subsystem}")
        
        # Удалить 'Unassigned' если он там есть (заменить на пустую строку)
        elif reaction.subsystem == 'Unassigned':
            reaction.subsystem = ''
            changed = True
            stats['removed_garbage'] += 1
        
        # Шаг 2: Нормализовать существующие подсистемы
        if reaction.subsystem and reaction.subsystem != '':
            subsystems = [s.strip() for s in reaction.subsystem.split(';')]
            normalized = []
            
            for sub in subsystems:
                # Пропустить мусорные категории
                if sub in SUBSYSTEMS_TO_REMOVE:
                    stats['removed_garbage'] += 1
                    changed = True
                    continue
                
                # Нормализовать название
                if sub in SUBSYSTEM_NORMALIZATION:
                    normalized.append(SUBSYSTEM_NORMALIZATION[sub])
                    stats['normalized_subsystems'] += 1
                    changed = True
                elif sub:  # не добавляем пустые строки
                    normalized.append(sub)
            
            # Убрать дубликаты, сохраняя порядок
            unique_normalized = list(dict.fromkeys(normalized))
            
            # Обновить subsystem
            if unique_normalized:  # если осталось что-то после фильтрации
                new_subsystem = '; '.join(unique_normalized)
            else:  # если все подсистемы были мусорными
                new_subsystem = ''
            
            if new_subsystem != reaction.subsystem:
                if changed and original_subsystem not in ['', 'Unassigned']:
                    print(f"  Нормализовано: {reaction.id}")
                    print(f"    ДО:  {original_subsystem}")
                    print(f"    ПОСЛЕ: {new_subsystem}")
                reaction.subsystem = new_subsystem
                stats['reactions_changed'] += 1
        
        if not changed:
            stats['unchanged'] += 1
    
    return stats


def analyze_subsystems(model, title=""):
    """Анализ подсистем в модели"""
    print(f"\n{'='*60}")
    print(f"{title}")
    print(f"{'='*60}")
    
    subs = set()
    unassigned = []
    empty = []
    
    for r in model.reactions:
        if r.subsystem == 'Unassigned':
            unassigned.append(r.id)
        elif not r.subsystem or r.subsystem == '':
            empty.append(r.id)
        elif r.subsystem:
            for s in r.subsystem.split(";"):
                subs.add(s.strip())
    
    print(f"\nСтатистика:")
    print(f"  Уникальных подсистем: {len(subs)}")
    print(f"  Реакций с 'Unassigned': {len(unassigned)}")
    print(f"  Реакций с пустой subsystem: {len(empty)}")
    print(f"  Всего реакций: {len(model.reactions)}")
    
    if unassigned:
        print(f"\n⚠ Реакции с 'Unassigned' ({len(unassigned)}):")
        for rid in unassigned[:10]:
            r = model.reactions.get_by_id(rid)
            print(f"  - {rid}: {r.name}")
        if len(unassigned) > 10:
            print(f"  ... и еще {len(unassigned) - 10}")
    
    if empty:
        print(f"\n⚠ Реакции с пустой subsystem ({len(empty)}):")
        for rid in empty[:10]:
            r = model.reactions.get_by_id(rid)
            print(f"  - {rid}: {r.name}")
        if len(empty) > 10:
            print(f"  ... и еще {len(empty) - 10}")
    
    return subs


# ============ ОСНОВНОЙ КОД ============

print("Загрузка модели...")
model = cobra.io.load_json_model("iES1300_prepared.json")

# Анализ ДО
subs_before = analyze_subsystems(model, "ДО НОРМАЛИЗАЦИИ")

# Нормализация
print(f"\n{'='*60}")
print("ПРОЦЕСС НОРМАЛИЗАЦИИ")
print(f"{'='*60}\n")
stats = normalize_subsystems(model)

print(f"\n{'='*60}")
print("ИТОГОВАЯ СТАТИСТИКА")
print(f"{'='*60}")
print(f"  ✓ Исправлено пустых subsystem: {stats['fixed_empty']}")
print(f"  ✓ Нормализовано названий подсистем: {stats['normalized_subsystems']}")
print(f"  🗑 Удалено мусорных категорий: {stats['removed_garbage']}")
print(f"  ✓ Изменено реакций: {stats['reactions_changed']}")
print(f"  ○ Без изменений: {stats['unchanged']}")

# Анализ ПОСЛЕ
subs_after = analyze_subsystems(model, "ПОСЛЕ НОРМАЛИЗАЦИИ")

# Сравнение
print(f"\n{'='*60}")
print("СРАВНЕНИЕ")
print(f"{'='*60}")
print(f"Подсистем было: {len(subs_before)}")
print(f"Подсистем стало: {len(subs_after)}")
reduction = len(subs_before) - len(subs_after)
if len(subs_before) > 0:
    reduction_pct = 100 * reduction / len(subs_before)
    print(f"Сокращение: {reduction} ({reduction_pct:.1f}%)")

# Показать новые подсистемы (если появились)
new_subs = subs_after - subs_before
if new_subs:
    print(f"\n✨ Новые подсистемы ({len(new_subs)}):")
    for s in sorted(new_subs):
        print(f"  + {s}")

# Показать удаленные подсистемы
removed_subs = subs_before - subs_after
if removed_subs:
    print(f"\n🗑 Удаленные подсистемы ({len(removed_subs)}):")
    for s in sorted(removed_subs):
        print(f"  - {s}")

# Финальный список подсистем
print(f"\n{'='*60}")
print("ФИНАЛЬНЫЙ СПИСОК ПОДСИСТЕМ")
print(f"{'='*60}")
sorted_subs = sorted(list(subs_after))
for i, sub in enumerate(sorted_subs, 1):
    print(f"{i:3d}. {sub}")

# Сохранение
output_file = "iES1300_normalized.json"
print(f"\n{'='*60}")
print(f"Сохранение модели в {output_file}...")
cobra.io.save_json_model(model, output_file)
print("✓ Готово!")