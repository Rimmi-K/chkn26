import cobra
import riptide
import os
import logging
import re
import pandas as pd
import sys

# Настройка логгирования в файл
log_file = "riptide_log.txt"
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

# Класс-перенаправитель вывода
class StreamToFile:
    def __init__(self, filename):
        self.file = open(filename, "a", encoding="utf-8")

    def write(self, message):
        self.file.write(message)
        self.file.flush()

    def flush(self):
        self.file.flush()

    def close(self):
        self.file.close()


def set_tissue_specific_bounds(model, tissue):
    """
    Устанавливает tissue-specific границы exchange реакций на основе 
    моделей млекопитающих (Recon2/3, MMRN), адаптированных для кур.
    
    Parameters
    ----------
    model : cobra.Model
        Genome-scale metabolic model
    tissue : str
        Тип ткани: 'liver', 'breast', 'leg'
    
    Returns
    -------
    cobra.Model
        Модель с обновленными границами и проверенной feasibility
    """
    
    BASE_AA_UPTAKE = -0.05
    
    amino_acids = [
        'ala__L', 'arg__L', 'asn__L', 'asp__L', 'cys__L', 
        'gln__L', 'glu__L', 'gly', 'his__L', 'ile__L', 
        'leu__L', 'lys__L', 'met__L', 'phe__L', 'pro__L', 
        'ser__L', 'thr__L', 'trp__L', 'tyr__L', 'val__L'
    ]
    
    # Сброс всех AA к базовому уровню
    for aa in amino_acids:
        rxn_id = f"EX_{aa}[e]"
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).bounds = (BASE_AA_UPTAKE, 1000.0)
    
    # === TISSUE-SPECIFIC BOUNDS ===
    
    if tissue == 'liver':
        # Печень: gluconeogenesis hub, Cori cycle, nitrogen metabolism
        model.reactions.get_by_id("EX_lac__L[e]").bounds = (-20.0, 1000.0)
        model.reactions.get_by_id("EX_glc__D[e]").bounds = (-10.0, 1000.0)
        
        if "EX_glu__L[e]" in model.reactions:
            model.reactions.get_by_id("EX_glu__L[e]").bounds = (-2.0, 1000.0)
        
        if "EX_ala__L[e]" in model.reactions:
            model.reactions.get_by_id("EX_ala__L[e]").bounds = (-1.0, 1000.0)
        
        if "EX_o2[e]" in model.reactions:
            model.reactions.get_by_id("EX_o2[e]").bounds = (-50.0, 1000.0)
        
        print(f"  ✓ Liver bounds set: lac(-20), glc(-10), glu(-2), ala(-1), O2(-50)")
    
    elif tissue == 'breast':
        # ГРУДНАЯ МЫШЦА: glycolytic powerhouse, fast-twitch
        
        # КРИТИЧЕСКОЕ: высокий glucose demand для роста
        model.reactions.get_by_id("EX_glc__D[e]").bounds = (-50.0, 1000.0)
        
        # Lactate export (гликолиз → лактат)
        model.reactions.get_by_id("EX_lac__L[e]").bounds = (0.0, 20.0)
        
        # BCAA для protein synthesis
        bcaa_list = ['leu__L', 'ile__L', 'val__L']
        for bcaa in bcaa_list:
            rxn_id = f"EX_{bcaa}[e]"
            if rxn_id in model.reactions:
                model.reactions.get_by_id(rxn_id).bounds = (-0.5, 1000.0)
        
        if "EX_o2[e]" in model.reactions:
            model.reactions.get_by_id("EX_o2[e]").bounds = (-30.0, 1000.0)
        
        print(f"  ✓ Breast bounds set: glc(-50!), lac(0→20), BCAA(-0.5), O2(-30)")
    
    elif tissue == 'leg':
        # МЫШЦА НОГИ: oxidative + glycolytic mix, slow-twitch dominant
        
        model.reactions.get_by_id("EX_glc__D[e]").bounds = (-20.0, 1000.0)
        model.reactions.get_by_id("EX_lac__L[e]").bounds = (-5.0, 10.0)
        
        # BCAA
        bcaa_list = ['leu__L', 'ile__L', 'val__L']
        for bcaa in bcaa_list:
            rxn_id = f"EX_{bcaa}[e]"
            if rxn_id in model.reactions:
                model.reactions.get_by_id(rxn_id).bounds = (-0.3, 1000.0)
        
        # High oxidative capacity
        if "EX_o2[e]" in model.reactions:
            model.reactions.get_by_id("EX_o2[e]").bounds = (-50.0, 1000.0)
        
        # Fatty acid oxidation (если доступно)
        fa_reactions = ['EX_palmitate[e]', 'EX_stearate[e]', 'EX_oleate[e]']
        for fa_rxn in fa_reactions:
            if fa_rxn in model.reactions:
                model.reactions.get_by_id(fa_rxn).bounds = (-0.77, 1000.0)
        
        print(f"  ✓ Leg bounds set: glc(-20), lac(-5→10), BCAA(-0.3), O2(-50), FA(-0.77)")
    
    else:
        print(f"  ⚠ Unknown tissue '{tissue}', using default bounds")
    
    # === КРИТИЧЕСКАЯ ПРОВЕРКА: Оптимизация после установки границ ===
    print(f"\n  Testing model feasibility for {tissue}...")
    
    try:
        solution = model.optimize()
        
        if solution.status == 'optimal':
            biomass_flux = solution.objective_value
            print(f"  ✓ Model is feasible!")
            print(f"  ✓ Biomass objective flux: {biomass_flux:.6f} mmol/gDW/h")
            
            # Логирование для парсинга
            logging.info(f"  {tissue.upper()} post-bounds biomass: {biomass_flux:.6f}")
        else:
            print(f"  ✗ WARNING: Model optimization status = {solution.status}")
            print(f"  ✗ Model may be infeasible with current bounds!")
            logging.warning(f"  {tissue.upper()} infeasible after bounds: {solution.status}")
            
    except Exception as e:
        print(f"  ✗ ERROR during optimization: {e}")
        logging.error(f"  {tissue.upper()} optimization failed: {e}")
    
    print()  # пустая строка для читаемости
    
    return model


def process_group_with_fractions(group_name, model, slow_transcriptome_path, fast_transcriptome_path, output_dir="results"):
    os.makedirs(output_dir, exist_ok=True)

    logging.info(f'\n{"="*60}')
    logging.info(f'Загрузка данных для {group_name}')
    logging.info(f'{"="*60}')
    
    # === ВАЖНО: Установить tissue-specific bounds ПЕРЕД RIPTiDe ===
    logging.info(f'Установка tissue-specific bounds для {group_name}')
    model = set_tissue_specific_bounds(model, group_name)
    
    t_slow = riptide.read_transcription_file(slow_transcriptome_path, header=True, norm=False)
    t_fast = riptide.read_transcription_file(fast_transcriptome_path, header=True, norm=False)

    for fraction in [round(i * 0.05, 2) for i in range(2, 20)]:
        logging.info(f'\nОбработка {group_name} с fraction = {fraction}')
        
        # === ОБРАБОТКА SLOW ===
        logging.info(f'  Обработка group_{group_name}_slow с fraction = {fraction}')
        mcs_slow = riptide.contextualize(model=model, transcriptome=t_slow, fraction=fraction)

        output_slow_path = os.path.join(output_dir, f'mcs_slow_{group_name}_fraction_{fraction:.2f}.xml')
        cobra.io.write_sbml_model(mcs_slow.model, output_slow_path)
        logging.info(f'  Сохранено: {output_slow_path}')
        
        # === ОБРАБОТКА FAST ===
        logging.info(f'  Обработка group_{group_name}_fast с fraction = {fraction}')
        mcs_fast = riptide.contextualize(model=model, transcriptome=t_fast, fraction=fraction)

        output_fast_path = os.path.join(output_dir, f'mcs_fast_{group_name}_fraction_{fraction:.2f}.xml')
        cobra.io.write_sbml_model(mcs_fast.model, output_fast_path)
        logging.info(f'  Сохранено: {output_fast_path}')


def parse_log_to_table(log_path="riptide_log.txt", output_path="riptide_summary.csv"):
    entries = []
    current_group = None
    current_tissue = None
    current_fraction = None
    change = None

    with open(log_path, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()

            # Более гибкая регулярка для парсинга group, tissue и fraction
            match = re.search(r"group_(\w+)_(slow|fast).+fraction\s*=\s*([\d.]+)", line)
            if match:
                current_tissue = match.group(1)
                current_group = match.group(2)
                current_fraction = float(match.group(3))
                continue

            if "Flux through the objective DECREASED to" in line:
                change_match = re.search(r"\(([\d.]+)% change\)", line)
                if change_match:
                    change = float(change_match.group(1))

            if "Context-specific metabolism correlates with transcriptome" in line:
                corr_match = re.search(r"r=([-\d.]+), p=([\d.e-]+)", line)
                if corr_match:
                    r = float(corr_match.group(1))
                    p = float(corr_match.group(2))
                    entries.append([current_group, current_tissue, current_fraction, change, r, p])

    df = pd.DataFrame(entries, columns=["group", "tissue", "fraction", "change", "correlation(r)", "p-value"])
    df.to_csv(output_path, index=False)
    print(f"Таблица сохранена в {output_path}")


if __name__ == "__main__":
    sys.stdout = StreamToFile(log_file)
    sys.stderr = sys.stdout

    output_dir = 'data/models/context_specific'
    tissues = ['breast', 'leg', 'liver']
    transcript_path_template = './data/processed/transcriptomics/{tissue}_{group}.csv'
    json_model_path = 'data/models/iES1300_prepared.json'

    print("Загрузка модели...")
    model = cobra.io.load_json_model(json_model_path)

    print("Проверка модели...")
    solution = model.optimize()
    print(f"Biomass flux: {solution.objective_value}")

    for tissue in tissues:
        slow_path = transcript_path_template.format(tissue=tissue, group="slow")
        fast_path = transcript_path_template.format(tissue=tissue, group="fast")

        # ВАЖНО: Создать копию модели для каждой ткани
        tissue_model = model.copy()
        
        process_group_with_fractions(
            group_name=tissue,
            model=tissue_model,  # используем копию
            slow_transcriptome_path=slow_path,
            fast_transcriptome_path=fast_path,
            output_dir=output_dir
        )

    parse_log_to_table(log_path=log_file, output_path="riptide_summary.csv")
    
    if isinstance(sys.stdout, StreamToFile):
        sys.stdout.close()
    
    print("\n" + "="*60)
    print("ЗАВЕРШЕНО! Проверьте:")
    print("  - riptide_log.txt (детальный лог)")
    print("  - riptide_summary.csv (сводная таблица)")
    print("="*60)
