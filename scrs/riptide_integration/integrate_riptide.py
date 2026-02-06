import cobra
import riptide
import os
import logging
import re
import pandas as pd
import sys

log_file = "riptide_log.txt"
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

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

def process_group_with_fractions(group_name, model, slow_transcriptome_path, fast_transcriptome_path, output_dir):
    model.objective = "BIOMASS_maintenance"
    model.optimize()
    os.makedirs(output_dir, exist_ok=True)        
    t_slow = riptide.read_transcription_file(slow_transcriptome_path, header=True, norm=False)
    t_fast = riptide.read_transcription_file(fast_transcriptome_path, header=True, norm=False)

    for fraction in [round(i * 0.05, 2) for i in range(2, 20)]:
        logging.info(f'\nProcessing {group_name} с fraction = {fraction}')
        
        # === ОБРАБОТКА SLOW ===
        logging.info(f'  Processing group_{group_name}_slow с fraction = {fraction}')
        mcs_slow = riptide.contextualize(model=model, transcriptome=t_slow, fraction=fraction)

        output_slow_path = os.path.join(output_dir, f'mcs_slow_{group_name}_fraction_{fraction:.2f}.json')
        cobra.io.save_json_model(mcs_slow.model, output_slow_path)
        logging.info(f'  Saved: {output_slow_path}')
        
        # === ОБРАБОТКА FAST ===
        logging.info(f'  Processing group_{group_name}_fast с fraction = {fraction}')
        mcs_fast = riptide.contextualize(model=model, transcriptome=t_fast, fraction=fraction)
        
        output_fast_path = os.path.join(output_dir, f'mcs_fast_{group_name}_fraction_{fraction:.2f}.json')
        cobra.io.save_json_model(mcs_fast.model, output_fast_path)
        logging.info(f'  Saved: {output_fast_path}')


def parse_log_to_table(log_path="riptide_log.txt", output_path="riptide_summary.csv"):
    entries = []
    current = None

    group_re = re.compile(r"Processing\s+group_(\w+)_(slow|fast)\s+с\s+fraction\s*=\s*([\d.]+)")
    reactions_re = re.compile(r"Reactions pruned to (\d+) from (\d+) \(([\d.]+)% change\)")
    metabolites_re = re.compile(r"Metabolites pruned to (\d+) from (\d+) \(([\d.]+)% change\)")
    flux_re = re.compile(
        r"Flux through the objective DECREASED to ~?([-\d.eE+]+) from ~?([-\d.eE+]+) \(([\d.]+)% change\)"
    )
    corr_re = re.compile(r"Context-specific metabolism correlates with transcriptome \(r=([-\d.]+), p([<=>]?)([\d.eE-]+)")
    completed_re = re.compile(r"RIPTiDe completed in (\d+) seconds")
    saved_re = re.compile(r"Saved:\s*(.+)")

    def finalize_entry():
        if current is None:
            return
        entries.append(current.copy())

    with open(log_path, "r", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue

            match = group_re.search(line)
            if match:
                if current is not None:
                    finalize_entry()
                current = {
                    "tissue": match.group(1),
                    "group": match.group(2),
                    "fraction": float(match.group(3)),
                }
                continue

            if current is None:
                continue

            match = reactions_re.search(line)
            if match:
                current["reactions_pruned_to"] = int(match.group(1))
                current["reactions_pruned_from"] = int(match.group(2))
                current["reactions_change"] = float(match.group(3))
                continue

            match = metabolites_re.search(line)
            if match:
                current["metabolites_pruned_to"] = int(match.group(1))
                current["metabolites_pruned_from"] = int(match.group(2))
                current["metabolites_change"] = float(match.group(3))
                continue

            match = flux_re.search(line)
            if match:
                current["objective_after"] = float(match.group(1))
                current["objective_before"] = float(match.group(2))
                current["change"] = float(match.group(3))
                continue

            match = corr_re.search(line)
            if match:
                current["correlation(r)"] = float(match.group(1))
                current["p_relation"] = match.group(2) or "="
                current["p-value"] = float(match.group(3))
                continue

            match = completed_re.search(line)
            if match:
                current["seconds"] = int(match.group(1))
                continue

            match = saved_re.search(line)
            if match:
                current["output_path"] = match.group(1).strip()
                finalize_entry()
                current = None

    if current is not None:
        finalize_entry()

    columns = [
        "group",
        "tissue",
        "fraction",
        "change",
        "correlation(r)",
        "p-value",
        "p_relation",
        "objective_after",
        "objective_before",
        "reactions_pruned_to",
        "reactions_pruned_from",
        "reactions_change",
        "metabolites_pruned_to",
        "metabolites_pruned_from",
        "metabolites_change",
        "seconds",
        "output_path",
    ]
    df = pd.DataFrame(entries)
    df = df.reindex(columns=columns)
    df.to_csv(output_path, index=False)
    print(f"Таблица сохранена в {output_path}")


if __name__ == "__main__":
    sys.stdout = StreamToFile(log_file)
    sys.stderr = sys.stdout

    output_dir = 'data/models/context_specific'
    tissues = ['breast', 'leg', 'liver']
    transcript_path_template = './data/processed/transcriptomics/{tissue}_{group}.csv'
    model_path = 'data/models/curated_model55.json'

    try:
        for tissue in tissues:
            slow_path = transcript_path_template.format(tissue=tissue, group="slow")
            fast_path = transcript_path_template.format(tissue=tissue, group="fast")
            model = cobra.io.load_json_model(model_path)
            solution = model.optimize()

            if tissue == "liver":
                model.reactions.get_by_id("EX_lac__L_e").bounds = (-5, 0)

            process_group_with_fractions(
                group_name=tissue,
                model=model,
                slow_transcriptome_path=slow_path,
                fast_transcriptome_path=fast_path,
                output_dir=output_dir
            )

        parse_log_to_table(log_path=log_file, output_path="riptide_summary.csv")
    finally:
        if isinstance(sys.stdout, StreamToFile):
            sys.stdout.close()

        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)

