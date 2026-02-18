import re
import sys
import os
import logging

import cobra
import pandas as pd
import riptide

LOG_FILE = "riptide_log.txt"
logging.basicConfig(filename=LOG_FILE, level=logging.INFO, format='%(message)s')


class StreamToFile:
    def __init__(self, filename: str):
        self.file = open(filename, "a", encoding="utf-8")

    def write(self, message: str):
        self.file.write(message)
        self.file.flush()

    def flush(self):
        self.file.flush()

    def close(self):
        self.file.close()


def process_group_with_fractions(
    group_name: str,
    model,
    slow_transcriptome_path: str,
    fast_transcriptome_path: str,
    output_dir: str,
):
    model.objective = "BIOMASS_maintenance"
    model.optimize()
    os.makedirs(output_dir, exist_ok=True)

    t_slow = riptide.read_transcription_file(slow_transcriptome_path, header=True, norm=False)
    t_fast = riptide.read_transcription_file(fast_transcriptome_path, header=True, norm=False)

    for fraction in [round(i * 0.05, 2) for i in range(2, 20)]:
        for growth, transcriptome in (("slow", t_slow), ("fast", t_fast)):
            logging.info(f"\n  Processing group_{group_name}_{growth} with fraction={fraction}")
            mcs = riptide.contextualize(model=model, transcriptome=transcriptome, fraction=fraction)
            out_path = os.path.join(output_dir, f"mcs_{growth}_{group_name}_fraction_{fraction:.2f}.json")
            cobra.io.save_json_model(mcs.model, out_path)
            logging.info(f"  Saved: {out_path}")


def parse_log_to_table(log_path: str = LOG_FILE, output_path: str = "riptide_summary.csv") -> pd.DataFrame:
    group_re     = re.compile(r"Processing\s+group_(\w+)_(slow|fast)\s+with\s+fraction\s*=\s*([\d.]+)")
    reactions_re = re.compile(r"Reactions pruned to (\d+) from (\d+) \(([\d.]+)% change\)")
    metabolites_re = re.compile(r"Metabolites pruned to (\d+) from (\d+) \(([\d.]+)% change\)")
    flux_re      = re.compile(r"Flux through the objective DECREASED to ~?([-\d.eE+]+) from ~?([-\d.eE+]+) \(([\d.]+)% change\)")
    corr_re      = re.compile(r"Context-specific metabolism correlates with transcriptome \(r=([-\d.]+), p([<=>]?)([\d.eE-]+)")
    completed_re = re.compile(r"RIPTiDe completed in (\d+) seconds")
    saved_re     = re.compile(r"Saved:\s*(.+)")

    entries, current = [], None

    def finalize():
        if current is not None:
            entries.append(current.copy())

    with open(log_path, "r", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue

            m = group_re.search(line)
            if m:
                finalize()
                current = {"tissue": m.group(1), "group": m.group(2), "fraction": float(m.group(3))}
                continue

            if current is None:
                continue

            m = reactions_re.search(line)
            if m:
                current.update(reactions_pruned_to=int(m.group(1)), reactions_pruned_from=int(m.group(2)), reactions_change=float(m.group(3)))
                continue

            m = metabolites_re.search(line)
            if m:
                current.update(metabolites_pruned_to=int(m.group(1)), metabolites_pruned_from=int(m.group(2)), metabolites_change=float(m.group(3)))
                continue

            m = flux_re.search(line)
            if m:
                current.update(objective_after=float(m.group(1)), objective_before=float(m.group(2)), change=float(m.group(3)))
                continue

            m = corr_re.search(line)
            if m:
                current.update(**{"correlation(r)": float(m.group(1)), "p_relation": m.group(2) or "=", "p-value": float(m.group(3))})
                continue

            m = completed_re.search(line)
            if m:
                current["seconds"] = int(m.group(1))
                continue

            m = saved_re.search(line)
            if m:
                current["output_path"] = m.group(1).strip()
                finalize()
                current = None

    finalize()

    columns = [
        "group", "tissue", "fraction", "change", "correlation(r)", "p-value", "p_relation",
        "objective_after", "objective_before",
        "reactions_pruned_to", "reactions_pruned_from", "reactions_change",
        "metabolites_pruned_to", "metabolites_pruned_from", "metabolites_change",
        "seconds", "output_path",
    ]
    df = pd.DataFrame(entries).reindex(columns=columns)
    df.to_csv(output_path, index=False)
    print(f"Summary saved: {output_path}")
    return df


if __name__ == "__main__":
    sys.stdout = StreamToFile(LOG_FILE)
    sys.stderr = sys.stdout

    output_dir = 'data/models/context_specific'
    transcript_template = './data/transcriptomics/{tissue}_{group}.csv'
    model_path = 'data/models/iES1470.json'

    try:
        for tissue in ('breast', 'leg', 'liver'):
            model = cobra.io.load_json_model(model_path)
            model.optimize()

            if tissue == "liver":
                model.reactions.get_by_id("EX_lac__L_e").bounds = (-5, 0)

            process_group_with_fractions(
                group_name=tissue,
                model=model,
                slow_transcriptome_path=transcript_template.format(tissue=tissue, group="slow"),
                fast_transcriptome_path=transcript_template.format(tissue=tissue, group="fast"),
                output_dir=output_dir,
            )

        parse_log_to_table(log_path=LOG_FILE, output_path="riptide_summary.csv")
    finally:
        if isinstance(sys.stdout, StreamToFile):
            sys.stdout.close()
        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)
