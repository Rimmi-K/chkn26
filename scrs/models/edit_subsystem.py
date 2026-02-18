import cobra

REACTION_SUBSYSTEM_FIXES = {
    'r0340g': 'Glycolysis / Gluconeogenesis',
    'r0384':  'Citric Acid Cycle',
    'r1391':  'Starch And Sucrose Metabolism',
    'r1392':  'Starch And Sucrose Metabolism',
    'r0510':  'Biosynthesis Of Unsaturated Fatty Acids',
    'r1174':  'Steroid Metabolism',
    'r1177':  'Steroid Metabolism',
    'r1167':  'Glycerophospholipid Metabolism',
    'r1479':  'Fatty Acid Oxidation',
    'R_group_phosphotase_1': 'Nucleotide Metabolism',
}

SUBSYSTEM_NORMALIZATION = {
    'Citrate Cycle (tca Cycle)':              'Citric Acid Cycle',
    'Glycolysis/gluconeogenesis':             'Glycolysis / Gluconeogenesis',

    'Glycine':                                'Glycine, Serine And Threonine Metabolism',
    'Serine':                                 'Glycine, Serine And Threonine Metabolism',
    'Alanine And Threonine Metabolism':       'Glycine, Serine And Threonine Metabolism',

    'Alanine And Aspartate Metabolism':       'Alanine, Aspartate And Glutamate Metabolism',
    'Glutamate Metabolism':                   'Alanine, Aspartate And Glutamate Metabolism',

    'Valine, Leucine And Isoleucine Degradation':  'Valine, Leucine And Isoleucine Metabolism',
    'Valine, Leucine And Isoleucine Biosynthesis': 'Valine, Leucine And Isoleucine Metabolism',
    'Valine':                                 'Valine, Leucine And Isoleucine Metabolism',
    'Leucine':                                'Valine, Leucine And Isoleucine Metabolism',
    'And Isoleucine Metabolism':              'Valine, Leucine And Isoleucine Metabolism',

    'Methionine And Cysteine Metabolism':     'Cysteine And Methionine Metabolism',
    'Cysteine Metabolism':                    'Cysteine And Methionine Metabolism',

    'Purine Synthesis':                       'Purine Metabolism',
    'Purine Catabolism':                      'Purine Metabolism',
    'Pyrimidine Synthesis':                   'Pyrimidine Metabolism',
    'Pyrimidine Catabolism':                  'Pyrimidine Metabolism',

    'Amino Sugar And Nucleotide Sugar Metabolism': 'Aminosugar Metabolism',
    'Biosynthesis Of Nucleotide Sugars':      'Biosynthesis Of Various Nucleotide Sugars',

    'Coa Synthesis':                          'Pantothenate And Coa Biosynthesis',
    'Coa Catabolism':                         'Pantothenate And Coa Biosynthesis',

    'Vitamin A Metabolism':                   'Retinol Metabolism',

    'Heme Synthesis':                         'Heme Metabolism',
    'Heme Degradation':                       'Heme Metabolism',
    'Heme':                                   'Heme Metabolism',

    'Fatty Acid Synthesis':                   'Fatty Acid Biosynthesis',
    'Fatty Acid Degradation':                 'Fatty Acid Oxidation',

    'Primary Bile Acid Biosynthesis':         'Bile Acid Synthesis',
    'Linoleic Acid Metabolism':               'Linoleate Metabolism',

    'Steroid Biosynthesis':                   'Steroid Metabolism',
    'Steroid Hormone Biosynthesis':           'Steroid Metabolism',

    'Lysine Biosynthesis':                    'Lysine Metabolism',
    'Lysine Degradation':                     'Lysine Metabolism',

    'Ubiquinone Synthesis':                   'Ubiquinone And Other Terpenoid-quinone Biosynthesis',
    'Sulfur Cycle':                           'Sulfur Metabolism',

    'Mitochondrial':                          'Mitochondrial Transport',
    'Peroxisomal':                            'Peroxisomal Transport',
    'Endoplasmic Reticular':                  'Endoplasmic Reticular Transport',
    'Lysosomal':                              'Lysosomal Transport',
    'Nuclear':                                'Nuclear Transport',
    'Golgi Apparatus':                        'Golgi Transport',
    'Endosomal':                              'Endosomal Transport',

    'Selenoamino Acid Metabolism':            'Selenium Metabolism',
    'Selenocompound Metabolism':              'Selenium Metabolism',
}

SUBSYSTEMS_TO_REMOVE = {
    'Metabolic Pathways',
    'Biosynthesis Of Secondary Metabolites',
    'Microbial Metabolism In Diverse Environments',
    'Miscellaneous',
    'Unassigned',
    'Exchange/demand Reaction',
    'Extracellular exchange',
    'Transport',
    'Streptomycin Biosynthesis',
    'Neomycin, Kanamycin And Gentamicin Biosynthesis',
    'Carbapenem Biosynthesis',
    'Monobactam Biosynthesis',
    'Mycolic Acid Biosynthesis',
    'Novobiocin Biosynthesis',
    'Teichoic Acid Biosynthesis',
    'Aflatoxin Biosynthesis',
    'Betalain Biosynthesis',
    'Glucosinolate Biosynthesis',
    'Tropane, Piperidine And Pyridine Alkaloid Biosynthesis',
    'Isoquinoline Alkaloid Biosynthesis',
    'Biosynthesis Of Various Plant Secondary Metabolites',
    'Styrene Degradation',
}


def normalize_subsystems(model) -> dict:
    stats = {'fixed_empty': 0, 'normalized': 0, 'removed': 0, 'changed': 0, 'unchanged': 0}

    for rxn in model.reactions:
        original = rxn.subsystem
        changed = False

        if not rxn.subsystem or rxn.subsystem == 'Unassigned':
            if rxn.id in REACTION_SUBSYSTEM_FIXES:
                rxn.subsystem = REACTION_SUBSYSTEM_FIXES[rxn.id]
                stats['fixed_empty'] += 1
                print(f"Fixed empty: {rxn.id} -> {rxn.subsystem}")
            else:
                rxn.subsystem = ''
                stats['removed'] += 1
            changed = True

        if rxn.subsystem:
            parts = [s.strip() for s in rxn.subsystem.split(';')]
            normalized = []
            for sub in parts:
                if sub in SUBSYSTEMS_TO_REMOVE:
                    stats['removed'] += 1
                    changed = True
                elif sub in SUBSYSTEM_NORMALIZATION:
                    normalized.append(SUBSYSTEM_NORMALIZATION[sub])
                    stats['normalized'] += 1
                    changed = True
                elif sub:
                    normalized.append(sub)

            new_subsystem = '; '.join(dict.fromkeys(normalized))
            if new_subsystem != rxn.subsystem:
                if changed and original not in ('', 'Unassigned'):
                    print(f"  {rxn.id}: [{original}] -> [{new_subsystem}]")
                rxn.subsystem = new_subsystem
                stats['changed'] += 1

        if not changed:
            stats['unchanged'] += 1

    return stats


def collect_subsystems(model) -> set:
    result = set()
    for rxn in model.reactions:
        if rxn.subsystem and rxn.subsystem not in ('Unassigned', ''):
            result.update(s.strip() for s in rxn.subsystem.split(';'))
    return result


def main():
    model = cobra.io.load_json_model("iES1300_prepared.json")

    subs_before = collect_subsystems(model)
    stats = normalize_subsystems(model)
    subs_after = collect_subsystems(model)

    added = subs_after - subs_before
    removed = subs_before - subs_after

    if added:
        print("\nNew subsystems:")
        for s in sorted(added):
            print(f"  + {s}")

    if removed:
        print("\nRemoved subsystems:")
        for s in sorted(removed):
            print(f"  - {s}")

    print(f"\nSubsystems: {len(subs_before)} -> {len(subs_after)}")
    print(f"Stats: {stats}")

    print("\nAll subsystems after normalization:")
    for i, sub in enumerate(sorted(subs_after), 1):
        print(f"{i:3d}. {sub}")

    cobra.io.save_json_model(model, "iES1300_normalized.json")


if __name__ == "__main__":
    main()
