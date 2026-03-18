#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import yaml

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.titlesize'] = 12

with open('input_parameters.yaml', 'r') as f:
    config = yaml.safe_load(f)

gsea_dir = Path("results/transcriptomics/gsea")
output_dir = Path("results/transcriptomics/gsea/plots")
output_dir.mkdir(parents=True, exist_ok=True)
tissues = [t['tissue'] for t in config['tissues']]

def load_gsea_results(tissue):
    """Load GSEA results for a tissue"""
    file_path = gsea_dir / f"{tissue}_gsea_results_stat.csv"
    if not file_path.exists():
        print(f"Warning: {file_path} not found")
        return None

    df = pd.read_csv(file_path)
    df['tissue'] = tissue
    return df

def prepare_nes_matrix(results_dict, padj_cutoff=0.25, qvalue_cutoff=0.25, min_tissues=1):
    """
    Prepare NES matrix for heatmap
    """
    all_pathways = {}

    for tissue, df in results_dict.items():
        if df is None:
            continue
        sig_df = df[(df['p.adjust'] < padj_cutoff) & (df['qvalue'] < qvalue_cutoff)].copy()

        for _, row in sig_df.iterrows():
            pathway_id = row['ID']
            pathway_desc = row['Description']

            if pathway_id not in all_pathways:
                all_pathways[pathway_id] = {
                    'description': pathway_desc,
                    'nes': {},
                    'padj': {},
                    'count': 0
                }

            all_pathways[pathway_id]['nes'][tissue] = row['NES']
            all_pathways[pathway_id]['padj'][tissue] = row['p.adjust']
            all_pathways[pathway_id]['count'] += 1

    filtered_pathways = {
        pid: data for pid, data in all_pathways.items()
        if data['count'] >= min_tissues
    }

    if not filtered_pathways:
        print(f"No pathways found in at least {min_tissues} tissues with padj < {padj_cutoff}")
        return None

    nes_matrix = pd.DataFrame(index=filtered_pathways.keys(), columns=tissues)

    for pathway_id, data in filtered_pathways.items():
        for tissue in tissues:
            nes_matrix.loc[pathway_id, tissue] = data['nes'].get(tissue, 0)

    nes_matrix = nes_matrix.astype(float)

    nes_matrix.index.name = 'GO_ID'
    descriptions = pd.Series(
        {pid: data['description'] for pid, data in filtered_pathways.items()},
        name='Description'
    )

    return nes_matrix, descriptions

def plot_nes_heatmap(nes_matrix, descriptions, output_file, figsize=(12, 16)):
    """
    Create NES heatmap
    """
    plot_data = nes_matrix
    labels = descriptions.loc[plot_data.index].apply(
        lambda x: x[:60] + '...' if len(x) > 60 else x
    )
    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        plot_data,
        cmap='RdBu_r',
        center=0,
        vmin=-3,
        vmax=3,
        yticklabels=labels,
        xticklabels=[t.capitalize() for t in plot_data.columns],
        cbar_kws={
            'label': 'NES',
            'orientation': 'horizontal',
            'shrink': 0.6,
            'aspect': 30,
            'pad': 0.15  
        },
        linewidths=0.5,
        linecolor='lightgray',
        ax=ax
    )
    ax.set_title('GSEA Enrichment Across Tissues', fontweight='bold', pad=10)
    ax.set_xlabel('Tissue', fontweight='bold')
    ax.set_ylabel('GO Biological Process', fontweight='bold')

    plt.tight_layout()
    plt.subplots_adjust(top=0.97, bottom=0.12) 
    plt.savefig(output_file.with_suffix('.tiff'), dpi=300, bbox_inches='tight',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.savefig(output_file.with_suffix('.eps'), format='eps', bbox_inches='tight')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file.with_suffix('.tiff')} (PLOS format)")

    plt.close()

def plot_clustered_heatmap(nes_matrix, descriptions, output_file, figsize=(14, 16)):
    """
    Create clustered NES heatmap with hierarchical clustering

    """
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import pdist

    plot_data = nes_matrix

    labels = descriptions.loc[plot_data.index].apply(
        lambda x: x[:60] + '...' if len(x) > 60 else x
    )

    row_linkage = linkage(pdist(plot_data, metric='euclidean'), method='ward')

    g = sns.clustermap(
        plot_data,
        cmap='RdBu_r',
        center=0,
        vmin=-3,
        vmax=3,
        row_linkage=row_linkage,
        col_cluster=False,  
        yticklabels=labels,
        xticklabels=[t.capitalize() for t in plot_data.columns],
        cbar_kws={
            'label': 'NES',
            'orientation': 'horizontal'
        },
        cbar_pos=(0.3, -0.05, 0.4, 0.02),
        linewidths=0.5,
        linecolor='lightgray',
        figsize=figsize,
        dendrogram_ratio=0.1
    )


    g.ax_heatmap.set_xlabel('Tissue', fontweight='bold')
    g.ax_heatmap.set_ylabel('GO Biological Process', fontweight='bold')


    plt.savefig(output_file.with_suffix('.tiff'), dpi=300, bbox_inches='tight',
                pil_kwargs={'compression': 'tiff_lzw'})

    plt.savefig(output_file.with_suffix('.eps'), format='eps', bbox_inches='tight')

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file.with_suffix('.tiff')} (PLOS format)")

    plt.close()

def prepare_nes_matrix_unfiltered(results_dict, min_abs_nes=0.5, min_tissues=2):
    """
    Prepare NES matrix WITHOUT filtering by padj/qvalue
    Instead, show all pathways with sufficient NES values

    """
    all_pathways = {}

    for tissue, df in results_dict.items():
        if df is None:
            continue

        for _, row in df.iterrows():
            pathway_id = row['ID']
            pathway_desc = row['Description']
            nes = row['NES']
            qval = row['qvalue']

            if abs(nes) > min_abs_nes:  
                if pathway_id not in all_pathways:
                    all_pathways[pathway_id] = {
                        'description': pathway_desc,
                        'nes': {},
                        'qvalue': {},
                        'count': 0
                    }

                all_pathways[pathway_id]['nes'][tissue] = nes
                all_pathways[pathway_id]['qvalue'][tissue] = qval
                all_pathways[pathway_id]['count'] += 1

    filtered_pathways = {
        pid: data for pid, data in all_pathways.items()
        if data['count'] >= min_tissues
    }

    if not filtered_pathways:
        print(f"No pathways found in at least {min_tissues} tissues with |NES| > {min_abs_nes}")
        return None
    tissues = list(results_dict.keys())
    nes_matrix = pd.DataFrame(index=filtered_pathways.keys(), columns=tissues, dtype=float)
    qvalue_matrix = pd.DataFrame(index=filtered_pathways.keys(), columns=tissues, dtype=float)

    for pathway_id, data in filtered_pathways.items():
        for tissue in tissues:
            nes_matrix.loc[pathway_id, tissue] = data['nes'].get(tissue, 0)
            qvalue_matrix.loc[pathway_id, tissue] = data['qvalue'].get(tissue, 1.0)

    descriptions = pd.Series(
        {pid: data['description'] for pid, data in filtered_pathways.items()},
        name='Description'
    )

    return nes_matrix, qvalue_matrix, descriptions

def create_significance_annotations(qvalue_matrix):
    """
    Create significance annotations based on qvalue
    """
    annot = qvalue_matrix.copy()

    # Convert qvalues to significance stars
    annot = annot.applymap(lambda x:
        '***' if x < 0.01 else
        '**' if x < 0.05 else
        '*' if x < 0.1 else
        ''
    )

    return annot

def plot_annotated_heatmap(nes_matrix, qvalue_matrix, descriptions, output_file, figsize=(7.5, 8.5)):
    """
    Create heatmap with significance annotations inside cells
    """
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import pdist

    labels = descriptions.loc[nes_matrix.index].apply(
        lambda x: x[:60] + '...' if len(x) > 60 else x
    )
    annot = create_significance_annotations(qvalue_matrix)
    row_linkage = linkage(pdist(nes_matrix, metric='euclidean'), method='ward')

    g = sns.clustermap(
        nes_matrix,
        cmap='RdBu_r',
        center=0,
        vmin=-3,
        vmax=3,
        row_linkage=row_linkage,
        col_cluster=False,
        yticklabels=labels,
        xticklabels=[t.capitalize() for t in nes_matrix.columns],
        annot=annot,
        fmt='s',
        annot_kws={'fontsize': 8, 'fontweight': 'bold'},
        cbar_kws={
            'label': 'NES',
            'orientation': 'horizontal'
        },
        cbar_pos=(0.3, -0.05, 0.4, 0.02),
        linewidths=0.5,
        linecolor='lightgray',
        figsize=figsize,
        dendrogram_ratio=0.1
    )

    g.fig.suptitle('GSEA Enrichment Across Tissues (Annotated)',
                   fontweight='bold', y=0.995)
    g.ax_heatmap.set_xlabel('Tissue', fontweight='bold')
    g.ax_heatmap.set_ylabel('GO Biological Process', fontweight='bold')

    legend_text = '*** q<0.01   ** q<0.05   * q<0.1'
    g.fig.text(0.5, 0.01, legend_text, ha='center', fontsize=9, style='italic')

    plt.savefig(output_file.with_suffix('.tiff'), dpi=300, bbox_inches='tight',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.savefig(output_file.with_suffix('.eps'), format='eps', bbox_inches='tight')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file.with_suffix('.tiff')} (annotated version)")

    plt.close()

def create_summary_stats(results_dict, output_file):
    """Create summary statistics table"""
    summary = []

    for tissue, df in results_dict.items():
        if df is None:
            continue

        sig_005 = (df['p.adjust'] < 0.05).sum()
        sig_025 = (df['p.adjust'] < 0.25).sum()
        total = len(df)

        summary.append({
            'Tissue': tissue.capitalize(),
            'Total pathways': total,
            'Significant (padj<0.05)': sig_005,
            'Significant (padj<0.25)': sig_025,
            'Upregulated (NES>0)': (df['NES'] > 0).sum(),
            'Downregulated (NES<0)': (df['NES'] < 0).sum()
        })

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output_file, index=False)
    print(f"\nSummary statistics:\n{summary_df.to_string(index=False)}\n")
    print(f"Saved: {output_file}")

def main():

    print("Loading GSEA results...")
    results_dict = {}
    for tissue in tissues:
        df = load_gsea_results(tissue)
        if df is not None:
            print(f"  {tissue}: {len(df)} pathways")
            results_dict[tissue] = df

    if not results_dict:
        print("Error: No GSEA results found!")
        return
    print("\nCreating summary statistics...")
    create_summary_stats(results_dict, output_dir / "gsea_summary_stats.csv")

    print("\nPreparing NES matrix...")
    print("  Filtering: p.adjust < 0.25 AND qvalue < 0.25")
    print("  Min tissues: 2")
    result = prepare_nes_matrix(results_dict, padj_cutoff=0.05, qvalue_cutoff=0.05, min_tissues=2)

    if result is None:
        return

    nes_matrix, descriptions = result
    print(f"  Total pathways: {len(nes_matrix)}")
    full_matrix = nes_matrix.copy()
    full_matrix['Description'] = descriptions
    full_matrix.to_csv(output_dir / "nes_matrix_all_pathways.csv")
    print(f"  Saved full NES matrix: {output_dir / 'nes_matrix_all_pathways.csv'}")

    print("\nCreating visualizations...")

    n_pathways = len(nes_matrix)
    width = 7.5
    height = min(8.5, max(4, n_pathways * 0.15))  
    print(f"  Figure size: {width} x {height} inches")
    print(f"  At 300 dpi: {int(width*300)} x {int(height*300)} pixels")

    plot_nes_heatmap(
        nes_matrix,
        descriptions,
        output_dir / "gsea_nes_heatmap_min2tissues.png",
        figsize=(width, height)
    )

    plot_clustered_heatmap(
        nes_matrix,
        descriptions,
        output_dir / "gsea_nes_heatmap_clustered_min2tissues.png",
        figsize=(width, height)
    )
    print("  Using |NES| > 0.5 filter instead of padj/qvalue")
    result_unfiltered = prepare_nes_matrix_unfiltered(results_dict, min_abs_nes=0.5, min_tissues=2)

    if result_unfiltered is not None:
        nes_matrix_unfilt, qvalue_matrix, descriptions_unfilt = result_unfiltered
        print(f"  Total pathways: {len(nes_matrix_unfilt)}")

        n_pathways_unfilt = len(nes_matrix_unfilt)
        height_unfilt = min(8.5, max(4, n_pathways_unfilt * 0.15))

        plot_annotated_heatmap(
            nes_matrix_unfilt,
            qvalue_matrix,
            descriptions_unfilt,
            output_dir / "gsea_nes_heatmap_annotated.png",
            figsize=(width, height_unfilt)
        )
    print(f"Results saved to: {output_dir}")

if __name__ == "__main__":
    main()
