import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.stats import fisher_exact
import sys

# === PLOT PARAMETER DICTIONARIES FOR EASY MODIFICATION ===
DOT_PLOT_PARAMS = {
    'figsize': (6, 6),
    'scatter_size_scale': 80,
    'colormap': 'Reds',
    'alpha': 0.8,
    'edgecolors': 'black',
    'linewidth': 0.5,
    'pvalue_threshold': 0.05,
    'line_color': 'gray',
    'line_style': '--',
    'line_alpha': 0.7,
    'xlabel': '-log10(p-value)',
    'ylabel': 'Phylum',
    'title': '',
    'xlabel_fontsize': 9,
    'xlabel_fontweight': 'bold',
    'ylabel_fontsize': 9,
    'ylabel_fontweight': 'bold',
    'title_fontsize': 9,
    'title_fontweight': 'bold',
    'title_pad': 20,
    'colorbar_shrink': 0.8,
    'colorbar_label': '% Multi-gene Clusters',
    'colorbar_rotation': 270,
    'colorbar_labelpad': 20,
    'colorbar_fontsize': 12,
    'legend_or_values': [1, 2, 4],
    'legend_title': 'Odds Ratio & Significance',
    'legend_loc': 'lower right',
    'legend_frameon': True,
    'legend_fontsize': 10,
    'legend_titlefontsize': 11,
    'grid_alpha': 0.3,
    'dpi': 300
}

DIVERGING_BAR_PARAMS = {
    'figsize': (3, 3),
    'colors': {
        'Significantly Enriched': '#d62728',    # Red
        'Significantly Depleted': '#1f77b4',   # Blue  
        'Not Significant': 'gray'
    },
    'alpha': 0.8,
    'edgecolor': 'black',
    'linewidth': 0.5,
    'reference_line_color': 'black',
    'reference_line_style': '-',
    'reference_line_width': 1,
    'reference_line_alpha': 0.8,
    'xlabel': 'log2(Odds Ratio)',
    'ylabel': 'Phylum',
    'title': '',
    'xlabel_fontsize': 9,
    'xlabel_fontweight': 'bold',
    'ylabel_fontsize': 9,
    'ylabel_fontweight': 'bold',
    'title_fontsize': 9,
    'title_fontweight': 'bold',
    'title_pad': 20,
    'pvalue_significance_cutoff': 0.05,
    'text_offset': 0.15,
    'text_fontsize': 8,
    'annotation_fontsize': 10,
    'grid_alpha': 0.3,
    'legend_loc': 'lower right',
    'legend_fontsize': 8,
    'dpi': 300
}

def create_dot_plot(plot_df, output_dot_svg, params):
    """Create dot plot with configurable parameters"""
    fig, ax = plt.subplots(figsize=params['figsize'])
    
    scatter = ax.scatter(plot_df['neg_log_p'], 
                        plot_df['Phylum'],
                        s=plot_df['odds_ratio'] * params['scatter_size_scale'],  
                        c=plot_df['percent_multi_gene_clusters'], 
                        cmap=params['colormap'], 
                        alpha=params['alpha'],
                        edgecolors=params['edgecolors'],
                        linewidth=params['linewidth'])
    
    ax.axvline(-np.log10(params['pvalue_threshold']), color=params['line_color'], 
               linestyle=params['line_style'], alpha=params['line_alpha'], 
               label=f'p = {params["pvalue_threshold"]}')
    
    ax.set_xlabel(params['xlabel'], fontsize=params['xlabel_fontsize'], 
                  fontweight=params['xlabel_fontweight'])
    ax.set_ylabel(params['ylabel'], fontsize=params['ylabel_fontsize'], 
                  fontweight=params['ylabel_fontweight'])
    ax.set_title(params['title'], fontsize=params['title_fontsize'], 
                 fontweight=params['title_fontweight'], pad=params['title_pad'])
    
    cbar = fig.colorbar(scatter, ax=ax, shrink=params['colorbar_shrink'])
    cbar.set_label(params['colorbar_label'], rotation=params['colorbar_rotation'], 
                   labelpad=params['colorbar_labelpad'], fontsize=params['colorbar_fontsize'])
    
    legend_elements = []
    for or_val in params['legend_or_values']:
        legend_elements.append(plt.scatter([], [], s=or_val*params['scatter_size_scale'], 
                                         c='lightgray', edgecolors=params['edgecolors'], 
                                         alpha=params['alpha'], label=f'OR = {or_val}'))
    
    legend_elements.append(plt.Line2D([0], [0], color=params['line_color'], 
                                    linestyle=params['line_style'], alpha=params['line_alpha'], 
                                    label=f'p = {params["pvalue_threshold"]}'))
    
    ax.legend(handles=legend_elements, title=params['legend_title'], 
             loc=params['legend_loc'], frameon=params['legend_frameon'], 
             fontsize=params['legend_fontsize'], title_fontsize=params['legend_titlefontsize'])
    
    ax.grid(True, alpha=params['grid_alpha'])
    ax.set_axisbelow(True)
    
    plt.tight_layout()
    plt.savefig(output_dot_svg, format='svg', dpi=params['dpi'], bbox_inches='tight')
    plt.close()

def create_diverging_bar_plot(bar_df, output_divbar_svg, params):
    """Create diverging bar plot with configurable parameters"""
    fig2, ax2 = plt.subplots(figsize=params['figsize'])
    
    colors = params['colors']
    
    # Create horizontal bars
    bars = ax2.barh(bar_df['Phylum'], bar_df['log2_odds_ratio'], 
                   color=[colors[cat] for cat in bar_df['significance_cat']],
                   alpha=params['alpha'], edgecolor=params['edgecolor'], 
                   linewidth=params['linewidth'])
    
    # Add only the main reference line at 0 (removed blue and red threshold lines)
    ax2.axvline(0, color=params['reference_line_color'], linestyle=params['reference_line_style'], 
                linewidth=params['reference_line_width'], alpha=params['reference_line_alpha'])
    
    # Customize axes
    ax2.set_xlabel(params['xlabel'], fontsize=params['xlabel_fontsize'], 
                   fontweight=params['xlabel_fontweight'])
    ax2.set_ylabel(params['ylabel'], fontsize=params['ylabel_fontsize'], 
                   fontweight=params['ylabel_fontweight'])
    ax2.set_title(params['title'], fontsize=params['title_fontsize'], 
                  fontweight=params['title_fontweight'], pad=params['title_pad'])
    
    # Add p-value annotations for significant results
    for i, (_, row) in enumerate(bar_df.iterrows()):
        if row['p_value'] < params['pvalue_significance_cutoff']:
            x_pos = row['log2_odds_ratio'] / 2
            alignment = 'center'
            
            ax2.text(x_pos, i, f"p = {row['p_value']:.1e}", 
                    va='center', ha=alignment, fontsize=params['text_fontsize'], fontweight='bold')
    
    # Create legend (removed threshold line references)
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor=colors['Significantly Enriched'], 
                     edgecolor=params['edgecolor'], label='Enriched'),
        plt.Rectangle((0,0),1,1, facecolor=colors['Significantly Depleted'], 
                     edgecolor=params['edgecolor'], label='Depleted'),
        plt.Rectangle((0,0),1,1, facecolor=colors['Not Significant'], 
                     edgecolor=params['edgecolor'], label='NS')
    ]
    
    ax2.legend(handles=legend_elements, loc=params['legend_loc'], frameon=True, 
               fontsize=params['legend_fontsize'])
    
    # Add biological interpretation annotations
    
    ax2.grid(True, alpha=params['grid_alpha'], axis='x')
    ax2.set_axisbelow(True)
    
    plt.tight_layout()
    plt.savefig(output_divbar_svg, format='svg', dpi=params['dpi'], bbox_inches='tight')
    plt.close()

def analyze_enrichment_with_dual_plots(file_path, output_tsv, output_dot_svg, output_divbar_svg):
    """
    Analyzes phylum enrichment and creates TSV results, dot plot, and diverging bar plot.
    """
    try:
        print("Step 1: Reading the combined data file...")
        df = pd.read_csv(file_path, sep=r'\s+', header=0, dtype=str, engine='python')
        
        if len(df.columns) != 5:
            print(f"Error: Expected 5 columns, found {len(df.columns)}")
            return
            
        df.columns = ["cluster", "genome", "gene", "First_Field", "Phylum"]
        print(f"Successfully loaded {len(df)} rows")
        
        print("Step 2: Counting genes per cluster and assigning phylum...")
        gene_counts = df.groupby('cluster')['gene'].nunique().rename('gene_count')
        phylum_mode = df.groupby('cluster')['Phylum'].apply(lambda x: x.mode()[0] if not x.mode().empty else "Unknown")
        cluster_summary_df = pd.concat([gene_counts, phylum_mode], axis=1).reset_index()
        
        print("Step 3: Categorizing clusters as single-gene or multi-gene...")
        cluster_summary_df['cluster_type'] = cluster_summary_df['gene_count'].apply(lambda x: 'multi-gene' if int(x) > 1 else 'single-gene')
        total_clusters = len(cluster_summary_df)
        multi_gene_clusters = len(cluster_summary_df[cluster_summary_df['cluster_type'] == 'multi-gene'])
        print(f"Total clusters: {total_clusters}")
        print(f"Multi-gene clusters: {multi_gene_clusters} ({multi_gene_clusters/total_clusters*100:.1f}%)")
        
        print("Step 4: Performing Fisher's Exact Test for each phylum...")
        unique_phyla = cluster_summary_df['Phylum'].unique()
        results = []
        total_multi_gene = len(cluster_summary_df[cluster_summary_df['cluster_type'] == 'multi-gene'])
        total_single_gene = len(cluster_summary_df[cluster_summary_df['cluster_type'] == 'single-gene'])
        
        for phylum in unique_phyla:
            if phylum == "Unknown":
                continue
            
            in_phylum = cluster_summary_df[cluster_summary_df['Phylum'] == phylum]
            a = len(in_phylum[in_phylum['cluster_type'] == 'multi-gene'])
            c = len(in_phylum[in_phylum['cluster_type'] == 'single-gene'])
            b = total_multi_gene - a
            d = total_single_gene - c
            
            if a + c == 0:  # Skip phyla with no clusters
                continue
                
            table = [[a, b], [c, d]]
            odds_ratio, p_value = fisher_exact(table)
            
            total_in_phylum = a + c
            percent_multi = (a / total_in_phylum * 100) if total_in_phylum > 0 else 0
            results.append({
                'Phylum': phylum,
                'p_value': p_value,
                'odds_ratio': odds_ratio,
                'percent_multi_gene_clusters': round(percent_multi, 2),
                'multi_gene_clusters_in_phylum': a,
                'single_gene_clusters_in_phylum': c,
                'total_clusters_in_phylum': total_in_phylum
            })
        
        print("Step 5: Saving TSV results...")
        results_df = pd.DataFrame(results).sort_values(by='p_value', ascending=True)
        
        # Save TSV results
        with open(output_tsv, 'w') as f:
            f.write("Phylum\tp_value\todds_ratio\tpercent_multi_gene_clusters\tmulti_gene_clusters_in_phylum\tsingle_gene_clusters_in_phylum\ttotal_clusters_in_phylum\n")
            for _, row in results_df.iterrows():
                f.write(f"{row['Phylum']}\t{row['p_value']}\t{row['odds_ratio']}\t{row['percent_multi_gene_clusters']}\t{row['multi_gene_clusters_in_phylum']}\t{row['single_gene_clusters_in_phylum']}\t{row['total_clusters_in_phylum']}\n")
        
        print("Step 6: Creating dot plot visualization...")
        
        # === DOT PLOT ===
        plot_df = results_df.copy()
        plot_df['p_value_plot'] = plot_df['p_value'].replace(0, 1e-300)
        plot_df['neg_log_p'] = -np.log10(plot_df['p_value_plot'])
        
        # Adjust figure height based on number of phyla
        dot_params = DOT_PLOT_PARAMS.copy()
        dot_params['figsize'] = (12, max(6, len(plot_df) * 0.5))
        
        create_dot_plot(plot_df, output_dot_svg, dot_params)
        
        print("Step 7: Creating diverging bar plot visualization...")
        
        # === DIVERGING BAR PLOT ===
        bar_df = results_df[results_df['p_value'] <= 0.1].copy()  # Filter for relevance
        
        # Handle edge cases in odds_ratio calculation
        bar_df['odds_ratio_safe'] = bar_df['odds_ratio'].replace([0, np.inf], [1e-10, 1e10])
        bar_df['log2_odds_ratio'] = np.log2(bar_df['odds_ratio_safe'])
        bar_df['log2_odds_ratio'] = bar_df['log2_odds_ratio'].clip(lower=-10, upper=10)
        bar_df = bar_df.sort_values('log2_odds_ratio', ascending=True)
        
        # Create significance categories
        bar_df['significance_cat'] = 'Not Significant'
        bar_df.loc[(bar_df['p_value'] < 0.05) & (bar_df['log2_odds_ratio'] > 0), 'significance_cat'] = 'Significantly Enriched'
        bar_df.loc[(bar_df['p_value'] < 0.05) & (bar_df['log2_odds_ratio'] < 0), 'significance_cat'] = 'Significantly Depleted'
        
        # Adjust figure height based on number of phyla
        bar_params = DIVERGING_BAR_PARAMS.copy()
        bar_params['figsize'] = (5, max(3, len(bar_df) * 0.6))
        
        create_diverging_bar_plot(bar_df, output_divbar_svg, bar_params)
        
        print(f"\nAnalysis complete!")
        print(f"TSV results saved to: {output_tsv}")
        print(f"Dot plot saved to: {output_dot_svg}")
        print(f"Diverging bar plot saved to: {output_divbar_svg}")
        
        print("\n--- Top Enriched Phyla for Multi-Gene Clusters ---")
        print(f"{'Phylum':<20} {'p-value':<12} {'Odds Ratio':<12} {'% Multi-gene':<12}")
        print("-" * 60)
        for _, row in results_df.head().iterrows():
            print(f"{row['Phylum']:<20} {row['p_value']:<12.2e} {row['odds_ratio']:<12.2f} {row['percent_multi_gene_clusters']:<12.1f}%")
            
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze phylum enrichment with dot plot and diverging bar plot.")
    parser.add_argument("combined_file", help="Path to the combined input file.")
    parser.add_argument("output_tsv", help="Path for the TSV results file.")
    parser.add_argument("output_dot_svg", help="Path for the dot plot SVG file.")
    parser.add_argument("output_divbar_svg", help="Path for the diverging bar plot SVG file.")
    
    args = parser.parse_args()
    
    try:
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend for servers
        from scipy.stats import fisher_exact
    except ImportError as e:
        print(f"Error: Required library not found. Please install: pip install scipy matplotlib")
        sys.exit(1)
        
    analyze_enrichment_with_dual_plots(args.combined_file, args.output_tsv, args.output_dot_svg, args.output_divbar_svg)
