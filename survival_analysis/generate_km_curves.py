#!/usr/bin/env python3
"""
Kaplan-Meier Survival Curve Generator

This script generates Kaplan-Meier survival curves stratified by cell abundance levels.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 12

def load_and_prepare_data(data_dir):
    """Load cell abundance and clinical data."""
    
    data_path = Path(data_dir)
    
    # Load cell abundance data
    print("Loading cell abundance data...")
    abundance_df = pd.read_csv(data_path / 'cell_abundance.csv')
    
    # Load clinical data
    print("Loading clinical data...")
    clinical_df = pd.read_csv(data_path / 'project_labels.csv')
    
    # Get patient-level abundance (pivot to wide format)
    abundance_wide = abundance_df.pivot(
        index='patient_id', 
        columns='cell_type', 
        values='abundance'
    ).reset_index()
    
    # Merge with clinical data
    merged_df = abundance_wide.merge(
        clinical_df[['ID', 'FU_DAYS', 'FU_MONTHS', 'FU_STATUS', 'stage_1_4', 'age_at_diagnosis']].drop_duplicates(),
        left_on='patient_id',
        right_on='ID',
        how='inner'
    )
    
    # Filter out patients with missing survival data
    merged_df = merged_df[merged_df['FU_DAYS'].notna() & merged_df['FU_STATUS'].notna()]
    
    # Convert days to months for easier interpretation
    merged_df['FU_MONTHS'] = merged_df['FU_DAYS'] / 30.44
    
    print(f"\nMerged data (All Stages): {len(merged_df)} patients with complete data")
    print(f"Events (deaths): {merged_df['FU_STATUS'].sum()}")
    print(f"Censored: {len(merged_df) - merged_df['FU_STATUS'].sum()}")
    
    return merged_df

def create_km_curve(data, cell_type, output_dir, median_split=True):
    """Create KM curve for a specific cell type."""
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Check if cell type exists in data
    if cell_type not in data.columns:
        print(f"WARNING: {cell_type} not found in data")
        return None
    
    # Remove patients with missing values for this cell type
    plot_data = data[data[cell_type].notna()].copy()
    
    if len(plot_data) == 0:
        print(f"WARNING: No valid data for {cell_type}")
        return None
    
    # Stratify by median or tertiles
    if median_split:
        threshold = plot_data[cell_type].median()
        plot_data['group'] = (plot_data[cell_type] >= threshold).astype(int)
        group_labels = ['Low', 'High']
    else:
        # Tertile split
        tertiles = plot_data[cell_type].quantile([0.33, 0.67])
        plot_data['group'] = pd.cut(
            plot_data[cell_type], 
            bins=[-np.inf, tertiles.iloc[0], tertiles.iloc[1], np.inf],
            labels=['Low', 'Medium', 'High']
        )
        group_labels = ['Low', 'Medium', 'High']
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Fit and plot KM curves for each group
    kmf = KaplanMeierFitter()
    
    for i, label in enumerate(group_labels):
        mask = plot_data['group'] == (i if median_split else label)
        if mask.sum() == 0:
            continue
            
        group_data = plot_data[mask]
        
        kmf.fit(
            durations=group_data['FU_MONTHS'],
            event_observed=group_data['FU_STATUS'],
            label=f'{label} ({len(group_data)} pts)'
        )
        
        kmf.plot_survival_function(ax=ax, ci_show=True)
    
    # Perform chi-squared test on event status between groups
    if median_split:
        from scipy.stats import chi2_contingency
        
        group_0 = plot_data[plot_data['group'] == 0]
        group_1 = plot_data[plot_data['group'] == 1]
        
        # Create contingency table: [dead, alive] for each group
        contingency_table = [
            [group_0['FU_STATUS'].sum(), len(group_0) - group_0['FU_STATUS'].sum()],
            [group_1['FU_STATUS'].sum(), len(group_1) - group_1['FU_STATUS'].sum()]
        ]
        
        chi2_stat, p_value, dof, expected = chi2_contingency(contingency_table)
        
        # Calculate median survival times
        kmf_low = KaplanMeierFitter()
        kmf_low.fit(group_0['FU_MONTHS'], group_0['FU_STATUS'])
        median_low = kmf_low.median_survival_time_
        
        kmf_high = KaplanMeierFitter()
        kmf_high.fit(group_1['FU_MONTHS'], group_1['FU_STATUS'])
        median_high = kmf_high.median_survival_time_
        
    else:
        # For tertiles, use pairwise comparison
        p_value = None
        median_low = None
        median_high = None
    
    # Format plot
    ax.set_xlabel('Time (months)', fontsize=12)
    ax.set_ylabel('Survival Probability', fontsize=12)
    ax.set_ylim([0, 1.05])
    
    # Clean cell type name for title
    clean_name = cell_type.replace('+', ' ')
    title = f'Kaplan-Meier Curve: {clean_name}'
    if p_value is not None:
        title += f'\n(Chi² p = {p_value:.4f})'
        if median_low is not None and median_high is not None:
            title += f'\nMedian: Low={median_low:.1f}mo, High={median_high:.1f}mo'
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    ax.legend(loc='best', frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    filename = cell_type.replace('+', '_').replace('/', '_') + '_KM.png'
    filepath = output_path / filename
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()
    
    p_str = f"{p_value:.4f}" if p_value is not None else "N/A"
    print(f"✓ Saved: {filename} (p={p_str})")
    
    return {
        'cell_type': cell_type,
        'n_patients': len(plot_data),
        'p_value': p_value,
        'median_abundance': plot_data[cell_type].median(),
        'mean_abundance': plot_data[cell_type].mean(),
        'median_survival_low': median_low if median_split else None,
        'median_survival_high': median_high if median_split else None
    }

def main():
    """Main execution function."""
    
    print("="*60)
    print("Kaplan-Meier Survival Curve Generator")
    print("="*60)
    
    data_dir = '/Volumes/ext_TD/01-project/project-life/results'
    output_dir = '/Volumes/ext_TD/01-project/project-life/results/km_curves'
    
    # Load data
    data = load_and_prepare_data(data_dir)
    
    # Get list of cell types (exclude ID and clinical columns, and Negative)
    exclude_cols = ['patient_id', 'ID', 'FU_DAYS', 'FU_MONTHS', 'FU_STATUS', 'stage_1_4', 'age_at_diagnosis', 'Negative']
    cell_types = [col for col in data.columns if col not in exclude_cols]
    
    print(f"\nGenerating KM curves for {len(cell_types)} cell types...")
    print("="*60)
    
    # Generate KM curves for all cell types
    results = []
    for cell_type in cell_types:
        result = create_km_curve(data, cell_type, output_dir, median_split=True)
        if result:
            results.append(result)
    
    # Create summary table
    if results:
        summary_df = pd.DataFrame(results)
        summary_df = summary_df.sort_values('p_value')
        
        # Save summary
        summary_path = Path(output_dir) / 'km_summary.csv'
        summary_df.to_csv(summary_path, index=False)
        
        print("\n" + "="*60)
        print("Top 10 Most Significant Cell Types:")
        print("="*60)
        print(summary_df.head(10)[['cell_type', 'p_value', 'n_patients']].to_string(index=False))
        
        print(f"\n✓ All KM curves saved to: {output_dir}")
        print(f"✓ Summary saved to: {summary_path}")
    
    print("\n" + "="*60)
    print("✓ KM curve generation complete!")
    print("="*60)

if __name__ == '__main__':
    main()
