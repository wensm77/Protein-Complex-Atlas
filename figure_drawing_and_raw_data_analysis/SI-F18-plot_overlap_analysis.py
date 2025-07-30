#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

# Set font to DejaVu Sans (Arial-like font available on Linux)
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 12

def plot_overlap_analysis(csv_file, threshold=0.5):
    """
    Plot pie chart and histogram for overlap rate analysis
    """
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Get overlap rates
    overlap_rates = df['overlap_rate']
    
    # Create PDF file
    pdf_filename = csv_file.replace('.csv', '_analysis.pdf')
    
    with PdfPages(pdf_filename) as pdf:
        # Create figure with subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        # Pie Chart
        above_threshold = (overlap_rates >= threshold).sum()
        below_threshold = (overlap_rates < threshold).sum()
        
        labels = [f'≥{threshold}', f'<{threshold}']
        sizes = [above_threshold, below_threshold]
        colors = ['#2E86AB', '#A23B72']
        
        ax1.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors, startangle=90)
        ax1.set_title(f'Overlap Rate Distribution\n(Threshold: {threshold})', fontweight='bold')
        
        # Histogram
        counts, bins, _ = ax2.hist(overlap_rates, bins=50, color='#2E86AB', alpha=1.0, edgecolor='black', linewidth=0.5)
        ax2.axvline(x=threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold = {threshold}')
        ax2.set_xlabel('Overlap Rate', fontweight='bold')
        ax2.set_ylabel('Frequency', fontweight='bold')
        ax2.set_title('Overlap Rate Histogram', fontweight='bold')
        ax2.legend()
        ax2.grid(False)  # Remove grid lines
        
        # Remove spines for cleaner look
        for spine in ax2.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.5)
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight', dpi=300)
        plt.close()
        
        # Create a summary statistics page
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.axis('off')
        
        # Calculate statistics
        stats_text = f"""
        Overlap Rate Analysis Summary
        
        Dataset: {csv_file.split('/')[-1]}
        Total samples: {len(overlap_rates):,}
        Threshold: {threshold}
        
        Statistics:
        • Mean: {overlap_rates.mean():.4f}
        • Median: {overlap_rates.median():.4f}
        • Standard Deviation: {overlap_rates.std():.4f}
        • Minimum: {overlap_rates.min():.4f}
        • Maximum: {overlap_rates.max():.4f}
        
        Threshold Analysis:
        • Samples ≥ {threshold}: {above_threshold:,} ({above_threshold/len(overlap_rates)*100:.1f}%)
        • Samples < {threshold}: {below_threshold:,} ({below_threshold/len(overlap_rates)*100:.1f}%)
        
        Quartiles:
        • Q1 (25th percentile): {overlap_rates.quantile(0.25):.4f}
        • Q2 (50th percentile): {overlap_rates.quantile(0.50):.4f}
        • Q3 (75th percentile): {overlap_rates.quantile(0.75):.4f}
        """
        
        ax.text(0.1, 0.9, stats_text, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', fontfamily='DejaVu Sans')
        
        pdf.savefig(fig, bbox_inches='tight', dpi=300)
        plt.close()
    
    # Export histogram data to CSV
    histogram_data = pd.DataFrame({
        'bin_start': bins[:-1],
        'bin_end': bins[1:],
        'bin_center': (bins[:-1] + bins[1:]) / 2,
        'frequency': counts,
        'cumulative_frequency': np.cumsum(counts),
        'percentage': counts / len(overlap_rates) * 100,
        'cumulative_percentage': np.cumsum(counts) / len(overlap_rates) * 100
    })
    
    histogram_csv_filename = csv_file.replace('.csv', '_histogram_data.csv')
    histogram_data.to_csv(histogram_csv_filename, index=False)
    
    # Export summary statistics to CSV
    summary_stats = pd.DataFrame({
        'metric': [
            'total_samples', 'mean', 'median', 'std', 'min', 'max',
            'q1_25th_percentile', 'q2_50th_percentile', 'q3_75th_percentile',
            'above_threshold_count', 'above_threshold_percentage',
            'below_threshold_count', 'below_threshold_percentage'
        ],
        'value': [
            len(overlap_rates), overlap_rates.mean(), overlap_rates.median(), 
            overlap_rates.std(), overlap_rates.min(), overlap_rates.max(),
            overlap_rates.quantile(0.25), overlap_rates.quantile(0.50), 
            overlap_rates.quantile(0.75),
            above_threshold, above_threshold/len(overlap_rates)*100,
            below_threshold, below_threshold/len(overlap_rates)*100
        ]
    })
    
    summary_csv_filename = csv_file.replace('.csv', '_summary_statistics.csv')
    summary_stats.to_csv(summary_csv_filename, index=False)
    
    print(f"Analysis plots saved to: {pdf_filename}")
    print(f"Histogram data exported to: {histogram_csv_filename}")
    print(f"Summary statistics exported to: {summary_csv_filename}")
    print(f"Total samples: {len(overlap_rates):,}")
    print(f"Samples with overlap rate ≥ {threshold}: {above_threshold:,} ({above_threshold/len(overlap_rates)*100:.1f}%)")
    print(f"Samples with overlap rate < {threshold}: {below_threshold:,} ({below_threshold/len(overlap_rates)*100:.1f}%)")

if __name__ == "__main__":
    csv_file = '/mnt/hdd/yecheng/eval_prot_dimers/virus_human_matches_with_overlap_20250723_183111.csv'
    plot_overlap_analysis(csv_file, threshold=0.5)