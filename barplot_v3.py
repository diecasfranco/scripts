#!/usr/bin/env python3

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import argparse
import numpy as np

def load_data(file_path):
    """Load OTU data from tab-separated file."""
    print(f"Loading data from {file_path}...")
    
    # Read the data, parse strings in numeric columns to numbers
    df = pd.read_csv(file_path, sep='\t')
    
    # Convert numeric columns to numeric type
    # First, identify metadata columns
    metadata_cols = ['OTU_ID', 'OTU_assignment', 'Taxonomy', 'Sequence']
    
    # Then convert all other columns that could be numeric
    for col in df.columns:
        if col not in metadata_cols:
            try:
                df[col] = pd.to_numeric(df[col])
            except (ValueError, TypeError):
                print(f"Warning: Could not convert column '{col}' to numeric type")
    
    return df

def preprocess_data(df):
    """Preprocess the OTU data."""
    print("Preprocessing data...")
    
    # Split taxonomy column into separate levels
    taxonomy_levels = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    df[taxonomy_levels] = df['Taxonomy'].str.split(',', expand=True)
    
    # If there are 8 levels in some rows, adjust accordingly
    expected_columns = set(['OTU_ID', 'OTU_assignment', 'Taxonomy', 'Sequence', 'Total'] + taxonomy_levels)
    if not expected_columns.issubset(df.columns):
        # Validate the existence of the column following 'Species'
        species_index = df.columns.get_loc('Species')
        if species_index + 1 < len(df.columns):
            additional_col = df.columns[species_index + 1]
            # Ensure the column contains relevant data
            if df[additional_col].notna().any():
                pass  # Proceed with the logic
            else:
                additional_col = None
        else:
            additional_col = None
        # Get the column name after Species (if it exists)
        additional_col = df.columns[df.columns.get_loc('Species') + 1]
        # For rows where the additional column has values, concatenate with Species
        mask = df[additional_col].notna()
        df.loc[mask, 'Species'] = df.loc[mask, 'Species'] + '_' + df.loc[mask, additional_col]
    
    # Identify sample columns (exclude metadata columns)
    # Note: 'Total' column should not be considered as a sample
    metadata_cols = ['OTU_ID', 'OTU_assignment', 'Taxonomy', 'Sequence', 'Total'] + taxonomy_levels
    sample_cols = [col for col in df.columns if col not in metadata_cols]
    
    # Double-check that all sample columns are numeric
    for col in sample_cols:
        if not pd.api.types.is_numeric_dtype(df[col]):
            print(f"Converting column '{col}' to numeric type")
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
    
    print(f"Found {len(sample_cols)} sample columns")
    
    return df, sample_cols, taxonomy_levels

def create_relative_abundance(df, sample_cols):
    """Calculate relative abundance for each sample."""
    print("Calculating relative abundance...")
    
    # Create a copy to avoid modifying the original dataframe
    abundance_df = df.copy()
    
    # Calculate relative abundance for each sample
    for sample in sample_cols:
        sample_total = abundance_df[sample].sum()
        if sample_total > 0:  # Avoid division by zero
            abundance_df[sample] = abundance_df[sample] / sample_total * 100
        else:
            print(f"Warning: Sample {sample} has zero total counts")
    
    return abundance_df

def group_by_taxonomy(df, sample_cols, taxonomy_level, group_small=False, threshold=1.0):
    """Group data by the selected taxonomy level."""
    print(f"Grouping by {taxonomy_level}...")
    
    # Group by the selected taxonomy level and sum abundances
    grouped_df = df.groupby(taxonomy_level)[sample_cols].sum().reset_index()
    
    # Sort by total abundance across all samples
    grouped_df['Total_Abundance'] = grouped_df[sample_cols].sum(axis=1)
    grouped_df = grouped_df.sort_values('Total_Abundance', ascending=False)
    
    # Optionally group small taxa into "Other" category
    if group_small:
        # Create a mask for taxa with low abundance
        small_taxa_mask = grouped_df['Total_Abundance'] / len(sample_cols) < threshold
        
        if small_taxa_mask.sum() > 1:  # Only group if there are multiple small taxa
            # Create "Other" row by summing all small taxa
            other_row = grouped_df[small_taxa_mask][sample_cols].sum()
            other_row[taxonomy_level] = "Other"
            other_row['Total_Abundance'] = other_row[sample_cols].sum()
            
            # Remove small taxa and add the "Other" row
            grouped_df = grouped_df[~small_taxa_mask]
            grouped_df = pd.concat([grouped_df, pd.DataFrame([other_row])], ignore_index=True)
            
            # Resort by total abundance
            grouped_df = grouped_df.sort_values('Total_Abundance', ascending=False)
    
    # Drop the helper column
    grouped_df = grouped_df.drop(columns=['Total_Abundance'])
    
    return grouped_df

def melt_for_plotting(df, sample_cols, taxonomy_level):
    """Reshape data for plotting."""
    print("Preparing data for plotting...")
    
    # Melt the dataframe to long format for plotting
    melted_df = pd.melt(
        df,
        id_vars=[taxonomy_level],
        value_vars=sample_cols,
        var_name='Sample',
        value_name='Relative_Abundance'
    )
    
    return melted_df

def create_interactive_plot(melted_df, taxonomy_level, top_n=15, color_palette=None):
    """Create an interactive bar plot with Plotly."""
    print("Creating interactive plot...")
    
    # Get the top N taxa by overall abundance
    top_taxa = melted_df.groupby(taxonomy_level)['Relative_Abundance'].sum().nlargest(top_n).index.tolist()
    
    # Filter to just the top taxa
    plot_df = melted_df[melted_df[taxonomy_level].isin(top_taxa)]
    
    # Check if there are any samples with zero reads
    zero_samples = plot_df.groupby('Sample')['Relative_Abundance'].sum() == 0
    if zero_samples.any():
        zero_sample_names = zero_samples[zero_samples].index.tolist()
        print(f"Warning: The following samples have zero reads: {', '.join(zero_sample_names)}")
    
    # Create interactive bar plot
    fig = px.bar(
        plot_df,
        x='Sample',
        y='Relative_Abundance',
        color=taxonomy_level,
        title=f'Relative Abundance by {taxonomy_level} (Top {len(top_taxa)})',
        labels={'Relative_Abundance': 'Relative Abundance (%)'},
        height=600,
        color_discrete_sequence=px.colors.qualitative.Bold if color_palette is None else color_palette
    )
    
    # Improve layout
    fig.update_layout(
        barmode='stack',
        xaxis_title='Sample',
        yaxis_title='Relative Abundance (%)',
        legend_title=taxonomy_level,
        hovermode='closest',
        # Make the plot more readable
        font=dict(size=12),
        legend=dict(
            font=dict(size=10),
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="right",
            x=1.15
        ),
        margin=dict(l=50, r=150, t=80, b=50)
    )
    
    # Add text labels to show exact percentages when hovering
    fig.update_traces(
        hovertemplate='<b>%{fullData.name}</b><br>Sample: %{x}<br>Abundance: %{y:.2f}%<extra></extra>'
    )
    
    return fig

def save_figure(fig, output_path):
    """Save the figure as an interactive HTML file."""
    print(f"Saving figure to {output_path}...")
    fig.write_html(output_path)
    print(f"Successfully saved visualization to {output_path}")

def main():
    parser = argparse.ArgumentParser(description='Create interactive OTU abundance visualization')
    parser.add_argument('--input', '-i', required=True, help='Path to OTU table (TSV format)')
    parser.add_argument('--output', '-o', default='otu_visualization.html', help='Output HTML file path')
    parser.add_argument('--taxonomy-level', '-t', default='Family', 
                        choices=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'],
                        help='Taxonomy level to visualize')
    parser.add_argument('--top-n', '-n', type=int, default=15, help='Number of top taxa to display')
    parser.add_argument('--group-small', '-g', action='store_true', 
                       help='Group small taxa (less than threshold %% abundance) into "Other" category')
    parser.add_argument('--threshold', type=float, default=1.0,
                       help='Threshold percentage for grouping small taxa (default: 1.0%%)')
    parser.add_argument('--sample-pattern', type=str, default=None,
                       help='Only include samples matching this pattern (e.g. "PM*" for all PM samples)')
    parser.add_argument('--color-palette', type=str, default='Bold',
                       choices=['Bold', 'Pastel', 'Vivid', 'Dark', 'Light', 'Prism'],
                       help='Color palette to use for the plot')
    
    args = parser.parse_args()
    
    # Load and process data
    df = load_data(args.input)
    df, sample_cols, taxonomy_levels = preprocess_data(df)
    
    # Filter samples by pattern if provided
    if args.sample_pattern:
        import fnmatch
        matching_samples = [col for col in sample_cols if fnmatch.fnmatch(col, args.sample_pattern)]
        if not matching_samples:
            print(f"Warning: No samples match the pattern '{args.sample_pattern}'")
            return
        sample_cols = matching_samples
        print(f"Filtered to {len(sample_cols)} samples matching pattern '{args.sample_pattern}'")
    
    # Confirm samples were found
    if not sample_cols:
        raise ValueError("No sample columns identified in the input file. Please check your data format.")
    
    # Calculate relative abundance
    abundance_df = create_relative_abundance(df, sample_cols)
    
    # Group by selected taxonomy level
    grouped_df = group_by_taxonomy(abundance_df, sample_cols, args.taxonomy_level, 
                                  args.group_small, args.threshold)
    
    # Prepare data for plotting
    melted_df = melt_for_plotting(grouped_df, sample_cols, args.taxonomy_level)
    
    # Select color palette
    color_palette = getattr(px.colors.qualitative, args.color_palette)
    
    # Create and save the plot
    fig = create_interactive_plot(melted_df, args.taxonomy_level, args.top_n, color_palette)
    save_figure(fig, args.output)

if __name__ == "__main__":
    main()