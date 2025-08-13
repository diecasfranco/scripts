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
    
    # If there are 8 levels in some rows, handle that case safely
    species_index = df.columns.get_loc('Species')
    if species_index + 1 < len(df.columns):
        additional_col = df.columns[species_index + 1]
        # For rows where the additional column has values, concatenate with Species
        mask = df[additional_col].notna()
        if mask.any():  # Only proceed if there are values in the additional column
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

def create_interactive_plot_with_level_selector(df, sample_cols, taxonomy_levels, top_n=15, color_palette=None, group_small=False, threshold=1.0):
    """Create an interactive plot with buttons to select taxonomy level."""
    print("Creating interactive plot with taxonomy level selection...")
    
    # Create a subplot figure for the plot
    fig = make_subplots(rows=1, cols=1)
    
    # Create buttons for each taxonomy level
    buttons = []
    
    # Store all traces in a dictionary, keyed by taxonomy level
    all_traces = {}
    max_num_traces = 0
    
    # Process each taxonomy level separately
    for level in taxonomy_levels:
        # Create a new relative abundance aggregation specifically for this level
        grouped_df = group_by_taxonomy(df, sample_cols, level, group_small, threshold)
        melted_df = melt_for_plotting(grouped_df, sample_cols, level)
        
        # Get top taxa for this level by overall abundance
        top_taxa = melted_df.groupby(level)['Relative_Abundance'].sum().nlargest(top_n).index.tolist()
        
        # Filter to just the top taxa
        plot_df = melted_df[melted_df[level].isin(top_taxa)]
        
        # Get all samples in a consistent order
        all_samples = sorted(plot_df['Sample'].unique())
        
        # Store traces for this level
        level_traces = []
        
        # For each taxon, create a separate bar
        color_idx = 0
        for taxon in top_taxa:
            taxon_df = plot_df[plot_df[level] == taxon]
            
            # Ensure we have data for all samples
            sample_data = pd.DataFrame({'Sample': all_samples})
            sample_data = sample_data.merge(
                taxon_df[['Sample', 'Relative_Abundance']], 
                on='Sample', 
                how='left'
            ).fillna(0)
            
            # Sort by sample name for consistency
            sample_data = sample_data.sort_values('Sample')
            
            # Create the trace
            trace = go.Bar(
                name=taxon,
                x=sample_data['Sample'],
                y=sample_data['Relative_Abundance'],
                marker_color=color_palette[color_idx % len(color_palette)],
                visible=False  # Initially all traces are hidden
            )
            level_traces.append(trace)
            color_idx += 1
            
            # Add the trace to the figure
            fig.add_trace(trace)
        
        # Store the traces for this level
        all_traces[level] = level_traces
        max_num_traces = max(max_num_traces, len(level_traces))
    
    # Make the first level visible by default
    default_level = taxonomy_levels[0]
    for trace in all_traces[default_level]:
        trace.visible = True
    
    # Create buttons for each taxonomy level
    for level in taxonomy_levels:
        # Create a list of booleans for trace visibility
        visibility = [False] * len(fig.data)
        
        # Calculate the index range for traces of this level
        start_idx = 0
        for prev_level in taxonomy_levels:
            if prev_level == level:
                break
            start_idx += len(all_traces[prev_level])
        
        end_idx = start_idx + len(all_traces[level])
        
        # Set visibility for traces of this level
        for i in range(start_idx, end_idx):
            visibility[i] = True
        
        # Create the button
        button = dict(
            args=[{"visible": visibility},
                  {"title": f"Relative Abundance by {level} (Top {len(all_traces[level])})",
                   "legend": {"title": level}}],
            label=level,
            method="update"
        )
        buttons.append(button)
    
    # Add buttons to the figure
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                active=0,
                x=0.5,
                y=1.15,
                xanchor="center",
                yanchor="top",
                buttons=buttons,
                bgcolor='white',
                bordercolor='gray',
                font=dict(size=12)
            )
        ]
    )
    
    # Add button menu title
    fig.update_layout(
        annotations=[
            dict(
                text="Taxonomy Level:",
                x=0.25,  # Position to the left of the buttons
                y=1.15,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=14)
            )
        ]
    )
    
    # Update layout for better readability
    fig.update_layout(
        title=f"Relative Abundance by {default_level} (Top {len(all_traces[default_level])})",
        barmode="stack",
        xaxis_title="Sample",
        yaxis_title="Relative Abundance (%)",
        legend_title=default_level,
        hovermode="closest",
        margin=dict(l=50, r=50, t=120, b=50),
        height=700,
        showlegend=True
    )
    
    # Add hover information
    fig.update_traces(
        hovertemplate="<b>%{fullData.name}</b><br>Sample: %{x}<br>Abundance: %{y:.2f}%<extra></extra>"
    )
    
    return fig

def save_figure(fig, output_path):
    """Save the figure as an interactive HTML file."""
    print(f"Saving figure to {output_path}...")
    fig.write_html(output_path)
    print(f"Successfully saved interactive visualization to {output_path}")

def main():
    parser = argparse.ArgumentParser(description='Create interactive OTU abundance visualization with taxonomy level selection')
    parser.add_argument('--input', '-i', required=True, help='Path to OTU table (TSV format)')
    parser.add_argument('--output', '-o', default='otu_visualization.html', help='Output HTML file path')
    parser.add_argument('--default-level', '-d', default='Family', 
                        choices=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'],
                        help='Default taxonomy level to show initially')
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
    
    # Reorder taxonomy levels to put the default level first
    if args.default_level in taxonomy_levels:
        taxonomy_levels.remove(args.default_level)
        taxonomy_levels.insert(0, args.default_level)
    
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
    
    # Select color palette
    color_palette = getattr(px.colors.qualitative, args.color_palette)
    
    # Create interactive plot with taxonomy level selection
    fig = create_interactive_plot_with_level_selector(
        abundance_df, 
        sample_cols, 
        taxonomy_levels,
        args.top_n, 
        color_palette,
        args.group_small,
        args.threshold
    )
    
    # Save the figure
    save_figure(fig, args.output)

if __name__ == "__main__":
    main()