#! /usr/bin/env python3
import os
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from dash import Dash, dcc, html
import argparse
from plotly.colors import qualitative

# Function to process a single table
def process_table(file_path, threshold=0.1):
    df = pd.read_csv(file_path, sep='\t')
    df.loc[df['fraction_total_reads'] < threshold, 'name'] = 'Others'
    df_agg = df.groupby('name', as_index=False).agg({
        'fraction_total_reads': 'sum'
    }).sort_values(by='fraction_total_reads', ascending=False)
    return df_agg

# Argument parser
parser = argparse.ArgumentParser(description='Process .bracken files and generate outputs.')
parser.add_argument('tables_dir', type=str, help='Path to the directory containing .bracken files')
parser.add_argument('--sample', type=str, help='Name of the specific sample to plot (without extension)')
parser.add_argument('--sort', type=str, choices=['alphabetical', 'abundance'], default='abundance', help='Sort order for the legend: alphabetical or abundance')
args = parser.parse_args()

# Directory containing the tables
tables_dir = args.tables_dir

# Process all tables and gather data for plotting
all_data = {}
for filename in os.listdir(tables_dir):
    if filename.endswith('.bracken'):
        file_path = os.path.join(tables_dir, filename)
        sample_name = os.path.splitext(filename)[0]
        df_agg = process_table(file_path)
        all_data[sample_name] = df_agg

# Find all unique names
all_names = set()
for data in all_data.values():
    all_names.update(data['name'])

# Sort all names alphabetically or by abundance
if args.sort == 'alphabetical':
    all_names = sorted(all_names)
else:
    name_totals = {name: 0 for name in all_names}
    for data in all_data.values():
        for _, row in data.iterrows():
            name_totals[row['name']] += row['fraction_total_reads']
    all_names = sorted(all_names, key=lambda x: name_totals[x], reverse=True)

# Create a DataFrame with all samples and all taxonomic names, filling missing values with 0
combined_df = pd.DataFrame(index=all_names)
for sample, data in all_data.items():
    combined_df[sample] = combined_df.index.map(lambda name: data.loc[data['name'] == name, 'fraction_total_reads'].sum() if name in data['name'].values else 0)

# Save the combined DataFrame as a TSV file
combined_df.to_csv('combined_data.tsv', sep='\t')

# Define color palette
colors = qualitative.Plotly  # You can choose a different color palette from Plotly's built-in ones or define a custom one

def create_plot(dataframe, sample_names, title, output_file):
    fig = go.Figure()

    for i, name in enumerate(dataframe.index):
        fig.add_trace(go.Bar(
            x=sample_names,
            y=dataframe.loc[name, sample_names],
            name=name,
            marker=dict(color=colors[i % len(colors)])
        ))

    fig.update_layout(
        barmode='stack',
        title=title,
        xaxis_title='Sample',
        yaxis_title='Fraction Total Reads',
        height=600,
    )

    # Save the plot as an HTML file
    pio.write_html(fig, file=output_file, auto_open=False)

# Plot a single sample or all samples
if args.sample:
    if args.sample in combined_df.columns:
        create_plot(combined_df[[args.sample]], [args.sample], f'Proportional Taxonomic Distribution for {args.sample}', f'{args.sample}_taxonomic_distribution.html')
    else:
        print(f"Sample {args.sample} not found.")
else:
    create_plot(combined_df, combined_df.columns, 'Proportional Taxonomic Distribution Across Multiple Samples', 'combined_taxonomic_distribution.html')

    # Save individual plots for each sample
    for sample in combined_df.columns:
        create_plot(combined_df[[sample]], [sample], f'Proportional Taxonomic Distribution for {sample}', f'{sample}_taxonomic_distribution.html')

# Dash Application
app = Dash(__name__)

app.layout = html.Div([
    html.H1("Taxonomic Distribution Stacked Bar Plot"),
    dcc.Graph(figure=fig if args.sample else go.Figure())
])

if __name__ == '__main__':
    app.run_server(debug=True)
