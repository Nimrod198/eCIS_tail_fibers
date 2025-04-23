import pandas as pd
import argparse
import plotly.express as px
import os
import colorsys

def parse_tblout(file_path, threshold=2.0):
    """Parse hmmscan tblout file with small hits grouping"""
    df = pd.read_csv(file_path, comment='#', sep='\s+', header=None,
                     usecols=[0, 2], names=['target', 'query'], engine='python')
    
    counts = df['target'].value_counts(normalize=True).mul(100).reset_index()
    counts.columns = ['Target', 'Percentage']
    
    # Group small hits and add sorting priority
    major_hits = counts[counts['Percentage'] >= threshold].copy()
    minor_hits = counts[counts['Percentage'] < threshold].copy()
    
    major_hits['sort_order'] = 1
    minor_hits['sort_order'] = 2
    
    if not minor_hits.empty:
        minor_sum = pd.DataFrame({
            'Target': ['<2%'],
            'Percentage': [minor_hits['Percentage'].sum()],
            'sort_order': [2]
        })
        result = pd.concat([major_hits, minor_sum])
    else:
        result = major_hits
    
    return result.round(2)

def generate_beige_shades(base_hex="#BFA784", n_shades=6):
    """Generate shades around base beige color"""
    base_rgb = tuple(int(base_hex.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
    base_hls = colorsys.rgb_to_hls(*[x/255 for x in base_rgb])
    
    shades = []
    for i in range(n_shades):
        # Maintain hue, adjust lightness/saturation
        h = base_hls[0]
        l = base_hls[1] * (1 - 0.08*i)  # Smaller lightness variation
        s = base_hls[2] * (1 - 0.03*i)  # Smaller saturation variation
        rgb = colorsys.hls_to_rgb(h, l, s)
        shades.append('#%02x%02x%02x' % (
            int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255)
        ))
    return shades

def create_final_plot(file_paths, output_file='final_plot.svg'):
    """Create plot with exact color specifications"""
    # Collect and combine data
    data = []
    for f in file_paths:
        file_name = os.path.splitext(os.path.basename(f))[0]
        df = parse_tblout(f)
        df['File'] = file_name
        data.append(df)
    
    combined_df = pd.concat(data).sort_values(['File', 'sort_order'])
    
    # Create color mapping
    color_map = {
        'Tail_P2_I': '#FFA500',     # Orange
        'Baseplate_J': '#0C8D90',   # Forest green
        'GPW_gp25': '#387654',      # Teal
        '<2%': '#FAF0D7'            # Cream
    }
    
    # Generate beige shades for other targets
    beige_shades = generate_beige_shades(n_shades=12)
    other_targets = [t for t in combined_df['Target'].unique() 
                    if t not in color_map and t != '<2%']
    
    color_map.update({target: beige_shades[i % len(beige_shades)] 
                     for i, target in enumerate(other_targets)})
    
    # Legend ordering
    legend_order = ['Tail_P2_I', 'Baseplate_J', 'GPW_gp25'] + \
                  sorted(other_targets) + \
                  ['<2%']
    
    # Create figure with adjusted dimensions
    fig = px.bar(
        combined_df,
        x='File',
        y='Percentage',
        color='Target',
        color_discrete_map=color_map,
        category_orders={'Target': legend_order},
        title='<b>Pfam Hit Distribution</b>',
        labels={'Percentage': 'Hit Percentage (%)', 'File': 'Sample'},
        hover_data={'Percentage': ':.2f%'},
        barmode='stack',
        width=300,  # Increased base width
        height=600  # Increased base height
    )

    # Formatting with proportional sizing
    fig.update_layout(
        font_family="Arial",
        font_size=12,
        title_font_size=14,
        title_x=0.15,  # Position title left of center
        title_y=0.98,
        legend=dict(
            title=None,
            orientation="v",
            yanchor="middle",
            xanchor="left",
            x=1.02,  # Move legend to right side
            y=0.5,
            itemwidth=40,
            font=dict(size=14)
        ),
        plot_bgcolor='white',
        margin=dict(l=80, r=150, t=80, b=120),  # Adjusted margins
        xaxis=dict(
            tickangle=25,
            tickfont=dict(size=15.5),
            title_font=dict(size=12)
        ),
        yaxis=dict(
            title_font=dict(size=15),
            tickfont=dict(size=13)
        ),
        bargap=0.35,
        uniformtext=dict(
            mode="hide",
            minsize=10  # Minimum text size for percentages
        )
    )

    # Bar text formatting
    fig.for_each_trace(lambda trace: trace.update(
        texttemplate=['%{y:.1f}%' if (trace.name != '<2%' and y > 5) else '' for y in trace.y],
        textposition='inside',
        textfont=dict(
            color='white',
            size=10,
            family='Arial'
        ),
        width=0.8,
        marker_line_width=0
    ))

    # Export with higher scale factor
    fig.write_image(
        output_file,
        format='svg',
        scale=2.0,  # Increased scale for better text clarity
        engine="kaleido"
    )
    print(f"Final plot saved to {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create final Pfam distribution plots')
    parser.add_argument('files', nargs='+', help='hmmscan tblout files')
    parser.add_argument('-o', '--output', default='final_plot.svg', 
                       help='Output SVG filename')
    args = parser.parse_args()
    
    create_final_plot(args.files, args.output)
