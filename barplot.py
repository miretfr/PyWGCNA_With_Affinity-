import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
from matplotlib import gridspec
import PyWGCNA

# Load your saved WGCNA
pyWGCNA_5xFAD = PyWGCNA.readWGCNA("5xFAD.p")

def barplotModuleEigenGene_fixed(self, moduleName, metadata=None, colorBar=None, show=True):
    
    
    #  Get metadata from datExpr.obs
    sampleInfo = self.datExpr.obs
    
    if metadata is None:
        metadata = list(sampleInfo.columns)
    
    # Check if module exists
    modules = self.datExpr.var['moduleColors'].unique()
    if moduleName not in modules:
        print(f"Module '{moduleName}' not found!")
        return None
    
    # Get module eigengene from datME
    ME_df = pd.DataFrame({'eigengeneExp': self.datME["ME" + moduleName].values})
    
    # Build combined dataframe
    df = ME_df.copy(deep=True)
    df['all'] = ''
    for m in metadata:
        df[m] = sampleInfo[m].values  # Use sampleInfo
        df['all'] = df['all'] + '_' + df[m].astype(str)
    df['all'] = df['all'].apply(lambda x: x[1:])
    
    # Create category dataframe
    cat = pd.DataFrame(pd.unique(df['all']), columns=['all'])
    cat[metadata] = cat['all'].str.split('_', expand=True)
    
    # Calculate mean and std
    ybar = df[['all', 'eigengeneExp']].groupby(['all']).mean()['eigengeneExp']
    ebar = df[['all', 'eigengeneExp']].groupby(['all']).std()['eigengeneExp']
    ybar = ybar.loc[cat['all']]
    ebar = ebar.loc[cat['all']]
    
    # Prepare dots
    dot = df[['all', 'eigengeneExp']].copy()
    ind_map = {cat.loc[i, 'all']: i for i in range(cat.shape[0])}
    dot['all'] = dot['all'].map(ind_map)
    
    # Set colorBar
    if colorBar is None:
        colorBar = metadata[-1]
    
    # FIXED: Map colors without converting to int
    palette_values = []
    for val in cat[colorBar]:
        if hasattr(self, 'metadataColors') and colorBar in self.metadataColors and val in self.metadataColors[colorBar]:
            palette_values.append(self.metadataColors[colorBar][val])
        else:
            palette_values.append('lightblue')
    
    # Calculate height ratios
    height_ratios = []
    for m in metadata:
        if hasattr(self, 'metadataColors') and isinstance(self.metadataColors.get(m), dict):
            height_ratios.append(len(self.metadataColors[m].keys()))
        else:
            height_ratios.append(1)
    height_ratios.reverse()
    
    # Create figure
    fig, axs = plt.subplots(nrows=2, ncols=2, 
                           figsize=(max(cat.shape[0] * 0.5, 15), len(metadata) * 4),
                           sharex='col', 
                           gridspec_kw={
                               'height_ratios': [len(metadata) * 0.4, len(metadata) * 0.6],
                               'width_ratios': [cat.shape[0], 3]
                           })
    
    # Create legend subplot
    gs = axs[0, 1].get_gridspec()
    for ax in axs[:, 1]:
        ax.remove()
    ax_legend = fig.add_subplot(gs[:, 1])
    ax_legend.axis('off')
    axs_legend = gridspec.GridSpecFromSubplotSpec(len(metadata), 1,
                                                   subplot_spec=ax_legend.get_subplotspec(),
                                                   height_ratios=height_ratios)
    
    # Plot metadata strips
    ind = list(range(cat.shape[0]))
    for m in metadata:
        handles = []
        x = ind
        y = np.repeat(3000 * metadata.index(m), len(ind))
        color = cat[m].values
        
        if hasattr(self, 'metadataColors') and isinstance(self.metadataColors.get(m), dict):
            # Map colors
            color_mapped = []
            for val in color:
                if val in self.metadataColors[m]:
                    color_mapped.append(self.metadataColors[m][val])
                    if val not in [h.get_label() for h in handles]:
                        patch = mpatches.Patch(color=self.metadataColors[m][val], label=val)
                        handles.append(patch)
                else:
                    color_mapped.append('lightgray')
            
            if m != colorBar:
                axs[0, 0].scatter(x, y, c=color_mapped, s=1600, marker='s')
            
            # Add legend
            ax_leg = plt.Subplot(fig, axs_legend[len(metadata) - 1 - metadata.index(m)])
            ax_leg.legend(title=m, handles=handles, loc='center', fontsize=10)
            ax_leg.axis('off')
            fig.add_subplot(ax_leg)
    
    # Top plot styling
    axs[0, 0].set_title(f"Module Eigengene for {moduleName}", size=20, fontweight="bold")
    axs[0, 0].set_ylim(-2000, np.max(y) + 2000 if len(y) > 0 else 5000)
    axs[0, 0].grid(False)
    axs[0, 0].axis('off')
    
    # Bottom barplot
    axs[1, 0].bar(ind, ybar, align='center', color=palette_values, edgecolor='black', linewidth=0.5)
    axs[1, 0].errorbar(ind, ybar, yerr=ebar, fmt="none", ecolor="red", capsize=3, linewidth=2)
    axs[1, 0].scatter(dot['all'], dot['eigengeneExp'], c='black', alpha=0.4, s=30)
    axs[1, 0].set_xticks(ind)
    axs[1, 0].set_xticklabels(cat['all'].values, rotation=90, ha='right', fontsize=8)
    axs[1, 0].set_ylabel('Module Eigengene Expression', fontsize=12, fontweight='bold')
    axs[1, 0].set_xlabel('Sample Groups', fontsize=12, fontweight='bold')
    axs[1, 0].axhline(y=0, color='gray', linestyle='--', linewidth=1)
    axs[1, 0].set_facecolor('white')
    axs[1, 0].spines['top'].set_visible(False)
    axs[1, 0].spines['right'].set_visible(False)
    
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.25, right=0.85)
    
    # Save
    output_file = f"module_barplot_eigengene_{moduleName}_fixed.pdf"
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f" Saved: {output_file}")
    
    if show:
        plt.show()
    else:
        plt.close(fig)
    
    return fig

# Now create barplots for your modules!
print("\n" + "="*60)
print("CREATING BARPLOTS")
print("="*60 + "\n")
# red","yellow","brown","turquoise","pink","green","greenyellow"
# Your actual modules from the output
for module in ["magenta"]:
    print(f"Creating barplot for {module}...")
    barplotModuleEigenGene_fixed(pyWGCNA_5xFAD, 
                                 moduleName=module, 
                                 metadata=['Age', 'Tissue', 'Sex', 'Genotype'],
                                 colorBar='Genotype',
                                 show=True)

print("\n Done")