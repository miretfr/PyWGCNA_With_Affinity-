import PyWGCNA
from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  
import matplotlib.pyplot as plt
from sklearn.cluster import AffinityPropagation
import time

OKGREEN = "\033[92m"
OKBLUE = "\033[94m"
ENDC = "\033[0m"
BOLD = "\033[1m"




def moduleEigengenes(expr, colors):
    """Calculate module eigengenes."""
    # Get unique colors/modules
    unique_colors = np.unique(colors)
    eigengenes = []
    average_expr_list = []

    for color in unique_colors:
        # Get boolean mask for genes in this module
        mask = np.array(colors) == color
        
        if np.sum(mask) == 0:
            # If no genes in this module, create zero vectors
            eigengene = np.zeros(expr.shape[1])
            avg_expr = np.zeros(expr.shape[1])
        elif np.sum(mask) == 1:
            # If only one gene, use its expression as eigengene
            gene_expr = expr.iloc[mask, :].values.flatten()
            eigengene = gene_expr
            avg_expr = gene_expr
        else:
            # Get expression data for genes in this module
            module_expr = expr.iloc[mask, :]
            
            # Calculate average expression
            avg_expr = module_expr.mean(axis=0).values
            
            # Calculate first principal component (eigengene)
            # Note: PCA expects samples as rows, genes as columns
            # So we need to transpose: samples x genes
            pca = PCA(n_components=1)
            eigengene_values = pca.fit_transform(module_expr.T)  # Transpose here
            
            # Flatten to 1D array
            eigengene = eigengene_values.flatten()

        eigengenes.append(eigengene)
        average_expr_list.append(avg_expr)

        # Convert to DataFrames
    eigengenes_df = pd.DataFrame(
        np.array(eigengenes).T,  # Transpose: samples x modules
        columns=[f'ME{color}' for color in unique_colors],
        index=expr.columns  # Sample names
    )
    
    average_expr_df = pd.DataFrame(
        np.array(average_expr_list).T,  # Transpose: samples x modules
        columns=[f'AE{color}' for color in unique_colors],
        index=expr.columns
    )
    MEList = {
        'eigengenes': eigengenes_df,
        'averageExpr': average_expr_df
      }  
    return MEList

def test_affinity_clustering(wgcna_obj):
    """
    Test Affinity Propagation clustering on subset of genes
    """
    print(f"Testing with genes...")
    
    # Get subset of data
    datExpr_subset = wgcna_obj.datExpr[:, :]
    print(f"Using {datExpr_subset.shape[1]} genes and {datExpr_subset.shape[0]} samples")
    
    # Pick soft threshold
    from PyWGCNA import WGCNA
    power, sft = WGCNA.pickSoftThreshold(
        datExpr_subset.to_df(), 
        RsquaredCut=wgcna_obj.RsquaredCut,
        MeanCut=wgcna_obj.MeanCut, 
        powerVector=wgcna_obj.powers,
        networkType=wgcna_obj.networkType
    )
    print(f"Selected power: {power}")
    
    # Calculate adjacency
    adjacency = WGCNA.adjacency(
        datExpr_subset.to_df(), 
        power=power, 
        adjacencyType=wgcna_obj.networkType
    )
    adjacency = pd.DataFrame(
        adjacency,
        columns=datExpr_subset.to_df().columns,
        index=datExpr_subset.to_df().columns
    )
    
    # Calculate TOM
    TOM = WGCNA.TOMsimilarity(adjacency.to_numpy(), TOMType=wgcna_obj.TOMType)
    TOM.columns = datExpr_subset.to_df().columns
    TOM.index = datExpr_subset.to_df().columns
    print(f"TOM shape: {TOM.shape}")
    
    # Affinity Propagation clustering
    print("Clustering with Affinity Propagation...")
    start_time = time.time()

    # Calculate preference from TOM distribution
    triu_indices = np.triu_indices_from(TOM.values, k=1)
    tom_similarities = TOM.values[triu_indices]
    #preference = np.percentile(tom_similarities, 1)  # Adjust 25 to 10, 50, etc.

   # print(f"Using preference={preference} (25th percentile of TOM similarities)")

    
    clust = AffinityPropagation(
        preference=-5,
        damping=0.9,
        affinity='precomputed',
        random_state=42,
        max_iter=200,
        convergence_iter=15,
        verbose=True
    )
    
    labels = clust.fit_predict(TOM.values)
    
    elapsed = time.time() - start_time
    print(f"\n✓ Clustering took {elapsed:.2f} seconds")
    print(f"✓ Found {len(np.unique(labels))} modules")
    print(f"✓ Module sizes: {np.bincount(labels)}")

    # Convert numeric labels into colors
    datExpr_subset.var['moduleLabels'] = labels
    datExpr_subset.var['dynamicColors'] = [f"module_{i}" for i in labels]

    # Calculate module eigengenes (CRITICAL!)
    print("Calculating module eigengenes...")
    from PyWGCNA import WGCNA
    
    # Convert numeric labels to color strings for eigengene calculation
    color_names = ['turquoise', 'blue', 'brown', 'yellow', 'green',
                   'red', 'black', 'pink', 'magenta', 'purple',
                   'greenyellow', 'tan', 'salmon', 'cyan', 'midnightblue',
                   'lightcyan', 'grey60', 'lightgreen', 'lightyellow', 'royalblue']
    
    # Map numeric labels to colors
    unique_labels = np.unique(labels)
    color_map = {i: color_names[i % len(color_names)] for i in range(len(np.unique(labels)))}
    module_colors = [color_map[label] for label in labels]
    datExpr_subset.var['moduleColors'] = module_colors
    
    # Calculate eigengenes
    MEList = moduleEigengenes(
        expr=datExpr_subset.to_df().T, 
        colors=module_colors
    )
    wgcna_obj.MEs = MEList['eigengenes']
    if 'MEgrey' in wgcna_obj.MEs.columns:
        wgcna_obj.MEs.drop(['MEgrey'], axis=1, inplace=True)
    
    # Also save to datME for compatibility
    wgcna_obj.datME = wgcna_obj.MEs.copy()
    # Remove grey module if it exists
    if 'MEgrey' in wgcna_obj.MEs.columns:
        wgcna_obj.MEs = wgcna_obj.MEs.drop(['MEgrey'], axis=1)
        wgcna_obj.datME = wgcna_obj.datME.drop(['MEgrey'], axis=1)
    
    # Also update the original datExpr with the subset
    wgcna_obj.datExpr = datExpr_subset
    
    print(f"✓ Calculated eigengenes for {wgcna_obj.MEs.shape[1]} modules")
    print(f"✓ Eigengenes shape: {wgcna_obj.MEs.shape}")
    
    print("\tDone testing WGCNA..\n")
    
    return labels, TOM

# Now call it:

pyWGCNA_5xFAD = PyWGCNA.readWGCNA("5xFAD.p")

pyWGCNA_5xFAD.findModules = test_affinity_clustering(pyWGCNA_5xFAD )
# After your clustering, check these:
print("Checking WGCNA object state:")
print(f"self.MEs is None: {pyWGCNA_5xFAD.MEs is None}")
print(f"self.datME is None: {pyWGCNA_5xFAD.datME is None}")

# Check if module colors were assigned
if hasattr(pyWGCNA_5xFAD.datExpr, 'var'):
    print(f"Columns in datExpr.var: {pyWGCNA_5xFAD.datExpr.var.columns.tolist()}")
    if 'moduleColors' in pyWGCNA_5xFAD.datExpr.var.columns:
        print(f"Number of genes with module assignments: {len(pyWGCNA_5xFAD.datExpr.var['moduleColors'])}")
        print(f"Unique modules: {pyWGCNA_5xFAD.datExpr.var['moduleColors'].unique()}")
print("done")
# add/update metadata in datExpr and geneExpr anndata
pyWGCNA_5xFAD.updateSampleInfo(path='5xFAD_paper/sampleInfo.csv', sep=',')
print("done")
# set color pallete for each group of metadata
pyWGCNA_5xFAD.setMetadataColor('Sex', {'Female': 'green',
                                       'Male': 'yellow'})
pyWGCNA_5xFAD.setMetadataColor('Genotype', {'5xFADWT': 'darkviolet',
                                            '5xFADHEMI': 'deeppink'})
pyWGCNA_5xFAD.setMetadataColor('Age', {'4mon': 'thistle',
                                       '8mon': 'plum',
                                       '12mon': 'violet',
                                       '18mon': 'purple'})
pyWGCNA_5xFAD.setMetadataColor('Tissue', {'Hippocampus': 'red',
                                          'Cortex': 'blue'})

geneList = PyWGCNA.getGeneList(dataset='mmusculus_gene_ensembl',
                               attributes=['ensembl_gene_id', 
                                           'external_gene_name', 
                                           'gene_biotype'],
                               maps=['gene_id', 'gene_name', 'gene_biotype'])
print("done")

# add/update genes info in datExpr and geneExpr anndata
pyWGCNA_5xFAD.updateGeneInfo(geneList)
print("done")



PyWGCNA.WGCNA.barplotModuleEigenGene = lambda *args, **kwargs: None

# Analysing results: 1.calculating module trait relationship 2.plotting module heatmap eigengene 3.finding GO term for each module
# pyWGCNA_5xFAD.analyseWGCNA()
import os; os.makedirs('figures', exist_ok=True); pyWGCNA_5xFAD.analyseWGCNA()
print("done")
# Saves the current WGCNA in pickle format with the .p extension
pyWGCNA_5xFAD.saveWGCNA()

# find top n hub genes based on connectivity in given module
pyWGCNA_5xFAD.top_n_hub_genes(moduleName="coral", n=10)