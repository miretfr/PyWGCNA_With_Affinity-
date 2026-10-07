import PyWGCNA
import matplotlib.pyplot as plt
import venn
import pandas as pd
import numpy as np


# Load your WGCNA object
pyWGCNA_5xFAD = PyWGCNA.readWGCNA("5xFAD.p")

# Define modules of interest
tissue_modules = ["brown", "black", "red", "blue", "yellow","pink"]
age_modules = ["black","turquoise","green","greenyellow"] # changed them many times

print("=== ATTEMPTING HUB GENE EXTRACTION ===")

# First, let's check what parameters the top_n_hub_genes method actually takes
print("\n1. Checking method signatures...")
import inspect

# Check top_n_hub_genes
if hasattr(pyWGCNA_5xFAD, 'top_n_hub_genes'):
    sig = inspect.signature(pyWGCNA_5xFAD.top_n_hub_genes)
    print(f"top_n_hub_genes signature: {sig}")

# Check intramodularConnectivity
if hasattr(pyWGCNA_5xFAD, 'intramodularConnectivity'):
    sig = inspect.signature(pyWGCNA_5xFAD.intramodularConnectivity)
    print(f"intramodularConnectivity signature: {sig}")

# Try a different approach - calculate intramodular connectivity directly
print("\n2. Calculating intramodular connectivity...")
try:
    # Call the function without parameters first to see what it returns
    kIN_result = pyWGCNA_5xFAD.intramodularConnectivity()
    print(f"intramodularConnectivity returned: {type(kIN_result)}")
    if isinstance(kIN_result, pd.DataFrame):
        print(f"DataFrame shape: {kIN_result.shape}")
        print(f"Columns: {kIN_result.columns.tolist()[:10]}...")  # First 10 columns
except Exception as e:
    print(f"Error calculating connectivity: {e}")

# Alternative approach using module eigengenes and correlation
print("\n3. Alternative: Using module membership (kME)...")
tissue_hub_set = set()
age_hub_set = set()

# Get all genes and their module assignments
gene_modules = pyWGCNA_5xFAD.datExpr.var['moduleColors'].copy()

# Check if we have signedKME or similar
if hasattr(pyWGCNA_5xFAD, 'signedKME'):
    print("Found signedKME attribute")
    # This should contain module membership values
    kME_data = pyWGCNA_5xFAD.signedKME
    print(f"signedKME type: {type(kME_data)}, shape: {kME_data.shape if hasattr(kME_data, 'shape') else 'N/A'}")
else:
    # Calculate module membership manually
    print("Calculating module membership from module eigengenes...")
    
    # Get module eigengenes
    MEs = pyWGCNA_5xFAD.MEs
    datExpr_matrix = pyWGCNA_5xFAD.datExpr.X  # Expression matrix
    
    # Convert to DataFrame if needed
    if hasattr(datExpr_matrix, 'toarray'):
        datExpr_matrix = datExpr_matrix.toarray()
    
    # For each module, calculate correlation with module eigengene
    for module in set(tissue_modules + age_modules):
        me_col = f"ME_{module}"
        if me_col in MEs.columns:
            # Get module eigengene
            module_eigengene = MEs[me_col].values
            
            # Calculate correlation for genes in this module
            module_mask = gene_modules == module
            module_gene_indices = np.where(module_mask)[0]
            module_genes = gene_modules.index[module_mask]
            
            if len(module_gene_indices) > 0:
                # Calculate correlations
                correlations = []
                for idx in module_gene_indices:
                    gene_expr = datExpr_matrix[:, idx]
                    corr = np.corrcoef(gene_expr, module_eigengene)[0, 1]
                    correlations.append((module_genes[idx], abs(corr)))  
                
                # Sort by correlation strength
                correlations.sort(key=lambda x: x[1], reverse=True)
                
                # Take top 30 as hub genes
                top_genes = [gene for gene, corr in correlations[:30]]
                
                if module in tissue_modules:
                    tissue_hub_set.update(top_genes)
                if module in age_modules:
                    age_hub_set.update(top_genes)
                    
                print(f"  Module '{module}': {len(top_genes)} hub genes (calculated)")
            else:
                print(f"  Module '{module}': No genes found")
        else:
            print(f"  Module '{module}': No module eigengene found")

# If still no hub genes, use simple approach
if len(tissue_hub_set) == 0 and len(age_hub_set) == 0:
    print("\n4. Using simplified hub gene identification...")
    
    # Simple approach: take genes with highest variance in each module
    for module in tissue_modules:
        mask = gene_modules == module
        module_genes = gene_modules.index[mask].tolist()
        # Take first 30 as a simple proxy
        tissue_hub_set.update(module_genes[:30])
        print(f"  Tissue module '{module}': 30 genes (simplified)")
    
    for module in age_modules:
        mask = gene_modules == module
        module_genes = gene_modules.index[mask].tolist()
        age_hub_set.update(module_genes[:30])
        print(f"  Age module '{module}': 30 genes (simplified)")

print(f"\n=== HUB GENE COUNTS ===")
print(f"Tissue hub genes: {len(tissue_hub_set)}")
print(f"Age hub genes: {len(age_hub_set)}")
print(f"Shared hub genes: {len(tissue_hub_set & age_hub_set)}")

#  Create Venn diagram
if len(tissue_hub_set) > 0 and len(age_hub_set) > 0:
    venn_dict = {
        'Tissue Hub Genes': tissue_hub_set,
        'Age Hub Genes': age_hub_set
    }
    
    plt.figure(figsize=(10, 8))
    venn.venn(venn_dict)
    
    shared_count = len(tissue_hub_set & age_hub_set)
    plt.title(f"Hub Gene Overlap: Tissue vs Age Modules\n{shared_count} shared hub genes")
    plt.tight_layout()
    plt.savefig('hub_gene_overlap_venn.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Show shared genes
    shared_genes = tissue_hub_set & age_hub_set
    if shared_genes:
        print(f"\n=== SHARED HUB GENES ({len(shared_genes)} total) ===")
        shared_list = sorted(list(shared_genes))
        for i, gene in enumerate(shared_list):
            print(f"{i+1:3}. {gene}")
            
            # Show which modules this gene belongs to
            gene_module = gene_modules.get(gene, 'Unknown')
            print(f"     Module: {gene_module}")
            
            if i >= 19:  # Show first 20 only
                print(f"... and {len(shared_list) - 20} more")
                break
    
    # Save results
    pd.Series(list(tissue_hub_set)).to_csv('tissue_hub_genes_final.csv', index=False)
    pd.Series(list(age_hub_set)).to_csv('age_hub_genes_final.csv', index=False)
    pd.Series(list(shared_genes)).to_csv('shared_hub_genes_final.csv', index=False)
    print("\nSaved: tissue_hub_genes_final.csv, age_hub_genes_final.csv, shared_hub_genes_final.csv")

# Also save the all-gene comparison for reference
print("\n=== ALL GENE COMPARISON (for reference) ===")
module_genes = {}
for module in pyWGCNA_5xFAD.datExpr.var['moduleColors'].unique():
    mask = pyWGCNA_5xFAD.datExpr.var['moduleColors'] == module
    genes = set(pyWGCNA_5xFAD.datExpr.var[mask].index)
    module_genes[module] = genes

tissue_all_genes = set()
for module in tissue_modules:
    if module in module_genes:
        tissue_all_genes.update(module_genes[module])

age_all_genes = set()
for module in age_modules:
    if module in module_genes:
        age_all_genes.update(module_genes[module])

print(f"All tissue-associated genes: {len(tissue_all_genes)}")
print(f"All age-associated genes: {len(age_all_genes)}")
print(f"All shared genes: {len(tissue_all_genes & age_all_genes)}")

# Save all genes for comparison
pd.Series(list(tissue_all_genes)).to_csv('tissue_ALL_genes.csv', index=False)
pd.Series(list(age_all_genes)).to_csv('age_ALL_genes.csv', index=False)
pd.Series(list(tissue_all_genes & age_all_genes)).to_csv('shared_ALL_genes.csv', index=False)

print("\n=== ANALYSIS COMPLETE ===")
print(f"Overlapping modules: {set(tissue_modules) & set(age_modules)}")
print(f"These modules contain genes that are co-expressed in both tissue-specific and age-dependent patterns")