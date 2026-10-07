import PyWGCNA


geneExp = '5xFAD_paper/expressionList.csv'
pyWGCNA_5xFAD = PyWGCNA.WGCNA(name='5xFAD', 
                              species='mus musculus', 
                              geneExpPath=geneExp, 
                              outputPath='',
                              save=True,
                              TPMcutoff=100.0)

# A class used to create gene expression anndata along data trait including both genes and samples information.
pyWGCNA_5xFAD.geneExpr.to_df().head(5)

# Preprocessing PyWGCNA object including removing obvious outlier on genes and samples (plot)
pyWGCNA_5xFAD.preprocess() 

#pyWGCNA_5xFAD.saveWGCNA()

# Clustering genes through original WGCNA pipeline: 
# 1.pick soft threshold 
# 2.calculating adjacency matrix 
# 3.calculating TOM similarity matrix 
# 4.cluster genes base of dissTOM (***)
# 5.merge similar cluster dynamically
pyWGCNA_5xFAD.findModules()

# add/update metadata in datExpr and geneExpr anndata
pyWGCNA_5xFAD.updateSampleInfo(path='5xFAD_paper/sampleInfo.csv', sep=',')

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


# add/update genes info in datExpr and geneExpr anndata
pyWGCNA_5xFAD.updateGeneInfo(geneList)

# bar plot of module eigen gene figure in given module , but with an unstopable error
# Just disable barplot with one line: 
PyWGCNA.WGCNA.barplotModuleEigenGene = lambda *args, **kwargs: None

# Analysing results: 1.calculating module trait relationship 2.plotting module heatmap eigengene 3.finding GO term for each module
pyWGCNA_5xFAD.analyseWGCNA()

# Saves the current WGCNA in pickle format with the .p extension
pyWGCNA_5xFAD.saveWGCNA()

# find top n hub genes based on connectivity in given module
pyWGCNA_5xFAD.top_n_hub_genes(moduleName="coral", n=10)