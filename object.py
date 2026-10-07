# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 02:28:28 2025

@author: miret
"""
import PyWGCNA
import numpy as np
import pandas as pd
from sklearn.cluster import AffinityPropagation
import time
from PyWGCNA import WGCNA
import matplotlib.pyplot as plt





OKGREEN = "\033[92m"
OKBLUE = "\033[94m"
ENDC = "\033[0m"
BOLD = "\033[1m"
import PyWGCNA
pyWGCNA_5xFAD = PyWGCNA.readWGCNA("5xFAD.p")

print(f"Raw expresion data along with information:\n {pyWGCNA_5xFAD.geneExpr}")

print(f"Raw sample inforamtions:")
pyWGCNA_5xFAD.geneExpr.obs.head(5)

print(f"Raw gene inforamtions:")
pyWGCNA_5xFAD.geneExpr.var.head(5)


print(f"Processed expresion data along with information:\n {pyWGCNA_5xFAD.datExpr}")

print(f"Raw sample inforamtions:")
pyWGCNA_5xFAD.datExpr.obs.head(5)

print(f"Raw gene inforamtions:")
pyWGCNA_5xFAD.datExpr.var.head(5)

geneList = PyWGCNA.getGeneList(dataset='mmusculus_gene_ensembl',
                               attributes=['ensembl_gene_id', 
                                           'external_gene_name', 
                                           'gene_biotype'],
                               maps=['gene_id', 'gene_name', 'gene_biotype'])


pyWGCNA_5xFAD.updateGeneInfo(geneList)


pyWGCNA_5xFAD.updateSampleInfo(path='5xFAD_paper/sampleInfo.csv', sep=',')


pyWGCNA_5xFAD.sft

pyWGCNA_5xFAD.adjacency

pyWGCNA_5xFAD.TOM

