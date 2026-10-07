# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 03:11:39 2025

@author: miret
"""
import matplotlib
matplotlib.use("Agg") 
import PyWGCNA
pyWGCNA_5xFAD = PyWGCNA.readWGCNA("5xFAD.p")
pyWGCNA_5xFAD.datExpr.var.head(5)

#pyWGCNA_5xFAD.CoexpressionModulePlot(modules=["purple"], numGenes=10, numConnections=100, minTOM=0)

modules = pyWGCNA_5xFAD.datExpr.var.moduleColors.unique().tolist()
pyWGCNA_5xFAD.CoexpressionModulePlot(modules=modules, numGenes=100, numConnections=1000, minTOM=0, file_name="all")