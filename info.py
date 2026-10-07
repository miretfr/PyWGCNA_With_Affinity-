import PyWGCNA



# Load the object
pyWGCNA_5xFAD = PyWGCNA.readWGCNA("5xFAD.p")

# Check what attributes exist
print("Available attributes:")
attrs = [attr for attr in dir(pyWGCNA_5xFAD) if not attr.startswith('_')]
for attr in attrs:
    print(f"  - {attr}")

# Check specifically for metadata
print("\n" + "="*60)
print("Looking for metadata:")
print("="*60)

if hasattr(pyWGCNA_5xFAD, 'datMeta'):
    print(" datMeta exists")
    print(f"   Shape: {pyWGCNA_5xFAD.datMeta.shape}")
    print(f"   Columns: {pyWGCNA_5xFAD.datMeta.columns.tolist()}")
elif hasattr(pyWGCNA_5xFAD, 'sampleInfo'):
    print(" sampleInfo exists")
    print(f"   Shape: {pyWGCNA_5xFAD.sampleInfo.shape}")
    print(f"   Columns: {pyWGCNA_5xFAD.sampleInfo.columns.tolist()}")
elif hasattr(pyWGCNA_5xFAD, 'datExpr'):
    print("  datExpr.obs exists (sample metadata)")
    print(f"   Shape: {pyWGCNA_5xFAD.datExpr.obs.shape}")
    print(f"   Columns: {pyWGCNA_5xFAD.datExpr.obs.columns.tolist()}")
else:
    print("❌ No metadata found!")

# Check what modules you have
print("\n" + "="*60)
print("Modules:")
print("="*60)
if hasattr(pyWGCNA_5xFAD, 'datExpr'):
    modules = pyWGCNA_5xFAD.datExpr.var['moduleColors'].value_counts()
    print(modules)