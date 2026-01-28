!pip install rdkit pandas
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

df = pd.read_csv("simeonov_actives.tsv", sep="\t")
print(f"Loaded {len(df)} rows.")

def sanitize_smiles(smiles):
    if pd.isna(smiles):
        return None, False
    try:
        # Create a molecule object
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, False
            
        # Sanitize
        clean_mol = rdMolStandardize.Cleanup(mol)
        clean_smiles = Chem.MolToSmiles(clean_mol, isomericSmiles=True)
        
        # Check
        orig_canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        is_changed = (orig_canonical != clean_smiles)
        
        return clean_smiles, is_changed
    except:
        return None, False

print("Sanitizing... this might take 1-2 minutes...")
results = df['SMILES'].apply(sanitize_smiles)

# Save results back to the dataframe
df['SMILES_sanitized'] = [res[0] for res in results]
df['Was_Changed'] = [res[1] for res in results]

print(f"Total SMILES changed: {df['Was_Changed'].sum()}")

number_changed = df['Was_Changed'].sum()

print(f"The answer is: {number_changed}")

# Save the new file
df.to_csv("simeonov_actives_sanitized.tsv", sep="\t", index=False)
print("Sanitization complete.")

