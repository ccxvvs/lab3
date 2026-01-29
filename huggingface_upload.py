from huggingface_hub import HfApi, login
from datasets import Dataset
import pandas as pd
import os

# --- PART 1: LOGIN ---
print("Logging in...")
login()

# --- SETUP ---
api = HfApi()
username = api.whoami()['name']
repo_id = f"{username}/Simeonov2008"
print(f"Working on repository: {repo_id}")

# --- UPLOAD DATA ---
print("Loading sanitized data...")
# Load the file created in Step 3
df = pd.read_csv("simeonov_actives_sanitized.tsv", sep="\t")
dataset = Dataset.from_pandas(df)

print(f"Pushing data to {repo_id}...")
dataset.push_to_hub(repo_id)

# --- CREATE & UPLOAD CARD ---
print("Creating Dataset Card...")

# Define the text (README.md)
readme_content = f"""---
license: mit
tags:
- chemistry
- biology
---

# Simeonov2008 Active Compounds

## 1. Description and Citations
This dataset contains active compounds identified in a Quantitative High-Throughput Screening (qHTS) study. The compounds were screened for Luciferase inhibition to identify false positives in bioluminescence assays.

**Citation:**
> Simeonov A, et al. "Interference with bioluminescence imaging cross-talk in quantitative high-throughput screening." *J Med Chem.* 2008; 51(8): 2363-2374.

## 2. Quick-Start Guide
Use the following Python code to load this dataset:

```python
from datasets import load_dataset

# Load the dataset
dataset = load_dataset("{repo_id}")"""
