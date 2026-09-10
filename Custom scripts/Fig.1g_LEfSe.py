import pandas as pd
import numpy as np
from scipy import stats
import os

abund_path = r"/path/to/abundance_table.txt"
meta_path = r"/path/to/metadata.xlsx"
out_dir = r"/path/to/output_directory"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

print("Reading Metadata...")
df_meta = pd.read_excel(meta_path)
df_meta = df_meta[df_meta['origin'].isin(['animal', 'plant'])]
df_meta['Sample_id'] = df_meta['Sample_id'].astype(str).str.strip()
sample_to_origin = dict(zip(df_meta['Sample_id'], df_meta['origin']))

print("Reading abundance table...")
try:
    df_abund = pd.read_csv(abund_path, sep='\t', encoding='utf-16')
except UnicodeError:
    print("Failed to read with UTF-16, trying GBK/GB18030...")
    df_abund = pd.read_csv(abund_path, sep='\t', encoding='gb18030')

df_abund.columns = [str(c).strip() for c in df_abund.columns]

sample_cols = [col for col in df_abund.columns if col in sample_to_origin]
print(f"Successfully matched {len(sample_cols)} samples.")

if len(sample_cols) == 0:
    print("❌ Warning: No samples matched! Please check if the sample names in TXT exactly match the Sample_id in Excel.")
else:
    df_genus = df_abund[['genus'] + sample_cols].groupby('genus').sum(numeric_only=True)

    sample_sums = df_genus.sum(axis=0)
    sample_sums = sample_sums.replace(0, 1) 
    df_genus_norm = df_genus.div(sample_sums, axis=1) * 1000000

    print("Performing differential statistical analysis...")
    results = []
    for genus in df_genus_norm.index:
        abundances = df_genus_norm.loc[genus]
        
        animal_abund = abundances[[s for s in sample_cols if sample_to_origin[s] == 'animal']].values
        plant_abund = abundances[[s for s in sample_cols if sample_to_origin[s] == 'plant']].values
        
        if np.sum(animal_abund) + np.sum(plant_abund) == 0:
            continue
            
        try:
            stat, p_value = stats.kruskal(animal_abund, plant_abund)
        except ValueError:
            p_value = np.nan
            
        mean_animal = np.mean(animal_abund)
        mean_plant = np.mean(plant_abund)
        
        if mean_animal > mean_plant:
            enriched = 'animal'
            diff = mean_animal - mean_plant
        elif mean_plant > mean_animal:
            enriched = 'plant'
            diff = mean_plant - mean_animal
        else:
            enriched = 'None'
            diff = 0
            
        lda_score = np.log10(diff) if diff > 1 else 0
        
        results.append({
            'Genus': genus,
            'Enriched_in': enriched,
            'p_value': p_value,
            'LDA_Score': lda_score,
            'Mean_animal_CPM': mean_animal,
            'Mean_plant_CPM': mean_plant
        })

    df_res = pd.DataFrame(results)
    if not df_res.empty:
        df_res = df_res.dropna(subset=['p_value'])
        df_res = df_res.sort_values('p_value')
        out_res_path = os.path.join(out_dir, "Predator_LEfSe_Calculation_Results.csv")
        df_res.to_csv(out_res_path, index=False, encoding='utf-8-sig')

    print("Generating official LEfSe format input file...")
    lefse_input = pd.DataFrame(index=['Origin_Class'] + list(df_genus.index), columns=sample_cols)

    for s in sample_cols:
        lefse_input.loc['Origin_Class', s] = sample_to_origin[s]

    lefse_input.loc[df_genus.index, sample_cols] = df_genus.loc[df_genus.index, sample_cols]

    out_lefse_path = os.path.join(out_dir, "Predator_LEfSe_Official_Input.txt")
    lefse_input.index.name = 'Sample_ID'
    lefse_input.to_csv(out_lefse_path, sep='\t', header=True, index=True)

    print(f"\n✅ Analysis complete! All results output to: {out_dir}")
