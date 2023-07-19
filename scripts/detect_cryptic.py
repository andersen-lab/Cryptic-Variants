import argparse
import pandas as pd
from outbreak_data import outbreak_data
import re

parser = argparse.ArgumentParser()

parser.add_argument("covar_file", help="Path to freyja covariants output")
parser.add_argument(
    "--max_clinical_count",
    help="Maximum number of GISAID sequences for a cluster to be saved",
    default=10,
)
parser.add_argument("--min_cluster_size", help="Minimum cluster size", default=2)
parser.add_argument("--min_site", help="Minimum genomic region", default=22556)
parser.add_argument("--max_site", help="Maximum genomic region", default=23156)

parser.add_argument("--location_id", help="Location id to query",
                    default='global')
parser.add_argument("-o", "--output", help="Output file",
                    default="cryptic_var.tsv")

args = parser.parse_args()

def sort_cluster(variant):
    if 'DEL' in variant:
        return int(variant.split('DEL')[0][:-1])
    else:
        return int(variant.split(':')[1][1:-2])
    
def extract_gene_aa_mutation(cluster):
    # Parse freyja covariants output
    cluster_final = []
    for variant in cluster:
        if ":" in variant:
            if "*" in variant or "INS" in variant:
                continue  # Ignore stop codons and insertions
            elif "DEL" in variant:
                if int(variant.split(")(")[0].split(",")[1]) % 3 != 0:
                    continue # Ignore single nt deletions
                else:
                    site = int(variant.split(',')[0][1:])
                    if (site >= int(args.min_site) and site <= int(args.max_site)):
                        cluster_final.append(variant.split(")(")[1][:-1]) # Deletion
            else:
                site = int(variant.split('(')[0][1:-1])
                if (site >= int(args.min_site) and site <= int(args.max_site)):
                    cluster_final.append(variant.split("(")[1][:-1]) # SNV

    if len(cluster_final) < int(args.min_cluster_size):
        return pd.NA
    
    cluster_final = list(set(cluster_final))
    cluster_final.sort(key=sort_cluster)

    return cluster_final

def get_clinical_data(cluster):
    mutations = ",".join(cluster)
    
    if args.location_id == 'global':
        query = f"mutations={mutations}"
    else:
        query = f"mutations={mutations}&location={args.location_id}"
        
    try:
        results = outbreak_data.get_outbreak_data(
            "genomics/mutations-by-lineage", argstring=query
        )["results"]
    except NameError as e:  # No clinical results found
        return pd.Series([0, pd.NA])

    results = list(results[mutations])

    count = 0
    lins_found = []
    for dict in results:
        count += dict["mutation_count"]
        lins_found.append(dict["pangolin_lineage"])

    return pd.Series([count, lins_found])

df = pd.read_csv(args.covar_file, sep="\t")
df = df.rename(columns={"Count": "WW_Count"})
df = df.drop(columns=["Coverage_start", "Coverage_end"])

# Extract gene-aa mutation from mutation cluster
df["Covariants"] = df["Covariants"].apply(lambda x: x.split(" "))
df["Covariants"] = df["Covariants"].apply(extract_gene_aa_mutation)

df = df.dropna()

# Aggregate by mutation cluster
df = df.groupby(df['Covariants'].astype(str))\
       .aggregate({'Covariants': 'first', 'WW_Count': 'sum'})
df = df[df["Covariants"].map(len) > 0]

# Get clinical data for each mutation cluster
success = False
try:
    df[["Clinical_Count", "Lineages"]] = df["Covariants"]\
                                        .apply(get_clinical_data)
    success = True
except ValueError:
    print("Empty covariants column found. Skipping this sample.")

if success:
    # Select clusters with clinical counts below threshold
    df = df.fillna('NA')
    df = df[df["Clinical_Count"] < int(args.max_clinical_count)]
    df = df.sort_values(by=["Clinical_Count"], ascending=True)

    # Save to file if there are cryptic variants present
    df.to_csv(args.output, sep="\t", index=False)
else:
    # Save empty dataframe
    colnames = ["Covariants", "WW_Count", "Clinical_Count", "Lineages"]
    pd.DataFrame(columns=colnames).to_csv(args.output, sep="\t", index=False)

print(f"Output saved to {args.output}")