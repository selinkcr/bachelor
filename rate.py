import pandas as pd
import ace_tools_open as tools; 

# Load the tab-separated file
# Load the tab-separated file
file_path = "mapping_21apr.txt"
df = pd.read_csv(file_path, sep="\t")

# Select relevant columns for growth rate analysis
growth_data = df[["taxon_oid","ncbi_taxon","genome_name", "Genus", "Species", "doubling_predicted"]]

# Sort by predicted doubling time (lower = faster growth)
growth_data_sorted = growth_data.sort_values("doubling_predicted")

# Save the sorted DataFrame to a new TSV file
growth_data_sorted.to_csv("growth_rates_sorted.tsv", sep="\t", index=False)

# Display to user (optional, can remove if just writing the file)
tools.display_dataframe_to_user(name="Bacterial Growth Rates", dataframe=growth_data_sorted)