import os

import pandas as pd

#create a function to read dataframes, concat two dataframes and delete the original files
def merge_dataframes(df1_path: str, df2_path: str, output_path: str):
    df1 = pd.read_csv(df1_path)
    df2 = pd.read_csv(df2_path)

    merged_df = pd.concat([df1, df2], ignore_index=True)
    merged_df.to_csv(output_path, index=False)

    print(f"Merged {df1_path} and {df2_path} into {output_path}")
    # Optionally delete original files
    os.remove(df1_path)
    os.remove(df2_path)
    
# automatically read files from output/dataframes and merge them into output/merged_data.csv
