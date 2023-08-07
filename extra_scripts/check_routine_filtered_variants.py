import json
import pandas as pd

from pathlib import Path

ROOT_DIR = Path(__file__).absolute().parents[1]

# Open the JSON with all the GM numbers and their file info
with open(
        ROOT_DIR.joinpath('resources', 'sample_file_IDs_outcome.json'),
        "r", encoding='utf8'
) as json_file:
        json_dict = json.load(json_file)

# For each sample and its info, get the X number and add a key with the X number

for gm_number, vcf_info in json_dict.items():
    file_name = vcf_info.get('filename')
    if file_name:
        if 'GM' in file_name.upper():
             x_number = file_name.split('-')[0]
        else:
            x_number = file_name.split('_')[0]
        json_dict[gm_number]['x_number'] = x_number
    else:
         json_dict[gm_number]['x_number'] = ''


with open(ROOT_DIR.joinpath('resources', 'total_included_variants.txt')) as f:
     lines = f.readlines()

# Group the two lines together which are the file name and the number of variants
stripped = [line.replace("\n"," ").strip() for line in lines]
grouped = [stripped[i:i+2] for i in range(0, len(stripped), 2)]

# Split the filename to get the X number with the number of included variants
new_list = []
for nested_list in grouped:
    file_name = nested_list[0]
    if 'GM' in file_name.upper():
        x_number = file_name.split('-')[0]
    else:
        x_number = file_name.split('_')[0]
    updated_list = [x_number, nested_list[1]]
    new_list.append(updated_list)

# Create df with two columns - x number and the number of variants remaining
# after routine filtering
df = pd.DataFrame(new_list, columns = ['x_number', 'routine_filter_variants'])

# Create df from the big JSON with each GM number and info about the VCF found
gm_df = pd.DataFrame.from_dict(json_dict, orient='index').reset_index()

# Subset to only have the GM number, X number and the report outcome
subset_gm = gm_df[["index", "x_number", "report_outcome"]]

# Merge to get one df with the GM number, X number, report outcome and no
# of variants left after routine filtering
result = pd.merge(subset_gm ,df , on='x_number', how='left')

# Remove NAs (the 4 samples with no VCF found)
new = result[result['routine_filter_variants'].notna()]
new = new.rename(columns={"index": "gm_number"})

# Write out to new CSV
new.to_csv('variants_with_routine_filters.csv', index=False)


# Find mean number of variants returned to scientists for positive cases
positive_cases = new.loc[new.report_outcome.isin(['MUT','UVAR'])]

positive_cases['routine_filter_variants'] = positive_cases['routine_filter_variants'].astype(int)
print(positive_cases['routine_filter_variants'].mean())

# Do same for NMD cases
nmd_cases = new.loc[new.report_outcome == 'NMD']
nmd_cases['routine_filter_variants'] = nmd_cases['routine_filter_variants'].astype(int)
print(nmd_cases['routine_filter_variants'].mean())
