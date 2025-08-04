# Only use this when figures are ready to be edited.
import pickle
import numpy as np
import pandas as pd
import sys
from numba import njit, prange
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')
sns.set(rc={'figure.figsize':(11.7,8.27)})


# Define a global mapping of modifications to ML_Tag indices
MODIFICATION_MAPPING = {
    "m6A": 0,          # ML_Tag1
    "other_mod": 1,          # ML_Tag2
    "another_mod": 2,    # ML_Tag3
    "psi": 3   # ML_Tag4
}

def load_pickle_file(file_path, debug=False):
    try:
        with open(file_path, "rb") as f:
            obj = pickle.load(f)
        if debug:
            print(f"[DEBUG] Successfully loaded pickle file: {file_path}")
        return obj
    except FileNotFoundError:
        print(f"[ERROR] File '{file_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] An error occurred while loading '{file_path}': {e}")
        sys.exit(1)

def prepare_data(read_arrays_dict, read_to_ref_dict, allowed_ref_names=None, debug=False):
    # If allowed_ref_names is provided, intersect them with what we actually have from read_to_ref_dict
    ref_names = sorted(set(read_to_ref_dict.values()).intersection(allowed_ref_names)) \
                if allowed_ref_names else sorted(set(read_to_ref_dict.values()))
    ref_name_to_idx = {name: idx for idx, name in enumerate(ref_names)}

    # Subtract 3 from the reference position if name ends with '_IVT', else 0
    subtract_amount = np.array([3 if name.endswith("_IVT") else 0 for name in ref_names], dtype=np.int16)

    all_read_arrays = []
    all_read_to_ref_idx = []
    
    # Determine the actual structure of the arrays (number of columns)
    if read_arrays_dict:
        sample_key = next(iter(read_arrays_dict))
        sample_array = read_arrays_dict[sample_key]
        array_length = len(sample_array) if np.isscalar(sample_array[0]) else sample_array.shape[1]
        if debug:
            print(f"[DEBUG] Sample array length/width: {array_length}")
            print(f"[DEBUG] Sample array: {sample_array[:5] if np.isscalar(sample_array[0]) else sample_array[:5, :]}")
    
    for read_id, ref_name in read_to_ref_dict.items():
        if allowed_ref_names and ref_name not in allowed_ref_names:
            continue
        arr = read_arrays_dict.get(read_id)
        if arr is not None:
            # Ensure arr is 2D
            if arr.ndim == 1:
                arr = arr.reshape(1, -1)
            
            # No need to pad or truncate, we'll handle variable column lengths
            all_read_arrays.append(arr)
            
            # Extend reference indices
            ref_idx = ref_name_to_idx[ref_name]
            all_read_to_ref_idx.extend([ref_idx] * arr.shape[0])

    # Convert to numpy arrays
    all_read_arrays = np.vstack(all_read_arrays) if all_read_arrays else np.empty((0, 1), dtype=np.float64)
    all_read_to_ref_idx = np.array(all_read_to_ref_idx, dtype=np.int32)
    num_ml_scores = all_read_arrays.shape[1] - 3  # Number of ML score columns

    if debug:
        print(f"[DEBUG] Ref names: {ref_names}")
        print(f"[DEBUG] Ref name to index mapping: {ref_name_to_idx}")
        print(f"[DEBUG] All Read Arrays Shape: {all_read_arrays.shape}")
        print(f"[DEBUG] All Read to Ref Indices Shape: {all_read_to_ref_idx.shape}")
        print(f"[DEBUG] Subtract Amount Shape: {subtract_amount.shape}")
        print(f"[DEBUG] Number of ML Scores: {num_ml_scores}")

    return all_read_arrays, all_read_to_ref_idx, subtract_amount, num_ml_scores, ref_names, ref_name_to_idx

@njit(parallel=False)
def gather_data(all_read_to_ref_idx, all_read_arrays, subtract_amount, num_ml_scores, target_ml_tag_idx=None):
    num_rows = all_read_arrays.shape[0]

    # Prepare output arrays
    ref_idx_list = np.empty(num_rows, dtype=np.uint8)
    ref_pos_list = np.empty(num_rows, dtype=np.int16)
    ml_score_list = np.empty(num_rows, dtype=np.uint8)
    condition_list = np.empty(num_rows, dtype=np.uint8)

    # Always use the last column if only one ML score, otherwise use target_ml_tag_idx
    if num_ml_scores == 1 or target_ml_tag_idx is None:
        # Last column is the ML score column
        ml_idx = all_read_arrays.shape[1] - 1
    else:
        # Use the specified ML score column based on target_ml_tag_idx
        ml_idx = 3 + target_ml_tag_idx

    count = 0
    for i in prange(num_rows):
        ref_idx = all_read_to_ref_idx[i]
        if ref_idx == -1:
            continue
        subtract = subtract_amount[ref_idx]
        ref_pos = int(all_read_arrays[i, 1] - subtract)
        
        # Get the ML score from the appropriate column
        # Make sure ml_idx is within bounds
        if ml_idx < all_read_arrays.shape[1]:
            score = all_read_arrays[i, ml_idx]
        else:
            continue

        if score != -1:
            ref_idx_list[count] = ref_idx
            ref_pos_list[count] = ref_pos
            ml_score_list[count] = score
            # condition: 1 for IVT (subtract==3), else 0
            condition_list[count] = 1 if subtract == 3 else 0
            count += 1

    return (
        ref_idx_list[:count],
        ref_pos_list[:count],
        ml_score_list[:count],
        condition_list[:count]
    )

def filter_gathered_data(ref_idx, ref_pos, ml_score, condition, positions_dict, ref_name_to_idx, debug=False):
    """
    Filter the gathered data down to only the rows where (ref_idx, ref_pos)
    matches the positions_dict allowed positions for each reference.
    """
    valid_indices = []
    for ref_name, allowed_positions in positions_dict.items():
        if ref_name in ref_name_to_idx:
            idx = ref_name_to_idx[ref_name]
            mask = (ref_idx == idx) & np.isin(ref_pos, allowed_positions)

            if debug:
                print(f"[DEBUG] Filtering for {ref_name}: Allowed positions: {allowed_positions}")
                print(f"[DEBUG] Rows before filtering: {np.sum(ref_idx == idx)}")
                unique_positions_before = np.unique(ref_pos[ref_idx == idx])
                print(f"[DEBUG] Unique ref_pos values before filtering: {unique_positions_before}")
                print(f"[DEBUG] Mask for {ref_name}: {mask}")
                print(f"[DEBUG] Rows after filtering: {np.sum(mask)}")

            valid_indices.append(np.where(mask)[0])

    if valid_indices:
        valid_indices = np.concatenate(valid_indices)
    else:
        valid_indices = np.array([], dtype=int)

    return (
        ref_idx[valid_indices],
        ref_pos[valid_indices],
        ml_score[valid_indices],
        condition[valid_indices]
    )

def convert_to_dataframe(ref_idx, ref_pos, ml_score, condition, ref_names_sorted):
    df = pd.DataFrame({
        'ref_idx': ref_idx,
        'ref_pos': ref_pos,
        'ML_score': ml_score,
        'condition': np.where(condition == 1, "IVT", "modified")
    })
    df['ref_name'] = df['ref_idx'].apply(lambda x: ref_names_sorted[x])
    return df[['ref_idx', 'ref_name', 'ref_pos', 'ML_score', 'condition']]

def process_data(main_pkl, ref_lookup_pkl, mod="m6A", positions_dict=None, debug=False):
    # Load pickled data
    read_arrays_dict = load_pickle_file(main_pkl, debug)
    read_to_ref_dict = load_pickle_file(ref_lookup_pkl, debug)

    # Check the structure of the data to determine if we have a single ML score or multiple
    if read_arrays_dict:
        sample_key = next(iter(read_arrays_dict))
        sample_array = read_arrays_dict[sample_key]
        
        # Determine array structure
        if np.isscalar(sample_array[0]):  # 1D array
            ml_array_length = len(sample_array)
            num_ml_scores_expected = ml_array_length - 3  # Assuming [read_pos, ref_pos, numerical_base, ML_scores...]
        else:  # 2D array
            ml_array_length = sample_array.shape[1]
            num_ml_scores_expected = ml_array_length - 3
            
        if debug:
            print(f"[DEBUG] Sample array length: {ml_array_length}")
            print(f"[DEBUG] Expected number of ML scores: {num_ml_scores_expected}")
            print(f"[DEBUG] Sample array: {sample_array[:5] if np.isscalar(sample_array[0]) else sample_array[:5, :]}")
    
    # Convert the mod to a target ML index
    target_ml_tag_idx = MODIFICATION_MAPPING.get(mod)
    
    # Check if we should use the last column (only one ML score)
    if num_ml_scores_expected == 1:
        if debug:
            print(f"[INFO] Only one ML score column detected. Using this column regardless of specified mod ({mod}).")
        target_ml_tag_idx = None  # Signal to use the last column
    elif target_ml_tag_idx is None:
        if debug:
            print(f"[WARNING] Invalid modification '{mod}'. Will use the last ML score column as fallback.")

    # Debug Checks
    if debug:
        print(f"[DEBUG] positions_dict keys: {sorted(positions_dict.keys())}")
        unique_names_in_read_to_ref = set(read_to_ref_dict.values())
        print(f"[DEBUG] Unique reference names in read_to_ref_dict: {sorted(unique_names_in_read_to_ref)}")
        print(f"[DEBUG] Number of read_arrays_dict entries: {len(read_arrays_dict)}")

    # Prepare arrays for processing
    all_read_arrays, all_read_to_ref_idx, subtract_amount, num_ml_scores, ref_names_sorted, ref_name_to_idx = prepare_data(
        read_arrays_dict,
        read_to_ref_dict,
        positions_dict.keys() if positions_dict else None,
        debug
    )

    # Gather data into single arrays (Numba-accelerated)
    ref_idx, ref_pos, ml_score, condition = gather_data(
        all_read_to_ref_idx, 
        all_read_arrays, 
        subtract_amount, 
        num_ml_scores, 
        target_ml_tag_idx
    )

    # Filter gathered data by allowed positions
    ref_idx, ref_pos, ml_score, condition = filter_gathered_data(
        ref_idx,
        ref_pos,
        ml_score,
        condition,
        positions_dict,
        ref_name_to_idx,
        debug
    )

    # Convert to a final DataFrame
    filtered_df = convert_to_dataframe(
        ref_idx,
        ref_pos,
        ml_score,
        condition,
        ref_names_sorted
    )

    if debug:
        print(f"[DEBUG] Final DataFrame shape: {filtered_df.shape}")
        print(filtered_df.head())

    # Check if DataFrame is empty and provide sample data for visualization if needed
    if filtered_df.empty:
        print("[WARNING] The resulting DataFrame is empty. Check your filtering criteria.")

    return filtered_df

# Usage
if __name__ == "__main__": 
    main_pickle = "/work/Genometechlab/andrew/RNA002_vs_RNA004/RNA004/Synthetics/Read_Filtering_and_Boxen_Plots/data/12_16_24_Mod_IVT_dorado_0.8.1_ino_m6A_no-trim_emit-moves_threshold_0_aligned_sorted_filtered_for_valid_reads_data_structure.pkl"
    ref_lookup_pickle = "/work/Genometechlab/andrew/RNA002_vs_RNA004/RNA004/Synthetics/Read_Filtering_and_Boxen_Plots/bam/12_16_24_Mod_IVT_dorado_0.8.1_ino_m6A_no-trim_emit-moves_threshold_0_aligned_sorted_filtered_for_valid_reads_read_id_to_reference.pkl"
    positions_dict_m6a = {
        'Both_High': [53, 56, 60, 64, 65],
        'Both_High_IVT': [53, 56, 60, 64, 65],
        'Dorado_High_m6Anet_Low': [56, 59, 60, 63, 64],
        'Dorado_High_m6Anet_Low_IVT': [56, 59, 60, 63, 64],
        'm6Anet_High_Dorado_Low': [54, 55, 60, 64, 68],
        'm6Anet_High_Dorado_Low_IVT': [54, 55, 60, 64, 68]
    }
    positions_dict_psi = {
        'psU_PhMh' : [51, 57, 60, 67, 69], 
        'psU_PhMh_IVT' : [51, 57, 60, 67, 69],
        'psU_PhMl' : [51, 57, 60, 62, 65], 
        'psU_PhMl_IVT' : [51, 57, 60, 62, 65],
        'psU_PlMh' : [48, 53, 60, 64, 82], 
        'psU_PlMh_IVT' : [48, 53, 60, 64, 82] 
    }

    filtered_df = process_data(
        main_pickle, 
        ref_lookup_pickle, 
        mod="m6A", 
        positions_dict=positions_dict_m6a, 
        debug=True
    )


# filepath
# filepath = "/scratch/stein.an/RNA004_Synthetics/analysis/data/12_16_24_filtered_data_updated.csv"

ref_names = {}

# df = pd.read_csv(filepath)
# print(df.keys())
df = filtered_df

# Step 1: Create 'group_key' by removing '_IVT' suffix
df['group_key'] = df['ref_name'].str.replace('_IVT$', '', regex=True)

# Step 2: Set seaborn style without grid
sns.set(style="white")  # Changed from "whitegrid" to "white"

# Step 3: Identify unique groups
groups = df['group_key'].unique()
num_groups = len(groups)
print(f"\nThe number of unique groups for plotting is: {num_groups}")  # QC Check

# Step 4: Create subplots: num_groups rows for groups, 1 column
fig, axes = plt.subplots(nrows=num_groups, ncols=1, figsize=(14, 5 * num_groups), sharex=False)

# If only one group, axes may not be an array
if num_groups == 1:
    axes = [axes]

for i, (ax, group) in enumerate(zip(axes, groups)):
    # Subset data for the current group
    group_data = df[df['group_key'] == group]
    
    # Create box plot
    sns.boxenplot(
        x='ref_pos',
        y='ML_score',
        hue='condition',
        data=group_data,
        palette='Set1',
        ax=ax,
        dodge=True,
        k_depth=10
    )
    
    # Define the threshold
    threshold = 255 * 0.70
    
    if i == 0:
        # Add horizontal line with label only in the first subplot
        line = ax.axhline(y=threshold, color='black', linestyle='--', linewidth=1.5, label='modkit threshold')
    else:
        # Add horizontal line without label in other subplots
        ax.axhline(y=threshold, color='black', linestyle='--', linewidth=1.5)
    
    # Customize the subplot
    ax.set_title(f'ML_score Box Plot for {group}', fontsize=16)
    ax.set_xlabel('Reference Position', fontsize=14)
    ax.set_ylabel('ML Score', fontsize=14)
    
    # Update the legend
    ax.legend(title='Condition', fontsize=12, title_fontsize=13)

plt.tight_layout()
plt.savefig("/work/Genometechlab/andrew/RNA002_vs_RNA004/RNA004/Synthetics/Read_Filtering_and_Boxen_Plots/plots/12_16_24_m6A_no_trim_valid_reads_boxen_plot.pdf")
