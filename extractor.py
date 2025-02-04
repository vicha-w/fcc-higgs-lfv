#!/work/app/modules/software/Python/3.9.6-GCCcore-11.2.0/bin/python
import ROOT
import os
import pandas as pd



from processes import BACKGROUND, MUTAU_SIGNAL, ETAU_SIGNAL, PURPLE_ETAU_E, PURPLE_MUTAU_E

ROOT.gInterpreter.ProcessLine('gErrorIgnoreLevel = kFatal;')

ALL_PROCESSES = {**BACKGROUND, **MUTAU_SIGNAL, **ETAU_SIGNAL, **PURPLE_ETAU_E, **PURPLE_MUTAU_E}


# Define the branches to extract
BRANCHES = {
    "Jet": ["PT", "Eta", "Phi", "BTag", "TauTag", "Charge", "Mass"],
    "Muon": ["PT", "Eta", "Phi", "Charge"],
    "Electron": ["PT", "Eta", "Phi", "Charge"],
    "MissingET": ["Phi", "MET"],
    "Event": ["Number"],
}

extracted_file_universal_path = "/work/project/physics/psriling/FCC/extraction/"

def extract_columns(input_file_path, output_file_path, branches_dict):
    # Open the input ROOT file
    input_file = ROOT.TFile.Open(input_file_path, "READ")
    if not input_file:
        print(f"Error: Could not open input file {input_file_path}")
        return

    # Create the output ROOT file
    output_file = ROOT.TFile.Open(output_file_path, "RECREATE")
    if not output_file:
        print(f"Error: Could not create output file {output_file_path}")
        return

    # Get the Delphes tree from the input file
    tree = input_file.Get("Delphes")
    if not tree:
        print("Error: Delphes tree not found in input file")
        return

    # Create a new tree in the output file
    output_tree = tree.CloneTree(0)

    # Enable only the desired branches
    tree.SetBranchStatus("*", 0)
    for key, branches in branches_dict.items():
        for branch in branches:
            full_branch_name = f"{key}.{branch}"
            tree.SetBranchStatus(full_branch_name, 1)
    tree.SetBranchStatus('Jet_size', 1)
    tree.SetBranchStatus('Muon_size', 1)
    tree.SetBranchStatus('Electron_size', 1)
    # Fill the new tree with the selected branches
    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        output_tree.Fill()

    # Write the new tree to the output file
    output_tree.Write()

    # Close the files
    input_file.Close()
    output_file.Close()

def get_files(path) -> list:
    files = []
    mapped = [path] if isinstance(path, str) else path  # Mapped to list
    for element in mapped:
        if os.path.isdir(element):
            files.extend(
                [f"{element}/{f}" for f in os.listdir(element) if f.endswith(".root")]
            )
        else:
            files.append(element)

    return files

def process_files(process_name, files):
    output_dir = os.path.join(extracted_file_universal_path, process_name)
    os.makedirs(output_dir, exist_ok=True)
    
    info_list = []
    processed_files = []
    existing_csv = False
    # read the info.csv to get the processed files
    if os.path.exists(os.path.join(output_dir, "info.csv")):
        old_info_df = pd.read_csv(os.path.join(output_dir, "info.csv"))
        processed_files = old_info_df["original_file"].tolist()
        existing_csv = True
    
    for idx, input_file in enumerate(files):
        print(f"Processing file {idx+1}/{len(files)}: {input_file}")
        
        if existing_csv and input_file in processed_files:
            print(f"File {input_file} already processed. Skipping...")
            continue
        
        original_size = os.path.getsize(input_file)
        output_file = os.path.join(output_dir, f"extracted_{idx}_{os.path.basename(input_file)}")
        
        extract_columns(input_file, output_file, BRANCHES)
        
        extracted_size = os.path.getsize(output_file)
        
        info_list.append({
            "ID": idx,
            "original_file": input_file,
            "original_size": original_size,
            "extracted_size": extracted_size
        })
    info_df = pd.DataFrame(info_list)
    info_df['reduction_pct'] = (1 - info_df['extracted_size'] / info_df['original_size']) * 100
    
    if existing_csv:
        info_df = pd.concat([old_info_df, info_df], ignore_index=True)
    info_df.to_csv(os.path.join(output_dir, "info.csv"), index=False)

# Example usage
# BACKGROUND = {
#     "ttbar_hvq": {
#         "path": [
#             "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/hvqMulti1M",        # 1M
#             "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/hvqMulti10M",       # 10M set1
#             "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/hvqMulti10M_set2",  # 10M set2
#         ],
#         "cross-section": 31429.120490154033, # pb - POWHEG
#     },
# }

# ttbar
process_name = "ttbar_hvq"
files = get_files(ALL_PROCESSES[process_name]["path"])
process_files(process_name, files)
