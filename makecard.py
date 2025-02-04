import ROOT
import os
import copy
import logging
from processes import MUTAU_SIGNAL, ETAU_SIGNAL, PURPLE_ETAU_E, PURPLE_MUTAU_E, BACKGROUND

SIGNAL_SCALE = 1
LUMI_SCALE = 1

# Systematic uncertainties

# Apply the lumi uncertainty to all processes
FIXED_UNCERTAINTIES = {
    "lumi": 0.025,
}

SIGNALS_GROUP = {**PURPLE_MUTAU_E, **PURPLE_ETAU_E}


# Apply the uncertainties to the background processes
UNCERTAINTIES = {
    "ttb_bg": {
        "ttbar_hvq": 0.13,
    },
    "WA_bg": {
        "WAToLNuA0123j": 0.11,
    },
    "DY_bg": {
        "DYToLLJets": 0.11,
        "DY2Jets": 0.11,
    },
}


# Test read root file and get the list of histograms name
def read_root_file(filename):
    f = ROOT.TFile.Open(filename)
    if f.IsZombie():
        print("Error: Cannot open file")
        return
    keys = f.GetListOfKeys()
    name_list = [key.GetName() for key in keys]
    f.Close()
    return name_list

def get_root_content(filepath: str) -> tuple:
    assert os.path.exists(filepath), f"File {filepath} does not exist."
    f = ROOT.TFile.Open(filepath)
    
    assert not f.IsZombie(), f"Error: Cannot open file {filepath}"
    
    hist_dict = {}
    for key in f.GetListOfKeys():
        hist = key.ReadObj()
        hist_dict[key.GetName()] = copy.deepcopy(hist)
    
    f.Close()
    return hist_dict    

def merge_hist(hist_dict: dict) -> ROOT.TH1:
    # Use .Add() to merge the histograms
    copied_hist_dict = copy.deepcopy(hist_dict)
    combined_hist = copied_hist_dict.popitem()[1]
    for hist in copied_hist_dict.values():
        combined_hist.Add(hist)
    return combined_hist


def get_info_from_name(name):  # Return, base name, category, jetgroup
    jetgroup = name[-2:]
    if len(name.split(f"_highmass_{jetgroup}")) > 1:
        return name.split(f"_highmass_{jetgroup}")[0], "highmass", jetgroup
    elif len(name.split(f"_lowmass_{jetgroup}")) > 1:
        return name.split(f"_lowmass_{jetgroup}")[0], "lowmass", jetgroup
    raise ValueError(f"Cannot determine category and jetgroup from name {name}")


def extract_list_backgrounds(name_list, category, jetgroup):
    backgrounds = []
    for name in name_list:
        base_name, cat, jg = get_info_from_name(name)
        if cat == category and jg == jetgroup:
            if base_name not in ["combined_background", "data_obs"] + list(
                SIGNALS_GROUP.keys()
            ):
                backgrounds.append(base_name)
    return backgrounds


# EXAMPLE_PATH:
background_dirs = {
    "ttbar_hvq": "..some_path",
    "WAToLNuA0123j": "..some_path",
}

signal_dir = {
    "mutauM200": "..some_path",
}

def make_datacard(background_dirs, signal_dir, cut_type, category, jetgroup, save_path="./datacard"):
    # cut_type: mutau, etau
    # category: lowmass, highmass
    # jetgroup: 0j, 1j
    
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("makecard")

    common_hist_name = f'weighted_{cut_type}_{category}_{jetgroup}'
    signal = list(signal_dir.keys())[0]
    signal_full_name = f"{signal}_{cut_type}_{category}_{jetgroup}"
    
    background_hist = {}
    for key, path in background_dirs.items():
        hist_dict = get_root_content(path)
        full_name = f"{key}_{cut_type}_{category}_{jetgroup}"
        background_hist[full_name] = hist_dict[common_hist_name]
    
    signal_hist = get_root_content(signal_dir[list(signal_dir.keys())[0]])[common_hist_name]
    combined_background = {f'combined_background_{cut_type}_{category}_{jetgroup}': merge_hist(background_hist)}
    
    # Save to save_path.root
    # Create save_path.root
    saved_root_name = save_path + ".root"
    f = ROOT.TFile(saved_root_name, "RECREATE")
    
    for key, hist in background_hist.items():
        hist.SetName(key)
        hist.Write()
        
    signal_hist.SetName(signal_full_name)
    signal_hist.Write()
    
    for key, hist in combined_background.items():
        hist.SetName(f"{key}_{cut_type}_{category}_{jetgroup}")
        hist.Write()
        
    f.Close()
    
    # Create datacard
    contents = []
    
    imax = 1 # number of channels
    jmax = len(background_hist) # number of backgrounds
    
    max_len = max([len(key) for key in list(background_hist.keys())] + [len(signal)]) + 1
    sep = "\t\t" * 2
    sep2 = "\t\t"
    
    # Write the header
    contents += [
        f"imax {imax} number of channels",
        f"jmax {jmax} number of backgrounds",
        "kmax * number of nuisance parameters",
        "-" * 100,
        f"shapes data_obs bin1 {saved_root_name} combined_background_{cut_type}_{category}_{jetgroup}",
        f"shapes * * {saved_root_name} $PROCESS",
        "-" * 100,
        "bin bin1",
        "observation -1",
        "-" * 100,
    ]
    # TODO: change to access weight_... instead of the name
    # Define the signal and background processes
    total_processes = list(background_hist.keys()) + [signal_full_name]
    bins = sep.join(["bin1".ljust(max_len)] * len(total_processes))
    processes = sep.join([p.ljust(max_len) for p in total_processes])
    process_idx = sep.join([str(i).ljust(max_len) for i in range(len(total_processes))])
    rates = sep.join(["-1".ljust(max_len)] * len(total_processes))
    contents += [
        "bin\t" + sep + bins,
        "process" + sep + processes,
        "process" + sep + process_idx,
        "rate" + sep + rates,
        "-" * 100,
    ]
    
    # Systematic uncertainties
    contents += ["# Systematic uncertainties"]
    # Apply the fixed uncertainties
    for name, value in FIXED_UNCERTAINTIES.items():
        contents.append(
            f"{name}{sep2}lnN{sep2}{sep.join([str(1+value).ljust(max_len)]*len(total_processes))}"
        )
        
    # Apply the uncertainties to the background processes
    for name, values in UNCERTAINTIES.items():
        for process, value in values.items():
            process_name = process + f"_{cut_type}_{category}_{jetgroup}"
            line = [
                "-".ljust(max_len) if process_name != p else str(1 + value).ljust(max_len)
                for p in total_processes
            ]
            contents.append(f"{name}{sep2}lnN{sep2}{sep.join(line)}")
            
    contents.append("* autoMCStats 0")
    contents.append(f"lumiscale rateParam * * {LUMI_SCALE}")
    contents.append(f"signalscale rateParam * {signal_full_name} {SIGNAL_SCALE}")
    
    final_contents = "\n".join(contents)
    
    with open(save_path + ".txt", "w") as f:
        f.write(final_contents)
    logger.info(f"Datacard is saved as {save_path}.txt")
    
    

# Test
all_masses = ['200', '300', '450', '450', '600', '750', '900']
all_categories = ['lowmass', 'lowmass', 'lowmass', 'highmass', 'highmass', 'highmass', 'highmass']
processes = {
    'mutau': {
        'mass': all_masses,
        'category': all_categories,
        'file': 'mutau_collinear_mass.root'
    },
    'etau': {
        'mass': all_masses,
        'category': all_categories,
        'file': 'etau_collinear_mass.root'
    },
}

backgrounds = {
    # "ttbar_hvq": "./outputs/ttbar_hvq/merged.root",
    "ttbar_hvq": "./merged_ttbar_hvq.root",
    # "WAToLNuA0123j": "some_path",
}

# make_datacard("collinear_mass.root", "mutau", "200", "lowmass", "0j")

# OLD
# for jetgroup in ["0j", "1j"]:
#     for cut_type, info in processes.items():
#         for mass, category in zip(info['mass'], info['category']):
#             make_datacard(info['file'], process, mass, category, jetgroup, save_path=f"{process}_{mass}_{category}_{jetgroup}.txt")


# NEW
signal_dir = {'purple_mutauM200': './merged_mutauM200.root'}

# test
make_datacard(backgrounds, signal_dir, 'mutau_e', 'lowmass', '0j', save_path="test_card")