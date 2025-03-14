import os
import sys
import logging
import time
import argparse
import ROOT
import pandas as pd
import numpy as np
from multiprocessing import Pool
from lhe_split import split_lhe

# ********** PATHS **********
MG5_PYTHIA_INTERFACE = "/work/app/pythia8/MGInterface/1.3/MG5aMC_PY8_interface"
DELPHES_PATH = "/work/home/psriling/MsProject/Delphes-3.5.0/DelphesHepMC2"
DELPEHS_CARD = "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/FCChh_I.tcl"


ROOT.gInterpreter.ProcessLine("gErrorIgnoreLevel = kFatal;")
os.system(
    "export ROOT_INCLUDE_PATH=/work/home/psriling/MsProject/Delphes-3.5.0/external:$ROOT_INCLUDE_PATH"
)


# ### BUDDHA BLESS THIS CODE ###

"""""" """""" """""" """""" """""" """""" """""" """
                           _
                        _ooOoo_
                       o8888888o
                       88" . "88
                       (| -_- |)
                       O\  =  /O
                    ____/`---'\____
                  .'  \\|     |//  `.
                 /  \\|||  :  |||//  \
                /  _||||| -:- |||||_  \
                |   | \\\  -  /'| |   |
                | \_|  `\`---'//  |_/ |
                \  .-\__ `-. -'__/-.  /
              ___`. .'  /--.--\  `. .'___
           ."" '<  `.___\_<|>_/___.' _> \"".
          | | :  `- \`. ;`. _/; .'/ /  .' ; |
          \  \ `-.   \_\_`. _.'_/_/  -' _.' /
===========`-.`___`-.__\ \___  /__.-'_.'_.-'================
""" """""" """""" """""" """""" """""" """""" """"""


def get_args():
    ##### ARGUMENTS #####
    parser = argparse.ArgumentParser(description="Split LHE files")
    parser.add_argument("--input_file"  ,   type=str,     help="Input LHE file")
    parser.add_argument("--njobs"       ,   type=int,     help="Number of jobs to split the LHE file into")
    parser.add_argument("--seed"        ,   type=int,     help="Seed for random number generator")
    parser.add_argument("--pythia8_card",   type=str,     help="Pythia8 card", default="pythia8_card.cmd")
    
    ##### FLAGS #####
    parser.add_argument("--postonly"    ,   action="store_true", help="Only run post-processing")
    parser.add_argument("--overwriteLHE",   action="store_true", help="Overwrite existing files")
    parser.add_argument("--overwriteHEPMC", action="store_true", help="Overwrite existing files")
    parser.add_argument("--overwriteROOT",  action="store_true", help="Overwrite existing files")
    return parser.parse_args()


def file_checking(args):
    '''Check if the input files are present'''
    # **** INPUT FILE, PYTHIA8 CARD ****
    assert os.path.exists(args.input_file), f"Input file {args.input_file} not found"
    assert os.path.exists(args.pythia8_card), f"Pythia8 card {args.pythia8_card} not found"
    
    # **** FLAGS ****
    target_LHE = [f"events_{i}.lhe.gz" for i in range(args.seed, args.seed + args.njobs)]
    target_HEPMC = [f"events_{i}.hepmc" for i in range(args.seed, args.seed + args.njobs)]
    target_ROOT = [f"delphes_output_{i}.root" for i in range(args.seed, args.seed + args.njobs)]
    
    avail_LHE = [f for f in os.listdir() if f.endswith(".lhe.gz")]
    avail_HEPMC = [f for f in os.listdir() if f.endswith(".hepmc")]
    availe_ROOT = [f for f in os.listdir() if f.endswith(".root")]
    
    # *** Check overlapping files ***
    overlap_LHE = set(target_LHE).intersection(set(avail_LHE))
    overlap_HEPMC = set(target_HEPMC).intersection(set(avail_HEPMC))
    overlap_ROOT = set(target_ROOT).intersection(set(availe_ROOT))
    
    if not args.overwriteLHE and len(overlap_LHE) > 0:
        assert set(target_LHE) == overlap_LHE, f"There are overlapping LHE files: {overlap_LHE}, but --overwriteLHE flag is not set"
    
    if not args.overwriteHEPMC and len(overlap_HEPMC) > 0:
        assert set(target_HEPMC) == overlap_HEPMC, f"There are overlapping HEPMC files: {overlap_HEPMC}, but --overwriteHEPMC flag is not set"
    
    if not args.overwriteROOT and len(overlap_ROOT) > 0:
        assert set(target_ROOT) == overlap_ROOT, f"There are overlapping ROOT files: {overlap_ROOT}, but --overwriteROOT flag is not set"

def check_status(logger, args):
    '''Check the status of the jobs'''
    files = os.listdir("run_status")
    n_jobs = args.njobs
    
    n_pythia_running = len([f for f in files if "pythia" in f])
    n_delphes_running = len([f for f in files if "delphes" in f])
    n_completed = len([f for f in files if "completed" in f])

    logger.info(
        f"\tRunning: {n_pythia_running} Pythia\t" + \
        f"{n_delphes_running} Delphes\t" + \
        f"Completed: {n_completed}/{n_jobs}"
    )


def Pythia8Delphes(input_args: dict) -> dict:
    seed = input_args["seed"]
    infile = input_args["input_file"]
    root_outfile = f"delphes_output_{seed}.root"
    hepmc_outfile = f"events_{seed}.hepmc"
    timer = {}

    logger = logging.getLogger(f"P8D_{seed}")
    logger.setLevel(logging.INFO)

    ##### INITIALIZE PYTHIA8 CARD #####
    original_pythia8_card = input_args["pythia8_card"]
    new_pythia8_card = f"pythia8_card_{seed}.cmd"
    
    os.system(f"cp {original_pythia8_card} {new_pythia8_card}")
    os.system(f"sed -i 's/$SEED/{seed}/g' {new_pythia8_card}")
    os.system(f"sed -i 's/$INFILE/{infile}/g' {new_pythia8_card}")
    os.system(f"sed -i 's/$OUTFILE/{hepmc_outfile}/g' {new_pythia8_card}")

    start_time = time.time()
    os.system(f"touch run_status/pythia_{seed}.running")
    check_status(logger, args)
    os.system(f"{MG5_PYTHIA_INTERFACE} {new_pythia8_card} >> log_{seed}.txt")
    os.system(f"rm run_status/pythia_{seed}.running")
    timer["Pythia"] = time.time() - start_time

    start_time = time.time()
    os.system(f"touch run_status/delphes_{seed}.running")
    check_status(logger, args)
    os.system(
        f"{DELPHES_PATH} {DELPEHS_CARD} {root_outfile} events_{seed}.hepmc >> log_{seed}.txt 2>&1" # 2>&1 redirects stderr to stdout (avoid printing progress bar)
    )
    os.system(f"rm run_status/delphes_{seed}.running")
    timer["Delphes"] = time.time() - start_time

    status = 0 if os.path.exists(root_outfile) else 1

    ##### POST-PROCESSING #####
    os.system(f"rm {new_pythia8_card}")
    os.system(f"mv log_{seed}.txt logs/")
    os.system(f"touch run_status/{seed}.completed")
    
    check_status(logger, args)

    return {"status": status, "time": timer}


def get_tot_evts():
    files = os.listdir("./")
    files = [f for f in files if f.endswith(".root")]

    tot_evts = 0
    for f in files:
        f = ROOT.TFile(f"{f}")
        tree = f.Get("Delphes")
        tot_evts += tree.GetBranch("Event_size").GetEntries()

    return tot_evts


def format_time(s) -> str:
    days, remainder  = divmod(s, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    return_str = f"{seconds:.2f}s"
    
    if minutes > 0: return_str = f"{minutes:.0f}m " + return_str
    if hours > 0:   return_str = f"{hours:.0f}h " + return_str
    if days > 0:    return_str = f"{days:.0f}d " + return_str
    
    return return_str


def get_cross_section(text):
    # Get the cross-section from the log file: target format: "Inclusive cross section: 2.17868639e-05  +-  3.08637916e-07 mb"
    lines = text.split("\n")
    for line in lines:
        if "Inclusive cross section" in line:
            num = line.split("cross section: ")[1]
            cross_section, unc = num.split("  +-  ")
            unc = unc.split()[0]
            return float(cross_section) * 1e9, float(unc) * 1e9
    return None, None


def get_gained_events(text):
    lines = text.split("\n")
    for line in lines:
        if "| sum" in line:
            num = line.split("|")[-3]
            spt_num = num.split()
            accepted, gained = spt_num[0], spt_num[-1]
            return int(accepted), int(gained)
    return None, None


def phys_summary(logger, args):
    # Print the summary of each seed in format: Seed: , Requested events: , Gained events: , Efficiency: ,cross-section:, uncertainty:
    ljust_space = 20
    seeds = np.arange(args.seed, args.seed + args.njobs)
    # print header
    print("\n\n")
    logger.info("=================== Physics Summary ===================")
    text = (
        "Seed".ljust(ljust_space)
        + "Requested events".ljust(ljust_space)
        + "Gained events".ljust(ljust_space)
    )
    text += (
        "Efficiency".ljust(ljust_space)
        + "Cross-section (pb)".ljust(ljust_space)
        + "Uncertainty (pb)".ljust(ljust_space)
    )
    logger.info(text)
    for seed in seeds:
        # Read from info.csv
        info = pd.read_csv("info.csv")
        try:
            requested_evts = info[info["seed"] == seed]["nevents"].values[0]
        except:
            print(f"Seed {seed} not found in info.csv, skipping...")
            continue

        with open(f"logs/log_{seed}.txt", "r") as f:
            log_text = f.read()
        cross_section, unc = get_cross_section(log_text)

        accepted, gained = get_gained_events(log_text)
        efficiency = gained / requested_evts

        text = (
            f"{seed}".ljust(ljust_space)
            + f"{requested_evts}".ljust(ljust_space)
            + f"{gained}".ljust(ljust_space)
        )
        text += (
            f"{efficiency:.4f}".ljust(ljust_space)
            + f"{cross_section:.4f}".ljust(ljust_space)
            + f"{unc:.4f}".ljust(ljust_space)
        )
        logger.info(text)


def main(args):
    file_checking(args)
    start_time = time.time()

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger("Main")

    os.system("mkdir logs")
    os.system("mkdir run_status")
    
    logger.info(f"Splitting {args.input_file} into {args.njobs} sub-files...")
    nevents = split_lhe(args.input_file, args.njobs)
    sub_nevents = nevents // args.njobs

    seeds = list(range(args.seed, args.seed + args.njobs))
    
    prefix_name = args.input_file.split(".lhe")[0]
    for i, seed in enumerate(seeds): os.rename(f"{prefix_name}_{i}.lhe.gz", f"events_{seed}.lhe.gz")

    input_args = []
    for seed in seeds:
        input_args.append(
            {
                "seed": seed,
                "input_file": f"events_{seed}.lhe.gz",
                "pythia8_card": args.pythia8_card,
                "njobs": args.njobs,
            }
        )

    logger.info(
        f"Split {args.input_file} with total {nevents} events into {args.njobs} sub-files with {sub_nevents} events each"
    )

    ##### SUBMIT JOBS #####
    logger.info(f"Submitting {args.njobs} jobs...")
    with Pool(args.njobs) as p:
        results = p.map(Pythia8Delphes, input_args)

    ##### CREATE INFO #####
    info = pd.DataFrame(
        {
            "seed": seeds,
            "nevents": [sub_nevents] * args.njobs,
            "status": [r["status"] for r in results],
            "out_file": [f"delphes_output_{seed}.root" for seed in seeds],
            "pythia": [r["time"]["Pythia"] for r in results],
            "delphes": [r["time"]["Delphes"] for r in results],
            "total": [r["time"]["Pythia"] + r["time"]["Delphes"] for r in results],
        }
    )
    info.to_csv("info.csv", index=False)
    
    for job in results:
        if job["status"] != 0: logger.error(f"Seed {job['seed']} failed. Check logs/log_{job['seed']}.txt for more information.")
    
    os.system("rm -r run_status")
    phys_summary(logger, args)
    
    ##### SUMMARY #####
    runtime = time.time() - start_time
    tot_pythia_time = sum(info["pythia"])
    tot_delphes_time = sum(info["delphes"])
    total_time = tot_pythia_time + tot_delphes_time
    
    
    logger.info(f"Total time taken: {format_time(runtime)}")
    logger.info(f"Time taken for each step:")
    logger.info(f"\tPythia8: {format_time(tot_pythia_time)} ({100 * tot_pythia_time / total_time:.2f}%)")
    logger.info(f"\tDelphes: {format_time(tot_delphes_time)} ({100 * tot_delphes_time / total_time:.2f}%)")
    
    # nevent/sec
    tot_evts = get_tot_evts()
    logger.info(f"Total events generated: {tot_evts}/{nevents}")
    logger.info(f"Efficiency: {tot_evts / nevents:.2f}")
    logger.info(f"Events generated per second: {tot_evts / runtime:.2f}")


if __name__ == "__main__":
    args = get_args()
    if args.postonly: phys_summary(logging.getLogger("Main"), args)
    else: main(args)
