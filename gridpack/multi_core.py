#!/work/app/modules/software/Python/3.9.6-GCCcore-11.2.0/bin/python

import os
import argparse
import time
import ROOT
import pandas as pd
import numpy as np
import logging
from multiprocessing import Pool


# Script to parallelize gridpack from Madgraph -> Pythia -> Delphes -> ROOT
# Input: nevents(total), seed (start number), njobs, PY8 card


# ********** PATHS **********
MG5_PYTHIA_INTERFACE = "/work/app/pythia8/MGInterface/1.3/MG5aMC_PY8_interface"
DELPHES_PATH = "/work/home/psriling/MsProject/Delphes-3.5.0/DelphesHepMC2"
DELPEHS_CARD = "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/FCChh_I.tcl"

ROOT.gInterpreter.ProcessLine("gErrorIgnoreLevel = kFatal;")
os.system(
    "export ROOT_INCLUDE_PATH=/work/home/psriling/MsProject/Delphes-3.5.0/external:$ROOT_INCLUDE_PATH"
)


def get_args():
    parser = argparse.ArgumentParser(
        description="Run gridpack from Madgraph -> Pythia -> Delphes -> ROOT"
    )
    parser.add_argument(
        "--nevents", type=int, default=1000, help="Number of events to generate"
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="Seed for random number generator"
    )
    parser.add_argument(
        "--njobs", type=int, default=1, help="Number of jobs to run in parallel"
    )
    parser.add_argument(
        "--noclean", action="store_true", help="Do not clean the intermediate files"
    )
    parser.add_argument("--verbose", action="store_true", help="Print verbose output")

    parser.add_argument(
        "--mg5parall", action="store_true", help="Run Madgraph in parallel"
    )

    # optional arguments: pythia_card, output_path
    parser.add_argument(
        "--pythia_card",
        type=str,
        default="pythia8_card.cmd",
        help="Path to the pythia card",
    )
    
    parser.add_argument(
        "--postonly", action="store_true", help="Run only post-processing (read cross-section)"
    )

    return parser.parse_args()


def run_gridpack(lhe_in, nevents, seed, pythia_card, out_file, args):
    runtime = {
        "Madgraph": 0,
        "Pythia": 0,
        "Delphes": 0,
    }

    # Create logger
    logging.basicConfig(level=logging.DEBUG)
    # Set format of logger to [TIME] [LEVEL] [logger] [MESSAGE]
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] [%(name)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    sub_logger = logging.getLogger(f"Subrun_S{seed} ")
    n_limit = 1
    if args.mg5parall:
        
        # Limiter, only initialize n directories at a time to avoid overloading in disk
        # Checking if there are directory of previous seeds, mg5_(seed-1)/madevent exists
        while not os.path.exists(f"mg5_{seed-n_limit}/run.sh"):
            
            if seed < args.seed + n_limit:
                break
            time.sleep(5)
        sub_logger.info(f"Initalizing Madgraph for seed {seed}")
        os.system(f"mkdir mg5_{seed}")
        os.system(f"cp madevent.tar.gz mg5_{seed}/")
        os.chdir(f"mg5_{seed}")
        os.system("tar -xzvf madevent.tar.gz > /dev/null")
        os.system("cp ../run.sh .")
        # Run the Madgraph
        sub_logger.info(f"Running Madgraph for seed {seed}")
        st = time.time()
        os.system(f"./run.sh {nevents} {seed} > log_{seed}.txt")
        # Exit the directory & clean up
        os.chdir("..")
        os.system(f"mv mg5_{seed}/events_{seed}.lhe.gz .")
        os.system(f"mv mg5_{seed}/log_{seed}.txt .")
        os.system(f"rm -r mg5_{seed}")
        runtime["Madgraph"] = time.time() - st

    # Preprocess
    pythia_card_seed = f"pythia_card_{seed}.cmd"
    with open(pythia_card, "r") as f:
        pythia_card_content = f.read()
    pythia_card_content = pythia_card_content.replace("$SEED", str(seed))
    pythia_card_content = pythia_card_content.replace("$INFILE", lhe_in)
    pythia_card_content = pythia_card_content.replace(
        "$OUTFILE", f"events_{seed}.hepmc"
    )

    with open(pythia_card_seed, "w") as f:
        f.write(pythia_card_content)

    # Pythia
    st = time.time()
    sub_logger.info(f"Running Pythia for seed {seed}")
    os.system(f"{MG5_PYTHIA_INTERFACE} {pythia_card_seed} >> log_{seed}.txt")
    runtime["Pythia"] = time.time() - st
    # Delphes
    st = time.time()
    sub_logger.info(f"Running Delphes for seed {seed}")
    os.system(
        f"{DELPHES_PATH} {DELPEHS_CARD} {out_file} events_{seed}.hepmc >> log_{seed}.txt"
    )
    runtime["Delphes"] = time.time() - st

    return {
        "tot_time": sum([v for k, v in runtime.items()]),
        "madgraph": runtime["Madgraph"],
        "pythia": runtime["Pythia"],
        "delphes": runtime["Delphes"],
        "seed": seed,
        "nevents": nevents,
        "out_file": out_file,
    }


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
    # Convert seconds to human readable format
    days, remainder = divmod(s, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, seconds = divmod(remainder, 60)
    # Construct string
    return_str = f"{seconds:.2f}s"
    if minutes > 0: return_str = f"{minutes:.0f}m " + return_str
    if hours > 0: return_str = f"{hours:.0f}h " + return_str
    if days > 0: return_str = f"{days:.0f}d " + return_str
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
    # Format:  | sum                                                |       50000      50000       4983 |   2.179e-05  2.762e-07 |
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
    text = "Seed".ljust(ljust_space) + "Requested events".ljust(ljust_space) + "Gained events".ljust(ljust_space)
    text += "Efficiency".ljust(ljust_space) + "Cross-section (pb)".ljust(ljust_space) + "Uncertainty (pb)".ljust(ljust_space)
    logger.info(text)
    for seed in seeds:
        requested_evts = args.nevents // args.njobs
        # Get the text in log file from logs/log_{seed}.txt
        with open(f"logs/log_{seed}.txt", "r") as f:
            log_text = f.read()
        # Get the cross-section from the log file: target format: "Inclusive cross section: 2.17868639e-05  +-  3.08637916e-07 mb"
        cross_section, unc = get_cross_section(log_text)
        # Get the gained events
        accepted, gained = get_gained_events(log_text)
        efficiency = gained / requested_evts
        
        text = f"{seed}".ljust(ljust_space) + f"{requested_evts}".ljust(ljust_space) + f"{gained}".ljust(ljust_space)
        text += f"{efficiency:.4f}".ljust(ljust_space) + f"{cross_section:.4f}".ljust(ljust_space) + f"{unc:.4f}".ljust(ljust_space)
        logger.info(text)
    
    
        
        


def parallen_run(args):
    # Create gridspace
    runtime = {
        "Madgraph": 0.1,
        "All": 0.1,
    }

    # Create logger
    logging.basicConfig(level=logging.DEBUG)
    
    logger = logging.getLogger("Main ")
    if not args.postonly:
        all_start = time.time()
        seeds = np.arange(args.seed, args.seed + args.njobs)
        nevents_per_job = args.nevents // args.njobs
        out_files = [f"delphes_output_{seed}.root" for seed in seeds]

        out_files = [os.path.abspath(f) for f in out_files]
        lhe_files = [f"events_{seed}.lhe.gz" for seed in seeds]

        info = pd.DataFrame()
        info["seed"] = seeds
        info["nevents"] = nevents_per_job
        info["out_file"] = [f"delphes_output_{seed}.root" for seed in seeds]
        info["runtime"] = 0.0
        info["madgraph"] = 0.0
        info["pythia"] = 0.0
        info["delphes"] = 0.0
        info.to_csv("info.csv", index=False)


        # Run Madgraph in serial
        if not args.mg5parall:
            st_mg = time.time()
            for idx, seed in enumerate(seeds):
                logger.info(f"Running Madgraph for seed {seed} [{idx + 1}/{len(seeds)}]")
                start_time = time.time()
                os.system(f"./run.sh {nevents_per_job} {seed} > log_{seed}.txt")
                logger.info(
                    f"\t- Runtime: {time.time() - start_time:.2f} s, {nevents_per_job / (time.time() - start_time):.2f} evts/s"
                )
            runtime["Madgraph"] = time.time() - st_mg
        else:
            os.system("tar -czvf madevent.tar.gz madevent > /dev/null")

        arg_list = [
            (lhe_file, nevents_per_job, seed, args.pythia_card, out_file, args)
            for lhe_file, seed, out_file in zip(lhe_files, seeds, out_files)
        ]

        # Run Pythia and Delphes in parallel
        with Pool(args.njobs) as p:
            out_dicts = p.starmap(run_gridpack, arg_list)

        info = pd.read_csv("info.csv")
        info["runtime"] = [d["tot_time"] for d in out_dicts]
        info["madgraph"] = [d["madgraph"] for d in out_dicts]
        info["pythia"] = [d["pythia"] for d in out_dicts]
        info["delphes"] = [d["delphes"] for d in out_dicts]
        info.to_csv("info.csv", index=False)


        # Clean up
        if not args.noclean:
            for seed in seeds:
                os.system(f"rm events_{seed}.lhe.gz")
                os.system(f"rm events_{seed}.hepmc")
                os.system(f"rm pythia_card_{seed}.cmd")
            os.system("mkdir logs")
            os.system("mv log_* logs/")

        if args.mg5parall:
            # Remove tar
            os.system("rm -r madevent.tar.gz")

        runtime["All"] = time.time() - all_start
    
    info = pd.read_csv("info.csv")
    
    tot_evts = get_tot_evts()

    sum_all_runtime = sum(info["runtime"])

    phys_summary(logger, args)
    print("\n")
    logger.info("=================== Performance Report ===================")

    logger.info(f"Expected events: {args.nevents}")
    logger.info(f"Gained events: {tot_evts}")
    logger.info(f"Efficiency: {tot_evts / args.nevents:.4f}")
    logger.info(f"Runtime: {format_time(runtime['All'])}")
    if not args.mg5parall:
        logger.info(f"\t- Madgraph Tot. {format_time(sum(info['madgraph']))}")
        groups = ["pythia", "delphes"]
    else:
        groups = ["madgraph", "pythia", "delphes"]

    ljust_space = 30
    for group in ["madgraph", "pythia", "delphes"]:
        pct = sum(info[group]) / sum_all_runtime * 100
        logger.info(
            # f"\t- {group.capitalize()} total, single-core: {format_time(sum(info[group]))} ({pct:.2f}%)"
            f"\t- {group.capitalize()} total, single-core:".ljust(ljust_space+10) + f"{format_time(sum(info[group]))}".ljust(ljust_space) + f"({pct:.2f}%)"
        )

    # Single and multi core
    print("\n")
    logger.info(
        # f"\t- Overall, single-core: {format_time(sum(info['runtime']))}, {tot_evts / sum(info['runtime']):.2f} evts/s"
        "\t- Overall, single-core:".ljust(ljust_space+10) + f"{format_time(sum(info['runtime']))}".ljust(ljust_space) + f"{tot_evts / sum(info['runtime']):.2f} evts/s"
    )
    logger.info(
        # f"\t- Overall, multi-cores: {format_time(runtime['All'])}, {tot_evts / runtime['All']:.2f} evts/s"
        "\t- Overall, multi-cores:".ljust(ljust_space+10) + f"{format_time(runtime['All'])}".ljust(ljust_space) + f"{tot_evts / runtime['All']:.2f} evts/s"
    )


if __name__ == "__main__":
    args = get_args()
    if args.seed is None:
        args.seed = np.random.randint(0, 1000)
    parallen_run(args)
