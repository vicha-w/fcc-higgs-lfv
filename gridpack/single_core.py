#!/work/app/modules/software/Python/3.9.6-GCCcore-11.2.0/bin/python

# Script to run gridpack from Madgraph -> Pythia -> Delphes -> ROOT
# Input to this script: nevents, seed, PY8 card, output path
import os
import argparse
import time


# ********** PATHS **********
MG5_PYTHIA_INTERFACE = "/work/app/pythia8/MGInterface/1.3/MG5aMC_PY8_interface"
DELPHES_PATH = "/work/home/psriling/MsProject/Delphes-3.5.0/DelphesHepMC2"
DELPEHS_CARD = "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/FCChh_I.tcl"


# ********** ARGUMENTS **********
def get_args():
    parser = argparse.ArgumentParser(
        description="Run gridpack from Madgraph -> Pythia -> Delphes -> ROOT"
    )
    parser.add_argument(
        "--nevents", type=int, default=1000, help="Number of events to generate"
    )
    parser.add_argument(
        "--seed", type=int, default=12345, help="Seed for random number generator"
    )
    parser.add_argument(
        "--noclean", action="store_true", help="Do not clean the intermediate files"
    )
    parser.add_argument("--verbose", action="store_true", help="Print verbose output")

    # optional arguments: pythia_card, output_path
    parser.add_argument(
        "--pythia_card",
        type=str,
        default="pythia8_card.cmd",
        help="Path to the pythia card",
    )
    parser.add_argument(
        "--output_path",
        type=str,
        default="delphes_output.root",
        help="Path to save the output ROOT file",
    )

    return parser.parse_args()



def run_gridpack(args):
    nevents = args.nevents
    seed = args.seed
    pythia_card = args.pythia_card
    output_path = args.output_path
    runtime = {
        "Madgraph": {"start": None, "end": None},
        "Pythia": {"start": None, "end": None},
        "Delphes": {"start": None, "end": None},
    }
    # Preprocess
    pythia_card_seed = f"pythia_card_{seed}.cmd"
    with open(pythia_card, "r") as f:
        pythia_card_content = f.read()
    pythia_card_content = pythia_card_content.replace("$SEED", str(seed))
    pythia_card_content = pythia_card_content.replace(
        "$INFILE", f"events_{seed}.lhe.gz"
    )
    pythia_card_content = pythia_card_content.replace(
        "$OUTFILE", f"events_{seed}.hepmc"
    )

    with open(pythia_card_seed, "w") as f:
        f.write(pythia_card_content)

    # Madgraph
    runtime["Madgraph"]["start"] = time.time()
    if args.verbose:
        print(f"Running Madgraph with {nevents} events and seed {seed}")
        os.system(f"./run.sh {nevents} {seed} | tee log_{seed}.txt")
    else:
        os.system(f"./run.sh {nevents} {seed} > log_{seed}.txt")
    runtime["Madgraph"]["end"] = time.time()

    # Pythia
    runtime["Pythia"]["start"] = time.time()
    if args.verbose:
        print("Running Pythia")
        os.system(f"{MG5_PYTHIA_INTERFACE} {pythia_card_seed} | tee -a log_{seed}.txt")
    else:
        os.system(f"{MG5_PYTHIA_INTERFACE} {pythia_card_seed} >> log_{seed}.txt")
    runtime["Pythia"]["end"] = time.time()

    # Delphes
    runtime["Delphes"]["start"] = time.time()
    if args.verbose:
        print("Running Delphes")
        os.system(
            f"{DELPHES_PATH} {DELPEHS_CARD} {output_path} events_{seed}.hepmc | tee -a log_{seed}.txt"
        )
    else:
        os.system(
            f"{DELPHES_PATH} {DELPEHS_CARD} {output_path} events_{seed}.hepmc >> log_{seed}.txt"
        )
    runtime["Delphes"]["end"] = time.time()

    if not args.noclean:
        os.system(f"rm {pythia_card_seed} events_{seed}.lhe.gz events_{seed}.hepmc")

    if args.verbose:
        print(f"Output saved at {output_path}")
        print("Runtime summary:")
        for step, times in runtime.items():
            print(
                f"{step}: {times['end'] - times['start']:.2f} s, {nevents / (times['end'] - times['start']):.2f} nevents/s"
            )

        # Total runtime
        total_runtime = sum(
            [times["end"] - times["start"] for times in runtime.values()]
        )
        print(
            f"Total runtime: {total_runtime:.2f} s, {nevents / total_runtime:.2f} nevents/s"
        )


if __name__ == "__main__":
    args = get_args()
    run_gridpack(args)
