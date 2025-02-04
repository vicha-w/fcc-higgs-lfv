#!/work/app/modules/software/Python/3.9.6-GCCcore-11.2.0/bin/python

import sys
import os
import argparse
import time
from processes import BACKGROUND, MUTAU_SIGNAL, ETAU_SIGNAL, PURPLE_ETAU_E, PURPLE_MUTAU_E


# Universal script to submit jobs, monitoring jobs, and post-processing

ALL_PROCESSES = {**BACKGROUND, **MUTAU_SIGNAL, **ETAU_SIGNAL, **PURPLE_ETAU_E, **PURPLE_MUTAU_E}


slurm_configs = {
    "--qos": "cu_htc",
    "--partition": "cpu",
    "--job-name": None,  # processes_name
    "--output": "sbatch.log",
    "--nodes": 1,
    "--ntasks": 1,
    "--cpus-per-task": None,
    "--mem": None,
}

envs_configs = {
    "module purge": "",
    "source /work/app/share_env/hepsw-gcc11p2-py3p9p9.sh": "",
}


def argpass():
    parser = argparse.ArgumentParser(
        description="Submit jobs, monitor jobs, and post-process the output"
    )
    parser.add_argument(
        "--mode", type=str, default="job_submit", help="Mode to run the script"
    )
    parser.add_argument(
        "--process", type=str, default=None, help="Process to run the script"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite the output directory", default=False
    )
    parser.add_argument("--outdir", type=str, default="auto", help="Output directory")
    parser.add_argument(
        "--postonly",
        action="store_true",
        help="Run post-processing only",
        default=False,
    )

    # job_submit ONLY
    parser.add_argument(
        "--ncpus", type=int, default=-1, help="Number of CPUs to use"
    )  # -1 = auto
    parser.add_argument(
        "--max_auto_cpus", type=int, default=20, help="Max number of CPUs to use"
    )
    
    # Experimental feature
    # extracted file path
    parser.add_argument(
        "--minimal", action="store_true", help="Use minimal file", default=False
    )
    
    args = parser.parse_args()

    return args


def get_files(args) -> list:

    if args.minimal: path = ALL_PROCESSES[args.process]["minimal"]
    else: path = ALL_PROCESSES[args.process]["path"]
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


def job_submit(args):
    def validate_args(args):
        
        assert args.process in ALL_PROCESSES.keys(), f"Process {args.process} is not available"

        outdir = args.outdir
        if outdir == "auto":
            outdir = f"outputs/{args.process}"
        
        if args.minimal: outdir += "_minimal"

        if os.path.exists(outdir):
            if args.overwrite:
                print(f"Overwriting the directory {outdir}")
                os.system(f"rm -rf {outdir}")
                
            else:
                print(f"Directory {outdir} already exists")
                sys.exit(1)

        args.outdir = outdir
        files = get_files(args)

        if args.ncpus == -1:
            args.ncpus = len(files) + 1  # +1 for monitoring script
            args.ncpus = min(
                args.ncpus, max(args.max_auto_cpus, 1)
            )  # Handle the case of 0 files

        return args

    def construct_slurm_script(args):
        # Construct SLURM script
        slurm_configs["--job-name"] = args.process
        slurm_configs["--cpus-per-task"] = args.ncpus
        slurm_configs["--mem"] = f"{args.ncpus * 4}G"

        slurm_script = "#!/bin/bash\n\n\n"  # Header
        for key, value in slurm_configs.items():
            slurm_script += f"#SBATCH {key}={value}\n"

        slurm_script += "\n"
        for key, value in envs_configs.items():
            slurm_script += f"{key}\n"

        slurm_script += "\n"
        slurm_script += (
            f"python -u pyinterface.py --mode job_monitor --process {args.process}"
        )
        
        if args.minimal: slurm_script += " --minimal"

        return slurm_script

    def setup_dir(args):
        # Create the output directory
        outdir = args.outdir
        os.makedirs(outdir)

        script_files = ["Delphes.C", "Delphes.h", "read-fcc-higgs-v2.cpp"]
        for script_file in script_files:
            os.system(f"cp {script_file} {outdir}")

        os.system(f"cp pyinterface.py {outdir}")
        os.system(f"cp processes.py {outdir}")

        # Write the SLURM script
        slurm_script = construct_slurm_script(args)
        with open(f"{outdir}/slurm.sh", "w") as f:
            f.write(slurm_script)

        return outdir

    """
        Setup the job submission
        - Create the directory: outputs/{processes_alias}
        - Copy the script file (for later validation)
            1. Delphes.C
            2. Delphes.h
            3. read-fcc-higgs-v2.cpp
            4. processes.py
        
        - Create the SLURM script
        - Copy this file (pyinterface.py) to the outputs/{processes_alias}
        # These step ensure the reproducibility of the results, version control of the script
    """

    validate_args(args)
    setup_dir(args)



def post_process(args):
    """
    - Real all .root files in current dir
    - Each root, open and get the histogram:
        - 'mutau_e_{catagory}_{jetgroup}' where jet = 0, 1, category = highmass, lowmass
        - 'etau_mu_{catagory}_{jetgroup}' where jet = 0, 1, category = highmass, lowmass
    - Perform counting every file and merge the counting into single .root

    """
    import copy
    # 30 ab^-1 / (N/ cross-section) = scale factor
    def compute_scale(cross_section, N, target_lumi):
        if N == -1:
            return None
        return target_lumi / (N / cross_section)

    # Mapping all files
    files = [f for f in os.listdir() if (f.endswith(".root") and f != "merged.root")]
    print(f"Files: {files}")

    # exclude the merged file
    import ROOT

    hist_dicts = {}
    for f in files:
        file = ROOT.TFile(f)

        keys = file.GetListOfKeys()
        for key in keys:
            hist = key.ReadObj()
            name = hist.GetName()
            if name not in hist_dicts:
                hist_dicts[name] = copy.deepcopy(hist)
            else:
                hist_dicts[name].Add(hist)

    # Print entries of merged histograms
    print("Merge results")
    ljust_space = 20
    for key, hist in hist_dicts.items():
        # print(f"histogram: {key}, Entries: {hist.GetEntries()}")
        print(f"histogram: {key.ljust(ljust_space)}, Entries: {int(hist.GetEntries())}")

    output = ROOT.TFile("merged.root", "RECREATE")
    for key, hist in hist_dicts.items():
        hist.Write()



    # ================== Weighting ==================
    tot_evts = hist_dicts["mutau_e_step00"].GetEntries()
    print(f"\nComputing scale factor with total events: {tot_evts}")
    process = args.process
    if process is None:
        import pandas as pd
        print("Process is not defined, reading from info.csv")
        df = pd.read_csv("info.csv")
        process = df["process"].unique()[0]
        print(f"Process: {process}")
        
    cross_section = ALL_PROCESSES[process]["cross-section"]
    target_lumi = 30
    target_lumi_pb = target_lumi * 1000
    scale = compute_scale(cross_section, tot_evts, target_lumi_pb) # 30 ab^-1
    print(f"\tTotal events: {int(tot_evts)}, \n\tCross-section: {cross_section} pb, \n\tTarget lumi: {target_lumi_pb} pb^-1, \n\tScale: {scale}")
    
    weighted_hist_dicts = {}
    for key, hist in hist_dicts.items():
        weighted_hist = hist.Clone(f"weighted_{key}")
        weighted_hist.Scale(scale)
        weighted_hist_dicts[key] = weighted_hist

    for key, hist in weighted_hist_dicts.items():
        hist.Write()

    output.Close()

def run_cut(file, out_file):
    command = (
        f'root -l -b -q "read-fcc-higgs-v2.cpp(\\"{file}\\", \\"{out_file}\\")"'
    )
    start_time = time.time()
    os.system(command)
    end_time = time.time()
    time_taken = end_time - start_time
    
    # Experiment, open info.csv and get the line with column file == file
    # Update the time_taken column with the time_taken
    import pandas as pd
    df = pd.read_csv("info.csv")
    df.loc[df["file"] == file, "time_taken"] = time_taken
    df.to_csv("info.csv", index=False)
    
    return {"file": file, "time_taken": time_taken}


def job_monitor(args):
    from multiprocessing import Pool
    import pandas as pd
    import ROOT

    def convert_time(s):
        return time.strftime("%H:%M:%S", time.gmtime(s))

    # Running from the output directory
    def get_cpu_usage():
        # Read slurm.sh
        # extract line with #SBATCH --cpus-per-task

        with open("slurm.sh", "r") as f:
            for line in f:
                if "--cpus-per-task" in line:
                    return int(line.split("=")[-1])

    def get_tot_evts():
        file = ROOT.TFile("merged.root")
        hist = file.Get("mutau_e_step00")
        return hist.GetEntries()

    def get_tot_passed():
        types = ["mutau_e", "etau_mu"]
        categories = ["highmass", "lowmass"]
        jetgroups = ["0j", "1j"]

        key_dict = {}

        file = ROOT.TFile("merged.root")
        for t in types:
            for c in categories:
                for j in jetgroups:
                    key = f"{t}_{c}_{j}"
                    hist = file.Get(key)
                    key_dict[key] = hist.GetEntries()

        return key_dict

    # Availiable cpus
    ncpus = get_cpu_usage()
    print(f"Used CPUs: {ncpus}")

    process = args.process
    files = get_files(args)
    print(f"List of files ({len(files)}):")
    for f in files: print(f"\t- {f}")

    pre_df = {
        "file": files,
        "out_file": [f"{process}_{n}.root" for n in range(len(files))],
        "process": [process] * len(files),
        "time_taken": [0] * len(files),
    }
    df = pd.DataFrame(pre_df)
    df.to_csv("info.csv", index=False)

    start_time = time.time()
    with Pool(ncpus) as p:
        out_dict = p.starmap(run_cut, df[["file", "out_file"]].values)

        for out in out_dict:
            file = out["file"]
            df.loc[df["file"] == file, "time_taken"] = out["time_taken"]

        df.to_csv("info.csv", index=False)
    tot_time = time.time() - start_time
    print("Done")

    post_process(args)
    # Summary of the job, time taken, total events, events/s
    tot_evts = get_tot_evts()
    tot_passed = get_tot_passed()
    print(f"\n\nTotal events: {tot_evts}")
    print("Time taken")
    print(f"\tSingle core, sum: {df['time_taken'].sum():.2f} s\t with {tot_evts / df['time_taken'].sum():.2f} events/s")
    print(f"\tMulti core, sum: {tot_time:.2f} s\t with {tot_evts / tot_time:.2f} events/s")
    print(f"\tEfficiency: {df['time_taken'].sum() / tot_time:.2f}")
    print("Total passed events:")
    for key, value in tot_passed.items():
        print(f"\t- {key}: {value}")


if __name__ == "__main__":
    args = argpass()
    if args.postonly:
        post_process(args)
    else:
        if args.mode == "job_submit":
            job_submit(args)
        elif args.mode == "job_monitor":
            job_monitor(args)
        elif args.mode == "post_process":
            post_process(args)
        else:
            print("Invalid mode")
            sys.exit(1)

