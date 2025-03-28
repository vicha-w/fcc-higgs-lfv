#!/work/app/modules/software/Python/3.9.6-GCCcore-11.2.0/bin/python

import os
import sys
import argparse
import time
import ROOT
import pandas as pd
import random
import numpy as np
import logging
from multiprocessing import Pool

from classes import FullChainRunner, LogHandler, format_time

ROOT.gInterpreter.ProcessLine("gErrorIgnoreLevel = kFatal;") # suppress ROOT TClass warnings

# ********** PATHS **********
MG5_PYTHIA_INTERFACE    = "/work/app/pythia8/MGInterface/1.3/MG5aMC_PY8_interface"
DELPHES_PATH            = "/work/home/psriling/MsProject/Delphes-3.5.0/DelphesHepMC2"
DELPEHS_CARD            = "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/FCChh_I.tcl"



def get_args():
    parser = argparse.ArgumentParser(description="Run gridpack from Madgraph -> Pythia -> Delphes -> ROOT")
    
    # ********** COMMON ARGUMENTS **********
    parser.add_argument("--seed",                   type=int,   default=None,               help="Random seed")
    parser.add_argument("--n_threads",              type=int,   default=1,                  help="Number of threads to use")
    parser.add_argument("--delphes_path",           type=str,   default=DELPHES_PATH,       help="Path to Delphes")
    parser.add_argument("--delphes_card",           type=str,   default=DELPEHS_CARD,       help="Path to the Delphes card")
    parser.add_argument("--pythia_card",            type=str,   default="pythia8_card.cmd", help="Path to the pythia card",)
    parser.add_argument("--mg5_pythia_interface",   type=str,   default=MG5_PYTHIA_INTERFACE,help="Path to the MG5-Pythia interface")
    parser.add_argument("--keep_lhe",               action="store_true",                    help="Keep splited LHE files (use_lhe mode) or keep LHE file (run_gridpack mode)")
    # ********** CASE I: RUN GRIDPACK **********
    parser.add_argument("--run_gridpack",           action="store_true",                    help="Run gridpack")
    parser.add_argument("--gridpack_path",          type=str,   default="gridpacks",        help="Path to the gridpacks")
    parser.add_argument("--n_events",               type=int,   default=1000,               help="Number of events to generate")
    # ********** CASE II: GIVEN LHE **********
    parser.add_argument("--use_lhe",                action="store_true",                    help="Use LHE file")
    parser.add_argument("--lhe_path",               type=str,   default="lhe",              help="Path to the LHE file")

    # ********** FEATURES TO ADD **********
    parser.add_argument("--example_inputs",         action="store_true",                    help="Show example inputs (both gridpack and LHE cases)")
    return parser.parse_args()


def print_example_inputs():
    out_str = '''
    Example inputs:
    Gridpack:
    python fullchain_interface.py \\
        --run_gridpack \\
        --gridpack_path GRIDPACK_PATH \\
        --seed SEED \\
        --n_events NEVENTS \\
        --n_threads NTHREADS
        # Optional: --keep_lhe
        
    LHE:
    python fullchain_interface.py \\
        --use_lhe \\
        --lhe_path LHE_PATH \\
        --seed SEED \\
        --n_threads NTHREADS
        # Optional: --keep_lhe
    '''
    print(out_str)
    # Exit the program
    sys.exit(0)

def run_job(args): # FOR ASSIGN TO MULTIPROCESSING
    runner = FullChainRunner(**args)
    runner.run_full_chain()
    return runner





class fullChainInterface:
    def __init__(self):
        self.args = get_args()
        
        if self.args.example_inputs: print_example_inputs()
        
        self.initialize_logging()
        
        if self.args.run_gridpack and self.args.use_lhe:
            raise ValueError("Cannot run gridpack and use LHE at the same time")
        if not self.args.run_gridpack and not self.args.use_lhe:
            raise ValueError("Must run gridpack or use LHE")
        
        # ********** SEEDS **********
        if self.args.seed is None: self.args.seed = random.randint(1,1000000)
        self.seeds = [self.args.seed + i for i in range(self.args.n_threads)]
        
        # ********** PRINTOUT CONFIG **********
        self.printout_config()
        
        # ********** INITIALIZATION **********
        if self.args.use_lhe                : self.init_lhe_chain()
        if self.args.run_gridpack           : self.init_gridpack_chain()
        if not os.path.exists("info.csv")   : self.init_info_csv()
        if not os.path.exists("run_status") : os.system("mkdir run_status")
        else                                : os.system("rm -r run_status/*")
        
        self.fullChainArgs = [self.get_fullchain_args(seed) for seed in self.seeds]


    def printout_config(self):
        # ********** PRINTOUT CONFIG **********
        out_str = '\n'
        out_str +=  f"{'='*20} CONFIG {'='*20}\n" + \
                    f"Number of threads: {self.args.n_threads}\n" + \
                    f"Seeds: {self.seeds}\n"
                    
        if self.args.run_gridpack:
            out_str += "Running mode: Gridpack\n" + \
                        f"Gridpack path: {self.args.gridpack_path}\n" + \
                        f"Requested events: {self.args.n_events}/ {self.args.n_events // self.args.n_threads} per thread\n"
        if self.args.use_lhe:
            nevents = os.popen(f"zgrep -c \"<event>\" {self.args.lhe_path}").read().strip()
            out_str += "Running mode: LHE\n" + \
                        f"LHE path: {self.args.lhe_path}\n" + \
                        f"Number of events: {nevents}/ {nevents // self.args.n_threads} per thread\n"
        # Paths        
        out_str +=  "\nPATHS:\n" + \
                    f"Delphes path: {self.args.delphes_path}\n" + \
                    f"Delphes card: {self.args.delphes_card}\n" + \
                    f"MG5-Pythia interface: {self.args.mg5_pythia_interface}\n" + \
                    f"Pythia card: {self.args.pythia_card}\n" 
        # Optional
        out_str +=  "\nOptional:\n" + \
                    f"Keep LHE: {self.args.keep_lhe}\n" + \
                    f"{'='*50}\n"
        self.logger.info(out_str)


    def initialize_logging(self):
        self.logger = logging.getLogger('.')
        self.logger.setLevel(logging.INFO)
    
    
    def split_lhe(self, input_file: str, num_parts: int) -> int:
        MG_CLASS_DIR  = "/work/home/psriling/MsProject/MG5_aMC_v3_5_3/madgraph/various"  # Directory containing lhe_parser.py
        MADGRAPH_DIR  = "/work/home/psriling/MsProject/MG5_aMC_v3_5_3/"
        
        sys.path.append(MADGRAPH_DIR)
        sys.path.append(MG_CLASS_DIR)
        
        from lhe_parser import EventFile
        
        old_name = input_file
        new_name = input_file.replace(".lhe", "")
        
        os.rename(old_name, new_name)
        lhe                 = EventFile(new_name)
        nevents             = len(lhe)
        EventFile.parsing   = False
        lhe.split(
            partition       = [nevents // num_parts for i in range(num_parts)], cwd=".", zip=True
        )
        os.rename(new_name, old_name)

        return nevents
        
    
    def init_lhe_chain(self):
        self.logger.info("Splitting LHE file")
        if not os.path.exists(self.args.lhe_path): 
            raise FileNotFoundError(f"LHE file {self.args.lhe_path} does not exist")

        if self.args.lhe_path.endswith(".gz"):
            os.system(f"gunzip {self.args.lhe_path}")
            self.lhe_path = self.args.lhe_path.replace(".gz", "")
        else:
            self.lhe_path = self.args.lhe_path
            
        self.nevents = self.split_lhe(self.lhe_path, self.args.n_threads) // self.args.n_threads
        prefix_name  = self.lhe_path.replace(".lhe", "")
        
        for i in range(self.args.n_threads):
            old_name = f"{prefix_name}_{i}.lhe.gz"
            new_name = f"events_{self.seeds[i]}.lhe.gz"
            os.rename(old_name, new_name)
            
        self.lhe_paths = [f"events_{seed}.lhe.gz" for seed in self.seeds]
        
    
    def init_gridpack_chain(self):
        self.logger.info("Extracting gridpack")
        if not os.path.exists(self.args.gridpack_path): 
            raise FileNotFoundError(f"Gridpack path {self.args.gridpack_path} does not exist")
        
        if self.args.gridpack_path.endswith(".tar.gz"):
            self.gridpack_path   = self.args.gridpack_path
            self.gridpack_prefix = self.args.gridpack_path.split(".tar.gz")[0]
        else:
            os.system(f"tar -czf {self.args.gridpack_path}.tar.gz {self.args.gridpack_path}")
            self.gridpack_path   = f"{self.args.gridpack_path}.tar.gz"
            self.gridpack_prefix = self.args.gridpack_path
        
        for i, seed in enumerate(self.seeds):
            self.logger.info(f"Extracting gridpack for seed {seed} [{i+1}/{len(self.seeds)}]")
            os.system(f"mkdir {self.gridpack_prefix}_{seed}")
            os.system(f"tar -xzf {self.gridpack_path} -C {self.gridpack_prefix}_{seed}")
        
        self.gridpack_paths = [f"./{self.gridpack_prefix}_{seed}/{self.gridpack_prefix}/" for seed in self.seeds]
        self.lhe_paths      = [f"{self.gridpack_prefix}_{seed}/events.lhe.gz" for seed in self.seeds]
        self.nevents        = self.args.n_events // self.args.n_threads # Number of events per thread
        
        if self.nevents * self.args.n_threads != self.args.n_events:
            self.logger.warning(f"Number of events {self.args.n_events} is not divisible by number of threads {self.args.n_threads}")
            self.logger.warning(f"Each thread will generate {self.nevents} events")
        if self.nevents < 1 : raise ValueError("Number of events per thread is less than 1")


    def init_info_csv(self):
        '''Initialize info.csv file'''
        if not os.path.exists("info.csv"):
            with open("info.csv", "w") as f:
                f.write("seed,expected_nevents,lhe_nevents,gained_nevents,runtime,madgraph,pythia,delphes\n")


    def write_info_csv(self, log_handler : LogHandler):
        '''Write to info.csv file'''
        # Expected columns: seed, expected_nevents, lhe_nevents, gained_nevents, runtime, madgraph, pythia, delphes
        with open("info.csv", "a") as f:
            for i, seed in enumerate(self.seeds):
                log_dict = log_handler.get_log_content(seed)
                runtime  = self.runners[i].runtime # Dictionary
                
                f.write(f"{seed},{self.nevents},"+
                        f"{log_dict['requested']},"+
                        f"{log_dict['accepted']},"+
                        f"{runtime['overall']},"+
                        f"{runtime['gridpack']},"+
                        f"{runtime['pythia']},"+
                        f"{runtime['delphes']}\n")


    def get_fullchain_args(self, seed):
        args = {
            "seed"              : seed,
            "run_gridpack"      : self.args.run_gridpack,
            "run_pythia"        : True,
            "run_delphes"       : True,
            "mg5_pythia_path"   : self.args.mg5_pythia_interface,
            "pythia_card"       : self.args.pythia_card,
            "delphes_card"      : self.args.delphes_card,
            "delhpes_path"      : self.args.delphes_path,
            "log_name"          : f"log_{seed}.log",
            "log_level"         : "INFO",
        }
        
        if self.args.run_gridpack:
            args["gridpack_path"]    = self.gridpack_paths[self.seeds.index(seed)]
            args["gridpack_nevents"] = self.nevents
            args["gridpack_verbose"] = False
            args["remove_gridpack"]  = True
            args["keep_lhe"]         = self.args.keep_lhe
            args["lhe_outfile"]      = f"events_{seed}.lhe.gz"
            
        if self.args.use_lhe:
            args["lhe_infile"]       = self.lhe_paths[self.seeds.index(seed)]
            
        return args
    
    
    def run_fullchain_parallel(self):
        self.logger.info("Running full chain in parallel")
        start_time = time.time()
        '''Run full chain in parallel'''
        with Pool(self.args.n_threads) as p:
            self.runners = p.map(run_job, self.fullChainArgs)
        self.runtime = time.time() - start_time
        self.logger.info("All jobs finished")
        self.cleanup()
        self.summarize_results()
    
    
    def cleanup(self):
        '''Cleanup after running full chain'''
        if self.args.run_gridpack:
            for seed in self.seeds:
                os.system(f"rm -r {self.gridpack_prefix}_{seed}")
        
        if not os.path.exists("logs"): os.system("mkdir logs")
        for seed in self.seeds:
            if os.path.exists(f"log_{seed}.log"):
                os.system(f"mv log_{seed}.log logs/")
            else:
                self.logger.warning(f"Log file log_{seed}.log does not exist")


 
    def summarize_results(self):
        log_handler = LogHandler(
            logs_dir = "logs",
            info_path= "info.csv",
        )

        self.write_info_csv(log_handler)
        ljust_len       = 20
        nevents         = self.nevents * self.args.n_threads
        total_accepted  = sum([log_handler.get_log_content(seed)["accepted"] for seed in self.seeds])
        efficiency      = total_accepted / self.args.n_events
        
        
        # Total runtime
        runner_list     = [runner.runtime for runner in self.runners]
        madgraph_total  = sum([runtime["gridpack"] for runtime in runner_list])
        pythia_total    = sum([runtime["pythia"] for runtime in runner_list])
        delphes_total   = sum([runtime["delphes"] for runtime in runner_list])
        
        # Overall runtime
        multi_runtime   = self.runtime
        single_runtime  = madgraph_total + pythia_total + delphes_total
        
        # Percentage of total runtime
        madgraph_percent = 100 * madgraph_total / single_runtime
        pythia_percent   = 100 * pythia_total / single_runtime
        delphes_percent  = 100 * delphes_total / single_runtime
        
        # Average runtime per thread
        madgraph_per_thread = madgraph_total / self.args.n_threads
        pythia_per_thread   = pythia_total / self.args.n_threads
        delphes_per_thread  = delphes_total / self.args.n_threads
        
        # Event rate
        multi_rate      = nevents / multi_runtime
        single_rate     = nevents / single_runtime
        
        # CPU efficiency
        cpu_efficiency  = single_runtime / multi_runtime
    
        out_str   = "\n"
        out_str  += log_handler.get_summarize(format="str") + \
                    f"\n{'='*20} Performance Report {'='*20}\n" + \
                    f"Expected events: {nevents}\n" + \
                    f"Gained events: {total_accepted}\n" + \
                    f"Efficiency: {efficiency:.4f}\n" + \
                    f"\n>Runtime (total, single-core)\n" + \
                    f"\t- Madgraph :\t" + f'{format_time(madgraph_total)}'.ljust(ljust_len) + f"({madgraph_percent:.2f}%)\n" + \
                    f"\t- Pythia   :\t" + f'{format_time(pythia_total)}'.ljust(ljust_len) + f"({pythia_percent:.2f}%)\n" + \
                    f"\t- Delphes  :\t" + f'{format_time(delphes_total)}'.ljust(ljust_len) + f"({delphes_percent:.2f}%)\n" + \
                    f"\n>Average runtime per thread:\n" + \
                    f"\t- Madgraph : \t{format_time(madgraph_per_thread)}\n" + \
                    f"\t- Pythia   : \t{format_time(pythia_per_thread)}\n" + \
                    f"\t- Delphes  : \t{format_time(delphes_per_thread)}\n" + \
                    f"\n>Overall runtime:\n" + \
                    f"\t- Single-core :\t" + f'{format_time(single_runtime)}'.ljust(ljust_len) + f"{single_rate:.2f} evts/s\n" + \
                    f"\t- Multi-cores :\t" + f'{format_time(multi_runtime)}'.ljust(ljust_len) + f"{multi_rate:.2f} evts/s\n" + \
                    f"\n>CPU efficiency: {cpu_efficiency:.2f}\n"
                    
        # Logging
        self.logger.info(out_str)
        
            
    
    
if __name__ == "__main__":
    fullChainInterface().run_fullchain_parallel()
    # Example usage:
    # python fullchain_interface.py --run_gridpack --gridpack_path <path_to_gridpack> --n_events 1000 --n_threads 5
    # Using existing LHE file:
    # python fullchain_interface.py --use_lhe --lhe_path <path_to_lhe> --n_threads 5


