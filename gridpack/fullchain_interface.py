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

from classes import FullChainRunner, format_time

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
    # ********** CASE I: RUN GRIDPACK **********
    parser.add_argument("--run_gridpack",           action="store_true",                    help="Run gridpack")
    parser.add_argument("--gridpack_path",          type=str,   default="gridpacks",        help="Path to the gridpacks")
    parser.add_argument("--n_events",               type=int,   default=1000,               help="Number of events to generate")
    # ********** CASE II: GIVEN LHE **********
    parser.add_argument("--use_lhe",                action="store_true",                    help="Use LHE file")
    parser.add_argument("--lhe_path",               type=str,   default="lhe",              help="Path to the LHE file")

    return parser.parse_args()



def run_job(args): # FOR ASSIGN TO MULTIPROCESSING
    runner = FullChainRunner(**args)
    runner.run_full_chain()
    return runner



class fullChainInterface:
    def __init__(self):
        self.args = get_args()
        self.initialize_logging()
        
        if self.args.run_gridpack and self.args.use_lhe:
            raise ValueError("Cannot run gridpack and use LHE at the same time")
        if not self.args.run_gridpack and not self.args.use_lhe:
            raise ValueError("Must run gridpack or use LHE")
        
        # ********** SEEDS **********
        if self.args.seed is None: self.args.seed = random.randint(1,1000000)
        self.seeds = [self.args.seed + i for i in range(self.args.n_threads)]
        self.logger.info(f"Seeds: {self.seeds}")
        
        # ********** INITIALIZATION **********
        if self.args.use_lhe        : self.init_LHE_chain()
        if self.args.run_gridpack   : self.init_gridpack_chain()
        
        self.fullChainArgs = [self.get_fullchain_args(seed) for seed in self.seeds]
        
        if not os.path.exists("run_status") : os.system("mkdir run_status")
        else                                : os.system("rm -r run_status/*")


    def initialize_logging(self):
        self.logger = logging.getLogger('Interface')
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
        
    
    def init_LHE_chain(self):
        if not os.path.exists(self.args.lhe_path): 
            raise FileNotFoundError(f"LHE file {self.args.lhe_path} does not exist")

        if self.args.lhe_path.endswith(".gz"):
            os.system(f"gunzip {self.args.lhe_path}")
            self.lhe_path = self.args.lhe_path.replace(".gz", "")
        else:
            self.lhe_path = self.args.lhe_path
            
        self.nevents = self.split_lhe(self.lhe_path, self.args.n_threads)
        prefix_name  = self.lhe_path.replace(".lhe", "")
        
        for i in range(self.args.n_threads):
            old_name = f"{prefix_name}_{i}.lhe.gz"
            new_name = f"events_{self.seeds[i]}.lhe.gz"
            os.rename(old_name, new_name)
            
        self.lhe_paths = [f"events_{seed}.lhe.gz" for seed in self.seeds]
        
    
    def init_gridpack_chain(self):
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
            args["gridpack_nevents"] = self.args.n_events
            args["gridpack_verbose"] = False
            args["LHE_outfile"]      = f"events_{seed}.lhe.gz"
            
        if self.args.use_lhe:
            args["LHE_infile"]       = self.lhe_paths[self.seeds.index(seed)]
            
        return args
    
    
    def run_fullchain_parallel(self):
        '''Run full chain in parallel'''
        with Pool(self.args.n_threads) as p:
            self.runners = p.map(run_job, self.fullChainArgs)
        
        self.logger.info("All jobs finished")
    
    # TODO: Implement method to clean up (mostly already done in FullChainRunner)
    '''
        Clean up:
        - Remove gridpacks
        - Remove LHE files (if split)
    '''
    
    # TODO: Implement method to summarize results
    '''
        Physics quantities to summarize:
        - Number of requested events
        - Number of events generated (gained after running full chain)
        - Efficiency (events generated / requested events)
        - Cross section
        
        Performance metrics:
        - Time taken for each step (single threaded sum)
        - Percentage of time spent in each step
        
    '''
    # def summarize_results(self):
    
    
if __name__ == "__main__":
    fullChainInterface().run_fullchain_parallel()
    # Example usage:
    # python fullchain_interface.py --run_gridpack --gridpack_path <path_to_gridpack> --n_events 1000 --n_threads 5
    # Using existing LHE file:
    # python fullchain_interface.py --use_lhe --lhe_path <path_to_lhe> --n_threads 5


