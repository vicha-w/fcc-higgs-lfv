import os
import sys
import logging
import subprocess
import random
import time

import pandas as pd


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



def format_time(s : float) -> str:
    '''Convert seconds to human-readable format'''
    days, remainder     = divmod(s,         86400)
    hours, remainder    = divmod(remainder, 3600)
    minutes, seconds    = divmod(remainder, 60)
    return_str          = f"{seconds:.2f}s"
    
    if minutes > 0  : return_str = f"{minutes:.0f}m " + return_str
    if hours > 0    : return_str = f"{hours:.0f}h "   + return_str
    if days > 0     : return_str = f"{days:.0f}d "    + return_str
    
    return return_str



class GridpackRunner:
    def __init__(
        self,
        nevents         : int  = None,
        seed            : int  = None,
        outfile         : str  = None,
        log_file        : str  = None,
        granularity     : int  = 1,
        gridpack_path   : str  = "./madevent",
        verbose         : bool = True,
    ):
        
        self.verbose        = verbose
        self.nevents        = nevents
        self.seed           = seed
        self.granularity    = granularity
        self.gridpack_path  = gridpack_path
        self.outfile        = outfile
        self.log_file       = log_file
        
        self.initialize_logging()

        if os.path.exists(gridpack_path):
            if not os.path.isdir(gridpack_path) and gridpack_path.endswith(".tar.gz"):
                self.tarball = True
            else:
                self.tarball = False
        else:
            raise FileNotFoundError(f"Gridpack path {gridpack_path} does not exist.")

        if outfile is None  : self.outfile = f"events_{seed}.lhe.gz"

    def dump(self):
        """Dump all attributes (for debugging)"""
        for attr in vars(self):
            print(f"{attr}: {getattr(self, attr)}")

    def initialize_logging(self):
        if self.verbose : level = logging.INFO
        else            : level = logging.CRITICAL
        self.logger             = logging.getLogger()
        self.logger.setLevel(level)

    def extract_tarball(self):
        if self.tarball:
            self.logger.debug(f"Extracting tarball {self.gridpack_path}")
            os.system(f"tar -xvf {self.gridpack_path}")
            self.gridpack_path = self.gridpack_path.replace(".tar.gz", "")
            self.logger.debug(f"Extracted to {self.gridpack_path}")

    def postprocess(self):
        event_file_path = f"{self.gridpack_path}/Events/GridRun_{self.seed}/unweighted_events.lhe.gz"
        
        if os.path.exists(event_file_path)  : os.rename(event_file_path, self.outfile)
        else:
            event_file_path = f"./Events/GridRun_{self.seed}/unweighted_events.lhe.gz"
            if os.path.exists(event_file_path):
                os.rename(event_file_path, self.outfile)
                os.system("rm -rf Events Cards P* *.dat randinit &> /dev/null")
                
        self.logger.debug(f"write ./{self.outfile}")

    def run(self):
        self.extract_tarball()
        
        self.logger.debug(f"Now generating {self.nevents} events with random seed {self.seed} and granularity {self.granularity}")
        base_cmd = f"{self.gridpack_path}/bin/gridrun {self.nevents} {self.seed} {self.granularity}"
        
        if self.verbose:
            if self.log_file is None    : gridrun_cmd = base_cmd
            else                        : gridrun_cmd = f"{base_cmd} | tee -a {self.log_file}"
        else:
            if self.log_file is None    : gridrun_cmd = f"{base_cmd} &> /dev/null"
            else                        : gridrun_cmd = f"{base_cmd} >> {self.log_file} 2>&1"
        
        self.logger.debug(f"Running command: {gridrun_cmd}")
        subprocess.run(gridrun_cmd, shell=True, check=True)

        self.postprocess()
        
        

class FullChainRunner:
    def __init__(
        self,
        seed: int = None,
        # === Gridpack ===
        run_gridpack      : bool = True,
        gridpack_verbose  : bool = True,
        keep_lhe          : bool = False,
        remove_gridpack   : bool = False,
        gridpack_nevents  : int  = None,
        lhe_outfile       : str  = None,
        gridpack_path     : str  = "./madevent",
        # === Pythia ===
        run_pythia        : bool = True,
        keep_hepmc        : bool = False,
        zip_hepmc         : bool = True,
        lhe_infile        : str  = None,
        hepmc_outfile     : str  = None,
        pythia_card       : str  = None,
        mg5_pythia_path   : str  = "/work/app/pythia8/MGInterface/1.3/MG5aMC_PY8_interface",
        # === Delphes ===
        run_delphes       : bool = True,
        hepmc_infile      : str  = None,
        root_outfile      : str  = None,
        delphes_card      : str  = "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/FCChh_I.tcl",
        delhpes_path      : str  = "/work/home/psriling/MsProject/Delphes-3.5.0/DelphesHepMC2",
        # === Others ===
        cleanup           : bool = True,
        save_log          : bool = True,
        log_name          : str  = None,
        log_level         : str  = "INFO", # DEBUG, INFO, WARNING, ERROR, CRITICAL
    ):
        # ****** Set all variables ******
        # === Gridpack ===
        self.run_gridpack       = run_gridpack
        self.lhe_outfile        = lhe_outfile
        self.gridpack_nevents   = gridpack_nevents
        self.gridpack_verbose   = gridpack_verbose
        self.gridpack_path      = gridpack_path
        self.keep_lhe           = keep_lhe
        self.remove_gridpack    = remove_gridpack
        # === Pythia ===
        self.run_pythia         = run_pythia
        self.lhe_infile         = lhe_infile
        self.hepmc_outfile      = hepmc_outfile
        self.mg5_pythia_path    = mg5_pythia_path
        self.pythia_card        = pythia_card
        self.zip_hepmc          = zip_hepmc
        self.keep_hepmc         = keep_hepmc
        # === Delphes ===
        self.hepmc_infile       = hepmc_infile
        self.root_outfile       = root_outfile
        self.delphes_card       = delphes_card
        self.delhpes_path       = delhpes_path
        self.run_delphes        = run_delphes
        # === Others ===
        self.cleanup            = cleanup
        self.save_log           = save_log
        self.log_name           = log_name
        self.log_level          = log_level
        self.runtime            = {}
        
        # ****** Initialize variables ******
        if seed is None         : self.seed = random.randint(0, 100000)
        else                    : self.seed = seed
        
        self.var_init_others()
        self.initialize_logging()
        if self.run_gridpack    : self.var_init_gridpack()
        if self.run_pythia      : self.var_init_pythia()
        if self.run_delphes     : self.var_init_delphes()
        
        # ****** Check for mismatches ******
        if self.run_gridpack    : self.lhe_mismatch = self.lhe_outfile != self.lhe_infile
        else                    : self.lhe_mismatch = False
        if self.run_pythia      : self.hepmc_mismatch = self.hepmc_outfile != self.hepmc_infile
        else                    : self.hepmc_mismatch = False


    def initialize_logging(self):
        self.logger         = logging.getLogger("FullChainRunner")
        log_level_mapping   = {
            "DEBUG"    : logging.DEBUG,
            "INFO"     : logging.INFO,
            "WARNING"  : logging.WARNING,
            "ERROR"    : logging.ERROR,
            "CRITICAL" : logging.CRITICAL,
        }

        if self.log_level not in log_level_mapping:
            raise ValueError(f"Invalid log level {self.log_level}")

        self.logger.setLevel(log_level_mapping[self.log_level])
            
            
    def var_init_gridpack(self):
        if self.gridpack_nevents is None : raise ValueError("gridpack_nevents must be specified")
        if self.lhe_outfile is None      : self.lhe_outfile = f"events_{self.seed}.lhe.gz"
        
        self.gridpack_runner = GridpackRunner(
            nevents         = self.gridpack_nevents,
            seed            = self.seed,
            outfile         = self.lhe_outfile,
            gridpack_path   = self.gridpack_path,
            verbose         = self.gridpack_verbose,
        )
        
        if self.save_log    : self.gridpack_runner.log_file = self.log_name
    
    
    def var_init_pythia(self):
        if self.lhe_infile is None:
            if self.run_gridpack        : self.lhe_infile = self.lhe_outfile
            else                        : raise ValueError("lhe_infile must be specified if not running gridpack")
            
        if self.hepmc_outfile is None   : self.hepmc_outfile = f"events_{self.seed}.hepmc.gz"
        
        for path in [self.mg5_pythia_path, self.pythia_card]:
            if not os.path.exists(path) : raise FileNotFoundError(f"Path {path} does not exist.")

        new_pythia_card                 = f"pythia_card_{self.seed}.dat"
        os.system(f"cp {self.pythia_card} {new_pythia_card}")
        self.pythia_card                = new_pythia_card
        
        with open(self.pythia_card, "r+") as f:
            contents = f.read()
            contents = contents.replace("$SEED",       str(self.seed))  \
                               .replace("$INFILE",     self.lhe_infile) \
                               .replace("$OUTFILE",    self.hepmc_outfile)
            f.seek(0)
            f.write(contents)
            f.truncate()

                
    def var_init_delphes(self):
        if self.hepmc_infile is None:
            if self.run_pythia          : self.hepmc_infile = self.hepmc_outfile
            else                        : raise ValueError("hepmc_infile must be specified if not running pythia")
        if self.root_outfile is None    : self.root_outfile = f"events_{self.seed}.root"
        
        for path in [self.delphes_card, self.delhpes_path]:
            if not os.path.exists(path) : raise FileNotFoundError(f"Path {path} does not exist.")
    
    
    def var_init_others(self):
        if self.save_log:
            if self.log_name is None            : self.log_name = f"log_{self.seed}.log"
            if os.path.exists(self.log_name)    : os.system(f"rm -f {self.log_name}")
        if os.path.exists("run_status")         : self.multicore = True
        else                                    : self.multicore = False
        

    def update_status(self, status : str) -> None:
        if self.multicore:
            os.system(f"rm -f run_status/{self.seed}.*")
            os.system(f"touch run_status/{self.seed}.{status}")
    
    
    def check_status(self):
        file_list       = os.listdir("run_status")
        groups          = set(file.split(".")[1] for file in file_list)
        group_counts    = {group: len([file for file in file_list if file.split(".")[1] == group]) for group in groups}
        sort_guide      = ["Gridpack", "Pythia", "Delphes"]
        out_str         = "Running: "

        for group in sort_guide:
            if group in group_counts:
                out_str += f"{group_counts[group]} {group}, "
        
        for group in group_counts:
            if group not in sort_guide and group != "Done":
                out_str += f"{group_counts[group]} {group}, "

        out_str     = out_str[:-2]
        done_count  = group_counts.get("Done", 0)  # Use .get() to avoid KeyError
        out_str     += f" Done: {done_count}/{len(file_list)}"

        self.logger.info(out_str)


    def run_gridpack_chain(self):
        self.update_status("Gridpack")
        if self.multicore           : self.check_status()
        start_time                  = time.time()
        self.gridpack_runner.run()
        self.runtime["gridpack"]    = time.time() - start_time
    
    
    def run_pythia_chain(self):
        self.update_status("Pythia")
        self.logger.debug(f"Running Pythia with seed {self.seed}")
        base_cmd                    = f"{self.mg5_pythia_path} {self.pythia_card}"
        start_time                  = time.time()
        if self.multicore           : self.check_status()
        if self.lhe_mismatch        : os.system(f"ln -s {self.lhe_outfile} {self.lhe_infile}")
        if self.save_log            : os.system(f"{base_cmd} >> {self.log_name} 2>&1")
        else                        : os.system(base_cmd)
        self.runtime["pythia"]      = time.time() - start_time
        
        
    def run_delphes_chain(self):
        self.update_status("Delphes")
        self.logger.debug(f"Running Delphes with seed {self.seed}")
        base_cmd                    = f"{self.delhpes_path} {self.delphes_card} {self.root_outfile} {self.hepmc_infile}"
        start_time                  = time.time()
        if self.multicore           : self.check_status()
        if self.hepmc_mismatch      : os.system(f"ln -s {self.hepmc_outfile} {self.hepmc_infile}")
        if self.save_log            : os.system(f"{base_cmd} >> {self.log_name} 2>&1")
        else                        : os.system(base_cmd)
        self.runtime["delphes"]     = time.time() - start_time
        self.update_status("Done")
    
    
    def cleanfiles(self):
        if self.keep_lhe: 
            if not os.path.exists("./LHE") : os.system("mkdir ./LHE")
        # Clean gridpack
        if self.run_gridpack:
            if not self.keep_lhe    : os.system(f"rm -rf {self.lhe_outfile}")
            else                    : os.system(f"mv {self.lhe_outfile} ./LHE/")
            if self.remove_gridpack : os.system(f"rm -rf {self.gridpack_path}")
        # Clean pythia
        if self.run_pythia:
            os.system(f"rm -rf {self.pythia_card}")
            if self.keep_hepmc: 
                if self.zip_hepmc   : os.system(f"gzip {self.hepmc_outfile} && mv {self.hepmc_outfile}.gz ./HEPMC/")
                else                : os.system(f"mv {self.hepmc_outfile} ./HEPMC/")
            else                    : os.system(f"rm -rf {self.hepmc_outfile}")
        # Clean delphes
        if self.run_delphes:
            if os.path.exists("./ROOT") : os.system(f"mv {self.root_outfile} ./ROOT/")
            else                        : os.system(f"mkdir ./ROOT && mv {self.root_outfile} ./ROOT/")
        
        if os.path.exists("djrs.dat") : os.system("rm djrs.dat")
        if os.path.exists("pts.dat")  : os.system("rm pts.dat")
        
        
    def run_full_chain(self):
        if self.run_gridpack    : self.run_gridpack_chain()
        if self.run_pythia      : self.run_pythia_chain()
        if self.run_delphes     : self.run_delphes_chain()
        if self.multicore       : self.check_status()
        if self.cleanup         : self.cleanfiles()
        
        self.runtime["overall"] = sum(self.runtime.values())
        self.logger.debug(f"Overall runtime: {format_time(self.runtime['overall'])}")
        for item in self.runtime.items():
            self.logger.debug(f"{item[0]}: {format_time(item[1])}")
        


'''
    Class signature:
    - Read the log file in the current directory (or given path if specified)
    - Read the info.csv file in the current directory
    - Summarize outputs of the jobs: gained event, efficiencies, cross-sections, etc.
    
    This supposed to be a sup-class for the FullChainRunner class, when handling multiple jobs.
'''
class LogHandler:
    def __init__(
        self,
        logs_paths  : str = None,
        info_path   : str = None,
    ):
        self.logs_paths     = logs_paths    # List of log file paths
        self.info_path      = info_path     # Path to info.csv file
        self.multiple_logs  = False         # Flag for multiple log files
        self.log_contents   = []            # Contain dictionaries of log contents: list[dict]
        
        if self.log_paths is None : self.log_paths = self.get_log_paths()
        if self.info_path is None : self.info_path = self.get_info_path()
        
    
    def get_log_paths(self):
        log_files = [file for file in os.listdir() if file.endswith(".log")]
        if len(log_files) == 0  : raise FileNotFoundError("No log files found in the current directory.")
        elif len(log_files) > 1 : self.multiple_logs = True
        return log_files
    
    
    def get_info_path(self):
        info_files = [file for file in os.listdir() if file == "info.csv"]
        if len(info_files) == 0  : raise FileNotFoundError("No info.csv file found in the current directory.")
        elif len(info_files) > 1 : raise FileNotFoundError("Multiple info.csv files found in the current directory.")
        return info_files[0]
    
    
    def extract_seed(self, log_path : str) -> int:
        '''Extract seed from the log file name'''
        seed = log_path.split("_")[1].split(".")[0]
        return int(seed)
    
    
    def extract_cross_section(self, lines : list) -> tuple:
        '''Extract cross section from the log file: "Inclusive cross section: 2.17868639e-05  +-  3.08637916e-07 mb"'''
        for line in lines:
            if "Inclusive cross section" in line:
                num = line.split("cross section: ")[1]
                cross_section, unc = num.split("  +-  ")
                unc = unc.split()[0]
                return float(cross_section) * 1e9, float(unc) * 1e9
        return None, None
    
    
    def extract_gained_events(self, lines : list) -> tuple:
        for line in lines:
            if "| sum" in line:
                num = line.split("|")[-3]
                spt_num = num.split()
                selected, accepted = spt_num[0], spt_num[-1]
                return int(selected), int(accepted)
        return None, None
    
    
    def extract_info_row(self, seed : int) -> dict:
        '''Read info.csv file and extract the row with the given seed'''
        '''Header: seed,nevents,out_file,runtime,madgraph,pythia,delphes'''
        df  = pd.read_csv(self.info_path)
        row = df[df["seed"] == seed]
        return row.to_dict()
            
    
    def read_log(self) -> list:
        self.log_contents = [] 
        for log_path in self.logs_paths:
            with open(log_path, "r") as f: lines = f.read().split("\n")

            cross_section, unc  = self.extract_cross_section(lines)
            selected, accepted  = self.extract_gained_events(lines)
            seed                = self.extract_seed(log_path)
            info_row            = self.extract_info_row(seed)
            
            log_dict = {
                "seed"          : seed,                 # Seed
                "cross_section" : cross_section,        # Cross section in pb
                "uncertainty"   : unc,                  # Uncertainty in pb
                "selected"      : selected,
                "accepted"      : accepted,
                "requested"     : info_row.get("nevents", None),
                "info_row"      : info_row,             # Info row from info.csv
            }
            self.log_contents.append(log_dict)
        return self.log_contents
    
    
    def get_log_contents(self) -> list:
        if len(self.log_contents) == 0 : self.read_log()
        return self.log_contents
    

    def get_summarize(self, format : str = "str") -> str: 
        '''Summarize the log contents into a multiple line of string
            if format == "str": return a string,
            if format == "list"" return a string with split lines (\n)
        '''
        log_contents    = self.get_log_contents()
        ljust_width     = 20
        out_str         = f"{'='*19} Physics Summary {'='*19}\n"
        header          = "Seed".ljust(ljust_width) + \
                          "Requested events".ljust(ljust_width) + \
                          "Gained events".ljust(ljust_width) + \
                          "Efficiency".ljust(ljust_width) + \
                          "Cross-section (pb)".ljust(ljust_width) + \
                          "Uncertainty (pb)".ljust(ljust_width) + "\n"
        out_str += header
        
        for log in log_contents:
            seed        = str(log.get("seed", "None")).ljust(ljust_width)
            requested   = str(log.get("requested", "None")).ljust(ljust_width)
            gained      = str(log.get("gained", "None")).ljust(ljust_width)
            efficiency  = str(log.get("efficiency", "None")).ljust(ljust_width)
            cross_sec   = str(log.get("cross_section", "None")).ljust(ljust_width)
            unc         = str(log.get("uncertainty", "None")).ljust(ljust_width)
            
            out_str += seed + requested + gained + efficiency + cross_sec + unc + "\n"
        
        if format == "str"    : return out_str
        elif format == "list" : return out_str.split("\n")
        else                  : raise ValueError("Invalid format")
        
'''
    Test
    1. Test running standalone gridpack (GridpackRunner)
    1.x. Test 4 cases of verbose and log file
    2. Test running full chain (FullChainRunner)
    2.1. Run full chain with verbose output
    2.2. Run full chain without verbose output (save to log file)
'''
class Tester:
    def __init__(
        self,
        seed: int = 1,
        gridpack_path = "./madevent",
        gridpack_nevents = 100,
        lhe_outfile = "test.lhe.gz",
    ):
        self.seed = seed
        self.gridpack_path = gridpack_path
        self.gridpack_nevents = gridpack_nevents
        self.lhe_outfile = lhe_outfile
        self.gridpack_log = "gridpack.log"
        
        self.initialize_logging()
    
    def initialize_logging(self):
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
    
    def test_gridpack(self):
        
        self.logger.info("Testing GridpackRunner")
        
        # case 1: verbose is True, log file is None
        self.logger.info("Case 1: verbose is True, log file is None")
        gridpack_runner = GridpackRunner(
            nevents         = self.gridpack_nevents,
            seed            = self.seed,
            outfile         = self.lhe_outfile,
            gridpack_path   = self.gridpack_path,
            verbose         = True,
        )
        gridpack_runner.run()
        
        # case 2: verbose is False, log file is None
        self.logger.info("Case 2: verbose is False, log file is None")
        gridpack_runner = GridpackRunner(
            nevents         = self.gridpack_nevents,
            seed            = self.seed,
            outfile         = self.lhe_outfile,
            gridpack_path   = self.gridpack_path,
            verbose         = False,
        )
        gridpack_runner.run()
        
        # case 3: verbose is True, log file is specified
        self.logger.info("Case 3: verbose is True, log file is specified")
        gridpack_runner = GridpackRunner(
            nevents         = self.gridpack_nevents,
            seed            = self.seed,
            outfile         = self.lhe_outfile,
            gridpack_path   = self.gridpack_path,
            verbose         = True,
            log_file        = 'case3.log',
        )
        gridpack_runner.run()
        
        # case 4: verbose is False, log file is specified
        self.logger.info("Case 4: verbose is False, log file is specified")
        gridpack_runner = GridpackRunner(
            nevents         = self.gridpack_nevents,
            seed            = self.seed,
            outfile         = self.lhe_outfile,
            gridpack_path   = self.gridpack_path,
            verbose         = False,
            log_file        = 'case4.log',
        )
        gridpack_runner.run()
        

    def test_full_chain(self):
        full_chain_runner = FullChainRunner(
            seed            = self.seed,
            run_gridpack    = True,
            gridpack_nevents= self.gridpack_nevents,
            lhe_outfile     = self.lhe_outfile,
            gridpack_path   = self.gridpack_path,
            gridpack_verbose= False,
            pythia_card     = "pythia8_card.cmd",
            run_pythia      = True,
            run_delphes     = True,
            log_name        = "full_chain.log",
            log_level       = "DEBUG",
        )
        full_chain_runner.run_full_chain()
        
tester = Tester()

# tester.test_gridpack()    # Passed
# tester.test_full_chain()    # Passed