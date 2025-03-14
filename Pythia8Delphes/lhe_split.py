import sys

class_path = "/work/home/psriling/MsProject/MG5_aMC_v3_5_3/madgraph/various"  # Directory containing lhe_parser.py
MADGRAPH_DIR = "/work/home/psriling/MsProject/MG5_aMC_v3_5_3/"
sys.path.append(MADGRAPH_DIR)
sys.path.append(class_path)
from lhe_parser import EventFile
import os


def split_lhe(input_file, num_parts):
    old_name = input_file
    new_name = input_file.replace(".lhe", "")

    os.rename(old_name, new_name)

    lhe = EventFile(new_name)
    EventFile.parsing = False
    nevents = len(lhe)

    lhe.split(
        partition=[nevents // num_parts for i in range(num_parts)], cwd=".", zip=True
    )

    os.rename(new_name, old_name)

    return nevents
