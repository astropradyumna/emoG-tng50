'''
This is a second attempt to evolve subhaolos from TNG50hydro 
'''
import numpy as np
import matplotlib.pyplot as plt 
import os 
import sys
os.environ["USE_LZMA"] = "0"
import pandas as pd
sys.path.append(os.path.abspath('/bigdata/saleslab/psadh003/tng50/dwarf_formation'))
from errani_plus_tng_subhalo import Subhalo
from tqdm import tqdm
import galpy
import IPython
import illustris_python as il
from matplotlib.backends.backend_pdf import PdfPages
from subhalo_profiles import ExponentialProfile, NFWProfile
import warnings
from populating_stars import *
from joblib import Parallel, delayed #This is to parallelize the code
from hydro_path_for_files import * #This is to get the paths for the files
from constants import * #This is to get the constants


# Suppress the lzma module warning
# warnings.filterwarnings("ignore", category=UserWarning, module="pandas.compat")
warnings.simplefilter(action='ignore', category=FutureWarning)

