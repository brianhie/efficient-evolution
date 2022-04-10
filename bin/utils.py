from collections import Counter
import datetime
from dateutil.parser import parse as dparse
import errno
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import random
import scipy
from scipy.sparse import csr_matrix, dok_matrix
import scipy.stats as ss
import seaborn as sns
import sys
import time
import torch
import warnings

from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from Bio import Seq, SeqIO

np.random.seed(1)
random.seed(1)

AAs = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
]

def tprint(string):
    string = str(string)
    sys.stdout.write(str(datetime.datetime.now()) + ' | ')
    sys.stdout.write(string + '\n')
    sys.stdout.flush()

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
        
def deep_mutational_scan(sequence, exclude_noop=True):
    for pos, wt in enumerate(sequence):
        for mt in AAs:
            if exclude_noop and wt == mt:
                continue
            yield (pos, wt, mt)
