from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import warnings
warnings.filterwarnings("ignore",
                        message='Matplotlib is building the font cache using fc-list. This may take a moment.')
warnings.filterwarnings("ignore",
                        message='The installed version of numexpr 2.4.4 is not supported in pandas and will be not be used')

import matplotlib
matplotlib.use('TkAgg')

import glob
import os
import shutil
import time
import sys
import datetime
import pickle
import docopt
import gzip
import socket

import numpy as np
from sklearn.decomposition import FastICA, PCA

import scipy
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy import interpolate


from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from astropy.io import fits as pf
import ephem
import pylightcurve as plc

import os
import sys
import glob
import time
import pickle
import shutil
import socket

if sys.version_info[0] > 2:
    from urllib.request import urlretrieve
else:
    from urllib import urlretrieve
    input = raw_input