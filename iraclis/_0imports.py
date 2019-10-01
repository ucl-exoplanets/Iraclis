from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys

if sys.version_info[0] > 2:
    from urllib.request import urlretrieve
else:
    from urllib import urlretrieve
    input = raw_input

import glob
import gzip
import time
import numpy as np
import scipy

import pickle
import shutil
import socket
import datetime

import docopt
import pylightcurve as plc

matplotlib = plc.matplotlib
plt = plc.plt
pf = plc.pf
griddata = plc.griddata
curve_fit = plc.curve_fit
interp1d = plc.interp1d
FastICA = plc.FastICA
PCA = plc.PCA

import warnings
from scipy.optimize import OptimizeWarning