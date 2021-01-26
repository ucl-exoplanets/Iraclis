
import os

__version__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__version__.txt')).read()

from .pipeline import *
