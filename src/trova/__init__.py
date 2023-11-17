import logging
import warnings

from .trova import *
from .functions import *
from ._version import get_versions 

__version__ = get_versions()['version']
__author__ = get_versions()['author']
__contact__ = get_versions()['contact']
__last_update__ = get_versions()['last_update']
del get_versions

try:
    	logging.lastResort
except AttributeError:
        logging.getLogger(__name__).addHandler(logging.StreamHandler())
