import math
from functools import cached_property
from pydantic import BaseModel, ConfigDict, Field
from functools import cached_property
from typing import Callable, Literal
import ctypes
from pathlib import Path
import numpy as np
from numpy.ctypeslib import ndpointer

L = {
    "S": [
        (0, 0, 0),
    ],
    "P": [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
    ],

    "D": [
        (2, 0, 0),
        (1, 0, 1),
        (0, 2, 0),
        (0, 1, 1),
    ],
    "F": [
        (3, 0, 0),
        (2, 1, 0),
        (2, 0, 1),
        (1, 2, 0),
        (1, 1, 1),
        (1, 0, 2),
        (0, 3, 0),
        (0, 2, 1),
        (0, 1, 2),
        (0, 0, 3),
    ],
}

CWD = Path(__file__).resolve().parent.parent
GTO_LIB = ctypes.CDLL(str(CWD/"c/gto.dylib"))

## Define the argument and return types for the gto function in the C library
GTO_LIB.gto.argtypes = [
    ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # X, Y, Z
    ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # x, y, z
    ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),  # nx, ny, nz
    ctypes.c_double                                     # alpha
]                                                     

GTO_LIB.gto.restype = ctypes.c_double


GTO_LIB.sto_ng.argtypes = [
    ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # X, Y, Z
    ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # x, y, z
    ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # cc
    ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # alpha
    ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),    # nx
    ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),    # ny
    ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),    # nz
    ctypes.c_int]                                       # n
GTO_LIB.sto_ng.restype = ctypes.c_double