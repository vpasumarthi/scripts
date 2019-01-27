#!/usr/bin/env python

# module file with functions relating to read and plot density of states data

import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType


