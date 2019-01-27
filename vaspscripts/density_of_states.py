#!/usr/bin/env python

# module file with functions relating to read and plot density of states data

import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType
from pymatgen.core.periodic_table import Element

def read_dos(src_path):
    # parsed object of vasprun.xml file
    dosrun = Vasprun(src_path / "vasprun.xml")
    # complete dos object which incorporates the total dos and all projected dos
    cdos = dosrun.complete_dos
    # total dos at the end of the run
    tdos = dosrun.tdos
    # fermi energy
    efermi = dosrun.efermi
    dos_data = {"dosrun" : dosrun,
                "cdos"   : cdos,
                "tdos"   : tdos,
                "efermi" : efermi}
    return dos_data

def get_element_spd_dos(dos_data, desired_orbitals):
    unique_elements = set([orbital.split("_")[0] for orbital in desired_orbitals])
    element_spd_dos_data = {}
    for element in unique_elements:
        atomic_number = Element(element).number
        element_spd_dos_data[element] = dos_data["cdos"].get_element_spd_dos(atomic_number)
    return element_spd_dos_data

