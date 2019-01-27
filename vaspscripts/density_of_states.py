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

def get_orbital_density_data(spd_dos, orbital_type, spin_type):
    if orbital_type == "s":
        orbital_projected_dos = spd_dos[OrbitalType.s]
    elif orbital_type = "p":
        orbital_projected_dos = spd_dos[OrbitalType.p]
    elif orbital_type = "d":
        orbital_projected_dos = spd_dos[OrbitalType.d]

    if spin_type == "up":
        orbial_density = orbital_projected_dos.densities[Spin.up]
    else
        orbial_density = -orbital_projected_dos.densities[Spin.down]

    return orbial_density

def plot_dos(dos_data, desired_orbitals):
    if desired_orbitals != ["total"]:
        element_spd_dos_data = get_element_spd_dos(dos_data, desired_orbitals)
    return None

