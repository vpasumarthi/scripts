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

def plot_dos(dos_data, desired_orbitals, dst_path, plot_properties):
    if desired_orbitals != ["total"]:
        element_spd_dos_data = get_element_spd_dos(dos_data, desired_orbitals)

    # setup matplotlib plot
    plt.switch_backend('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    font = {'family': plot_properties["font_family",
            'size': plot_properties["font_size"]}
    plt.rc('font', **font)

    if plot_properties["zero_at_fermi"] == "yes":
        energy_data = dos_data["tdos"].energies - dos_data["efermi"]
    else:
        energy_data = dos_data["tdos"].energies

    for index, orbital in enumerate(desired_orbitals):
        element = orbital.split("_")[0]
        orbital_type = orbital.split("_")[1][-1]
        
        if plot_properties["spin_type"] == "up":
            spin_type = "up"
            ax.plot(energy_data,
                    get_orbital_density_data(element_spd_dos_data[element], orbital_type, spin_type),
                    color=plot_properties["color_list"][index],
                    label=orbital,
                    lw=plot_properties["line_width"])
        elif plot_properties["spin_type"] == "down":
            spin_type = "down"
            ax.plot(energy_data,
                    get_orbital_density_data(element_spd_dos_data[element], orbital_type, spin_type),
                    color=plot_properties["color_list"][index],
                    label=orbital,
                    lw=plot_properties["line_width"])
        elif plot_properties["spin_type"] == "both":
            spin_type = "up"
            ax.plot(energy_data,
                    get_orbital_density_data(element_spd_dos_data[element], orbital_type, spin_type),
                    color=plot_properties["color_list"][index],
                    label=orbital,
                    lw=plot_properties["line_width"])
            spin_type = "down"
            ax.plot(energy_data,
                    get_orbital_density_data(element_spd_dos_data[element], orbital_type, spin_type),
                    color=plot_properties["color_list"][index],
                    lw=plot_properties["line_width"])
    
    if plot_properties["zero_at_fermi"] == "yes":
        ax.xlabel(plot_properties["x_label_fermi0"]
    else:
        ax.xlabel(plot_properties["x_label"]
    ax.ylabel(plot_properties["y_label"])
    ax.legend(prop={'size': plot_properties["legend_font_size"]})
    plt.savefig(dst_path / plot_properties["output_file_name"] / "." / plot_properties["output_file_type"], format=plot_properties["output_file_type"])
    return None

