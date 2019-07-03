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
        if element == "total":
            element_spd_dos_data[element] = ''
        else:
            atomic_number = Element(element).number
            element_spd_dos_data[element] = dos_data["cdos"].get_element_spd_dos(atomic_number)
    return element_spd_dos_data

def get_orbital_density_data(spd_dos, orbital_type, spin_type):
    if orbital_type == "s":
        orbital_projected_dos = spd_dos[OrbitalType.s]
    elif orbital_type == "p":
        orbital_projected_dos = spd_dos[OrbitalType.p]
    elif orbital_type == "d":
        orbital_projected_dos = spd_dos[OrbitalType.d]

    if spin_type == "up":
        orbial_density = orbital_projected_dos.densities[Spin.up]
    else:
        orbial_density = -orbital_projected_dos.densities[Spin.down]

    return orbial_density

def plot_orbital_projected_element_dos(dos_data, desired_orbitals, dst_path, plot_properties):
    if desired_orbitals != ["total"]:
        element_spd_dos_data = get_element_spd_dos(dos_data, desired_orbitals)

    # setup matplotlib plot
    plt.switch_backend('Agg')
    font = {'family': plot_properties["font_family"],
            'size': plot_properties["font_size"]}
    plt.rc('font', **font)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    raw_energy_data = dos_data["tdos"].energies
    if plot_properties["reference_orbital"] == "fermi0":
        energy_data = raw_energy_data - dos_data["efermi"]
    elif plot_properties["reference_orbital"] == "raw":
        energy_data = raw_energy_data
    else:
        orbital = plot_properties["reference_orbital"]
        element = orbital.split("_")[0]
        orbital_type = orbital.split("_")[1][-1]
        orbital_density_data = get_orbital_density_data(element_spd_dos_data[element], orbital_type, plot_properties["spin_type"])
        reference_energy_level = raw_energy_data[np.argmax(orbital_density_data)]
        energy_data = raw_energy_data - reference_energy_level

    for index, orbital in enumerate(desired_orbitals):
        if orbital == "total":
            if plot_properties["spin_type"] == "up":
                ax.plot(energy_data,
                        dos_data["tdos"].densities[Spin.up],
                        color=plot_properties["color_list"][index],
                        label=orbital,
                        lw=plot_properties["line_width"])
            elif plot_properties["spin_type"] == "down":
                ax.plot(energy_data,
                        -dos_data["tdos"].densities[Spin.down],
                        color=plot_properties["color_list"][index],
                        label=orbital,
                        lw=plot_properties["line_width"])
            elif plot_properties["spin_type"] == "both":
                ax.plot(energy_data,
                        dos_data["tdos"].densities[Spin.up],
                        color=plot_properties["color_list"][index],
                        label=orbital,
                        lw=plot_properties["line_width"])
                ax.plot(energy_data,
                        -dos_data["tdos"].densities[Spin.down],
                        color=plot_properties["color_list"][index],
                        lw=plot_properties["line_width"])
            
        else:
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

    x_axis_lims = plot_properties["x_axis_lims"]
    y_axis_lims = plot_properties["y_axis_lims"]
    if len(x_axis_lims):
        ax.set_xlim(x_axis_lims[0], x_axis_lims[1])
    if len(y_axis_lims):
        ax.set_ylim(y_axis_lims[0], y_axis_lims[1])

    if plot_properties["indicate_vbm"] == "yes":
        ylim = ax.get_ylim()
        if plot_properties["reference_orbital"] == "fermi0":
            fermi_data = [0] * 2
        elif plot_properties["reference_orbital"] == "raw":
            fermi_data = [dos_data["efermi"]] * 2
        ax.plot(fermi_data, [ylim[0], ylim[1]], '--', color='black', lw=plot_properties["line_width"])
        ax.set_ylim(ylim[0], ylim[1])
    
    if plot_properties["reference_orbital"] == "fermi0":
        ax.set_xlabel(r"$E - E_f$ (eV)")
    elif plot_properties["reference_orbital"] == "raw":
        ax.set_xlabel("Energy (eV)")
    else:
        ax.set_xlabel(r"$E - E_shift$ (eV)")
    ax.set_ylabel(plot_properties["y_label"])
    ax.set_title(plot_properties["title"])
    ax.legend(prop={'size': plot_properties["legend_font_size"]})
    output_path = str(dst_path / f'{plot_properties["output_file_name"]}.{plot_properties["output_file_type"]}')
    plt.savefig(output_path, format=plot_properties["output_file_type"], dpi=plot_properties["dpi"])
    return None

def plot_site_spd_dos(dos_data, site_index, dst_path, plot_properties):
    if site_index == "total":
        spd_dos_data = dos_data["cdos"].get_spd_dos()
    else:
        site = dos_data["dosrun"].initial_structure[site_index]
        spd_dos_data = dos_data["cdos"].get_site_spd_dos(site)

    # setup matplotlib plot
    plt.switch_backend('Agg')
    font = {'family': plot_properties["font_family"],
            'size': plot_properties["font_size"]}
    plt.rc('font', **font)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if plot_properties["reference_orbital"] == "fermi0":
        energy_data = dos_data["tdos"].energies - dos_data["efermi"]
    elif plot_properties["reference_orbital"] == "raw":
        energy_data = dos_data["tdos"].energies

    base_orbital_types = ["s", "p", "d", "f"]
    num_orbitals = len(spd_dos_data)
    for orbital_index in range(num_orbitals):
        orbital_type = base_orbital_types[orbital_index]
        if plot_properties["spin_type"] == "up":
            spin_type = "up"
            ax.plot(energy_data,
                    get_orbital_density_data(spd_dos_data, orbital_type, spin_type),
                    color=plot_properties["color_list"][orbital_index],
                    label=orbital_type,
                    lw=plot_properties["line_width"])
        elif plot_properties["spin_type"] == "down":
            spin_type = "down"
            ax.plot(energy_data,
                    get_orbital_density_data(spd_dos_data, orbital_type, spin_type),
                    color=plot_properties["color_list"][orbital_index],
                    label=orbital_type,
                    lw=plot_properties["line_width"])
        elif plot_properties["spin_type"] == "both":
            spin_type = "up"
            ax.plot(energy_data,
                    get_orbital_density_data(spd_dos_data, orbital_type, spin_type),
                    color=plot_properties["color_list"][orbital_index],
                    label=orbital_type,
                    lw=plot_properties["line_width"])
            spin_type = "down"
            ax.plot(energy_data,
                    get_orbital_density_data(spd_dos_data, orbital_type, spin_type),
                    color=plot_properties["color_list"][orbital_index],
                    lw=plot_properties["line_width"])
        
    x_axis_lims = plot_properties["x_axis_lims"]
    y_axis_lims = plot_properties["y_axis_lims"]
    if len(x_axis_lims):
        ax.set_xlim(x_axis_lims[0], x_axis_lims[1])
    if len(y_axis_lims):
        ax.set_ylim(y_axis_lims[0], y_axis_lims[1])

    if plot_properties["indicate_vbm"] == "yes":
        ylim = ax.get_ylim()
        if plot_properties["reference_orbital"] == "fermi0":
            fermi_data = [0] * 2
        elif plot_properties["reference_orbital"] == "raw":
            fermi_data = [dos_data["efermi"]] * 2
        ax.plot(fermi_data, [ylim[0], ylim[1]], '--', color='black', lw=plot_properties["line_width"])
        ax.set_ylim(ylim[0], ylim[1])
    
    if plot_properties["reference_orbital"] == "fermi0":
        ax.set_xlabel(r"$E - E_f$ (eV)")
    elif plot_properties["reference_orbital"] == "raw":
        ax.set_xlabel("Energy (eV)")
    else:
        ax.set_xlabel(r"$E - E_shift$ (eV)")
    ax.set_ylabel(plot_properties["y_label"])
    ax.set_title(plot_properties["title"])
    ax.legend(prop={'size': plot_properties["legend_font_size"]})
    output_path = str(dst_path / f'{plot_properties["output_file_name"]}.{plot_properties["output_file_type"]}')
    plt.savefig(output_path, format=plot_properties["output_file_type"], dpi=plot_properties["dpi"])
    return None
