#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


def plot_energy_diagram(input_data_file_path, color_info, marker_info,
                        linestyle_info, marker_size, font_size, shift):
    plt.switch_backend('Agg')
    energy_data = np.loadtxt(input_data_file_path)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    num_data = len(energy_data)
    rc_data = np.linspace(0, 1, num_data)
    ax.plot(rc_data[0:int(num_data / 2) + 1 + (shift if shift < 0 else 0)],
            energy_data[0:int(num_data / 2) + 1 + (shift if shift < 0 else 0)],
            c=color_info, marker=marker_info, ls=linestyle_info,
            markersize=marker_size)
    ax.plot(rc_data[int(num_data / 2) + (shift if shift > 0 else 0):num_data],
            energy_data[int(num_data / 2) + (shift if shift > 0 else 0):num_data],
            c=color_info, marker=marker_info, ls=linestyle_info,
            markersize=marker_size)
    ax.set_xlabel('Reaction Coordinate', fontsize=font_size)
    ax.set_ylabel('Energy (eV)', fontsize=font_size)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.show()
    figure_path = input_data_file_path.with_suffix('.png')
    plt.savefig(str(figure_path))
