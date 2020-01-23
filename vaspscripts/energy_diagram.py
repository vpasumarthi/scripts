# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def plot_energy_diagram(input_data_file_path, color_info, marker_info,
                        linestyle_info, marker_size, font_size, shift,
                        num_forecast_periods, cubic_spline):
    plt.switch_backend('Agg')
    energy_data = np.loadtxt(input_data_file_path)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    num_data = len(energy_data)
    rc_data = np.linspace(0, 1, num_data)
    if cubic_spline:
        spline_data = interp1d(rc_data, energy_data, kind='cubic')
        x_new = np.linspace(rc_data[0], rc_data[-1], num=41, endpoint=True)
        ax.plot(x_new, spline_data(x_new), c=color_info, ls=linestyle_info)
        ax.plot(rc_data,
                energy_data,
                c=color_info, marker=marker_info, ls='None', markersize=marker_size)
        
    if not cubic_spline:
        ax.plot(rc_data[0:int(num_data / 2) + 1 + (shift if shift < 0 else 0)],
                energy_data[0:int(num_data / 2) + 1 + (shift if shift < 0 else 0)],
                c=color_info, marker=marker_info, ls=linestyle_info,
                markersize=marker_size)

    if num_forecast_periods > 0:
        poly = np.polyfit(rc_data[0:int(num_data / 2) + 1 + (shift if shift < 0 else 0)],
                          energy_data[0:int(num_data / 2) + 1 + (shift if shift < 0 else 0)], deg=2)
        y_int  = np.polyval(poly, rc_data[0:int(num_data / 2) + 1 + (shift if shift < 0 else 0) + num_forecast_periods])
        ax.plot(rc_data[0:int(num_data / 2) + 1 + (shift if shift < 0 else 0) + num_forecast_periods],
                y_int, c=color_info, marker=marker_info, ls='--',
                markersize=marker_size)
    if not cubic_spline:
        ax.plot(rc_data[int(num_data / 2) + (shift if shift > 0 else 0):num_data],
                energy_data[int(num_data / 2) + (shift if shift > 0 else 0):num_data],
                c=color_info, marker=marker_info, ls=linestyle_info,
                markersize=marker_size)
    if num_forecast_periods < 0:
        poly = np.polyfit(rc_data[int(num_data / 2) + (shift if shift > 0 else 0):num_data],
                          energy_data[int(num_data / 2) + (shift if shift > 0 else 0):num_data], deg=2)
        y_int  = np.polyval(poly, rc_data[int(num_data / 2) + (shift if shift > 0 else 0) + num_forecast_periods:num_data])
        ax.plot(rc_data[int(num_data / 2) + (shift if shift > 0 else 0) + num_forecast_periods:num_data],
                y_int, c=color_info, marker=marker_info, ls='--',
                markersize=marker_size)
    ax.set_xlabel('Reaction Coordinate', fontsize=font_size)
    ax.set_ylabel('Energy (eV)', fontsize=font_size)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.show()
    figure_path = input_data_file_path.with_suffix('.png')
    plt.savefig(str(figure_path), dpi=600)
