# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def lennard_jones(x, a, b):
    return a*((b/x)**12 - (b/x)**6)

def main():
    src_path = Path.cwd()
    data_file_name = 'stabilization_energies.dat'
    
    data_file_path = src_path / data_file_name
    input_data = np.loadtxt(data_file_name)
    shell_data = input_data[:, 0]
    w_stabilization_data = input_data[:, 1]
    mo_stabilization_data = input_data[:, 2]
    # tol = 0 suggests to ignore zero-shell index
    # tol = <small value like 1.00E-04 asks to consider zero-shell index energy value at tol value>
    tol = 0
    if tol != 0:
        tol_shell_data = np.copy(shell_data)
        tol_shell_data[0] = tol
        w_params = curve_fit(lennard_jones, tol_shell_data, w_stabilization_data)
    else:
        w_params = curve_fit(lennard_jones, shell_data[1:], w_stabilization_data[1:])
    plt.switch_backend('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title(f'Curve fitting of shell-wise stabilization energy around W dopant')
    ax.plot(shell_data, w_stabilization_data, 'o', mfc='#607c8e', mec='black')
    if tol != 0:
        ax.plot(shell_data, lennard_jones(tol_shell_data, w_params[0][0], w_params[0][1]), 'b-')
    else:
        ax.plot(shell_data[1:], lennard_jones(shell_data[1:], w_params[0][0], w_params[0][1]), 'b-')
    ax.set_xlabel('Shell number')
    ax.set_xlabel('Stabilization energy (eV)')
    ax.text(7, 0.5, f'eps={w_params[0][0]:.4f}, rm={w_params[0][1]:.4f}')
    figure_name = f'W-dopant_curve_fitting_stabilization.png'
    figure_path = src_path / figure_name
    plt.savefig(str(figure_path))

    if tol != 0:
        mo_params = curve_fit(lennard_jones, tol_shell_data, mo_stabilization_data)
    else:
        mo_params = curve_fit(lennard_jones, shell_data[1:], mo_stabilization_data[1:])
    plt.switch_backend('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title(f'Curve fitting of shell-wise stabilization energy around Mo dopant')
    ax.plot(shell_data, mo_stabilization_data, 'o', mfc='#607c8e', mec='black')
    if tol != 0:
        ax.plot(shell_data, lennard_jones(tol_shell_data, mo_params[0][0], mo_params[0][1]), 'b-')
    else:
        ax.plot(shell_data[1:], lennard_jones(shell_data[1:], mo_params[0][0], mo_params[0][1]), 'b-')
    ax.set_xlabel('Shell number')
    ax.set_xlabel('Stabilization energy (eV)')
    ax.text(7, 0.06, f'eps={mo_params[0][0]:.4f}, rm={mo_params[0][1]:.4f}')
    figure_name = f'Mo-dopant_curve_fitting_stabilization.png'
    figure_path = src_path / figure_name
    plt.savefig(str(figure_path))
    
if __name__ == '__main__':
    main()
