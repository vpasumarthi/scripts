#!/usr/bin/env python

import ase.io.vasp


def identify_atom_pairs(src_file_path):
    cell = ase.io.vasp.read_vasp(src_file_path)
    return None

