import re
import numpy as np
import statistics
import argparse
import pathlib
import itertools

import harmonic_analytical

T = 17.8  # Kelvin
nparticles = 3

def read_ipi_output(filename):
    """ Reads an i-PI output file and returns a dictionary with the properties in a tidy order. """

    f = open(filename, "r")

    regex = re.compile(".*column *([0-9]*) *--> ([^ {]*)")

    fields = [];
    cols = []
    for line in f:
        if line[0] == "#":
            match = regex.match(line)
            if match is None:
                print("Malformed comment line: ", line)
                continue  # TODO: was error
            fields.append(match.group(2))
            cols.append(slice(int(match.group(1)) - 1, int(match.group(1))))
        else:
            break  # done with header
    f.close()

    columns = {}
    raw = np.loadtxt(filename)
    for i, c in enumerate(fields):
        columns[c] = raw[:, cols[i]].T
        if columns[c].shape[0] == 1:
            columns[c].shape = columns[c].shape[1]
    return columns


def analyze_mean_energy(infile):
    o = read_ipi_output(infile)
    num_points_to_drop = 100
    avg_kinetic = -statistics.mean(o['virial_fq'][num_points_to_drop:])
    avg_potential = statistics.mean(o['potential'][num_points_to_drop:])
    return avg_kinetic + avg_potential

def harmonic_analytical_energy(nparticles, temp):
    return harmonic_analytical.analytical_energy(nparticles, temp)

def get_average_energy(infiles):
    avg_energy = analyze_mean_energy(infiles[0])
    analytical_energy = harmonic_analytical_energy(nparticles, T)

    print("Avg energy:", avg_energy)
    print("Analytical:", analytical_energy)
    return avg_energy


def main():
    parser = argparse.ArgumentParser(description='analyze boson convergence results')

    parser.add_argument('infile', type=str, nargs='*',
                        help='path to i-PI data.out output file')
    args = parser.parse_args()

    # analyze_mean_energy(args.infile[0])
    get_average_energy(args.infile)


if __name__ == "__main__":
    main()