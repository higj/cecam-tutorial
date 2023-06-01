import re
import math
import numpy as np
import statistics
import argparse
import os
import pathlib
import itertools
#from create_converged_job import job_params_from_id

T = 17.4  ## JACOB THE KING
nparticles = 3 ## JACOB ADD TO ARGUMENTS
infile = 'C:/Users/Jacob/Desktop/bosons_converged_3_32_17.8_1/#data.out#0#'
def read_ipi_output(filename):
    """ Reads an i-PI output file and returns a dictionary with the properties in a tidy order. """
    
    f = open(filename, "r")
    
    regex = re.compile(".*column *([0-9]*) *--> ([^ {]*)")
    
    fields = []; cols = []
    for line in f:
        if line[0] == "#":
            match = regex.match(line)
            if match is None:
                print("Malformed comment line: ", line)
                continue # TODO: was error
            fields.append(match.group(2))
            cols.append(slice(int(match.group(1))-1,int(match.group(1))))
        else:
            break # done with header
    f.close()
    
    columns = {}
    raw = np.loadtxt(filename)
    for i, c in enumerate(fields):
        columns[c] = raw[:,cols[i]].T
        if columns[c].shape[0] == 1:
            columns[c].shape = columns[c].shape[1]
    return columns

def analyze_mean_energy(infile):
    o = read_ipi_output(infile)
    num_points_to_drop = 2001
    kin_est = -o['virial_fq'][num_points_to_drop:]
    pot_est = o['potential'][num_points_to_drop:]

    avg_kinetic = -statistics.mean(o['virial_fq'][num_points_to_drop:])
    avg_potential = statistics.mean(o['potential'][num_points_to_drop:])
    return avg_kinetic + avg_potential, kin_est, pot_est



def job_id_from_filepath(infile):
    return pathlib.Path(infile).parts[-2]

def groupby(lst, key_func):
    return itertools.groupby(sorted(lst, key=key_func), key_func)

def group_and_analyze_files(infiles):
    avg_energy,  kin_est, pot_est = analyze_mean_energy(infiles[0])

    print("Avg energy:", avg_energy)
    return avg_energy,  kin_est, pot_est

def main():
    parser = argparse.ArgumentParser(description='analyze boson convergence results')

    parser.add_argument('infile', type=str, nargs='*',
                        help='path to i-PI data.out output file')
    args = parser.parse_args()
    
    #analyze_mean_energy(args.infile[0])
    avg_energy,  kin_est, pot_est = group_and_analyze_files([infile])
    return  avg_energy,  kin_est, pot_est



if __name__ == "__main__":
    main()
