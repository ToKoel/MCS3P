import sys
import argparse
import json
import glob
import os

def splitArg(argument):
    if "," in argument:
        splitargs = argument.split(",")
        return splitargs
    else:
        return (0.0,0.0,0.0)

def generate(args):
    arguments = vars(args)
    exchangeConstants = splitArg(arguments["exchangeconstants"])
    latticeParameters = splitArg(arguments["latticepar"])
    
    del arguments["exchangeconstants"]
    del arguments["latticepar"]
    
    arguments["Fe_TT"] = float(exchangeConstants[0])
    arguments["Fe_OO"] = float(exchangeConstants[1])
    arguments["Fe_TO"] = float(exchangeConstants[2])
    arguments["lattice_a"] = float(latticeParameters[0])
    arguments["lattice_b"] = float(latticeParameters[1])
    arguments["lattice_c"] = float(latticeParameters[2])
    
    if arguments["number"] == 0:
        arguments["structure_file"] = glob.glob(os.path.join(arguments["structure_file_path"], "D*_template"))[0]
    else:
        arguments["structure_file"] = glob.glob(os.path.join(arguments["structure_file_path"], "D*_{}".format(arguments["number"])))[0]
    
    with open(os.path.join(arguments["json_file_path"], "settings_{}".format(arguments["number"])), "w") as outfile:
        json.dump(arguments, outfile, indent=4)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--measurement', type=str, default="None")
    parser.add_argument('--output', type=str, default="")
    parser.add_argument('--structure_file_path', type=str, default="")
    parser.add_argument('--json_file_path', type=str, default="")
    parser.add_argument('--num_particles', type=int, default=0)
    parser.add_argument('--particle_size', type=int, default=0)
    parser.add_argument('--meas_field', type=float, default=0.0)
    parser.add_argument('--temperature', type=float, default=0.0)
    parser.add_argument('--dipole', type=str, default="")
    parser.add_argument('--steps', type=float, default=0.0)
    parser.add_argument('--av_steps', type=int, default=0)
    parser.add_argument('--num_or', type=int, default=0)
    parser.add_argument('--cool_field', type=float, default=0.0)
    parser.add_argument('--Tupper', type=float, default=0.0)
    parser.add_argument('--Tlower', type=float, default=0.0)
    parser.add_argument('--T_step', type=float, default=0.0)
    parser.add_argument('--alpha', type=float, default=0.0)
    parser.add_argument('--beta', type=float, default=0.0)
    parser.add_argument('--gamma', type=float, default=0.0)
    parser.add_argument('--macrocell_size', type=float, default=0.0)
    parser.add_argument('--anisotropy_constant', type=float, default=0.0)
    parser.add_argument('--ZFC', action='store_true')
    parser.add_argument('--FC', action='store_true')
    parser.add_argument('--APB', action='store_true')
    parser.add_argument('--Bupper', type=float, default=0.0)
    parser.add_argument('--Blower', type=float, default=0.0)
    parser.add_argument('--Bstep', type=float, default=0.0)
    parser.add_argument('--start_temperature', type=float, default=0.0)
    parser.add_argument('--cooling_steps', type=float, default=0.0)
    parser.add_argument('--sigma', type=float, default=0.0)
    parser.add_argument('--latticepar', type=str, default="")
    parser.add_argument('--exchangeconstants', type=str, default="")
    parser.add_argument('--APB_constant', type=float, default=0.0)
    parser.add_argument('--nearest_neighbour_distance', type=float, default=0.2505)
    parser.add_argument('--number', type=int)
    parser.add_argument('--seed', type=int)
    args = parser.parse_args()
    generate(args)

