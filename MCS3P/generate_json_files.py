import sys
import argparse
import json

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
    
    arguments["Fe_TT"] = exchangeConstants[0]
    arguments["Fe_OO"] = exchangeConstants[1]
    arguments["Fe_TO"] = exchangeConstants[2]
    arguments["lattice_a"] = latticeParameters[0]
    arguments["lattice_b"] = latticeParameters[1]
    arguments["lattice_c"] = latticeParameters[2]
    
    with open("json_test.json", "w") as outfile:
        json.dump(arguments, outfile, indent=4)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--measurement', type=str, default="None")
    parser.add_argument('--output', type=str, default="")
    parser.add_argument('--structure_path', type=str, default="")
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
    args = parser.parse_args()
    generate(args)

