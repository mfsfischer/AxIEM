#!/usr/bin/env python3
"""Inserts b-factor values into a PDB file.

b-factor values are provided through a file with chain id and sequence id.
A copy of the PDB file containing the b-factor values is written out.
"""

from argparse import ArgumentParser
from collections import defaultdict

# create the command line parser
parser = ArgumentParser(
    """This script inserts b-factor values into a PDB file. The PDB file name
    and the name of the file containing the b-factor values need to be given.
    An optional default value can be provided, which will be assigned to
    residues for which no b-factor value was provided."""
)
parser.add_argument(
    "-p",
    "--pdb",
    dest="pdb_file",
    required=True,
    help="""The PDB file the b-factor values shall be inserted into. This file
    won\"t be changed by this script."""
)
parser.add_argument(
    "-b",
    "--b_factors",
    dest="b_factors_file",
    required=True,
    help="""A file containing the b-factor values that should be inserted into
    "the PDB file. The file format is expected as
    <chain_id> <sequence_id> <b-factor>."""
)
parser.add_argument(
    "-o",
    "--output_file",
    dest="output_file",
    required=True,
    help="""The output PDB file produced by this script. It will be identical
    to the input PDB file but will contain the provided b-factors."""
)
parser.add_argument(
    "-d",
    "--default_value",
    dest="default_value",
    default="NaN",
    help="""The default b-factor value used if no b-factor value was given for
    certain residues."""
)

# parse the command line
args = parser.parse_args()


def create_bfac_dict(bfactor_file_name):
    """Parses the provided b-factor file and returns a dictionary that maps
    each b-factor value to a chain id and a sequence id.

    @param bfactor_file_name: The name of file containing the b-factor values.
    """

    # initialize the dictionary, which is a dictionary of dictionaries
    bfac_dict = defaultdict(dict)

    # parse the b-factor file
    with open(bfactor_file_name) as bfac_file:
        for line in bfac_file:
            line = line.split()
            chain_id = line[0]
            seq_id = line[1]
            bfac = float(line[2])
            bfac_dict[chain_id][seq_id] = bfac

    return bfac_dict


def insert_bfac_into_pdb(pdb, bfac_dict, default_value):
    """Inserts the b-factor values provided in the dictionary into the provided
    PDB string and returns a copy of the PDB string with the inserted b-factor
    values. If no b-factor value was provided for certain residues,
    the provided default value is inserted.

    @param pdb: The PDB into which the b-factor values should be inserted.

    @param bfac_dict: Dictionary to map residues to b-factor values.

    @param default_value: Default b-factor value for residues, which are not in
    the dictionary.
    """

    new_pdb = []
    for line in pdb:
        if line[0:6] == "ATOM  ":
            chain_id = line[20:22].strip()
            seq_id = line[23:26].strip()
            if chain_id in bfac_dict.keys() and seq_id in bfac_dict[chain_id].keys():
                new_pdb.append("%s%6.2F%s" % (line[:60], bfac_dict[chain_id][seq_id],
                                              line[66:]))
            else:
                new_pdb.append("%s%6s%s" % (line[:60], default_value, line[66:]))
        elif line[0:6] == "HETATM":
            new_pdb.append("%s%6s%s" % (line[:60], default_value, line[66:]))
        else:
            new_pdb.append(line)

    return new_pdb


def main():
    """Reads in the PDB file and writes out a PDB file that contains the provided
    b-factor values.
    """

    # create the dictionary to map b-factor values to residues
    dictionary = create_bfac_dict(args.b_factors_file)

    # read in the provided PDB file
    with open(args.pdb_file) as pdb_file:
        pdb = pdb_file.readlines()

    # insert the b-factor values into the PDB string and write out the PDB file
    insert_bfac_into_pdb(pdb, dictionary, args.default_value)
    new_pdb = "".join(insert_bfac_into_pdb(pdb, dictionary, args.default_value))
    with open(args.output_file, "w") as output_file:
        output_file.write(new_pdb)

# call main function
main()
