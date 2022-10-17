#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:19:28 2020

@author: marionsauer
"""

#!/usr/bin/python

import math, numpy, sys
from Bio.PDB.PDBParser import PDBParser
from argparse import ArgumentParser
numpy.set_printoptions( threshold = numpy.inf, 
                       formatter={'float': '{: 0.4f}'.format, 
                                  'all':lambda x: '{}\n'.format(x)})

# create the command line parser
parser = ArgumentParser(
    """Assigns a score quantifying changes in a protein ensemble's Cb-Cb distance
(contact) map for each residue as contact proximity 
deviation (Sauer et al., 2020) and relative surface exposure as neighbor vector
(Durham et al. 2009)

Requires at least two PDB files that contain identical number of residues, with 
each aligned position representing congruous positions within an ensemble.

@ Returns an array of size nx4, with n being the number of aligned residues, 
as columns as <residue contact promimity RMSD value> <residue Neighbor Vector> <PDB> <Res #>.
    """
)
parser.add_argument(
    "-l",
    "--pdbs_list",
    nargs="*",
    dest="pdbs",
    default=[sys.stdin],
    help="""A txt file containing a list of PDB file paths. At least two PDB files 
    should be listed within file"""
)
parser.add_argument(
    "-p",
    "--predictors",
    dest="predictors",
    help="""Output file name for residues' contact proximity deviations and Neighbor Vector scores"""
)
# Parse the command line
args = parser.parse_args()
# Parse PDB files
pdb_parser = PDBParser()


def getNeighborWeights(structure):
    """
    Get all Cb-Cb distance vectors within a PDB structure 
    and evaluate neighbor count weights for each distance
    @ Returns neighbor count weights and Cb-Cb distance vectors 
    for a single PDB
    """
    residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
    x = 0
    y = 0
    n = len(residues)
    wts = numpy.zeros( [n,n] )
    dist_uvec = numpy.zeros( [n,n,3] )
    for x in range( len(residues) ):
        cb_distance = 0.0
        # Check for glycine for reference residue
        # If glycine, use virtual C-beta atom coordinate
        if residues[x].get_resname()=='GLY':
            cb1 = numpy.array( getVirtualCbeta( residues[x] ) )
        else:
            cb1 = numpy.array( residues[x]["CB"].get_coord() )
        for y in range( len(residues) ):
            # For when the reference and comparison residue are identical
            if y == x:
                wts[x][y] = 0.0
                # Excluding neighbor weights of each residue to itself,
                # calculate for only half of matrix since diagonal is equivalent
            if y > x:
                # Check for glycine for comparison residue 
                if residues[y].get_resname()=='GLY':
                    cb2 = numpy.array( getVirtualCbeta( residues[y] ) )
                else:
                    cb2 = numpy.array( residues[y]["CB"].get_coord() )
                cb_distance = numpy.linalg.norm(cb2 - cb1)
                if cb_distance <= 4.0:
                    wts[x][y] = 1.0
                    wts[y][x] = 1.0
                if cb_distance > 4.0 and cb_distance <= 12.8:
                    distance_bound_ratio  =  ( ( cb_distance - 4.0 ) / (12.8 - 4.0 ) ) * math.pi
                    sigmoidal_proximity = 0.5 * ( math.cos( distance_bound_ratio ) + 1 )
                    wts[x][y] = sigmoidal_proximity
                    wts[y][x] = sigmoidal_proximity
                else:
                    wts[x][y] = 0.0
                    wts[y][x] = 0.0
                # Find unit vectors for each Cb-Cb distance vector
                dist_vec = numpy.subtract( cb2, cb1 )
                if cb_distance == 0:
                    dist_uvec[x][y] = 0.0
                else:
                    dist_uvec[x][y] = numpy.divide( dist_vec, cb_distance )
                    dist_uvec[y][x] = dist_uvec[x][y]
    return wts, dist_uvec

def countNumResidues(structure):
    residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
    return len(residues)

def getVirtualCbeta( glycine ):
    """
    Gets a glycine's N, C, and CA coordinates as vectors to find a rotation
    matrix that rotates N -120 degrees along the CA-C vector, with the origin
    being the CA coordinates
    """
    n = glycine["N"].get_coord()
    c = glycine["C"].get_coord()
    ca = glycine["CA"].get_coord()
    # Place origin at C-alpha atom
    n = n - ca
    c = c - ca
    # Find rotation matrix that rotates n -120 degrees along the 
    # ca-c vector
    theta = -120.0/180.0 * math.pi
    rot = rotaxis( theta, c )
    # Apply rotation to ca-n vector
    cb_origin = numpy.dot( n, rot )
    return cb_origin + ca

def rotaxis( theta , vec ):
    """Calculate left multiplying rotation matrix
    Taken from Bio.PDB.vectors.rotaxis2m since rotaxis2m vector was being
    assigned as a numpy.ndarray vectors object not as a rotaxis2m function

    Parameters
    ----------
    theta : type float. The rotation angle.
    vec : type array. The rotation axis.

    Returns
    -------
    rot : The rotation matrix, a 3x3 array.

    """
    vec = normalize(vec)
    c = numpy.cos(theta)
    s = numpy.sin(theta)
    t = 1 - c
    x, y, z = numpy.array( vec )
    rot = numpy.zeros((3,3))
    # 1st row
    rot[0, 0] = t * x * x + c
    rot[0, 1] = t * x * y - s * z
    rot[0, 2] = t * x * z + s * y
    # 2nd row
    rot[1, 0] = t * x * y + s * z
    rot[1, 1] = t * y * y + c
    rot[1, 2] = t * y * z - s * x
    # 3rd row
    rot[2, 0] = t * x * z - s * y
    rot[2, 1] = t * y * z + s * x
    rot[2, 2] = t * z * z + c
    return rot

def normalize( vec ):
    norm = numpy.sqrt(sum(vec * vec))
    vec_normalized = vec / norm
    return vec_normalized

def getNeighborVector( vects, wts ):
    """
    Returns NeighborVector calculated as 
    = || [Sum unit vector of cb2-cb1 distance * weight] / Sum of weights ||
    where a value of 0 indicates a highly buried residue 
    and a value of 1 indicates a highly exposed residue
    """

    # Get weighted distance unit vectors and their sums of each residue (numerators) 
    # using only half of matrix
    vec_weighted = numpy.zeros( [vects.shape[0], vects.shape[1], vects.shape[2]] )
    for x in range( vects.shape[0] ):
        for y in range( vects.shape[1] ):
            if y > x:
                for c in range( vects.shape[2] ):
                    coord = vects[x][y][c]
                    wt = wts[x][y]
                    vec_weighted[x][y][c] = coord * wt
                    vec_weighted[y][x][c] = vec_weighted[x][y][c]
    residue_vec_wt_sums = numpy.sum( vec_weighted, axis = 0 )
    residue_wt_sums = numpy.sum( wts, axis = 0 )

    # Norm of the average weighted distance unit vector
    residue_nv_un_norm = numpy.zeros( [residue_vec_wt_sums.shape[0], residue_vec_wt_sums.shape[1]] )
    residue_nv = numpy.zeros( [residue_vec_wt_sums.shape[0]] )
    for x in range( residue_vec_wt_sums.shape[0] ):
        for y in range( residue_vec_wt_sums.shape[1] ):
            residue_nv_un_norm[x][y] = numpy.true_divide( residue_vec_wt_sums[x][y], residue_wt_sums[x] )
        residue_nv[x] = numpy.linalg.norm( residue_nv_un_norm[x] )
    return residue_nv
        

def getContactProximityRMSD( wts_array ):
    """
                Calculates the variation of each CB-CB distance calculated
                within each PDB structure ensemble. All PDBs must be of equal length with 
                congruent, or matching, residues.
                
    Returns: The deviation of all possible CB-CB proximity scores as summed for 
                each residue - meant as a metric to illustrate which residues change their
                local chemical environment to achieve larger conformation rearrangements
                """
    proximity_variation = numpy.var(wts_array, 0, dtype = numpy.float64)
    residue_sums = numpy.sum( proximity_variation, axis = 0, dtype = numpy.float64)
    return numpy.sqrt( residue_sums )

def main():
    # Get length of protein
    models = ""
    n = 0
    for file in args.pdbs:
        models=file
        if file is sys.stdin:
            file = file.fileno()
        with open(file) as p:
            pdb = p.readline().strip()
            structure = pdb_parser.get_structure("apdb", pdb)
            n = countNumResidues( structure )

    # Get number of pdb models
    count = 0
    for file in args.pdbs:
        if file is sys.stdin:
            file = file.fileno()
        with open(file) as p:
            count = len(p.readlines(  ))

    cp_ensemble = numpy.empty( [count,n,n] ) # Stores Cb-Cb weights for ensemble of PDBs
    vect_ensemble = numpy.empty( [count,n,n,3] ) # Stores Cb-Cb distance vectors for ensemble of PDBs
    
    # Read in PDB files from list file to calculate weights + distance vectors
    print("\n")
    print("Calculating Neighbor Weights and distance vectors for ensemble...")
    pdb_list = []
    with open(models) as file:
        p = 0
        for pdb in file:
            pdb = pdb.strip()
            pdb_list.append( pdb )
            structure = pdb_parser.get_structure("apdb", pdb)
            # Check to make sure PDB files have same number of residuess
            # Eventually should check if residues align instead of matching length
            num_res = countNumResidues( structure )
            print(pdb, "has",num_res,"residues")
            if num_res != n:
                print("Aborting - all input PDB files must contain identical number of residues.")
                sys.exit(-1)
            cp_ensemble[p], vect_ensemble[p] = getNeighborWeights( structure )
            p += 1

    print("Calculating Neighbor Vector scores...")
    neighbor_vec = numpy.empty( [count, n] )
    for pdb in range(count):
        neighbor_vec[pdb] = getNeighborVector( vect_ensemble[pdb], cp_ensemble[pdb] )
        
    print("Calculating Neighbor Weight RMSD of each Cb-Cb distance...")
    contact_proximity_rmsd_scores = getContactProximityRMSD( cp_ensemble )
    
    print("Done. Writing output.")
    
    with open(args.predictors, "w") as output_file:
        # Print output with PDB labels
        # Eventually add other predictors to generate matrix
        for pdb in range( len(pdb_list) ):
            residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
            for p in range( 0, n):
                output_file.write( str(contact_proximity_rmsd_scores[p]) + " " + str(neighbor_vec[pdb][p]) + " " + pdb_list[pdb] + " " + str(residues[p].get_id()[1]) + "\n")
                
#invoke main
main()
