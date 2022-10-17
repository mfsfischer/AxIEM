#!/usr/bin/env python3
"""
Performs leave-one-out analyses of epitope mapping models using a range of epitope radii
@author marionfischer
"""
# Parsing files
from argparse import ArgumentParser
from Bio.PDB.PDBParser import PDBParser
from Bio import BiopythonWarning

# Basic calculations and statistical analysis
from collections import Counter, defaultdict
from itertools import chain
import math
import numpy
import pandas
import sys

# Mute Biopython PDB parsing warnings when reading Rosetta output files
import warnings
warnings.simplefilter('ignore', BiopythonWarning)

# Create command line parser
parser = ArgumentParser(
    """Step 1 in AxIEM benchmark:
    Calculates per-residue features CP (contact proximity) and NV (Neighbor Vector) using at least two PDB structures
    of identical protein length. PDBs should contain similar/same virus protein sequence.
    """
)
parser.add_argument(
    "--data",
    dest="data",
    required=True,
    help="""INPUT text file - see documentation for formatting"""
)
parser.add_argument(
    "--features",
    dest="features",
    help="""OUTPUT file to write feature sets - see documentation for formatting"""
)

# Parse the command line
args = parser.parse_args()
# Parse PDB file
pdbparser = PDBParser()

def getPDBdata_byVirus(dataset):
    """
    Parameters
    ----------
    dataset : txt file of annotated data (e.g. AxIEM.data)

    Returns
    -------
    viruses_cb_distances : dict of dict containing numpy arrays size (n,n) of each virus's PDB's contact map distances 
    viruses_CP_var       : dict of dict containing numpy arrays size (n,n) of each virus's PDB's contact proximity variation scores
    viruses_NV           : dict of dict containing numpy arrays size (n,n) of each virus's PDB's neighbor vector scores
    proteins : list of str; virus ids within dataset
    pdb_lens : list of int; number of residues within each PDB for each virus
    """
    data = pandas.read_csv(dataset, delimiter=' ')
    viruses = data.iloc[:,0].to_list()
    pdbs = data.iloc[:,1].to_list()
    classifiers = data.iloc[:,2].to_list()
    resno = data.iloc[:,3].to_list()
    reu = data.iloc[:,4].to_list()
    viruses_cb_distances = defaultdict(lambda : defaultdict([]))
    viruses_CP = defaultdict(lambda : defaultdict([]))
    viruses_NV = defaultdict(lambda : defaultdict([]))
    features = []
    pdblen = 0
    for i in range(len(viruses)):
        if not viruses[i] in viruses_cb_distances:
            viruses_cb_distances[viruses[i]] = {}
            viruses_CP[viruses[i]] = {}
            viruses_NV[viruses[i]] = {}
            if not pdbs[i] in viruses_cb_distances[viruses[i]]:
                pdb = "pdb_structures/{}".format(pdbs[i])
                structure = pdbparser.get_structure("apdb", pdb)
                residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
                print("{} structure {} with {} residues".format(viruses[i], pdbs[i], len(residues)))
                pdblen = len(residues)
                dist, cb_wts, cb_uvec = getCbInfo(residues)
                viruses_cb_distances[viruses[i]][pdbs[i]] = dist
                viruses_CP[viruses[i]][pdbs[i]] = cb_wts
                viruses_NV[viruses[i]][pdbs[i]] = getNeighborVector( cb_wts, cb_uvec)
                
        else:
            if not pdbs[i] in viruses_cb_distances[viruses[i]]:
                pdb = "pdb_structures/{}".format(pdbs[i])
                structure = pdbparser.get_structure("apdb", pdb)
                residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
                print("{} structure {} with {} residues".format(viruses[i], pdbs[i], len(residues)))
                if len(residues)!=pdblen:
                    print("Aborting... all PDB files for the {} fusion protein must contain identical number of residues".format(viruses[i]))
                    sys.exit(-1)
                dist, cb_wts, cb_uvec = getCbInfo(residues)
                viruses_cb_distances[viruses[i]][pdbs[i]] = dist
                viruses_CP[viruses[i]][pdbs[i]] = cb_wts
                viruses_NV[viruses[i]][pdbs[i]] = getNeighborVector(cb_wts, cb_uvec)
    viruses_CP_var = defaultdict(lambda : [])
    for virus, pdb in viruses_CP.items():
        virus_cp = []
        i = 0
        for p, cp in pdb.items():
            virus_cp.append(numpy.array(cp))
            i += 1
        if not virus in viruses_CP_var:
            viruses_CP_var[virus] = getContactProximityVariation(numpy.asarray(virus_cp))
    return viruses_cb_distances, viruses_CP_var, viruses_NV, classifiers, resno, reu

def getCbInfo(residues):
    """
    Parameters
    ----------
    residues : list of len(n) of Biopython PDB structure's residue entities

    Returns
    -------
    cb_distances : numpy array size (n,n) of all C-beta atoms' pairwise distances
    cb_wts       : numpy array size (n,n) of single PDB's contact map probability weights
                   with contact likely at or less than 4.0 Angstroms, not likely greater
                   than 12.8 Angstroms, and with a cosine-weighted probability between
                   4.0 and 12.8 Angstroms. (Durham et al., 2009; Sauer et al., 2020)
    cb_uvec      : numpy array size (n,n) of each Cb-Cb distance unit vector
    """
    u = 12.8
    l = 4.0
    cb_distances = numpy.zeros((len(residues),len(residues)))
    cb_coords = numpy.zeros((len(residues),3))
    cb_wts = numpy.zeros((len(residues),len(residues)))
    cb_uvec = numpy.zeros((len(residues),len(residues),3))
    for x in range(len(residues)):
        cb1 = []
        if residues[x].get_resname()=='GLY':
            cb1 = numpy.array( getVirtualCbeta(residues[x]) )
        else:
            cb1 = numpy.array( residues[x]["CB"].get_coord() )
        cb_coords[x] = cb1
        for y in range(len(residues)):
            # Cb contact map is symmetric - only compute half of matrix
            cb2 = []
            if y == x:
                cb_distances[x][y] = 0.0
                cb_wts[x][y] = 1.0
                cb_uvec[x][y] = 0.0
            if y > x:
                if residues[y].get_resname()=='GLY':
                    cb2 = numpy.array( getVirtualCbeta(residues[y]) )
                else:
                    cb2 = numpy.array( residues[y]["CB"].get_coord() )
                cb_distances[x][y] = numpy.linalg.norm(cb2 - cb1)
                cb_distances[y][x] = cb_distances[x][y]
                if cb_distances[x][y] <= l:
                    cb_wts[x][y] = 1.0
                    cb_wts[y][x] = 1.0
                if cb_distances[x][y] > l and cb_distances[x][y] <= u:
                    distance_bound_ratio  =  ( ( cb_distances[x][y] - l ) / (u - l ) ) * math.pi
                    sigmoidal_proximity = 0.5 * ( math.cos( distance_bound_ratio ) + 1 )
                    cb_wts[x][y] = sigmoidal_proximity
                    cb_wts[y][x] = sigmoidal_proximity
                else:
                    cb_wts[x][y] = 0.0
                    cb_wts[y][x] = 0.0
                vec = numpy.subtract(cb2, cb1)
                cb_uvec[x][y] = numpy.divide(vec, cb_distances[x][y])
                cb_uvec[y][x] = cb_uvec[x][y]
    return cb_distances, cb_wts, cb_uvec

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
    assigned as a numpy.ndarray vectors object not as a rotaxis2m parameter

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

def getNeighborVector( wts, uvec ):
    """
    Returns NeighborVector calculated as 
    = || Sum [unit vector of cb2-cb1 distance * cb2-cb1 weight] / Sum [cb2-cb1 weights] ||
    where a value of 0 indicates a highly buried residue 
    and a value of 1 indicates a highly exposed residue

    Parameters
    ----------
    wts   : numpy array size (n,n) of Cb-Cb contact map probabilities
    vects : numpy array size (n,n,3) of Cb-Cb unit vectors

    Returns
    -------
    nv : numpy array size (n) of per-residue neighbor vector feature values
    """

    # Get weighted distance unit vectors and their sums of each residue (numerators) 
    # using only half of matrix
    vec_weighted = numpy.zeros( [uvec.shape[0], uvec.shape[1], uvec.shape[2]] )
    for x in range( uvec.shape[0] ):
        for y in range( uvec.shape[1] ):
            if y > x:
                for c in range( uvec.shape[2] ):
                    coord = uvec[x][y][c]
                    vec_weighted[x][y][c] = coord * wts[x][y]
                    vec_weighted[y][x][c] = vec_weighted[x][y][c]
    vec_wt_sums = numpy.sum( vec_weighted, axis = 0 )
    wt_sums = numpy.sum( wts, axis = 0 )

    # Norm of the average weighted distance unit vector
    nv_un_norm = numpy.zeros( [vec_wt_sums.shape[0], vec_wt_sums.shape[1]] )
    nv = numpy.zeros( [vec_wt_sums.shape[0]] )
    for x in range( vec_wt_sums.shape[0] ):
        for y in range( vec_wt_sums.shape[1] ):
            nv_un_norm[x][y] = numpy.true_divide( vec_wt_sums[x][y], wt_sums[x] )
            nv[x] = numpy.linalg.norm( nv_un_norm[x] )
    return nv

def getContactProximityVariation( wts_array ):
    """
    Calculates the variation of each CB-CB distance calculated
    within each PDB structure ensemble. All PDBs must be of equal length with 
    congruent, or matching, residues.

    Parameters
    ----------
    wts_array : numpy array size (i,n,n) of i number of PDB contact proximity weights
                
    Returns
    -------
    residue_sums : The summed variation of all possible CB-CB proximity scores for each residue 
                - meant as a metric to illustrate which residues change their local
                chemical environment to achieve larger conformation rearrangements
    """
    proximity_variation = numpy.var(wts_array, 0, dtype = numpy.float64)
    residue_sums = numpy.sum( proximity_variation, axis = 0, dtype = numpy.float64)
    return residue_sums


def main():

    # -------Set up input and output-------
    print("\nCalculating distance maps, contact proximity variation, and neighbor vector feature values for:")
    
    # -------Set up initial feature calculations--------
    # Get contact map distances to shorten NS_u calculations
    # Calculate contact proximity (CP) and neighbor vector (NV) features for each protein ensemble and protein/PDB respectively
    CBdist_byVirusPDB, CPvar_byVirus, NV_byVirusPDB, classifiers, resno, reu = getPDBdata_byVirus(args.data)

    # -------Write output for Step 2--------
    print("Writing contact maps as 'pdb_contact_maps/PDB.pdbmaps'")
    for virus, pdb in CBdist_byVirusPDB.items():
        for p, dists in pdb.items():
            with open('{}maps'.format(str(p)), 'w') as maps:
                for i in range(len(dists)):
                    for j in range(len(dists)):
                        if j==len(dists)-1:
                            maps.write(str(dists[i][j])+"\n")
                        else:
                            maps.write(str(dists[i][j])+" ")
    print("Writing features to {}".format(args.features))
    with open(args.features, 'w') as feat:
        feat.write("Virus PDB Classifier ResNo REU CPvar NV\n")
        j = 0
        for v1, v2 in zip(CPvar_byVirus, NV_byVirusPDB):
            cp = CPvar_byVirus[v1]
            for pdb, nv in NV_byVirusPDB[v2].items():
                for i in range(len(nv)):
                    line =[v1, pdb, classifiers[j], resno[j], reu[j], cp[i], nv[i]]
                    feat.write(" ".join(map(str, line)))
                    feat.write("\n")
                    j += 1
        
    print("Done")

if __name__ == "__main__":
    main()
