import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Compute Hubbard V values in a Quantum Espresso structure file.")
parser.add_argument("--input_file", "-i", type=str, default="conf.qe")
parser.add_argument("--fe2_file", type=str, default="Fe2_V.txt")
parser.add_argument("--fe3_file", type=str, default="Fe3_V.txt")
parser.add_argument("--output_file", "-o", type=str, default="V.txt")
parser.add_argument("--log_file", type=str, default="log.txt")

args = parser.parse_args()

def read_qe_structure(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    cell_start = lines.index("CELL_PARAMETERS (angstrom)\n") + 1
    cell = np.array([list(map(float, lines[cell_start + i].split())) for i in range(3)])
    
    atom_start = lines.index("ATOMIC_POSITIONS (crystal)\n") + 1
    atoms = []
    for i in range(atom_start, len(lines)):
        parts = lines[i].split()
        if len(parts) < 4:
            break
        atoms.append((parts[0], np.dot(np.array(list(map(float, parts[1:]))), cell)))
    
    return cell, atoms

def distance_PBC(atom1_pos, atom2_pos, cell):
    """ Calculate the minimum distance between two atoms in a periodic system """
    min_dist = float('inf')
    for shift in np.array([[i, j, k] for i in [-1, 0, 1] for j in [-1, 0, 1] for k in [-1, 0, 1]]):
        shifted_atom1_pos = atom1_pos + np.dot(shift, cell)
        dist = np.linalg.norm(atom2_pos - shifted_atom1_pos)
        min_dist = min(min_dist, dist)
    return min_dist

def find_nearest_oxygen(fe_position, o_atoms, cell):
    """ Find the 8 nearest oxygen atoms to a given iron atom """
    distances = []
    for o_index, o_symbol, o_position in o_atoms:
        dist = distance_PBC(fe_position, o_position, cell)
        distances.append((o_index, dist))
    
    distances.sort(key=lambda x: x[1])
    return distances[:8]

def read_v_values(file_path):
    """ Read the Hubbard V values from a file """
    with open(file_path, 'r') as f:
        values = [float(line.split()[-1]) for line in f.readlines()[:8]]
    return values

def process_structure(structure_file):
    """ Process the structure file and write the Hubbard V values """
    cell, atoms = read_qe_structure(structure_file)
    
    fe2_atoms = [(i+1, symbol, pos) for i, (symbol, pos) in enumerate(atoms) if symbol == 'Fe2']
    fe3_atoms = [(i+1, symbol, pos) for i, (symbol, pos) in enumerate(atoms) if symbol == 'Fe3']
    o_atoms = [(i+1, symbol, pos) for i, (symbol, pos) in enumerate(atoms) if symbol == 'O']
    
    v_file = open(args.output_file, "w")
    log_file = open(args.log_file, "w")
    log_file.write("\tFe_i\tO_i\td (A)\tV (eV)\n")

    for fe_type, fe_atoms, v_file_name in [("Fe2", fe2_atoms, args.fe2_file), ("Fe3", fe3_atoms, args.fe3_file)]:
        v_values = read_v_values(v_file_name)
        
        for fe_index, fe_symbol, fe_position in fe_atoms:
            nearest_oxygens = find_nearest_oxygen(fe_position, o_atoms, cell)
            
            for (o_index, distance), v_value in zip(nearest_oxygens, v_values):
                v_file.write(f"V {fe_type}-3d O-2p {fe_index} {o_index} {v_value:.2f}\n")
                log_file.write(f"{fe_type}\t{fe_index}\t{o_index}\t{distance:.3f}\t{v_value:.2f}\n")
    
    v_file.close()
    log_file.close()

process_structure(args.input_file)
