# README

## Quantum Espresso Hubbard V Computation Script

This script processes a Quantum Espresso (QE) structure file, calculates the nearest oxygen neighbors for ferrous iron (Fe2) and ferric iron (Fe3) atoms considering periodic boundary conditions, and extracts Hubbard V values from provided files.

### Features
- Reads a QE structure file to extract cell parameters and atomic positions.
- Converts atomic positions from fractional (crystal) coordinates to Cartesian coordinates.
- Computes distances between Fe and O atoms while accounting for periodic boundary conditions.
- Identifies the 8 nearest O atoms for each Fe2 and Fe3 atom.
- Reads Hubbard V values from specified files (`Fe2_V.txt` and `Fe3_V.txt`).
- Outputs the results in a formatted `V.txt` file and logs detailed computations in `log.txt`.

## Requirements
- Python 3.x
- NumPy

## Usage
Run the script from the command line with:
```
python script.py --input_file conf.qe --fe2_file Fe2_V.txt --fe3_file Fe3_V.txt --output_file V.txt --log_file log.txt
```

### Command-Line Arguments
| Argument        | Description                          | Default |
|---------------|----------------------------------|---------|
| `--input_file` | Quantum Espresso input file      | `conf.qe` |
| `--fe2_file`   | File containing Fe2 V values    | `Fe2_V.txt` |
| `--fe3_file`   | File containing Fe3 V values    | `Fe3_V.txt` |
| `--output_file` | Output file for computed values | `V.txt` |
| `--log_file`   | Log file for computation details | `log.txt` |

## Output Format
The script generates `V.txt` with the format:
```
V Fe2-3d O-2p Fe_index O_index V_value
V Fe3-3d O-2p Fe_index O_index V_value
```
where `Fe_index` and `O_index` are the atom indices, and `V_value` is the Hubbard V parameter.

The log file (`log.txt`) contains detailed calculations including Fe-O distances:
```
Fe_type    Fe_index    O_index    Distance(A)    V_value(eV)
```

## Notes
- The script assumes periodic boundary conditions when calculating distances.
- Ensure the input structure file follows the Quantum Espresso format.
- Only the first 8 values from `Fe2_V.txt` and `Fe3_V.txt` are used for each Fe atom.

## Example
The inputs and results of an example run can be found in `example/`.
