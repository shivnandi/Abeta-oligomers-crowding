# Amyloid-beta-oligomers-study-with-Molecular-Crowders


# README for Research Code Repository

## Overview
This GitHub repository contains all the codes used to calculate the properties in the study titled:

"Impact of Crowder Size, Hydrophobicity, and Hydration on the Structure of Amyloid-β Oligomers"

The repository includes scripts and programs for computing key molecular properties from molecular dynamics (MD) simulations.

## Included Codes
1. **Tetrahedral Order Calculation** ('qtet.f90' and related scripts)
   - Computes the tetrahedral order parameter for water molecules in the system.
2. **Mean Squared Displacement (MSD) Calculation** ('msd.final')
   - Calculates the MSD of water molecules in MD simulations.
3. **Preferential Binding Analysis** ('pref-binding.f90')
   - Determines the preferential binding coefficient of crowders around the Aβ(16-22) oligomers.
4. **Secondary Structure Analysis** ('secondary_structure.f90')
   - Evaluates the secondary structure distribution based on data from VMD's STRIDE algorithm.

## Additional Files and Scripts
- **Bash Scripts** ('bash_qtet.sh')
  - Automate the execution of Fortran programs and pre-process input files.


## Compilation and Execution
All Fortran codes can be compiled using 'gfortran', e.g.,

    gfortran filename.f90 -o output_name
    ./output_name

Ensure required input files are present before execution.

## License
This repository is released under the **MIT License**. You are free to use, modify, and distribute the code with proper attribution.

## Citation  
If you use these codes in your research, please cite the study:  
> **"Impact of Crowder Size, Hydrophobicity, and Hydration on the Structure of Amyloid-β Oligomers"**  
> *Shivnandi, Divya Nayar*  

For any issues or contributions, feel free to submit a pull request or raise an issue in the repository.  



