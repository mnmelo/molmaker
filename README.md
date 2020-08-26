molmaker.py version 1.5.0-26-08-2020 by Manuel Melo

This tool creates a .gro from an .itp file. It works by randomly scattering coordinates along a linear stretch and then performing an evil minimization as VdW and charges are faded in (using the free energy code). As you already guess, it's totally useless for proteins unless you want a linear segment (in which case it works pretty well!).

Additionally, molmaker.py will likely not preserve your chiral centers unless you protect them in your topology using some sort of dihedral potential/restraint. Alternatively you might want to hand-correct each center using other tools and then energy-minimizing.

Check the -h flag for more details. Please report bugs in the GitHub project (https://github.com/mnmelo/molmaker).
