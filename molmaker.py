#!/usr/bin/env python3

###  Copyright 2016-2020 Manuel N. Melo ###
#########################################################################
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################


## Template for the .mdp file #############
#   edit to suit your needs               #
#   but leave the %s fields alone!        #
#   (you can set them with option flags)  #
###########################################

mdp_template = """
define                   = -DFLEXIBLE
integrator               = %s
emtol                    = 0.1
nsteps                   = %s
nstxout                  = 1
nstcgsteep               = 50
nstlist                  = 1
ns-type                  = grid 
pbc                      = no
epsilon_r                = %s
lincs_order              = 8
lincs_iter               = 2
continuation             = yes
free-energy              = yes
init-lambda              = 0
delta-lambda             = 0.002
couple-moltype           = system
couple-lambda0           = none
couple-lambda1           = vdw-q
couple-intramol          = yes
sc-power                 = 2
sc-alpha                 = 1.5
sc-sigma                 = 0.2
"""

mdp_template_gmx5 = """
; The following line(s) will be appended if running gmx >= 5.0
cutoff-scheme            = group
"""

## OK, no touchy from here onwards ########
###########################################

import sys
import os
import re
import subprocess
import random
import argparse
import shutil
from pathlib import Path
import copy

# Functions
############

def find_exec(names, error='error', extra=''):
    if isinstance(names, str):
        names = (names,)
    for name in names:
        proc = subprocess.Popen(['which', name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        path = proc.communicate()[0].decode('utf-8')
        if path:
            path = Path(path.strip()).resolve()
            return path
    else:
        if error in ('error', 'warning'):
            namelist = ', '.join(names)
            msg = f'{error}: can\'t find {namelist} in $PATH. {extra}'
            if error == 'error':
                sys.exit(msg)
            print(msg, file=sys.stderr)
        return False

def convert_val(val):
    try:
        try:
            return int(val)
        except ValueError:
            return float(val)
    except ValueError:
        return str(val)

def clean(line):
    line = line.split(';')[0]
    line = line.strip()
    if not line:
        return None
    line = line.strip('[]')
    line = line.strip()
    return [convert_val(val) for val in line.split()]


# Classes
##########

class ProperFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """A hackish class to get proper help format from argparse.

    """
    def __init__(self, *args, **kwargs):
        super(ProperFormatter, self).__init__(*args, **kwargs)

class ITP:
    def __init__(self, val_lists):
        self.topology = {}
        self._parse(val_lists)

    @property
    def name(self):
        return self.topology['moleculetype'][0][0]

    @name.setter
    def name(self, newname):
        self.topology['moleculetype'][0][0] = str(newname)

    def _parse(self, val_lists):
        if val_lists[0][0] != 'moleculetype':
            raise TypeError('Topology doesn\'t start with "moleculetype"')
        curr_section = None
        for vals in val_lists:
            if len(vals) == 1:
                self.topology[vals[0]] = []
                curr_section = self.topology[vals[0]]
            else:
                curr_section.append(vals)
        # Make the order predictable
        for vals in self.topology.values():
            vals[:] = sorted(vals)


class MolMaker(argparse.ArgumentParser):
    def __init__(self, formatter_class=ProperFormatter, *args, **kwargs):
        argparse.ArgumentParser.__init__(self, *args, formatter_class=formatter_class, **kwargs)

    def checkopts(self, args=sys.argv):
        self.values = self.parse_args(args)
        #
        # Input options checking and assignment of defaults not done by argparse.
        self.itp = self.values.itp
        if not self.itp:
            self.itp = Path(input('Topology file: '))

        if not self.itp.exists():
            sys.exit('Error: can\'t find the specified topology file.')

        self.gro = self.values.gro
        if not self.gro:
            self.gro = self.itp.with_suffix('.gro')
        self.name = self.gro.stem
        if not self.values.tmpprefix:
            sys.exit('Error: -temp-prefix cannot be an empty string.')
        self.basename = self.gro.with_name(self.values.tmpprefix + self.name)

        self.values.fuzzy = max(self.values.fuzzy, 0)

        # Are we running gromacs >=5?
        self.gmx = find_exec(('gmx', 'gmx_d', 'gmx_mpi'), error=None)

        if self.gmx:
            self.grompp = [self.gmx, "grompp"]
            self.mdrun = [self.gmx, "mdrun"]
            self.trjconv = [self.gmx, "trjconv"]
            self.editconf = [self.gmx, "editconf"]
        else:
            # No gmx executable. Let's try the gromacs 4 names.
            self.grompp = [find_exec('grompp')]
            self.mdrun = [find_exec('mdrun')]
            self.trjconv = [find_exec('trjconv', error='warning', extra='Will not center output coordinates nor output minimization trajectory (if requested).')]
            self.editconf = [find_exec('editconf', error='warning', extra='Will not be able to process single-atom topologies.')]

        # Version checking
        self.gmxversion = subprocess.Popen(self.grompp + ["-version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.gmxversion = self.gmxversion.communicate()[0].decode('utf-8')
        self.gmxversion = float(re.search('\D(\d+\.\d+)(\.\d+)?\D', self.gmxversion ).groups()[0])
        if self.gmxversion < 4.5:
            sys.exit('Error: GMX version must be 4.5 or higher.')
        elif self.gmxversion >= 2020:
            sys.exit('Error: GMX versions 2020 and greater are not supported.')

        if not self.values.ff:
            try:
                self.gmxlibdir = Path(os.environ['GMXLIB'])
            except KeyError:
                self.gmxlibdir = self.grompp[0].parents[1]/'share/gromacs/top'
            if not self.gmxlibdir.exists():
                sys.exit('Error: forcefield not specified and I can\'t find the standard GMX forcefield tree.')

            ffs = [top.name for top in self.gmxlibdir.iterdir() if top.suffix == '.ff']
            if not ffs:
                sys.exit('Error: forcefield not specified and I can\'t find any in $GMXLIB.')
            for i, ff in enumerate(ffs):
                print('%3i: %s' % (i, ff[:-3]))
            chosenff = int(input('Forcefield to use (number only): '))
            chosenff = ffs[chosenff]
            self.values.ff = self.gmxlibdir/chosenff
    
    def getmol(self):
        molecules = {}
        curr_molecule = []
        for line in open(self.itp):
            line = clean(line)
            if not line:
                continue
            if len(line) == 1 and 'moleculetype' in line:
                # starting a new molecule. let's process the current one, if it exists
                if curr_molecule:
                    itp = ITP(curr_molecule)
                    molecules[itp.name] = itp
                    curr_molecule = []
            curr_molecule.append(line)
        # The final molecule
        itp = ITP(curr_molecule)
        molecules[itp.name] = itp

        if not molecules:
            sys.exit(f'Error: Wasn\'t able to identify any molecules in {self.itp}.')

        if len(molecules) == 1:
            self.molecule = molecules[list(molecules.keys())[0]]
        else:
            if self.values.mol is not None:
                self.mol = molecules[self.values.mol]
            else:
                molnames = '\n'.join(molecules.keys())
                print(molnames)
                self.molecule = molecules[input('Which molecule to convert? ')]
        self.molname = self.molecule.name

    def createitp(self):
        if not self.values.keep_constraints:
            newbonds = copy.deepcopy(self.molecule.topology['constraints'])
            del self.molecule.topology['constraints']
            for bond in newbonds:
                bond.append(1000000)
            bonds = self.molecule.topology.setdefault('bonds', [])
            bonds.extend(newbonds)

        self.run_itp = self.basename.with_suffix('.itp')

        with open(self.run_itp, 'w') as itp_out:
            for directive, lines in self.molecule.topology.items():
                print(f'[ {directive} ]', file=itp_out)
                for line in lines:
                    text = '  '.join([str(v) for v in line])
                    print(text, file=itp_out)

    def creategro(self):
        self.atoms = len(self.molecule.topology['atoms'])
        
        self.init_gro = self.basename.with_suffix('.gro')
        with open(self.init_gro, 'w') as gro:
            print(self.molname, file=gro)
            print(self.atoms, file=gro)
            for i in range(self.atoms):
                print('%5d%-5s%5s%5d%8.3f%8.3f%8.3f'
                      % (i+1,"DUM","DUM",i+1,
                         0.05*i,
                         self.values.fuzzy*random.random()+0.025*self.atoms-0.5*self.values.fuzzy,
                         self.values.fuzzy*random.random()+0.025*self.atoms-0.5*self.values.fuzzy),
                      file=gro)
            box_len = max(0.2*self.atoms, 2.0)
            print('%f %f %f' % (box_len, box_len, box_len), file=gro)
    
    def createtop(self):
        template = '''#include "%s"
#include "%s"
[ system ]
%s
[ molecules ]
%s   1
'''
        self.top = self.basename.with_suffix('.top')
        top = open(self.top, 'w')
        top.write(template % (self.values.ff.resolve(), self.run_itp.resolve(), self.molname, self.molname))
        top.close()

    def createmdp(self):
        template = mdp_template % (self.values.intg, str(self.values.nsteps), str(self.values.eps))
        if self.gmx: #GMX >= 5.0
            template += mdp_template_gmx5
        self.mdp = self.basename.with_suffix('.mdp')
        mdp = open (self.mdp, 'w')
        mdp.write(template)
        if self.values.xmdp:
            mdp.write(''.join(open(self.values.xmdp).readlines()))
        mdp.close()

    def minimize(self):
        self.tpr = self.basename.with_suffix('.tpr')
        self.mdout = self.basename.with_suffix('.out.mdp')
        self.pplog = self.basename.with_suffix('.pp.log')
        self.deffnm = self.basename.with_name(self.basename.name + '_md')
        self.mdgro = self.deffnm.with_suffix('.gro')
        self.mdtrr = self.deffnm.with_suffix('.trr')
        self.mdlog = self.deffnm.with_suffix('.log')
        self.mdedr = self.deffnm.with_suffix('.edr')

        ppargs = self.grompp + f'-f {self.mdp} -p {self.top} -c {self.init_gro} -maxwarn 3 -po {self.mdout} -o {self.tpr}'.split()
        mdargs = self.mdrun + f'-s {self.tpr} -nt 1 -deffnm {self.deffnm} -cpt 0'.split()

        with open(self.pplog, 'w') as pplog:
            md_env = os.environ.copy()
            md_env['GMX_SUPPRESS_DUMP'] = '1'
            md_env['GMX_MAXBACKUP'] = '-1'
            pp = subprocess.call(ppargs, stdout=pplog, stderr=pplog, env=md_env)
        if pp:
            sys.exit(f'grompp error: check {self.pplog}')
        with open(self.mdlog, 'w') as mdlog:
            if self.atoms > 1:
                md = subprocess.call(mdargs, stdout=mdlog, stderr=mdlog, env=md_env)
                if md:
                    sys.exit(f'mdrun error: check {self.mdlog}')
            else:
                # The minimization gives all 'nan' in this case, let's just replace
                #  with the .tpr structure
                if not self.editconf[0]:
                    sys.exit(f'Cannot handle single-atom topologies because no \'editconf\' executable was found.')
                ecargs = self.editconf + f'-f {self.tpr} -o {self.deffnm}.gro'.split()
                ec = subprocess.call(ecargs, stdout=mdlog, stderr=mdlog, env=md_env)
                if ec:
                    sys.exit(f'Error creating single-atom structure: check {self.mdlog}')

    def cleanup(self):
        self.outtpr = self.gro.with_suffix('.tpr')
        self.outtrr = self.gro.with_suffix('.trr')
        if self.trjconv[0]:
            tcargs = self.trjconv + f'-f {self.mdgro} -o {self.gro} -s {self.tpr} -center -pbc mol'.split()
            trjc = subprocess.Popen(tcargs, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            trjc.communicate(b'0\n0\n')
            tc = trjc.returncode
            if tc:
                sys.stderr.write('Warning: could not center the final structure.\n')
                shutil.copy(self.mdgro, self.gro)
            if self.values.traj:
                tcargs = self.trjconv + f'-f {self.mdtrr} -o {self.outrr} -s {self.tpr} -center -pbc mol'.split()
                trjc = subprocess.Popen(tcargs, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                trjc.communicate(b'0\n0\n')
                tc = trjc.returncode
                if tc:
                    sys.stderr.write("Warning: could not center the trajectory.\n")
        else:
            shutil.copy(self.mdgro, self.gro)
        if self.values.traj:
            shutil.copy(self.tpr, self.outtpr)
            
        if not self.values.keep:
            for fname in [self.run_itp, self.init_gro, self.top, self.mdp,
                          self.mdout, self.pplog, self.tpr, self.mdgro,
                          self.mdtrr, self.mdlog, self.mdedr]:
                try:
                    fname.unlink()
                except FileNotFoundError:
                    pass

    def justdoit(self):
        self.getmol()
        self.createitp()
        self.creategro()
        self.createmdp()
        self.createtop()
        self.minimize()
        self.cleanup()

##########
# Main
##########
def main(argv=None):
    if not argv:
        argv = sys.argv[1:]
    #
    # Option parsing and 'documentation'.
    desc ="""
******************************************
molmaker.py attempts to generate a structure from a GROMACS .itp and a forcefield, by minimizing initially semi-random coordinates while slowly fading in nonbonded interactions.

It is hopeless to use this for molecules with many degrees of freedom, especially if they're somewhat coiled or branched (try increasing the fuzziness in these cases). It does work out well for creating linear stretches of protein.

The script will always try to do the RightThing(TM), but it'll need help at times, namely a) if you're using a custom forcefield not present under the $GMXLIB path, and b) if your molecule requires exotic .mdp directives to minimize properly.

The used VdW and Coulombic interactions are cutoff at 1nm. If potential shapes/cuttoffs really are important, the .mdp template is readily editable at the top of the script.

If you're using charges and a forcefield with implicit screening (such as Martini) you'll want to set -eps to something else than the default.

By default, constraints in the topology are converted to harmonic bonds with 1e6 kJ/mol/nm^2 force constant.

The minimal .mdp used for minimization is the following:

%s
%s
******************************************

where %%intg%%, %%nsteps%% and %%eps%% are the values specified with -intg, -nsteps and -eps.

Finally, the script is not that clever that it can read #included files from the .top/.itp you provide. Make sure that at least the [ moleculetype ] and [ atoms ] directives are in the file you supply. If they appear more than once you'll get to choose the molecule you want, or you can preempt that with -mol.

Version 2.1.0-30-06-2020 by Manuel Melo (m.n.melo@itqb.unl.pt)""" % (mdp_template%("%intg%","%nsteps%","%eps%"), mdp_template_gmx5)

    builder = MolMaker(description=desc)
    builder.add_argument('-i', dest='itp', type=Path, help='The target molecule topolgy file. Default is to ask for it if none given.')
    builder.add_argument('-x', dest='xmdp', type=Path, help='File with extra directives to include in the minimization .mdp.')
    builder.add_argument('-ff', dest='ff', type=Path, help='The forcefield itp. Default is to ask for it if none given. You\'ll most likely want to use this flag if using the MARTINI CG forcefield.')
    builder.add_argument('-o', dest='gro', type=Path, help='The output .gro file. Default is to infer from the topology file.')
    builder.add_argument('-mol', default=None, help='Which molecule to convert, in case there are multiple in the topology file.')
    builder.add_argument('-temp-prefix', dest='tmpprefix', default='.molmk_', help='The prefix of the temp files.')
    builder.add_argument('-fuzzy', type=float, dest='fuzzy', default=2, help='Fuzziness of the randomly placed initial coordinates.')
    builder.add_argument('-keep', action='store_true', dest='keep', help='Whether to keep all temporary output files (note that these will be hidden by default anyway).')
    builder.add_argument('-traj', action='store_true', dest='traj', help='Whether to save the minimization trajectory and a .tpr file.')
    builder.add_argument('-intg', dest="intg", default='cg', help='The energy-minimization integrator to use.')
    builder.add_argument('-nsteps', type=int, dest='nsteps', default=5000, help='Maximum number of steps in minimization (might be useful to limit this for CG minimizations).')
    builder.add_argument('-eps', type=float, dest='eps', default=1.0, help='The relative dielectric constant.')
    builder.add_argument('-constr', dest='keep_constraints', action='store_true', help='Whether to keep constraints as such, instead of converting to harmonic bonds.')
    
    builder.checkopts(args=argv)
    builder.justdoit()

if __name__ == '__main__':
    sys.exit(main())


