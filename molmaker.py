#!/usr/bin/env python3

###  Copyright 2016-2021 Manuel N. Melo ###
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
import warnings

# Constants
############
COMMENT = 1
PREPROCESSING = 2
comment_pp_chars = {';': COMMENT,
                    '#': PREPROCESSING}

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

def is_comment_or_preprocessing(val):
    '''Returns 1 if a comment, 2 if preprocessing, 0 if neither.'''
    if not val:
        return 0
    try:
        while True:
            try:
                return comment_pp_chars[val]
            except KeyError:
                pass
            oldval, val = val, val[0]
            if val is oldval:
                return 0
    except TypeError:
        return 0

def convert_val(val):
    try:
        try:
            return int(val)
        except ValueError:
            return float(val)
    except ValueError:
        return str(val)

def clean(line, keep_preprocessing=True, keep_comments=True):
    line = line.strip()
    comment_or_pp = is_comment_or_preprocessing(line)

    if ((comment_or_pp is COMMENT and not keep_comments) or
        (comment_or_pp is PREPROCESSING and not keep_preprocessing)):
        return None
    if not keep_comments:
        line = line.split(';')[0]
        line = line.strip()
    if not line:
        return None

    # Let's deal with preprocessing/comments here and pass them (mostly) untouched
    if comment_or_pp:
        return [line]
    # remove possible spaces around directive names
    if line.startswith('['):
        line = ''.join(line.split())
    # comments midway through the line get a bit more special treatment
    try:
        line, comment = line.split(';', maxsplit=1)
        comment = ';' + comment
        vals = line.split() + [comment]
    except ValueError:
        vals = line.split()
    return [convert_val(val) for val in vals]


# Classes
##########

class ProperFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """A hackish class to get proper help format from argparse.

    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class ITP:
    at_ids = {'moleculetype': 0,
              'atoms': 1,
              'bonds': 2,
              'pairs': 2,
              'pairs_nb': 2,
              'angles': 3,
              'dihedrals': 4,
              'exclusions': None,
              'constraints': 2,
              'settles': 1,
              'virtual_sites2': 3,
              'virtual_sites3': 4,
              'virtual_sites4': 5,
              'virtual_sitesn': 0,
              'position_restraints': 1,
              'distance_restraints': 2,
              'dihedral_restraints': 4,
              'orientation_restraints': 2,
              'angle_restraints': 4,
              'angle_restraints_z': 2,}

    def __init__(self, val_lists, order_sections=None):
        self.order_sections = order_sections
        try:
            # Constructor from existing topology dictionary
            val_lists['moleculetype']
            self.topology = val_lists
        except TypeError:
            self.topology = {}
            self._parse(val_lists)

    @property
    def name(self):
        return self.topology['moleculetype'][0][0]

    @name.setter
    def name(self, newname):
        self.topology['moleculetype'][0][0] = str(newname)

    def _parse(self, val_lists):
        curr_section = pre_section = self.topology['pre_section'] = []
        for vals in val_lists:
            first_str = str(vals[0])
            if first_str.startswith('['): # new section
                section = first_str[1:-1]
                if curr_section is pre_section:
                    if section != 'moleculetype':
                        raise TypeError('Topology doesn\'t start with directive "moleculetype"')
                elif section == 'moleculetype':
                    raise TypeError('Directive "moleculetype" is not first in topology.')
                curr_section = self.topology.setdefault(section, [])
            else:
                if is_comment_or_preprocessing(vals):
                    if self.order_sections:
                        raise ValueError('Section ordering was requested, but '
                                         'incompatible topology lines are '
                                         'present (either comment-only lines or '
                                         'preprocessing directives)')
                    self.order_sections = False
                elif curr_section is pre_section:
                    raise ValueError('Non-comment/preprocessing line before '
                                     'directive "moleculetype".')
                curr_section.append(vals)

        if not pre_section:
            del self.topology['pre_section']

        # Make the order predictable
        if self.order_sections:
            for section, vals in self.topology.items():
                if section not in ('pre_section', 'atoms'):
                    vals[:] = sorted(vals)

    def write(self, outfile):
        try:
            OUT = open(outfile, 'w')
        except TypeError:
            OUT = outfile
        with OUT:
            for section, lines in self.topology.items():
                if section != 'pre_section':
                    print(f'[ {section} ]', file=OUT)
                for line in lines:
                    if is_comment_or_preprocessing(line):
                        print(line[0], file=OUT)
                        continue
                    at_strs = ' '.join([f'{v:<3}' for v in line[:self.at_ids[section]]])
                    if self.at_ids[section] is None:  # exclusions' special case
                        print(f'  {at_strs}', file=OUT)
                        continue
                    xtra_strs = ' '.join([f'{v:<5}' for v in line[self.at_ids[section]:]])
                    print(f'  {at_strs} {xtra_strs}', file=OUT)
                print(file=OUT)
            print(file=OUT)

    def renum(self, delta):
        '''Adds delta to each atom index in the topology and returns a new ITP
        with the renumbering.
           Does not change the ITP object in-place!
        '''
        new_top = copy.deepcopy(self.topology)
        for directive, lines in new_top.items():
            if directive == 'virtual_sitesn':
                warnings.warn('This topology uses virtual_sitesn, which cannot '
                              'be automatically renumbered. Please correct '
                              'those afterwards.')
            n_ats = self.at_ids[directive]
            if n_ats == 0: # moleculetype case
                continue
            for line in lines:
                if is_comment_or_preprocessing(line):
                    continue
                new_vals = [v + delta for v in line[:n_ats]]
                line[:n_ats] = new_vals
        return ITP(new_top, order_sections=False)

    @classmethod
    def merge(cls, *itps, newname=None):
        new_itp = ITP(itps[0].topology)
        for itp in itps[1:]:
            for directive, lines in itp.topology.items():
                if directive == 'moleculetype':
                    continue
                try:
                    new_itp.topology[directive].extend(lines)
                except KeyError:
                    new_itp.topology[directive] = copy.deepcopy(lines)
        if newname is not None:
            new_itp.topology['moleculetype'][0][0] = newname
        return new_itp


class ITPFile:
    def __init__(self, itp_path, keep_preprocessing=True, keep_comments=False):
        self.itp_path = Path(itp_path)
        order_sections = not (keep_preprocessing or keep_comments)
        self.lines = []
        for line in open(itp_path):
            line = clean(line,
                         keep_preprocessing=keep_preprocessing,
                         keep_comments=keep_comments)
            if line:
                self.lines.append(line)

        linetypes = []
        mol_starts = []
        for lineno, line in enumerate(self.lines):
            val = str(line[0])
            if val == '[moleculetype]':
                linetypes.append('molstart')
                mol_starts.append(lineno)
            elif val.startswith('['):
                linetypes.append('section')
            elif val.startswith(';'):
                linetypes.append('comment')
            elif val.startswith('#end'):
                linetypes.append('directive_pre')
            elif val.startswith('#'):
                linetypes.append('directive')
            else:
                linetypes.append('val')

        for mol_idx, start_lineno in enumerate(mol_starts):
            if not start_lineno:  # start_lineno == 0 breaks the reverse slice below
                continue
            for count, prev_type in enumerate(linetypes[start_lineno-1::-1]):
                if prev_type in ('directive_pre', 'val'):
                    mol_starts[mol_idx] -= count
                    break
            else:
                mol_starts[mol_idx] = 0

        mol_starts.append(None)  # So that the zip trick works for the last slice
        self.molecules = {}
        for a, b in zip(mol_starts, mol_starts[1:]):
            itp = ITP(self.lines[slice(a, b)], order_sections=order_sections)
            self.molecules[itp.name] = itp

    def write(self, outfile):
        with open(outfile, 'w') as OUT:
            for molecule in self.molecules.values():
                molecule.write(outfile)


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
            self.gro = Path(self.itp.with_suffix('.gro').name)
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
            self.values.ff = self.gmxlibdir/chosenff/'forcefield.itp'
    
    def getmol(self):
        mol_file = ITPFile(self.itp, keep_preprocessing=False, keep_comments=False)
        molecules = mol_file.molecules

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
        if ('constraints' in self.molecule.topology
                and not self.values.keep_constraints):
            newbonds = copy.deepcopy(self.molecule.topology['constraints'])
            del self.molecule.topology['constraints']
            for bond in newbonds:
                try:
                    bond.append(1000000)
                except AttributeError:
                    print(newbonds)
                    raise
            bonds = self.molecule.topology.setdefault('bonds', [])
            bonds.extend(newbonds)

        self.run_itp = self.basename.with_suffix('.itp')
        self.molecule.write(self.run_itp)

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
                tcargs = self.trjconv + f'-f {self.mdtrr} -o {self.outtrr} -s {self.tpr} -center -pbc mol'.split()
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


