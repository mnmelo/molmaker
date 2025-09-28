#!/usr/bin/env python3
#-*- coding: utf-8 -*-

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
#   but leave the {} fields alone!        #
#   (you can set them with option flags)  #
###########################################

mdp_template = '''\
define            = -DFLEXIBLE
integrator        = {integrator}
emtol             = 0.1
nsteps            = {nsteps}
nstxout           = {nstxout}
nstcgsteep        = 50

rlist             = {rlist}
rcoulomb          = {rcoulomb}
rvdw              = {rvdw}
epsilon_r         = {eps}

lincs_order       = 8
lincs_iter        = 2
continuation      = yes

free-energy       = yes
init-lambda       = 0
delta-lambda      = {delta_lambda}
couple-moltype    = system
couple-lambda0    = none
couple-lambda1    = vdw-q
couple-intramol   = yes
sc-function       = gapsys
'''

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
        proc = subprocess.Popen(['which', name],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
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

def sanitize_filename(fname):
    return re.sub('[.\/]', '_', fname)

def is_comment_or_preprocessing(val):
    '''Returns 1 if a comment, 2 if preprocessing, 0 if neither.'''
    if not val:
        return 0

    # Unpack the vals to get the very first character
    while True:
        try:
            oldval, val = val, val[0]
        except (TypeError, IndexError):
            break
        if val is oldval:
            break

    return comment_pp_chars.get(val, 0)

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

    # We deal with preprocessing/comments here and pass them (mostly) untouched
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

def get_format_lens(lines, nats=None, minlen=3):
    formats = []
    for line in lines:
        if is_comment_or_preprocessing(line):
            continue
        for i, val in enumerate(line):
            if is_comment_or_preprocessing(val):
                break
            val_len = len(str(val))
            try:
                formats[i] = max(formats[i], val_len)
            except IndexError:
                formats.append(max(val_len, minlen))
    if nats is None:
        # Cool/hackish way to do an infinite iterator without
        # having to import itertools
        return iter(lambda:max(formats), 0)
    if nats > 1:
        formats[:nats] = [max(formats[:nats])] * nats
    return formats


# Classes
##########

class ProperFormatter(argparse.RawTextHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    '''A hackish class to get proper help format from argparse.'''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class ITP:
    at_ids = {'pre_section': 0,
              'post_section': 0,
              'moleculetype': 0,
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
              'virtual_sitesn': None,
              'position_restraints': 1,
              'distance_restraints': 2,
              'dihedral_restraints': 4,
              'orientation_restraints': 2,
              'angle_restraints': 4,
              'angle_restraints_z': 2,}

    def __init__(self, val_lists=None, moleculetype='default', nrexc=3, order_sections=None, parent_collection=None):
        '''Creates an ITP representation as a dict of directive sections.

        Each section is a list of lists of line values.
        The first argument should either be a list of cleaned lines, or a
        dict of sections to make a copy from.'''

        self.order_sections = order_sections
        if val_lists is None:
            self.topology = {'moleculetype': [[moleculetype, nrexc]]}
        else:
            try:
                # Constructor from existing topology dictionary
                val_lists['moleculetype']
                self.topology = copy.deepcopy(val_lists)
            except TypeError:
                self.topology = {}
                self._parse(val_lists)
        # We register with the parent collection, if given
        if parent_collection is not None:
            if ('pre_section' not in self.topology
                    and parent_collection.molecules):
                post_section = parent_collection[-1].topology.pop('post_section',
                                                                  None)
                if post_section is not None:
                    self.topology['pre_section'] = post_section
            self.register(parent_collection)

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
                        raise TypeError('Topology doesn\'t start with '
                                        'directive "moleculetype"')
                elif section == 'moleculetype':
                    raise TypeError('Directive "moleculetype" is not first '
                                    'in topology.')
                curr_section = self.topology.setdefault(section, [])
            else:
                if is_comment_or_preprocessing(vals):
                    if self.order_sections:
                        raise ValueError('Section ordering was requested, but '
                                         'incompatible topology lines are '
                                         'present (either comment-only lines '
                                         'or preprocessing directives)')
                elif curr_section is pre_section:
                    raise ValueError(f'Non-comment/preprocessing line {vals} '
                                     'before directive "moleculetype".')
                curr_section.append(vals)
    
        post_section = self.topology['post_section'] = []
        for line in list(self.topology.values())[-1][::-1]:
            if is_comment_or_preprocessing(line):
                post_section.append(line)
            else:
                break
        if post_section:
            post_section[:] = post_section[::-1]
        else:
            del self.topology['post_section']

        if not pre_section:
            del self.topology['pre_section']

        # Make the order predictable
        if self.order_sections:
            for section, vals in self.topology.items():
                if section not in ('pre_section', 'post_section', 'atoms'):
                    vals[:] = sorted(vals)

    def write(self, outfile):
        own_out = True
        try:
            OUT = open(outfile, 'w')
        except TypeError:
            OUT = outfile
            own_out = False # To know whether to close it ourselves

        for section, lines in self.topology.items():
            if lines:
                format_lens = get_format_lens(lines, self.at_ids[section])
                if section not in ('pre_section', 'post_section'):
                    print(f'[ {section} ]', file=OUT)

            for line in lines:
                if is_comment_or_preprocessing(line):
                    print(line[0], file=OUT)
                    continue
                line_str = ''
                for val, fmt_len in zip(line, format_lens):
                    line_str += f'  {val:<{fmt_len}}'
                print(line_str, file=OUT)
            if is_comment_or_preprocessing(line) != COMMENT: 
                print(file=OUT)
        print(file=OUT)
        if own_out:
            OUT.close()

    def reindex(self, new_index, just_renumber=False):
        '''Reorders an itp and/or renumbers each atom interaction index.

        New index list must be 1-based.
        Can possibly also remove atoms, if not all indices are passed. 
        Does not change the ITP object in-place!
        '''
        if just_renumber:
            convert_dict = {old_id: new_id for old_id, new_id in enumerate(new_index, 1)}
        else:
            convert_dict = {old_id: new_id for new_id, old_id in enumerate(new_index, 1)}

        new_top = copy.deepcopy(self.topology)
        for directive, lines in new_top.items():
            n_ats = self.at_ids[directive]
            if n_ats == 0: # moleculetype case
                continue
            new_directive = []
            for line in lines:
                try:
                    new_line = line.copy() 
                    if is_comment_or_preprocessing(new_line):
                        pass
                    elif directive == 'virtual_sitesn':
                        change_ndx = list(range(len(new_line)))
                        if new_line[1] == 3:
                            del change_ndx[1::2]
                        else:
                            del change_ndx[1]
                        for c_ndx in change_ndx:
                            new_line[c_ndx] = convert_dict[new_line[c_ndx]]
                    else:
                        new_vals = [convert_dict[v] for v in new_line[:n_ats]]
                        new_line[:n_ats] = new_vals
                except KeyError:
                    continue
                new_directive.append(new_line)
            lines[:] = new_directive

        # Reorder chgrps
        #cgrps = [at_line[5] for at_line in new_top['atoms']
        #         if not is_comment_or_preprocessing(at_line)]
        #new_idxs = list(range(1, len(cgrps) + 1))
        #new_cgrps = []
        #old_cgrp = None
        #for cgrp in cgrps:
        #    if cgrp == old_cgrp:
        #        new_cgrps.append(new_cgrps[-1])
        #    else:
        #        new_cgrps.append(new_idxs.pop(0))
        #    old_cgrp = cgrp

        for at_line in new_top['atoms']:
            if not is_comment_or_preprocessing(at_line):
                #at_line[5] = new_cgrps.pop(0)
                at_line[5] = at_line[0]

        return ITP(new_top, order_sections=False)

    def renum(self, at_delta=0, res_delta=0):
        '''Adds delta to each atom/residue index and returns a renumbered ITP

        Does not change the ITP object in-place!
        '''
        new_index = [at_line[0] + at_delta for at_line in self.topology['atoms']
                     if not is_comment_or_preprocessing(at_line)]
        new_top = self.reindex(new_index, just_renumber=True)

        if res_delta:
            for at_line in new_top['atoms']:
                if not is_comment_or_preprocessing(at_line):
                    at_line[2] += res_delta

        return new_top

    def register(self, parent_collection):
        self.parent_collection = parent_collection
        self.parent_collection[self.name] = self

    def copy(self):
        return ITP(self.topology)

    @classmethod
    def merge(cls, *itps, newname=None):
        new_itp = ITP(itps[0].topology)
        for i, itp in enumerate(itps[1:]):
            new_itp.append(itp)
        if newname is not None:
            new_itp.topology['moleculetype'][0][0] = newname
        return new_itp

    def append(self, other):
        other = other.renum(self.n_atoms, self.n_res)
        for directive, lines in other.topology.items():
            if directive == 'moleculetype':
                continue
            try:
                self.topology[directive].extend(copy.deepcopy(lines))
            except KeyError:
                self.topology[directive] = copy.deepcopy(lines)

    def __getitem__(self, key):
        return self.topology[key]

    @property
    def n_atoms(self):
        if not 'atoms' in self.topology:
            return 0
        return len([at_line for at_line in self.topology['atoms']
                    if not is_comment_or_preprocessing(at_line)])

    @property
    def n_res(self):
        if not 'atoms' in self.topology:
            return 0
        return len(set([at_line[2] for at_line in self.topology['atoms']
                        if not is_comment_or_preprocessing(at_line)]))


class ITPFile:
    def __init__(self, itp_path, keep_preprocessing=True, keep_comments=False,
                 ITPclass=ITP):
        '''Parses an .itp containing one or more molecules.

        Optionally keeps preprocessing stuf and/or comments.
        Constructs a dictionary of ITP instances; the actual used ITP class can
        be overriden as an argument, to allow for customization/subclassing.

        The ITPFile instance is always passed as the first argument to ITPclass
        initialization.
        '''
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
            # start_lineno == 0 breaks the reverse slice below
            if not start_lineno:
                continue
            for count, prev_type in enumerate(linetypes[start_lineno-1::-1]):
                if prev_type in ('val'):
                    mol_starts[mol_idx] -= count
                    break
            else:
                mol_starts[mol_idx] = 0

        # So that the zip trick works for the last slice
        mol_starts.append(None)
        self.molecules = {}
        for a, b in zip(mol_starts, mol_starts[1:]):
            ITPclass(self.lines[slice(a, b)],
                     order_sections=order_sections,
                     parent_collection=self)
            # The molecule self-registers, so no need to do it here

    def write(self, outfile):
        with open(outfile, 'w') as OUT:
            for molecule in self.molecules.values():
                molecule.write(OUT)

    @property
    def pre_section(self):
        return self[0].topology.get('pre_section')

    @property
    def post_section(self):
        return self[-1].topology.get('post_section')

    def __getitem__(self, key):
        try:
            return self.molecules[key]
        except KeyError as err:
            # Let's assume it's int indexing
            try:
                return list(self.molecules.values())[key]
            except (TypeError, KeyError):
                raise(err)

    def __setitem__(self, key, val):
        self.molecules[key] = val


class MolMaker(argparse.ArgumentParser):

    top_template = '''\
#include "{ff_itp}"
#include "{mol_itp}"
[ system ]
{molname}
[ molecules ]
{molname}   1
'''

    def __init__(self, formatter_class=ProperFormatter, *args, **kwargs):
        argparse.ArgumentParser.__init__(self,
                                         *args,
                                         formatter_class=formatter_class,
                                         **kwargs)

    def checkopts1(self, args=sys.argv):
        self.values = self.parse_args(args)
        #
        # Options checking and assignment of defaults not done by argparse.
        self.itp = self.values.itp
        if not self.itp:
            self.itp = Path(input('Topology file: '))

        if not self.itp.exists():
            sys.exit('Error: can\'t find the specified topology file.')

        if not self.values.tmpprefix:
            sys.exit('Error: -temp-prefix cannot be an empty string.')

        self.values.fuzzy = max(self.values.fuzzy, 0)

        # Are we running gromacs >=5?
        self.gmx = find_exec(('gmx', 'gmx_d', 'gmx_mpi', 'gmx_mpi_d'),
                             error=None)

        self.grompp = [self.gmx, "grompp"]
        self.mdrun = [self.gmx, "mdrun"]
        self.trjconv = [self.gmx, "trjconv"]
        self.editconf = [self.gmx, "editconf"]

        ## Version checking
        #self.gmxversion = subprocess.Popen(self.grompp + ["-version"],
        #                                   stdout=subprocess.PIPE,
        #                                   stderr=subprocess.PIPE)
        #self.gmxversion = self.gmxversion.communicate()[0].decode('utf-8')
        #self.gmxversion = float(re.search('\D(\d+\.\d+)(\.\d+)?\D',
        #                                  self.gmxversion ).groups()[0])
        #if self.gmxversion < 4.5:
        #    sys.exit('Error: GMX version must be 4.5 or higher.')

        if not self.values.ff:
            try:
                self.gmxlibdir = Path(os.environ['GMXLIB'])
            except KeyError:
                self.gmxlibdir = self.grompp[0].parents[1]/'share/gromacs/top'
            if not self.gmxlibdir.exists():
                sys.exit('Error: forcefield not specified and I can\'t find '
                         'the standard GMX forcefield tree.')

            ffs = [top.name for top in self.gmxlibdir.iterdir()
                   if top.suffix == '.ff']
            if not ffs:
                sys.exit('Error: forcefield not specified and I can\'t find '
                         'any in $GMXLIB.')
            for i, ff in enumerate(ffs):
                print('%3i: %s' % (i, ff[:-3]))
            chosenff = int(input('Forcefield to use (number only): '))
            chosenff = ffs[chosenff]
            self.values.ff = self.gmxlibdir/chosenff/'forcefield.itp'

    def checkopts2(self):
        self.gro = self.values.gro
        if not self.gro:
            self.gro = Path(self.molname).with_suffix('.gro')
        self.name = sanitize_filename(self.gro.stem)
        self.basename = self.gro.with_name(self.values.tmpprefix + self.name)

    def getmol(self):
        mol_file = ITPFile(self.itp, keep_preprocessing=True,
                           keep_comments=False)
        molecules = mol_file.molecules

        if not molecules:
            sys.exit(f'Error: Wasn\'t able to identify any molecules '
                      'in {self.itp}.')

        if len(molecules) == 1:
            self.molecule = molecules[list(molecules.keys())[0]]
        else:
            if self.values.mol is not None:
                self.molecule = molecules[self.values.mol]
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
                bond.append(1000000)

            if 'bonds' not in self.molecule.topology:
                # must add bonds, but ensure exclusions come last
                self.molecule.topology['bonds'] = []
                if 'exclusions' in self.molecule.topology:
                    exc = self.molecule.topology.pop('exclusions')
                    self.molecule.topology['exclusions'] = exc

            self.molecule.topology['bonds'].extend(newbonds)

        self.run_itp = self.basename.with_suffix('.itp')
        self.molecule.write(self.run_itp)

    def creategro(self):
        self.atoms = len(self.molecule.topology['atoms'])

        self.init_gro = self.basename.with_suffix('.gro')
        with open(self.init_gro, 'w') as gro:
            print(self.molname, file=gro)
            print(self.atoms, file=gro)
            xs = []
            ys = []
            zs = []
            for i in range(self.atoms):
                x = 0.05*i
                y = self.values.fuzzy * (random.random() - .5) + 0.025 * self.atoms 
                z = self.values.fuzzy * (random.random() - .5) + 0.025 * self.atoms 
                xs.append(x)
                ys.append(y)
                zs.append(z)
                print('%5d%-5s%5s%5d%8.3f%8.3f%8.3f'
                      % (i+1,"DUM","DUM",i+1,x,y,z),file=gro)
            box_len = max([max(dim) - min(dim) for dim in (xs, ys, zs)]) * 2
            box_len = max(box_len, 5.0)
            print('%f %f %f' % (box_len, box_len, box_len), file=gro)
        return box_len

    def createtop(self):
        self.top = self.basename.with_suffix('.top')
        with open(self.top, 'w') as top:
            top.write(self.top_template.format(ff_itp=self.values.ff.resolve(),
                                               mol_itp=self.run_itp.resolve(),
                                               molname=self.molname))

    def createmdp(self, cutoff=1.2):
        template = mdp_template.format(integrator=self.values.intg,
                                       nsteps=self.values.nsteps,
                                       nstxout=int(self.values.traj or self.values.keep),
                                       eps=self.values.eps,
                                       rlist=cutoff,
                                       rcoulomb=cutoff,
                                       rvdw=cutoff,
                                       delta_lambda=1/self.values.nsteps)
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

        ppargs = self.grompp + (f'-f {self.mdp} -p {self.top} '
                                f'-c {self.init_gro} -maxwarn 3 '
                                f'-po {self.mdout} -o {self.tpr}'.split())
        mdargs = self.mdrun + (f'-s {self.tpr} -nt 1 '
                               f'-deffnm {self.deffnm} -cpt 0'.split())

        with open(self.pplog, 'w') as pplog:
            md_env = os.environ.copy()
            md_env['GMX_SUPPRESS_DUMP'] = '1'
            md_env['GMX_MAXBACKUP'] = '-1'
            pp = subprocess.call(ppargs, stdout=pplog, stderr=pplog,
                                 env=md_env)
        if pp:
            sys.exit(f'grompp error: check {self.pplog}')
        with open(self.mdlog, 'w') as mdlog:
            if self.atoms > 1:
                md = subprocess.call(mdargs, stdout=mdlog, stderr=mdlog,
                                     env=md_env)
                if md:
                    sys.exit(f'mdrun error: check {self.mdlog}')
            else:
                # The minimization gives all 'nan' in this case, let's just
                #  replace with the .tpr structure
                if not self.editconf[0]:
                    sys.exit('Cannot handle single-atom topologies because '
                             'no \'editconf\' executable was found.')
                ecargs = self.editconf + (f'-f {self.tpr} '
                                          f'-o {self.deffnm}.gro'.split())
                ec = subprocess.call(ecargs, stdout=mdlog, stderr=mdlog,
                                     env=md_env)
                if ec:
                    sys.exit('Error creating single-atom structure: '
                             f'check {self.mdlog}')

    def cleanup(self):
        self.outtpr = self.gro.with_suffix('.tpr')
        self.outtrr = self.gro.with_suffix('.trr')
        if self.trjconv[0]:
            tcargs = self.trjconv + (f'-f {self.mdgro} -o {self.gro} '
                                     f'-s {self.tpr} -center -pbc mol'.split())
            trjc = subprocess.Popen(tcargs,
                                    stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            trjc.communicate(b'0\n0\n')
            tc = trjc.returncode
            if tc:
                print('Warning: could not center '
                      'the final structure.', file=sys.stderr)
                shutil.copy(self.mdgro, self.gro)
            if self.values.traj:
                tcargs = self.trjconv + (f'-f {self.mdtrr} -o {self.outtrr} '
                                         f'-s {self.tpr} -center '
                                         f'-pbc mol'.split())
                trjc = subprocess.Popen(tcargs,
                                        stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
                trjc.communicate(b'0\n0\n')
                tc = trjc.returncode
                if tc:
                    print('Warning: could not center '
                          'the trajectory.', file=sys.stderr)
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

    def justdoit(self, args):
        self.checkopts1(args)
        self.getmol()
        self.checkopts2()
        self.createitp()
        box_len = self.creategro()
        self.createmdp(cutoff=box_len/2.15)
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
    desc ='''\
******************************************
molmaker.py attempts to generate a structure from a GROMACS .itp and a
forcefield, by minimizing initially semi-random coordinates while slowly fading
in nonbonded interactions. It uses the free-energy code for this soft-core
potential fading, which means that your molecule should be compatible with
such potential de-coupling.

It is hopeless to use this for molecules with many degrees of freedom,
especially if they\'re somewhat coiled or branched (try increasing the
fuzziness in these cases). It does work out well for creating linear stretches
of protein.

The script will try to do the RightThing™, but it\'ll need help at times,
namely a) if you\'re using a custom forcefield not present under the $GMXLIB
path, and b) if your molecule requires exotic .mdp directives to minimize
properly.

The used VdW and Coulombic interactions are cutoff at 1nm. If potential
shapes/cuttoffs really are important, the .mdp template is readily editable at
the top of the script.

If you\'re using charges and a forcefield with implicit screening (such as
Martini) you\'ll want to set -eps to something else than the default.

By default, constraints in the topology are converted to harmonic bonds with
1e6 kJ/mol/nm^2 force constant.

The minimal .mdp used for minimization is the following:

{mdp}
******************************************

where %%intg%%, %%nsteps%% and %%eps%% are the values specified with -intg,
-nsteps and -eps.

Finally, the script is not that clever that it can read #included files from
the .top/.itp you provide. Make sure that at least the [ moleculetype ] and
[ atoms ] directives are in the file you supply. If they appear more than once
you\'ll get to choose the molecule you want, or you can preempt that with -mol.

Version 3.0.0-05-09-2025 by Manuel Melo (m.n.melo@itqb.unl.pt)'''.format(
        mdp=mdp_template.format(integrator='%intg%',
                                nsteps='%nsteps%',
                                nstxout='*0, or 1 if -keep/-traj is passed*',
                                eps='%eps%',
                                rlist='*dependent on system size*',
                                rcoulomb='*dependent on system size*',
                                rvdw='*dependent on system size*',
                                delta_lambda='1/%nsteps%'))

    builder = MolMaker(description=desc)
    builder.add_argument('-i', dest='itp', type=Path,
                         help='The target molecule topolgy file. Default is '
                              'to ask for it if none given.')
    builder.add_argument('-x', dest='xmdp', type=Path,
                         help='File with extra directives to include in the '
                              'minimization .mdp.')
    builder.add_argument('-ff', dest='ff', type=Path,
                         help='The forcefield itp. Default is to ask for it '
                              'if none given. You\'ll most likely want to use '
                              'this flag if using the MARTINI CG forcefield.')
    builder.add_argument('-o', dest='gro', type=Path,
                         help='The output .gro file. Default is to infer from '
                              'the topology file.')
    builder.add_argument('-mol', default=None,
                         help='Which molecule to convert, in case there are '
                              'multiple in the topology file.')
    builder.add_argument('-temp-prefix', dest='tmpprefix', default='.molmk_',
                         help='The prefix of the temp files.')
    builder.add_argument('-fuzzy', type=float, dest='fuzzy', default=2,
                         help='Fuzziness of the randomly placed '
                              'initial coordinates.')
    builder.add_argument('-traj', action='store_true', dest='traj',
                         help='Whether to save the minimization trajectory '
                              'and a .tpr file.')
    builder.add_argument('-keep', action='store_true', dest='keep',
                         help='Whether to keep all temporary output files '
                              '(note that these will be hidden by default '
                              'anyway).')
    builder.add_argument('-intg', dest="intg", default='cg',
                         help='The energy-minimization integrator to use.')
    builder.add_argument('-nsteps', type=int, dest='nsteps', default=5000,
                         help='Maximum number of steps in minimization '
                              '(might be useful to limit this for '
                              'CG minimizations).')
    builder.add_argument('-eps', type=float, dest='eps', default=1.0,
                         help='The relative dielectric constant.')
    builder.add_argument('-constr', dest='keep_constraints',
                         action='store_true',
                         help='Whether to keep constraints as such, '
                              'instead of converting to harmonic bonds.')

    builder.justdoit(args=argv)

if __name__ == '__main__':
    sys.exit(main())
