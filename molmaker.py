#!/usr/bin/env python

###  Copyright 2016 Manuel N. Melo ###
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

# Functions
############

def find_exec(name, error='error', extra=""):
    proc = subprocess.Popen(["which", name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    path = proc.communicate()[0]
    if path:
        path = path.strip()
        if os.path.islink(path):
            path = os.readlink(path)
        return path
    else:
        if error in ('error', 'warning'):
            sys.stderr.write("%s: can't find %s in $PATH.\n" % (error, name))
            if extra:
                sys.stderr.write("%s\n" % (extra,))
            if error == 'error':
                sys.exit()
        return False


def deduplicate_mdp(mdp_string):
    """Only keep the last definition of each option.

    Take the MDP file as a string, and returns the new MDP file as a string.
    Only the last definition of each option is kept. Comments and empty lines
    are lost in the process.
    """
    # Save the options in a dictionary
    mdp_dict = {}
    for line in mdp_string.split('\n'):
        equal_pos = line.find('=')
        if equal_pos != -1 and not line.lstrip().startswith(';'):
            key = line[:equal_pos].strip()
            value = line[equal_pos + 1:].strip()
            mdp_dict[key] = value
    # Format the options back to MDP format
    key_length = max(len(key) for key in mdp_dict.keys())
    template = '{key:<{length}} = {value}'
    lines = (template.format(key=key, value=value, length=key_length)
             for key, value in mdp_dict.items())
    return '\n'.join(lines)


# Classes
##########

class ProperFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """A hackish class to get proper help format from argparse.

    """
    def __init__(self, *args, **kwargs):
        super(ProperFormatter, self).__init__(*args, **kwargs)

class MolMaker(argparse.ArgumentParser):
    def __init__(self, formatter_class=ProperFormatter, *args, **kwargs):
        argparse.ArgumentParser.__init__(self, *args, formatter_class=formatter_class, **kwargs)

    def checkopts(self, args=sys.argv):
        self.values = self.parse_args(args)
        #
        # Input options checking and assignment of defaults not done by optparse.
        if not self.values.itp:
            self.values.itp = raw_input ("Topology file: ")
        if not os.path.exists(self.values.itp):
            sys.stderr.write("Error: can't find the specified topology file.\n")
            sys.exit()

        if not self.values.gro:
            self.values.gro = re.sub("(\.itp|.top)?$", ".gro", os.path.basename(self.values.itp))
        self.name = re.match("(.*)\.gro$", self.values.gro).groups()[0]
        self.basename = self.values.tmpprefix+self.name
            
        if self.values.fuzzy < 0:
            self.values.fuzzy = 0

        # Are we running gromacs 5?
        self.gmx = find_exec('gmx', error=None)
        if self.gmx:
            self.grompp = [self.gmx, "grompp"]
            self.mdrun = [self.gmx, "mdrun"]
            self.trjconv = [self.gmx, "trjconv"]
        else:
            # No gmx executable. Let's try the gromacs 4 names.
            self.grompp = [find_exec('grompp')]
            self.mdrun = [find_exec('mdrun')]
            self.trjconv = [find_exec('trjconv', error='warning', extra="Will not center output coordinates nor output minimization trajectory (if requested).\n")]

        # Version checking
        self.gmxversion = subprocess.Popen(self.grompp + ["-version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.gmxversion = self.gmxversion.communicate()[0]
        self.gmxversion = float(re.search("(\d+\.\d+)\.\d+\D", self.gmxversion ).groups()[0])
        if self.gmxversion < 4.5:
            sys.stderr.write("Error: GMX version must be 4.5 or higher.\n")
            sys.exit()

        if not self.values.itp:
            self.values.itp = raw_input ("Topology file: ")
        
        if not self.values.ff:
            try:
                self.gmxlibdir = os.environ["GMXLIB"] 
            except KeyError:
                self.gmxlibdir = re.sub("/bin/.*$", "/share/gromacs/top", self.grompp[0])
                if not os.path.exists(self.gmxlibdir):
                    sys.stderr.write("Error: forcefield not specified and I can't find $GMXLIB.\n")
                    sys.exit()
            tops = os.listdir(self.gmxlibdir)
            ffs = []
            for i in tops:
                if re.search('\.ff$', i):
                    ffs.append(i)
            if not ffs:
                sys.stderr.write("Error: forcefield not specified and I can't find any in $GMXLIB.\n")
                sys.exit()
            for ff,i in enumerate(ffs):
                sys.stdout.write("%3i: %s\n" % (ff, i[:-3]))
            chosenff = int(raw_input ("Forcefield to use (number only): "))
            chosenff = ffs[chosenff % len(ffs)]
            self.values.ff = "%s/%s/forcefield.itp" % (self.gmxlibdir,chosenff)

    
    def getmolname(self):
        itp=open(self.values.itp, "r")
        self.molname = ""
        while not re.match("\s*\[\s*moleculetype\s*\]" , itp.readline()):
            pass
        while not self.molname:
            self.molname = re.match("\s*(\w+)\s+" , itp.readline())
        self.molname = self.molname.groups()[0]
        itp.close()

    def creategro(self):
        itp=open(self.values.itp, "r")
        line = ""
        atoms = 0
        while not re.match("\s*\[\s*atoms\s*\]" , itp.readline()):
            pass
        while not re.match("\s*\[\s*" , line):
            if re.match("\s*\d",line):
                atoms += 1
            line = itp.readline()
        itp.close()
        
        self.gro = "%s.gro" % (self.basename)
        gro = open (self.gro,"w")
        gro.write("%s\n" % (self.molname))
        gro.write("%i\n" % (atoms))
        for i in range(0,atoms):
            gro.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (i+1,"DUM","DUM",i+1,0.05*i,self.values.fuzzy*random.random()+0.025*atoms-0.5*self.values.fuzzy,self.values.fuzzy*random.random()+0.025*atoms-0.5*self.values.fuzzy))
        gro.write("%f %f %f" % (0.2*atoms, 0.2*atoms, 0.2*atoms))
        gro.close()
    
    def createtop(self):
        template = """#include \"%s\"
#include \"%s\"
[ system ]
%s
[ molecules ]
%s   1
"""
        self.top = "%s.top" % (self.basename)
        top = open (self.top,"w")
        top.write(template % (self.values.ff, self.values.itp, self.molname, self.molname))
        top.close()

    def createmdp(self):
        template = mdp_template % (self.values.intg, str(self.values.nsteps), str(self.values.eps))
        if self.gmx: #GMX >= 5.0
            template += mdp_template_gmx5
        self.mdp = "%s.mdp" % (self.basename)
        with open (self.mdp,"w") as mdp:
            if self.values.xmdp:
                template += "".join(open(self.values.xmdp).readlines())
            template = deduplicate_mdp(template)
            mdp.write(template)

    def minimize(self):
        self.tpr = "%s.tpr" % (self.basename)
        self.mdout = "%s.out.mdp" % (self.basename)
        self.pplog = "%s.pp.log" % (self.basename)
        self.deffnm = "%s.md" % (self.basename)
        self.mdgro = "%s.gro" % (self.deffnm)
        self.mdtrr = "%s.trr" % (self.deffnm)
        self.mdlog = "%s.log" % (self.deffnm)
        self.mdedr = "%s.edr" % (self.deffnm)

        ppargs = self.grompp + ["-f",self.mdp, "-p",self.top, "-c",self.gro, "-maxwarn","3", "-po",self.mdout, "-o",self.tpr]
        mdargs = self.mdrun + ["-s",self.tpr, "-nt","1", "-deffnm",self.deffnm, "-cpt","0"]

        with open(self.pplog, "w") as pplog:
            md_env = os.environ.copy()
            md_env["GMX_SUPPRESS_DUMP"] = '1'
            md_env["GMX_MAXBACKUP"] = '-1'
            pp = subprocess.call(ppargs, stdout=pplog, stderr=pplog, env=md_env)
        if pp:
            sys.stderr.write("grompp error: check %s\n" % (self.pplog))
            sys.exit()
        md = subprocess.call(mdargs, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=md_env)
        if md:
            sys.stderr.write("mdrun error: check %s\n" % (self.mdlog))
            sys.exit()

    def cleanup(self):
        self.outtpr = "%s.tpr"%(self.name)
        self.outtrr = "%s.trr"%(self.name)
        if self.trjconv[0]:
            tcargs = self.trjconv + ["-f",self.mdgro, "-o",self.values.gro, "-s",self.tpr, "-center", "-pbc","mol"]
            trjc = subprocess.Popen(tcargs, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            trjc.communicate("0\n0\n")
            tc = trjc.returncode
            if tc:
                sys.stderr.write("Warning: could not center the final structure.\n")
                shutil.copy(self.mdgro, self.values.gro)
            if self.values.traj:
                tcargs = self.trjconv + ["-f",self.mdtrr, "-o",self.outtrr, "-s",self.tpr, "-center", "-pbc","mol"]
                trjc = subprocess.Popen(tcargs, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                trjc.communicate("0\n0\n")
                tc = trjc.returncode
                if tc:
                    sys.stderr.write("Warning: could not center the trajectory.\n")
        else:
            shutil.copy(self.mdgro, self.values.gro)
        if self.values.traj:
            shutil.copy(self.tpr,self.outtpr)
            
        if not self.values.keep:
            os.remove(self.gro)
            os.remove(self.top)
            os.remove(self.mdp)
            os.remove(self.mdout)
            os.remove(self.pplog)
            os.remove(self.tpr)
            os.remove(self.mdgro)
            os.remove(self.mdtrr)
            os.remove(self.mdlog)
            os.remove(self.mdedr)

    def justdoit(self):
        self.getmolname()
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
The minimal .mdp used for minimization is the following:

%s
%s
******************************************

where %%intg%%, %%nsteps%% and %%eps%% are the values specified with -intg, -nsteps and -eps.

Finally, the script is not that clever that it can read #included files from the .top/.itp you provide. Make sure that at least the [ moleculetype ] and [ atoms ] directives are in the file you supply (If they appear more than once you'll get the first molecule in the file).
%%prog [options]\n
Version 1.4.0-15-02-2016 by Manuel Melo (m.n.melo@rug.nl)""" % (mdp_template%("%intg%","%nsteps%","%eps%"), mdp_template_gmx5)

    builder = MolMaker(description=desc)
    builder.add_argument("-i", dest="itp", help="The target molecule topolgy file. Default is to ask for it if none given.")
    builder.add_argument("-x", dest="xmdp", help="File with extra directives to include in the minimization .mdp.")
    builder.add_argument("-ff", dest="ff", help="The forcefield itp. Default is to ask for it if none given. You'll most likely want to use this flag if using the MARTINI CG forcefield.")
    builder.add_argument("-o", dest="gro", help="The output .gro file. Default is to infer from the topology file.")
    builder.add_argument("-temp-prefix", dest="tmpprefix", default=".molmk_", help="The prefix of the temp files.")
    builder.add_argument("-fuzzy", type=float, dest="fuzzy", default=2, help="Fuzziness of the randomly placed initial coordinates.")
    builder.add_argument("-keep", action="store_true", dest="keep", help="Whether to keep all temporary output files (note that these will be hidden by default anyway).")
    builder.add_argument("-traj", action="store_true", dest="traj", default=False, help="Whether to save the minimization trajectory and a .tpr file.")
    builder.add_argument("-intg", dest="intg", default="cg", help="The energy-minimization integrator to use.")
    builder.add_argument("-nsteps", type=int, dest="nsteps", default=5000, help="Maximum number of steps in minimization (might be useful to limit this for CG minimizations).")
    builder.add_argument("-eps", type=float, dest="eps", default=1.0, help="The relative dielectric constant.")
    
    builder.checkopts(args=argv)
    builder.justdoit()

if __name__ == "__main__":
    sys.exit(main())


