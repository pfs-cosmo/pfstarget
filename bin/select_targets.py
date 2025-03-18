#!/usr/bin/env python
import os, sys
import glob 
import numpy as np
from astropy.table import Table

from pfstarget import cuts as Cuts
from pfstarget import extinction as E

from argparse import ArgumentParser
ap = ArgumentParser(description='Generate PFS targets from HSC tract files')
ap.add_argument("tractsdir",
                help="tract file or root directory with tract files")
ap.add_argument("dest",
                help="Output target selection directory")
ap.add_argument("--dust", type=str,
                help='galactic extinction dust model [defaults to desi dust map]',
                default='desi')
ns = ap.parse_args()

infiles = [] 
if os.path.isfile(ns.tractsdir): infiles = glob.glob(ns.tractsdir) 
elif os.path.isdir(ns.tractsdir): infiles = glob.glob('%s/*' % ns.tractsdir) 
else: raise ValueError("no tract files found") 

if len(infiles) == 0:
    raise ValueError("no tract files found") 
    sys.exit(1)

# output file name 
if os.path.isfile(ns.dest): 
    fout = ns.dest
elif os.path.isdir(ns.dest): 
    fout = os.path.join(ns.dest, f'pfs_target.dust_{ns.dust}.fits') 
else: 
    raise ValueError('specify output directory or filename')  

# loop through tract files 
# and select PFS cosmology targets
targets = [] 
for infile in infiles: 
    # read tract file 
    tract = Table.read(infile)

    # preprocess tract file (using specified galactic extinction dust model) 
    _hsc = Cuts._prepare_hsc(tract, dust_extinction=ns.dust)

    # apply PFS cosmology target selection 
    is_pfscosmo = Cuts.isCosmology(_hsc)
    
    # targets 
    targs = _hsc[is_pfscosmo]
    targets.append(targs) 

# combine all targets 
targets = Table(np.concatenate(targets)) 

# write to file 
targets.write(fout) 
