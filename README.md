# pfstarget
This package contains scripts for the PFS Cosmology Survey target selection.

## Installation 
```bash
git clone https://github.com/pfs-cosmo/pfstarget.git
cd pfstarget
pip install -e . 
```

## getting started

1. Download tracts from HSC using the script in `bin/hsc/sql`. You will need a
   [STARS](https://stars2.naoj.hawaii.edu/) account for this. This is the
   imaging catalog that will be used for target selection. 

```bash 
python3 hscReleaseQuery.py s23b_wide_pfs.sql -D --user username -r s23
```

2. Apply target selection. You can do this in python by manually loading the
   tracts.  

```python
from astropy.table import Table
from pfstarget import cuts as Cuts

tract = Table.read('your_tract_file_name.fits') 

is_pfs_cosmo = Cuts.isCosmology(tract)

targets = tract[is_pfs_cosmo]
```

Or, you can use the scripts in `bin/` (COMING SOON). 


## Contribution
If you'd like to contribute, please do so through forking, as described in https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project
