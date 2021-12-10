
This folder contains a quick and easy example: generate models for the
barnase-barstar complex and compute all descriptors. 

In the ppdb directory there are the files describing the complex: template and
sequence.  The template (`1brs-cleaned.pdb`) is derived from PDB 1BRS. The
chains A and D have been extracted and used as is. The sequence in `ppdb.seq`
is directly from 1BRS.  The file `ppdb.txt` just list the name of the complex
(matching the name in the sequence file) the number of chains in the ligand and
receptor, the binding affinity (just a placeholder in this case), and the name
of the template file.

The file `config-ppdg.ini` helps configure ppdg, indicating where to find all
executables or scripts. Not all of these are needed, it depends on which
descriptors you want to compute.  The python script `compute_descriptors.py`
read the protein-protein info from the ppdb directory and compute the
descriptors. You can choose which protocol to generate the complex you want to
use, how many models to generate, which descriptors to compute and how (serial,
paralled, parsl) to compute them. `descriptors-all.json` contains sample outputs.

