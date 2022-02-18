#!/usr/bin/env python3

import os
import urllib.request
import logging
log = logging.getLogger(__name__)

PDB_protein = ['ALA', 'CYS', 'GLU', 'ASP', 'GLY',
               'PHE', 'ILE', 'HIS', 'LYS', 'MET',
               'LEU', 'ASN', 'GLN', 'PRO', 'SER',
               'ARG', 'THR', 'TRP', 'VAL', 'TYR',
               'HSD' ]
PDB_dna = ['DA', 'DC', 'DG', 'DT', 'DI' ]
PDB_rna = ['A', 'C', 'G', 'U', 'I' ]
PDB_water = ['HOH']
PDB_ions = ['MG', 'CA', 'ZN', 'CL', 'SO4', 'FES', 'PBM', 'PO4', 'MN', 'ALF', 'CD', 'AF3', 'IOD', 'NA', 'Y1', 'CU', 'NI', 'SR', 'ACT', 'NO3', 'FLC']
PDB_sugars = ['MAN', 'NAG', 'BMA', 'FUC', 'GAL', 'NDG', 'GLA', 'FUL']
PDB_ligands = ['MPD', '3AA', 'PON', 'GOL', 'GNP', 'COT', 'ANP', 'AMP', 'ATP', 'GDP', 'GTP', 'FAD', 'NAD', 'YBT', 'REA', 'B3P', 'PG4', 'PLP', 'HEM', 'GVE', 'TAD', 'GSP', 'MES', 'TXP', 'DKA', 'NPS', 'IPA', 'BTB', 'EDO']
PDB_others = ['ACE']

#AA_3to1 = {
#    'ALA':'A', 'CYS':'C', 'GLU':'E', 'ASP':'D', 'GLY':'G',
#    'PHE':'F', 'ILE':'I', 'HIS':'H', 'LYS':'K', 'MET':'M',
#    'LEU':'L', 'ASN':'N', 'GLN':'Q', 'PRO':'P', 'SER':'S',
#    'ARG':'R', 'THR':'T', 'TRP':'W', 'VAL':'V', 'TYR':'Y',
#    'XXX':'X',
#    'A':'A', 'C':'C', 'T':'T', 'U':'U', 'G':'G',
#    'HSD':'H',
#    }
#AA_1to3 = { v:k for k,v in AA_3to1.items() }


# MN  : MANGANESE (II) ION
# CD  : CADMIUM ION
# FES : FE2/S2 (INORGANIC) CLUSTER
# PBM : TRIMETHYL LEAD ION
# PO4 : PHOSPHATE ION
# ALF : TETRAFLUOROALUMINATE ION
# AF3 : ALUMINUM FLUORIDE
# IOD : IODIDE ION
# NA  : SODIUM ION
# Y1  : YTTRIUM ION
# CU  : COPPER (II) ION
# NI  : NICKEL (II) ION
# SR  : STRONTIUM ION
# ACT : ACETATE ION
# NO3 : NITRATE ION
# FLC : CITRATE ANION

# MPD : (4S)-2-METHYL-2,4-PENTANEDIOL
# 3AA : 3-AMINOPYRIDINE-ADENINE DINUCLEOTIDE PHOSPHATE
# PON : IMIDO DIPHOSPHATE
# GOL : GLYCEROL
# ANP : PHOSPHOAMINOPHOSPHONIC ACID-ADENYLATE ESTER
# GNP : PHOSPHOAMINOPHOSPHONIC ACID-GUANYLATE ESTER
# COT : COA-S-ACETYL TRYPTAMINE
# YBT : BIS-(2-HYDROXYETHYL)AMINO-TRIS(HYDROXYMETHYL)METHANE YTTRIUM
# REA : RETINOIC ACID
# B3P : 2-[3-(2-HYDROXY-1,1-DIHYDROXYMETHYL-ETHYLAMINO)-PROPYLAMINO]-2-HYDROXYMETHYL-PROPANE-1,3-DIOL
# PG4 : TETRAETHYLENE GLYCOL
# PLP : VITAMIN B6 PHOSPHATE
# HEM : HEME
# GVE : METHYL 4-AMINOBUTANOATE
# TAD : BETA-METHYLENE-THIAZOLE-4-CARBOXYAMIDE-ADENINE DINUCLEOTIDE 
# GSP : 5'-GUANOSINE-DIPHOSPHATE-MONOTHIOPHOSPHATE
# MES : 2-(N-MORPHOLINO)-ETHANESULFONIC ACID
# TXP : 1,4,5,6-TETRAHYDRONICOTINAMIDE ADENINE DINUCLEOTIDE PHOSPHATE
# DKA : DECANOIC ACID
# NPS : (2S)-2-(6-METHOXYNAPHTHALEN-2-YL)PROPANOIC ACID
# IPA : ISOPROPYL ALCOHOL
# BTB : 2-[BIS-(2-HYDROXY-ETHYL)-AMINO]-2-HYDROXYMETHYL-PROPANE-1,3-DIOL
# EDO : ETHYLENE GLYCOL

PDB_protein_mod = dict()
PDB_protein_mod['HIC'] = ['HIS', '4-METHYL-HISTIDINE']
PDB_protein_mod['SEP'] = ['SER', 'PHOSPHOSERINE']
PDB_protein_mod['TPO'] = ['THR', 'PHOSPHOTHREONINE']
PDB_protein_mod['TYS'] = ['TYR', 'O-SULFO-L-TYROSINE']
PDB_protein_mod['MSE'] = ['MET', 'SELENOMETHIONINE']
PDB_protein_mod['DDE'] = ['HIS', 'HISTIDINE DERIVATIVE']
PDB_protein_mod['PCA'] = ['GLU', 'PYROGLUTAMIC ACID']
PDB_protein_mod['LLP'] = ['LYS', "LLP N'-PYRIDOXYL-LYSINE-5'-MONOPHOSPHATE"]
PDB_protein_mod['TRQ'] = ['TRP', "TRP-TRP QUINONE"]
PDB_protein_mod['GLX'] = ['GLN', "Glutamine/glutamic-acid"]


class Atom():
    def __init__(self, pdbstring=None, cifstring=None):
        if pdbstring:
            self.pdbstring(pdbstring)
        if cifstring:
            self.cifstring(cifstring)

    def pdbstring(self, line):
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            raise ValueError('Invalid input line!')
        self.name = line[13-1:16].strip()
        self.resname = line[18-1:21].strip()
        self.chain = line[22-1].strip()
        self.resnum = line[23-1:27].strip()
        self.x = float(line[31-1:38])
        self.y = float(line[39-1:46])
        self.z = float(line[47-1:54])
        self.occupancy = float(line[55-1:60])
        self.beta = float(line[61-1:66])
        self.segid = line[73-1:77].strip()

    def cifstring(self, line):
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            raise ValueError('Invalid input line!')
        items = line.split()
        self.name = items[3].replace('"', '')
        self.resname = items[5]
        self.chain = items[18] if len(items[18])==1 else ''
        self.resnum = items[16]
        self.x = float(items[10])
        self.y = float(items[11])
        self.z = float(items[12])
        self.occupancy = float(items[13])
        self.beta = float(items[14])
        self.segid = items[18]


class Pdb():

    def __init__(self, pdb=None):
        self.atoms = list()
        self.seqres = dict()
        if pdb:
            self.read(pdb)

    def __str__(self):
        return "PDB containing %d atoms" % (len(self.atoms))

    def __add__(self, other):
        newpdb = Pdb()
        newpdb.atoms = self.atoms.copy()
        newpdb.atoms.extend(other.atoms)
        newpdb.seqres = self.seqres.copy()
        newpdb.seqres.update(other.seqres)
        return newpdb

    def read(self, fname):
        if fname[-3:]=='pdb':
            with open(fname, 'r') as fp:
                for line in fp.readlines():
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        self.atoms.append(Atom(pdbstring=line))
                    if line.startswith("SEQRES"):
                        self.parse_seqres(line)
                    if line.startswith('END'): ## Stop reading after finding END or ENDMDL
                        break
        elif fname[-3:]=='cif':
            with open(fname, 'r') as fp:
                for line in fp.readlines():
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        self.atoms.append(Atom(cifstring=line))
        else:
            raise ValueError('Cannot read file %s. Allowed filetypes are pdb and cif.' % (fname))

    def convert_res_to_one(self, res):
        '''
            Convert a 3 letter resname to a 1 letter code.
            Keep only protein and DNA/RNA.
            Convert modified amino acids to normal ones.
            Skip all the rest.
        '''
        AA_3to1 = {
            'ALA':'A', 'CYS':'C', 'GLU':'E', 'ASP':'D', 'GLY':'G',
            'PHE':'F', 'ILE':'I', 'HIS':'H', 'LYS':'K', 'MET':'M',
            'LEU':'L', 'ASN':'N', 'GLN':'Q', 'PRO':'P', 'SER':'S',
            'ARG':'R', 'THR':'T', 'TRP':'W', 'VAL':'V', 'TYR':'Y',
            'HSD':'H',
            'DA':'A', 'DC':'C', 'DG':'G', 'DT':'T', 'DU':'U',
            'A':'A' , 'C':'C' , 'G':'G' , 'T':'T' , 'U':'U' ,
            }
        # Normal aminoacids and nucleotides
        if res in PDB_protein+PDB_dna+PDB_rna:
            return AA_3to1[res]

        # Modified aminoacids
        if res in PDB_protein_mod.keys():
            log.info('Found residue %s which I assume is %s and treat is as a %s' % 
                (res, PDB_protein_mod[res][1], PDB_protein_mod[res][0]))
            return AA_3to1[PDB_protein_mod[res][0]]

        # Stuff to skip (water, ions, sugars, ligands, other stuff...)
        if res in PDB_ions+PDB_water+PDB_sugars+PDB_ligands+PDB_others:
            return None

        # ???
        log.warning('Found residue %s, but I do not know what it is. Skipping it.' % (res))
        return None

    def parse_seqres(self, line):
        '''
            Read the sequence from the SEQRES field in PDB.
        '''
        if not line.startswith('SEQRES'):
            log.error('Expected line starting with SEQRES. Got <%s> instead.' % (line))
            return
        chain = line[11]
        seq = ''
        for res in line[19:].split():
            one = self.convert_res_to_one(res)
            if one:
                seq += one
        if not chain in self.seqres.keys():
            self.seqres[chain] = ''
        self.seqres[chain] += seq

    def get_sequence(self):
        '''
            Parse the list of atoms in the Pdb and return the aminoacid sequence.
        '''
        seq = ''
        prev=None
        prevC=self.atoms[0].chain
        for atom in self.atoms:
            curr = atom.resname+atom.resnum+atom.chain
            if curr==prev:
                continue
            else:
                if atom.chain != prevC:
                    seq += '/'
                    prevC = atom.chain
                one = self.convert_res_to_one(atom.resname)
                if one:
                    seq += one
                prev = curr
        return seq

    def write(self, fname, original_resnumber=False):
        if len(self.atoms)<1:
            log.warning("No Atom to write! This Pdb is empty!")
            return
        with open(fname, 'w') as fp:
            res_number = 0
            prev_res = 0
            prev_resname = self.atoms[0].resname
            prev_chain = self.atoms[0].chain
            prev_segid = self.atoms[0].segid
            atom_number = 1
            for atom in self.atoms:
                if atom.chain+atom.segid != prev_chain+prev_segid:
                    fp.write("%-6s%5d %4s %-4s%1s%4d    \n" % 
                        ("TER", atom_number%100000, '', prev_resname, prev_chain, res_number))
                    atom_number += 1
                    prev_chain = atom.chain
                    prev_segid = atom.segid
                if atom.resnum!=prev_res:
                    res_number += 1
                    prev_res = atom.resnum
                    prev_resname = atom.resname
                if len(atom.name)<4:
                    if atom.name[0].isdigit():
                        aname = "%-3s " % (atom.name)
                    else:
                        aname = " %-3s" % (atom.name)
                else:
                    aname = atom.name
                if original_resnumber:
                    fp.write("%6s%5d %4s %-4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s    \n" % 
                        ("ATOM  ", atom_number%100000, aname, atom.resname, atom.chain, atom.resnum, 
                        atom.x, atom.y, atom.z, atom.occupancy, atom.beta, atom.segid))
                else:
                    fp.write("%6s%5d %4s %-4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s    \n" % 
                        ("ATOM  ", atom_number%100000, aname, atom.resname, atom.chain, res_number, 
                        atom.x, atom.y, atom.z, atom.occupancy, atom.beta, atom.segid))
                atom_number += 1
            fp.write("END ")

    def set_occupancy(self, occupancy):
        for atom in self.atoms:
            atom.occupancy = occupancy

    def set_beta(self, beta):
        for atom in self.atoms:
            atom.beta = beta

    def set_beta_residue(self, beta):
        prev = None
        ndx = -1
        for atom in self.atoms:
            curr = atom.resname+atom.resnum+atom.chain
            if curr==prev:
                atom.beta = beta[ndx]
            else:
                ndx += 1
                atom.beta = beta[ndx]
                prev = curr
        print(ndx, len(beta))

    def set_chain(self, chain):
        if chain and len(chain)==1:
            for atom in self.atoms:
                atom.chain = chain
        else:
            log.error("Chain id is not one letter! It is <%s>" % (chain))

    def set_segid(self, segid):
        if segid and len(segid)<=4:
            for atom in self.atoms:
                atom.segid = segid
        else:
            log.error("SegID id is not <= 4 letter! It is <%s>" % (segid))

    def chain2segid(self):
        for atom in self.atoms:
            atom.segid = atom.chain

    def segid2chain(self):
        for atom in self.atoms:
            if len(atom.segid)>0:
                atom.chain = atom.segid[0]
            else:
                atom.chain = ""

    def extract(self, chain=None):
        new = Pdb()
        for atom in self.atoms:
            if atom.chain in chain:
                new.atoms.append(atom)
        return new

    def split_by_chain(self):
        chains = {}
        for atom in self.atoms:
            ch = atom.chain
            if ch not in chains:
                chains[ch] = Pdb()
            chains[ch].atoms.append(atom)
        for ch in chains.keys():
            if ch in self.seqres:
                chains[ch].seqres[ch] = self.seqres[ch]
        return chains

    def split_by_segid(self, waterions=False):
        chains = {}
        for atom in self.atoms:
            if waterions:
                if atom.resname in ['HOH', 'TIP3', 'MG']:
                    atom.segid += 'WT'
            ch = atom.segid
            if ch not in chains:
                chains[ch] = Pdb()
            chains[ch].atoms.append(atom)
        return chains

    def remove_hydrogens(self):
        old = self.atoms
        self.atoms = []
        for atom in old:
            if (atom.name.find('HG')>=0 or 
                atom.name.find('HE')>=0 or
                atom.name.find('HH')>=0 or
                atom.name.find('HD')>=0 or
                atom.name.find('HT')>=0 or
                atom.name.find('HA')>=0 or
                atom.name.find('HB')>=0 or
                atom.name.find('HN')>=0 or
                atom.name.find('HZ')>=0 or
                atom.name in ['H', 'H1', 'H2', 'H3']):
                continue
            else:
                self.atoms.append(atom)

    def remove_ions(self):
        old = self.atoms
        self.atoms = []
        for atom in old:
            if atom.resname in PDB_ions:
                continue
            else:
                self.atoms.append(atom)

    def remove_water_ions(self):
        old = self.atoms
        self.atoms = []
        for atom in old:
            if atom.resname in PDB_ions+PDB_water:
                continue
            else:
                self.atoms.append(atom)

    def keep_protein_only(self):
        old = self.atoms
        self.atoms = []
        for atom in old:
            if atom.resname in PDB_protein:
                self.atoms.append(atom)

    def extract_ca(self):
        old = self.atoms
        self.atoms = []
        for atom in old:
            if atom.name=='CA':
                self.atoms.append(atom)

    def fix4charmm(self):
        for atom in self.atoms:
            # Proteins:
            if atom.name=='1H'  : atom.name = 'HT1'
            if atom.name=='2H'  : atom.name = 'HT2'
            if atom.name=='3H'  : atom.name = 'HT3'
            if atom.name=='H2'  : atom.name = 'HT2'
            if atom.name=='H3'  : atom.name = 'HT3'
            if atom.name=='OXT' : atom.name = 'OT2'
            if atom.name=='OCT1': atom.name = 'OT1'
            if atom.name=='OTC2': atom.name = 'OT2'
            if atom.name=='H'   : atom.name = 'HN'
            if atom.name=='HA3' : atom.name = 'HA1'
            if atom.resname!='ALA' and atom.name=='HB3': atom.name = 'HB1'
            if atom.name=='HG3' and atom.resname in ['GLN', 'LYS', 'ARG', 'GLU', 'PRO'] : atom.name = 'HG1'
            if atom.name=='HG' and atom.resname in ['SER', 'CYS']:  atom.name = 'HG1'
            if atom.resname=='LYS' and atom.name=='HE3': atom.name = 'HE1'
            if atom.resname=='LYS' and atom.name=='HD3': atom.name = 'HD1'
            if atom.resname=='ARG' and atom.name=='HD3': atom.name = 'HD1'
            if atom.resname=='PRO' and atom.name=='HD3': atom.name = 'HD1'
            if atom.resname=='ILE' and atom.name=='CD1': atom.name = 'CD'
            if atom.resname=='ILE' and atom.name=='HD11':atom.name = 'HD1'
            if atom.resname=='ILE' and atom.name=='HD12':atom.name = 'HD2'
            if atom.resname=='ILE' and atom.name=='HD13':atom.name = 'HD3'
            if atom.resname=='CYS' and atom.name=='HG' : atom.name = 'HG3'
            if atom.resname=='SER' and atom.name=='HG' : atom.name = 'HG3'
            if atom.resname=='HIS': atom.resname = 'HSD'
            # Water:
            if atom.resname=='HOH': atom.resname='TIP3'
            if atom.name=='O' and atom.resname=='TIP3': atom.name='OH2'
            # Nucleic:
            if atom.resname=='A'  : atom.resname='ADE'
            if atom.resname=='G'  : atom.resname='GUA'
            if atom.resname=='C'  : atom.resname='CYT'
            if atom.resname=='T'  : atom.resname='THY'
            if atom.resname=='U'  : atom.resname='URA'
            # Phosphate
            if atom.name=='OP1'  : atom.name='O1P'
            if atom.name=='OP2'  : atom.name='O2P'

    def make_standard(self):
        for atom in self.atoms:
            # Proteins atoms:
            if atom.name=='HT1' : atom.name = 'H'
            if atom.name=='HT2' : atom.name = 'H2'
            if atom.name=='HT3' : atom.name = 'H3'
            if atom.name=='OT1' : atom.name = 'O'
            if atom.name=='OT2' : atom.name = 'OXT'
            if atom.name=='OTC2': atom.name = 'OXT'
            if atom.name=='HN'  : atom.name = 'H'
            if atom.name=='HA1' : atom.name = 'HA3'
            if atom.resname!='ALA' and atom.name=='HB1': atom.name = 'HB3'
            if atom.name=='HG1':
                if atom.resname in ['GLN', 'LYS', 'ARG', 'GLU', 'PRO']:
                    atom.name = 'HG3'
                elif atom.resname in ['SER', 'CYS']:
                    atom.name = 'HG'
            if atom.resname=='LYS' and atom.name=='HE1': atom.name = 'HE3'
            if atom.resname=='LYS' and atom.name=='HD1': atom.name = 'HD3'
            if atom.resname=='ARG' and atom.name=='HD1': atom.name = 'HD3'
            if atom.resname=='PRO' and atom.name=='HD1': atom.name = 'HD3'
            if atom.resname=='ILE' and atom.name=='CD':  atom.name = 'CD1'
            if atom.resname=='ILE' and atom.name=='HD1':atom.name = 'HD11'
            if atom.resname=='ILE' and atom.name=='HD2':atom.name = 'HD12'
            if atom.resname=='ILE' and atom.name=='HD3':atom.name = 'HD13'
            if atom.resname=='CYS' and atom.name=='HG3':atom.name = 'HG'
            if atom.resname=='SER' and atom.name=='HG3':atom.name = 'HG'
            if atom.resname=='HSD': atom.resname = 'HIS'

#def pdb2charmm(res, aname, atom_map):
#    res = res.strip()
#    aname = aname.strip()
#    ndx = atom_map[res]['PDB'].index(aname)
#    new_aname = atom_map[res]['XPLOR'][ndx]
#    return res, new_aname

#def _make_resmap():
#    res_map = dict()
#    res_map['PDB'] = [
#    res_map['C36'] = [


#def _read_atomnom(fname='atom_nom.tbl'):
#    """
#        Read atom name table from
#        Data from http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl
#    """
#    atom_map = dict()
#    for k,v in AA_3to1.items():
#        atom_map[k] = dict()
#        atom_map[k]['PDB'] = list()
#        atom_map[k]['XPLOR'] = list()
#    with open(fname, 'r') as fp:
#        for line in fp.readlines():
#            if len(line.strip())==0 or line[0]=='#':
#                continue
#            r1, _, bmrb, sc, pdb, ucsf, msi, xplor, sybyl, midas, diana = line.split('\t')
#            atom_map[AA_1to3[r1]]['PDB'].append(pdb)
#            atom_map[AA_1to3[r1]]['XPLOR'].append(xplor)
#    return atom_map


