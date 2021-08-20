#!/usr/bin/env python2.7
import sys
import amino_acids
import Bio
from Bio.Align.Applications import ClustalwCommandline as cw
import warnings
from multiprocessing import Pool
import argparse
import os
from Bio import AlignIO
import rosettaScore_beta as rsb
from Bio.PDB import PDBExceptions
from rosettautil.protein import util, pdbStat
from rosettautil import resfile
import glob

clustalw_exe = r"/sb/apps/clustalw/Linux2/i686/v1.83/clustalw"

warnings.simplefilter('ignore', PDBExceptions.PDBConstructionWarning)

if sys.version_info < (2, 7):
    raise Exception("You must use python2.7 to run this")

# usage
usage = "%prog [options] -n native.pdb *pdbs"
parser = argparse.ArgumentParser(prog="score_vs_rmsd_full.py", description="A Rosetta score vs. rmsd script. Pass a resfile to superimpose and calculate score vs. rmsd of just the listed residues. \
	Two output files will be generated. One file will be the whole model aligned. The other will be \
    	just the resfile aligned.")

parser.add_argument("pdbs", metavar="*pdb", nargs=argparse.REMAINDER, help="The pdb files that you want superimposed")
parser.add_argument("--native", "-n", dest="native", help="Native structure")
parser.add_argument("--out", "-o", dest="table", help="A prefix for the output file", default="score_vs_rmsd")
parser.add_argument("--res", "-r", dest="residues", help="The res file to use for localized alignment.", default="")
args = parser.parse_args()

if len(args.pdbs) < 1:
    parser.error("specify at least 1 protein to compare to native")

def get_fasta_from_pdb(pdb):
    holder = []
    for i in pdb.get_residues():
        holder.append((i,amino_acids.longer_names[i.get_resname()]))
    return holder

def get_pairwise_alignment(native,target):
    alignment_dict = {}
    fasta1 = get_fasta_from_pdb(native)
    fasta2 = get_fasta_from_pdb(target)

    with open("clustal_tmp_in.fasta",'w') as f:
        f.write(">native\n")
        for letter1 in fasta1:
            f.write(letter1[1])
        f.write("\n>target\n")
        for letter2 in fasta2:
            f.write(letter2[1])
    cline = cw(clustalw_exe,infile="clustal_tmp_in.fasta")
    assert os.path.isfile(clustalw_exe), "Clustal W binary is missing"
    stdout,stderr = cline()
    alignment = AlignIO.read('clustal_tmp_in.aln','clustal')
    alignment_dict['native'] = str(alignment[0,:].seq)
    alignment_dict['target'] = str(alignment[1,:].seq)

    native_str = str(alignment[0,:].seq)
    target_str = str(alignment[1,:].seq)

    print("Native: " + native_str + " (" + str(len(fasta1)) + ")")
    print("Target: " + target_str + " (" + str(len(fasta2)) + ")")

    use_native = []
    #nat_len = 0
    use_target = []
    #tar_len = 0
    for (native_res, target_res) in zip(native_str, target_str):
        if target_res<>"-" and native_res<>"-":
            use_native.append(True)
            use_target.append(True)
            #nat_len += 1
            #tar_len += 1
        elif target_res=="-" and native_res<>"-":
            use_native.append(False)
        elif target_res<>"-" and native_res=="-":
            use_target.append(False)
        else:
            #Do nothing
            continue

    #print "Native: " + str(nat_len)
    #print "Target: " + str(tar_len)

    return [use_native, use_target]



def main(args):
    tag_list = []
    for i in args.pdbs:
     # stats is a dictionary
        stats = find_rmsd(i)
        tag_list.append(stats)
    make_table(tag_list, args.table)
   


def make_table(name_score_rmsd, out_name):
    additional_headers = name_score_rmsd[0].values()[0]["pose_score"].keys()
    header = ["MODEL", "CA_RMSD"] + additional_headers + ["\n"]
    with open(out_name + "_align_all_model.tsv", 'w') as f:
        f.write("\t".join(header))
        for entity in name_score_rmsd:
            for model_entity in entity:
                model = model_entity
                model_rms = entity[model]["rms_all"]
                pose_score = entity[model]["pose_score"]
                f.write(model + "\t")
                f.write(str(model_rms) + "\t")
                #f.write("{0:.2f}\t".format(all_rmsds["ca_rmsd"]))
                #f.write("{0:.2f}\t".format(all_rmsds["bb_rmsd"]))
                #f.write("{0:.2f}\t".format(all_rmsds["all_rmsd"]))
                for scores in pose_score:
                    f.write("{0:.2f}\t".format(pose_score[scores]))
                #f.write(str(pose_score) + "\t")
                f.write("\n")
    if args.residues:
        #print "Again, resfiles aren't supported by this script!"
        additional_headers = name_score_rmsd[0].values()[0]["residue_scores"].keys()
        header = ["MODEL", "CA_RMSD"] + additional_headers + ["\n"]
        with open(out_name + "_align_by_residue.tsv", 'w') as f:
            f.write("\t".join(header))
            for entity in name_score_rmsd:
                for model_entity in entity:
                    model = model_entity
                    residue_rmsds = entity[model]["rms_res"]
                    residue_scores = entity[model]["residue_scores"]
                    f.write(model + "\t")
                    f.write(str(residue_rmsds) + "\t")
                    # f.write("{0:.2f}\t".format(residue_rmsds["ca_rmsd"]))
                    # f.write("{0:.2f}\t".format(residue_rmsds["bb_rmsd"]))
                    # f.write("{0:.2f}\t".format(residue_rmsds["all_rmsd"]))
                    for scores in residue_scores:
                        f.write("{0:.2f}\t".format(residue_scores[scores]))
                    f.write("\n")

def find_rmsd(model):
    name_score_rmsd = {}
    native = args.native
    rms_all = []
    rms_res = []
    residue_scores = {}
    
    # Start the parser
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
     
    # Get the structures
    ref_structure = pdb_parser.get_structure("reference", native)
    sample_structure = pdb_parser.get_structure("sample", model)
     
    # Use the first model in the pdb-files for alignment
    # Change the number 0 if you want to align to another structure
    ref_model    = ref_structure[0]
    sample_model = sample_structure[0]
     
    # Make a list of the atoms (in the structures) you wish to align.
    # In all cases:
    print "File: " + model
    pwa = get_pairwise_alignment(util.load_pdb(args.native), util.load_pdb(model))
    ref_atoms = []
    sample_atoms = []  

    # When we have a resfile:
    if args.residues:
        residues = resfile.Residue_File(args.residues).get_designed_entities()
        #print residues
        ref_res_atoms = []
        sample_res_atoms = []

    # Iterate of all chains in the model in order to find all residues
    n = 0
    for ref_chain in ref_model:
      # Iterate of all residues in each model in order to find proper atoms
      for ref_res in ref_chain:
        # Check if residue number ( .get_id() ) is in the list
        if pwa[0][n]:
          # Append CA atom to list
          ref_atoms.append(ref_res['CA'])
        # If we have a resfile, look to see if this residue is also in that list:
        if args.residues:
        	#print (ref_chain.id, ref_res.get_id()[1])
            if (ref_chain.id, ref_res.get_id()[1]) in residues:
                ref_res_atoms.append(ref_res['CA'])
        n += 1
    
    n = 0
    for sample_chain in sample_model:
      # Iterate of all residues in each model in order to find proper atoms
      for sample_res in sample_chain:
        # Check if residue number ( .get_id() ) is in the list
        if pwa[1][n]:
          # Append CA atom to list
          sample_atoms.append(sample_res['CA'])
        # If we have a resfile, look to see if this residue is also in that list:
        if args.residues:
            if (sample_chain.id, sample_res.get_id()[1]) in residues:
                sample_res_atoms.append(sample_res['CA'])
        n += 1

    #for (rres, sres) in zip(ref_atoms, sample_atoms):
    #    print str(rres.parent) + " vs " + str(sres.parent)

    # Now we initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
     
    # Print RMSD:
    rms_all = super_imposer.rms
    print "CA_RMSD: " + str(rms_all)
    #sys.exit()

    # will return three rmsds aligned by everything in decoy
    # give it an empty residue file and it will autmatically calculate all rmsds
    # rms_all = pdbStat.calculate_all_superpositions(util.load_pdb(native), util.load_pdb(model), [], args.debug,pwa)
    # print "RMSD: " + str(rms_all)
    
    score_table = rsb.ScoreTable(model)
    # dictionary with all the pose scores
    pose_score = score_table.get_pose_all_scores()
    if args.residues:
        #print "This code does not support res files yet!"
        #rms_res = 0.0
        #residue_scores = 0
        # dictionary with all the residue scores
        
        si2 = Bio.PDB.Superimposer()
        si2.set_atoms(ref_res_atoms, sample_res_atoms)
        si2.apply(sample_model.get_atoms())
        rms_res = si2.rms
        print "RES CA_RMSD: " + str(rms_res)

        #rms_res = pdbStat.calculate_all_superpositions(util.load_pdb(native), util.load_pdb(model), residues, args.debug)
        for i in residues:
            score_table_per_residue = score_table.get_all_score_terms(chain=i[0], pdbres=int(i[1]),)
            for i in score_table_per_residue.iterkeys():
                try:
                    residue_scores[i] += score_table_per_residue[i]
                except KeyError:
                    residue_scores[i] = score_table_per_residue[i]

    name_score_rmsd[model] = {"rms_all": rms_all, "rms_res": rms_res,
                              "pose_score": pose_score, "residue_scores": residue_scores}
    return name_score_rmsd


if __name__ == "__main__":
    main(args)
