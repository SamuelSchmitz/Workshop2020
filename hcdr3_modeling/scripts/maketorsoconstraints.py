#!python2
import sys
import argparse

usage = "%prog [options] -b -s 95 -e 105"
parser = argparse.ArgumentParser(prog="maketorsoconstraints.py", description="A little script to generate an antibody torso constraint file.")

parser.add_argument("--bulged", "-b", dest="isBulged", action="store_true", default=False, help="Generate bulged constraints.")
parser.add_argument("--start", "-s", dest="intStart", help="Torso starting residue.")
parser.add_argument("--end", "-e", dest="intEnd", help="Torso ending residue.")
args = parser.parse_args()

def main(args):
	if args.isBulged:
		#print "Writing bulged constraints file."
		print "Dihedral C " + str(int(args.intStart)-1) + " N " + args.intStart + " CA " + args.intStart + " C " + args.intStart + " CIRCULARHARMONIC -2.52274631284 0.16437124259"
		print "Dihedral C " + args.intStart + " N " + str(int(args.intStart)+1) + " CA " + str(int(args.intStart)+1) + " C " + str(int(args.intStart)+1) + " CIRCULARHARMONIC -1.76431095224 0.377634488226"
		print "Dihedral C " + str(int(args.intStart)+1) + " N " + str(int(args.intStart)+2) + " CA " + str(int(args.intStart)+2) + " C " + str(int(args.intStart)+2) + " CIRCULARHARMONIC -1.87264994529 0.552906731961"
		print "Dihedral C " + str(int(args.intEnd)-4) + " N " + str(int(args.intEnd)-3) + " CA " + str(int(args.intEnd)-3) + " C " + str(int(args.intEnd)-3) + " CIRCULARHARMONIC -2.11225785093 0.846788046792"
		print "Dihedral C " + str(int(args.intEnd)-3) + " N " + str(int(args.intEnd)-2) + " CA " + str(int(args.intEnd)-2) + " C " + str(int(args.intEnd)-2) + " CIRCULARHARMONIC -1.65408729106 0.602333675851"
		print "Dihedral C " + str(int(args.intEnd)-2) + " N " + str(int(args.intEnd)-1) + " CA " + str(int(args.intEnd)-1) + " C " + str(int(args.intEnd)-1) + " CIRCULARHARMONIC -1.51771540627 0.306977026918"
		print "Dihedral C " + str(int(args.intEnd)-1) + " N " + args.intEnd + " CA " + args.intEnd + " C " + args.intEnd + " CIRCULARHARMONIC -2.19175703236 0.248404020469"
		print "Dihedral N " + args.intStart + " CA " + args.intStart + " C " + args.intStart + " N " + str(int(args.intStart)+1) + " CIRCULARHARMONIC 2.58567239175 0.211900697875"
		print "Dihedral N " + str(int(args.intStart)+1) + " CA " + str(int(args.intStart)+1) + " C " + str(int(args.intStart)+1) + " N " + str(int(args.intStart)+2) + " CIRCULARHARMONIC 2.46937990303 0.221695520757"
		print "Dihedral N " + str(int(args.intStart)+2) + " CA " + str(int(args.intStart)+2) + " C " + str(int(args.intStart)+2) + " N " + str(int(args.intStart)+3) + " CIRCULARHARMONIC 2.39435565136 0.572186942591"
		print "Dihedral N " + str(int(args.intEnd)-2) + " CA " + str(int(args.intEnd)-2) + " C " + str(int(args.intEnd)-2) + " N " + str(int(args.intEnd)-1) + " CIRCULARHARMONIC 1.70230342428 0.451352627783"
		print "Dihedral N " + str(int(args.intEnd)-1) + " CA " + str(int(args.intEnd)-1) + " C " + str(int(args.intEnd)-1) + " N " + args.intEnd + " CIRCULARHARMONIC -0.529477998154 0.449317109408"
		print "Dihedral N " + args.intEnd + " CA " + args.intEnd + " C " + args.intEnd + " N " + str(int(args.intEnd)+1) + " CIRCULARHARMONIC 2.33736630416 0.17136710291"
	else:
		print "This script does not support generation of non-bulged torso constraints at this time."
if __name__ == "__main__":
	main(args)
