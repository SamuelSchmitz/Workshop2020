#!/bin/env python3

from __future__ import print_function

import os
import re
import math

from Bio.PDB.PDBParser import PDBParser
from Bio.SeqUtils import seq1

sctablex = re.compile("([A-Z]{3})_?[A-Z]?:?([A-Za-z]*)_([0-9]+)_?([A-Za-z]?[A-Za-z]?[A-Za-z]?)")

def substract_scoretables(table1, table2, args):
    import re

    if len(table1) != len(table2):
        raise Exception("Scoring tables must have the same y dimension (%s vs %s)" % (len(table1), len(table2)))

    rowk1 = list( k for k,v in table1)
    rowk2 = list( k for k,v in table2)

    mutkeynames = list()
    # Row integrity check
    for ki, (k1, k2) in enumerate(zip(rowk1, rowk2)):
        m1 = sctablex.match(k1)
        m2 = sctablex.match(k2)
        keyname = k1

        if m1 and m2:
            resi = int(m1.group(3))
            resj = int(m2.group(3))
            if resi != resj:
                raise Exception("Row %s contains residue scores for different positions %s and %s" % (ki, resi, resj))
            # Give mutated positions a different key
            resn1 = m1.group(1)
            resn2 = m2.group(1)
            if resn1 != resn2:
                keyname = "%s_%s_%s" % (resn2, resi, resn1)
        elif k1 != k2:
            raise Exception("Keys '%s' and '%s' do not match and are not residue entries" % (k1, k2))
        mutkeynames.append(keyname)

    # Substract
    result = []
    for rowi in range(len(table1)):
        row1, cols1 = table1[rowi]
        row2, cols2 = table2[rowi]
        kset1 = set(cols1.keys())
        kset2 = set(cols2.keys())
        wl = list(kset1.intersection(kset2))
        kdiff = kset1.difference(kset2)
        if len(kdiff):
            print("Discarding %s columns to reduce scoretables to common features: %s" % (len(kdiff), ", ".join(kdiff)))

        rowdat = dict()
        for col1, val1 in cols1.items():
            if not col1 in wl:
                continue

            kn = mutkeynames[rowi]
            val2 = cols2[col1]
            try:
                rowdat[col1] = val1 - val2
            except:
                print("Unable to substract scoretable entries %s.%s: '%s' - '%s'. Taking first" % (row1, col1, val1, val2))
                rowdat[col1] = val1
        result.append( (kn, rowdat) )
    return result

def scoretable2dataframe(scoretable):
    from Bio.SeqUtils import seq1

    from pandas import DataFrame

    data = list()
    for ii, (key, val) in enumerate(scoretable):
        match = sctablex.match(key)
        line = val
        line["key"] = key
        line["resi"] = None
        line["resn"] = None
        if match:
            resi = int(match.group(3))
            resn = seq1(match.group(1))
            tag = match.group(2)
            resm = seq1(match.group(4))
            #if resm:
            #    print(resn, resi, resm)
            line["resm"] = resm
            line["resi"] = resi
            line["resn"] = resn
            line["resname"] = "%s%s%s" % (resn, resi, resm)

        data.append(line)

    frame = DataFrame(data)
    return frame

def scoretable_crop_post_cterminals(table):
    cropped = list()
    for rowk, cols in table:
        cropped.append( (rowk,cols) )
        match = sctablex.match(rowk)
        if match:
            tag = match.group(2)
            if tag == "CtermProteinFull":
                break

    return cropped


def pdb2scoretable(pdbfile, args):
    table = list()
    try:
        with open(pdbfile, 'r') as fh:
            line_istable = False
            header = None
            data = []
            for line in fh:
                line = line.strip().split()
                if not line:
                    continue

                if line[0] == "#BEGIN_POSE_ENERGIES_TABLE":
                    line_istable = True
                elif line[0] == "#END_POSE_ENERGIES_TABLE":
                    line_istable = False
                    break
                elif line_istable:
                    if line[0] == "label":
                        header = line
                    else:
                        data.append(line)

        for row in data:
            resdat = dict(zip(header,row))
            resn = resdat['label']
            del resdat['label']
            for label, value in resdat.items():
                if value == "NA":
                    continue
                score = float(value)
                resdat[label] = score
            table.append( (resn, resdat) )
    except Exception as excptn:
        print("%s: Unable to read Rosetta score table\n%s" % (pdbfile, excptn))
        return None

    return table

def pdb2seq(pdb_path, args):
    pdbinfos = []
    pdb_name = os.path.dirname(os.path.basename(pdb_path))
    pdbparser = PDBParser(PERMISSIVE=1)
    structure = None

    try:
        structure = pdbparser.get_structure(pdb_name, pdb_path)
    except:
        util.perror("%s: Unable to read PDB" % pdb_path, args)
        return [None]

    posenum = 1
    for model in structure.get_models():
        for chain in model.get_chains():
            sequence = ""
            for res in chain.get_residues():
                if res.get_full_id()[3][0] == " ":
                    sequence += seq1(res.get_resname())
            id = "%s%s %s" % (structure.id, model.id or "", chain.id)
            pdbinfo = {"file" : pdb_path,
                       "seq" : sequence,
                       "chain" : chain.id,
                       "posenum" : posenum}
            pdbinfos.append(pdbinfo)
            posenum += len(sequence)
    return pdbinfos

class VariantsByChain:
    def run(self, args, data):
        newdata = dict()

        assigned = False
        while not assigned:
            removeme = list()
            for qdati in data:
                queryi = qdati["name"]
                pdbi = qdati["pdb"]
                for qdatj in data:
                    queryj = qdatj["name"]
                    pdbj = qdatj["pdb"]
                    if pdbi["chain"] == pdbj["chain"] and queryi != queryj:
                        print("Assigning variant '%s' to reference '%s'" % (queryj, queryi))
                        newdata[queryi] = qdati
                        newdata[queryi]['variants'] = {
                            queryj : qdatj
                        }
                        removeme.append(queryj)
                        break
                if removeme:
                    removeme.append(queryi)
                    break
            data = [ e for e in data if e["name"] not in removeme ]
            assigned = not bool(data)
            if not assigned and not removeme:
                raise Exception("Unable to find variant for '%s'" % data[0]["name"])

        return newdata

class AAREUfromPDB:
    def run(self, args, data):
        if data:
            util.pwarn("%s is meant to initialize the data" %self.id(), args)
        data = list()

        for pdb in args["pdbs"]:
            pdbinfos = []
            scoretables = []
            if os.path.isdir(pdb):
                print("Reading directory '%s'..." % (pdb))
                pdbinfos += pdbs2seq(pdb, args)
                scoretables += pdbs2scoretables(pdb, args)
            else:
                pdbinfos += pdb2seq(pdb, args)
                scoretables += [pdb2scoretable(pdb, args)] * len(pdbinfos)

            print("Reading %s pdb contents from '%s'..." % (len(pdbinfos), pdb))

            for filei, pdbinfo in enumerate(pdbinfos):
                if not pdbinfo:
                    continue

                pdbfile = pdbinfo["file"]
                pdbname = os.path.splitext(os.path.basename(pdbfile))[0]
                seq = pdbinfo["seq"]

                name = "%s:%s:%s" % (pdbname, pdbinfo["chain"], pdbinfo["posenum"])
                if name in data:
                    print("Entry '%s' already exists. Skipping..." % name)
                    continue

                if not scoretables[filei]:
                    print("Entry '%s' without scoretable. Skipping..." % name)
                    continue

                if args["crop_post_cterminals"]:
                    set0 = set( k for k,v in scoretables[filei] )
                    scoretables[filei] = scoretable_crop_post_cterminals(scoretables[filei])
                    set1 = set( k for k,v in scoretables[filei] )
                    diff = list(set0.difference(set1))
                    if len(diff):
                        print("Cropped %s post C-terminal(s): %s" % (len(diff), ",".join(diff)))

                REU = scoretables[filei]
                if args["omit_per_residue_scores"]:
                    REU = [ (k,v) for (k,v) in scoretables[filei] if k == "pose" ]

                data.append( {
                    "name" : name,
                    "aa" : seq,
                    "pdb" : pdbinfo,
                    "REU" : REU
                } )

        return data

    @staticmethod
    def create_parser():
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("--omit_per_residue_scores", action="store_true",
                            help="Does not store per residue scoring information but only pose scores")
        parser.add_argument("--crop_post_cterminals", action="store_true",
                            help="If entries exist after the C terminus (like VRT dummies), remove them.")
        parser.add_argument("-p", "--pdbs", nargs="+", required=True,
                            help="List of PDB files or directories containting PDBs.")
        return parser

class VariantsREUPerResidueBreakdown:
    def run(self, args, data):
        title = "%s" % args["title"]
        for queryi, qdata in data.items():
            reference = qdata["aa"]
            refscores = qdata["REU"]
            variants = qdata.get('variants', None)

            rposenum = qdata.get("pdb", {}).get("posenum", None)
            rchainseq = qdata.get("pdb", {}).get("seq", None)
            rchainrange = None
            if rposenum == None or rchainseq == None:
                print("%s: No PDB information found. Unable to check for chain boundaries!" % (queryi))
            else:
                rchainrange = (rposenum-1, rposenum+len(rchainseq)-1)

            if variants:
                for varianti, vdata in variants.items():
                    vposenum = vdata.get("pdb", {}).get("posenum", None)
                    vchainseq = vdata.get("pdb", {}).get("seq", None)
                    vchainrange = None
                    if vposenum == None or vchainseq == None:
                        print("%s: No PDB information for variant %s found. Unable to check for chain boundaries!" % (queryi, varianti))
                    else:
                        vchainrange = (vposenum-1, vposenum+len(vchainseq)-1)

                    if vchainrange != None and rchainrange != None:
                        if vchainrange[0] != rchainrange[0] or vchainrange[1] != rchainrange[1]:
                            print("%s: Incompatible chain ranges. Is the variant '%s' a different protein and not comparable?" % (queryi, varianti))

                    outfile = "%s_%s_%s_%s.svg" % (args["output"] or self.id(), queryi, varianti, title)
                    outfile = outfile.replace(" ", "_")
                    scores = vdata["REU"]
                    try:
                        self.plotPerResidueREU(refscores, scores, outfile, rchainrange, args)
                    except Exception as excptn:
                        print("%s: Unable to plot variant '%s': %s" % (queryi, varianti, excptn))
                        #import traceback
                        #traceback.print_exc()
            else:
                outfile = "%s_%s_%s.svg" % (args["output"] or self.id(), queryi, title)
                outfile = outfile.replace(" ", "_")
                self.plotPerResidueREU(None, refscores, outfile, rchainrange, args)


        return data

    @classmethod
    def plotPerResidueREU(cls, refscores, scores, outfile, chainrange, args):
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from pandas import DataFrame
        import seaborn as sns
        sns.set()

        frame = None

        if refscores:
            table = substract_scoretables(scores, refscores, args)
            frame = scoretable2dataframe(table)

        else:
            frame = scoretable2dataframe(scores)

        #mutant_resi = []
        #mutants = []
        droprows = []
        pos = []
        neg = []

        posenum = 0
        #print(chainrange)
        for rowi, index in frame.iterrows():
            droprow = False
            inboundary = True
            if chainrange != None:
                inboundary = posenum >= chainrange[0] and posenum <= chainrange[1]
            if index["resn"] != None and index["resi"] != None:
                # It appears to be a valid residue, increase posenum!
                posenum += 1

            if index["resn"] == None or index["resi"] == None or not inboundary:
                droprow = True
            else:
                if args["deadzone"] != None and math.fabs(index["total"]) <= args["deadzone"]:
                    if not args["highlight_mutations"] or not index["resm"]:
                        droprow = True
                if args["resi"]:
                    droprow = not int(index["resi"]) in args["resi"]
            #if index["resm"] and not droprow:
            #    mutants.append(rowi-len(droprows))
            #    print(mutants[-1])
            #    mutant_resi.append(int(index["resi"]))
            #if not droprow:
            #    print(posenum, index["resi"], index["resn"], index["resm"], droprow, inboundary)

            if droprow:
                droprows.append(rowi)
            else:
                if index["total"] > 0:
                    pos.append(int(index["resi"]))
                elif index["total"] < 0:
                    neg.append(int(index["resi"]))
        frame = frame.drop(droprows)
#        if not frame.empty:
#            raise Exception("Nothing to plot!")

        col_key = "Score Term"
        col_val = "Score"
        xcol = "total"
        hue = None
        if args["breakdown"] != None:
            tmp = frame.set_index("key")
            cols_keep = ["resi", "resm", "resn", "resname"]
            cols_collapse = set(tmp.keys()) - set(cols_keep)

            breakdown = args["breakdown"]
            deadzone = None
            if breakdown:
                if len(breakdown) == 1:
                    try:
                        deadzone = float(breakdown[0])
                    except:
                        pass
                if deadzone == None:
                    cols_collapse = cols_collapse.intersection(set(breakdown))

            rows = list()
            for rowi, index in tmp.iterrows():
                for collapse in cols_collapse:
                    cells = dict()
                    for keep in cols_keep:
                        cells[keep] = index[keep]
                    cells[col_key] = collapse
                    cells[col_val] = index[collapse]

                    if deadzone == None or math.fabs(cells[col_val]) > deadzone:
                        rows.append(cells)
            frame = DataFrame(rows)
            hue = col_key
            xcol = col_val

#        if not frame.empty:
#            raise Exception("Nothing to plot!")

        mutants = []
        mutant_resi = []
        plotrow = 0

        try:
            for ii, (rowi, index) in enumerate(frame.drop_duplicates(subset=["resi", "resn", "resm", "resname"]).iterrows()):
                if index["resm"]:
                    mutants.append(ii)
                    mutant_resi.append(int(index["resi"]))
        except ValueError as error:
            pass
            # No duplicates

        print("%s improved residue(s): %s" % (len(neg), ",".join(map(str, sorted(neg)))))
        print("%s worsened residue(s): %s" % (len(pos), ",".join(map(str, sorted(pos)))))
        print("%s mutation(s): %s" % (len(mutant_resi), ",".join(map(str, sorted(mutant_resi)))))

        title = "%s" % args["title"]
        plt.clf()

        fontsize_pt = plt.rcParams['ytick.labelsize']
        dpi = 52.25
        matrix_height_pt = fontsize_pt * frame.shape[0]
        if not matrix_height_pt:
            raise Exception("Empty plot!")
        matrix_height_in = matrix_height_pt / dpi
        top_margin = 24./matrix_height_pt  # in percentage of the figure heigh
        bottom_margin = 24./matrix_height_pt  # in percentage of the figure height
        figure_height = matrix_height_in / (1 - top_margin - bottom_margin)
        fig, ax = plt.subplots(
            figsize=(8,figure_height),
            gridspec_kw=dict(top=1-top_margin, bottom=bottom_margin,
                             right=0.70))

        sns.set(font_scale=1.)
        ax = sns.barplot(data=frame, x=xcol, y="resname",
                         palette=sns.color_palette("tab20", 20) if hue else None,
                         hue=hue, hue_order=sorted(set(frame[hue])) if hue else None)
        ax.set_xlabel("%sPer Residue Score [REU]" % ("Relative " if refscores else ""))
        ax.set_ylabel("Residue")
        ax.tick_params(axis='both', labeltop=True)
        for mutant in mutants:
            ax.get_yticklabels()[mutant].set_color("red")
        plt.suptitle(title)

        if args["breakdown"]:
            from matplotlib import ticker
            minorLocator = ticker.AutoMinorLocator(2)
            ax.yaxis.set_minor_locator(minorLocator)
            #plt.tick_params(which='both', width=2)
            plt.grid(b=True, which='minor', color='w',
                    linestyle='-', linewidth=10)

            leglow = plt.legend(bbox_to_anchor=(1.05, 0), loc=3,
                                borderaxespad=0.5, frameon=True,
                                title=col_key, fancybox=True)
            leglow.get_frame().set_facecolor('none')
            leglow.get_frame().set_linewidth(1.0)
            legupp = plt.legend(bbox_to_anchor=(1.05, 1), loc=2,
                                borderaxespad=0.5, frameon=True,
                                title=col_key, fancybox=True)
            legupp.get_frame().set_facecolor('none')
            legupp.get_frame().set_linewidth(1.0)

            if len(set(frame["resname"])) > 8:
                plt.gca().add_artist(leglow)
        else:
            ax.grid(True)

        print("Saving '%s'..." % outfile)
        plt.savefig(outfile, dpi=600)

    @staticmethod
    def create_parser():
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("--title", default="",
                            help="Title to add to the output image")
        parser.add_argument("--deadzone", type=float,
                            help="Hide bars with absolute values <= provided positive float value")
        parser.add_argument("--highlight_mutations", action="store_true",
                            help="Highlight and show mutated residues of variants even if they are within the deadzone")
        parser.add_argument("--breakdown", nargs="*", default=None,
                            help="Display scoring term contributions. Accepts either a list of scoring term names or a single deadzone cutoff (positive float) to automatically hide and display individual scoring terms")
        parser.add_argument("--resi", nargs="+", type=int,
                            help="List of residues to be printed even if they are within the deadzone")
        return parser

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input pose")
    parser.add_argument("design", help="ROSETTA design")
    #parser.add_argument("-o", "--output", required=True, help="Output file")
    parser.add_argument("-r", "--resi", required=True, nargs="+", type=int,
                        help="Pose numbers to plot")
    args = parser.parse_args()

    args = {
        "title" : "",
        "pdbs" : [args.input, args.design],
        "crop_post_cterminals" : True,
        "omit_per_residue_scores" : False,
        "highlight_mutations" : True,
        "output" : "PerResidueEnergies",
        "breakdown" : None,
        "deadzone" : 0.0,
        "resi" : args.resi
    }

    data = dict()
    reader = AAREUfromPDB()
    relative = VariantsByChain()
    plotter = VariantsREUPerResidueBreakdown()

    data = reader.run(args, data)
    data = relative.run(args, data)
    data = plotter.run(args, data)

    args["output"] = "PerResidueEnergiesBreakdown"
    args["breakdown"] = [0.1]
    plotter.run(args, data)

    #from IPython import embed; embed()
