#! /usr/bin/env python2.7

usage = '''
Usage: compare_design_to_control_multistate.py PDBID control.fasc design.fasc

Plots metrics of a designed protein-protein interface vs a native interface in a multistate design problem
'''

import matplotlib
matplotlib.use('Agg')
import sys
import matplotlib.pyplot as plt
plt.style.use('ggplot')

if len(sys.argv) < 4:
	print(usage)
	exit()

script, pdbid, control_sc, design_sc = sys.argv

header_line = open( control_sc ).readlines()[1].split()
score_index = header_line.index('total_score')
ddg_index = header_line.index('dG_separated')
dens_index = header_line.index('dG_separated/dSASAx100')

boxplot_args = { 'showcaps':False, 'sym':'ko', 'whiskerprops':{'color':'black', 'ls':'-'}, 
    'boxprops':{'color':'black'}, 'medianprops':{'color':'blue'}, 'whis':'range' }

control_scores = []
control_ddgs = []
control_dens = []
for line in open( control_sc ).readlines()[2:]:
	split = line.split()
	if pdbid in split[-1]:
		control_scores.append( float(split[score_index]) )
		control_ddgs.append( float(split[ddg_index]) )
		control_dens.append( float(split[dens_index]) )

header_line = open( design_sc ).readlines()[1].split()
score_index = header_line.index('total_score')
ddg_index = header_line.index('dG_separated')
dens_index = header_line.index('dG_separated/dSASAx100')

design_scores = []
design_ddgs = []
design_dens = []
for line in open( design_sc ).readlines()[2:]:
	split = line.split()
	if pdbid in split[-1]:
		design_scores.append( float(split[score_index]) )
		design_ddgs.append( float(split[ddg_index]) )
		design_dens.append( float(split[dens_index]) )


plt.figure(1)
plt.boxplot( [control_scores, design_scores], \
	labels=('Control', 'Design'), **boxplot_args )
plt.ylabel( 'Score (REU)' )
plt.title(pdbid+' control vs design - Score')
plt.savefig( pdbid+'_control_vs_design_score.png' )

plt.figure(2)
plt.boxplot( [[i for i in control_ddgs], [i for i in design_ddgs]], \
	labels=('Control', 'Design'), **boxplot_args )
plt.ylabel( 'Binding energy (REU)' )
plt.title(pdbid+' control vs design - Binding energy')
plt.savefig( pdbid+'_control_vs_design_ddg.png' )

plt.figure(3)
plt.boxplot( [control_dens, design_dens], \
	labels=('Control', 'Design'), **boxplot_args )
plt.ylabel( 'Binding density (REU/A2)' )
plt.title(pdbid+' control vs design - Binding density')
plt.savefig( pdbid+'_control_vs_design_density.png' )

