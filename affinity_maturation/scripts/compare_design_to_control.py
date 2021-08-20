#! /usr/bin/env python2.7

usage = '''
Usage: compare_design_to_control.py control.fasc design.fasc

Plots metrics of a designed protein-protein interface vs a native interface
'''

import sys
import matplotlib.pyplot as plt
plt.style.use('ggplot')

if len(sys.argv) < 3:
	print(usage)
	exit()

script, control_sc, design_sc = sys.argv

header_line = open( control_sc ).readlines()[1].split()
score_index = header_line.index('total_score')
ddg_index = header_line.index('dG_separated')
dens_index = header_line.index('dG_separated/dSASAx100')

boxplot_args = { 'showcaps':False, 'sym':'ko', 'whiskerprops':{'color':'black', 'ls':'-'}, 
    'boxprops':{'color':'black'}, 'medianprops':{'color':'blue'}, 'whis':'range' }

control_scores, control_ddgs, control_dens = zip(* [ (float(line.split()[score_index]), float(line.split()[ddg_index]), \
	float(line.split()[dens_index])) \
	for line in open( control_sc ).readlines()[2:] ] )
design_scores, design_ddgs, design_dens = zip(* [ (float(line.split()[score_index]), float(line.split()[ddg_index]), \
	float(line.split()[dens_index])) \
	for line in open( design_sc ).readlines()[2:] ] )

plt.figure(1)
plt.boxplot( [control_scores, design_scores], \
	labels=('Control', 'Design'), **boxplot_args )
plt.ylabel( 'Score (REU)' )
plt.title('Control vs design - Score')
plt.savefig( 'control_vs_design_score.png' )

plt.figure(2)
plt.boxplot( [control_ddgs, design_ddgs], \
	labels=('Control', 'Design'), **boxplot_args )
plt.ylabel( 'Binding energy (REU)' )
plt.title('Control vs design - Binding energy')
plt.savefig( 'control_vs_design_ddg.png' )

plt.figure(3)
plt.boxplot( [control_dens, design_dens], \
	labels=('Control', 'Design'), **boxplot_args )
plt.ylabel( 'Binding density (REU/A2)' )
plt.title('Control vs design - Binding density')
plt.savefig( 'control_vs_design_density.png' )

