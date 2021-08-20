~/rosetta_workshop/rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease @design.options \
	-parser:protocol 3UBQ_design_control.xml -out:suffix _control -scorefile 3UBQ_control.fasc -s 3UBQ_relax.pdb > design_control.log &
~/rosetta_workshop/rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease @design.options \
        -parser:protocol 4HKX_design_control.xml -out:suffix _control -scorefile 4HKX_control.fasc -s 4HKX_relax.pdb > design_control.log &
