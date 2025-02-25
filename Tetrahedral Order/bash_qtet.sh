#!/bin/bash

# defining path and file names
tpr="prod.tpr"
gro="prod.gro"
xtc="traj.xtc"
exe="gmx"
r="0.5"
out="qtet_water_ave.dat"
frame="50000"

# loop over number of frames
for i in $(seq 1 $frame) ;do
	time=$((i * 10)) 
	$exe trjconv -s $tpr -f $xtc -o frame.gro -b $time -e $time <<EOF
0
EOF
	$exe select -f frame.gro -s frame.gro -select 'name OW and within 0.5 of group protein' -on temp1.ndx
	$exe select -f frame.gro -s frame.gro -select 'name C O N and group protein' -on temp2.ndx
	cat temp1.ndx temp2.ndx > temp3.ndx
	$exe make_ndx -f frame.gro -n temp3.ndx -o final_index.ndx << EOF
0 | 1
name 2 heavy
q
EOF
	rm temp*
	$exe trjconv -f frame.gro -s frame.gro  -n final_index.ndx -o heavy.gro << EOF
	
2
EOF
	num_mol=$(sed -n '2p' "heavy.gro")
	echo $num_mol
	sed -i "s/XXX/${num_mol}/g" qtet.f90	
	gfortran qtet.f90
	./a.out
	qtet_value=$(cat fort.100)	
	echo " $i: $qtet_value" >> $out
	rm heavy.gro
	rm final_index.ndx
	echo $i
	rm frame.gro
	sed -i "s/${num_mol}/XXX/g" qtet.f90
	rm fort.100	
done


