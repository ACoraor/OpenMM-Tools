#!/bin/bash

COUNT=0
f=force_0pN
#for f in force_*/; do
	cd $f
	for m in */; do
		cd $m
		#Now in each ptm
		for r in {0..2}; do
			cd $r/processing
			echo "On ${f}/${m}/${r}"
			#echo `pwd`
			#sleep 20
			for n in {0..35}; do
				if [ -f "../files/log.lammps.${n}.prev" ]; then
					#Previous files exist. Check if they're up-to-date.
					var1=`tail -n 1 ../files/log.lammps.${n}`
					var2=`tail -n 1 ../files/log.lammps.${n}.prev`
				
					if [ "$var1" = "$var2" ]; then
						#Prev includes current run.
						cp ../files/log.lammps.${n}.prev log.lammps.${n}.all
					else
						#Prev doesn't include current run.
						cat ../files/log.lammps.${n}.prev ../files/log.lammps.${n} > log.lammps.${n}.all
					fi
				else
					cat ../files/log.lammps.${n}.prev ../files/log.lammps.${n} > log.lammps.${n}.all
				fi
			done

			if [ -f "../files/log.lammps.prev" ]; then
				#Assume prev doesn't include current run
				cat ../files/log.lammps.prev ../files/log.lammps > log.lammps
			fi
			python ../../../../make_log.py 

			let COUNT=COUNT+1
		    if [ $((${COUNT} % 20)) -eq 0 ]; then
		        echo "Waiting at ${COUNT}"
		        wait
		    fi
			cd ../..
		done
		cd ..
	done
	cd ..
#done
