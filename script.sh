#!/bin/bash

#run_time=3
output_dir="y0_y25_s3_ec_ics_fbetahash_5min_5_iter"
mkdir -p $output_dir
for run_time in 5 
do
	
	#mkdir -p "$run_time"min
	for i in {1..1}
	do
		for j in {0..4}
		do
			for graph_dens in syeast0 syeast0_complement syeast25 syeast25_complement
			do
				for network_type in "$graph_dens"_node_del_100_add_edges_0.3 
				do					
					base_network="$network_type"_iter_"$j" 
					
					for network in "$graph_dens"
					do		
						for measure in ec ics fbetahash
						do
							result_path="$PWD"/"$output_dir"/"$run_time"min/"$measure"/"$base_network"_"$network"_"$i"
							mkdir -p "$result_path"		
							./sana -g1 "$base_network" -g2 "$network" -t $run_time -s3 0 -"$measure" 1 -truealignment true_alignments/"$base_network".align
							mv -v "$PWD"/sana.out "$result_path"/sana_"$measure".out
							mv -v "$PWD"/sana.align "$result_path"/sana_"$measure".align
						done
		
						result_path="$PWD"/"$output_dir"/"$run_time"min/s3/"$base_network"_"$network"_"$i"
						mkdir -p "$result_path"		
						./sana -g1 "$base_network" -g2 "$network" -t $run_time -truealignment true_alignments/"$base_network".align
						mv -v "$PWD"/sana.out "$result_path"/sana_s3.out
						mv -v "$PWD"/sana.align "$result_path"/sana_s3.align
						
						#for beta in 0.1 0.33 1 3 10
						#do
						#	result_path="$PWD"/"$output_dir"/"$run_time"min/"F$beta"/"$base_network"_"$network"_"$i"
						#	mkdir -p "$result_path"	
						#	./sana -g1 "$base_network" -g2 "$network" -t $run_time -s3 0 -fbeta 1 -b "$beta" -truealignment true_alignments/"$base_network".align
						#	mv -v "$PWD"/sana.out "$result_path"/"sana_f$beta.out"
						#	mv -v "$PWD"/sana.align "$result_path"/"sana_f$beta.align"
						#						
						#done
					done				
				done  
			done
		done	
	done	
done

 
