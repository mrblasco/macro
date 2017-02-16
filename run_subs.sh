#!/bin/bash

base_dir=`pwd`

data_dir="data"
file_eco="$base_dir/$data_dir/economies_final.csv"
file_shocks="$base_dir/$data_dir/shocks_system.csv"

submissions_num=`find ./ -maxdepth 1 -name "sub*" | wc -l`

sub_bin_dir="bin"
sub_out_dir="$base_dir/output"
if [ ! -d $sub_out_dir ]; then mkdir $sub_out_dir; fi

if [ "$1" = "clean" ]; then
	echo "*** Cleaning ***"
	if [ -d $sub_out_dir ]; then rm -rf $sub_out_dir; fi
	exit
fi

sub_list=(`seq 1 $submissions_num`)
if [ "$1" != "" ]; then
	sub_list=($1)
fi

time_tot0=`date +%s%N | cut -b1-13`
for ii in `seq 0 $(( ${#sub_list[@]} - 1 ))`; do
	i=`printf '%03g' ${sub_list[$ii]}`
	echo "*** $i ***"
	
	# Submission info
	sub_dir="submission_$i"
	sub_lang=`tail -n 2 $sub_dir/competitor.info | head -n 1`
	sub_handle=`tail -n 6 $sub_dir/competitor.info | head -n 1`
	sub_bin_name="submission_$i"
	sub_out_name="$sub_out_dir/submission_$i.csv"
	
	if [ "$sub_lang" = "C++" ] || [ "$sub_lang" = "C#" ] || [ "$sub_lang" = "Python" ]; then
		echo "C++"
		cd $sub_bin_dir
		./$sub_bin_name.run $file_eco $file_shocks "$sub_out_name"
		cd $base_dir
	
	elif [ "$sub_lang" = "Java" ]; then
		echo "Java"
		cd $sub_bin_dir
		java -jar "$sub_bin_name.jar" -eco $file_eco -shocks $file_shocks -out "$sub_out_name"
		cd $base_dir
	
	elif [ "$sub_lang" = "MATLAB" ]; then
		echo "MATLAB"
		cd "$base_dir/$sub_dir"
		matlab -nojvm -r "run_simu('$file_eco', '$file_shocks', '$sub_out_name'); quit"
		cd $base_dir

	else
		echo "Unsupported"
	fi
done
time_tot1=`date +%s%N | cut -b1-13`
time_tot_delta=$(( $time_tot1 - $time_tot0 ))
echo "Elapsed time: $time_tot_delta ms"

