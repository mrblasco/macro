#!/bin/bash

base_dir=`pwd`

submissions_num=`find ./ -maxdepth 1 -name "sub*" | wc -l`

src_dir="src"
sub_bin_dir="bin"
if [ ! -d $sub_bin_dir ]; then mkdir $sub_bin_dir; fi

if [ "$1" = "clean" ]; then
	echo "*** Cleaning ***"
	if [ -d $sub_bin_dir ]; then rm -rf $sub_bin_dir; fi
	rm submission*/run_simu.m
	exit
fi

sub_start=1
sub_end=$submissions_num
if [ "$1" != "" ]; then
	sub_start=$1
	sub_end=$1
fi

for i in `seq -f '%03g' $sub_start $sub_end`; do
	echo "*** $i ***"
	sub_dir="submission_$i"
	sub_lang=`tail -n 2 $sub_dir/competitor.info | head -n 1`
	sub_handle=`tail -n 6 $sub_dir/competitor.info | head -n 1`
	sub_bin_name="submission_$i"
	if [ "$sub_lang" = "C++" ] || [ "$sub_lang" = "C#" ] || [ "$sub_lang" = "Python" ]; then
		echo "C++"
		build_dir="build_cpp"
		if [ -d $build_dir ]; then rm -rf $build_dir; fi
		mkdir $build_dir
		simu_file="$build_dir/run_simu.cpp"
		cat $sub_dir/submission_$i.cpp > $simu_file
		echo "" >> $simu_file
		cat $src_dir/run_simu.cpp >> $simu_file
		cd $build_dir
		g++ -O2 -std=c++11 run_simu.cpp -o "$sub_bin_name.run"
		cd $base_dir
		mv "$build_dir/$sub_bin_name.run" "$sub_bin_dir/$sub_bin_name.run"
		rm -rf $build_dir
	elif [ "$sub_lang" = "Java" ]; then
		echo "Java"
		build_dir="build_java"
		rm -rf $build_dir; mkdir $build_dir
		cp $src_dir/NationSaveTester.java $build_dir/
		cp $sub_dir/submission_$i.java $build_dir/NationSave.java
		cd $build_dir
		javac NationSaveTester.java
		jar cfe "$sub_bin_name.jar" NationSaveTester *.class
		chmod +x "$sub_bin_name.jar"
		cd $base_dir
		mv "$build_dir/$sub_bin_name.jar" "$sub_bin_dir/$sub_bin_name.jar"
		rm -rf $build_dir
	elif [ "$sub_lang" = "MATLAB" ]; then
		echo "MATLAB"
		cp -f "$src_dir/run_simu.m" "$sub_dir/run_simu.m"
	else
		echo "Unsupported"
	fi
done
