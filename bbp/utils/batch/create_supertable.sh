#!/bin/bash

for method in gp sdsu ucsb exsim song irikura_recipe_m1 irikura_recipe_m2; do
    cmd="combine_bbp_data.py ${method}-results-${1}.txt"
    for event in ch ar whittier nps tottori niigata nr lomap landers chuetsu parkfield sansimeon iwate; do
	cmd="${cmd} ${2}/${method}-${event}-50r"
    done;
    $cmd
done;

for method in gp sdsu ucsb exsim song irikura_recipe_m1 irikura_recipe_m2; do
    report_bbp_data.py ${method}-results-${1}.txt > ${method}-part-a-table-${1}.txt
done;
