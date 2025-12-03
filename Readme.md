
1  nohup ./00_Pipeline_launch.sh 1 3-5 > nohup_01_3-5_pmat.out 2>&1 &
2  nohup ./00_Pipeline_launch.sh 1 --part2 6-10 > nohup_01_6-10_part2.out 2>&1 &
3  nohup ./00_Pipeline_launch.sh 2 3 > nohup_02_Sum_mut_loci_3.out 2>&1 &
4  nohup ./00_Pipeline_launch.sh 3 3-5 > nohup_03_Merging_CTGA_3-5.out 2>&1 &
5  nohup ./00_Pipeline_launch.sh 4 --merge 5 3-5 > nohup_04_Rider_processing_merge5_3-5.out 2>&1 &