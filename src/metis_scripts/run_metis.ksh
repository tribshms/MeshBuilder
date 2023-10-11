#!/usr/bin/ksh

################################################################################
# Script used to create a *.reach file starting from a connectivity.meshb file
# using METIS. The *.reach file is needed for parallel runs in tRIBS.
#
# Usage ksh ./run_metis.ksh NumberOfNodes OPT_Part basename
#
# where
# 
# NumberOfNodes is the number of nodes that will be used to run the 
#                 parallel simulation in tRIBS. 
#
# OPT_Part      is the option used to select the kind of partitioning. 
#               It can be: 
#		1 --> SurfaceFlow (SF)
#		2 --> Surface-Subsurface Flow (SSF)
#		3 --> Surface-Subsurface Flow with Headwaters (SSF-H).
#
# basename        is the the name that indicates your simulation in tRIBS.
#
#
# Example: ksh ./run_metis.ksh 4 1 Baron_Fork
################################################################################

 
if [[ $# -ne 3 ]];then
	print "\n Error! \n"
   	print "After ksh ./run_metis.ksh, you must indicate 3 options:"
	print " 1st: number of nodes"
	print " 2nd: the option used to select the kind of partitioning"
	print " 3rd: a basename for tRIBS simulations"
	print "\n"
	print "Example: ksh ./run_metis.ksh 4 1 Baron_Fork\n"
   exit
fi

if [[ $2 -ge 4 ]] || [[ $2 -le 0 ]];then
	print "\n Error! \n"
   	print "Option 2 specifies the kind of partitioning and must be equal "
	print "to either 1 (SurfaceFlow), 2 (Surface-Subsurface Flow), "
	print "or 3 (Surface-Subsurface Flow with Headwaters)!"
	print "\n"
	print "Example: ksh ./run_metis.ksh 4 1 Baron_Fork\n"
   exit
fi

if [[ $2 -eq 1 ]];then
	perl connectivity2metis.pl connectivity.meshb $3_metis_flow.dat 1 0 0 0 0 
	# The above instruction creates a file for SF partitioning called:  basename_metis_flow.dat 

	./kmetis $3_metis_flow.dat $1 > flow_$3nodes_$1_stats.out
	# The above instruction creates a file called:  basename_metis_flow.dat.part.NumberOfNodes 
	# The report of the instruction is outputted in flow_basename_NumberOfNodesnodes_stats.out

	perl metis2tribs.pl $3_metis_flow.dat.part.$1 connectivity.meshb  $3_flow_$1nodes.reach > flow_$3_$1nodes_stats.out
	# The above instruction creates a file called:  basename_flow_NumberOfNodesnodes.reach
	# This is the file that nwill be used for parallel tRIBS simulation with NumberOfNodes nodes.
	# The report of the instruction is outputted in flow_basename_NumberOfNodesnodes_stats.out
fi

if [[ $2 -eq 2 ]];then
	perl connectivity2metis.pl connectivity.meshb $3_metis_nconn.dat 1 0 0 0 2
	# The above instruction creates a file for SSF partitioning called:  basename_metis_nconn.dat

	./kmetis $3_metis_nconn.dat $1 > nconn_$3nodes_$1_stats.out
	# The above instruction creates a file called:  basename_metis_nconn.dat.part.NumberOfNodes 
	# The report of the instruction is outputted in nconn_basename_NumberOfNodesnodes_stats.out

	perl metis2tribs.pl $3_metis_nconn.dat.part.$1 connectivity.meshb  $3_nconn_$1nodes.reach > nconn_$3_$1nodes_stats.out
	# The above instruction creates a file called:  basename_nconn_NumberOfNodesnodes.reach
	# This is the file that nwill be used for parallel tRIBS simulation with NumberOfNodes nodes.
	# The report of the instruction is outputted in nconn_basename_NumberOfNodesnodes_stats.out
fi

if [[ $2 -eq 3 ]];then
	perl connectivity2metis.pl connectivity.meshb $3_metis_upnconn.dat 1 1 0 0 2
	# The above instruction creates a file for SSF-H partitioning called:  basename_metis_upnconn.dat

	./kmetis $3_metis_upnconn.dat $1 > upnconn_$3nodes_$1_stats.out
	# The above instruction creates a file called:  basename_metis_upnconn.dat.part.NumberOfNodes 
	# The report of the instruction is outputted in upnconn_basename_NumberOfNodesnodes_stats.out

	perl metis2tribs.pl $3_metis_upnconn.dat.part.$1 connectivity.meshb  $3_upnconn_$1nodes.reach > upnconn_$3_$1nodes_stats.out
	# The above instruction creates a file called:  basename_upnconn_NumberOfNodesnodes.reach
	# This is the file that nwill be used for parallel tRIBS simulation with NumberOfNodes nodes.
	# The report of the instruction is outputted in upnconn_basename_NumberOfNodesnodes_stats.out
fi




#perl collectPartitionInfo.pl t_$3.connectivity.$1 connectivity.meshb >  flow_part_$3_$1nodes.txt
# The above instruction creates a file called:  flow_part_basename_NumberOfNodesnodes.txt
# containing some statistics about the partitioning about numbers of nodes, reaches, 
# downstream partitions, and flux partitions.


