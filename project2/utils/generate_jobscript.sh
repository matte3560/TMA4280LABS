#!/bin/bash

# Check arguments
if [ $# -lt 2 ]; then
	printf "Usage:\n"
	printf "  $0 jobname command [processes] [threads]\n\n"
	printf "Arguments:\n"
	printf "  jobname:    the name of the job\n"
	printf "  command:    the command to be executed by the job script\n"
	printf "  processes:  optional, number of proccesses to use for the job\n"
	printf "  threads:    optional, number of threads per process\n"
	exit 0
fi


# Initialize vars used to generate jobscript
NAME=$1
FILE="${NAME}_job.sh"
COMMAND=$2
NUM_PROCS=0
NUM_THR=0
if [ $# -ge 3 ]; then
	NUM_PROCS=$3
fi
if [ $# -ge 4 ]; then
	NUM_THR=$4
fi

# Output to jobscript
printf "Writing jobscript to ${FILE}\n"
(
	printf "#!/bin/bash\n\n"
	printf "#PBS -N ${NAME}\n"
	printf "#PBS -q training\n"
	printf "#PBS -W group_list=imf_lille-tma4280\n"
	printf "#PBS -l walltime=00:05:00\n"

	# Determine number of nodes
	NODES=1
	if [ $NUM_PROCS -gt 20 ]; then
		NODES=2
	else
		if [ $(expr $NUM_PROCS \* $NUM_THR) -gt 20 ]; then
			NODES=2
		fi
	fi

	# Request resources
	printf "#PBS -l select=${NODES}:ncpus=20"
	if [ $NUM_PROCS -ne 0 ]; then
		if [ $NODES -eq 1 ]; then
			printf ":mpiprocs=${NUM_PROCS}"
		else
			printf ":mpiprocs=$(expr ${NUM_PROCS} / 2)"
		fi
	fi
	# Append number of threads if provided and end the line
	if [ $NUM_THR -ne 0 ]; then
		printf ":ompthreads=${NUM_THR}\n\n"
	else
		printf "\n\n"
	fi

	# Change into dir job was submitted from and load modules
	printf 'cd $PBS_O_WORKDIR'
	printf "\n"
	printf "module load gcc\n"
	printf "module load openmpi\n"
	printf "module load openblas\n\n"

	# Insert the actual command
	if [ $NUM_PROCS -eq 0 ]; then
		printf "exec ${COMMAND}\n"
	else
		printf "mpirun ${COMMAND}\n"
	fi
) > $FILE
# Make file executable
chmod +x $FILE
