
print_general_slurm_env_vars() {
    # Display allocated resources
    echo "JOB INFO"
    echo -e "\t Job name: $SLURM_JOB_NAME"
    echo -e "\t Job ID: $SLURM_JOB_ID"
    echo -e "\t Job node list: $SLURM_JOB_NODELIST"
    echo -e "\t Node list: $SLURM_NODELIST"
    echo -e "\t Submit dir: $SLURM_SUBMIT_DIR"

    echo -e "\nRESOURCES"
    echo -e "\tNumber of ntasks: $SLURM_NTASKS"
    echo -e "\tCPUs per task: $SLURM_CPUS_PER_TASK"
    echo -e "\tTotal memory allocated to this job: $SLURM_MEM_PER_NODE"
    echo -e "\tTotal memory per CPU: $SLURM_MEM_PER_CPU\n\n"

}

print_job_array_stuff() {
    # Echo job array stuff
    echo -e "\tARRAY INFO"
    echo -e "\tjob ID: $SLURM_ARRAY_JOB_ID" ## will be set to the first job ID of the array.
    echo -e "\ttask ID: $SLURM_ARRAY_TASK_ID" ## will be set to the job array index value.
    echo -e "\ttask count: $SLURM_ARRAY_TASK_COUNT" ## will be set to the number of tasks in the job array.
    echo -e "\ttask max: $SLURM_ARRAY_TASK_MAX" ## will be set to the highest job array index value.
    echo -e "\ttask min: $SLURM_ARRAY_TASK_MIN\n\n" ## will be set to the lowest job array index value.
}

is_slurm_array_job() {
    if [[ -n "$SLURM_ARRAY_TASK_ID" ]]; then
        echo -e "Running as a job array (task ID: $SLURM_ARRAY_TASK_ID)\n\n"
        return 0
    else
        echo -e "Not running as a job array\n\n"
        return 1
    fi
}

redirect_output() {
    # if is_slurm_array_job; then
    # OUTDIR="logs"
    OUTDIR=$1
    LOGFILE="${OUTDIR}/${SLURM_JOB_ID}.log"
    ERRFILE="${OUTDIR}/${SLURM_JOB_ID}.err"
    exec > "$LOGFILE"
    exec 2> "$ERRFILE"
    # fi
}

remove_empty_file() {
    # Check if file is empty
    if [ ! -s "$1" ]; then
        echo "$1 is empty. Deleting..."
        rm "$1"  # Delete the empty file
    else
        echo "$1 is not empty."
    fi
}


move_db() {
    cp $1 $2
}