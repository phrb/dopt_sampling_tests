#! /bin/bash

echo "Deploying orio_experiments image"

kadeploy3 -f ${OAR_NODE_FILE} -e debian9-x64-big -k

echo "Waiting for deply to finish"
sleep 10

echo "Launching jobs"

JOB_SCRIPT="/home/pbruel/dopt_sampling_tests/src/run.sh"
taktuk -c "ssh" -l root -f <( uniq ${OAR_FILE_NODES} ) broadcast exec [ ${JOB_SCRIPT} ]
