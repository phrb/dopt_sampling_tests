#! /bin/bash

echo "Installing R packages"

Rscript -e 'install.packages(c("AlgDesign", "rsm", "uuid", "dplyr", repos="https://cran.rstudio.com")'

CLONE_TARGET="/root/dopt_sampling_tests"

echo "Updating target data directory"

if [ -d "$CLONE_TARGET" ]; then
    git -C ${CLONE_TARGET} pull
else
    git clone --depth 1 https://github.com/phrb/dopt_sampling_tests.git
fi

USR="pbruel"
USR_TARGET="/home/${USR}/dopt_sampling_tests/data/results/"
NODE_NAME="xeon_e5_2630_v3_$(uname -n | cut -d. -f1)"

APP_TARGET="/root/dopt_sampling_tests/src/"
cd $APP_TARGET

Rscript run_experiments.r

NODE_NAME="xeon_e5_2630_v3_$(uname -n | cut -d. -f1)_$(date +%s)"

mkdir -p $NODE_NAME
mv ${APP_TARGET}/${NODE_NAME}_* /tmp/

su ${USR} -c "mkdir -p ${USR_TARGET}"
su ${USR} -c "mv /tmp/${NODE_NAME}_* ${USR_TARGET}"