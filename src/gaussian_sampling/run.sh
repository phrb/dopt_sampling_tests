#! /bin/bash

echo "Installing Julia"

wget https://julialang-s3.julialang.org/bin/linux/x64/1.1/julia-1.1.0-linux-x86_64.tar.gz
tar xvf julia-1.1.0-linux-x86_64.tar.gz

JULIA_PROGRAM="julia-1.1.0/bin/julia"

echo "Installing Julia Packages"

$JULIA_PROGRAM -e 'using Pkg; Pkg.add(["Distributions", "DataFrames", "CSV", "RCall"]);'

echo "Installing R packages"

Rscript -e 'install.packages(c("AlgDesign", "rsm", "uuid", "dplyr"), repos="https://cran.rstudio.com")'

CLONE_TARGET="/root/dopt_sampling_tests"

echo "Updating target data directory"

if [ -d "$CLONE_TARGET" ]; then
    git -C ${CLONE_TARGET} pull
else
    git clone --depth 1 https://github.com/phrb/dopt_sampling_tests.git
fi

USR="pbruel"
USR_TARGET="/home/${USR}/dopt_sampling_tests/data/gaussian_sampling/"

APP_TARGET="/root/dopt_sampling_tests/src/gaussian_sampling/"
cd $APP_TARGET

OUTPUT_FILE="results.csv"

echo "Starting script run"

Rscript run_experiments.r

echo "Copying files to main machine"

NODE_NAME="$(uname -n | cut -d. -f1)_$(date +%s)"

mkdir -p $NODE_NAME
mv $OUTPUT_FILE $NODE_NAME
mv ${APP_TARGET}/${NODE_NAME} /tmp/

su ${USR} -c "mkdir -p ${USR_TARGET}"
su ${USR} -c "cp -r /tmp/${NODE_NAME} ${USR_TARGET}"
