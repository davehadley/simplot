
startdir=$(pwd);

if [ -z $SIMPLOTROOT ]; then
  scriptdir=$(dirname $(readlink -f ${BASH_SOURCE[0]:-$0}));
  simplotdir=$(readlink -f $scriptdir);
  export SIMPLOTROOT="${simplotdir}"
fi

SIMPLOT_SRC_DIR=${SIMPLOTROOT}
SIMPLOT_BIN_DIR=${SIMPLOTROOT}/bin
export PYTHONPATH=${SIMPLOT_SRC_DIR}:${PYTHONPATH}
export PATH=${SIMPLOT_BIN_DIR}:${PATH}
export LD_LIBRARY_PATH=${SIMPLOT_SRC_DIR}/simplot/rootprob3pp/Prob3++.20121225:${LD_LIBRARY_PATH}


