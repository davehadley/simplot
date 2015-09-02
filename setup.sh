
startdir=$(pwd);

if [ -z $SIMPLOTROOT ]; then
  scriptdir=$(dirname $(readlink -f $0));
  simplotdir=$(readlink -f $scriptdir);
  export SIMPLOTROOT="${simplotdir}"
fi

SIMPLOT_SRC_DIR=${SIMPLOTROOT}
SIMPLOT_BIN_DIR=${SIMPLOTROOT}/bin
export PYTHONPATH=${SIMPLOT_SRC_DIR}:${PYTHONPATH}
export PATH=${SIMPLOT_BIN_DIR}:${PATH}
#echo $SIMPLOTROOT

