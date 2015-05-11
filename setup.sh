
startdir=$(pwd);

if [ -z $SIMPLOTROOT ]; then
  scriptdir=$(dirname $(readlink -f $0));
  simplotdir=$(readlink -f $scriptdir);
  export SIMPLOTROOT="${simplotdir}"
fi

SIMPLOT_SRC_DIR=${SIMPLOTROOT}
export PYTHONPATH=${SIMPLOT_SRC_DIR}:${PYTHONPATH}
#echo $SIMPLOTROOT

