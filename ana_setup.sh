#!/bin/sh


export ANA_BASE_DIR=$PWD
export ANA_DATA_DIR=${ANA_BASE_DIR}/data
export ANA_SRC_DIR=${ANA_BASE_DIR}/src
export ANA_BIN_DIR=${ANA_BASE_DIR}/bin
export ANA_RUN_DIR=${ANA_BASE_DIR}/run
export ANA_LIST_DIR=${ANA_BASE_DIR}/list
export ANA_HIST_DIR=${ANA_BASE_DIR}/hist
export ANA_PLOT_DIR=${ANA_BASE_DIR}/plot
export ANA_PS_DIR=${ANA_BASE_DIR}/ps

if [ ! -d $ANA_DATA_DIR ]; then
  echo "Directory $ANA_DATA_DIR does not exist! Creating..."
  mkdir data
fi
if [ ! -d $ANA_BIN_DIR ]; then
  echo "Directory $ANA_BIN_DIR does not exist! Creating..."
  mkdir bin
fi
if [ ! -d $ANA_RUN_DIR ]; then
  echo "Directory $ANA_RUN_DIR does not exist! Creating..."
  mkdir run
fi
if [ ! -d $ANA_LIST_DIR ]; then
  echo "Directory $ANA_LIST_DIR does not exist! Creating..."
  mkdir list
fi
if [ ! -d $ANA_HIST_DIR ]; then
  echo "Directory $ANA_HIST_DIR does not exist! Creating..."
  mkdir hist
fi
if [ ! -d $ANA_PLOT_DIR ]; then
  echo "Directory $ANA_PLOT_DIR does not exist! Creating..."
  mkdir plot
fi
if [ ! -d $ANA_PS_DIR ]; then
  echo "Directory $ANA_PS_DIR does not exist! Creating..."
  mkdir ps
fi

printenv | grep ANA_

