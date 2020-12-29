#!/bin/env bash

ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Export TAX4PROC bin
export TAX4PROC=$ROOT
export TAX4PROC_BIN=$TAX4PROC/bin
export PATH=$TAX4PROC_BIN:$PATH

# Export data directories
export DATA_ROOT_DIR=$HOME/data
export BRTAX4_DATA_ROOT=$DATA_ROOT_DIR/brtax4
export BRTAX4_DATA_RAW=$DATA_ROOT_DIR/brtax4/raw
export BRTAX4_DATA_PASS5=$DATA_ROOT_DIR/brtax4/pass5/data
export BRTAX4_MC_ROOT=$DATA_ROOT_DIR/brtax4/mc
export BRTAX4_HYBRID_ROOT=$DATA_ROOT_DIR/hyb_brtax4

export MDTAX4_DATA_ROOT=$DATA_ROOT_DIR/mdtax4
export MDTAX4_DATA_RAW=$DATA_ROOT_DIR/mdtax4/raw
export MDTAX4_DATA_PASS5=$DATA_ROOT_DIR/mdtax4/pass5/data
export MDTAX4_MC_ROOT=$DATA_ROOT_DIR/mdtax4/mc
export MDTAX4_HYBRID_ROOT=$DATA_ROOT_DIR/hyb_mdtax4

export SDTAX4_DATA_ROOT=$DATA_ROOT_DIR/sdtax4/
