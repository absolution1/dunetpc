function echo-info(){
    echo "INFO   : $@" >&2
}

function echo-error(){
    echo "ERROR  : $@" >&2
}

function usage(){
    echo-info "`basename $0` <runnum> <version> <qualifiers> [debug]"
    echo-info "`basename $0` 1234 v04_30_03 e9"
}

function parse_args(){
    if [ $# -lt 3 ];then
	usage
	exit 1
    fi
    RUN=$1
    VERSION=$2
    QUALIFIERS=$3

    if [ $# -eq 4 ];then
	DEBUG=1
    else
	DEBUG=0
    fi

    RELEASE_DIR=/home/lbnedaq/nearline/nearline_test_release_${VERSION}
    SCRIPT_PATH=${RELEASE_DIR}/srcs/dunetpc/dune/NearlineMonitor/scripts
    PRODUCTS_DIR=${RELEASE_DIR}/localProducts_larsoft_${VERSION}_${QUALIFIERS}_prof
    OUTPUT_PATH=/lbne/data2/users/lbnedaq/nearline_evd/${VERSION}

    echo-info "RUN $RUN"
    echo-info "RELEASE_DIR $RELEASE_DIR"
    echo-info "SCRIPT_PATH $SCRIPT_PATH"
    echo-info "PRODUCTS_DIR $PRODUCTS_DIR"
    echo-info "OUTPUT_PATH $OUTPUT_PATH"
}

function find_run_file(){
    RUN_ZEROS=`printf %06i $RUN`
    echo-info "RUN_ZEROS $RUN_ZEROS"
    TRANSFERRED_DIR=/data/lbnedaq/data/transferred_files
    NUM_FILES=`find ${TRANSFERRED_DIR}/lbne_r${RUN_ZEROS}_sr??_*.root 2>/dev/null | wc -l`

    if [ $NUM_FILES -eq 0 ];then
	echo-error "Failed to find RUN $RUN_ZEROS"
	exit 1
    elif [ $NUM_FILES -gt 1 ];then
	echo-error "Found more than one file for RUN $RUN_ZEROS : `find ${TRANSFERRED_DIR}/lbne_r${RUN_ZEROS}_sr??_*.root 2>/dev/null`"
	exit 1
    fi
    FILE_FULL_PATH=`find ${TRANSFERRED_DIR}/lbne_r${RUN_ZEROS}_sr??_*.root 2>/dev/null`
    echo-info "FILE_FULL_PATH $FILE_FULL_PATH"
}

function process_file(){
    FILE=`basename $FILE_FULL_PATH`
    bigrun=${FILE:6:3}
    run=${FILE:6:6}
    subrun=${FILE:15:2}
    RunDir=${OUTPUT_PATH}/$bigrun/$run
    mkdir -p $RunDir
    OLD_PWD=$PWD
    cd $SCRIPT_PATH
    echo-info "./ProcessSingleFileEVD.sh $RunDir $FILE_FULL_PATH $PRODUCTS_DIR"
    if [ $DEBUG -eq 1 ];then
	./ProcessSingleFileEVD.sh $RunDir $FILE_FULL_PATH $PRODUCTS_DIR
    else
	nohup ./ProcessSingleFileEVD.sh $RunDir $FILE_FULL_PATH $PRODUCTS_DIR >> /dev/null 2>&1 &
    fi
    cd $OLD_PWD
}

parse_args "$@"
find_run_file
process_file

