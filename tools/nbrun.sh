#!/bin/bash

function usage () {
    echo "USAGE: ${0} <path/to/notebook.ipynb> [-batch] [some variables, used in notebook in format: variable1=value1 variable2=value2]"
    exit 1
}

test -z ${1} && usage

SCRIPT_PATH=${1}

shift

if [ "$1" == "-batch" ]; then
    BATCH=yes
    shift
else
    BATCH=no
fi


cd $(dirname ${SCRIPT_PATH})

TMPFILE=`mktemp .nbrunXXXXXXXX.ipynb`

cp -f $(basename ${SCRIPT_PATH}) $TMPFILE

for i in "$@"; do
    VAR=$(echo $i | cut -d'=' -f1)
    VAL=$(echo $i | cut -d'=' -f2)
    # echo "var:" $VAR and "val:" $VAL
    sed -i "s|\"${VAR}[[:blank:]]*=.*\\\\n\"|\"${VAR}=${VAL}\\\\n\"|g" $TMPFILE
done

jupyter nbconvert ${TMPFILE} --stdout --to python | sed '/get_ipython/d' > ${TMPFILE}.py

echo  >> ${TMPFILE}.py

if [ "$BATCH" == "no" ]; then
    echo 'input("Please press RETURN to exit ")' >> ${TMPFILE}.py
fi

python ${TMPFILE}.py

rm -f $TMPFILE ${TMPFILE}.py
