#!/bin/bash

ME=$(realpath $0)
MYDIR=$(dirname $ME)
TMPDIR="$(mktemp -d /tmp/.picopic.XXXX)"

cd $MYDIR && cd ../../
make distclean 2>&1 >/dev/null
cp -rp . ${TMPDIR}

cd ${TMPDIR}
git reset --hard 2>&1 >/dev/null

SUM=""

for var in `seq 5`; do
    make distclean 2>&1 >/dev/null
    ./autogen.sh 2>&1 >/dev/null
    ./configure 2>&1 >/dev/null
    make test 2>&1 >/dev/null
    SUM+=$(make test 2>&1 | grep 'Execution time is' | cut -d' ' -f4 | head -n1)
    if [ ! $var == 5 ]; then
	SUM+="+"
    fi

    git checkout HEAD~1 2>&1 >/dev/null
done

MIDDLE=$(echo \(${SUM}\)/$(echo ${SUM} | sed 's/\+/\ /g' | wc -w) | bc)

FIRST=$(echo ${SUM} | cut -d'+' -f1)

if [ "$(echo ${FIRST} > ${MIDDLE} | bc)" == 0 ]; then
    echo FAIL
else
    echo OK
fi

rm -rf ${TMPDIR}
