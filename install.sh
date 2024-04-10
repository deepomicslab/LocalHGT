#!/bin/bash

make



## below is to avoid the bug: No module named _sysconfigdata_x86_64_conda_cos7_linux_gnu
if [ -f "${CONDA_PREFIX}/lib/python3.7/_sysconfigdata_x86_64_conda_cos6_linux_gnu.py" ];then

cp ${CONDA_PREFIX}/lib/python3.7/_sysconfigdata_x86_64_conda_cos6_linux_gnu.py ${CONDA_PREFIX}/lib/python3.7/_sysconfigdata_x86_64_conda_cos7_linux_gnu.py

fi

python setup.py install
