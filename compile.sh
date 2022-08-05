#! /bin/sh
cd cUtil/samtools
make libbam.so
cd ../
#cd cUtil/
python setup.py build_ext --inplace
cd ../solve/
python setup.py build_ext --inplace

