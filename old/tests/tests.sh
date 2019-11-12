#!/usr/bin/env bash
dirname=$(dirname $0)
set -u
echo -e "\tLavReader ::: R19S5-trunc.lav:"
PYTHONPATH=$dirname/.. $dirname/unit.test.py LavReader -l $dirname/R19S5-trunc.lav | diff -s - $dirname/R19S5-trunc.lav.out
#TODO: check new output after fixing conversion bug in 78aeb53
# echo -e "\tLavReader.convert ::: R34S1-part.lav:"
# PYTHONPATH=$dirname/.. $dirname/unit.test.py LavReader.convert -l $dirname/R34S1-part.lav | diff -s - $dirname/R34S1-part.lav.convert.out
echo -e "\tFastaLineGenerator.extract ::: R26S11.fa:"
PYTHONPATH=$dirname/.. $dirname/unit.test.py FastaLineGenerator.extract -C $dirname/R26S11.fa.extract.in -f $dirname/R26S11.fa | diff -s - $dirname/R26S11.fa.extract.out
