#!/usr/bin/env bash
dirname=$(dirname $0)
set -ue
echo -e "\tLavReader/R19S5-trunc.lav:"
PYTHONPATH=.. $dirname/unit.test.py LavReader $dirname/R19S5-trunc.lav | diff -s - $dirname/R19S5-trunc.lav.out
echo -e "\tLavReader.convert/R34S1-part.lav:"
PYTHONPATH=.. $dirname/unit.test.py LavReader.convert $dirname/R34S1-part.lav | diff -s - $dirname/R34S1-part.lav.convert.out
