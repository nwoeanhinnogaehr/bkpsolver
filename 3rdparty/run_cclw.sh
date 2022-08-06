#!/bin/bash
python 3rdparty/Run_Bilevel_KP.py > out
sed -i 1d out
cat out
rm out
