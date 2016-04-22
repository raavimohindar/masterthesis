#! /bin/sh


set -x
set -e
a=top
echo " compiling $a"

latex --interaction=  $a.tex

dvips -o $a.ps $a.dvi

rm *.aux *.dvi *.log

