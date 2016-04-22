#! /bin/sh


set -x
set -e
a=get
echo " compiling $a"

latex --interaction=  $a.tex

dvips -o $a.ps $a.dvi

ps2epsi $a.ps $a.eps

rm *.aux *.dvi *.log *.ps


