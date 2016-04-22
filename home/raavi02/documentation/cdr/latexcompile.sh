#! /bin/sh


set -x
set -e
a=c1_mu
echo " compiling $a"

latex --interaction=  $a.tex

dvips -o $a.ps $a.dvi

ps2epsi $a.ps $a.eps

rm *.aux *.dvi *.log *.ps


