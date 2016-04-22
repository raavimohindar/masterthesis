#! /bin/sh

latex slides
dvips -Z -V -t a4 -o slides.ps slides


pstops "12:8U@.25(20cm,08.0cm)+4U@.25(13cm,08.0cm)+0U@.25(6cm,08.0cm)+9U@.25(20cm,15.0cm)+5U@.25(13cm,15.0cm)+1U@.25(6cm,15.0cm)+10U@.25(20cm,22.0cm)+6U@.25(13cm,22.0cm)+2U@.25(6cm,22.0cm)+11U@.25(20cm,29.0cm)+7U@.25(13cm,29.0cm)+3U@.25(6cm,29.0cm)" slides.ps > slides12up.ps

pstops "4:2U@.5(20.1cm,15.7cm)+0U@.5(10.5cm,15.7cm)+1U@.5(10.5cm,29.5cm)+3U@.5(20.1cm,29.5.0cm)" slides.ps > slides4up.ps

pstops "2:0U@.5(10.5cm,15.7cm)+1U@.5(20.1cm,15.7cm)" slides.ps > slides2up.ps

rm -f *~ *.bck *.bak *.log *.aux *.dvi *.toc *.bbl *.blg *.lof *.ind *.ilg *.idx

ps2pdf slides.ps slides.pdf
