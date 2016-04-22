#!/bin/sh
# farbiger Prompt
#################
export PS1="\[\033[1;7;34m\]\u@\h\[\033[27m\] [\w]\[\033[39;22;0m\]: "
alias kill="kill -9"
alias wo="ps -aux|grep mildner|grep mozilla|grep usr"
alias rm="rm -i"
alias xload="xload -hl red -fg green -bg black"

if [ $OSTYPE == linux-gnu ]; then
  alias ls="ls --color=never -h"
  alias du="du -h"
fi
GARNOME=/ant/huge/users/student/raavi02/gnome
     PATH=$GARNOME/bin:$GARNOME/sbin:$PATH
     LD_LIBRARY_PATH=$GARNOME/lib:$LD_LIBRARY_PATH
     PYTHONPATH=$GARNOME/lib/python2.3/site-packages
     PKG_CONFIG_PATH=$GARNOME/lib/pkgconfig:/usr/lib/pkgconfig
     XDG_DATA_DIRS=$GARNOME/share
     XDG_CONFIG_DIRS=$GARNOME/etc/xdg
     GDK_USE_XFT=1
     export PATH LD_LIBRARY_PATH PYTHONPATH PKG_CONFIG_PATH GDK_USE_XFT XDG_DATA_DIRS XDG_CONFIG_DIRS
