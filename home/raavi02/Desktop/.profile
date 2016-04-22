
#echo ".profile"
export DOT_PROFILE_READ=1



export OSTYPE





unset PATH
unset LIBRARY_PATH
unset LD_LIBRARY_PATH
unset LD_RUN_PATH
unset MANPATH
unset INFOPATH
unset INCLUDE_PATH
unset C_INCLUDE_PATH
unset CPLUS_INCLUDE_PATH

# Einstellungen, damit das System ueberhaupt funktioniert
export PATH=/bin:/sbin
export LD_LIBRARY_PATH=/lib
export PGSQLDIR=/usr
export RINASDIR=/home/rinas/Research/soft

# export KDEDIR=/usr/local/kde-3.1.1__deleteme_on_29.09.2003
# export KDEDIR=/usr/local/kde-3.1.5
# export KDEDIR=/usr/local/kde
# export KDEHOME=~/.kdetest
# export KDETMP=/var/tmp
# export QTDIR=$KDEDIR
# /home/rinas/Research/Algorithmen/DRM/soft
# $KDEDIR
 
BASEDIRS="
          $RINASDIR/$OSTYPE
          $RINASDIR

          /home/rinas/Research/soft_add
          /home/rinas/Research/soft_add/$OSTYPE
 
	  /home/rinas/Allgemeines
	  
	  /home/rinas/Buecher

	  /opt/sfw
	  /usr/ccs
	  
	  /usr/local/pgsql
          /usr/local
	  /usr
	  /usr/dt
	  /usr/X11R6" 


#          /usr/local/lib/gcc-lib/sparc-sun-solaris2.8/2.95.3
#          /usr/local/lib/gcc-lib/sparc-sun-solaris2.9/3.2


for BASEDIR in $BASEDIRS; do
  if [ -d ${BASEDIR} ]; then

    if [ -d ${BASEDIR}/bin ]; then
      export PATH="${PATH}:${BASEDIR}/bin"
    fi
    
    if [ -d ${BASEDIR}/sbin ]; then
      export PATH="${PATH}:${BASEDIR}/sbin"
    fi
    
    if [ -d ${BASEDIR}/lib ]; then
      # eine dusselige Ausnahme, damit matlab funktioniert!
      if [ ! ${BASEDIR} == "/usr/X11R6" ]; then
        export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${BASEDIR}/lib" # see man ld.so: the search path for libs at __run time__
        export LIBRARY_PATH="${LIBRARY_PATH}:${BASEDIR}/lib"       # see info gcc: the search path for libraries at __compile time__
        export LD_RUN_PATH="${LD_RUN_PATH}:${BASEDIR}/lib"         # see man ld: the search path path to libs at __link time__
      fi
    fi  

    if [ -d ${BASEDIR}/share/man ]; then
      export MANPATH="${MANPATH}:${BASEDIR}/share/man"
    fi  

    if [ -d ${BASEDIR}/man ]; then
      export MANPATH="${MANPATH}:${BASEDIR}/man"
    fi  
    
    if [ -d ${BASEDIR}/info ]; then
      export INFOPATH="${INFOPATH}:${BASEDIR}/info"
    fi  

    if [ -d ${BASEDIR}/include ]; then
      export INCLUDE_PATH="${INCLUDE_PATH}:${BASEDIR}/include"             # not used
      export C_INCLUDE_PATH="${C_INCLUDE_PATH}:${BASEDIR}/include"         # see info gcc: search path for include files     
      export CPLUS_INCLUDE_PATH="${CPLUS_INCLUDE_PATH}:${BASEDIR}/include" # see info gcc: search path for include files
    fi  
  fi
done

# ueberzaehlige Doppelpunkte entfernen
export MANPATH="${MANPATH:1}"
export INFOPATH="${INFOPATH:1}"
export INCLUDE_PATH="${INCLUDE_PATH:1}"


export PATH="${PATH}:."



export CXXFLAGS=-fpermissive


if [ -x /usr/bin/lesspipe ]; then
  eval $(lesspipe)
fi  

source ~/.bashrc
