#########################################################

dnl AC_EXTRACT_REGEX(list, regEx)
dnl variable $regExResult returns the first entry in 'list' that matches the 
dnl regular expression 'regEx', using the 'expr' utility
dnl $regExResult = "" if nothing is found
AC_DEFUN(AC_EXTRACT_REGEX, 
[
    regExResult=""
    if test "$1" != ""; then
        for i in $1; do
            regExResult=`expr "$i" : "$2"`
            if test "$regExResult" != ""; then
                break
            fi
        done
    fi
])

#########################################################

dnl AC_FIND_PACKAGE(packageName, packageLib, packageInc, packageComment)
dnl defines with_packageName=yes/no
dnl         with_packageNamelib=<path>/empty if not found
dnl         with_packageNameinc=<path>/empty if not found
dnl example:
dnl     AC_FIND_PACKAGE(tiff, tiff, tiff.h, support import/export of tiff images)
AC_DEFUN(AC_FIND_PACKAGE, 
[
    AC_ARG_WITH([$1], [  --with-$1=dir : $4. 
                    if dir='yes': $1 package files will be searched for in 
                       some standard directories.
                    if dir is a directory: $1 package files will be searched for   
                       below dir using 'find'.
                    alternatively, you can specify:], ,)
    AC_ARG_WITH([$1lib], [    --with-$1lib=dir : the $1 package's lib directory], ,)
    AC_ARG_WITH([$1inc], [    --with-$1inc=dir : the $1 package's include directory
                         ], ,)

    if test "$libext" = ""; then
        if test "$CYGWIN" = "yes" || test "$MINGW32" = "yes"; then
          libext="lib"
        else
          libext="a"
        fi
        libext=a
    fi
    
    if test "$solibext" = ""; then
        if test "$CYGWIN" = "yes" || test "$MINGW32" = "yes"; then
          solibext="dll"
        else
          solibext="so"
        fi
    fi  
    if test "$CYGWIN" = "yes" || test "$MINGW32" = "yes"; then
      libpre=""
    else
      libpre="lib"
    fi

    
    if test ${with_[$1]:-""} = "" -a ${with_[$1]lib:-""} = "" -a ${with_[$1]inc:-""} = ""; then
        with_[$1]="no"
    fi
    
    if test ${with_[$1]:-""} != "no"; then
        AC_MSG_CHECKING([for lib$2 ])
        dirs=""
        if test "$with_[$1]lib" != ""; then
            dirs=$with_[$1]lib
        elif test "$with_[$1]" != "yes"; then
            dirs=$with_[$1]
        else
            dirs="/usr/local/lib /usr/local/gnu/lib /usr/local/[$1] /opt/lib /opt/gnu/lib /opt/[$1] /usr/lib"
        fi
        found=""
        for i in $dirs; do
            if test -d $i; then
                found="$found "`find $i -name "${libpre}[$2].$solibext" -print 2> /dev/null; find $i -name "${libpre}[$2].$libext" -print 2> /dev/null`
            fi
        done

        AC_EXTRACT_REGEX($found, \(.*lib\)/${libpre}$2\.$solibext)
        if test "$regExResult" = ""; then
            AC_EXTRACT_REGEX($found, \(.*lib\)/${libpre}$2\.$libext)
        fi
        if test "$regExResult" = ""; then
            AC_EXTRACT_REGEX($found, \(.*\)/${libpre}$2\.$solibext)
        fi
        if test "$regExResult" = ""; then
            AC_EXTRACT_REGEX($found, \(.*\)/${libpre}$2\.$libext)
        fi
        if test "$regExResult" = ""; then
            with_[$1]lib=""
            AC_MSG_RESULT("not found in $dirs")
        else
            with_[$1]lib=$regExResult
            AC_MSG_RESULT($with_[$1]lib)
        fi

        AC_MSG_CHECKING([for $3 ])
        dirs=""
        if test "$with_[$1]inc" != ""; then
            dirs=$with_[$1]inc
        elif test "$with_[$1]" != "yes"; then
            dirs=$with_[$1]
        else
            dirs="/usr/local/include /usr/local/gnu/include /usr/local/[$1] /opt/include /opt/gnu/include /opt/[$1] /usr/include"
        fi
        found=""
        for i in $dirs; do
            if test -d $i; then
                found="$found "`find $i -name patsubst([$3], .*/, ) -print 2> /dev/null`
            fi
        done
        AC_EXTRACT_REGEX($found, \(.*include\)/patsubst([$3], \., \\.))
        if test "$regExResult" = ""; then
            AC_EXTRACT_REGEX($found, \(.*\)/patsubst([$3], \., \\.))
        fi
        if test "$regExResult" = ""; then
            with_[$1]inc=""
            AC_MSG_RESULT("not found in $dirs")
        else
            with_[$1]inc=$regExResult
            AC_MSG_RESULT($with_[$1]inc)
        fi
        
        if test "$with_[$1]lib" = "" -o "$with_[$1]inc" = ""; then
            with_[$1]="no"
            AC_MSG_WARN(  Configuring without [$1] support)
        else
            with_[$1]="yes"
        fi
    fi
    
])

pushdef([VIGRA_PROG_INSTALL],
[
  dnl our own version, testing for the -p flag (--preserve-timestamp).
  dnl we first have to save if the user specified INSTALL as the
  dnl autoconf AC_PROG_INSTALL overwrites INSTALL:
  test -n "$INSTALL" && vigra_save_INSTALL_given=$INSTALL
  AC_PROG_INSTALL

  if test -z "$vigra_save_INSTALL_given" ; then
    # user hasn't overwritten INSTALL, autoconf found one for us
    # now we'll test if it supports the -p flag
    AC_MSG_CHECKING(for -p flag to install)
    rm -f confinst.$$.* > /dev/null 2>&1
    echo "Testtest" > confinst.$$.orig
    ac_res=no
    if ${INSTALL} -p confinst.$$.orig confinst.$$.new > /dev/null 2>&1 ; then
      if test -f confinst.$$.new ; then
        # OK, -p seems to do no harm to install
        INSTALL="${INSTALL} -p"
        ac_res=yes
      fi
    fi
    rm -f confinst.$$.*
    AC_MSG_RESULT($ac_res)
  fi
])dnl
