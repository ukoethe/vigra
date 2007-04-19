# redefine _AC_INIT_DEFAULTS in order to add 'ac_default_docdir=...'
# --------------------------------------------------------------
# Values which defaults can be set from `configure.ac'.
# `/bin/machine' is used in `glibcbug'.  The others are used in config.*
m4_define([_AC_INIT_DEFAULTS],
[m4_divert_push([DEFAULTS])dnl

# Name of the host.
# hostname on some systems (SVR3.2, Linux) returns a bogus exit status,
# so uname gets run too.
ac_hostname=`(hostname || uname -n) 2>/dev/null | sed 1q`

exec AS_MESSAGE_FD>&1

#
# Initializations.
#
ac_default_prefix=/usr/local
ac_default_docdir='${prefix}/doc/vigra'
ac_config_libobj_dir=.
cross_compiling=no
subdirs=
MFLAGS=
MAKEFLAGS=
AC_SUBST([SHELL], [${CONFIG_SHELL-/bin/sh}])dnl
AC_SUBST([PATH_SEPARATOR])dnl

# Maximum number of lines to put in a shell here document.
# This variable seems obsolete.  It should probably be removed, and
# only ac_max_sed_lines should be used.
: ${ac_max_here_lines=38}

# Identity of this package.
AC_SUBST([PACKAGE_NAME],
	 [m4_ifdef([AC_PACKAGE_NAME],      ['AC_PACKAGE_NAME'])])dnl
AC_SUBST([PACKAGE_TARNAME],
	 [m4_ifdef([AC_PACKAGE_TARNAME],   ['AC_PACKAGE_TARNAME'])])dnl
AC_SUBST([PACKAGE_VERSION],
	 [m4_ifdef([AC_PACKAGE_VERSION],   ['AC_PACKAGE_VERSION'])])dnl
AC_SUBST([PACKAGE_STRING],
	 [m4_ifdef([AC_PACKAGE_STRING],    ['AC_PACKAGE_STRING'])])dnl
AC_SUBST([PACKAGE_BUGREPORT],
	 [m4_ifdef([AC_PACKAGE_BUGREPORT], ['AC_PACKAGE_BUGREPORT'])])dnl

m4_divert_pop([DEFAULTS])dnl
m4_wrap([m4_divert_text([DEFAULTS],
[ac_subst_vars='m4_ifdef([_AC_SUBST_VARS],  [m4_defn([_AC_SUBST_VARS])])'
ac_subst_files='m4_ifdef([_AC_SUBST_FILES], [m4_defn([_AC_SUBST_FILES])])'])])dnl
])# _AC_INIT_DEFAULTS

#########################################################

# AC_DOCDIR_DEFAULT(DOCDIR)
# -------------------------
AC_DEFUN([VIGRA_DOCDIR_DEFAULT],
[m4_divert_text([DEFAULTS], [ac_default_docdir=$1])])

#########################################################

# redefine _AC_INIT_HELP in order to get better control over
# default settings for output directories
# ----------------------------------------------------------
# Handle the `configure --help' message.
m4_define([_AC_INIT_HELP],
[m4_divert_push([HELP_BEGIN])dnl

#
# Report the --help message.
#
if test "$ac_init_help" = "long"; then
  # Omit some internal or obsolete options to make the list less imposing.
  # This message is too long to be a string in the A/UX 3.1 sh.
  cat <<_ACEOF
\`configure' configures m4_ifset([AC_PACKAGE_STRING],
                        [AC_PACKAGE_STRING],
                        [this package]) to adapt to many kinds of systems.

Usage: $[0] [[OPTION]]... [[VAR=VALUE]]...

[To assign environment variables (e.g., CC, CFLAGS...), specify them as
VAR=VALUE.  See below for descriptions of some of the useful variables.

Defaults for the options are specified in brackets.

Configuration:
  -h, --help              display this help and exit
      --help=short        display options specific to this package
      --help=recursive    display the short help of all the included packages
  -V, --version           display version information and exit
  -q, --quiet, --silent   do not print \`checking...' messages
      --cache-file=FILE   cache test results in FILE [disabled]
  -C, --config-cache      alias for \`--cache-file=config.cache'
  -n, --no-create         do not create output files
      --srcdir=DIR        find the sources in DIR [configure dir or \`..']

_ACEOF

  cat <<_ACEOF
Installation directories:
  --prefix=PREFIX         install architecture-independent files in PREFIX
                          [$ac_default_prefix]
  --exec-prefix=EPREFIX   install architecture-dependent files in EPREFIX
                          [PREFIX]

By default, \`make install' will install all the files in
\`$ac_default_prefix/bin', \`$ac_default_prefix/lib' etc.  You can specify
an installation prefix other than \`$ac_default_prefix' using \`--prefix',
for instance \`--prefix=\$HOME'.

For better control, use the options below.

Fine tuning of the installation directories:
  --bindir=DIR           user executables [EPREFIX/bin]
  --libdir=DIR           object code libraries [EPREFIX/lib]
  --includedir=DIR       C header files [PREFIX/include]
  --docdir=DIR           html documentation [$ac_default_docdir]
_ACEOF

  cat <<\_ACEOF]
m4_divert_pop([HELP_BEGIN])dnl
dnl The order of the diversions here is
dnl - HELP_BEGIN
dnl   which may be prolongated by extra generic options such as with X or
dnl   AC_ARG_PROGRAM.  Displayed only in long --help.
dnl
dnl - HELP_CANON
dnl   Support for cross compilation (--build, --host and --target).
dnl   Display only in long --help.
dnl
dnl - HELP_ENABLE
dnl   which starts with the trailer of the HELP_BEGIN, HELP_CANON section,
dnl   then implements the header of the non generic options.
dnl
dnl - HELP_WITH
dnl
dnl - HELP_VAR
dnl
dnl - HELP_VAR_END
dnl
dnl - HELP_END
dnl   initialized below, in which we dump the trailer (handling of the
dnl   recursion for instance).
m4_divert_push([HELP_ENABLE])dnl
_ACEOF
fi

if test -n "$ac_init_help"; then
m4_ifset([AC_PACKAGE_STRING],
[  case $ac_init_help in
     short | recursive ) echo "Configuration of AC_PACKAGE_STRING:";;
   esac])
  cat <<\_ACEOF
m4_divert_pop([HELP_ENABLE])dnl
m4_divert_push([HELP_END])dnl
m4_ifset([AC_PACKAGE_BUGREPORT], [
Report bugs to <AC_PACKAGE_BUGREPORT>.])
_ACEOF
fi

if test "$ac_init_help" = "recursive"; then
  # If there are subdirs, report their specific --help.
  ac_popdir=`pwd`
  for ac_dir in : $ac_subdirs_all; do test "x$ac_dir" = x: && continue
    test -d $ac_dir || continue
    _AC_SRCPATHS(["$ac_dir"])
    cd $ac_dir
    # Check for guested configure; otherwise get Cygnus style configure.
    if test -f $ac_srcdir/configure.gnu; then
      echo
      $SHELL $ac_srcdir/configure.gnu  --help=recursive
    elif test -f $ac_srcdir/configure; then
      echo
      $SHELL $ac_srcdir/configure  --help=recursive
    elif test -f $ac_srcdir/configure.ac ||
           test -f $ac_srcdir/configure.in; then
      echo
      $ac_configure --help
    else
      AC_MSG_WARN([no configuration information is in $ac_dir])
    fi
    cd $ac_popdir
  done
fi

test -n "$ac_init_help" && exit 0
m4_divert_pop([HELP_END])dnl
])# _AC_INIT_HELP

#########################################################


# redefine _AC_INIT_PARSE_ARGS to add 'docdir'
# --------------------------------------------
m4_define([_AC_INIT_PARSE_ARGS],
[m4_divert_push([PARSE_ARGS])dnl

# Initialize some variables set by options.
ac_init_help=
ac_init_version=false
# The variables have the same names as the options, with
# dashes changed to underlines.
cache_file=/dev/null
AC_SUBST(exec_prefix, NONE)dnl
no_create=
no_recursion=
AC_SUBST(prefix, NONE)dnl
program_prefix=NONE
program_suffix=NONE
AC_SUBST(program_transform_name, [s,x,x,])dnl
silent=
site=
srcdir=
verbose=
x_includes=NONE
x_libraries=NONE

# Installation directory options.
# These are left unexpanded so users can "make install exec_prefix=/foo"
# and all the variables that are supposed to be based on exec_prefix
# by default will actually change.
# Use braces instead of parens because sh, perl, etc. also accept them.
AC_SUBST([bindir],         ['${exec_prefix}/bin'])dnl
AC_SUBST([libdir],         ['${exec_prefix}/lib'])dnl
AC_SUBST([includedir],     ['${prefix}/include'])dnl
AC_SUBST([docdir],         [$ac_default_docdir])dnl

ac_prev=
for ac_option
do
  # If the previous option needs an argument, assign it.
  if test -n "$ac_prev"; then
    eval "$ac_prev=\$ac_option"
    ac_prev=
    continue
  fi

  ac_optarg=`expr "x$ac_option" : 'x[[^=]]*=\(.*\)'`

  # Accept the important Cygnus configure options, so we can diagnose typos.

  case $ac_option in

  -bindir | --bindir | --bindi | --bind | --bin | --bi)
    ac_prev=bindir ;;
  -bindir=* | --bindir=* | --bindi=* | --bind=* | --bin=* | --bi=*)
    bindir=$ac_optarg ;;

  -build | --build | --buil | --bui | --bu)
    ac_prev=build_alias ;;
  -build=* | --build=* | --buil=* | --bui=* | --bu=*)
    build_alias=$ac_optarg ;;

  -cache-file | --cache-file | --cache-fil | --cache-fi \
  | --cache-f | --cache- | --cache | --cach | --cac | --ca | --c)
    ac_prev=cache_file ;;
  -cache-file=* | --cache-file=* | --cache-fil=* | --cache-fi=* \
  | --cache-f=* | --cache-=* | --cache=* | --cach=* | --cac=* | --ca=* | --c=*)
    cache_file=$ac_optarg ;;

  --config-cache | -C)
    cache_file=config.cache ;;

  -disable-* | --disable-*)
    ac_feature=`expr "x$ac_option" : 'x-*disable-\(.*\)'`
    # Reject names that are not valid shell variable names.
    expr "x$ac_feature" : "[.*[^-_$as_cr_alnum]]" >/dev/null &&
      AC_MSG_ERROR([invalid feature name: $ac_feature])
    ac_feature=`echo $ac_feature | sed 's/-/_/g'`
    eval "enable_$ac_feature=no" ;;

  -docdir | --docdir | --docdi | --docd | --doc | --do)
    ac_prev=docdir ;;
  -docdir=* | --docdir=* | --docdi=* | --docd=* | --doc=* | --do=* )
    docdir=$ac_optarg ;;

  -enable-* | --enable-*)
    ac_feature=`expr "x$ac_option" : 'x-*enable-\([[^=]]*\)'`
    # Reject names that are not valid shell variable names.
    expr "x$ac_feature" : "[.*[^-_$as_cr_alnum]]" >/dev/null &&
      AC_MSG_ERROR([invalid feature name: $ac_feature])
    ac_feature=`echo $ac_feature | sed 's/-/_/g'`
    case $ac_option in
      *=*) ac_optarg=`echo "$ac_optarg" | sed "s/'/'\\\\\\\\''/g"`;;
      *) ac_optarg=yes ;;
    esac
    eval "enable_$ac_feature='$ac_optarg'" ;;

  -exec-prefix | --exec_prefix | --exec-prefix | --exec-prefi \
  | --exec-pref | --exec-pre | --exec-pr | --exec-p | --exec- \
  | --exec | --exe | --ex)
    ac_prev=exec_prefix ;;
  -exec-prefix=* | --exec_prefix=* | --exec-prefix=* | --exec-prefi=* \
  | --exec-pref=* | --exec-pre=* | --exec-pr=* | --exec-p=* | --exec-=* \
  | --exec=* | --exe=* | --ex=*)
    exec_prefix=$ac_optarg ;;

  -gas | --gas | --ga | --g)
    # Obsolete; use --with-gas.
    with_gas=yes ;;

  -help | --help | --hel | --he | -h)
    ac_init_help=long ;;
  -help=r* | --help=r* | --hel=r* | --he=r* | -hr*)
    ac_init_help=recursive ;;
  -help=s* | --help=s* | --hel=s* | --he=s* | -hs*)
    ac_init_help=short ;;

  -host | --host | --hos | --ho)
    ac_prev=host_alias ;;
  -host=* | --host=* | --hos=* | --ho=*)
    host_alias=$ac_optarg ;;

  -includedir | --includedir | --includedi | --included | --include \
  | --includ | --inclu | --incl | --inc)
    ac_prev=includedir ;;
  -includedir=* | --includedir=* | --includedi=* | --included=* | --include=* \
  | --includ=* | --inclu=* | --incl=* | --inc=*)
    includedir=$ac_optarg ;;

  -libdir | --libdir | --libdi | --libd)
    ac_prev=libdir ;;
  -libdir=* | --libdir=* | --libdi=* | --libd=*)
    libdir=$ac_optarg ;;

  -nfp | --nfp | --nf)
    # Obsolete; use --without-fp.
    with_fp=no ;;

  -no-create | --no-create | --no-creat | --no-crea | --no-cre \
  | --no-cr | --no-c | -n)
    no_create=yes ;;

  -no-recursion | --no-recursion | --no-recursio | --no-recursi \
  | --no-recurs | --no-recur | --no-recu | --no-rec | --no-re | --no-r)
    no_recursion=yes ;;

  -prefix | --prefix | --prefi | --pref | --pre | --pr | --p)
    ac_prev=prefix ;;
  -prefix=* | --prefix=* | --prefi=* | --pref=* | --pre=* | --pr=* | --p=*)
    prefix=$ac_optarg ;;

  -program-prefix | --program-prefix | --program-prefi | --program-pref \
  | --program-pre | --program-pr | --program-p)
    ac_prev=program_prefix ;;
  -program-prefix=* | --program-prefix=* | --program-prefi=* \
  | --program-pref=* | --program-pre=* | --program-pr=* | --program-p=*)
    program_prefix=$ac_optarg ;;

  -program-suffix | --program-suffix | --program-suffi | --program-suff \
  | --program-suf | --program-su | --program-s)
    ac_prev=program_suffix ;;
  -program-suffix=* | --program-suffix=* | --program-suffi=* \
  | --program-suff=* | --program-suf=* | --program-su=* | --program-s=*)
    program_suffix=$ac_optarg ;;

  -program-transform-name | --program-transform-name \
  | --program-transform-nam | --program-transform-na \
  | --program-transform-n | --program-transform- \
  | --program-transform | --program-transfor \
  | --program-transfo | --program-transf \
  | --program-trans | --program-tran \
  | --progr-tra | --program-tr | --program-t)
    ac_prev=program_transform_name ;;
  -program-transform-name=* | --program-transform-name=* \
  | --program-transform-nam=* | --program-transform-na=* \
  | --program-transform-n=* | --program-transform-=* \
  | --program-transform=* | --program-transfor=* \
  | --program-transfo=* | --program-transf=* \
  | --program-trans=* | --program-tran=* \
  | --progr-tra=* | --program-tr=* | --program-t=*)
    program_transform_name=$ac_optarg ;;

  -q | -quiet | --quiet | --quie | --qui | --qu | --q \
  | -silent | --silent | --silen | --sile | --sil)
    silent=yes ;;

  -site | --site | --sit)
    ac_prev=site ;;
  -site=* | --site=* | --sit=*)
    site=$ac_optarg ;;

  -srcdir | --srcdir | --srcdi | --srcd | --src | --sr)
    ac_prev=srcdir ;;
  -srcdir=* | --srcdir=* | --srcdi=* | --srcd=* | --src=* | --sr=*)
    srcdir=$ac_optarg ;;

  -target | --target | --targe | --targ | --tar | --ta | --t)
    ac_prev=target_alias ;;
  -target=* | --target=* | --targe=* | --targ=* | --tar=* | --ta=* | --t=*)
    target_alias=$ac_optarg ;;

  -v | -verbose | --verbose | --verbos | --verbo | --verb)
    verbose=yes ;;

  -version | --version | --versio | --versi | --vers | -V)
    ac_init_version=: ;;

  -with-* | --with-*)
    ac_package=`expr "x$ac_option" : 'x-*with-\([[^=]]*\)'`
    # Reject names that are not valid shell variable names.
    expr "x$ac_package" : "([.*[^-_$as_cr_alnum]])" >/dev/null &&
      AC_MSG_ERROR([invalid package name: $ac_package])
    ac_package=`echo $ac_package| sed 's/-/_/g'`
    case $ac_option in
      *=*) ac_optarg=`echo "$ac_optarg" | sed "s/'/'\\\\\\\\''/g"`;;
      *) ac_optarg=yes ;;
    esac
    eval "with_$ac_package='$ac_optarg'" ;;

  -without-* | --without-*)
    ac_package=`expr "x$ac_option" : 'x-*without-\(.*\)'`
    # Reject names that are not valid shell variable names.
    expr "x$ac_package" : "([.*[^-_$as_cr_alnum]])" >/dev/null &&
      AC_MSG_ERROR([invalid package name: $ac_package])
    ac_package=`echo $ac_package | sed 's/-/_/g'`
    eval "with_$ac_package=no" ;;

  --x)
    # Obsolete; use --with-x.
    with_x=yes ;;

  -x-includes | --x-includes | --x-include | --x-includ | --x-inclu \
  | --x-incl | --x-inc | --x-in | --x-i)
    ac_prev=x_includes ;;
  -x-includes=* | --x-includes=* | --x-include=* | --x-includ=* | --x-inclu=* \
  | --x-incl=* | --x-inc=* | --x-in=* | --x-i=*)
    x_includes=$ac_optarg ;;

  -x-libraries | --x-libraries | --x-librarie | --x-librari \
  | --x-librar | --x-libra | --x-libr | --x-lib | --x-li | --x-l)
    ac_prev=x_libraries ;;
  -x-libraries=* | --x-libraries=* | --x-librarie=* | --x-librari=* \
  | --x-librar=* | --x-libra=* | --x-libr=* | --x-lib=* | --x-li=* | --x-l=*)
    x_libraries=$ac_optarg ;;

  -*) AC_MSG_ERROR([unrecognized option: $ac_option
Try `$[0] --help' for more information.])
    ;;

  *=*)
    ac_envvar=`expr "x$ac_option" : 'x\([[^=]]*\)='`
    # Reject names that are not valid shell variable names.
    expr "x$ac_envvar" : "[.*[^_$as_cr_alnum]]" >/dev/null &&
      AC_MSG_ERROR([invalid variable name: $ac_envvar])
    ac_optarg=`echo "$ac_optarg" | sed "s/'/'\\\\\\\\''/g"`
    eval "$ac_envvar='$ac_optarg'"
    export $ac_envvar ;;

  *)
    # FIXME: should be removed in autoconf 3.0.
    AC_MSG_WARN([you should use --build, --host, --target])
    expr "x$ac_option" : "[.*[^-._$as_cr_alnum]]" >/dev/null &&
      AC_MSG_WARN([invalid host type: $ac_option])
    : ${build_alias=$ac_option} ${host_alias=$ac_option} ${target_alias=$ac_option}
    ;;

  esac
done

if test -n "$ac_prev"; then
  ac_option=--`echo $ac_prev | sed 's/_/-/g'`
  AC_MSG_ERROR([missing argument to $ac_option])
fi

# Be sure to have absolute paths.
for ac_var in exec_prefix prefix
do
  eval ac_val=$`echo $ac_var`
  case $ac_val in
    [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
    *)  AC_MSG_ERROR([expected an absolute directory name for --$ac_var: $ac_val]);;
  esac
done

# Be sure to have absolute paths.
for ac_var in bindir libdir includedir docdir
do
  eval ac_val=$`echo $ac_var`
  case $ac_val in
    [[\\/$]]* | ?:[[\\/]]* ) ;;
    *)  AC_MSG_ERROR([expected an absolute directory name for --$ac_var: $ac_val]);;
  esac
done

# There might be people who depend on the old broken behavior: `$host'
# used to hold the argument of --host etc.
# FIXME: To remove some day.
build=$build_alias
host=$host_alias
target=$target_alias

# FIXME: To remove some day.
if test "x$host_alias" != x; then
  if test "x$build_alias" = x; then
    cross_compiling=maybe
    AC_MSG_WARN([If you wanted to set the --build type, don't use --host.
    If a cross compiler is detected then cross compile mode will be used.])
  elif test "x$build_alias" != "x$host_alias"; then
    cross_compiling=yes
  fi
fi

ac_tool_prefix=
test -n "$host_alias" && ac_tool_prefix=$host_alias-

test "$silent" = yes && exec AS_MESSAGE_FD>/dev/null

m4_divert_pop([PARSE_ARGS])dnl
])# _AC_INIT_PARSE_ARGS

dnl-------------------------------------------------------------------
dnl      VIGRA external package configuration macros
dnl-------------------------------------------------------------------

dnl VIGRA_FIND_PATH(srcvariable, targetvariable, regexPattern)
dnl example:
dnl     VIGRA_FIND_PATH(CPPFLAGS, INCLUDEPATH, [-I\(.*\)])
AC_DEFUN([VIGRA_FIND_PATH],
[
for VIGRA_FIND_PATH_i in $[$1];  do
    VIGRA_FIND_PATH_r=`expr "$VIGRA_FIND_PATH_i" : "$3"`
    if test -n "$VIGRA_FIND_PATH_r"; then
        [$2]="$[$2] $VIGRA_FIND_PATH_r"
    fi
done
])dnl VIGRA_FIND_PATH

dnl VIGRA_DECLARE_PACKAGE(packageName, comment, defaultSetting)
dnl declare command line switches and comment for the given package
dnl example:
dnl     VIGRA_FIND_PATH(tiff, [import TIFF images], "yes")
AC_DEFUN([VIGRA_DECLARE_PACKAGE],
[
    if test "x$3" = "x"; then
        with_[$1]default="yes"
    else
        with_[$1]default="[$3]"
    fi
    # $1 is the package name, $2 the description, $3 the default setting
    AC_ARG_WITH([$1], [
  --with-$1
  --with-$1=dir
  --without-$1
  --with-$1lib=dir
  --with-$1inc=dir
      $2.
      default: --with-$1=$3
      if --with-$1 or --with-$1=yes is given: $1 package files will be
         searched for by 'pkg-config' and in some standard directories (the default).
      if --with-$1=dir is given, and dir is a directory: $1 library files
         will be looked up in 'dir'.
      if --with-$1=no or --without-$1 is given: $1 package will
         not be used.
      alternatively, you can specify:], ,[with_[$1]="$with_[$1]default"])
    AC_ARG_WITH([$1lib], [        --with-$1lib=dir : the $1 package's lib directory], ,)
    AC_ARG_WITH([$1inc], [        --with-$1inc=dir : the $1 package's include directory], ,)

    # bring the user's response into a standard form
    # if any flag is "no", the main flag should be "no", the others empty
    if test "x$with_[$1]" = "xno" -o "x$with_[$1]lib" = "xno" -o "x$with_[$1]inc" = "xno"; then
        with_[$1]="no"
        with_[$1]lib=""
        with_[$1]inc=""
    else
        # here, the main flag is either "yes" or a path
        # treat "yes" in the auxilliary flags as ""
        if test "x$with_[$1]lib" = "xyes"; then
            with_[$1]lib=""
        fi
        if test "x$with_[$1]inc" = "xyes"; then
            with_[$1]inc=""
        fi
        # use a path given in the main flag as the library path
        if test "$with_[$1]" != "yes" -a "$with_[$1]" != "$with_[$1]default" -a "x$with_[$1]lib" = "x"; then
            with_[$1]lib="$with_[$1]"
        fi
        # remember in the main flag whether an explicit path was given
        if test "x$with_[$1]lib" != "x" -o "x$with_[$1]inc" != "x"; then
            with_[$1]="explicit"
        else
            with_[$1]="search"
        fi
    fi
    # now $with_[$1] is either "no", "search", or "explicit"
    # $with_[$1]lib and $with_[$1]inc are either "" or a path
])dnl VIGRA_DECLARE_PACKAGE

dnl VIGRA_FIND_PACKAGE_PKGCONFIG(packageName, pkgconfigName, includeFile)
dnl search for package using pkg-config (or $PKGCONFIG)
dnl on success, $with_packageName is "yes", and $with_package_I, $with_package_L, $with_package_l are set
dnl             otherwise $with_packageName is unchanged (i.e. "search")
dnl example:
dnl     VIGRA_FIND_PACKAGE_PKGCONFIG(fftw, fftw3, "fftw3.h")
AC_DEFUN([VIGRA_FIND_PACKAGE_PKGCONFIG],
[
    if test "x$PKGCONFIG" != "x"; then
        AC_MSG_CHECKING([if pkg-config knows about $1 ])
        if $PKGCONFIG --exists "$2"; then
            # save the state
            SAVECPPFLAGS="$CPPFLAGS"
            SAVELDFLAGS="$LDFLAGS"
            SAVELIBS="$LIBS"
            AC_LANG_SAVE
            AC_LANG_CPLUSPLUS

            # $1 is the package, $2 the pkg-config name, $3 the include name
            CPPFLAGS=`$PKGCONFIG --cflags-only-I "$2"`
            LDFLAGS=`$PKGCONFIG --libs-only-L "$2"`
            LIBS=`$PKGCONFIG --libs-only-l "$2"`
            AC_TRY_LINK(
                    [#include <stdio.h>  /* for jpeglib.h */
                     #include <$3>
], [],
                    [with_[$1]="yes"
                     with_[$1]_I="$CPPFLAGS"
                     with_[$1]_L="$LDFLAGS"
                     with_[$1]_l="$LIBS"
                    ],
                    [with_[$1]_I=""
                     with_[$1]_L=""
                     with_[$1]_l=""
                    ])
            CPPFLAGS="$CPPFLAGS $SAVECPPFLAGS"
            LDFLAGS="$LDFLAGS $SAVELDFLAGS"
            LIBS="$LIBS $SAVELIBS"
            AC_LANG_RESTORE
        fi
        if test "x$with_[$1]" = "xyes"; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no (check \$PKG_CONFIG_PATH)])
        fi
    fi
])dnl VIGRA_FIND_PACKAGE_PKGCONFIG

dnl VIGRA_FIND_LIBRARY_EXPLICIT(packageName, libraryName, path)
dnl search for a library in the given path
dnl example:
dnl     VIGRA_FIND_LIBRARY_EXPLICIT(zlib, z, /usr/local/lib)
AC_DEFUN([VIGRA_FIND_LIBRARY_EXPLICIT],
[
    # save the state
    SAVELDFLAGS="$LDFLAGS"
    SAVELIBS="$LIBS"
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS

    # $1 is the package, $2 the library name, $3 the path
    if test "x$3" != "x" ; then
        LDFLAGS="-L$3 $SAVELDFLAGS"
        LIBS="-l[$2] $SAVELIBS"
        AC_TRY_LINK([],
                    [],
            [with_[$1]lib="$3"
             with_[$1]_L="$LDFLAGS"
             with_[$1]_l="$LIBS"],
            with_[$1]lib="")
    else
        LIBS="-l[$2] $SAVELIBS"
        AC_TRY_LINK([],
                    [],
            [with_[$1]lib="in default path"
             with_[$1]_L="$LDFLAGS"
             with_[$1]_l="$LIBS"],
            with_[$1]lib="")
    fi
    # restore the state
    AC_LANG_RESTORE
    LDFLAGS="$SAVELDFLAGS"
    LIBS="$SAVELIBS"
])dnl VIGRA_FIND_LIBRARY_EXPLICIT

dnl VIGRA_FIND_INCLUDE_EXPLICIT(packageName, includeName, path)
dnl search for an include in the given path
dnl example:
dnl     VIGRA_FIND_INCLUDE_EXPLICIT(jpeg, jpeglib.h, /usr/local/include)
AC_DEFUN([VIGRA_FIND_INCLUDE_EXPLICIT],
[
    # save the state
    SAVECPPFLAGS="$CPPFLAGS"
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS

    # $1 is the package, $2 the include name, $3 the path
    if test "x$3" != "x" ; then
        CPPFLAGS="-I$3 $SAVECPPFLAGS"
        AC_TRY_COMPILE(
                [#include <stdio.h>  /* for jpeglib.h */
                 #include <$2>
], [],
            [with_[$1]inc="$3"
             with_[$1]_I="$CPPFLAGS"],
             with_[$1]inc="")
    else
         AC_TRY_COMPILE(
                [#include <stdio.h>  /* for jpeglib.h */
                 #include <$2>
], [],
            [with_[$1]inc="in default path"
             with_[$1]_I="$CPPFLAGS"],
             with_[$1]inc="")
    fi
    # restore the state
    AC_LANG_RESTORE
    CPPFLAGS="$SAVECPPFLAGS"
])dnl VIGRA_FIND_INCLUDE_EXPLICIT

dnl VIGRA_FIND_LIBRARY(packageName, libraryName)
dnl search for a library in a number of paths
dnl if $with_packageName is "search" the following paths are tried
dnl              "", all in $LIBSEARCHPATH, /usr/local/lib/packageName, /usr/local/lib
dnl if $with_packageName is "explicit" the path in $with_packageNamelib is tried
dnl on success, $with_packageNamelib is either "in default path" or contains the path
dnl             otherwise it is empty
dnl example:
dnl     VIGRA_FIND_LIBRARY(jpeg, jpeg)
AC_DEFUN([VIGRA_FIND_LIBRARY],
[
    AC_MSG_CHECKING([for lib$2 ])

    # if we got an explicit path, that's the only one we check
    if test "x$with_[$1]" = "xexplicit" -a "x$with_[$1]lib" != "x"; then
        VIGRA_FIND_LIBRARY_EXPLICIT($1, $2, "$with_[$1]lib")
    else
        # otherwise check whether the library is found without explicit path
        # specification (e.g. because it is in /usr/lib)
        if test "x$with_[$1]lib" = "x" ; then
            VIGRA_FIND_LIBRARY_EXPLICIT($1, $2, "")
        fi
        # if this didn't work check whether the library path was in the
        # original $LDFLAGS or a few standard directories
        if test "x$with_[$1]lib" = "x" ; then
            for i in $LIBSEARCHPATH /usr/local/lib/[$1] /usr/local/lib; do
                VIGRA_FIND_LIBRARY_EXPLICIT($1, $2, "$i")
                if test "x$with_[$1]lib" != "x" ; then break; fi
            done
        fi
    fi

    # report the results back
    if test "x$with_[$1]lib" = "x"; then
        AC_MSG_RESULT([not found])
        with_[$1]="no"
    else
        AC_MSG_RESULT($with_[$1]lib)
        if test "$with_[$1]lib" = "in default path"; then
            with_[$1]lib=""
        fi
    fi
])dnl VIGRA_FIND_LIBRARY

dnl VIGRA_FIND_INCLUDE(packageName, includeName, pathPostfix)
dnl search for an include in a number of paths
dnl if $with_packageName is "search" the following paths are tried
dnl              "", include directories in the same directory tree as $with_packageNamelib,
dnl              all in $LIBSEARCHPATH, /usr/local/include/packageName, /usr/local/include
dnl if $with_packageName is "explicit" the path in $with_packageNamenc is tried
dnl on success, $with_packageNameinc is either "in default path" or contains the path
dnl             otherwise it is empty
dnl example:
dnl     VIGRA_FIND_INCLUDE(jpeg, jpeglib.h, )
AC_DEFUN([VIGRA_FIND_INCLUDE],
[
    AC_MSG_CHECKING([for $2 ])

    # first check whether we got an explicit path
    if test "x$with_[$1]" = "xexplicit" -a "x$with_[$1]inc" != "x"; then
        VIGRA_FIND_INCLUDE_EXPLICIT($1, $2, "$with_[$1]inc")
    else
        # otherwise check whether the library is found without explicit path
        # specification (e.g. because it is in /usr/lib)
        if test "x$with_[$1]inc" = "x" ; then
            VIGRA_FIND_INCLUDE_EXPLICIT($1, $2, "")
        fi
        # if this didn't work check whether the library path is in the
        # same tree as the library path, or among the directories in the
        # original $CPPFLAGS or a few standard directories
        if test -n "$3"; then
            PATHPOSTFIX="/$3"
        else
            PATHPOSTFIX=""
        fi
        EXTRAINCSEARCHPATH=""
        i="$with_[$1]lib"
        while test "x$i" != "x"; do
            if test -d "${i}/include"; then
                EXTRAINCSEARCHPATH="$EXTRAINCSEARCHPATH ${i}/include"
            fi
            j=`expr "$i" : "\(.*\)/"`
            i=$j
        done
        if test "x$with_[$1]inc" = "x" ; then
            for i in $EXTRAINCSEARCHPATH $INCSEARCHPATH /usr/local/include/[$1] /usr/local/include; do
                 VIGRA_FIND_INCLUDE_EXPLICIT($1, $2, "$i$PATHPOSTFIX")
                if test "x$with_[$1]inc" != "x" ; then break; fi
            done
        fi
    fi

    # report the results back
    if test "x$with_[$1]inc" = "x"; then
        AC_MSG_RESULT([not found])
        with_[$1]="no"
        with_[$1]lib=""
    else
        AC_MSG_RESULT($with_[$1]inc)
        if test "$with_[$1]inc" = "in default path"; then
            with_[$1]inc=""
        fi
    fi
])dnl VIGRA_FIND_INCLUDE

dnl VIGRA_FIND_PACKAGE(packageName, libraryName, includeName)
dnl call VIGRA_FIND_LIBARY(packageName, libraryName) and VIGRA_FIND_INCLUDE(packageName, includeName)
dnl on success, $with_packageName is "yes", and $with_package_I, $with_package_L, $with_package_l are set
dnl             otherwise $with_packageName is "no"
dnl example:
dnl     VIGRA_FIND_PACKAGE(jpeg, jpeg, jpeglib.h)
AC_DEFUN([VIGRA_FIND_PACKAGE],
[
    if test ${with_[$1]:-""} != "no"; then
        VIGRA_FIND_LIBRARY($1, $2)
        if test ${with_[$1]:-""} != "no"; then
            VIGRA_FIND_INCLUDE($1, $3)
            if test ${with_[$1]:-""} != "no"; then
                with_[$1]="yes"
            fi
        fi
    fi

])dnl

#########################################################

AC_DEFUN([VIGRA_PROG_INSTALL],
[
  dnl our own version, testing for the -p flag (--preserve-timestamp).
  dnl we first have to save if the user specified INSTALL as the
  dnl autoconf AC_PROG_INSTALL overwrites INSTALL:
  test -n "$INSTALL" && vigra_save_INSTALL_given=$INSTALL
  AC_PROG_INSTALL

  INSTALL_MKDIR="$INSTALL -d"
  AC_SUBST(INSTALL_MKDIR)

  if test -z "$vigra_save_INSTALL_given" ; then
    # user hasn't overwritten INSTALL, autoconf found one for us
    # now we'll test if it supports the -p flag
    AC_MSG_CHECKING(whether $INSTALL accepts the -p flag)
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
