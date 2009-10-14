if test "x$2" != "x" ; then
    PATH=$2${PATH} ;
fi
$1 || { rm $1; exit 1; }
