PROG=$1
shift
until [ -z "$1" ]  # Until all parameters used up . . .
do
   PATH=$1:${PATH}
   shift
done
$PROG || { rm touch testsuccess.cxx; exit 1; }
