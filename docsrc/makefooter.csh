#!/bin/csh -f

setenv LC_TIME C
set today = `date +"%b %d %Y"`
cat footer.dxx | nawk '$1=="TODAY_TODAY_TODAY"  {print "VIGRA", version, "(", today,")"; next;} {print;}' today="$today" version="$1" > footer.html
