FILES = algorithm \
        deque \
	functional \
	hash_map \
	hash_set \
	iostream \
	iterator \
	list \
	map \
	memory \
	new \
	numeric \
	queue \
	set \
	stack \
	stdexcept \
	string \
	strstream \
	typeinfo \
	utility \
	vector

all::
	rm -f $(FILES)
	echo "#include <algo.h>" > algorithm
	echo "#include <deque.h>" > deque
	echo "#include <function.h>" > functional
	echo "#include <hash_map.h>" > hash_map
	echo "#include <hash_set.h>" > hash_set
	echo "#include <iostream.h>" > iostream
	echo "#include <iterator.h>" > iterator
	echo "#include <list.h>" > list
	(echo "#include <map.h>"; echo "#include <multimap.h>") > map
	(echo "#include <defalloc.h>"; \
         echo "#include <tempbuf.h>"; \
         echo "#include <iterator.h>"; \
         echo "#include <algo.h>") > memory
	echo "#include <new.h>" > new
	echo "#include <algo.h>" > numeric
	echo "#include <stack.h>" > queue
	(echo "#include <set.h>"; echo "#include <multiset.h>") > set
	echo "#include <stack.h>" > stack
	echo "#include <stdexcept.h>" > stdexcept
	echo "#include <mstring.h>" > string
	echo "#include <strstream.h>" > strstream
	echo "#include <typeinfo.h>" > typeinfo
	(echo "#include <pair.h>"; echo "#include <function.h>") > utility
	echo "#include <vector.h>" > vector
