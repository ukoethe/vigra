#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <stdexcept>
#include <vigra/hdf5impex.hxx>
#include <vigra/array_vector.hxx>

using namespace std;
using namespace vigra;

int main(int argc, char** argv)
{
	string line;
	ifstream a(argv[1], ios::in);
	ArrayVector<ArrayVector<string> > token;
	char delim;

	map<int, ArrayVector<double> > 	num_vecs;
	map<int, ArrayVector<string> > 	str_vecs; 
	map<int, string> 		  		col_nmes;
	string 							src_name = "crappy";
	string							hdf5_name = "asd.hdf5";

	if(argc < 3)
		delim = '|';
	else
		delim = argv[2][0];
	while(!(getline(a, line).eof()))
	{
		if(line.find(',') != string::npos && delim == '|')
			delim = ',';
		else if(delim == '|')
			delim = ' ';
		token.push_back(ArrayVector<string>());
		ArrayVector<string> & c_line = token.back();
		while(line.size() != 0)
		{
			if(line.find(delim) != string::npos)
			{
				c_line.push_back(line.substr(0, line.find(delim)));
				line = line.substr(line.find(delim)+1, line.size()); 
			}
			else
			{
				c_line.push_back(line);
				line = "";
			}
		}
		if(c_line.size() == 0)
		{
			token.pop_back();
		}
		else if(c_line.size() != token[0].size())
		{
			std::runtime_error("unequal columns in row");
		}
	}
	int rowCt = token.size();
	int colCt = token[0].size(); 
	
	ArrayVector<int> aie(size_t(colCt), 0); 

	for(int ii = 0; ii < colCt; ++ii)
	{
		ArrayVector<string> & c_line = token[0]; 
		istringstream sin(c_line[ii]);
		double in;
		sin >> in;
		if(sin.fail())
		{
			aie[ii] = 1;
			str_vecs[ii] = ArrayVector<string>(1, c_line[ii]); 
		}
		else
		{
			aie[ii] = 0; 
			num_vecs[ii] = ArrayVector<double>(1, in);
		}
	}

	for(int ii = 0; ii < colCt; ++ii)
	{
		if(aie[ii])
		{
			for(int jj = 1; jj < rowCt; ++jj)
			{
				str_vecs[ii].push_back(token[jj][ii]);
			}
		}
		else
		{
			for(int jj = 1; jj < rowCt; ++jj)
			{
				ArrayVector<string> & c_line = token[jj];
				istringstream sin(c_line[ii]);
				double in;
				sin >> in;
				num_vecs[ii].push_back(in);
			}
		}
	}
	{
		map<int, ArrayVector<double> >::iterator iter;
		for(iter = num_vecs.begin(); iter != num_vecs.end(); ++iter)
		{
			MultiArrayView<2, double> stuff(MultiArrayShape<2>::type((iter->second).size(),1), 
											(iter->second).data());
			ostringstream sout;
			sout << src_name << "/" << iter->first;
			writeToHDF5File(hdf5_name.c_str(), sout.str().c_str(), stuff);
		}
	}
	{
		map<int, ArrayVector<string> >::iterator iter;
		for(iter = str_vecs.begin(); iter != str_vecs.end(); ++iter)
		{
			MultiArrayView<1, std::string> stuff(MultiArrayShape<1>::type((iter->second).size()), 
											(iter->second).data());
			ostringstream sout;
			sout << src_name << "/" << iter->first;
			std::cerr << sout.str().c_str() <<std::endl;
			writeToHDF5File(hdf5_name.c_str(), sout.str().c_str(), stuff);
		}
	}
	return 0;
}
