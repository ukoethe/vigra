#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <stdexcept>
#include <vigra/hdf5impex.hxx>
#include <vigra/array_vector.hxx>
#include <tclap/CmdLine.h>
#include <vigra/hdf5impex.hxx>
#include <set>
using namespace std;
using namespace vigra;

enum IDS_ {NUMERAL, LITERAL, ORDINAL, CATEGORICAL, IFSTRING, IN, INOUT, OUT};
void ARFF2HDF5(std::string arff, std::string hdf5)
{
}

template<class In>
void writeHDF5attr(std::string filename, std::string path, In in)
{
}

void tokenize(std::string  str, ArrayVector<std::string> & in, std::string delim = ",")
{
    if(str.size() < delim.size() || str.rfind(delim) != str.size() - delim.size())
        str = str+delim;

    while(str.size() != 0)

    {
        if(str.find(delim) != string::npos)
        {
            in.push_back(str.substr(0, str.find(delim)));
            str = str.substr(str.find(delim)+delim.size(), str.size()); 
        }
    }
}
template <class T, class C>
void writeColAttr(std::string name, std::string path, MultiArrayView<2, T, C> ar, IDS_ type, bool in)
{
    if(ar.size() > 0)
    {
    ArrayVector<std::string> col_type;
    if(in)
        col_type.push_back("RESPONSE");
    else
        col_type.push_back("INPUT");
	writeHDF5attr(name, path+".type", col_type);
	std::set<T> a;
	if(type == CATEGORICAL)
	{
        for(int ii = 0 ; ii < ar.size() ; ++ii)
		    a.insert(ar[ii]);
		ArrayVector<T> cat(a.begin(), a.end());
		writeHDF5attr(name, path+".domain", cat);
    }
	else
	{
		ArrayVector<std::string> cat(1, "NUMERAL");
		writeHDF5attr(name, path+".domain", cat);
	}
    }
}



void processCSVname(std::string in, 
					std::string& name, 
					int & beg, 
					int & end, 
					IDS_ & type,
					bool & in_)
{
	ArrayVector<std::string>  tokens;
	tokenize(in, tokens, "::");

	if(tokens.size() < 4)
		in_ = true;
	else if(tokens[3] == "IN")
		in_ = true;
	else
		in_ = false; 

	if(tokens.size() < 3)
		type = IFSTRING;
	else if(tokens[2] == "ORDINAL" || tokens[2] == "")
		type = ORDINAL;
	else 
		type = CATEGORICAL;

	if(tokens.size() < 2)
		end = -1;
	else
	{
		ArrayVector<std::string>  tokensloc;
		tokenize(tokens[1], tokensloc, ":");
		beg = atoi(tokensloc[0].c_str());
		if(tokensloc.size() < 2)
			end = beg +1;
		else
			end = atoi(tokensloc[1].c_str());
	}
	name = tokens[0];

}



struct DlmInfo
{
	std::string 				csvFile;

	ArrayVector<IDS_>			type_;
    std::string					delim_;

	int							doubleCount_;
	int							stringCount_;
	int							rowCount_;
	
	int							h_count_;
	int							t_count_;
	ArrayVector<std::string>	names_;
	
	int colCount()
	{
		return doubleCount_ + stringCount_;
	}
	DlmInfo(std::string name, std::string delim = ",", 
			int headercount		= 0, 
			int trailercount	= -1)
	:
		csvFile(name),
		delim_(delim),
		doubleCount_(0),
		stringCount_(0),
		rowCount_(0),
		h_count_(headercount),
		t_count_(0)
	{
		string 		line;
		ifstream 	fin(name.c_str(), ios::in);
		int			count = 0;

		//get the first line!
		while(!getline(fin, line).eof() && count <= h_count_)
			++count;

		ArrayVector<string> tokens;
		tokenize(line, tokens, delim);
		ArrayVector<std::string>::iterator iter;
		for(iter = tokens.begin(); iter != tokens.end(); ++iter)
		{
            std::istringstream sin(*iter);
			double in;  
			sin >> in;
			if(sin.fail())
			{
				type_.push_back(LITERAL);
				++stringCount_;
			}
			else
			{
				type_.push_back(NUMERAL);
				++doubleCount_;
			}
		}
		while(!getline(fin, line).eof())
			++count;
        if(line == "")
            t_count_ = count - trailercount +1;
        else
		    t_count_ = count - trailercount;
		rowCount_ = t_count_ - h_count_;
		fin.close();
	}
};

template<class T, class C1, class C2>
void dlmread(DlmInfo & info, 
			 MultiArrayView<2, T, C1> & dblMult, 
			 MultiArrayView<2, std::string, C2> strMult = MultiArrayView<2, std::string, C2>())
{
	if(info.stringCount_ != 0)
	vigra_precondition(dblMult.shape(0) == strMult.shape(0),
					   "row Mismatch");
	vigra_precondition(dblMult.shape(0) == info.rowCount_,
					   "row_Mismatch2");
	vigra_precondition(dblMult.shape(1) == info.doubleCount_,
					   "dbl col mismatch");
	if(info.stringCount_ != 0)
	vigra_precondition(strMult.shape(1) == info.stringCount_,
					   "str col mismatch");


	ifstream 	fin(info.csvFile.c_str(), ios::in);
	int			count = 0;
	std::string line;
	//get the first line!
	while(count < info.h_count_ && !getline(fin, line).eof())
		++count;

	for(int ii = 0; ii < info.rowCount_; ++ii)
	{
		getline(fin, line);
		ArrayVector<std::string> tokens;
		tokenize(line, tokens, info.delim_);
		int str_jj = 0, dbl_jj= 0;
		for(int jj = 0; jj < info.colCount(); ++jj)
		{
			if(info.type_[jj] == NUMERAL)
			{
                std::istringstream sin(tokens[jj]);
                double a;
				sin >> a;
                dblMult(ii, dbl_jj) = a;
				++dbl_jj;
			}
			else
			{
				strMult(ii, str_jj) = tokens[jj];
				++str_jj;
			}
		}
	}
}


void CSV2HDF5(std::string csv, 
              std::string hdf5, 
              std::string path,
              std::vector<string> const & names, 
              int headerCt = 0, 
              int trailerCt = 0,
              std::string delim = ",")
{
	typedef MultiArrayShape<2>::type Shp;
	
	DlmInfo info(csv,delim, headerCt, trailerCt);
	MultiArray<2, double> dblArr(Shp(info.rowCount_, info.doubleCount_));
	MultiArray<2, string> strArr(Shp(info.rowCount_, info.stringCount_), "");
	dlmread(info, dblArr, strArr);
    int beg, end;
	std::string r_name;
	IDS_ type;
	bool in_, first = true;
	if(names.size() != 0)
	{
		for(int ii = 0; ii < int(names.size()); ++ii)
		{
			processCSVname(names[ii], r_name, beg, end, type, in_);
			if(end <= dblArr.shape(1))
			{
				writeToHDF5File(hdf5.c_str(), (path +"/"+ r_name).c_str(), dblArr.subarray(Shp(0, beg),
															   Shp(info.rowCount_, end)));
				if(type == IFSTRING)
					type = ORDINAL;
				if(in_ == INOUT)
					in_ = IN;
				writeColAttr(hdf5, 
							 path+"/"+ r_name, 
							 dblArr.subarray(Shp(0, beg), Shp(info.rowCount_, end)),
							 type, 
							 in_);
			}
			else
			{
				//writeToHDF5File(hdf5.c_str(), (path + r_name).c_str(), strArr.subarray(Shp(0, beg),
				//											   Shp(info.rowCount_, end)));
				if(type == IFSTRING)
					type = CATEGORICAL;
				if( in_ == INOUT && first)
				{
					in_ = OUT;
					first = false;
				}
				else if( in_ == INOUT)
					in_ = IN;

				writeColAttr(hdf5, 
							 path+ r_name, 
							 strArr.subarray(Shp(0, beg -info.doubleCount_), Shp(info.rowCount_, end- info.doubleCount_)),
							 type, 
							 in_);
			}
		}
	}
	else
	{
		writeToHDF5File(hdf5.c_str(), (path+"doubleCols").c_str(), dblArr);
		writeColAttr(hdf5, path + "doubleCols",dblArr, ORDINAL, IN);
		writeToHDF5File(hdf5.c_str(), (path+"stringCols").c_str(), strArr);
		writeColAttr(hdf5, path + "doubleCols",strArr, CATEGORICAL, IN);
	}
}
int main(int argc, char** argv)
{
    TCLAP::CmdLine cmd("Command description message", ' ', "0.9");
    TCLAP::SwitchArg isARFF("","isARFF","input is ARFF file", cmd, false);
    TCLAP::SwitchArg isCSV("","isCSV","input is CSV file", cmd, true);
    
    TCLAP::ValueArg<std::string> input("i","input","inputfile", true,"a.in",    "string");   
    cmd.add(input);
    TCLAP::ValueArg<std::string> output("o","output","outputfile", true,"a.hdf5","string");   
    cmd.add(output);
    TCLAP::ValueArg<std::string> delim("d","delim","csv delimiter", false, ",", "string");   
    cmd.add(delim);
    TCLAP::ValueArg<std::string> src("s","source","source", false, "dataset", "string");   
    cmd.add(src);
    TCLAP::ValueArg<int> header("t","topper","header lines", false, 0, "int");   
    cmd.add(header);
    TCLAP::ValueArg<int> footer("f","footer","footer lines", false, 0, "int");   
    cmd.add(footer);
    TCLAP::MultiArg<std::string> names("n", "names", "column names", false, "string");
    cmd.add(names); 
    cmd.parse( argc, argv );


	if(isARFF.getValue())
		ARFF2HDF5(input.getValue(), output.getValue());
	else if(isCSV.getValue())
		CSV2HDF5(input.getValue(), 
                 output.getValue(),
                 src.getValue(),
                 names.getValue(), 
                 header.getValue(), 
                 footer.getValue(), 
                 delim.getValue());
	return 0;
}
