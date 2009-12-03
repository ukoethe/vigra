#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <stdexcept>
#include <vigra/array_vector.hxx>
#include <tclap/CmdLine.h>
#include "../../include/vigra/hdf5impex.hxx"
#include <set>
using namespace std;
using namespace vigra;

/** various enums used
 */
enum IDS_ {NUMERAL, LITERAL, ORDINAL, NOMINAL, IFSTRING, INPUT, OUTPUT,TAG, DUNNO};
void ARFF2HDF5(std::string arff, std::string hdf5)
{
}

std::string stripstr(std::string in)
{
    std::string out = "";
    for(int ii = 0; ii < int(in.size()); ++ii)
    {
        if(out.size() == 0  && (in[ii] == ' ' || in [ii] == '\t'))
            continue;
        else if(*out.rbegin() == ' ' && (in[ii] == ' ' || in [ii] == '\t'))
            continue;
        else if(in[ii] == '\t' || in[ii] == ' ')
            out = out + ' ';
        else
            out = out + in[ii];
    }
    return out;
}

/** tokenize a string given delimiter string
 */
void tokenize(std::string  str, ArrayVector<std::string> & in, std::string delim = ",")
{
    if(delim == "w")
    {
        str = stripstr(str);
        delim = " ";
    }
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

bool isInteger(std::string a)
{
    return false;
}
bool isInteger(double a)
{
    return double(int(a)) == a;
}
/** write the attributes of a dataset column
 */
template <class T, class C>
void writeColAttr(std::string filename, 
                  std::string path, 
                  MultiArrayView<2, T, C> ar, 
                  IDS_ type, 
                  IDS_ in, 
                  int number)
{
    //write precedence number
    ArrayVector<int> order(size_t(1), number);
    writeHDF5Attr(filename, path + ".order" , order);

    //write column type
    ArrayVector<std::string> col_type;
    if(in == INPUT)
        col_type.push_back("INPUT");
    else if (in == OUTPUT)
        col_type.push_back("OUTPUT");
    else
    {
		ArrayVector<std::string> tmp(1, "na");
        writeHDF5Attr(filename, path+".scale", tmp);
        writeHDF5Attr(filename, path+".domain", tmp);
        col_type.push_back("TAGS");
        writeHDF5Attr(filename, path+".type", col_type);
        return;
    }
    writeHDF5Attr(filename, path+".type", col_type);

    //write domain and scale of the columns
	if(type == NOMINAL)
	{
		ArrayVector<std::string> scale(1, "NOMINAL");
        writeHDF5Attr(filename, path+".scale", scale);

	    std::set<T> a;
        for(int ii = 0 ; ii < ar.size() ; ++ii)
		    a.insert(ar[ii]);
		ArrayVector<T> cat(a.begin(), a.end());
        writeHDF5Attr(filename, path+".domain", cat);
    }
	else
	{
		ArrayVector<std::string> scale(1, "ORDINAL");
        writeHDF5Attr(filename, path+".scale", scale);

        bool isInt = true;
        for(int ii = 0; ii < ar.size(); ++ii)
            isInt = isInt && isInteger(ar[ii]);
        ArrayVector<std::string> domain;
        if(isInt)
            domain.push_back("INTEGER");
        else
            domain.push_back("REAL");
        writeHDF5Attr(filename, path + ".domain", domain);
	}
}


/** processes the name string delivered by cmd
 */
void processCSVname(std::string in, 
					std::string& name,
                    int last_col_index,
					int & beg, 
					int & end, 
					IDS_ & type,
					IDS_ & in_)
{
	ArrayVector<std::string>  tokens;
	tokenize(in, tokens, "::");
    
    name = tokens[0];
    type = DUNNO;
    in_  = DUNNO;
    beg  = last_col_index + 1;
    end  = beg ;

    for(int ii = 1; ii < int(tokens.size()); ++ii)
    {
        if      (tokens[ii] == "IN")
            in_ = INPUT;
        else if (tokens[ii] == "OUT")
            in_ = OUTPUT;
        else if (tokens[ii] == "TAG")
            in_ = TAG;
        else if (tokens[ii] == "NOMINAL")
            type = NOMINAL;
        else if (tokens[ii] == "ORDINAL")
            type = ORDINAL;
        else if (tokens[ii] == "all")
            beg = -1;
        else
        { 
            double tmp;
            std::istringstream sin(tokens[ii]);
            sin >> tmp;
            if(sin.fail() && tokens[ii].find(":") == string::npos)
                throw std::runtime_error("could not parse name");
		    ArrayVector<std::string>  tokensloc;
		    tokenize(tokens[ii], tokensloc, ":");
	    	beg = atoi(tokensloc[0].c_str());
		    if(tokensloc.size() < 2)
			    end = beg;
		    else
			    end = atoi(tokensloc[1].c_str());
        }
    }
}



struct DlmInfo
{
	std::string 				csvFile;

	ArrayVector<IDS_>			type_;
    ArrayVector<int>            offsets_;
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
                offsets_.push_back(stringCount_);
				++stringCount_;
			}
			else
			{
				type_.push_back(NUMERAL);
                offsets_.push_back(doubleCount_);
				++doubleCount_;
			}
		}
        bool found_blank_line = false;
		while(!getline(fin, line).eof())
        {
            if(stripstr(line) == "")
            {
                found_blank_line = true;
                break;
            }
			++count;
        }
        if(found_blank_line)
            t_count_ = count;
        else
		    t_count_ = count - trailercount;
		rowCount_ = t_count_ - h_count_;
		fin.close();
	}
};

/**read a delimmited table into a double matrix and an optional
 * string matrix
 */
template<class T, class C1, class C2>
void dlmread(DlmInfo & info, 
			 MultiArrayView<2, T, C1> & dblMult, 
			 MultiArrayView<2, std::string, C2> strMult 
                        = MultiArrayView<2, std::string, C2>())
{
    //precondtion checking
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

	//get information for the seperated file and read
    std::cerr << "getting file info" <<std::endl;
	DlmInfo info(csv,delim, headerCt, trailerCt);
	MultiArray<2, double> dblArr(Shp(info.rowCount_, info.doubleCount_));
	MultiArray<2, string> strArr(Shp(info.rowCount_, info.stringCount_), "");
    std::cerr << "loading text file "<<info.doubleCount_ <<" "<< info.rowCount_<< " " <<  info.stringCount_ << std::endl;
	dlmread(info, dblArr, strArr);



    int beg = 0, end = 0, last_end = -1, last_beg = 0;
	std::string r_name;
	IDS_ type, in_;
    bool first_string_column = true;

    //user has supplied column data.
	if(names.size() != 0)
	{
		for(int ii = 0; ii < int(names.size()); ++ii)
		{
            std::cerr << "writing " << names[ii] << std::endl;
			processCSVname(names[ii], r_name, last_end, beg, end, type, in_);
            if(beg == -1)
            {
                beg = 0;
                end = dblArr.shape(1);
            }
            else
            {
                last_end = end;
                last_beg = beg; 
                end = end - beg;
                beg = info.offsets_[beg];
                end = beg + end + 1;
            }
            std::string col_path = path + "/" + r_name;
			if(info.type_[last_beg] == NUMERAL)
			{
				writeToHDF5File(hdf5.c_str(), 
                                col_path.c_str(), 
                                dblArr.subarray(Shp(0, beg),
									            Shp(info.rowCount_, end)));
				if(type == DUNNO)
					type = ORDINAL;
				if(in_ == DUNNO)
					in_ = INPUT;
				writeColAttr(hdf5.c_str(), 
							 col_path.c_str(),
							 dblArr.subarray(Shp(0, beg), 
                                             Shp(info.rowCount_, end)),
							 type, 
							 in_,
                             ii);
			}
			else
			{
                if(end - beg == 1)
				    writeToHDF5File(hdf5.c_str(), 
                                    col_path.c_str(), 
                                    strArr.subarray(Shp(0, beg),
								    				Shp(info.rowCount_, end)).bindOuter(0));
                else
				    writeToHDF5File(hdf5.c_str(), 
                                    col_path.c_str(), 
                                    strArr.subarray(Shp(0, beg),
							    					Shp(info.rowCount_, end)));
				if(type == DUNNO)
					type = NOMINAL;
				if( in_ == DUNNO && first_string_column)
				{
					in_ = OUTPUT;
					first_string_column = false;
				}
				else if(in_ == DUNNO)
					in_ = TAG;

				writeColAttr(hdf5.c_str(), 
							 col_path.c_str(), 
							 strArr.subarray(Shp(0, beg), 
                                             Shp(info.rowCount_, end)),
							 type, 
							 in_,
                             ii);
			}
        std::cerr << "done " << names[ii] << std::endl;
		}
	}
	else
	{
        if(dblArr.size() > 0)
        {
		writeToHDF5File(hdf5.c_str(), 
                        (path+"/doubleCols").c_str(), 
                        dblArr);
		writeColAttr(hdf5, 
                     path + "/doubleCols",
                     dblArr, 
                     ORDINAL, 
                     INPUT,
                     0);
        }
        if(strArr.size() > 0)
        {
		writeToHDF5File(hdf5.c_str(), 
                        (path+"/stringCols").c_str(), 
                        strArr);
		writeColAttr(hdf5, 
                     path + "/stringCols",
                     strArr, 
                     DUNNO, 
                     TAG,
                     1);
        }
	}
}


void INFO2HDF5(std::string info, 
              std::string hdf5, 
              std::string name)
{

	ifstream 	fin(info.c_str(), ios::in);
	std::string line;
    MultiArray<1, std::string> text_file(MultiArrayShape<1>::type(1), "");

	//get the first line!
	while(!getline(fin, line).eof())
		text_file[0] = text_file[0] +"\n" + line;

    name = "_" + name;
    writeToHDF5File(hdf5.c_str(), name.c_str(),  text_file);
}
int main(int argc, char** argv)
{
    TCLAP::CmdLine cmd("Command description message", ' ', "0.9");
    TCLAP::SwitchArg isARFF("","isARFF","input is ARFF file", cmd, false);
    TCLAP::SwitchArg isCSV("","isCSV","input is CSV file", cmd, false);
    TCLAP::SwitchArg isSparse("","isSparse","input is Sparse csv file", cmd, false);
    TCLAP::SwitchArg isINFO("","isINFO","input is txt information", cmd, false);
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
    bool iscsv = isCSV.getValue();
    if(iscsv == false && !isARFF.getValue() && !isINFO.getValue())
        iscsv = true;
    
    ifstream fin(input.getValue().c_str());
    if(fin.fail())
    {
        throw std::runtime_error("could not find inputfile");
    }
    else
        fin.close();
    std::cerr<< "loading file :" << input.getValue() <<std::endl;
	if(isARFF.getValue())
		ARFF2HDF5(input.getValue(), output.getValue());
    if(isSparse.getValue())
        std::cerr << "Sparse Matrix representation currently not supported" <<std::endl;
	else if(iscsv)
		CSV2HDF5(input.getValue(), 
                 output.getValue(),
                 src.getValue(),
                 names.getValue(), 
                 header.getValue(), 
                 footer.getValue(), 
                 delim.getValue());
    else if(isINFO.getValue())
        INFO2HDF5(input.getValue(),
                  output.getValue(),
                  src.getValue());
	return 0;
}
