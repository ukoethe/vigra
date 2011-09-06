/************************************************************************/
/*                                                                      */
/*       Copyright 2010 by Joachim Schleicher and Ullrich Koethe        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/


/*                                                          
 *  Opens an Andor .sif file as MultiImageView.             
 *  The width, height and number of images are extracted    
 *  from the ASCII encoded variable length header.          
 *                                                          
 *  Based on the Java-Code from                             
 *  http://rsb.info.nih.gov/ij/plugins/open-sif.html        
 *  written by                                              
 *  L. Stirling Churchman (stirling at stanford.edu)        
 *  Philippe Carl (pcarl at uni-muenster.de)                
 *  Yoshiyuki Arai (arai at phys1.med.osaka-u.ac.jp)        
 *                                                          
 *  Currently tested SIF versions: 4.16.12005.0             
 *                                 4.16.30001.0             
 *                                 4. 6.    3.0             
 *                                 4. 6.    0.0             
*/



#include "vigra/sifImport.hxx"
#include "vigra/utilities.hxx"
#include "byteorder.hxx"

namespace vigra {

// private helper functions
namespace helper {
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

class BadConversion : public std::runtime_error {
 public:
   BadConversion(std::string const& s)
     : std::runtime_error(s)
     { }
};
 
inline double convertToDouble(std::string const& s) {
   std::istringstream i(s);
   double x;
   if (!(i >> x))
     throw BadConversion("convertToDouble(\"" + s + "\")");
   return x;
} 

inline int convertToInt(std::string const& s) {
   std::istringstream i(s);
   int x;
   if (!(i >> x))
     throw BadConversion("convertToDouble(\"" + s + "\")");
   return x;
} 

}// namespace helper
 
/********************************************************/
/*                                                      */
/*                   SIFImportInfo                      */
/*                                                      */
/********************************************************/

int SIFImportInfo::width() const {    return m_width; }
int SIFImportInfo::height() const {    return m_height; }
int SIFImportInfo::stacksize() const {    return m_stacksize; }
std::ptrdiff_t SIFImportInfo::getOffset() const {    return m_offset; }
const char * SIFImportInfo::getFileName() const  {    return m_filename;    }    




SIFImportInfo::SIFImportInfo(const char* filename) :
    m_filename(filename), m_offset(0), headerlen(32)
{

    // some initialisations
    left = 1, right = 512, bottom = 1, top = 512;
    xbin = 1, ybin = 1;
    xres = 1, yres = 1;
    
    std::ifstream siffile (filename);
    if( !siffile.is_open() )
    {
        std::string msg("Unable to open file '");
        msg += filename;
        msg += "'. File not found.";
        vigra_precondition(0, msg.c_str());
    }    
    int spool = 0;    //If Data is Spooled one, spool is set to 1
    for(int i=0;i<headerlen+spool+1;i++) {
        std::string str;
        getline(siffile, str);
#ifdef DEBUG
        std::std::cout << str << std::std::endl;
#endif
        if(i==0) {
            vigra_precondition(str=="Andor Technology Multi-Channel File", "The file is not a valid sif File.");
        }
        if(i==2) { // Extract parameters such as temperature, exposureTime and so on. 
            std::vector<std::string> tokens = helper::split(str, ' ');
            // TODO
            //~ d = new Date(Long.parseLong(tokens[4])*1000); // Date is recorded as seconds counted from 1970.1.1 00:00:00
            temperature1 = helper::convertToDouble(tokens[5]);
            exposureTime = tokens[12];
            cycleTime = tokens[13];
            readout = 1/helper::convertToDouble(tokens[18])/1e+6;
            EMGain = tokens[21];
            verticalShiftSpeed = tokens[41];
            preAmpGain = tokens[43];
            temperature2 = helper::convertToDouble(tokens[47]);
            if(tokens.size() > 57) {
                version = tokens[54]+"."+tokens[55]+"."+tokens[56]+"."+tokens[57];
            }
            if(temperature1 == -999) 
                // If the temperature is unstable, temperature1 value is -999 and unstable temperature value is recored in temperature2
                temperature = asString(temperature2) + " (Unstable)";
            else
                temperature = asString(temperature1);
        }
        if(i==3) {
            model = str; // Model of EMCCD camera
        }
        if(i==5) {
            originalFilename = str; // Read original filename
        }
        if(i==7) { // If the Data is spooled one, "spool" value is set to 1
            std::vector<std::string> tokens = helper::split(str, ' ');
            if(tokens.size() >= 1 && tokens[0]=="Spooled") {
                spool=1;
            }
        }
        if(i > 7 && i < headerlen-12) 
        {
            if(str.size() == 17 && 
                 str.substr(0,6)=="65539 " && 
                 str[6] == 0x01 && str[7] == 0x20 && str[8] == 0x00) 
            { // seems to be always ten lines before the dimensions-line 
                 headerlen = i+12;
            }
        }
        if(i==(headerlen-2)) { // Read size of stack (frame length)
            std::string str2 = str.substr ( 0, 12 );
            if(str2 == "Pixel number") str = str.substr(12); // drop "Pixel number" as first letters
            std::vector<std::string> tokens = helper::split(str, ' ');
            vigra_precondition(tokens.size() >= 6, "format error. Not able to read stacksize.");
            yres = helper::convertToInt(tokens[2]);
            xres = helper::convertToInt(tokens[3]);
            m_stacksize = helper::convertToInt(tokens[5]);
            
        }
        if(i==(headerlen-1)) { // Read information about the size and bin
            std::vector<std::string> tokens = helper::split(str, ' ');
            vigra_precondition(tokens.size() >= 7, "format error. Not able to read image dimensions.");
            left = helper::convertToInt(tokens[1]);
            top = helper::convertToInt(tokens[2]);
            right = helper::convertToInt(tokens[3]);
            bottom = helper::convertToInt(tokens[4]);
            xbin = helper::convertToInt(tokens[5]);
            ybin = helper::convertToInt(tokens[6]);
        }
    }
    
    // determine filesize
    siffile.seekg (0, std::ios::end); // goto the end
    filesize = siffile.tellg();
    siffile.seekg (0, std::ios::beg);  // goto the beginning
    filesize -= siffile.tellg();

    // Estimate the offset value (header length)
    for (int i = 0; i < (headerlen+m_stacksize); i++){
        while(siffile.get() != 10) {
            m_offset++;
        }
        m_offset++;
    }
    if(siffile.get() == '0' && siffile.get() == 10) { // Newer sif version
        m_offset += 2;
    }
    siffile.close();
    
    // Calc the width and the height value of the ROI
    m_width = right-left+1;
    mod = m_width % xbin;
    m_width = (m_width-mod)/ybin;
    m_height = top-bottom+1;
    mod = m_height % ybin;
    m_height = (m_height-mod)/xbin;


    size_t data_size = (size_t)m_width * m_height * 4 * m_stacksize;
    vigra_precondition(m_offset + data_size + 8 == filesize, "error reading sif file: data with header should be equal to filesize. ");
    


}

// this function only works for MultiArrayView<3, float> so we don't use a template here.
void readSIF(const SIFImportInfo &info, MultiArrayView<3, float, UnstridedArrayTag> array) 
{
    vigra_precondition(sizeof(float) == 4, "SIF files can only be read into MultiArrayView<float32>. On your machine a float has more than 4 bytes.");
    float * memblock =  array.data();        // here we assume that MultiArray hat float32 values as the sif raw data!!

    std::ifstream file (info.getFileName(), std::ios::in|std::ios::binary);
    vigra_precondition(file.is_open(), "Unable to open sif file");

    byteorder bo = byteorder("little endian");  // SIF file is little-endian

    std::ptrdiff_t pos = file.tellg();        // pointer to beginning of the file
    file.seekg(pos+info.getOffset());
    read_array( file, bo, memblock, info.width()*info.height()*info.stacksize() );
    file.close();

    return;
    
}


std::ostream& operator<<(std::ostream& os, const SIFImportInfo& info)
{
    // output
    os << "\n" <<
        "SIF Image Information: " <<
        "\nOriginal Filename:\t" << info.originalFilename <<
        "\nDate and Time:\t" << info.d <<
        "\nSoftware Version:\t" << info.version <<
        "\nCamera Model:\t\t\t" << info.model <<
        "\nTemperature (C):\t\t"        << info.temperature <<
        "\nExposure Time (s):\t\t"    << info.exposureTime <<
        "\nCycle Time (s):\t\t\t" << info.cycleTime <<
        "\nPixel Readout Rate (MHz):\t" << info.readout <<
        "\nHorizontal Camera Resolution:\t"    << info.xres <<
        "\nVertical Camera Resolution:\t"        << info.yres <<
        "\nImage width:\t\t"          << info.width() <<
        "\nImage Height:\t\t"        << info.height() <<
        "\nHorizontal Binning:\t"    << info.xbin <<
        "\nVertical Binning:\t"        << info.ybin <<
        "\nEM Gain level:\t"        << info.EMGain <<
        "\nVertical Shift Speed (s):\t" << info.verticalShiftSpeed <<
        "\nPre-Amplifier Gain:\t" << info.preAmpGain <<
        "\nStacksize: \t\t\t" << info.m_stacksize <<
        "\nFilesize: \t\t\t" << info.filesize <<
        "\nOffset to Image Data: \t" << info.m_offset <<
        "\n";    
        
    return os;
}

} // namespace vigra
