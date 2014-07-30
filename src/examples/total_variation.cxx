/************************************************************************/
/*                                                                      */
/*                 Copyright 2012 by Frank Lenzen &                     */
/*                                           Ullrich Koethe             */
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

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <vigra/multi_array.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/impex.hxx>
#include <vigra/tv_filter.hxx>

#ifdef _MSC_VER
#define strcasecmp _stricmp
#endif

using namespace vigra; 


// Read a string from file <file>.
void read_string(FILE *file,const char *name, char* out){
    char dummy[80];
    int s=fscanf(file,"%s %s ",dummy,out);
    if (s==EOF || s<2){
        std::cout<<"Could not read from file.\n";
        exit(0); 
    }
    else if (strcasecmp(name,dummy)!=0){
        std::cout<<"Found parameter "<<dummy<<" instead of "<<name<<"\n";
        exit(0);
    }
    else
    std::cout<<"Parameter "<<name<<" = "<<out<<std::endl;
}

// read a flota-value from file <file>
float read_value(FILE *file,const char *name){
    char dummy[80];
    char dummy2[80];
    float f=0;

    int s=fscanf(file,"%s %s ",dummy,dummy2); 


    if (s==EOF || s<2){
        std::cout<<"Could not read from file.\n";
        exit(0); 
    }
    else if (strcasecmp(name,dummy)!=0){
        std::cout<<"Found parameter "<<dummy<<" instead of "<<name<<"\n";
        exit(0);
    }
    else{
        f=atof(dummy2);
        std::cout<<"Parameter "<<name<<" = "<<f<<std::endl;
    }
    return f;
}


int main(int argc, char ** argv)
{

    using namespace multi_math;

    if(argc <2)
    {
        std::cout << "Usage: " << argv[0] << " parameterfile" << std::endl;
        std::cout << "(supported formats for images: " << impexListFormats() << ")" << std::endl;
        
        return 1;
    }

    try
    {   
        char infile[80],outfile[80];
        double alpha0=0.01*sqrt(255.0),beta0=0.1*sqrt(255.0),sigma=0.1,rho=.5,K=10,eps=0,gamma0=0;
        int mode,inner_steps=200,outer_steps=5,write_steps=0;
        
        FILE *file=fopen(argv[1],"r");         // open parameter file 
        if (!file){
            std::cout<<"Cannot open file "<<argv[1]<<std::endl;
            return 1;
        }
        read_string(file,"input_file",infile);    // read name of input file (an image)
        read_string(file,"output_file",outfile);  // read name of output file (an image)
        mode=(int)read_value(file,"mode");        // read mode:  1- standard Total Variation
        //             2- Anisotropic Total Variation
        //             3- Anisotropic Total Variation with Higher Order term
        
        // further parameters depend on mode					      
        switch(mode){
        case 1:
        case 2:
            inner_steps=1000;                 // read number of iteration steps
            alpha0=read_value(file,"alpha");  // read alpha
            eps=read_value(file,"epsilon");   // read eps
            break;
        case 3:
            outer_steps=(int)read_value(file,"outer_steps"); // read number of outer iteration steps
            inner_steps=(int)read_value(file,"inner_steps"); // read number of inner iteration steps
            alpha0=read_value(file,"alpha");  // read alpha
            beta0=read_value(file,"beta");    // read beta 
            sigma=read_value(file,"sigma");   // read sigma
            rho=read_value(file,"rho");       // read rho
            K=read_value(file,"edge_factor"); // read parameter K (edge sensitivity)
            write_steps=(int)read_value(file,"write_outer_steps"); //read flag for writing intermediate result after each outer iteration
            break;  
            
        case 4:
            outer_steps=(int)read_value(file,"outer_steps"); // read number of outer iteration steps
            inner_steps=(int)read_value(file,"inner_steps"); // read number of inner iteration steps
            alpha0=read_value(file,"alpha"); // read alpha
            beta0=read_value(file,"beta");   // read beta 
            gamma0=read_value(file,"gamma"); // read gamma
            sigma=read_value(file,"sigma");  // read sigma
            rho=read_value(file,"rho");      // read rho
            K=read_value(file,"edge_factor");// read parameter K (edge sensitivity)
            write_steps=(int)read_value(file,"write_outer_steps");//read flag for writing intermediate result after each outer iteration
            break;
        default:
            std::cout<<"Unknown mode "<<mode<<std::endl;
        }
        int stretch=(int)read_value(file,"histogram_stretch");        // flag for histogram stretching: 0 -no, 1 -yes
        
        
        
        ImageImportInfo info(infile); // get information about input image
        
        
        
        if(info.isGrayscale())
        {
            MultiArray<2,double> data(info.shape());
            MultiArray<2,double> out(info.shape()); 
            MultiArray<2,double> weight(info.shape());  
            
            
            for (int y=0;y<data.shape(1);y++){
                for (int x=0;x<data.shape(0);x++){
                    weight(x,y)=1;                         // set weight to 1 (needed for anisotropicTotalVariationFilter and
                    // higherOrderTotalVariationFilter
                }
            }
            importImage(info, destImage(data));        // read image data
            data=data*(1/255.);                        // scale to [0,1] (smoothing parameters are tuned to this scaling)
            
            //------------------------------------
            //     Choose specific TV Filter
            //------------------------------------
            switch(mode){
            case 1:
                std::cout<<"Standard TV filter"<<std::endl;
                totalVariationFilter(data,out,alpha0,inner_steps,eps);
                break;
            case 2: 
                std::cout<<"Weighted TV filter"<<std::endl;
                totalVariationFilter(data,out,weight,alpha0,inner_steps,eps);
                break;
            case 3:{
                    std::cout<<"Anisotropic TV filter"<<std::endl;
                    
                    MultiArray<2,double> phi(info.shape());
                    MultiArray<2,double> alpha(info.shape());
                    MultiArray<2,double> beta(info.shape()); 
                    
                    out=data;   //'data' serves as initial value
                    
                    for (int i=0;i<outer_steps;i++){   // Outer loop to update anisotropic data 
                        std::cout<<"outer step "<<i<<"\n";
                        
                        getAnisotropy(out,phi,alpha,beta,alpha0,beta0,sigma,rho,K);  // get anisotropic data
                        anisotropicTotalVariationFilter(data,weight,phi,alpha,beta,out,inner_steps); //perform smoothing 
                        
                        if(write_steps){
                            char dummy[80];
                            sprintf(dummy,"output_step_%03d.png",i);
                            std::cout<<"Writing temp file\n";
                            if (stretch)
                            exportImage(out, ImageExportInfo(dummy)); //exportImage performes histogramm stretching
                            else{
                                MultiArray<2,unsigned char>  buffer(info.shape());
                                buffer=max(min(out*255+0.5,0.),255.);                              //scaling back to [0,255], clipping, cast to uint8
                                exportImage(buffer, ImageExportInfo(dummy));
                            }  
                        }
                    }
                    break;
                }
                
            case 4:{
                    std::cout<<"Second order TV filter"<<std::endl;
                    
                    MultiArray<2,double> phi(info.shape());
                    MultiArray<2,double> alpha(info.shape());
                    MultiArray<2,double> beta(info.shape());  
                    MultiArray<2,double> gamma(info.shape()); 
                    MultiArray<2,double> xedges(info.shape()); 
                    MultiArray<2,double> yedges(info.shape()); 
                    
                    for (int y=0;y<data.shape(1);y++){  // set data needed for second order term
                        for (int x=0;x<data.shape(0);x++){
                            gamma(x,y)=gamma0;
                            xedges(x,y)=1;
                            yedges(x,y)=1;
                        }
                    }
                    
                    out=data;  //'data' serves as initial value
                    
                    for (int i=0;i<outer_steps;i++){// Outer loop to update anisotropic data 
                        std::cout<<"outer step "<<i<<"\n";
                        
                        getAnisotropy(out,phi,alpha,beta,alpha0,beta0,sigma,rho,K); // get anisotropic data
                        secondOrderTotalVariationFilter(data,weight,phi,alpha,beta,gamma,xedges,yedges,out,inner_steps);//perform smoothing 
                        
                        if(write_steps){
                            char dummy[80];
                            sprintf(dummy,"output_step_%03d.png",i);
                            std::cout<<"Writing temp file\n"; 
                            if (stretch)
                            exportImage(out, ImageExportInfo(dummy));
                            else{
                                MultiArray<2,unsigned char>  buffer(info.shape());
                                buffer=max(min(out*255+0.5,0.),255.);                              //scaling back to [0,255], clipping, cast to uint8
                                exportImage(buffer, ImageExportInfo(dummy));
                            }  
                        }
                    }
                }
            }
            //-------------------------------------
            if (stretch)
            exportImage(out, ImageExportInfo(outfile));  //write result to file
            else{
                MultiArray<2,unsigned char>  buffer(info.shape());
                buffer=min(max(out*255.+0.5,0.),255.);                              //scaling back to [0,255], clipping, cast to uint8
                exportImage(buffer, ImageExportInfo(outfile));
            }
        }
        else
        {
            std::cout<<"Color images are currently not supported !\n";
        }
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}
