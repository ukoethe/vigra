/************************************************************************/
/*                                                                      */
/*               Copyright 2002 by Ullrich Koethe                       */
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

#include <fstream>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include "vigra/config.hxx"
#include "vigra/sized_int.hxx"
#include "error.hxx"
#include "byteorder.hxx"
#include "void_vector.hxx"
#include "gif.hxx"

#define BitSet(byte,bit)  (((byte) & (bit)) == (bit))

namespace vigra {

namespace {

    int read_data_block(std::ifstream & stream, void_vector<UInt8> & data)
    {
        int count;

        count = stream.get();
        if(!stream.good())
            return -1;
        if(count == 0)
            return 0;
        data.resize(count);
        stream.read( reinterpret_cast< char * >(data.begin()), count);
        if(!stream.good())
            return -1;
        return count;
    }

    struct ColorCluster
    {
        UInt8 cmin[3], cmax[3];
        std::vector<UInt8 *> entries;
        mutable int largest_dim, largest_diff;

        typedef UInt8 rgb[3];

        struct ColorSorter
        {
            int dim;

            ColorSorter(int d)
            : dim(d)
            {}

            bool operator()(UInt8 * l, UInt8 * r) const
            {
                return l[dim] < r[dim];
            }
        };

        ColorCluster()
        : largest_dim(-1)
        {
            reset_minmax();
        }

        void add(UInt8 * entry)
        {
            entries.push_back(entry);

            update_minmax(entry);
        }

        void reset_minmax()
        {
            for(int i=0; i<3; ++i)
            {
                cmin[i] = 255;
                cmax[i] = 0;
            }

            largest_dim = -1;
        }

        void update_minmax(UInt8 * entry)
        {
            for(int i=0; i<3; ++i)
            {
                if(entry[i] < cmin[i])
                    cmin[i] = entry[i];
                if(cmax[i] < entry[i])
                    cmax[i] = entry[i];
            }

            largest_dim = -1;
        }

        void update_largest() const
        {
            if(largest_dim >= 0)
                return;
            largest_diff = cmax[0] - cmin[0];
            largest_dim = 0;

            for(int i=1; i<3; ++i)
            {
                if(largest_diff < cmax[i] - cmin[i])
                {
                    largest_dim = i;
                    largest_diff = cmax[i] - cmin[i];
                }
            }
        }

        bool operator<(ColorCluster const & o) const
        {
            update_largest();
            o.update_largest();
            return largest_diff < o.largest_diff;
        }

        void split(ColorCluster & o)
        {
            update_largest();
            std::sort(entries.begin(), entries.end(), ColorSorter(largest_dim));

            std::vector<UInt8 *> old_list;
            old_list.swap(entries);
            reset_minmax();

            UInt32 i = 0;
            for(; i<old_list.size()/2; ++i)
            {
                add(old_list[i]);
            }
            for(; i<old_list.size(); ++i)
            {
                o.add(old_list[i]);
            }
        }

        void average(UInt8 * color) const
        {
            UIntBiggest r = 0, g = 0, b = 0;

            for(size_t i=0; i<entries.size(); ++i)
            {
                r += entries[i][0];
                g += entries[i][1];
                b += entries[i][2];
            }

            color[0] = (UInt8)(r / entries.size());
            color[1] = (UInt8)(g / entries.size());
            color[2] = (UInt8)(b / entries.size());
        }

        size_t size() const
        {
            return entries.size();
        }
    };


    void find_color_clusters(void_vector<UInt8> & data,
            std::vector<ColorCluster> & clusters, void_vector<UInt8> & colors)
    {
        size_t count = clusters.size();
        size_t size = data.size() / 3;
        size_t i, current;
        for(i=0; i<size; ++i)
        {
            clusters[0].add(data.begin()+3*i);
        }

        for(current = 1; current < count; ++current)
        {
            size_t largest_index = 0;
            for(i=1; i<current; ++i)
            {
                if(clusters[largest_index] < clusters[i])
                {
                    largest_index = i;
                }
            }
            if(clusters[largest_index].size() == 1)
                break;
            clusters[largest_index].split(clusters[current]);
        }

        for(i=0; i<count; ++i)
        {
            if(clusters[i].size() == 0)
            {
                colors[3*i] = colors[3*i+1] = colors[3*i+2] = 0;
            }
            else
            {
                clusters[i].average(colors.begin() + 3*i);
            }
        }
    }

    void find_color_indices(void_vector<UInt8> & data,
           std::vector<ColorCluster> & clusters, void_vector<UInt8> & indices)
    {
        size_t count = clusters.size();
        UInt8 * base = data.begin();

        size_t i;
        for(i=0; i<count; ++i)
        {
            for(size_t j=0; j<clusters[i].size(); ++j)
            {
                size_t offset = (clusters[i].entries[j] - base) / 3;
                indices[offset] = static_cast<UInt8>(i);
            }
        }
    }
} // anonymous namespace

    CodecDesc GIFCodecFactory::getCodecDesc() const
    {
        CodecDesc desc;

        // init file type
        desc.fileType = "GIF";

        // init pixel types
        desc.pixelTypes.resize(1);
        desc.pixelTypes[0] = "UINT8";

        // init compression types
        desc.compressionTypes.resize(0);

        // init magic strings
        desc.magicStrings.resize(1);
        desc.magicStrings[0].resize(4);
        desc.magicStrings[0][0] = 'G';
        desc.magicStrings[0][1] = 'I';
        desc.magicStrings[0][2] = 'F';
        desc.magicStrings[0][3] = '8';

        // init file extensions
        desc.fileExtensions.resize(1);
        desc.fileExtensions[0] = "gif";

        desc.bandNumbers.resize(2);
        desc.bandNumbers[0] = 1;
        desc.bandNumbers[1] = 3;

        return desc;
    }

    std::auto_ptr<Decoder> GIFCodecFactory::getDecoder() const
    {
        return std::auto_ptr<Decoder>( new GIFDecoder() );
    }

    std::auto_ptr<Encoder> GIFCodecFactory::getEncoder() const
    {
        return std::auto_ptr<Encoder>( new GIFEncoder() );
    }

    struct GIFHeader
    {
        // attributes

        UInt16 width, height, maplength;
        UInt8  bits_per_pixel;
        bool  global_colormap, interlace;

        // methods

        void global_from_stream( std::ifstream & stream, const byteorder & bo );
        bool local_from_stream( std::ifstream & stream, const byteorder & bo );
        void global_to_stream( std::ofstream & stream, const byteorder & bo );
        void local_to_stream( std::ofstream & stream, const byteorder & bo );
    };

    void GIFHeader::global_from_stream( std::ifstream & stream, const byteorder & bo )
    {
        UInt8 flag, c, background;
        read_field( stream, bo, width );
        read_field( stream, bo, height );
        read_field( stream, bo, flag );
        read_field( stream, bo, background );
        read_field( stream, bo, c );
        global_colormap = BitSet(flag, 0x80);
        if(global_colormap)
        {
            bits_per_pixel = (flag & 0x07)+1;
            maplength = 3*( 1 << bits_per_pixel);
        }
    }

    void GIFHeader::global_to_stream( std::ofstream & stream, const byteorder & bo )
    {
        write_field( stream, bo, width );
        write_field( stream, bo, height );
        write_field( stream, bo, (UInt8)0xf7 );
        write_field( stream, bo, (UInt8)0 );  // background
        write_field( stream, bo, (UInt8)0 );  // must be zero
    }

    bool GIFHeader::local_from_stream( std::ifstream & stream, const byteorder & bo )
    {
        UInt8 c, flag;
        for ( ; ; )
        {
            c = stream.get();
            if(!stream.good() || c == ';')
                return false;
            if(c == '!')
            {
                void_vector<UInt8> extensions;

                // read and ignore extension data
                read_field( stream, bo, c );
                while (read_data_block(stream, extensions) > 0) /* empty */;
            }
            if(c == ',')
                break;
        }

        UInt16 x,y;

        read_field( stream, bo, x );
        read_field( stream, bo, y );
        read_field( stream, bo, width );
        read_field( stream, bo, height );
        read_field( stream, bo, flag );
        interlace=BitSet(flag,0x40);
        if(BitSet(flag,0x80))
        {
            global_colormap = false;
            bits_per_pixel = (flag & 0x07)+1;
            maplength = 3*( 1 << bits_per_pixel);
        }
        return true;
    }

    void GIFHeader::local_to_stream( std::ofstream & stream, const byteorder & bo )
    {
        write_field( stream, bo, ',' );
        write_field( stream, bo, (UInt16)0 ); // x
        write_field( stream, bo, (UInt16)0 ); // y
        write_field( stream, bo, width );
        write_field( stream, bo, height );
        write_field( stream, bo, (UInt8)0); // use global colormap, no interlace
   }

    struct GIFDecoderImpl
    {
        // attributes

        GIFHeader header;
        std::ifstream stream;
        byteorder bo;
        void_vector< UInt8 > maps, bands;
        UInt32 components;
        UInt8 * scanline;

        // methods

        void decodeGIF();

        // ctor

        GIFDecoderImpl( const std::string & filename );
    };

    GIFDecoderImpl::GIFDecoderImpl( const std::string & filename )
#ifdef VIGRA_NEED_BIN_STREAMS
        : stream( filename.c_str(), std::ios::binary ),
#else
        : stream( filename.c_str() ),
#endif
          bo("little endian"),
          maps(0),
          bands(0),
          scanline(0)
    {
        if(!stream.good())
        {
            std::string msg("Unable to open file '");
            msg += filename;
            msg += "'.";
            vigra_precondition(0, msg.c_str());
        }

        // read the magic number
        char buf[6];
        read_array( stream, bo, buf, 6 );
        std::string magic(6, (std::string::value_type)0);

        std::copy(buf, buf + 6, magic.begin());
        vigra_precondition( magic == "GIF87a" || magic == "GIF89a",
                            "the stored magic number is invalid" );

        // read the header
        header.global_from_stream( stream, bo );

        // read the global color map, if there is one
        if (header.global_colormap)
        {
            // read the maps
            maps.resize(header.maplength);
            read_array( stream, bo, maps.data(), header.maplength );
        }

        if(!header.local_from_stream( stream, bo ))
        {
            std::string msg("Unable to read file '");
            msg += filename;
            msg += "'.";
            vigra_precondition(0, msg.c_str());
        }

        // read the local color map, if there is one
        if (!header.global_colormap)
        {
            // read the maps
            maps.resize(header.maplength);
            read_array( stream, bo, maps.data(), header.maplength );
        }

        // check if image is Gray or RGB
        int i=0;
        components = 1;
        for(; i < header.maplength/3; ++i)
        {
            if(maps[3*i] != maps[3*i+1] || maps[3*i] != maps[3*i+2])
            {
                components = 3;
                break;
            }
        }
    }

    void GIFDecoderImpl::decodeGIF()
    {
        #define MaxStackSize  4096
        #define NullCode  (-1)

        int
            available,
            clear,
            code_mask,
            code_size,
            end_of_information,
            in_code,
            old_code;

        register int
            bits,
            code,
            count;

        register unsigned long
        datum;

        void_vector<Int16> prefix(MaxStackSize);
        void_vector<UInt8> suffix(MaxStackSize);
        void_vector<UInt8> pixel_stack(MaxStackSize+1);
        void_vector<UInt8> packet(256);
        void_vector<UInt16> indices(header.width*header.height);

        register UInt8 *c;
        register UInt16 *p = indices.begin();

        UInt8
        data_size,
        first,
        *top_stack;

        /*
        Initialize GIF data stream decoder.
        */
        data_size = stream.get();
        clear=1 << data_size;
        end_of_information=clear+1;
        available=clear+2;
        old_code=NullCode;
        code_size=data_size+1;
        code_mask=(1 << code_size)-1;
        for (code=0; code < clear; code++)
        {
            prefix[code]=0;
            suffix[code]=code;
        }
        /*
        Decode GIF pixel stream.
        */
        datum=0;
        bits=0;
        c=0;
        count=0;
        first=0;
        top_stack=pixel_stack.begin();
        while (p < indices.end())
        {
            if (top_stack == pixel_stack.begin())
            {
                if (bits < code_size)
                {
                    /*
                      Load bytes until there is enough bits for a code.
                    */
                    if (count == 0)
                    {
                        /*
                          Read a new data block.
                        */
                        count=read_data_block(stream, packet);
                        if (count <= 0)
                          break;
                        c=packet.begin();
                    }
                    datum+=(*c) << bits;
                    bits+=8;
                    c++;
                    count--;
                    continue;
                }
                /*
                  Get the next code.
                */
                code=datum & code_mask;
                datum>>=code_size;
                bits-=code_size;
                /*
                  Interpret the code
                */
                if ((code > available) || (code == end_of_information))
                  break;
                if (code == clear)
                {
                    /*
                      Reset decoder.
                    */
                    code_size=data_size+1;
                    code_mask=(1 << code_size)-1;
                    available=clear+2;
                    old_code=NullCode;
                    continue;
                }
                if (old_code == NullCode)
                {
                    *top_stack++=suffix[code];
                    old_code=code;
                    first=code;
                    continue;
                }
                in_code=code;
                if (code == available)
                {
                    *top_stack++=first;
                    code=old_code;
                }
                while (code > clear)
                {
                  *top_stack++=suffix[code];
                  code=prefix[code];
                }
                first=suffix[code];
                /*
                  Add a new string to the string table,
                */
                if (available >= MaxStackSize)
                  break;
                *top_stack++=first;
                prefix[available]=old_code;
                suffix[available]=first;
                available++;
                if (((available & code_mask) == 0) && (available < MaxStackSize))
                {
                    code_size++;
                    code_mask+=available;
                }
                old_code=in_code;
            }
            /*
              Pop a pixel off the pixel stack.
            */
            top_stack--;
            *p++ =(UInt16) *top_stack;
        }

        // decode intelaced image
        if (header.interlace)
        {
            void_vector<UInt16> non_interlaced(header.width*header.height);

            int pass, x, y;

            register UInt16 *q;

            static int
              interlace_rate[4] = { 8, 8, 4, 2 },
              interlace_start[4] = { 0, 4, 2, 1 };

            p=indices.begin();
            for (pass=0; pass < 4; pass++)
            {
              y=interlace_start[pass];
              while (y < header.height)
              {
                q=non_interlaced.begin()+(y*header.width);
                for (x=0; x < header.width; x++)
                {
                  *q=(*p);
                  p++;
                  q++;
                }
                y+=interlace_rate[pass];
              }
            }

            swap_void_vector( indices, non_interlaced );
            header.interlace = false;
        }

        // apply colormap
        bands.resize(header.width*header.height*components);
        for(int i=0; i<header.width*header.height; ++i)
        {
            if(components == 1)
            {
                bands[i] = maps[3*indices[i]];
            }
            else
            {
                bands[3*i] = maps[3*indices[i]];
                bands[3*i+1] = maps[3*indices[i]+1];
                bands[3*i+2] = maps[3*indices[i]+2];
            }
        }
    }

    void GIFDecoder::init( const std::string & filename )
    {
        pimpl = new GIFDecoderImpl( filename );
    }

    GIFDecoder::~GIFDecoder()
    {
        delete pimpl;
    }

    std::string GIFDecoder::getFileType() const
    {
        return "GIF";
    }

    unsigned int GIFDecoder::getWidth() const
    {
        return pimpl->header.width;
    }

    unsigned int GIFDecoder::getHeight() const
    {
        return pimpl->header.height;
    }

    unsigned int GIFDecoder::getNumBands() const
    {
        return pimpl->components;
    }

    std::string GIFDecoder::getPixelType() const
    {
        return "UINT8";
    }

    unsigned int GIFDecoder::getOffset() const
    {
        return pimpl->components;
    }

    const void * GIFDecoder::currentScanlineOfBand( unsigned int band ) const
    {
        return pimpl->scanline + band;
    }

    void GIFDecoder::nextScanline()
    {
        if (pimpl->scanline)
            pimpl->scanline += getWidth()*getNumBands();
        else
        {
            pimpl->decodeGIF();
            pimpl->scanline = pimpl->bands.begin();
        }
    }

    void GIFDecoder::close() {}
    void GIFDecoder::abort() {}

    struct GIFEncoderImpl
    {
        // attributes

        GIFHeader header;
        std::ofstream stream;
        byteorder bo;
        void_vector< UInt8 > bands;
        void_vector< UInt8 > maps;
        void_vector< UInt8 > indices;
        UInt32 components;
        UInt8 *scanline;
        bool finalized;

        // methods

        void finalize();
        void writeHeader();
        void writeColormap();
        void writeImageData();
        void reduceTo256Colors();
        void outputEncodedData(void_vector< UInt8 > &);

        // ctor

        GIFEncoderImpl( const std::string & filename );
    };

    GIFEncoderImpl::GIFEncoderImpl( const std::string & filename )
#ifdef VIGRA_NEED_BIN_STREAMS
        : stream( filename.c_str(), std::ios::binary ),
#else
        : stream( filename.c_str() ),
#endif
          bo("little endian"),
          bands(0),
          maps(0),
          indices(0),
          scanline(0),
          finalized(false)
    {
        if(!stream.good())
        {
            std::string msg("Unable to open file '");
            msg += filename;
            msg += "'.";
            vigra_precondition(0, msg.c_str());
        }
        // write the magic number
        write_array( stream, bo, "GIF87a", 6 );
    }

    void GIFEncoderImpl::finalize()
    {
        // color depth
        vigra_precondition( components == 1 || components == 3,
                            "number of bands is not supported" );
    }

    void GIFEncoderImpl::writeHeader()
    {
        // write the header
        header.global_to_stream( stream, bo );
        writeColormap();
        header.local_to_stream( stream, bo );
    }

    void GIFEncoderImpl::writeColormap()
    {
        write_array( stream, bo, maps.data(), header.maplength );
    }

    void GIFEncoderImpl::writeImageData()
    {
        stream.put(header.bits_per_pixel); // code size
        if(components == 3)
        {
            outputEncodedData(indices);
        }
        else
        {
            outputEncodedData(bands);
        }
        stream.put(0);   // end of raster stream
        stream.put(';'); // GIF terminator
    }

    void GIFEncoderImpl::reduceTo256Colors()
    {
        header.bits_per_pixel = 8;
        header.maplength = 3*256;

        maps.resize(header.maplength);
        if(components == 3)
        {
            std::vector<ColorCluster> clusters(256);
            find_color_clusters(bands, clusters, maps);
            indices.resize(header.width*header.height);
            find_color_indices(bands, clusters, indices);
        }
        else
        {
            for(int i=0; i<256; ++i)
            {
                maps[3*i] = maps[3*i+1] = maps[3*i+2] = i;
            }
        }
    }

    void GIFEncoderImpl::outputEncodedData(void_vector<UInt8> & indices)
    {
        #define MaxCode(number_bits)  ((1 << (number_bits))-1)
        #define MaxHashTable  5003
        #define MaxGIFBits  12
        #if defined(HasLZW)
        #define MaxGIFTable  (1 << MaxGIFBits)
        #else
        #define MaxGIFTable  max_code
        #endif
        #define GIFOutputCode(code) \
        { \
          /*  \
            Emit a code. \
          */ \
          if (bits > 0) \
            datum|=((long) code << bits); \
          else \
            datum=(long) code; \
          bits+=number_bits; \
          while (bits >= 8) \
          { \
            /*  \
              Add a character to current packet. \
            */ \
            packet[byte_count++]=(UInt8) (datum & 0xff); \
            if (byte_count >= 254) \
              { \
                stream.put(byte_count); \
                stream.write(reinterpret_cast< char * >(packet.begin()),byte_count); \
                byte_count=0; \
              } \
            datum>>=8; \
            bits-=8; \
          } \
          if (free_code > max_code)  \
            { \
              number_bits++; \
              if (number_bits == MaxGIFBits) \
                max_code=MaxGIFTable; \
              else \
                max_code=MaxCode(number_bits); \
            } \
        }

        int
          bits,
          byte_count,
          number_bits,
          data_size = header.bits_per_pixel;
        UInt32  i;

        long
          datum;

        register int k;

        register UInt8 *p;

        void_vector<Int16> hash_code(MaxHashTable);
        void_vector<Int16> hash_prefix(MaxHashTable);
        void_vector<Int16> hash_suffix(MaxHashTable);

        Int16
          clear_code,
          end_of_information_code,
          free_code,
          index,
          max_code,
          waiting_code;

        void_vector<UInt8> packet(256);

        /*
          Initialize GIF encoder.
        */
        number_bits=data_size+1;
        max_code=MaxCode(number_bits);
        clear_code=((Int16) 1 << data_size);
        end_of_information_code=clear_code+1;
        free_code=clear_code+2;
        byte_count=0;
        datum=0;
        bits=0;
        for (i=0; i < MaxHashTable; i++)
          hash_code[i]=0;
        GIFOutputCode(clear_code);
        /*
          Encode pixels.
        */
        p=indices.begin();
        waiting_code=*p;
        for (i=0; i < indices.size(); i++)
        {
          if(i > 0)
          {
            /*
              Probe hash table.
            */
            index=*p & 0xff;
            k=(int) ((int) index << (MaxGIFBits-8))+waiting_code;
            if (k >= MaxHashTable)
              k-=MaxHashTable;
            GIFOutputCode(waiting_code);
            if (free_code < MaxGIFTable)
            {
                hash_code[k]=free_code++;
                hash_prefix[k]=waiting_code;
                hash_suffix[k]=index;
            }
            else
            {
                /*
                  Fill the hash table with empty entries.
                */
                for (k=0; k < MaxHashTable; k++)
                  hash_code[k]=0;
                /*
                  Reset compressor and issue a clear code.
                */
                free_code=clear_code+2;
                GIFOutputCode(clear_code);
                number_bits=data_size+1;
                max_code=MaxCode(number_bits);
            }
            waiting_code=index;
          }
          p++;
        }
        /*
          Flush out the buffered code.
        */
        GIFOutputCode(waiting_code);
        GIFOutputCode(end_of_information_code);
        if (bits > 0)
        {
            /*
              Add a character to current packet.
            */
            packet[byte_count++]=(UInt8) (datum & 0xff);
            if (byte_count >= 254)
            {
                stream.put(byte_count);
                stream.write(reinterpret_cast< char * >(packet.begin()),byte_count);
                byte_count=0;
            }
        }
        /*
          Flush accumulated data.
        */
        if (byte_count > 0)
        {
                stream.put(byte_count);
                stream.write(reinterpret_cast< char * >(packet.begin()),byte_count);
        }
    }

    void GIFEncoder::init( const std::string & filename )
    {
        pimpl = new GIFEncoderImpl(filename);
    }

    GIFEncoder::~GIFEncoder()
    {
        delete pimpl;
    }

    std::string GIFEncoder::getFileType() const
    {
        return "GIF";
    }

    void GIFEncoder::setWidth( unsigned int width )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->header.width = width;
    }

    void GIFEncoder::setHeight( unsigned int height )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->header.height = height;
    }

    void GIFEncoder::setNumBands( unsigned int numBands )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->components = numBands;
    }

    void GIFEncoder::setCompressionType( const std::string & /* comp */, int /* quality */)
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
    }

    void GIFEncoder::setPixelType( const std::string & pixeltype )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        vigra_precondition( pixeltype == "UINT8",
                            "GIFEncoder::setPixelType(): "
                            "GIF raster supports only the UINT8 pixeltype" );
    }

    unsigned int GIFEncoder::getOffset() const
    {
        return pimpl->components;
    }

    void GIFEncoder::finalizeSettings()
    {
        pimpl->finalize();
        pimpl->finalized = true;
    }

    void * GIFEncoder::currentScanlineOfBand( unsigned int band )
    {
        if (!pimpl->scanline)
        {
            pimpl->bands.resize(pimpl->header.width*pimpl->header.height*pimpl->components);
            pimpl->scanline = pimpl->bands.begin();
        }
         return pimpl->scanline + band;
    }

    void GIFEncoder::nextScanline()
    {
        pimpl->scanline += pimpl->header.width*pimpl->components;
    }

    void GIFEncoder::close()
    {
        pimpl->reduceTo256Colors();
        pimpl->writeHeader();
        pimpl->writeImageData();
    }

    void GIFEncoder::abort() {}
}
