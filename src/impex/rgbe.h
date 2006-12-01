#ifndef _H_RGBE
#define _H_RGBE
/* THIS CODE CARRIES NO GUARANTEE OF USABILITY OR FITNESS FOR ANY PURPOSE.
 * WHILE THE AUTHORS HAVE TRIED TO ENSURE THE PROGRAM WORKS CORRECTLY,
 * IT IS STRICTLY USE AT YOUR OWN RISK.  */

/* utility for reading and writing Ward's rgbe image format.
   See rgbe.txt file for more details.
*/

/* changes by Pablo d'Angelo <pablo.dangelo@web.de>:
 * Added vigra_ prefix to all exported symbols
 */

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int valid;            /* indicate which fields are valid */
  char programtype[16]; /* listed at beginning of file to identify it 
                         * after "#?".  defaults to "RGBE" */ 
  float gamma;          /* image has already been gamma corrected with 
                         * given gamma.  defaults to 1.0 (no correction) */
  float exposure;       /* a value of 1.0 in an image corresponds to
			 * <exposure> watts/steradian/m^2. 
			 * defaults to 1.0 */
} vigra_rgbe_header_info;

/* flags indicating which fields in an rgbe_header_info are valid */
#define VIGRA_RGBE_VALID_PROGRAMTYPE 0x01
#define VIGRA_RGBE_VALID_GAMMA       0x02
#define VIGRA_RGBE_VALID_EXPOSURE    0x04

/* return codes for rgbe routines */
#define VIGRA_RGBE_RETURN_SUCCESS 0
#define VIGRA_RGBE_RETURN_FAILURE -1

/* read or write headers */
/* you may set rgbe_header_info to null if you want to */
int VIGRA_RGBE_WriteHeader(FILE *fp, int width, int height, vigra_rgbe_header_info *info);
int VIGRA_RGBE_ReadHeader(FILE *fp, int *width, int *height, vigra_rgbe_header_info *info);

/* read or write pixels */
/* can read or write pixels in chunks of any size including single pixels*/
int VIGRA_RGBE_WritePixels(FILE *fp, float *data, int numpixels);
int VIGRA_RGBE_ReadPixels(FILE *fp, float *data, int numpixels);

/* read or write run length encoded files */
/* must be called to read or write whole scanlines */
int VIGRA_RGBE_WritePixels_RLE(FILE *fp, float *data, int scanline_width,
			 int num_scanlines);
int VIGRA_RGBE_ReadPixels_RLE(FILE *fp, float *data, int scanline_width,
			int num_scanlines);

int VIGRA_RGBE_ReadPixels_Raw_RLE(FILE *fp, unsigned char *data, int scanline_width,
            int num_scanlines);

#ifdef _CPLUSPLUS
/* define if your compiler understands inline commands */
#define INLINE inline
#else
#define INLINE
#endif

INLINE void VIGRA_float2rgbe(unsigned char rgbe[4], float red, float green, float blue);
INLINE void VIGRA_rgbe2float(float *red, float *green, float *blue, unsigned char rgbe[4]);

#ifdef __cplusplus
}
#endif

#endif /* _H_RGBE */



