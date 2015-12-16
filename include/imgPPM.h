#pragma warning( disable: 4996 )


/* Modes this code works for (returned in *mode for ReadPPM(), values to    */
/*    pass into mode in WritePPM()).  Note that PBMs are handled in a very  */
/*    primitive way.                                                        */
#define PPM_ASCII    3
#define PPM_RAW      6
#define PGM_ASCII    2
#define PGM_RAW      5
#define PBM_ASCII    1
#define PBM_RAW      4

/* info about types of images (Max # of types above) */
#define PPM_MAX      7

/* define return codes for WritePPM() */
#ifndef GFXIO_ERRORS
#define GFXIO_ERRORS
    #define GFXIO_OK            0
    #define GFXIO_OPENERROR     1
    #define GFXIO_BADFILE       2
    #define GFXIO_UNSUPPORTED   3
#endif

/* Reads the PPM/PGM/PBM from the file 'f'                                  */
/*    Returns:  A pointer to a character array storing the image.  Data is  */
/*              laid out how it is arranged in the PPM file.  1st byte is   */
/*              the red component of the upper left pixel, 2nd byte is the  */
/*              green component of the ul-pixel, then the blue component,   */
/*              the next 3 bytes contain the rgb values of the pixel to the */
/*              right, and pixels are stored in scan-line order.            */
/*    The value stored in *mode is one of the modes above (e.g., PPM_ASCII) */
/*    The values stored in *w and *h are the image width & height           */             
unsigned char *ReadPPM( char *f, int *mode, int *width, int *height, bool invertY=false );

/* Writes a PPM/PGM/PBM to the file 'f'                                     */
/*    Returns:  One of the error codes from above or GFXIO_OK               */
/*    Input:  'f', the filename to write to                                 */
/*            mode, the way to store the data (e.g., PPM_ASCII or PGM_RAW)  */
/*            width, the image width                                        */
/*            height, the image height                                      */
/*            ptr, a pointer to an unsigned character / unsigned byte array */
/*                 NOTE:  ptr should have length 3*w*h no matter what       */
/*                 format the data is to be written as (even b&w images)!   */
int WritePPM( char *f, int mode, int width, int height, unsigned char *ptr );




