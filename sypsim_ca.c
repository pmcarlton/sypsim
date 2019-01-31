/* sypsim_ca -- redoing as completely discrete cellular automaton 20190131pmc
 *
 * 1 2D array that holds:
 *  level 1: SYP-1p
 *  level 2: SYP-1np
 *  level 3: PLK_LEVEL (level of target PLK protein)
 *  level 4: PP1_LEVEL (level of PP1/LAB-1)
 *
 *  under the model that:
 *  1. SYP-1p causes PLK_LEVEL to increase 
 *  2. movement of SYP-1p and SYP-1np is decreased by PLK_LEVEL
 *  3. SYP-1p is increased near a crossover
 *  4. PP1_LEVEL is decreased by PLK_LEVEL
 *  5. PLK_LEVEL is decreased by PP1_LEVEL
 *  6. At each timestep, SYP-1p and SYP-1np move identically:
 *      - probability "off" of being removed from chromosome
 *      - probability "stay" (== (1-off)) of staying on chromosome
 *      - probability stay*"sticky" of staying in the center
 *      - probability (stay-(stay*sticky))/2 of moving left or right
 *  7. PLK_LEVEL reduces "off" and increases "sticky" 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "gd.h"
#include <unistd.h>


//setting up the XS matrix
#define ROWS 4
#define SYPP 0 //row 0 for SYP-1phos
#define SYPN 1 //row 1 for SYP-1nonphos
#define PLKL 2 //row 2 for 'PLK-target' level
#define PP1L 3 //row 3 for PP1 level
#define QUANTUM 200 //number of columns in the xs matrix
#define CO 50 //crossover position; TODO need way for multiple COs
#define PLKMAX 3 //matrix values can be between 0 and this level...start small
#define PP1MAX 3
#define SYPMAX 3

//triggers for things at a certain run number:
#define PLKRUN 500 //when does PLK-2 start working
#define PLKEND 5000
#define COPHOSRUN 1000 //when does SYP-1 phos at CO start
#define COPHOSEND 2850
#define RUNMAX 10000 //when does program end

#define IMGWIDTH 800
#define IMGSPACING 50

#define rnd ( (0.0 + rand() ) / (RAND_MAX + 1.0) )


/************************************************************************
 * Internal function declarations
 ***********************************************************************/

void writeimage(syp_t *syps, int nsyps_current); 
void writekymorow(syp_t *syps, int nsyps_current, gdImagePtr kymo, int xs[ROWS][QUANTUM]); 
void writekymopng(gdImagePtr kymo); 
void print_ppmat (); 
void xs_zero (int xs[ROWS][QUANTUM]); 
void plk_step  (int xs[ROWS][QUANTUM]); 
void pp1_step  (int xs[ROWS][QUANTUM]); 
int  syp_place (int xs[ROWS][QUANTUM]); 
int  syp_step (int xs[ROWS][QUANTUM]); 
/************************************************************************
 * Global vars
 ***********************************************************************/
int run;


/*******************************************
 * other functions
 * *************************/

int syp_place(int* xs) {
    int li;
    int tries=0;
    int phos;
    int successflag=1;
    phos=((rnd<PHOSGEN)?0:1);
    while(successflag | (tries>=MAXPLACETRIES)) {
        li=(int)(rnd*(QUANTUM)); //position of placement -- totally random for now
        if(syp[SYPP][li]+syp[SYPN][li] < SYPMAX) {
        syp[phos][li]++;
        successflag=0;
        }
        tries++;
    }
    return(0==successflag); //returns 1 if successfully placed, 0 if unsuccessful
}
   
void pp1_step(int* xs) {
    int li;
    for (li=0;li<QUANTUM;li++) {
        xs[PP1L][li]++;
        xs[PP1L][li]-=(xs[PLKL][li]);
        if(xs[PP1L][li]>PP1MAX) {
            xs[PP1L][li]=PP1MAX;
        }
        xs[PP1L][li]+=1;
        if(xs[PP1L][li]>PP1MAX) {
            xs[PP1L][li]=PP1MAX;
        }
    }
} 

void plk_step(int* xs) {
    int li;
    for (li=0;li<QUANTUM;li++) {
        xs[PLKL][li]+=xs[SYPP][li];
        if(xs[PLKL][li]>PLKMAX) {
            xs[PLKL][li]=PLKMAX;
        }
        xs[PLKL][li]-=(xs[PP1L][li]>1);
        if(xs[PLKL][li]<0) {
            xs[PLKL][li]=0;
        }
    }
}

/************************************************************************
 * Main function
 ***********************************************************************/
int
main (int argc, char **argv)
{
  int nsyps_current = 0;
  int initsyps = 300;
  int syptotal = 0;
  int syp_add = 1;
  int li, li2;
  int xs[ROWS][QUANTUM]; //the chromosome data
  int xsu[ROWS][QUANTUM];//copy of the data for updating

  gdImagePtr kymo;

  srand (time (NULL));

  kymo = gdImageCreate(4*IMGWIDTH+3*IMGSPACING+QUANTUM, RUNMAX);
 
  xs_zero(xs);  //zero out the two XS-containing matrices
  xs_zero(xsu);

  for (li=0;li<initsyps;li++) {
      syptotal=syptotal + syp_place(xs); //adds 1 if placement successful
  }

  for (run = 0; run < RUNMAX; run++)
    {
        if(run>PLKRUN) {
            plk_step(xs);
        }
        if(run>PP1RUN) {
            pp1_step(xs);
        }

      // go through "free" pool of SYPs and place on the SC 
      // probability of placing should be proportional to concentration

      bound_syps = syp_step(xs); //TODO syp_step, sum up all bound syps@end and return
      free_syps = syptotal-bound_syps;
      free_reg = free_syps;

      prob=KON_BASE*(float)((free_reg+0.0)/(SYPMAX));

      for (li=0;li<free_reg;li++) {
          if(rnd<prob) {
              if(syp_place(xs)) {
                  free_syps--;
              }
          }
      }


      writekymorow(syps, nsyps_current, kymo, plk);
    }
 writeimage (syps, nsyps_current);
 writekymopng(kymo);

 return (0);
}

int syp_step(int xs[ROWS][QUANTUM]) {


}
/* zero out an 'xs' matrix */
void xs_zero (int xs[ROWS][QUANTUM]) {
    int li, li2;
for(li=0;li<ROWS;li++) {
for (li2=0;li2<QUANTUM;li2++) {xs[li][li2]=0;} // initialize chromosome with zeros
}
}

void writekymopng(gdImagePtr kymo) {
  char pngfilename[120];
  FILE *pngout;

  sprintf(pngfilename,"kymo.png");
  pngout = fopen(pngfilename, "wb");

  /* Output the image to the disk file in PNG format. */
  gdImagePng(kymo, pngout);

  /* Close the files. */
  fclose(pngout);

  /* Destroy the image in memory. */
  gdImageDestroy(kymo);

}

void writeimage(syp_t *syps, int nsyps_current) {
    static int imgnum=0;
    int li;
    int px;
  /* Declare the image */
  gdImagePtr im;
  /* Declare output files */
  FILE *pngout;
  /* Declare color indexes */
  int white;
  int black;
  int green;
  int magenta;
  char pngfilename[120];

  imgnum+=1;
  /* Allocate the image: IMGWIDTH pixels across by 64 pixels tall */
  im = gdImageCreate(IMGWIDTH, 64);

  /* Allocate the color black (red, green and blue all minimum).
    Since this is the first color in a new image, it will
    be the background color. */
  black = gdImageColorAllocate(im, 0, 0, 0);
  white = gdImageColorAllocate(im, 255, 255, 255);
  magenta = gdImageColorAllocate(im, 255, 0, 255);
  green = gdImageColorAllocate(im, 0, 255, 0);


for (li=0;li<nsyps_current;li++) {
    if(syps[li].pos > 0) {
    px=(IMGWIDTH*syps[li].pos);
    if(syps[li].phos) {
    gdImageLine(im, px, 0, px, 64, green); //SYP-1-phos is green
    } else {
        gdImageLine(im,px,0,px,64,magenta); //SYP-1-np is magenta
    }
    }
}

  /* Open a file for writing. "wb" means "write binary", important
    under MSDOS, harmless under Unix. */
  sprintf(pngfilename,"output_%.5d.png",imgnum);
  pngout = fopen(pngfilename, "wb");

  /* Output the image to the disk file in PNG format. */
  gdImagePng(im, pngout);

  /* Close the files. */
  fclose(pngout);

  /* Destroy the image in memory. */
  gdImageDestroy(im);
}
 
