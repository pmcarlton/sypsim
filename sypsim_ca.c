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
#include <unistd.h>
#include <tiffio.h>

#define MIN(x,y) ( (x>y) ? (y) : (x) )
#define MAX(x,y) ( (x<y) ? (y) : (x) )

//setting up the XS matrix
#define ROWS 4
#define SYPP 0 //row 0 for SYP-1phos
#define SYPN 1 //row 1 for SYP-1nonphos
#define PLKL 2 //row 2 for 'PLK-target' level
#define PP1L 3 //row 3 for PP1 level
#define QUANTUM 200 //number of columns in the xs matrix
#define CO 30 //crossover position; TODO need way for multiple COs
#define PLKMAX 4 //matrix values can be between 0 and this level...start small
#define PP1MAX 30
#define SYPMAX 10
#define POOLMAX 1000
#define INITSYPS 500
#define KON_BASE 1.0

//triggers for things at a certain run number:
#define TIFFWRITERUN 50
#define PP1RUN 1000
#define PP1END 80000
#define PLKRUN 1000 //when does PLK-2 start working
#define PLKEND 80000
#define PHOSGEN 0.3
#define COPHOSRUN 1000 //when does SYP-1 phos at CO start
#define COPHOSEND 50000
#define RUNMAX 100000 //when does program end
#define MAXPLACETRIES 20 //how many times a SYP tries to get on the SC before giving up

#define rnd ( (0.0 + rand() ) / (RAND_MAX + 1.0) )


/************************************************************************
 * Internal function declarations
 ***********************************************************************/

void writekymorow(int xs[ROWS][QUANTUM]); 
void print_ppmat (); 
void xs_zero (int xs[ROWS][QUANTUM]); 
void plk_step  (int xs[ROWS][QUANTUM]); 
void pp1_step  (int xs[ROWS][QUANTUM]); 
int  syp_place (int xs[ROWS][QUANTUM]); 
int  syp_step (int xs[ROWS][QUANTUM]); 
float get_stickyprob(int plk);
float get_offprob(int plk);
TIFF* open_tiff(); 
void scaleto8bit(int xs[ROWS][QUANTUM], char Image[ROWS*QUANTUM]);

/************************************************************************
 * Global vars
 ***********************************************************************/
int run;


/*******************************************
 * other functions
 * *************************/

int syp_place(int xs[ROWS][QUANTUM]) {
    int li;
    int tries=0;
    int phos;
    int successflag=1;
    phos=((rnd<PHOSGEN)?0:1);
    while(successflag | (tries>=MAXPLACETRIES)) {
        li=(int)(rnd*(QUANTUM)); //position of placement -- totally random for now
        if(xs[SYPP][li]+xs[SYPN][li] < SYPMAX) {
        xs[phos][li]++;
        successflag=0;
        }
        tries++;
    }
    return(0==successflag); //returns 1 if successfully placed, 0 if unsuccessful
}
   
void
pp1_step (int xs[ROWS][QUANTUM])
{
  int li;
  for (li = 0; li < QUANTUM; li++)
    {
      if (rnd < .1)
	{
	  xs[PP1L][li]++;
	  xs[PP1L][li] -= (xs[PLKL][li]);
	  if (xs[PP1L][li] < 0)
	    {
	      xs[PP1L][li] = 0;
	    }
	  if (xs[PP1L][li] > PP1MAX)
	    {
	      xs[PP1L][li] = PP1MAX;
	    }
	}
    }
}

void
plk_step (int xs[ROWS][QUANTUM])
{
  int li;
  for (li = 0; li < QUANTUM; li++)
    {
      if (rnd < .5)
	{
	  xs[PLKL][li] += xs[SYPP][li];
	  if (xs[PLKL][li] > PLKMAX)
	    {
	      xs[PLKL][li] = PLKMAX;
	    }
	  xs[PLKL][li] -= 1;
	  xs[PLKL][li] -= (xs[PP1L][li]>0);
	  if (xs[PLKL][li] < 0)
	    {
	      xs[PLKL][li] = 0;
	    }
	}
    }
}

/************************************************************************
 * Main function
 ***********************************************************************/
int
main (int argc, char **argv)
{
  int initsyps = INITSYPS;
  int syptotal = 0;
  int li, li2, bound_syps, free_syps, free_reg;
  int xs[ROWS][QUANTUM]; //the chromosome data
  int xsu[ROWS][QUANTUM];//copy of the data for updating
  float prob;
  char xsChar[ROWS*QUANTUM]; //for writing to image
  tdata_t buf;
  uint32 row=0;
  TIFF* tif;

  srand (time (NULL));

  
  tif = open_tiff();
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

    if(syptotal<POOLMAX) {
        syptotal=syptotal+syp_place(xs); //syp keeps growing throughout prophase
    }

if (run > COPHOSRUN && run < COPHOSEND)
  {
      for(li=(-1);li<=1;li++) {
          li2=CO+li;
    xs[SYPP][li2] += xs[SYPN][li2];
    xs[SYPN][li2] = 0;
    if (xs[SYPP][li2] > SYPMAX)
      {
	xs[SYPN][li2] = SYPMAX - xs[SYPN][li2];
	xs[SYPP][li2] = SYPMAX;
      }
  }
  }

      bound_syps = syp_step(xs); //TODO syp_step, sum up all bound syps@end and return
      free_syps = syptotal-bound_syps;
      free_reg = free_syps;

      prob=KON_BASE*(float)((free_reg+0.0)/(POOLMAX));

      for (li=0;li<free_reg;li++) {
          if(rnd<prob) {
              if(syp_place(xs)) {
                  free_syps--;
              }
          }
      }

    if(0==(run % (TIFFWRITERUN))) {
    scaleto8bit(xs,xsChar);
    buf=&xsChar;
    TIFFWriteScanline(tif, buf, row, 0);
    row=row+1;
      }
  }

TIFFClose (tif);
return (0);
}

void writekymorow(int xs[ROWS][QUANTUM]) {
int li,li2;
for(li=0;li<ROWS;li++) {
    for(li2=0;li2<QUANTUM;li2++) {
        printf("%i ",xs[li][li2]);
    }
}
printf("\n");
}


int syp_step(int xs[ROWS][QUANTUM]) {
int li,li2,syp,acc,r,i;
int xs2[ROWS][QUANTUM];
float off,stay,sticky,lateral,test;

xs_zero(xs2);

for (li=0;li<QUANTUM;li++) {
    for(r=0;r<2;r++) {
    syp=xs[r][li];
    if(syp) { 
        for (i=0;i<syp;i++) {
            off=get_offprob(xs[PLKL][li]);
            stay=1-off;
            sticky=stay*get_stickyprob(xs[PLKL][li]);
            lateral=(stay-sticky);
            test=rnd;
            if(test<lateral) { // move left or right
                 li2=(rnd<0.5?li+1:li-1); //equal prob. of left or right
                 if((li2==CO) | (li2<0) | (li2==QUANTUM)) {li2=li;} //forbidden moves
            xs2[r][li2]++;
              } 
            else if (test<stay) { xs2[r][li]++; }//stay in center
        }
    }
    }
}

//re-make xs from xs2
for(li=0;li<2;li++) {
    for (li2=0;li2<QUANTUM;li2++) {
        xs[li][li2]=xs2[li][li2];
    }
}

//sum up at the end:
acc=0;
for (li=0;li<QUANTUM;li++) {
acc += xs[SYPP][li];
acc += xs[SYPN][li];
}

return(acc);
}


float get_offprob(int plk) {
    return(0.0001*(1-(plk/PLKMAX)));
//return(1.0/(20.0+plk*8));
}

float get_stickyprob(int plk) {
return(0);
//return(1.0/2*(2+(plk)));
}

/* zero out an 'xs' matrix */
void xs_zero (int xs[ROWS][QUANTUM]) {
    int li, li2;
for(li=0;li<ROWS;li++) {
for (li2=0;li2<QUANTUM;li2++) {xs[li][li2]=0;} // initialize chromosome with zeros
}
}

TIFF* open_tiff() {

        TIFF* tif = TIFFOpen("rpm.tif", "w");
        TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField (tif, TIFFTAG_ARTIST, "Pete Carlton");
	TIFFSetField (tif, TIFFTAG_IMAGELENGTH, RUNMAX/TIFFWRITERUN);
	TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, QUANTUM*ROWS);
	TIFFSetField (tif, TIFFTAG_IMAGEDEPTH, 1);
	TIFFSetField (tif, TIFFTAG_PLANARCONFIG, 1);	
	TIFFSetField (tif, TIFFTAG_PHOTOMETRIC, 1);
        return(tif);
}
	

		           
void scaleto8bit(int xs[ROWS][QUANTUM], char Image[ROWS*QUANTUM])
{
int x,y,li=0;
float scl;
const int rowscl[4]={SYPMAX,SYPMAX,PLKMAX,PP1MAX};
/*double correction;*/
/*max -= min;*/
/*correction = 255 / max;*/

for (x=0;x<ROWS;x++) {
scl=255.0 / rowscl[x];
for (y=0;y<QUANTUM;y++) {
Image[li] = (char)(scl*xs[x][y]); li++;
}
}
Image[0]=(char)255;
Image[1*QUANTUM]=(char)255;
Image[2*QUANTUM]=(char)255;
Image[3*QUANTUM]=(char)255;
Image[4*QUANTUM-1]=(char)255;
}


