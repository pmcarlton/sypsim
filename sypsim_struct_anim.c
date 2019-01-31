/************************************************************************
 * sypsim_struct.c
 * 2017-12-08 Peter M. Carlton carlton.petermark.3v@kyoto-u.ac.jp
 * 2019-01-21 ....still working on it... 
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * 20190129 -- seems to be working functionally now, but need to find 
 *          parameters that can work or else this model is crap
 *          ...also still no sure things are werkin lit they shoud
 *
 ***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "gd.h"
#include <unistd.h>

#define PPIND 9
#define IMGWIDTH 800
#define IMGSPACING 50
#define CO 0.25
#define STEPNUMERATOR 10.0
#define QUANTUM 200
#define SYPMAX 1000
#define PLKUP 1.00
#define PLKDOWN 0.05
#define PLKLEAVEMOD 100
#define PLKMAX 60
#define PLKRUN 500
#define COPHOSRUN 1000
#define COPHOSEND 2850
#define RUNMAX 14800
#define rnd ( (0.0 + rand() ) / (RAND_MAX + 1.0) )


/************************************************************************
 * Internal data types
 ***********************************************************************/
typedef struct syp_t
{
  float pos;			// -1 for "free", [0..1] for chromosome 1 pos, 
                                //[1..2] for chromosome 2 pos, etc
  int phos;			// 0 or 1 for non-phos or phos
  int id;			// a unique identifier for each molecule
  int gen;			// generation of birth
} syp_t;

/************************************************************************
 * Internal function declarations
 ***********************************************************************/

void print_usage (const char *program_name);
int syp_init (syp_t * syps, int nsyps_current, int nsyps_new);
void print_syps (syp_t * syps, int nsyps_current);
void plk_step (syp_t * syps, int nsyps_current, int* plk);
void syp_step (syp_t * syps, int nsyps_current, int* plk);
float newpos (float pos, int phos, int* plk);
int newphos (float pos, int phos, int* plk);
void ppmat_init (float *ppmat);
void writeimage(syp_t *syps, int nsyps_current); 
void writekymorow(syp_t *syps, int nsyps_current, gdImagePtr kymo, int* plk); 
void writekymopng(gdImagePtr kymo); 
void print_ppmat (); 

/************************************************************************
 * Global vars
 ***********************************************************************/
float ppmat[9];
int run;

/************************************************************************
 * Main function
 ***********************************************************************/
int
main (int argc, char **argv)
{
  int nsyps_current = 0;
  int initsyps = 300;
  int syp_add = 1;
  int li, li2;
  int plk[QUANTUM];
  syp_t syps[SYPMAX];
  gdImagePtr kymo;

  srand (time (NULL));

  kymo = gdImageCreate(4*IMGWIDTH+3*IMGSPACING+QUANTUM, RUNMAX);
 
  nsyps_current = syp_init (syps, nsyps_current, initsyps);

//  ppmat_init (ppmat); - the pos/phos matrix
//  these are constants that determine probability at any step of :
//  - movement (based on whether subunit is phosphorylated or not, and bound or free)
//  - phosphorylation (based on whether subunit is phosphorylated or not, and bound or free)

  ppmat[0] = .0001;	//phos - free or step
  ppmat[1] = ppmat[0];	//nonphos - free or step
  /*ppmat[1] = .001;	//nonphos - free or step*/
  ppmat[2] = .05;	//free - free or random position

  ppmat[3] = 0.0;        //free phos - dephos if gt
  ppmat[4] = 1;        //free nonphos - phos if gt
  ppmat[5] = 1;         //near-CO phos - dephos if gt
  ppmat[6] = 0;         //near-CO nonphos - phos if gt
  ppmat[7] = .99; //other bound phos - dephos if gt
  ppmat[8] = 1;         //other bound nonphos - phos if gt

for (li=0;li<QUANTUM;li++) {plk[li]=0;} // initialize PLK with zeros

  while (nsyps_current < SYPMAX) { nsyps_current = syp_init (syps, nsyps_current, syp_add); }
  writeimage (syps, nsyps_current);

//  print_ppmat();
  for (run = 0; run < RUNMAX; run++)
    {
        if(run>PLKRUN) {
            plk_step(syps, nsyps_current, plk);
        }
/*        if(run==1) {*/
            /*print_syps(syps, nsyps_current);*/
/*        }*/

      syp_step (syps, nsyps_current, plk);
      /*printf("plk: "); for (li=0;li<QUANTUM;li++) {printf("%i ",plk[li]);} printf("\n");*/
  
  //    writeimage (syps, nsyps_current);
      writekymorow(syps, nsyps_current, kymo, plk);
    }
// print_syps (syps, nsyps_current);
 writeimage (syps, nsyps_current);
 writekymopng(kymo);

 return (0);
}

void
print_usage (const char *program_name)
{
  fprintf (stdout, "Usage: %s ... how to use ... \n", program_name);
}


int
syp_init (syp_t * syps, int nsyps_current, int nsyps_new)
{
  static int gen = 0;
  static int id = 0;
  int li;
  nsyps_new += nsyps_current;
  for (li = nsyps_current; li < nsyps_new; li++)
    {
      syps[li].id = id++;
      syps[li].phos = (rnd < 0.5) ? 0 : 1;
      syps[li].gen = gen;
      syps[li].pos = rnd;
    }
  gen++;
  return (nsyps_new);
}

void print_ppmat () {
printf("________PPMAT:________\n");
printf("phos - free or step: %.5f \n",ppmat[0]);
printf("nonphos - free or step: %.5f \n",ppmat[1]);
printf("free - free or random: %.5f \n",ppmat[2]);
printf("______________________\n");
}

void
print_syps (syp_t * syps, int nsyps_current)
{
  int li;
  for (li = 0; li < nsyps_current; li++)
    {
      printf
	("Syp id %.4i from generation %i is at pos %.3f and phos-state %i.\n",
	 syps[li].id, syps[li].gen, syps[li].pos, syps[li].phos);
    }
}

void
plk_step (syp_t * syps, int nsyps_current, int* plk)
{
  int li, phos, qpos;
  float pos;

  for (li = 0; li < nsyps_current; li++)
    {
      pos = syps[li].pos;
      phos = syps[li].phos;
      if(pos>0) {
          if(phos) {
              qpos=(int)(pos*(QUANTUM));
              /*printf("qpos %.3f ",qpos);*/
          plk[qpos] +=PLKUP;
      }
    }

  }

  for (li=0; li < QUANTUM; li++) {
      if(plk[li]>PLKMAX) {plk[li]=PLKMAX;}
      plk[li]-=(PLKDOWN*rnd);
      if(plk[li]<0) {plk[li]=0;}
  }
}

void
syp_step (syp_t * syps, int nsyps_current, int* plk)
{
  int li, phos;
  float pos;
  for (li = 0; li < nsyps_current; li++)
    {
      pos = syps[li].pos;
      phos = syps[li].phos;
      phos = newphos (pos, phos, plk);
      pos = newpos (pos, phos, plk);
      syps[li].pos = pos;
      syps[li].phos = phos;
    }
}

float
newpos (float pos, int phos, int* plk)
{
  int plk_level;
  int posindex;
  const float stepsize=(STEPNUMERATOR/QUANTUM);
  float truestepsize;
  int qpos;
  float pos_orig;
  static int pflag=0;

  pos_orig=pos;

  if (pos < 0)
    {
      posindex = 2;
      pos = ((rnd < ppmat[posindex]) ? (-1) : rnd);
//      if(pos<CO & pos>0) {printf("%.2f %.2f %i __ ",CO, pos,run);}
      return (pos);
    }
  else { posindex = phos ? 0 : 1; }

  qpos=(int)(pos*(QUANTUM));
  plk_level=plk[qpos];
  if(plk_level < 0) {printf("jacked!%.3f %i %i",pos,qpos,plk_level);}

  truestepsize=(stepsize/(1+(.01*plk_level)));

  if(truestepsize>stepsize) {printf("2. jacked since truestepsize is %.3f \n",truestepsize);}

  /*pos = (((rnd) < ppmat[posindex]) ? (-1) : (pos + ((rnd * truestepsize) - (truestepsize/2))));
   * this was bad -- since plk level reduced step size but had no effect on 
   * whether a SYP would come off the chromosome or not. Hard to do using the 
   * ppmat[] method -- 20190129pmc*/

  pos = (((rnd) < ppmat[posindex]) ? ((rnd)<(1/(PLKLEAVEMOD+plk_level))?-1:pos) : (pos + ((rnd * truestepsize) - (truestepsize/2))));

  if(pos==(-1)) {return pos;} /* needed to prevent accumulation solely to the right of CO, 20190129pmc*/
  if ((pos <= 0.0) & (pos > (-stepsize)) ) { pos = pos_orig; }
  if (pos >= 1.0) { pos = pos_orig; }
  if(pos_orig<CO & pos>CO) {pos=pos_orig;}
  if(pos_orig>CO & pos<CO) {pos=pos_orig;}
 return (pos);
}

int
newphos (float pos, int phos, int* plk)
{
  int phosindex;
  float codist;
  int qpos,plk_level;
  const float stepsize=(1.0/QUANTUM);

  if (pos < 0) { 
      phosindex = (phos ? 3 : 4); 
  }
  else {

    qpos=(int)(pos*(QUANTUM));
    if(qpos==QUANTUM) {printf("wakkojakkac!");}
    plk_level=plk[qpos];
    codist=pos-(CO); codist=(codist > 0 ? codist : (-codist));
    if ((run > COPHOSRUN & run<COPHOSEND) & (codist<(2*stepsize))) {
        phosindex = ( phos ? 5 : 6); 
    }
    else { 
        phosindex = (phos? 7 : 8); 
    }
  }
  phos = ((rnd*(1/(1+plk_level)) < ppmat[phosindex]) ? (phos) : (!phos));
  /*return ((rnd)<.5?0:1); [> testing <]*/
  return (phos);
}

void
ppmat_init (float *ppmat)
{
  int li;
  for (li = 0; li < PPIND; li++) { ppmat[li] = rnd; }
}

void writekymorow(syp_t *syps, int nsyps_current, gdImagePtr kymo, int* plk) {
    int li,px;
    static int white,black,green,magenta;
    static int rownum=0;
    static int notAllocated=1;

    if(notAllocated) {
    black=gdImageColorAllocate(kymo,0,0,0);
    white=gdImageColorAllocate(kymo,255,255,255);
    magenta=gdImageColorAllocate(kymo,255,0,255);
    green=gdImageColorAllocate(kymo,0,255,0);
    notAllocated=0;
    }

    for(li=0;li<nsyps_current;li++) {
        if(syps[li].pos>0) {
            px=IMGWIDTH*syps[li].pos;
            if(syps[li].phos) {
                gdImageSetPixel(kymo, px, rownum, green);
                gdImageSetPixel(kymo, px+IMGWIDTH+IMGSPACING, rownum, white);
            } else {
                gdImageSetPixel(kymo, px, rownum, magenta);
                gdImageSetPixel(kymo, px+(2*IMGWIDTH+2*IMGSPACING), rownum, white);
            }
            gdImageSetPixel(kymo, px+(3*IMGWIDTH+3*IMGSPACING), rownum, white);
        }
    }
    for(li=0;li<QUANTUM;li++) {
        px=li+(4*IMGWIDTH+4*IMGSPACING);
        gdImageSetPixel(kymo, px, rownum, (plk[li]>0 ? white : black));
    }
        
    rownum++;
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
