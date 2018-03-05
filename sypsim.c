/************************************************************************
 * sypsim.c
 * 2017-12-08 Peter M. Carlton pcarlton@icems.kyoto-u.ac.jp
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
 ***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "gd.h"
#include <unistd.h>

#define CO 0.75
#define SYPMAX 1000
#define RUNMAX 100000
#define PPIND 12
#define STEPSIZE 0.00
#define HSTEPSIZE 0.00
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
void syp_step (syp_t * syps, int nsyps_current);
float newpos (float pos, int phos, float co);
int newphos (float pos, int phos, float co);
void ppmat_init (float *ppmat);
int writeimage(syp_t *syps, int nsyps_current); 

/************************************************************************
 * Global vars
 ***********************************************************************/
float ppmat[PPIND];

/************************************************************************
 * Main function
 ***********************************************************************/
int
main (int argc, char **argv)
{
  int nsyps_current = 0;
  int initsyps = 20;
  int syp_add = 5;
  int li, li2, run;
  syp_t syps[SYPMAX];

  srand (time (NULL));

  nsyps_current = syp_init (syps, nsyps_current, initsyps);
//  ppmat_init (ppmat);
//  these are 12 constants that determine probability at any step of :
//  - movement (based on whether subunit is phosphorylated or not, and at long 
//  arm, short arm, or free)
//  - phos-to-nonphos transition (based on whether subunit is phosphorylated or 
//  not, and at long arm, short arm, or free)
  ppmat[0] = .9;		//long arm phos - free or step
  ppmat[1] = .1;		//long arm nonphos - free or step
  ppmat[2] = .1;		//short arm phos - free or step
  ppmat[3] = .9;		//short arm nonphos - free or step
  ppmat[4] = 0.1;		//free phos - free or random position
  ppmat[5] = 0.1;		//free nonphos - free or random position
  ppmat[6] = 1;		//long arm phos - dephosphorylate if gt
  ppmat[7] = 1;		//long arm nonphos - phosphorylate if gt
  ppmat[8] = 1;		//short arm phos - dephosphorylate if gt
  ppmat[9] = 1;		//short arm nonphos - phosphorylate if gt
  ppmat[10] = 1;	//free phos - dephosphorylate if gt
  ppmat[11] = 1;	//free nonphos - phosphorylate if gt

  for (run = 0; run < RUNMAX; run++)
    {
      if ((nsyps_current + syp_add) <= SYPMAX) { 
          nsyps_current = syp_init (syps, nsyps_current, syp_add); }
      syp_step (syps, nsyps_current);
    }
  print_syps (syps, nsyps_current);
  writeimage (syps, nsyps_current);
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
syp_step (syp_t * syps, int nsyps_current)
{
  int li, phos;
  float pos;
  const float co = CO;		//CO position
  for (li = 0; li < nsyps_current; li++)
    {
      pos = syps[li].pos;
      phos = syps[li].phos;
      pos = newpos (pos, phos, co);
      phos = newphos (pos, phos, co);
      syps[li].pos = pos;
      syps[li].phos = phos;
    }
}

float
newpos (float pos, int phos, float co)
{
  float pos_orig;
  int posindex;
  pos_orig=pos;
  if (pos < 0)
    {
      posindex = (phos ? 4 : 5);
      pos = ((rnd < ppmat[posindex]) ? (-1) : rnd);
      return (pos);
    }
  else if (pos > co) { posindex = phos ? 2 : 3; }
  else { posindex = phos ? 0 : 1; }
  pos = ((rnd < ppmat[posindex]) ? (-1) : (pos + ((rnd * STEPSIZE) - (HSTEPSIZE))));
  if (pos < 0) { pos = (-1); }
  if (pos > 1) { pos = (-1); }
  return (pos);
}

int
newphos (float pos, int phos, float co)
{
  int phosindex;
  if (pos < 0) { phosindex = (phos ? 10 : 11); }
  else if (pos > co) { phosindex = (phos ? 8 : 9); }
  else { phosindex = (phos ? 6 : 7); }
  phos = ((rnd < ppmat[phosindex]) ? (phos) : (!phos));
  return (phos);
}

void
ppmat_init (float *ppmat)
{
  int li;
  for (li = 0; li < PPIND; li++) { ppmat[li] = rnd; }
}

int writeimage(syp_t *syps, int nsyps_current) {
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
  int red;

  /* Allocate the image: 640 pixels across by 64 pixels tall */
  im = gdImageCreate(640, 64);

  /* Allocate the color black (red, green and blue all minimum).
    Since this is the first color in a new image, it will
    be the background color. */
  black = gdImageColorAllocate(im, 0, 0, 0);

  /* Allocate the color white (red, green and blue all maximum). */
  white = gdImageColorAllocate(im, 255, 255, 255);
  red = gdImageColorAllocate(im, 255, 0, 0);
  green = gdImageColorAllocate(im, 0, 255, 0);


for (li=0;li<nsyps_current;li++) {
    px=(640*syps[li].pos);
    if(syps[li].phos) {
    gdImageLine(im, px, 0, px, 64, green); //SYP-1-phos is green
    } else {
        gdImageLine(im,px,0,px,64,red); //SYP-1-np is red
    }
}

  /* Open a file for writing. "wb" means "write binary", important
    under MSDOS, harmless under Unix. */
  pngout = fopen("test.png", "wb");

  /* Output the image to the disk file in PNG format. */
  gdImagePng(im, pngout);

  /* Close the files. */
  fclose(pngout);

  /* Destroy the image in memory. */
  gdImageDestroy(im);
}

