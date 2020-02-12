/*Public Gamut Mapping Algorithm (PGMA) C Source Code - Version 2.1
Author: Ján Morovic <j.morovic@derby.ac.uk>
Release date: 26 July 2002

pgma_v2_1.h - header of PGMA package

2002 © Ján Morovic

See readme_pgma_v2_1.txt for details*/

#include <stdio.h>
#include <math.h>


#define SECTORS 16 /*number of divisions in alpha and theta [for GBD]*/
#define INTRITHLD 0.00001  /*difference between area of triangle and sum of subtriangles allowed*/
#define MSIZE 40 /*dimension of matrix*/


struct XYAR
{
	double x;
	double y;
	double alpha;
	double r;
};

struct XY
{
	double x;
	double y;
};


struct LAB
{
	int type;	/*-1 empty, 1 measured, 2 interpolated, 3 modelled, 4 arbitrary gamut point, 5 false gamut point, 6 gpoint C averaged*/
	double L;
	double a;
	double b;
	double C;
	double alpha;
	double theta;
	double r;
};

struct XYZ
{
	double x;
	double y;
	double z;
};

/*LINE defines a line using the parametric form:
x = a1 + t*u1
y = a2 + t*u2
z = a3 + t*u3*/
struct LINE
{
	double a1;
	double a2;
	double a3;
	double u1;
	double u2;
	double u3;
};

/*LINE2D defines a line using the parametric form:
x = a1 + t*u1
y = a2 + t*u2*/
struct LINE2D
{
	double a1;
	double a2;
	double u1;
	double u2;
};

/*PLANE defines a line using the parametric form:
x = b1 + r*v1 + s*w1
y = b2 + r*v2 + s*w2
z = b3 + r*v3 + s*w3  */

struct PLANE
{
	double b1;
	double b2;
	double b3;
	double v1;
	double v2;
	double v3;
	double w1;
	double w2;
	double w3;
};

/*function declarations*/

/*from pgma_gbcalc_v2_1.c*/
struct PLANE hueplane();
void gb_at_hue();
int line_line_isect_2d();
int line_hgb_intersects();
struct LAB hgb_structxyz2lab();
double dist_cl();
void calc_gb();
void interp_gb();
struct LAB findcusp();

/*pgma_a_v2_1.c*/
int eqlab();
struct LAB structxyz2lab();
struct XYZ structlab2xyz();
int solvegem();
void swapr();
void arctan();
struct XYAR structlab2xyar();

/*pgma_a_vgeom_v2_1.c*/
struct PLANE points2plane();
struct LINE points2line();
int testintria();
struct XYZ ppintercept();
void spherical();
void spher2ortho();
struct XYZ cpointll();
struct XYZ poline();
struct XYZ closestpl();
double xyzdist();
void ch2ab();
double triaarea();
struct PLANE orthoplane();
struct XYZ vectorprod();
double vectorsize();

/*pgma_gmas_v2_1.c*/
void init_sgck();
void sgck();
double C_dept_L();
void init_hp_minde();
void  hp_minde();
/**/