/*Public Gamut Mapping Algorithm (PGMA) C Source Code - Version 2.1
Author: Ján Morovic <j.morovic@derby.ac.uk>
Release date: 26 July 2002

pgma_a_v2_1.c - auxiliary functions used in PGMA

2002 © Ján Morovic

See readme_pgma_v2_1.txt for details*/


#include "pgma_v2_1.h"

/*compares two LABs (a and b) and returns 1 if they are the same [at 0.000001 threshold]*/
int eqlab(	struct LAB a,
			struct LAB b)
{
	int eq;

	eq=0;

	if(fabs(a.a-b.a)<0.000001 && fabs(a.b-b.b)<0.000001 && fabs(a.L-b.L)<0.000001)
		eq=1;

	return eq;
}

/*converts XYZ struct to LAB struct
centre		centre of spherical coordiantes
xyz			struct to be converted
returns converted LAB struct*/
struct LAB structxyz2lab(	struct XYZ centre,
							struct XYZ xyz)
{
	struct LAB temp;

	temp.a=xyz.x;
	temp.b=xyz.y;
	temp.L=xyz.z;
	temp.C=sqrt(temp.a*temp.a + temp.b*temp.b);
	spherical(centre,temp.a,temp.b,temp.L,&(temp.alpha),&(temp.theta),&(temp.r));

	temp.type=6;
	return temp;
}

/*converts LAB struct to XYZ struct
lab			struct to be converted
returns converted XYZ struct*/
struct XYZ structlab2xyz(struct LAB lab)
{
	struct XYZ temp;

	temp.x=lab.a;
	temp.y=lab.b;
	temp.z=lab.L;

	return temp;
}

/*Solves augumented matrix using Gaussian elimination method with back-substitution
matrix		matrix containing n x n and n x res parts
n 			n x n part of augumented matrix to be reduced,
res 		number of result vectors*/
int solvegem(	double matrix[MSIZE][MSIZE],
				int n,
				int res)
{
	int test;
	int i,j,k;
	double temp;

	/*test matrix for solvability*/

	for (i=0; i<n; i++)					/*test for zero row in nxn part*/
	{
		test=0;
		for (j=0; j<n; j++)
 			if (matrix[i][j]==0.00) test++;
 		if (test==n)
 		{
 			return 0;
 		}
 	}

 	for (i=0; i<n-1; i++)					/*test for linear dependence in nxn part*/
 	{
		for (j=i+1; j<n; j++)
		{
			test=0;
			for (k=0; k<n; k++)
				if ((matrix[i][k])==(matrix[j][k]*matrix[i][0]/matrix[j][0])) test++;
		}
 		if (test==n)
 		{
 			return 0;
 		}
 	}

	for (i=0; i<n; i++)					/*reduces matrix to echleon form using Gaussian elimination method*/
	{
		k=i+1;
		while (matrix[i][i]==0.0000)
		{
			swapr(matrix,i,k,(n+res));
			k++;
		}

		temp=matrix[i][i];
		for	(j=i; j<n+res; j++)
			matrix[i][j]=(matrix[i][j]/temp);

		for (j=i+1; j<n; j++)
		{
			temp=-matrix[j][i];
			for (k=i; k<n+res; k++)
				matrix[j][k]=matrix[j][k]+temp*matrix[i][k];
		}
	}

	for (i=n-1; i>0; i--)							/*backsubstitution*/
		for (j=i-1; j>-1; j--)
		{
			temp=-matrix[j][i];
			for (k=i; k<n+res; k++)
				matrix[j][k]=matrix[j][k]+temp*matrix[i][k];
		}

	return 1;
}

/*swaps rows a and b in matrix - l is the width of the matrix*/
void swapr(	double matrix[MSIZE][MSIZE],
			int a,
			int b,
			int l)		/*swaps rows a and b in matrix*/
{
	double temp;
	int i;

	for (i=0; i<l; i++)
	{
		temp=matrix[a][i];
		matrix[a][i]=matrix[b][i];
		matrix[b][i]=temp;
	}
}

/*calculates angle of line between [0,0] and [a,b] in degrees ranging from 0 to 360*/
void arctan(	double a,
				double b,
				double *h)
{
	double pi;

	pi=acos(-1.0);

	if (a==0.0 && b==0.0)
		*h=0;
	if (a==0.0)
	{
		if (b>0.0)
			*h=90;
		else if (b<0.0)
			*h=270;
	}
	else
	{
		*h=atan(b/a)*180/pi;
		if (a<0.0)
			*h=*h+180;
		else if (b<0.0)
			*h=*h+360;
	}
	if (b==0.0)
	{
		if (a>0.0)
			*h=0;
		else if (a<0.0)
			*h=180;
	}
}

struct XYAR structlab2xyar(	struct XYAR centrec,
							struct LAB lab)
{
	struct XYAR xy;

	xy.x=lab.C;
	xy.y=lab.L;
	arctan((xy.x-centrec.x),(xy.y-centrec.y),&(xy.alpha));
	xy.r=sqrt((xy.y-centrec.y)*(xy.y-centrec.y) + (xy.x-centrec.x)*(xy.x-centrec.x));

	return xy;
}

