/*Public Gamut Mapping Algorithm (PGMA) C Source Code - Version 2.1
Author: Ján Morovic <j.morovic@derby.ac.uk>
Release date: 26 July 2002

pgma_test_v2_1.c - test use of gamut boundary calculation and gamut mapping functions
in pgma_gbcalc_v2_1.c and pgma_gmas_v2_1.c

2002 © Ján Morovic & Pei-Li Sun

See readme_pgma_v2_1.txt for details*/


#include "pgma_v2_1.h"

void savegamut();

double gdata[3][4500];
struct LAB ogamut[SECTORS][SECTORS], rgamut[SECTORS][SECTORS], ogamutc[SECTORS][SECTORS];

void main(void)
{
	char fname[100];
	FILE *ofile;
	FILE *rfile;
	long int n,i,w,h,gma_i;
	struct LAB otest,rtest;
	unsigned char L,a,b;
	double ominL,omaxL,rminL,rmaxL;
	struct XYZ centre;
	double S_lut[1001];
	int m;


	extern double gdata[3][4500];
	extern struct LAB ogamut[SECTORS][SECTORS], rgamut[SECTORS][SECTORS], ogamutc[SECTORS][SECTORS];



	centre.x=centre.y=0.0, centre.z=50.0;

	printf("Test of public GMA code (2000 © Ján Morovic)\n\n");

	/*calculate original gamut boundary*/
	printf("Enter original gamut data file name: ");
	scanf("%s",fname);
	printf("Enter number of colours: ");
	scanf("%d",&n);

	ofile=fopen(fname,"r");

	for(i=0;i<n;i++)
		fscanf(ofile, "%lf %lf %lf", &gdata[0][i], &gdata[1][i], &gdata[2][i]);

	fclose(ofile);

	calc_gb(gdata,n,ogamut);

	savegamut(ogamut,"ogamut.sg");
	/**/

	/*calculate reproduction gamut boundary*/
	printf("Enter reproduction gamut data file name: ");
	scanf("%s",fname);
	printf("Enter number of colours: ");
	scanf("%d",&n);

	ofile=fopen(fname,"r");

	for(i=0;i<n;i++)
		fscanf(ofile, "%lf %lf %lf", &gdata[0][i], &gdata[1][i], &gdata[2][i]);

	fclose(ofile);

	calc_gb(gdata,n,rgamut);

	savegamut(rgamut,"rgamut.sg");
	/**/

	/*transform image*/
	printf("Enter original raw LAB image file name: ");
	scanf("%s",fname);
	printf("Enter width: ");
	scanf("%d",&w);
	printf("Enter height: ");
	scanf("%d",&h);

	ofile=fopen(fname,"rb");

	printf("Enter reproduction raw LAB image file name: ");
	scanf("%s",fname);

	rfile=fopen(fname,"wb");


	printf("Available GMAs:\n\n");
	printf("1   HPMINDE\n");
	printf("2   SGCK\n\n");

	printf("Enter choice: ");
	scanf("%d",&gma_i);

	if(gma_i==1)
		init_hp_minde(rgamut,&rminL,&rmaxL);
	else if(gma_i==2)
	{
		printf("Chose number of points (<=1000) used in the S-curve LUT:  ");
		scanf("%d",&m);
		init_sgck(ogamut,rgamut,ogamutc,&ominL,&omaxL,&rminL,&rmaxL,S_lut,m);
	}
	else
		printf("Error!\n");

	for(i=0;i<(w*h);i++)
	{
		fscanf(ofile, "%c%c%c", &L, &a, &b);

		otest.L=((double) L)/2.55;
		otest.a=((double) a)-128.0;
		otest.b=((double) b)-128.0;
		otest.C=sqrt(otest.a*otest.a + otest.b*otest.b);

		spherical(centre,otest.a,otest.b,otest.L,&(otest.alpha),&(otest.theta),&(otest.r));

		if(gma_i==1)
			hp_minde(rgamut,rminL,rmaxL,otest,&rtest);
		else if(gma_i==2)
			sgck(ogamutc,rgamut,ominL,omaxL,otest,&rtest,S_lut,m);

		rtest.L=(rtest.L<0.0) ? 0.0 : ((rtest.L>100.0) ? 100.0 : rtest.L);
		rtest.a=(rtest.a<-128.0) ? -128.0 : ((rtest.a>127.0) ? 127.0 : rtest.a);
		rtest.b=(rtest.b<-128.0) ? -128.0 : ((rtest.b>127.0) ? 127.0 : rtest.b);

		L=(unsigned char) floor(rtest.L*2.55);
		a=(unsigned char) floor(rtest.a+128.0);
		b=(unsigned char) floor(rtest.b+128.0);

		fprintf(rfile, "%c%c%c", L, a, b);

		if(i%(w*h/10)==0)
			printf("• %.2f per cent\n",100.0*i/(w*h));

	}

	fclose(ofile);
	fclose(rfile);
	/**/

}

/*saves information in gamut boundary descriptor to a file*/
void savegamut(	struct LAB gamut[SECTORS][SECTORS],
				char filen[32])
{
	FILE *outfile;
	int i,j;

	outfile=fopen(filen, "w");

	for (j=0; j<SECTORS; j++)
	{
		for (i=0; i<SECTORS; i++)
		{
			fprintf(outfile, "%d %f %f %f %f %f %f %f\n",	gamut[j][i].type,
															gamut[j][i].L,
															gamut[j][i].a,
															gamut[j][i].b,
															gamut[j][i].C,
															gamut[j][i].alpha,
															gamut[j][i].theta,
															gamut[j][i].r);
		}
	}
	fclose(outfile);

}
