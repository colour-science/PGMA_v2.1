/*Public Gamut Mapping Algorithm (PGMA) C Source Code - Version 2.1
Author: Ján Morovic <j.morovic@derby.ac.uk>
Release date: 26 July 2002

pgma_gmas_v2_1.c - gamut mapping algorithms

2002 © Ján Morovic & Pei-Li Sun

See readme_pgma_v2_1.txt for details*/

#include "pgma_v2_1.h"

/*initialisation for Sigmoid_Gaussian_Cusp_Knee (SGCK) algorithm
ogamut		original GBD
rgamut		reproduction GBD
ogamutc		lightness-compressed original GBD
ominL		minimum lightness of original gamut
omaxL		maximum lightness of original gamut
rminL		minimum lightness of reproduction gamut
rmaxL		maximum lightness of reproduction gamut
S_lut		sigmoidal lightness rescaling look-up-table
m			the number of points used in the discrete LUT*/
void init_sgck(	struct LAB ogamut[SECTORS][SECTORS],
				struct LAB rgamut[SECTORS][SECTORS],
				struct LAB ogamutc[SECTORS][SECTORS],
				double *ominL,
				double *omaxL,
				double *rminL,
				double *rmaxL,
				double S_lut[1001],
				int m)
{
	int i,j;
	struct XYZ centre;
	double w,pi=acos(-1.0),Xo,sigma,Si;

	/*look-up-table for Xo and sigma (normal lightness class)
	refer to Journal of Electronic Imaging (1999) Vol. 8(4):391,Table 3*/
	double table[4][2]= {{53.7,43.0},{56.8,40.0},{58.2,35.0},{60.6,34.5}};

	centre.x=centre.y=0.0, centre.z=50.0;

	*ominL=200.0;
	*rminL=200.0;
	*omaxL=-200.0;
	*rmaxL=-200.0;

	/*find lightenss ranges*/
	for (j=0; j<SECTORS; j++)
		for (i=0; i<SECTORS; i++)
		{
			ogamutc[i][j]=ogamut[i][j];

			if(ogamut[i][j].L<*ominL)
				*ominL=ogamut[i][j].L;
			if(ogamut[i][j].L>*omaxL)
				*omaxL=ogamut[i][j].L;

			if(rgamut[i][j].L<*rminL)
				*rminL=rgamut[i][j].L;
			if(rgamut[i][j].L>*rmaxL)
				*rmaxL=rgamut[i][j].L;
		}
	/**/


	/*calculate Xo and sigma for sigmoidal lightness rescaling*/
	i= (int)floor(*rminL/5)-1;
	w= (*rminL-i*5-5)/5.0;
	if (i>2) {i=2;w=1;}
	if (i<0) {i=0;w=0;}
	Xo=    (1.0-w)*table[i][0] +w*table[i+1][0];
	sigma= (1.0-w)*table[i][1] +w*table[i+1][1];
	/**/


	/*create S_lut for sigmoidal lightness rescaling*/
	for (i=0;i<=m;i++)
	{
		Si=1.0/(sqrt(2.0*pi)*sigma)*pow(exp(1.0),-(100.0/m*i-Xo)*(100.0/m*i-Xo)/(2.0*sigma*sigma));
		if (i==0) S_lut[0]=Si;
		else S_lut[i]=S_lut[i-1]+Si;
	}

	for (i=1;i<=m;i++)
		S_lut[i]=(S_lut[i]-S_lut[0])/(S_lut[m]-S_lut[0])*(*rmaxL-*rminL)+*rminL;

	S_lut[0]=*rminL;
	/**/


	/*compress original gamut goundary descriptor ogamut*/
	for (j=0; j<SECTORS; j++)
		for (i=0; i<SECTORS; i++)
		{
			ogamutc[i][j]=ogamut[i][j];

			/*normalise to [0,100]*/
			ogamutc[i][j].L= 100.0 -(*omaxL-ogamut[i][j].L) *100.0/(*omaxL-*ominL);

			/*sigmoidal lightness rescaling*/
			if(ogamutc[i][j].L<100.0)
			{
				w= floor(ogamutc[i][j].L*m/100.0);
				ogamutc[i][j].L= (1.0 -ogamutc[i][j].L*m/100.0 +w)*S_lut[(int)w]
								+(ogamutc[i][j].L*m/100.0 -w)*S_lut[(int)w +1];
			}
			else
				ogamutc[i][j].L=S_lut[m];

			/*Chroma dependent lightness rescaling*/
			ogamutc[i][j].L= C_dept_L(ogamut[i][j].C,ogamut[i][j].L,ogamutc[i][j].L);

			spherical(centre,ogamutc[i][j].a,ogamutc[i][j].b,ogamutc[i][j].L,&(ogamutc[i][j].alpha),&(ogamutc[i][j].theta),&(ogamutc[i][j].r));
		}
}

/*SGCK gamut mapping algorithm
ogamut		original GBD [lightness compressed]
rgamut		reproduction gamut
ominL		minimum lightness of original gamut
omaxL		maximum lightness of original gamut
rminL		minimum lightness of reproduction gamut
rmaxL		maximum lightness of reproduction gamut
otest		colour to be gamut-mapped
rtest		resulting gamut-mapped colour
S_lut		sigmoidal lightness rescaling look-up-table
m			the number of points used in the discrete LUT*/
void sgck(	struct LAB ogamut[SECTORS][SECTORS],
			struct LAB rgamut[SECTORS][SECTORS],
			double ominL,
			double omaxL,
			struct LAB otest,
			struct LAB *rtest,
			double S_lut[1001],
			int m)
{
	struct LAB gbho[SECTORS+2],gbhr[SECTORS+2];
	struct LINE2D gbhlo[SECTORS+2],gbhlr[SECTORS+2],mline;
	struct XYZ reso[2],resr[2],centre;
	struct LAB cuspo;
	double w,t;

	centre.x=centre.y=0.0, centre.z=50.0;

	*rtest=otest;

	/*normalise to [0,100]*/
	if(otest.L<ominL) otest.L=ominL;
	if(otest.L>omaxL) otest.L=omaxL;
	(*rtest).L= 100.0-(omaxL-otest.L)*100.0/(omaxL-ominL);

	/*sigmoidal lightness rescaling*/
	if((*rtest).L<100.0)
	{
		w= floor((*rtest).L*m/100.0);
		(*rtest).L= (1.0 -(*rtest).L*m/100.0 +w)*S_lut[(int)w]
					+((*rtest).L*m/100.0 -w)*S_lut[(int)w +1];
	}
	else
		(*rtest).L=S_lut[m];

	/*chroma dependent lightness mapping*/
	(*rtest).L=C_dept_L(otest.C,otest.L,(*rtest).L);

	/*find gamut boundary at given colour's hue*/
	gb_at_hue(otest.alpha,ogamut,gbho,gbhlo,centre);
	gb_at_hue((*rtest).alpha,rgamut,gbhr,gbhlr,centre);

	cuspo=findcusp(rgamut,otest.alpha);

	mline.a1=0.0;
	mline.a2=cuspo.L;
	mline.u1=(*rtest).C - mline.a1;
	mline.u2=(*rtest).L - mline.a2;

	line_hgb_intersects(gbhlo,mline,reso);
	line_hgb_intersects(gbhlr,mline,resr);

	/*chroma copression*/
	/*knee-point is set at 90 perc. of destination-gamut range*/
	if(mline.u1>(0.9*resr[0].x))
		t=(0.1*resr[0].x*(mline.u1-0.9*resr[0].x)/(reso[0].x-0.9*resr[0].x)+0.9*resr[0].x)/mline.u1;
	else t=1.0;
	if (resr[0].x>reso[0].x) t=1.0;

	(*rtest).C=mline.a1+t*mline.u1;
	(*rtest).L=mline.a2+t*mline.u2;
	(*rtest).alpha=otest.alpha;

	ch2ab((*rtest).C,(*rtest).alpha,&((*rtest).a),&((*rtest).b));
	spherical(centre,(*rtest).a,(*rtest).b,(*rtest).L,&((*rtest).alpha),&((*rtest).theta),&((*rtest).r));

}


/*1D chroma–dependent linear compression - based on GCUSP
o			original value
omin		minimum value of original range
omax		maximum value of original range
rmin		minimum value of reproduction range
rmax		maximum value of reproduction range
returns linearly compressed reproduction value*/
double C_dept_L(	double C,
					double o,
					double s)
{
	double p;

	p=1.0-sqrt(pow(C,3.0)/(pow(C,3.0)+500000.0));
	return ((1.0-p)*o + p*s);
}

/*initialisation for HPMINDE algorithm
rgamut		reproduction GBD
rminL		minimum lightness of reproduction gamut
rmaxL		maximum lightness of reproduction gamut*/
void init_hp_minde(	struct LAB rgamut[SECTORS][SECTORS],
					double *rminL,
					double *rmaxL)
{
	int i,j;
	struct XYZ centre;

	centre.x=centre.y=0.0, centre.z=50.0;

	*rminL=200.0;
	*rmaxL=-200.0;

	/*find lightenss ranges*/
	for (j=0; j<SECTORS; j++)
		for (i=0; i<SECTORS; i++)
		{
			if(rgamut[i][j].L<*rminL)
				*rminL=rgamut[i][j].L;
			if(rgamut[i][j].L>*rmaxL)
				*rmaxL=rgamut[i][j].L;
		}
	/**/
}

/*hue preserving minimum ∆E gamut mapping algorithm
rgamut		reproduction gamut
otest		colour to be gamut-mapped
rtest		resulting gamut-mapped colour*/
void  hp_minde(	struct LAB rgamut[SECTORS][SECTORS],
				double rminL,
				double rmaxL,
				struct LAB otest,
				struct LAB *rtest
			)
{
	struct XYZ centre,tempi,cl1,cl2,res[2];
	struct LAB gbh[SECTORS+2],rg;
	struct LINE2D gbhl[SECTORS+2],orthol,mline;
    int i,in,in2;
    double mdist1,mdist2,dist;

	centre.x=centre.y=0.0; centre.z=50.0;

	gb_at_hue(otest.alpha,rgamut,gbh,gbhl,centre);

    /*find whether point is in gamut*/
    mline.a1=0.0;
	mline.a2=(rmaxL+rminL)/2.0;
	mline.u1=otest.C - mline.a1;
	mline.u2=otest.L - mline.a2;

        if(fabs(mline.u1)+fabs(mline.u2)<0.00001)
            return;
        else
        {
            line_hgb_intersects(gbhl,mline,res);

            rg=hgb_structxyz2lab(centre,res[0],otest.alpha);
        }

	if(rg.C>otest.C && otest.L<=rmaxL && otest.L>=rminL)
		(*rtest)=otest;
	else	/*if otest is outside the reproduction gamut*/
	{
		mdist1=mdist2=999999.9;
	    in2=1;
	    /*find orthogonal intersection that is inside segment*/
	    for(i=0;i<(SECTORS+2);i++)
	    {
	    	orthol.a1=otest.C;
	    	orthol.a2=otest.L;
	    	orthol.u1=gbhl[i].u2;
	    	orthol.u2=-gbhl[i].u1;

	    	in=line_line_isect_2d(gbhl[i],orthol,&(tempi));

	    	dist=sqrt((otest.C-tempi.x)*(otest.C-tempi.x)+(otest.L-tempi.y)*(otest.L-tempi.y));

	    	if(in==2 && mdist1>dist)
	    	{
	    		in2=2;
	    		mdist1=dist;
	    		cl1=tempi;
	    	}
	    }

		/*check whether cl1 is min∆E*/
    	for(i=0;i<(SECTORS+2);i++)
    	{
    		dist=sqrt((otest.C-gbh[i].C)*(otest.C-gbh[i].C)+(otest.L-gbh[i].L)*(otest.L-gbh[i].L));

    		if(mdist2>dist)
	    	{
	    		mdist2=dist;
	    		cl2.x=gbh[i].C;
	    		cl2.y=gbh[i].L;
			}
		}

		if(mdist1<mdist2)
			(*rtest)=hgb_structxyz2lab(centre,cl1,otest.alpha);
		else
			(*rtest)=hgb_structxyz2lab(centre,cl2,otest.alpha);
	}
}
