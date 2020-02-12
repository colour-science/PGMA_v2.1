/*Public Gamut Mapping Algorithm (PGMA) C Source Code - Version 2.1
Author: Ján Morovic <j.morovic@derby.ac.uk>
Release date: 26 July 2002

pgma_gbcalc_v2_1.c - functions for calculating gamut boundaries

2002 © Ján Morovic

See readme_pgma_v2_1.txt for details*/


#include "pgma_v2_1.h"

/*Implementation of segment maxima gamut boundary description method
gbdata		contains a b L data of samples from gamut
n			number of samples in gbdata
gamut		resulting GBD*/
void calc_gb(	double gbdata[3][4500],
				int n,
				struct LAB gamut[SECTORS][SECTORS])
{
	int alpha,theta;
	struct XYZ centre;
	int i,j;
	struct LAB test;

	for (j=0; j<SECTORS; j++)
		for (i=0; i<SECTORS; i++)
			gamut[j][i].type=-1;

	centre.x=centre.y=0.0, centre.z=50.0;

	/*finds those colours from "gbdata" that have largest radia in each segment*/
	for(i=0;i<n;i++)
	{
		test.a=gbdata[0][i];
		test.b=gbdata[1][i];
		test.L=gbdata[2][i];
		test.C=sqrt(test.a*test.a + test.b*test.b);

		spherical(centre,test.a,test.b,test.L,&(test.alpha),&(test.theta),&(test.r));

		alpha = (int) floor(test.alpha/(360.0/SECTORS));
		theta = (int) floor(test.theta/(180.0/SECTORS));
		if(alpha==SECTORS)
			alpha=SECTORS-1;
		if(theta==SECTORS)
			theta=SECTORS-1;

		if (gamut[theta][alpha].type==-1)
		{
			gamut[theta][alpha]=test;
			gamut[theta][alpha].type=3;
		}
		else if (gamut[theta][alpha].r<test.r)
		{
			gamut[theta][alpha]=test;
			gamut[theta][alpha].type=3;
		}

	}

	/*if any segments are left empty, values for them are interpolated*/
	interp_gb(gamut);

}

/*interpolates missing points in GBD
gamut		GBD to be interpolated*/
void interp_gb(struct LAB gamut[SECTORS][SECTORS])
{
	int i,j,k,l,m,n,success,intria;
	struct LAB templ,gpoint,closel;
	struct LINE ray,edge;
	struct PLANE ro;
	double diff,dist,mindist;
	struct XYZ temp,closez[SECTORS],ep,lp;
	struct XYZ centre;
	struct LAB close[SECTORS],closef[SECTORS];

	centre.x=centre.y=0.0, centre.z=50.0;

	j=0;
	for(i=0;i<SECTORS;i++)
	{
		if(gamut[j][i].type==-1)
		{
			for(n=0;n<SECTORS;n++)
			{
				close[n].type=-1;
				closef[n].type=-1;
			}
			n=0;
			m=i;
			if(m<0)
				m=SECTORS-1;
			for(k=1;k<(SECTORS+1);k++)
			{
				l=j+k;
				if(l>(SECTORS-1))
				{
					k=0;
					l=j+k;
					m=m-1;
					if(m<0)
						m=SECTORS-1;

					if(m==i)
						k=(SECTORS+1);
				}

				if(gamut[l][m].type!=-1)
				{
					close[n]=gamut[l][m];
					n++;

					k=-1;
					m=m-1;
					if(m<0)
						m=SECTORS-1;

					if(m==i)
						k=(SECTORS+1);
				}
			}


			templ.alpha=(i+0.5)*360.0/SECTORS;
			templ.theta=(j+0.5)*180.0/SECTORS;
			templ.r=50.0;

			spher2ortho(centre,templ.alpha,templ.theta,templ.r,&(templ.a),&(templ.b),&(templ.L));

			ray.a1=centre.x;
			ray.a2=centre.y;
			ray.a3=centre.z;
			ray.u1=templ.a-ray.a1;
			ray.u2=templ.b-ray.a2;
			ray.u3=templ.L-ray.a3;

			l=0;
			for(k=0;k<SECTORS;k++)
				if(close[k].type!=-1)
				{
					closef[l]=close[k];
					closez[l]=structlab2xyz(close[k]);
					l++;
				}

			closel.r=0.0;

			for(k=0;k<(l-1);k++)
				for(m=k+1;m<l;m++)
				{
					edge=points2line(closez[k],closez[m]);

					temp=cpointll(ray,edge);
					templ=structxyz2lab(centre,temp);

					if(templ.r>closel.r && templ.theta>=(j*180.0/SECTORS)
										&& templ.theta<=((j+1)*180.0/SECTORS)
										&& templ.alpha>=(i*360.0/SECTORS)
										&& templ.alpha<=((i+1)*360.0/SECTORS))
						closel=templ;
				}

			gamut[j][i]=closel;
			gamut[j][i].C=sqrt(gamut[j][i].a*gamut[j][i].a + gamut[j][i].b*gamut[j][i].b);
			gamut[j][i].type=2;


			if(gamut[j][i].type==2)
			{
				if(!(gamut[j][i].alpha>=(i*(360.0/SECTORS)) && gamut[j][i].alpha<=((i+1)*(360.0/SECTORS))))
					printf("Wrong hue angle - interpolate\n");
				if(!(gamut[j][i].theta>=(j*(180.0/SECTORS)) && gamut[j][i].theta<=((j+1)*(180.0/SECTORS))))
					printf("Wrong theta - interpolate\n");
			}
		}
	}

	j=SECTORS-1;
	for(i=0;i<SECTORS;i++)
	{
		if(gamut[j][i].type==-1)
		{
			for(n=0;n<SECTORS;n++)
			{
				close[n].type=-1;
				closef[n].type=-1;
			}
			n=0;
			m=i;
			if(m<0)
				m=SECTORS-1;
			for(k=1;k<(SECTORS+1);k++)
			{
				l=j-k;
				if(l<0)
				{
					k=0;
					l=j-k;
					m=m-1;
					if(m<0)
						m=SECTORS-1;

					if(m==i)
						k=(SECTORS+1);
				}

				if(gamut[l][m].type!=-1)
				{
					close[n]=gamut[l][m];
					n++;

					k=-1;
					m=m-1;
					if(m<0)
						m=SECTORS-1;

					if(m==i)
						k=(SECTORS+1);
				}
			}


			templ.alpha=(i+0.5)*360.0/SECTORS;
			templ.theta=(j+0.5)*180.0/SECTORS;
			templ.r=50.0;

			spher2ortho(centre,templ.alpha,templ.theta,templ.r,&(templ.a),&(templ.b),&(templ.L));

			ray.a1=centre.x;
			ray.a2=centre.y;
			ray.a3=centre.z;
			ray.u1=templ.a-ray.a1;
			ray.u2=templ.b-ray.a2;
			ray.u3=templ.L-ray.a3;

			l=0;
			for(k=0;k<SECTORS;k++)
				if(close[k].type!=-1)
				{
					closef[l]=close[k];
					closez[l]=structlab2xyz(close[k]);
					l++;
				}

			closel.r=0.0;

			for(k=0;k<(l-1);k++)
				for(m=k+1;m<l;m++)
				{
					edge=points2line(closez[k],closez[m]);

					temp=cpointll(ray,edge);
					templ=structxyz2lab(centre,temp);

					if(templ.r>closel.r && templ.theta>=(j*180.0/SECTORS)
										&& templ.theta<=((j+1)*180.0/SECTORS)
										&& templ.alpha>=(i*360.0/SECTORS)
										&& templ.alpha<=((i+1)*360.0/SECTORS))
						closel=templ;
				}

			gamut[j][i]=closel;
			gamut[j][i].C=sqrt(gamut[j][i].a*gamut[j][i].a + gamut[j][i].b*gamut[j][i].b);
			gamut[j][i].type=2;


			if(gamut[j][i].type==2)
			{
				if(!(gamut[j][i].alpha>=(i*(360.0/SECTORS)) && gamut[j][i].alpha<=((i+1)*(360.0/SECTORS))))
					printf("Wrong hue angle - interpolate\n");
				if(!(gamut[j][i].theta>=(j*(180.0/SECTORS)) && gamut[j][i].theta<=((j+1)*(180.0/SECTORS))))
					printf("Wrong theta - interpolate\n");
			}
		}
	}

	for(j=1;j<(SECTORS-1);j++)
	{
		for(i=0;i<SECTORS;i++)
		{
			if(gamut[j][i].type==-1)
			{
				if(j!=0 && j!=(SECTORS-1))
				{
					for(n=0;n<SECTORS;n++)
					{
						close[n].type=-1;
						closef[n].type=-1;
					}
					n=0;
					m=i-1;
					if(m<0)
						m=SECTORS-1;
					for(k=1;k<SECTORS;k++)
					{
						l=j+k;
						if(l>(SECTORS-1))
						{
							m=m-1;
							if(m<0)
								m=SECTORS-1;
							if(m==(i+1))
							{
								n=n+1;
								l=j-n;
							}

							k=1;
							l=j+k;
						}

						if(gamut[l][m].type!=-1)
						{
							close[0]=gamut[l][m];
							k=SECTORS;
						}
					}

					n=0;
					m=i+1;
					if(m>=SECTORS)
						m=0;
					for(k=1;k<SECTORS;k++)
					{
						l=j+k;
						if(l>(SECTORS-1))
						{
							m=m+1;
							if(m>=SECTORS)
								m=0;
							if(m==(i-1))
							{
								n=n+1;
								l=j-n;
							}

							k=1;
							l=j+k;
						}

						if(gamut[l][m].type!=-1)
						{
							close[1]=gamut[l][m];
							k=SECTORS;
						}
					}

					n=0;
					m=i-1;
					if(m<0)
						m=SECTORS-1;
					for(k=1;k<SECTORS;k++)
					{
						l=j-k;
						if(l<0)
						{
							m=m-1;
							if(m<0)
								m=SECTORS-1;
							if(m==(i+1))
							{
								n=n-1;
								l=j+n;
							}
							k=1;
							l=j-k;
						}

						if(gamut[l][m].type!=-1)
						{
							close[2]=gamut[l][m];
							k=SECTORS;
						}
					}

					n=0;
					m=i+1;
					if(m>=SECTORS)
						m=0;
					for(k=1;k<SECTORS;k++)
					{
						l=j-k;
						if(l<0)
						{
							m=m+1;
							if(m>=SECTORS)
								m=0;
							if(m==(i-1))
							{
								n=n-1;
								l=j+n;
							}
							k=1;
							l=j-k;
						}

						if(gamut[l][m].type!=-1)
						{
							close[3]=gamut[l][m];
							k=SECTORS;
						}
					}
				}

				/**/
				templ.alpha=(i+0.5)*360.0/SECTORS;
				templ.theta=(j+0.5)*180.0/SECTORS;
				templ.r=50.0;

				spher2ortho(centre,templ.alpha,templ.theta,templ.r,&(templ.a),&(templ.b),&(templ.L));

				ray.a1=centre.x;
				ray.a2=centre.y;
				ray.a3=centre.z;
				ray.u1=templ.a-ray.a1;
				ray.u2=templ.b-ray.a2;
				ray.u3=templ.L-ray.a3;

				l=0;
				for(k=0;k<4;k++)
					if(close[k].type!=-1)
					{
						closef[l]=close[k];
						closez[l]=structlab2xyz(close[k]);
						l++;
					}

				if(l==4)
				{
					intria=0;

					if(!intria)
					{
						ro.b1=closef[0].a;
						ro.b2=closef[0].b;
						ro.b3=closef[0].L;
						ro.v1=closef[1].a - ro.b1;
						ro.v2=closef[1].b - ro.b2;
						ro.v3=closef[1].L - ro.b3;
						ro.w1=closef[2].a - ro.b1;
						ro.w2=closef[2].b - ro.b2;
						ro.w3=closef[2].L - ro.b3;

						temp=ppintercept(ray,ro,&success);

						gpoint.a=temp.x;
						gpoint.b=temp.y;
						gpoint.L=temp.z;

						spherical(centre,gpoint.a,gpoint.b,gpoint.L,&(gpoint.alpha),&(gpoint.theta),&(gpoint.r));

						gpoint.C=sqrt(temp.x*temp.x+temp.y*temp.y);
						gpoint.type=2;

						intria=testintria(closef[0],closef[1],closef[2],gpoint,&diff);

						if(intria==1 && success==1)
						{
							gamut[j][i]=gpoint;
							gamut[j][i].C=sqrt(gamut[j][i].a*gamut[j][i].a + gamut[j][i].b*gamut[j][i].b);
						}
					}
					if(!intria)
					{
						ro.b1=closef[1].a;
						ro.b2=closef[1].b;
						ro.b3=closef[1].L;
						ro.v1=closef[3].a - ro.b1;
						ro.v2=closef[3].b - ro.b2;
						ro.v3=closef[3].L - ro.b3;
						ro.w1=closef[2].a - ro.b1;
						ro.w2=closef[2].b - ro.b2;
						ro.w3=closef[2].L - ro.b3;

						temp=ppintercept(ray,ro,&success);

						gpoint.a=temp.x;
						gpoint.b=temp.y;
						gpoint.L=temp.z;

						spherical(centre,gpoint.a,gpoint.b,gpoint.L,&(gpoint.alpha),&(gpoint.theta),&(gpoint.r));

						gpoint.C=sqrt(temp.x*temp.x+temp.y*temp.y);
						gpoint.type=2;

						intria=testintria(closef[1],closef[3],closef[2],gpoint,&diff);

						if(intria==1 && success==1)
						{
							gamut[j][i]=gpoint;
							gamut[j][i].C=sqrt(gamut[j][i].a*gamut[j][i].a + gamut[j][i].b*gamut[j][i].b);
						}
					}

				}
				else if(l==3)
				{
					intria=0;

					if(!intria)
					{
						ro.b1=closef[0].a;
						ro.b2=closef[0].b;
						ro.b3=closef[0].L;
						ro.v1=closef[1].a - ro.b1;
						ro.v2=closef[1].b - ro.b2;
						ro.v3=closef[1].L - ro.b3;
						ro.w1=closef[2].a - ro.b1;
						ro.w2=closef[2].b - ro.b2;
						ro.w3=closef[2].L - ro.b3;

						temp=ppintercept(ray,ro,&success);

						gpoint.a=temp.x;
						gpoint.b=temp.y;
						gpoint.L=temp.z;

						spherical(centre,gpoint.a,gpoint.b,gpoint.L,&(gpoint.alpha),&(gpoint.theta),&(gpoint.r));

						gpoint.C=sqrt(temp.x*temp.x+temp.y*temp.y);
						gpoint.type=2;

						intria=testintria(closef[0],closef[1],closef[2],gpoint,&diff);

						if(intria==1 && success==1)
						{
							gamut[j][i]=gpoint;
							gamut[j][i].C=sqrt(gamut[j][i].a*gamut[j][i].a + gamut[j][i].b*gamut[j][i].b);
						}
					}

				}
				else if(l==2)
				{
					edge.a1=closez[0].x;
					edge.a2=closez[0].y;
					edge.a3=closez[0].z;
					edge.u1=closez[1].x-edge.a1;
					edge.u2=closez[1].y-edge.a2;
					edge.u3=closez[1].z-edge.a3;

					ep=poline(edge,0.0);
					temp=closestpl(ray,ep);
					mindist=xyzdist(ep.x,ep.y,ep.z,temp.x,temp.y,temp.z);

					for(k=1;k<101;k++)
					{
						ep=poline(edge,(double) k/100.0);
						lp=closestpl(ray,ep);
						dist=xyzdist(lp.x,lp.y,lp.z,ep.x,ep.y,ep.z);
						if(dist<mindist)
						{
							mindist=dist;
							temp=lp;
						}
					}

					gamut[j][i]=structxyz2lab(centre,temp);
					gamut[j][i].type=2;
				}

				if(gamut[j][i].type==2)
				{
					if(!(gamut[j][i].alpha>=(i*(360.0/SECTORS)) && gamut[j][i].alpha<=((i+1)*(360.0/SECTORS))))
						printf("Wrong hue angle - interpolate\n");
					if(!(gamut[j][i].theta>=(j*(180.0/SECTORS)) && gamut[j][i].theta<=((j+1)*(180.0/SECTORS))))
						printf("Wrong theta - interpolate\n");
				}

			}
		}
	}
}


/*returns gamut at hue angle "h" with SECTORS +2 points where first and last points are
on L* axis
h			hue angle at which 2D boundary is to be found
gamut		GBD from which 2D boundary is to be found
gbh			2D boundary points
gbhl		lines between neighbouring points in gbh
centre		orthogonal coordinates of centre of spherical coordinates*/
void gb_at_hue(	double h,
				struct LAB gamut[SECTORS][SECTORS],
				struct LAB gbh[SECTORS+2],
				struct LINE2D gbhl[SECTORS+2],
				struct XYZ centre)
{
	struct PLANE ro;
	struct LINE l;
	struct XYZ a,b,c,g;
	int i,j,found;
	double hue,diff;
	struct LAB temp,gtemp;

	hue=h;

	while(hue>360.0)
		hue-=360.0;

	ro=hueplane(hue);

	/*find points 1 to SECTORS*/
	for(i=0;i<SECTORS;i++)
	{
		found=0;
		for(j=0;j<(SECTORS-1);j++)
		{
			if(gamut[i][j].alpha<=hue && gamut[i][j+1].alpha>=hue)
			{
				a=structlab2xyz(gamut[i][j]);
				b=structlab2xyz(gamut[i][j+1]);
				found=1;
				break;
			}
		}

		if(!found)
		{
			a=structlab2xyz(gamut[i][SECTORS-1]);
			b=structlab2xyz(gamut[i][0]);
		}

		l=points2line(a,b);

		g=ppintercept(l,ro,&found);

		gbh[i+1]=structxyz2lab(centre,g);

	}

	/*find point 0*/
	a.x=a.y=0.0;
	a.z=50.0;
	b.x=b.y=0.0;
	b.z=150.0;
	l=points2line(a,b);
	temp.L=0.0;

	for(i=0;i<SECTORS;i++)
	{
		if(gamut[0][i].L>temp.L)
			temp=gamut[0][i];
	}

	a=structlab2xyz(temp);
	found=0;
	for(i=0;i<(SECTORS-1);i++)
	{
		if(!(eqlab(temp,gamut[0][i]) || eqlab(temp,gamut[0][i+1])))
		{
			b=structlab2xyz(gamut[0][i]);
			c=structlab2xyz(gamut[0][i+1]);

			ro=points2plane(a,b,c);
			g=ppintercept(l,ro,&found);
			gtemp=structxyz2lab(centre,g);

			found=testintria(temp,gamut[0][i],gamut[0][i+1],gtemp,&diff);

			if(found)
				break;
		}
	}

	if(!found)
	{
		b=structlab2xyz(gamut[0][0]);
		c=structlab2xyz(gamut[0][SECTORS-1]);

		ro=points2plane(a,b,c);
		g=ppintercept(l,ro,&found);
		gtemp=structxyz2lab(centre,g);

		found=testintria(temp,gamut[0][0],gamut[0][SECTORS-1],gtemp,&diff);

	}

	gbh[0]=gtemp;

	/*find point SECTORS+1*/
	a.x=a.y=0.0;
	a.z=50.0;
	b.x=b.y=0.0;
	b.z=-50.0;
	l=points2line(a,b);
	temp.L=100.0;

	for(i=0;i<SECTORS;i++)
	{
		if(gamut[SECTORS-1][i].L<temp.L)
			temp=gamut[SECTORS-1][i];
	}

	a=structlab2xyz(temp);
	found=0;
	for(i=0;i<(SECTORS-1);i++)
	{
		if(!(eqlab(temp,gamut[SECTORS-1][i]) || eqlab(temp,gamut[SECTORS-1][i+1])))
		{
			b=structlab2xyz(gamut[SECTORS-1][i]);
			c=structlab2xyz(gamut[SECTORS-1][i+1]);

			ro=points2plane(a,b,c);
			g=ppintercept(l,ro,&found);
			gtemp=structxyz2lab(centre,g);

			found=testintria(temp,gamut[SECTORS-1][i],gamut[SECTORS-1][i+1],gtemp,&diff);

			if(found)
				break;
		}
	}

	if(!found)
	{
		b=structlab2xyz(gamut[SECTORS-1][0]);
		c=structlab2xyz(gamut[SECTORS-1][SECTORS-1]);

		ro=points2plane(a,b,c);
		g=ppintercept(l,ro,&found);
		gtemp=structxyz2lab(centre,g);

		found=testintria(temp,gamut[SECTORS-1][0],gamut[SECTORS-1][SECTORS-1],gtemp,&diff);

	}

	gbh[SECTORS+1]=gtemp;

	/*calculate parametric line functions*/

	for(i=0;i<(SECTORS+1);i++)
	{
		gbhl[i].a1=gbh[i].C;
		gbhl[i].a2=gbh[i].L;
		gbhl[i].u1=gbh[i+1].C - gbhl[i].a1;
		gbhl[i].u2=gbh[i+1].L - gbhl[i].a2;
	}

	gbhl[SECTORS+1].a1=gbh[SECTORS+1].C;
	gbhl[SECTORS+1].a2=gbh[SECTORS+1].L;
	gbhl[SECTORS+1].u1=gbh[0].C - gbhl[SECTORS+1].a1;
	gbhl[SECTORS+1].u2=gbh[0].L - gbhl[SECTORS+1].a2;
}


/*returns two points resulting from the intersection of a line and the 2D GB
hgbl		lines determining 2D gamut boundary
l			line to be intersected with hgbl
res			result of intersection
NOTE: the values in the XYZ struct are not 3D orthogonal coordinates but the first two values
are 2D coordinates in the constant-hue-angle plane and the third is the parameter of the
intersection point along the vecor l. This use of the XYZ struct is specific and exclusive
to some functions in this file.*/
int line_hgb_intersects(	struct LINE2D hgbl[SECTORS+2],
							struct LINE2D l,
							struct XYZ res[2])
{
	int i,j,n;
	struct XYZ tres[10],ares[SECTORS+2];
	struct XYZ tint;
	double cdiff,diff=9999999999.9;

	n=0;
	for(i=0;i<(SECTORS+2);i++)
	{
		ares[i].z=999999999.0;

		j=line_line_isect_2d(hgbl[i],l,&tint);

		if(j==2 && tint.x>=-0.01 && tint.y>=-0.01)
		{
			tres[n]=tint;
			n++;
		}
		else if (j==1)
			ares[i]=tint;

	}

	if(n==0)
	{
		for(i=0;i<(SECTORS+2);i++)
		{
			cdiff=dist_cl(ares[i].x,ares[i].y,hgbl[i].a1,hgbl[i].a2);
			if(cdiff<diff)
				res[0]=ares[i], diff=cdiff;
		}

		n=1;

		return n;
	}

	if(n>1)
	{
		for(i=0;i<(n-1);i++)
		{
			for(j=i+1;j<n;j++)
			{
				if(tres[j].z>tres[i].z)
				{
					tint.x=tres[j].x, tint.y=tres[j].y, tint.z=tres[j].z;
					tres[j].x=tres[i].x, tres[j].y=tres[i].y, tres[j].z=tres[i].z;
					tres[i].x=tint.x, tres[i].y=tint.y, tres[i].z=tint.z;
				}
			}
		}
		res[0].x=tres[0].x, res[0].y=tres[0].y, res[0].z=tres[0].z;
		res[1].x=tres[n-1].x, res[1].y=tres[n-1].y, res[1].z=tres[n-1].z;

	}
	else
	{
		res[0].x=tres[0].x, res[0].y=tres[0].y, res[0].z=tres[0].z;
		res[1]=res[0];
	}

	if(n==1)
	{
		res[1]=res[0];
		return n;
	}
	else
		return 2;
}

/*intersection of two lines in 2D
l			first line
m			second line
p			coordinates of intersection [special use of XYZ struct]
returns 0 if lines are linearly dependent, 1 if they intersect and 2 if the intersection
is within the l segment*/
int line_line_isect_2d(	struct LINE2D l,
						struct LINE2D m,
						struct XYZ *p)
{
	double tl,tm;
	double matrix[MSIZE][MSIZE];

	matrix[0][0]=l.u1;
	matrix[0][1]=-m.u1;
	matrix[0][2]=m.a1-l.a1;
	matrix[1][0]=l.u2;
	matrix[1][1]=-m.u2;
	matrix[1][2]=m.a2-l.a2;

	if(solvegem(matrix,2,1)==0)
		return 0;

	tl=matrix[0][2];
	tm=matrix[1][2];

	(*p).x=m.a1+tm*m.u1;
	(*p).y=m.a2+tm*m.u2;
	(*p).z=tm;

	if((tl<=1.001) && (tl>=-0.001)) //1.01 and -0.01 are used instead of 1 and 0 due to computational accuracy
		return 2;
	else
		return 1;

}

/*returns plane which is orthogonal to xy plane and at the angle "h" from the x axis
"h" is in degrees*/
struct PLANE hueplane(double h)
{
	struct PLANE ro;
	struct XYZ a,b,c;
	double hue;

	double pi;

	pi=acos(-1.0);

	hue=h;

	while(hue>360.0)
		hue-=360.0;

	a.x=a.y=a.z=0.0;
	b.x=b.y=0.0;
	b.z=10.0;

	if(hue<90.0 || hue >270.0)
	{
		c.z=0.0;
		c.x=10.0;
		c.y=c.x*tan(hue*pi/180.0);
	}
	else if(hue==90.0)
	{
		c.z=0.0;
		c.x=0.0;
		c.y=10.0;
	}
	else if(hue==270.0)
	{
		c.z=0.0;
		c.x=0.0;
		c.y=-10.0;
	}
	else
	{
		c.z=0.0;
		c.x=-10.0;
		c.y=c.x*tan(hue*pi/180.0);
	}

	ro=points2plane(a,b,c);

	return ro;
}


/*converts special XYZ struct to LAB
centre		centre of spherical coordinates
point		special XYZ struct to be converted
hue			hue angle*/
struct LAB hgb_structxyz2lab(	struct XYZ centre,
								struct XYZ point,
								double hue)
{
	struct LAB temp;

	temp.C=point.x;
	temp.L=point.y;
	temp.alpha=hue;
	ch2ab(temp.C,temp.alpha,&(temp.a),&(temp.b));
	spherical(centre,temp.a,temp.b,temp.L,&(temp.alpha),&(temp.theta),&(temp.r));

	return temp;
}

/*2D euclidean distance*/
double dist_cl(	double c1,
				double l1,
				double c2,
				double l2)
{
	double tmp;
	tmp=sqrt((c2-c1)*(c2-c1) + (l2-l1)*(l2-l1));
	return tmp;
}

/*finds cusp of gamut at hue angle h*/
struct LAB findcusp(	struct LAB gamut[SECTORS][SECTORS],
						double h)
{

	struct LAB gbh[SECTORS+2],cusp;
	struct LINE2D gbhl[SECTORS+2];
	struct XYZ centre;
	int i;

	centre.x=centre.y=0.0, centre.z=50.0;

	gb_at_hue(h,gamut,gbh,gbhl,centre);

	cusp.C=-1.0;

	for(i=0;i<(SECTORS+2);i++)
		if(gbh[i].C>cusp.C)
			cusp=gbh[i];

	return cusp;
}
