/*Public Gamut Mapping Algorithm (PGMA) C Source Code - Version 2.1
Author: Ján Morovic <j.morovic@derby.ac.uk>
Release date: 26 July 2002

pgma_a_vgeom_v2_1.c - auxiliary vector geometry functions used in PGMA

2002 © Ján Morovic

See readme_pgma_v2_1.txt for details*/


#include "pgma_v2_1.h"

/*returns plane determined by three points a, b and c*/
struct PLANE points2plane(	struct XYZ a,
							struct XYZ b,
							struct XYZ c)
{
	struct PLANE temp;

	temp.b1=a.x;
	temp.b2=a.y;
	temp.b3=a.z;
	temp.v1=b.x-a.x;
	temp.v2=b.y-a.y;
	temp.v3=b.z-a.z;
	temp.w1=c.x-a.x;
	temp.w2=c.y-a.y;
	temp.w3=c.z-a.z;

	return temp;
}

/*returns line connecting two points a and b*/
struct LINE points2line(	struct XYZ a,
							struct XYZ b)
{
	struct LINE temp;

	temp.a1=a.x;
	temp.a2=a.y;
	temp.a3=a.z;
	temp.u1=b.x-a.x;
	temp.u2=b.y-a.y;
	temp.u3=b.z-a.z;

	return temp;
}

/*tests whether point t is in triangle abc and diff returns the difference in area between
triangle abc and the sum of areas of triangles tbc, atc and abt
the function returns 1 if the difference is below the threshold INTRITHLD ans 0 otherwise*/
int testintria(	struct LAB a,
				struct LAB b,
				struct LAB c,
				struct LAB t,
				double *diff)
{
	double area[3],tarea;

	tarea=triaarea(a,b,c);

	area[0]=triaarea(t,b,c);
	area[1]=triaarea(a,t,c);
	area[2]=triaarea(a,b,t);

	*diff=area[0]+area[1]+area[2]-tarea;

	if(*diff<INTRITHLD)
		return 1;
	else
		return 0;

}

/*returns coords of intercept between line and plane [success contains 0 if there is lin. dependence*/
struct XYZ ppintercept(	struct LINE line,
						struct PLANE plane,
						int *success)
{
	double m[MSIZE][MSIZE];
	struct XYZ temp;
	int n,i,j,k,test;

	m[0][0]=plane.v1;
	m[0][1]=plane.w1;
	m[0][2]=-line.u1;
	m[0][3]=line.a1-plane.b1;

	m[1][0]=plane.v2;
	m[1][1]=plane.w2;
	m[1][2]=-line.u2;
	m[1][3]=line.a2-plane.b2;

	m[2][0]=plane.v3;
	m[2][1]=plane.w3;
	m[2][2]=-line.u3;
	m[2][3]=line.a3-plane.b3;

	n=3;

	for (i=0; i<n-1; i++)					/*test for linear dependence in nxn part*/
 	{
		for (j=i+1; j<n; j++)
		{
			test=0;
			for (k=0; k<n; k++)
				if ((m[i][k])==(m[j][k]*m[i][0]/m[j][0]))
					test++;
		}
 		if (test==n) /*m is linearly dependent*/
 		{
 			temp.x=0.0;
 			temp.y=0.0;
 			temp.z=50.0;
 			*success=0;
 			return temp;
 		}
 	}


	solvegem(m,3,1);

	temp=poline(line,m[2][3]);

	*success=1;
	return temp;
}

/*calculates alpha theta and radius of XYZ coords
alpha goes from 0–360 and theta from 0–180, whereby 0 corresponds to 90 and 180 to 270
this has been done by using L* as the x axis and the radius in the xy plane as the y axis,
which is always positive*/
void spherical(	struct XYZ centre,
				double x,
				double y,
				double z,
				double *alpha,
				double *theta,
				double *r)
{
	double ia,ib,iL,a,b;

	ia=x-centre.x;
	ib=y-centre.y;
	iL=z-centre.z;

	arctan(ia,ib,alpha);

	a=sqrt(ia*ia + ib*ib);
	b=iL;

	arctan(b,a,theta);

	*r=sqrt(iL*iL + ia*ia + ib*ib);
}

/*converts spherical coordinates to orthogonal ones*/
void spher2ortho(	struct XYZ centre,
					double alpha,
					double theta,
					double r,
					double *x,
					double *y,
					double *z)
{
	double ta,tt,at,bt,lt,pi;
	struct XYZ centre2;

	pi=acos(-1.0);

	theta=90.0-theta;
	if(fabs(theta)<0.0001)
		theta=0.0001;


	centre2=centre;
	ta=tan(alpha*pi/180.0)*tan(alpha*pi/180.0);
	tt=tan(theta*pi/180.0)*tan(theta*pi/180.0);

	at=r*r/(ta + 1 + tt + tt*ta);
	bt=ta*at;
	lt=r*r-at-bt;

	if(theta==0.0)
		lt=0.0;

	at=sqrt(at);
	bt=sqrt(bt);
	lt=sqrt(lt);

	if(alpha==90.0 || alpha==270.0)
	{
		at=0.0;
		bt=fabs(r*(cos(theta*pi/180.0)));
		lt=fabs(r*(sin(theta*pi/180.0)));
	}

	if(alpha>=90.0 && alpha<=270.0)
		at=-at;
	if(alpha>=180.0 && alpha<=360.0)
		bt=-bt;
	if(theta<0.0)
		lt=-lt;

	*x=at+centre.x;
	*y=bt+centre.y;
	*z=lt+centre.z;

	if(theta==90.0)
	{
		*x=centre.x;
		*y=centre.y;
		*z=centre.z+r;
	}
	else if(theta==-90.0)
	{
		*x=centre.x;
		*y=centre.y;
		*z=centre.z-r;
	}

}

/*returns closest point on line1 to line2*/
struct XYZ cpointll(	struct LINE line1,
						struct LINE line2)
{
	struct XYZ mdtemp;
	double dist,mindist;
	struct XYZ ep,lp;
	int i;

	ep=poline(line2,0.0);
	mdtemp=closestpl(line1,ep);
	mindist=xyzdist(ep.x,ep.y,ep.z,mdtemp.x,mdtemp.y,mdtemp.z);

	for(i=1;i<301;i++)
	{
		ep=poline(line2,(double) i/300);
		lp=closestpl(line1,ep);
		dist=xyzdist(lp.x,lp.y,lp.z,ep.x,ep.y,ep.z);
		if(dist<mindist)
		{
			mindist=dist;
			mdtemp=lp;
		}
	}

	return mdtemp;

}

/*returns XYZ coords of point on line with a given t parameter value*/
struct XYZ poline(	struct LINE line,
					double t)
{
	struct XYZ temp;

	temp.x=line.a1+t*line.u1;
	temp.y=line.a2+t*line.u2;
	temp.z=line.a3+t*line.u3;

	return temp;
}

/*returns point on "line" closest to "point"*/
struct XYZ closestpl(	struct LINE line,
						struct XYZ point)
{
	struct PLANE ro;
	struct XYZ temp;
	int success;

	ro=orthoplane(line,point);

	temp=ppintercept(line,ro,&success);

	if(!success)/*should never happen ...*/
		printf("no intercept [closest point on line]\n");
	return temp;
}

/*euclidean distance between 2 points in 3D*/
double xyzdist(	double x1,
				double y1,
				double z1,
				double x2,
				double y2,
				double z2)
{
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}

/*CIELAB chroma and hue angle to a* and b**/
void ch2ab(	double c,
			double h,
			double *a,
			double *b)
{
	double pi;

	pi=acos(-1.0);

	*a=sqrt(c*c/(1+tan(h*pi/180.0)*tan(h*pi/180.0)));
	*b=sqrt(tan(h*pi/180.0)*tan(h*pi/180.0)*(*a)*(*a));

	if(h==90.0 || h==270.0)
	{
		*a=0.0;
		*b=c;
	}

	if(h>=90.0 && h<=270.0)
		*a=-*a;
	if(h>=180.0 && h<=360.0)
		*b=-*b;
}

/*returns area of triangle abc*/
double triaarea(	struct LAB a,
					struct LAB b,
					struct LAB c)
{
	struct XYZ u,v,w;

	u.x=b.a-a.a;
	u.y=b.b-a.b;
	u.z=b.L-a.L;

	v.x=c.a-a.a;
	v.y=c.b-a.b;
	v.z=c.L-a.L;

	w=vectorprod(u,v);

	return vectorsize(w)/2;

}

/*returns plane which is orthogonal to "line" and contains "point"*/
struct PLANE orthoplane(	struct LINE line,
							struct XYZ point)
{
	struct PLANE ro;

	/*ro is orthogonal to line and contains point*/
	ro.b1=point.x;
	ro.b2=point.y;
	ro.b3=point.z;

	if(line.u1==0.0)
	{
		if(line.u2==0.0)
		{
			ro.v1=10.0;
			ro.v2=0.0;
			ro.v3=0.0;
			ro.w1=10.0;
			ro.w2=10.0;
			ro.w3=0.0;
		}
		else if(line.u3==0.0)
		{
			ro.v1=10.0;
			ro.v2=0.0;
			ro.v3=0.0;
			ro.w1=10.0;
			ro.w2=0.0;
			ro.w3=10.0;
		}
		else
		{
			ro.v1=10.0;
			ro.v2=0.0;
			ro.v3=0.0;
			ro.w1=0.0;
			ro.w2=10.0;
			ro.w3=-(10*line.u2/line.u3);
		}
	}
	else if(line.u2==0.0)
	{
		if(line.u3==0.0)
		{
			ro.v1=0.0;
			ro.v2=10.0;
			ro.v3=0.0;
			ro.w1=0.0;
			ro.w2=10.0;
			ro.w3=10.0;
		}
		else
		{
			ro.v1=0.0;
			ro.v2=10.0;
			ro.v3=0.0;
			ro.w1=10.0;
			ro.w2=0.0;
			ro.w3=-(10*line.u1/line.u3);
		}
	}
	else if(line.u3==0.0)
	{
		ro.v1=0.0;
		ro.v2=0.0;
		ro.v3=10.0;
		ro.w1=10.0;
		ro.w2=-(10*line.u1/line.u2);
		ro.w3=10.0;
	}
	else
	{
		ro.v1=0.0;
		ro.v2=10.0;
		ro.v3=-(10*line.u2/line.u3);
		ro.w1=10.0;
		ro.w2=0.0;
		ro.w3=-(10*line.u1/line.u3);
	}

	return ro;
}

/*returns vector which is the vector product of two input vectors*/
struct XYZ vectorprod(	struct XYZ u,
						struct XYZ v)
{
	struct XYZ temp;

	temp.x=u.y*v.z-v.y*u.z;
	temp.y=u.z*v.x-v.z*u.x;
	temp.z=u.x*v.y-v.x*u.y;

	return temp;
}


/*returns size of vector*/
double vectorsize(struct XYZ u)
{
	return sqrt(u.x*u.x+u.y*u.y+u.z*u.z);
}