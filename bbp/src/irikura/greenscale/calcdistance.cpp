#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define E2 0.0066943800229
#define PAI 3.14159265358979323846264
#define RADIUS 6378137.0

void SetXYZ(double dLon,
			double dLat,
			double dHgt,
			double& dX, 
			double& dY, 
			double& dZ )
{
  double N;
  const double a = RADIUS;
         N=a/sqrt(1 - E2 * pow(sin(dLat), 2.0));
         dX = (N + dHgt) * cos(dLat) * cos(dLon);
         dY = (N + dHgt) * cos(dLat) * sin(dLon);
         dZ = (N * (1 - E2) + dHgt) * sin(dLat);
  return;
}

void Distance(double dLonA, 
			  double dLatA,
			  double dHgtA, 
			  double dLonB,
			  double dLatB,
			  double dHgtB,
			  double &dRet)
{
  double dXA, dYA, dZA, dXB, dYB, dZB;
  double dLonRadA = dLonA * PAI/180;
  double dLatRadA = dLatA * PAI/180;
  double dLonRadB = dLonB * PAI/180;
  double dLatRadB = dLatB * PAI/180;

  SetXYZ(dLonRadA, dLatRadA, dHgtA, dXA, dYA, dZA);
  SetXYZ(dLonRadB, dLatRadB, dHgtB, dXB, dYB, dZB);
  
  dRet=sqrt(pow(dXA-dXB,2.0)+pow(dYA-dYB,2.0)+pow(dZA-dZB,2.0));

  return;
}
