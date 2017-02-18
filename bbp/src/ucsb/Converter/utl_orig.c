#include <stdio.h>
#include <string.h>
#include <math.h>

#define	E_RAD	6378.163
#define E_FLAT  298.26
#define DRAD	1.7453292e-2
#define DRLT	9.9330647e-1

static double	orig_lon;		/* origin longitude */
static double	orig_lat;		/* origin latitude  */
static double	lat_fac;		/* conversion factor for latitude in km */
static double	lon_fac;		/* conversion factor for longitude in km */
static double	snr;			/* sin of rotation angle */
static double	csr;			/* cos of rotation angle */

void utl_orig(double lon,double lat,double rota)
{
   double	dlt1;
   double	dlt2;
   double	del;
   double	radius;


   /* convert everything to minutes */
   lat = 60.0 * lat;
   lon = 60.0 * lon;
   orig_lon = lon;
   orig_lat = lat;

   /* latitude */
   dlt1 = atan(DRLT * tan((double)lat * DRAD/60.0));
   dlt2 = atan(DRLT * tan(((double)lat +1.0) * DRAD/60.0));
   del  = dlt2 - dlt1;
   radius = E_RAD * (1.0 - (sin(dlt1)*sin(dlt1) / E_FLAT));
   lat_fac = radius * del;

   /* longitude */
   del = acos(1.0 - (1.0 - cos(DRAD/60.0)) * cos(dlt1) * cos(dlt1));
   dlt2 = radius * del;
   lon_fac = dlt2 / cos(dlt1);

   /* rotation */
   snr = sin((double)rota * DRAD);
   csr = cos((double)rota * DRAD);
}

void utl_cart(double lon,double lat,double *x,double *y)
{
    double	tmp;
    double	tmp_x, tmp_y;

    lat = 60.0 * lat;
    lon = 60.0 * lon;

    tmp_x = lon - orig_lon;
    tmp_y = lat - orig_lat;

    tmp   = atan(DRLT * tan(DRAD * (lat+orig_lat)/120.0));
    tmp_x = (double)tmp_x * lon_fac * cos(tmp);    
    tmp_y = (double)tmp_y * lat_fac;
    
    /* rotation */
    *x = csr*tmp_x - snr*tmp_y;
    *y = csr*tmp_y + snr*tmp_x;

}


void utl_lonlat(double x,double y,double *lon,double *lat)
{
    double	tmp_x;
    double	tmp_y;
    double	tmp;

    tmp_x = snr*y + csr*x;
    tmp_y = csr*y - snr*x;

    tmp_y = tmp_y/lat_fac;
    tmp_y += orig_lat;

    tmp = atan(DRLT * tan(DRAD * (tmp_y+orig_lat)/120.0));
    tmp_x = tmp_x / (lon_fac * cos(tmp));
    tmp_x += orig_lon;
	

    *lon = tmp_x/60.0;
    *lat = tmp_y/60.0;
}

/*
#ifdef DEBUG
   
main()
{
    double	lon;
    double	lat;
    double	rota;
    double	slon,slat,x,y;

    rota = 0.0;

    printf("Enter Longitude and Latitude of the Origin (f, f)->\n");
    scanf("%f, %f",&lon,&lat);

    printf("Enter the rotation angle (f)->\n");
    scanf("%f",&rota);

    utl_orig(lon,lat,rota);

    printf("LAT_FAC: %lf km   LON_FAL: %lf km    SNR: %lf    CSR: %lf\n", lat_fac, lon_fac, snr, csr);

    printf("Enter station coordinates(Lon, Lat) ->\n");
    scanf("%f, %f", &slon, &slat);

    utl_cart(slon,slat,&x,&y);
    printf("Relative station coordinates-> x: %f    y: %f\n",x,y);
    
    do
    {
    printf("Enter station coordinates(Lon, Lat) ->\n");
    scanf("%f, %f", &slon, &slat);
    utl_cart(slon,slat,&x,&y);
    printf("Relative station coordinates-> x: %f   y: %f computed for(Lon: %f, Lat: %f)\n",x,y,slon,slat);
    } while (x>-1000.0);

}

#endif
*/
