#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
float *real, *imag, amp;
float rmean, imean, rsig, isig;
int i, n;

long seed = 0;
float fsigma = 1.0;
float fmean = 0.0;

scanf("%d",&n);

real = (float *)check_malloc(n*sizeof(float));
imag = (float *)check_malloc(n*sizeof(float));

amp = 0.0;
for(i=0;i<n;i++)
   {
   real[i] = gaussian_rand(&fsigma,&fmean,&seed);
   imag[i] = gaussian_rand(&fsigma,&fmean,&seed);

   rmean = rmean + real[i];
   imean = imean + imag[i];

   amp = amp + (real[i]*real[i] + imag[i]*imag[i]);
   }

rmean = rmean/(float)(n);
imean = imean/(float)(n);

rsig = 0.0;
isig = 0.0;
for(i=0;i<n;i++)
   {
   rsig = rsig + (real[i] - rmean)*(real[i] - rmean);
   isig = isig + (imag[i] - imean)*(imag[i] - imean);
   }
rsig = sqrt(rsig/(float)(n-1));
isig = sqrt(isig/(float)(n-1));

printf("rmean= %f, rsig= %f\n",rmean,rsig);
printf("imean= %f, isig= %f\n",imean,isig);

printf("avg. amp= %f\n",sqrt(amp/(float)(n)));
}
