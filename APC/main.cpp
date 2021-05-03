#include <stdlib.h>
#include "apn.h"


int main(int argc, char** argv)
{
	printf("ARBITRARY PRECISION COMPUTING\n\n");

	apn a;
	apn b;
	uint32_t c = 0xFFFFFFFF;
	
	for(u32 k=0;k<MAX_LENGTH;k++)//MAX_LENGTH
	{ 
		//a[k] = 0xFFFFFFFF;
		b[k] = 0x00000000;
	}
	//a[0] = 0xFFFFFFFF;
	
	a[5] = 0x00000001;
	a[4] = 0x80000000;
	
	printf("a: ");
	a.printh();
	printf("\n");



	printf("a: ");
	a.printdf(5*32);
	printf("\n");



	printf("invsqrt2: ");
	invsqrt(a).printdf(5);
	printf("\n");

	
	//b[0] = 0x00000007;
	//b[0] = 0xFFFFFFFF;

	//printf("\n%d\n", leaddigitpos(a));

	

	int y = 1;
	apn res;
	while (y--)
	{
		//res=(a/b);
		//free(res.data);
	}
	
	//printf("\n");
	free(a.data);
	free(b.data);
}