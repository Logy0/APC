#include "apn.h"
#pragma warning(disable : 4996)

apn::apn()
{
	data = (uint32_t*)malloc(sizeof(uint32_t) * MAX_LENGTH);
	if (data == NULL) exit(1);

	for (u32 k = 0;k < MAX_LENGTH;k++)
	{
		this->data[k]=0;
	}
}
apn::apn(u32* initval)
{
	data = (uint32_t*)malloc(sizeof(uint32_t) * MAX_LENGTH);
	if (data == NULL) exit(1);
	for (u32 k = 0;k < MAX_LENGTH;k++)
	{
		data[k] = initval[k];
	}
}
apn::apn(u32 initval)
{
	data = (uint32_t*)malloc(sizeof(uint32_t) * MAX_LENGTH);
	if (data == NULL) exit(1);
	for (u32 k = 0;k < MAX_LENGTH;k++)
	{
		data[k] = initval;
	}
}
/*

apn::~apn()
{
	//std::cout << data << std::endl;
	//free(data);
	
}
*/

bool apn::isNull()
{
	u32 k = leaddigitpos(*this);
	if(k>0 || data[0]!=0)
		return false;
	return true;
}
void apn::operator=(apn v)
{
	for (u32 k = 0;k < MAX_LENGTH;k++)
	{
		data[k] = v[k];
	}
}
void apn::zero()
{
	for (u32 k = 0;k < MAX_LENGTH;k++)
		data[k] = 0;
}
void apn::print() //TODO: fix printf uint32_t limit
{

	for (u32 k = 0;k < MAX_LENGTH;k++)
	{
		printf("%0*u ",10, this->data[MAX_LENGTH-k-1]);
	}
	printf("\n");
}
void apn::print(u32 n)
{
	u32 p=std::min(n,(u32)MAX_LENGTH);
	for (int k = 0;k < p;k++)
	{
		printf("%0*u ",10 ,this->data[p - k - 1]);
	}
	printf("\n");
}
void apn::printh()
{
	for (int k = 0;k < MAX_LENGTH;k++)
	{
		printf("%0*X ",8 ,this->data[MAX_LENGTH-k-1]);
	}
	printf("\n");
}
void apn::printb()
{
	for (int k = 0;k < MAX_LENGTH;k++)
	{
		for (u16 i = 32;i > 0;i--)printf("%d", (this->data[MAX_LENGTH - k - 1] & (1 << (i-1)))>>(i-1));
		printf(" ");
	}
	printf("\n");
}void apn::printb(u32 n)
{
	u32 p = std::min(n, (u32)MAX_LENGTH);
	for (int k = 0;k < p;k++)
	{
		for (u16 i = 32;i > 0;i--)printf("%d", (this->data[MAX_LENGTH - k - 1] & (1 << (i-1))) >> (i-1));
		printf(" ");
	}
	printf("\n");
}
void apn::printh(u32 n)
{
	u32 p = std::min(n, (u32)MAX_LENGTH);
	for (int k = 0;k < p;k++)
	{
		printf("%0*X ",8 ,this->data[p-k-1]);
	}
	printf("\n");
}
/*
void apn::printd()
{
	for (u32 k = 0;k <MAX_LENGTH;k++)// MAX_LENGTH
	{
		u32 num = data[k];

		u32 temp1,temp2=0;
		bool flag = true;
		u32 res = 0;
		printf("te1:\n");
		while (flag)
		{
			temp1 = num & 0x0000FFFF; //low 16bits
			temp2 = (num & 0xFFFF0000) >> 16; //high 16 bits

			for (u32 i = 0;i < 16;i++)
			{
				temp2 <<= 1;
				if ((temp2 & 0xFFFF0000) >= 0xA0000)
				{
					temp2 -= 0xA0000;
					res |= 1 << (31 - i);
				}
			}

			temp1 += temp2 & 0xFFFF0000; //include remainder from high bits

			for (u32 i = 0;i < 16;i++)
			{
				temp1 <<= 1;
				if ((temp1 & 0xFFFF0000) >= 0xA0000)
				{
					temp1 -= 0xA0000;
					res |= 1 << (15 - i);
				}

			}

			temp1 >>= 16;

			printf("%d\n", temp1);
			if (res == 0)
				flag = false;

			num = res;
			res = 0;
		}

	}
}
*/
void apn::printd()
{
	apn divd(this->data);
	const int mlen = MAX_LENGTH * log10_2 * 32+1;
	char *output;
	output = (char*)malloc(sizeof(char) * mlen);
	if (output == NULL) exit(1);
	apn tempdiv;
	u32 step = mlen / 20;
	u32 num;
	bool flag = true;
	u32 y = 0;
	while (flag)
	{ 
		tempdiv.zero();
		if(0&& y%step==0)
			printf("ite:%d\n", y);
		u32 temp1 = 0; //stores low bits
		u32 temp2 = 0; //stores high bits
		for (u32 k = 0;k < MAX_LENGTH;k++)
		{
			num = divd[MAX_LENGTH - k - 1];

			temp2 = (num & 0xFFFF0000) >> 16; //high 16 bits
			temp2 += (temp1 & 0xFFFF0000); //add previous remainder (which should be less than 0xA)
			temp1 = num & 0x0000FFFF; //low 16bits
			for (u32 i = 0;i < 16;i++)
			{
				temp2 <<= 1;
				if ((temp2 & 0xFFFF0000) >= 0xA0000)
				{
					temp2 -= 0xA0000;
					tempdiv[MAX_LENGTH - k - 1] |= 1 << (31 - i);
				}
			}

			temp1 += temp2 & 0xFFFF0000; //include remainder from high bits

			for (u32 i = 0;i < 16;i++)
			{
				temp1 <<= 1;
				if ((temp1 & 0xFFFF0000) >= 0xA0000)
				{
					temp1 -= 0xA0000;
					tempdiv[MAX_LENGTH - k - 1] |= 1 << (15 - i);
				}
			}
		}

		output[mlen-y-1] = 48 + (temp1 >> 16);
		if (tempdiv.isNull())
			flag = false;
		divd = tempdiv;
		y++;
	}

	output += mlen - y;
	if (y < 500)
	{
		for (u32 k = 0;k < y;k++)
			std::cout << output[k];
	}
	else
	{
		FILE* f = fopen("bignum.txt", "wb");
		fwrite(output, sizeof(char), y, f);
		fclose(f);
	}
}
void apn::printd(u32 n)
{
	apn divd(this->data);
	const u32 mlen = MAX_LENGTH * log10_2 * 32 + 1;
	char* output;
	output = (char*)malloc(sizeof(char) * mlen);
	if (output == NULL) exit(1);
	u32 leadpos = leaddigitpos(*this) / 32+1;  //------------------
	u32 step = mlen / 20;
	u32 temp1, temp2,num;
	apn tempdiv;
	bool flag = true;
	u32 y = 0;
	//printf("\nmlen:%d\n", mlen);
	while (flag && y<n)
	{
		tempdiv.zero();
		if (0 && y % step == 0)
			printf("ite:%d\n", y);
		temp1 = 0; //stores low bits
		temp2 = 0; //stores high bits
		for (u32 k = 0;k < MAX_LENGTH;k++)
		{
			num = divd[MAX_LENGTH - k - 1];

			temp2 = (num & 0xFFFF0000) >> 16; //high 16 bits
			temp2 += (temp1 & 0xFFFF0000); //add previous remainder (which should be less than 0xA)
			temp1 = num & 0x0000FFFF; //low 16bits
			for (u32 i = 0;i < 16;i++)
			{
				temp2 <<= 1;
				if ((temp2 & 0xFFFF0000) >= 0xA0000)
				{
					temp2 -= 0xA0000;
					tempdiv[MAX_LENGTH - k - 1] |= 1 << (31 - i);
				}
			}

			temp1 += temp2 & 0xFFFF0000; //include remainder from high bits

			for (u32 i = 0;i < 16;i++)
			{
				temp1 <<= 1;
				if ((temp1 & 0xFFFF0000) >= 0xA0000)
				{
					temp1 -= 0xA0000;
					tempdiv[MAX_LENGTH - k - 1] |= 1 << (15 - i);
				}
			}
		}

		output[mlen - y - 1] = 48 + (temp1 >> 16);
		if (tempdiv.isNull())
			flag = false;
		divd = tempdiv;
		y++;

	}

	output += mlen - y;
	if (y < 500)
	{
		for (u32 k = 0;k < y;k++)
			std::cout << output[k];
	}
	else
	{
		FILE* f = fopen("bignum.txt", "wb");
		fwrite(output, sizeof(char), y, f);
		fclose(f);
	}
	
}


void apn::printdf(u32 comapos)
{
	apn temp3;
	temp3[0] = 0xA;
	
	if (comapos < 0 || comapos>32 * MAX_LENGTH - 1)exit(1);
	u32 index = comapos / 32;
	u32 bit = comapos % 32;//0-31
	//printf("index: %u\n", index);
	//printf("bit: %u\n", bit);

	apn temp1(this->data);//integer part
	apn temp2(this->data);//float part
	for (u32 k = 0;k < MAX_LENGTH - 1;k++) //not so efficient
	{
		temp1[k] >>= bit; //LSB format
		temp1[k] |= (((1 << bit) - 1) & temp1[k + 1]) << (32 - bit);
	}

	temp1[MAX_LENGTH - 1] >>= bit; //LSB format
	
	temp1.Rshift(index);
	/*
	printf("temp1:");
	temp1.printh();	
	*/
	for (u32 k = index+1 ;k < MAX_LENGTH;k++) temp2[k] = 0;
	temp2[index] &= ((1 << bit) - 1);
	//printf("temp2:");
	//temp2.printh();
	//printf("\n");
	temp1.printd();
	printf(".");
	u32 res = 0;
	if (index == MAX_LENGTH - 1)
	{
		u32 excess = temp2[0];
		//temp2.Rshift(1);

		while (!temp2.isNull())
		{
			temp2 = temp2 * temp3;
			temp2.printh();
			res = ((~((1 << bit) - 1)) & temp2[index]);
			printf("\n%u\n", res>>bit);
			temp2[index] &= ((1 << bit) - 1);
		}
	}
	else
	{
		while (!temp2.isNull())
		{
			temp2 = temp2 * temp3;
			res = (((1 << bit) - 1) & temp2[index + 1]) | ((~((1 << bit) - 1)) & temp2[index]);
			printf("%u", res>>bit);
			temp2[index] &= ((1<< bit) - 1);
			temp2[index+1] &= ~((1<<bit) - 1);			
		}
	}
}
u32 longdiv_10(apn u) //TO DO
{
	return 0;
} //TO DO

uint32_t leaddigitpos(apn& u)
{
	uint32_t pos = MAX_LENGTH - 1;
	while(pos>0 && u[pos] == 0)
	{
		pos--;
	}
	uint32_t i = 31;
	while ((u[pos] & (1 << i)) == 0 && i > 0)
		i--;
	return 32 * pos + i;
}

apn operator+(apn u, apn v)
{
	apn res(u);
	bool carry_next = 0;
	bool cr, cv, cc;
	for (int k = 0;k < MAX_LENGTH;k++)
	{
		cr = res.data[k] >> 31;
		cv = v.data[k] >> 31;
		res.data[k]&= ~(1 << 31);
		v.data[k]&=~(1<<31);
		res.data[k] += uint32_t(v.data[k]+carry_next);
		cc = res.data[k] >> 31;
		carry_next = cr&&cv || cr&&cc || cc&&cv;
		if(k==MAX_LENGTH-1 && carry_next)
			printf("+ operator OVERFLOW\n");
		if (cr && cv && (!cc) || cr && cc && (!cv) || cc && cv && (!cr))
		{
			res.data[k] ^= 1 << 31;
		}	
		if (cr && !(cv || cc) || cc && !(cr || cv) || cv && !(cc || cr))
		{
			res.data[k] |= 1 << 31;
		}		
	}
	return res;
}
apn operator-(apn u, apn v) //u-v
{
	apn res;
	bool borrow_next = 0;

	for (int k = 0;k < MAX_LENGTH;k++)
	{

		if (u[k] > 0)
		{
			if (u[k] - borrow_next >= v[k])
			{
				res[k] = u[k] - borrow_next - v[k];
				borrow_next = 0;
			}
			else
			{
				res[k] = u[k] - borrow_next - v[k];
				borrow_next = 1;
			}
		}
		else
		{
			if (borrow_next != 0 || v[k] != 0)
			{
				res[k] = u[k] - borrow_next - v[k];
				borrow_next = 1;
			}
			else
			{
				borrow_next = 0;
			}			
		}
		if (k == (MAX_LENGTH - 1) && borrow_next)
			printf("- operator UNDERFLOW\n");
	}
	return res;

}
//pass pointer *u reference to avoid many copies of res
apn operator*(apn u, apn v)
{
	apn res;
	uint32_t highu,highv, lowu,lowv; //MAYBE LOOK INTO USING UINT16_T and casting the products to UINT32_T
	uint32_t p1,p2, p3,p4; 
	for (int k = 0;k < MAX_LENGTH;k++)
	{
		
		for (int i = 0;i < MAX_LENGTH;i++)
		{
			highu = (u.data[k] & 0xFFFF0000) >> 16;
			highv = (v.data[i] & 0xFFFF0000) >> 16;
			lowu = (u.data[k] & 0x0000FFFF);
			lowv = (v.data[i] & 0x0000FFFF);
			p1 = lowu * lowv;
			p2 = lowu * highv;
			p3 = highu * lowv;
			p4 = highu * highv;

			if (k + i < MAX_LENGTH)
			{
				res = offsetadd(res, p1, k + i);
			}
			else
			{
				if (p1 != 0)
				{
					printf("* operator OVERFLOW %d\n", k + i);
					//throw std::overflow_error("MULTIPLICATION OV");
					exit(1);
				}
			}
			if ((k + i + 1 < MAX_LENGTH))
			{
				res = offsetaddhalf(res, p2, k + i); // offsetaddhalf shifts index by 1 -> k+i+1
				res = offsetaddhalf(res, p3, k + i);
				res = offsetadd(res, p4, k + i + 1);
			}
			else
			{
				if ((p2 & 0xFFFF0000) >> 16 != 0 || (p3 & 0xFFFF0000) >> 16 != 0 || p4 != 0)
				{
					printf("* operator OVERFLOW %d\n", k + i+1);
					//throw std::overflow_error("MULTIPLICATION OV");
					exit(1);
				}
			}			
		}
	}
	return res;
}
uint32_t& apn::operator[](u32 k)
{
	return data[k];
}
bool operator==(apn u, apn v)
{
	for (u32 k = 0;k < MAX_LENGTH;k++)
	{
		if (u[k] != v[k])
			return false;
	}
	return true;
}
void apn::Rshift(u32 n)
{

	for (u32 k = n;k < MAX_LENGTH;k++)
	{
		data[k - n] = data[k];
	}
	for (u32 k = MAX_LENGTH-n;k < MAX_LENGTH;k++)
	{
		data[k] = 0;
	}
}
void apn::Lshift(u32 n)
{
	//later
} //TO DO

apn offsetadd(apn u, uint32_t v, u32 n)
{
	apn res(u);
	bool carry_next = 0;
	bool cr, cv, cc;
	int k = n;
	if (n >= MAX_LENGTH || n < 0)
		exit(1);
	do {
		cr = res.data[k] >> 31;
		cv = v >> 31;
		res.data[k] &= ~(1 << 31);
		v &= ~(1 << 31);
		res.data[k] += uint32_t(v + carry_next);
		cc = res.data[k] >> 31;
		carry_next = cr && cv || cr && cc || cc && cv;

		if (cr && cv && (!cc) || cr && cc && (!cv) || cc && cv && (!cr))
		{
			res.data[k] ^= 1 << 31;
		}
		if (cr && !(cv || cc) || cc && !(cr || cv) || cv && !(cc || cr))
		{
			res.data[k] |= 1 << 31;
		}
		v = 0;
		k++;
	} while (carry_next != 0 and k < MAX_LENGTH);

	return res;
}
apn offsetaddhalf(apn u, uint32_t v, u32 n)//separates a uint32_t into 2 highbytes and 2 lowbytes then add them to u;
{
	uint32_t high = (v & 0xFFFF0000) >> 16;
	uint32_t low = (v & 0x0000FFFF) << 16;
	return offsetadd(offsetadd(u, low, n), high, n + 1);
}

apn def_nrc1(u32 n)//n is an index of data
{
	if (n >= MAX_LENGTH || n < 0)
		exit(1);
	apn nrc1;
	for (u32 k = 0;k < n;k++)
	{
		nrc1[k] = NEWTONRH_4817[1];		
	}
	nrc1[n] = NEWTONRH_4817[0];
	return nrc1;
}
apn def_nrc2(u32 n)//n is an index of data
{
	if (n >= MAX_LENGTH || n < 0)
		exit(1);
	apn nrc2;
	for (u32 k = 0;k < n ;k++)
	{
		nrc2[k] = NEWTONRH_3217[1];
	}
	nrc2[n] = NEWTONRH_3217[0];
	return nrc2;
}
apn operator/(apn u, apn v)
{
	u32 k = ceil(leaddigitpos(v)/32.0);
	printf("Shift factor :%d\n", k);
	apn res;
	apn temp;
	//apn nrc1 = def_nrc1(2);
	//apn nrc2 = def_nrc2(2);
	apn one;
	one[k] = 2;
	
	printf("v: ");
	v.printh(10);
	uint32_t idx = 1;
	printf("one:   ");
	one.printh(10);

	res[5] = 0xFFFF0000;
	//res[1] = 0x00000000;
	printf("res0:  ");
	res.printh(10);
	
	while ( idx<200 && !(temp==res))//NEWTONRH_LEN && MAX_LENGTH * 32>newton_rhapson_lut[5]
	{
		temp = res;
		
		res =  res*(one - res*v);
		res.Rshift(k);
		printf("res%d: ",idx);
		//res.printh(10);
		printf("diff%d: ",idx);
		(temp - res).printh(10);
		
		idx++;
	}

	/*
	for (int k = 0;k < idx;k++)
	{
		void;
	}
	*/
	return res;
}
apn invsqrt(apn u)
{
	u32 k = ceil(leaddigitpos(u) / 32.0);
	//printf("Shift factor :%d\n", k);
	apn res;
	apn temp;
	apn one;
	one[k] = 1;
	one[k-1] = 0xA0000000;
	
	uint32_t idx = 1;

	/*printf("one: ");
	one.printdf(k * 32);
	printf("\n");
	*/
	res[5] = 0xFFFF0000;

	/*printf("res0: ");
	res.printh(10);
	printf("\n");
	*/
	while (idx < 200 && !(temp == res))//NEWTONRH_LEN && MAX_LENGTH * 32>newton_rhapson_lut[5]
	{
		temp = res;

		res = res * (one - res* res);
		res.Rshift(k);
		printf("res%d: ", idx);
		//res.printh(10);
		printf("diff%d: ", idx);
		(temp - res).printh(10);

		idx++;
	}

	/*
	for (int k = 0;k < idx;k++)
	{
		void;
	}
	*/
	return res;
}

apn pi()
{
	apn pi;


}