#pragma once
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <math.h>



#define MAX_LENGTH 50
#define NEWTONRH_LEN 10
#define u32 uint32_t
#define u16 uint16_t
#define log10_2 0.30103

class apn
{
	public:
		u32 *data=nullptr;
		apn();
		//~apn();
		apn(u32 *initval);
		apn(u32 initval);
		void print();
		void print(u32 n);
		void printh();
		void printh(u32 n);
		void printb();
		void printb(u32 n);
		void printd();
		void printd(u32 n);
		void printdf(u32 comapos);
		uint32_t& operator[](u32 k);
		void operator=(apn v);
		void Rshift(u32 n);
		void Lshift(u32 n);
		bool isNull();
		void zero();
};

apn pi();

u32 leaddigitpos(apn& u);

apn offsetadd(apn u, u32 v, u32 n);
apn offsetaddhalf(apn u, u32 v, u32 n);
apn def_nrc1(u32 n);
apn def_nrc2(u32 n);

apn operator+(apn u, apn v);
apn operator-(apn u, apn v);
apn operator*(apn u, apn v);
apn operator/(apn u, apn v);
bool operator==(apn u, apn v);
apn invsqrt(apn u);


const uint32_t newton_rhapson_lut[28] = { 0 ,1 ,2 ,4 ,8 ,16 ,32 ,65 ,130 ,261 ,523 ,1046 ,2092 ,4185 ,8371 ,16742 ,33484 ,66968 ,133937 ,267875 ,535751 ,1071503 ,2143007 ,4286015 ,8572030 ,17144061 ,34288123 , 68576246 };
const uint32_t NEWTONRH_4817[2] = { 0x00000002,0xD2D2D2D2 };// 48/17 << 0+32 BAAAD NO YOU IDIOT (NEED 0.5<=D<2)
const uint32_t NEWTONRH_3217[2] = { 0x00000001,0xE1E1E1E1 }; //32/17 << 0+32 BAAAD NO YOU IDIOT
