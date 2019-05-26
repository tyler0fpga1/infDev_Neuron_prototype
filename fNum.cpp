

#pragma  once
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <chrono>
#include <random>
#include "fNum.h"
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#pragma once
double  seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937 Num(seed);
volatile float size[6]; //(this is really 5:0... because c++.)
unsigned int fNumCounter = 0;
volatile double compVal = 0;
unsigned int countToTen = 0;
/****NUMBER GENERATION****************************************/

float fRand() {

	fNumCounter++;

	size[0] = Num() - fNumCounter; //random Number 1. Basically a random seed
	size[1] = 0;
	size[2] = 0;
	size[3] = 0;
	size[4] = 0;
	size[5] = 0; // return val (-1 to +1)

	while (size[5] == 0) //dont return 0.
	{
		countToTen++;
		if (countToTen >= 10) {
			countToTen = 0;
			size[1] = 0;
		}
		size[1] = +Num();
		if (size[0] > size[1])
		{
			size[2] = size[0] / size[1];
			size[1] = size[0];
		}
		size[1] = Num();
		if (size[0] > size[1])
		{
			size[3] = size[0] / size[1];
			size[1] = size[0];
		}
		//asign pos/neg
		if ((size[0] == size[1])) //this will only be true after both if statements above run
		{
			compVal = Num();
			if (size[2] > compVal) { size[4] = size[3] / size[2]; }
			if (size[2] < compVal) { size[4] = size[2] / size[3]; }
			if ((size[4] > -1) & (size[4] < 1) & (size[4] != 0))
			{
				while (size[4] < .01) { size[4] = size[4] * 11; } //scaler need to adda const typedef var to it
				//for (long unsigned int counter = 0x0; counter > 0xFFFF; counter++) 
				//.95 is the min and max this func returns. need to add const type def var instead of .95
				if (compVal < size[0]) { if ((size[4] < .95) & (size[4] > -.95)) size[5] = size[4]; }
				if (compVal > size[0]) { if ((size[4] < .95) & (size[4] > -.95)) size[5] = ((-1) * size[4]); }
			}

		}

		size[2] = 0;
		size[3] = 0;
		size[4] = 0;

	}
	//found the perfect random number. Reset and return number.
	size[1] = 0;
	return size[5];
};
//....................

