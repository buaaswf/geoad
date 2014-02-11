#pragma once
#include "vol_math_filter_Interface.h"
#include "vol_math_RawImage.h"
#include <stdio.h>
#include <math.h> 
#include <conio.h>
#include <stdlib.h>
#include<iostream>
#include<vector>
using namespace std;
//using std::vector;
class Eigen;
class AnistropicMatrix
{
public:
	AnistropicMatrix(void);
	~AnistropicMatrix(void);
	bool EigenValuesAndEigenVectors(PIXTYPE [],Eigen &);

};

class Eigen
{
public:
	//int msize;
	vector <double> eigenvalues;
	vector<vector <double>> eigenvectors;
public:
	Eigen(void);
	~Eigen();
};
class Jacobi
{
	public:
		vector <double> jacobim;

};

