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
class Jacobi;
class AnistropicMatrix
{
public:
	AnistropicMatrix(void);
	~AnistropicMatrix(void);
	bool EigenValuesAndEigenVectors(PIXTYPE [],Eigen &);
	bool ComputeJacobiMatrix(Raw *src,vector<Jacobi> jacobim);

};

/**
 \brief	Eigen.store the eigen for every point
 */

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

/**
 \brief	Jacobi.store the jacobi matrix for every point
 */

class Jacobi
{
	public:
		vector <double> jacobimarix;
		
};

