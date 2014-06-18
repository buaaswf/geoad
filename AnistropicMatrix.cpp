#include "AnistropicMatrix.h"


AnistropicMatrix::AnistropicMatrix(void)
{
}


AnistropicMatrix::~AnistropicMatrix(void)
{

}
Eigen::Eigen()
{
	vector <double> eigenvalues(3);
	vector<vector <double>> eigenvectors(3,vector<double>(3));
}
bool AnistropicMatrix::EigenValuesAndEigenVectors( PIXTYPE inmatrix[],Eigen & eigen )
{
		int msize=3;
		vector<vector<double> > A(msize, vector<double>(msize));
		vector<vector<double> > V(msize, vector<double>(msize));
		
		
	for (int i=0;i<msize;i++)
	{
		for (int j=0 ; j <msize;j++)
		{
			A[i][j] = inmatrix[i*msize+j];
		}
	}






	for(int iv=0; iv < msize; iv++)
	{
			for(int iv2 = 0;iv2 < msize; iv2++)
			{
				if(iv==iv2){
					V[iv][iv2]=1;
				}
				else
				{
					V[iv][iv2]=0;
				}
			}
		}
		double *eigsv=new double[msize];
		for(int ieigsv=0;ieigsv<msize;ieigsv++)
		{
			eigsv[ieigsv]=0;
		}


		double epsl=0.0001;
		int maxt=10;
		int n=msize;

		int success=0; // 函数返回值
		double tao, t, cn, sn; // 临时变量

		
		
		double maxa;

		//写特征向量到C:\\tezhengxiangliang.txt
		FILE *fp2;
		if((fp2=fopen("F:\\tezhengxiangliang.txt","w"))==NULL)
		{

			printf("无法建立！");
			exit(0);
		}
		//写特征值到C:\\tezhengzhi.txt
		FILE *fp4;
		if((fp4=fopen("F:\\tezhengzhi.txt","w"))==NULL)
		{

			printf("无法建立！");
			exit(0);
		}


		for (int it=0; it< maxt; it++)
		{
			maxa=0;
			for (int p=0; p<n-1; p++)
			{
				for (int q=p+1; q<n; q++)
				{
					if (fabs(A[p][q])>maxa) // 记录非对角线元素最大值
					{
						maxa=fabs(A[p][q]);
					}
					//cout<<"fabs(A[p][q]"<<A[p][q]<<endl;
					if (fabs(A[p][q])>epsl) // 非对角线元素非0时才执行Jacobi变换
					{
						// 计算Givens旋转矩阵的重要元素:cos(theta), 即cn, sin(theta), 即sn
						
						tao=0.5*(A[q][q]-A[p][p])/A[p][q];

						if (tao>=0) // t为方程 t^2 + 2*t*tao - 1 = 0的根, 取绝对值较小的根为t
						{
							t=-tao+sqrt(1+tao*tao);
						} 
						else
						{
							t=-tao-sqrt(1+tao*tao);
						}
						cn=1/sqrt(1+t*t);
						sn=t*cn;

						// Givens旋转矩阵之转置左乘A, 即更新A的p行和q行
						for (int j=0; j<n; j++)
						{
							double apj=A[p][j];
							double aqj=A[q][j];
							A[p][j]=cn*apj-sn*aqj;
							A[q][j]=sn*apj+cn*aqj;
						}

						// Givens旋转矩阵右乘A, 即更新A的p列和q列
						for (int i=0; i<n; i++)
						{
							double aip=A[i][p];
							double aiq=A[i][q];
							A[i][p]=cn*aip-sn*aiq;
							A[i][q]=sn*aip+cn*aiq;
						}

						// 更新特征向量存储矩阵V, V=J0×J1×J2...×Jit, 所以每次只更新V的p, q两列
						for (int i2=0; i2<n; i2++)
						{
							double vip=V[i2][p];
							double viq=V[i2][q];
							V[i2][p]=cn*vip-sn*viq;
							V[i2][q]=sn*vip+cn*viq;
						}
					} 

				}  

			} 


			if (maxa<epsl) // 非对角线元素已小于收敛标准，迭代结束
			{
				// 特征值向量
				for (int j2=0; j2<n; j2++)
				{
					eigsv[j2]=A[j2][j2];
					//fprintf(fp2, "%f ",eigsv[j2]);
					//fprintf(fp2, "%s ","##");
				}
				// "对特征值向量从大到小进行排序, 并调整特征向量顺序 (直接插入法)"


				
				double* tmp=new double[n];
				for (int j=1; j<n; j++)
				{
					int i=j;
					double a=eigsv[j];
					for (int k=0; k<n; k++)
					{
						tmp[k]=V[k][j];
					}
					while(i>0 && eigsv[i-1]<a){
						eigsv[i]=eigsv[i-1];
						for (int k=0; k<n; k++)
						{
							V[k][i]=V[k][i-1];
						}
						i--;
					}
					eigsv[i]=a;
					
					for (int k2=0; k2<n; k2++)
					{
						V[k2][i]=tmp[k2];
						
					}    
				}
				for (int i=0;i<n;i++)
				{
					eigen.eigenvalues.push_back(eigsv[i]);
				}
				delete[] tmp;

				cout<<"Hello World!"<<"非对角线元素已小于收敛标准，迭代结束"<<endl;
				for(int ivc=0;ivc<msize;ivc++){
					for(int ivc2=0;ivc2<msize;ivc2++)
					{
						//cout<<V[ivc][ivc2]<<" ";
						fprintf(fp2, "%f%s",V[ivc][ivc2],"\n");
						
					}
					fprintf(fp2, "%s","#\n");
					//cout<<endl;
					eigen.eigenvectors=V;
				}
				fclose(fp2);

				for(int ieigsvc=0; ieigsvc < msize; ieigsvc++){

					fprintf(fp4, "%f%s",eigsv[ieigsvc],"\n");
					//cout<<endl;
				}
				fclose(fp2);
				return success=1;
			}
		}

		//---------------------------------------------------------------------------------------------------------
		//---------------------------------------------------------------------------------------------------------

		//jacbobi_loop(A,V,eigsv,0.0,100,10);
		//	printf("Hello World!\n"+c);
		printf("Hello World!\n");

		return success;


}

/**
 \brief	Calculates the jacobi matrix.

 \param [in,out]	src	If non-null, source for the.
 \param	jacobim		   	The jacobim.

 \return	true if it succeeds, false if it fails.
 */

bool AnistropicMatrix::ComputeJacobiMatrix( Raw *src,vector <Jacobi> jacobim )
{
	for (int i= 0; i< src->getXsize(); i++ )
	{
		for (int j = 0; j < src->getYsize(); j++ )
		{
			for (int k = 0; k < src->getZsize(); k++)
			{
				Jacobi * point = new  Jacobi();
				//9 data of matrix
				// IXX,IXY,IXZ
				// IYX,IYY,IYZ
				// IZX,IZY,IZZ
				// IXX
				point->jacobimarix.push_back(src->get(i+1,j,k)+src->get(i-1,j,k)-2*src->get(i,j,k));
				//IXY
				point->jacobimarix.push_back(src->get(i+1,j,k)+src->get(i-1,j,k)-2*src->get(i,j,k));
				jacobim.push_back(*point);

			}
		}
	}
	return true;
}


