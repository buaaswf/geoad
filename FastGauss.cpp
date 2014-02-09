#include "FastGauss.h"


FastGauss::FastGauss(void)
{
	windowsize = 3;
}


FastGauss::~FastGauss(void)
{
}
Raw * FastGauss:: fg(Raw *data,Raw *ret)
{
	Raw *src= new Raw(*ret,true);
	Raw *temp = new Raw(*data);
	for (int i = 1; i <src->getXsize()-1; i++)
	{
		for ( int j =0; j < src->getYsize() ; j++)
		{
			for ( int k = 0; k <src->getZsize() ; k++)
			{
				PIXTYPE val = temp->get(i-1,j,k)+2*temp->get(i,j,k)+temp->get(i+1,j,k);
				src->put(i,j,k,val*0.25);
			}
		}
	}
	temp = new Raw(*src);
	for (int i = 0; i <src->getXsize(); i++)
	{
		for ( int j = 1; j < src->getYsize()-1 ; j++)
		{
			for ( int k = 0; k <src->getZsize() ; k++)
			{
				PIXTYPE val = temp->get(i,j-1,k)+2*temp->get(i,j,k)+temp->get(i,j+1,k);
				src->put(i,j,k,val*0.25);
			}
		}
	}
	temp = new Raw(*src);
	for (int i = 0; i <src->getXsize(); i++)
	{
		for ( int j =0; j < src->getYsize() ; j++)
		{
			for ( int k = 1; k <src->getZsize() -1; k++)
			{
				PIXTYPE val = temp->get(i,j,k-1)+2*temp->get(i,j,k)+temp->get(i,j,k+1);
				src->put(i,j,k,val*0.25);
			}
		}
	}
	return src;


}