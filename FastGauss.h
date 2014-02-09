#pragma once

/************************************************************************/
/*fast gauss filter will be improved in the followed two aspects  
1\no need to compute the number use the integer filter window 1 2 1 
2\one dimensional for x,y,z directions 
*/                                                                   
/************************************************************************/
#include "vol_math_filter_Interface.h"
#include "vol_math_RawImage.h"
class FastGauss
{
private:
	Raw *src;
	size_t windowsize ;
public:
	FastGauss(void);
	~FastGauss(void);
	Raw* fg(Raw *data,Raw *ret);

};

