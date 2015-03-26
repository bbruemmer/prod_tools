#include <oxstd.h>
#include "prod_tools.h"

translog(const data, const flag, const asNames, ...)
/* Calculates the quadratic form of a given set of variables
 *  data: mat of vars (vars are columns)
 *  vIdx: index of cols of which the quadratic form is to be constructed
 *  flag: int to adjust options (0: nothing, 1 negative of all vars)
 * asNames: array of names to manipulate
*/
{
        decl out, cX, cTL, asVar, i, j, vIdx, args = va_arglist();
        cX = columns(data);
            //check for constants that are in vIdx
            // check if vIdx is given
        if (sizeof(args) > 0)
			vIdx = args[0];
		else
			vIdx = range(0,(cX-1));
            // check if names are given
        if (any(!varc(data[][vIdx]))){
            println("Constant in data matrix (col) ",
                vecindex(varc(data[][vIdx]),0));
            return 0;
        }
        cTL = cX+columns(vIdx)*(columns(vIdx)+1)/2;
        asVar = new array[cX];
        if (sizeof(asNames[0]) == 0) // empty array
        {
            for (i=0;i<cX;i++) asVar[i]=sprint("Var",i+1);
        }
        else asVar=asNames[0];
        asNames[0] = new array[cX];
        for (i=0;i<cX;i++) asNames[0][i]=sprint(asVar[i]);
        for (i=0;i<columns(vIdx);i++)
            asNames[0]~=sprint(".5*",asVar[vIdx[i]],"^2");
        for (i=0;i<columns(vIdx);i++)
            for (j = i; j < columns(vIdx)-1; j++)
                asNames[0] ~= sprint(asVar[vIdx[i]],"*",asVar[vIdx[j+1]]);

        out=data~0.5*(data[][vIdx] .* data[][vIdx]);
        for (i = 0 ; i < columns(vIdx)-1 ; i++)
        {
          out~=(data[][vIdx[i+1]:vIdx[columns(vIdx)-1]].* data[][vIdx[i]]);
        }
        //print("done.\n");
        if (flag){
            out= -out;
            for (i=0;i<sizeof(asNames[0]);i++)
                asNames[0][i]=sprint("-",asNames[0][i]);
        }
        return out;
}

db_expand(const asNames){
    decl asExpanded ={},i;
//  asExpanded = new array(3*(sizeof(asNames)));
    for (i=0;i<(sizeof(asNames));++i)
        asExpanded ~= {asNames[i]}~{0,0};
    return asExpanded;
}
	  