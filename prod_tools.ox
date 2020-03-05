#include <oxstd.h>
#include "prod_tools.h"

translog(const data, const flag, const asNames, ...)
/* Calculates the quadratic form of a given set of variables
 *  data: mat of vars (vars are columns)
 *  flag: int to adjust options (0: nothing, 1 negative of all vars)
 *  asNames: array of names to manipulate
 *  (Optional) vIdx: vector with col indices from which to construct the quadratic form
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
             // check if any constants in data
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
			println("The whole X-Matrix has been multiplied by -1!");
//            for (i=0;i<sizeof(asNames[0]);i++)
//                asNames[0][i]=sprint("-",asNames[0][i]);
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

fPlotPartial(const sXname, const oM, const m_iTl,const dWidth, ...){
/* Plots the dependent variable against the named sXname variable, also for flexible functional form
 *  sXname: Name of variable as used in .Select(X_VAR, ...)
 *  oM: instance of database (or derived) class, eg fob=new Sfamb();
 *  m_iTl: int to indicate function type
 *    0: only first order terms
 *    1: translog on all X_Vars
 *    >2: translog on first m_iTL vars
 *    <0: as above, negative sign indicates that all regressors are negative (output distance functions)
 * dWidth: between 0 and 1; middle dWidth% of sXName are plotted on the abscissa
 * bLn: TRUE for logarithmic functional forms (default: TRUE)
*/
//get Parameter names
decl vMeanVars, vPar, i, vIdx, maxParIdx, sTmp, sTmp1, sTmp2, bSign=1,
	iOffset=0, asNames, asNamesX, asNameY, args = va_arglist(), bLn=TRUE;

  if (sizeof(args) > 1){
		eprint("Too many arguments"); exit(9);
	}

  if (sizeof(args) > 0)
		bLn = args[0];
 

  asNames=oM.GetParNames();
  maxParIdx = strfind(asNames,"ln{\sigma_v}");
  // some error checking
  if (maxParIdx == -1){
    eprint("Error with finding sigma_v in the list of parameter names, exiting...");
	exit(1);
  }
  if (strfind(asNames,sXname) == -1){
    eprint("Error with finding ",sXname," in the list of parameter names, exiting...");
	exit(1);
  }

  vPar=oM.GetPar()[:maxParIdx-1];
  oM.GetGroupNames(0,&asNamesX);
  oM.GetGroupNames(1,&asNameY);
  if (m_iTl < 0){
    i = fabs(m_iTl);
	i = sizeof(asNamesX) - (i+1)*(i+2)/2	+i;
  	asNamesX=asNames[:i];
  }

  // if squares and interactions are not constructed from all vars, an offset is required later
  if (m_iTl < 0) bSign=-1;
  if (bSign * m_iTl > 1) iOffset = sizeof(asNamesX)-(bSign*m_iTl)-1;
//  print(iOffset); exit(0);
  decl vTmp=vPar[0];
  //first order terms
  for (i=1;i<sizeof(asNamesX);++i)
  	 vTmp += asNames[i] == sXname
	 	? oM.GetVar(asNames[i]) * vPar[i]
		: meanc(oM.GetVar(asNames[i])) * vPar[i];

  if ((bSign*m_iTl) >0 ){
  // square terms
	for (i = sizeof(asNamesX); i < 2*sizeof(asNamesX)-iOffset-1; ++i){
	  sTmp=asNames[i];
	  sTmp= sTmp[strfind(sTmp,".5*")+3:strfind(sTmp,"^")-1];println(sTmp);
  	  vTmp += sTmp == sXname
	 	? 0.5 * oM.GetVar(sTmp).^2 * vPar[i]
		: 0.5 * meanc(oM.GetVar(sTmp)).^2 * vPar[i];
	}
  // interaction terms    
  	for (i = 2*sizeof(asNamesX)-iOffset-1; i < maxParIdx; ++i){
      sTmp=asNames[i];//println(sTmp);
	  sTmp1 =sTmp[:strfind(sTmp,"*")-1];
	  sTmp2 = sTmp[strfind(sTmp,"*")+1:];
	  if (sTmp1 == sXname)
	    vTmp += oM.GetVar(sTmp1).* meanc(oM.GetVar(sTmp2)) * vPar[i];
	  else if (sTmp2 == sXname)
		  vTmp +=	oM.GetVar(sTmp2).* meanc(oM.GetVar(sTmp1)) * vPar[i];
	    else vTmp += meanc(oM.GetVar(sTmp1)).* meanc(oM.GetVar(sTmp2)) * vPar[i];
    }
  }
//	print(meanc(oM.GetVar("lnQ")~vTmp));
	// for smooth lines, plot data need to be sorted ...
	decl mOut;
	if (bLn) mOut = exp(vTmp) ~ exp(bSign*oM.GetVar(sXname));
	  else
	  	mOut = vTmp~(bSign*oM.GetVar(sXname));
	if (strfind(sXname,asNameY[0])> -1) mOut[][1] = mOut[][1] .* mOut[][0];
	mOut = sortbyc(mOut,1);
	//trim; plot only the middle dWidth %
	print(0.5*(1-dWidth)*rows(mOut)," to ", 0.5*(1+dWidth)*rows(mOut));
	mOut = mOut[0.5*(1-dWidth)*rows(mOut):0.5*(1+dWidth)*rows(mOut)][];
	println(" equal to ", rows(mOut), " observations");
	DrawXMatrix(0,mOut[][0]',"Output",mOut[][1]',sXname,2,2);
	DrawTitle(0,"Partial frontier plot for exp("+sXname+")");
	ShowDrawWindow();
//  print(asX);
return 1;
}
