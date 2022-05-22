#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
struct MDstruct
{   int         iRANK,              		iSIZE;
    int       	iF_BASEelem,        		iF_ADDelem,				iB_BASEelem,       			iB_ADDelem;  
    int       	*ivF_OFFSET,        		*ivF_START,				*ivB_OFFSET,        		*ivB_START;
	int         *ivF_ModOFFs,				*ivB_ModOFFs;

    int       	iRunNum,            		iSeedStart,				iUseProp,           		iUsePSeis;
    int       	iFPNum,             		iBPNum,					iFVNum,             		iBVNum; 
    int       	iFSegNum,          		 	iBSegNum,				iFUsdGrd,					iAlsoBoundForPostSeis;
    int       	iMaxIterat,      			iEQcntr,      			iMaxMRFlgth,				iMaxSTFlgth;				
    int      	iChgBtwEQs,					iGlobTTmax,				iWritePos,					iUseVpVs;
   
    float     	fLegLgth,          			fUnitSlip,				fDeltT,             		fCutStrss;
	float 		fFltLegs,					fBndLegs,               fMinMag4Prop,				fCutStFrac;

	float       fAftrSlipTime,				fISeisStep,         	fRecLgth,					fHealFact;
	float       fDeepRelaxTime,				fPSeis_Step;

    float       fg,							fMedDense,     			fAddNrmStrs,				fVp;
    float       fPoisson,           		fLambda,        		fShearMod,					fVs;
    float       fISeisTStp,         		fTimeYears,             fVpVsRatio,					fMeanStiffness;   

    char        cInputName[512];
};
//------------------------------------------------------------------
struct SGstruct
{	float     	*fvSlipRate,				*fvSlipRake,			*fvStrsRate;
	int         *ivFricLawUSED;
    float       *fvRefStaFr,				*fvRefDynFr,			*fvRefCritDist;
	float       *fvRefStaFr_vari,			*fvRefDynFr_vari,		*fvRefCritDist_vari;       
};
//------------------------------------------------------------------
struct TRstruct
{	int     	*ivFL_StabT,				*ivFL_Activated,		*ivFL_Ptch_t0,				*ivFL_FricLaw;	

    float    	*fvFL_RefNrmStrs,			*fvFL_Area;
	float		*fvFL_SelfStiffStk,			*fvFL_SelfStiffDip,		*fvFL_MeanSelfStiff;
	float		*fvBL_SelfStiffStk,			*fvBL_SelfStiffDip,		*fvBL_SelfStiffOpn;

    float    	*fvFL_RefStaFric,           *fvFL_RefDynFric,    	*fvFL_RefDcVal;
    float    	*fvFL_RefStaFric_vari,      *fvFL_RefDynFric_vari,	*fvFL_RefDcVal_vari;
    float    	*fvFL_StaFric,              *fvFL_DynFric,			*fvFL_CurFric;              
	float       *fvFL_B4_Fric,			    *fvFL_TempRefFric,		*fvFL_CurDcVal;

    float    	*fvFL_PSeis_T0_F,          	*fvFL_AccumSlp;
	float       *fvBL_PSeis_T0_S,			*fvBL_PSeis_T0_N;
	float    	*fvFL_SlipRate_temp,      	*fvFL_SlipRake_temp;
	 
	int      	*ivFG_SegID_temp,			*ivFG_FltID_temp,		*ivFG_Flagged_temp;
    int    		*ivFG_V1_temp,              *ivFG_V2_temp,        	*ivFG_V3_temp;
	
    int    		*ivBG_V1_temp,              *ivBG_V2_temp,       	*ivBG_V3_temp;
    float    	*fvFG_CentE_temp,          	*fvFG_CentN_temp,    	*fvFG_CentZ_temp;
    float    	*fvBG_CentE_temp,           *fvBG_CentN_temp,    	*fvBG_CentZ_temp;

	gsl_matrix_int	    *imFGL_TTP,         *imFGL_NextP;
	gsl_matrix_int    	*imFGL_TTS,			*imFGL_NextS; 
    gsl_matrix_float    *fmFGL_SrcRcvH,    	*fmFGL_SrcRcvV,      	*fmFGL_SrcRcvN;    

	gsl_vector_float	*fvFL_StrsRateStk, 	*fvFL_StrsRateDip;
	gsl_vector_float	*fvFL_CurStrsH,		*fvFL_CurStrsV,			*fvFL_CurStrsN;
	gsl_vector_float    *fvFL_B4_StrsH,		*fvFL_B4_StrsV,		    *fvFL_B4_StrsN;    
	gsl_vector_float	*fvBL_CurStrsH,		*fvBL_CurStrsV,			*fvBL_CurStrsN;
};
//------------------------------------------------------------------
struct VTstruct
{	float   	*fvFG_VlX_temp,             *fvFG_VlY_temp;
	float  		*fvFG_PosE_temp,          	*fvFG_PosN_temp,      	*fvFG_PosZ_temp;
	float   	*fvBG_PosE_temp,           	*fvBG_PosN_temp,      	*fvBG_PosZ_temp;
};
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void LoadInput(struct MDstruct *MD, struct SGstruct *SG, struct TRstruct *TR, struct VTstruct *VT)
{	
    //------------------------------------------------------------------
	gsl_rng *fRandN; // this is pretty much straight from the GSL reference, the default RNG has good performance, so no need to change
    const gsl_rng_type *RandT;
	gsl_rng_env_setup(); 
	RandT   = gsl_rng_default;
    fRandN  = gsl_rng_alloc(RandT);
	unsigned long RSeed =(unsigned long)(MD->iRANK + MD->iSeedStart);   
    gsl_rng_set(fRandN, RSeed);
	//-----------------------------------------------------------------
	int		iTemp0,             	i,            			j,					iGlobPos;
	int    	iFPnumT,				iFVnumT,				iBPnumT,			iBVnumT;
    float   fTemp0,					fTemp1,					fTemp2,				fTemp3;
    char    ctempVals[512],         cAppend[512];
	char    cFileName1[512], 		cFileName2[512],		cFileName3[512],	cFileName4[512];
    FILE    *fp1,               	*fp2,           		*fp3; 
	//-------------------------------------
	//reading/loading the remaining information from txt and dat files.. this step is done by each rank individually; also, some values are read locally (every rank only reads the information relevant for it) while other values are read globally
    strcpy(cFileName1,MD->cInputName);          strcat(cFileName1,"_Summary_RoughStrength.txt"); 
    //------------------
	if      (MD->iRunNum == 0)			{	strcpy(cFileName2,MD->cInputName);      strcat(cFileName2,"_FLTtrig.dat");    																												 		}
	else								{	strcpy(cFileName2,MD->cInputName);      strcat(cFileName2,"_");               	sprintf(cAppend, "%d",MD->iRunNum); strcat(cFileName2,cAppend);     	strcat(cFileName2,"_Roughn.dat");			}
	//------------------
	strcpy(cFileName3,MD->cInputName);      strcat(cFileName3,"_");             	sprintf(cAppend, "%d",MD->iRunNum); 	strcat(cFileName3,cAppend);     	strcat(cFileName3,"_Strgth.dat"); 
	//------------------
	strcpy(cFileName4,MD->cInputName);      strcat(cFileName4,"_BNDtrig.dat");
    //-----------------------------------------------------------------
    if ((fp1 = fopen(cFileName1,"r")) == NULL)          { 	printf("Error -cant open *.flt file. LoadInputParameter function...\n");      exit(10);     }
    
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %e",&MD->fMedDense);                                          }
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %e",&MD->fAddNrmStrs);    MD->fAddNrmStrs *= -1.0;            }// keep at MPa but make compression negative ==> in MATLAB and the input files, normal stress is still positive for compression => is switched here to compression == negative convention
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %e",&MD->fShearMod);      MD->fShearMod   *=  1.0E+9;         }// convert from GPa to Pa (for actual computation of K-matrix and also for magnitude calculation, I need Pa and not MPa)
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %e",&MD->fPoisson);                                           }
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %d",&MD->iChgBtwEQs);                                         }
     
    MD->fLambda = (2.0*MD->fShearMod*MD->fPoisson)/(1.0-2.0*MD->fPoisson); //elastic parameters => for K-matrix calculation, remember that the code uses a linear elastic half-space
    fTemp0      = (MD->fMedDense > 0.0) ? MD->fMedDense : 2700.0; // if zero density is used for depth gradient (e.g., to have a depth-independent normal stress) then I still!! need to define a density for the wave propagation speed.... => is done here...
    MD->fVp     = sqrtf((MD->fLambda +2.0*MD->fShearMod)/fTemp0); // in m/s
    MD->fVs     = sqrtf(MD->fShearMod/fTemp0); // in m/s 
    //-------------------------------------------   
    for (i = 0; i < MD->iFSegNum; i++) //first, read the segment information (from the txt file)
    { 	for (j = 0; j < 21;    j++)                    
		{	if (fgets(ctempVals, 512, fp1) != NULL)     {   					} //these rows are skipped...	
        }
        if (fgets(ctempVals, 512, fp1) != NULL)         {   sscanf(ctempVals,"%*s %d", &SG->ivFricLawUSED[i]);          }                        
    	if (fgets(ctempVals, 512, fp1) != NULL)         {   sscanf(ctempVals,"%*s %e", &SG->fvRefStaFr[i]);             }
        if (fgets(ctempVals, 512, fp1) != NULL)         {   sscanf(ctempVals,"%*s %e", &SG->fvRefStaFr_vari[i]);        }
        if (fgets(ctempVals, 512, fp1) != NULL)         {   sscanf(ctempVals,"%*s %e", &SG->fvRefDynFr[i]);             }
        if (fgets(ctempVals, 512, fp1) != NULL)         {   sscanf(ctempVals,"%*s %e", &SG->fvRefDynFr_vari[i]);        }
        if (fgets(ctempVals, 512, fp1) != NULL)         {   sscanf(ctempVals,"%*s %e", &SG->fvRefCritDist[i]);          }
        if (fgets(ctempVals, 512, fp1) != NULL)         {   sscanf(ctempVals,"%*s %e", &SG->fvRefCritDist_vari[i]);     }      
        for (j = 0; j < 12;    j++)                     
		{ 	if (fgets(ctempVals, 512, fp1) != NULL)     {            	        } //these rows are skipped...
        }	
        if (SG->ivFricLawUSED[i] == 1)                   
		{	SG->fvRefCritDist[i]      = 0.0;  //set critical slip distance to zero when classic static/dynamic friction is used
			SG->fvRefCritDist_vari[i] = 0.0;    
    }	}
    fclose(fp1);
	//-----------------------------------------------------------------
	if ((fp2 = fopen(cFileName2,"rb"))     == NULL)     {   printf("Error -cant open %s LoadInputParameter function...\n",cFileName2);      exit(10);     }
    if ((fp3 = fopen(cFileName3,"rb"))     == NULL)     {   printf("Error -cant open %s LoadInputParameter function...\n",cFileName3);      exit(10);     }    
	//-----------------------------------------------------------------
    for (i = 1; i < MD->iFUsdGrd; i++)                  // all grids before the used one are tossed => read into "dummmy", means I read them normally but also overwrite them until the correct grid is processed
    {   if (fread(&iFPnumT,   sizeof(int),1,fp2) != 1)  {   exit(10);  	}    // this is currently read PatchNum  
        if (fread(&iFVnumT,   sizeof(int),1,fp2) != 1)  {   exit(10);  	}    // this is currently read VertexNum 
        fseek(fp2, (5L*sizeof(int)*(long)iFPnumT+5L*sizeof(float)*(long)iFVnumT +3L*sizeof(float)*(long)iFPnumT), SEEK_CUR); // skip over the values that I don't need to read
        fseek(fp3, (1L*sizeof(int)*1L           +2L*sizeof(float)*(long)iFPnumT +1L*sizeof(  int)*(long)iFPnumT), SEEK_CUR); // skip over the values that I don't need to read 
		if (fread(&iTemp0, sizeof(int), 1, fp2) != 1)	{	exit(10);	} 	       // this is the number of different faults (not sections) that have been defined
        fseek(fp2, (1L*sizeof(float)*(long)iTemp0), SEEK_CUR); // skip over the roughness values (RMS values) of the imported grids
    }
    if (fread(&iFPnumT, sizeof(int),1,fp2) != 1)  		{	exit(10);	}// this is currently read PatchNum 
    if (fread(&iFVnumT, sizeof(int),1,fp2) != 1) 		{	exit(10);	}// this is currently read VertexNum
 	//-----------------------------------------------------
    if ((MD->iFPNum != iFPnumT) || (MD->iFVNum != iFVnumT))            		{  	fprintf(stdout,"Patch or vertex number not matching: patchum1 %d  patchum2 %d       vertexnum1 %d    vertexnum2 %d\n", MD->iFPNum, iFPnumT, MD->iFVNum, iFVnumT);        exit(10);        			}
    //-----------------------------------------------------
	if (fread(TR->ivFG_V1_temp,    sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	} //reading the geometric information about the faults
    if (fread(TR->ivFG_V2_temp,    sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	}
    if (fread(TR->ivFG_V3_temp,    sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	}
	if (fread(TR->ivFG_SegID_temp, sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	}
    if (fread(TR->ivFG_FltID_temp, sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)	    {	exit(10);	}    
	//-------------------------------------
	if (fread(VT->fvFG_VlX_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)	    {	exit(10);	}
	if (fread(VT->fvFG_VlY_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)	    {	exit(10);	}
    //-------------------------------------
    if (fread(VT->fvFG_PosE_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)		{	exit(10);	}
    if (fread(VT->fvFG_PosN_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)		{	exit(10);	}
    if (fread(VT->fvFG_PosZ_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)		{	exit(10);	}
    //----------------------------------------
    if (fread(TR->fvFG_CentE_temp, sizeof(float),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	}
    if (fread(TR->fvFG_CentN_temp, sizeof(float),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	}   
    if (fread(TR->fvFG_CentZ_temp, sizeof(float),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	}
	//------------------------------------- 	
    fclose(fp2);
    //-----------------------------------------------------------------
    fseek(fp3, (1L*sizeof(int)),       SEEK_CUR); // contains the patch number again -> skipped
    
	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
	if (fread(TR->fvFL_RefStaFric, sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])	{	exit(10);	}
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
 
	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if(fread(TR->fvFL_RefDynFric, sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])	{	exit(10);	}
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum

	fseek(fp3, (1L*sizeof(  int)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
	if (fread(TR->ivFL_StabT, sizeof( int),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])			{	exit(10);	}
	fseek(fp3, (1L*sizeof(  int)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
	//-------------------------------------
    fclose(fp3);
    //-----------------------------------------------------------------
	if (MD->iBPNum > 0)
    {   if ((fp1 = fopen(cFileName4,"rb")) == NULL)   {   printf("Error -cant open  %s LoadInputParameter function...\n",cFileName4);     exit(10);     }

   		if (fread(&iBPnumT, sizeof(int),1,fp1) != 1)	{	exit(10);	}         // this is currently read PatchNum for boundary 
        if (fread(&iBVnumT, sizeof(int),1,fp1) != 1)	{	exit(10);	}         // this is currently read VertexNum for boundary
        //-------------------------------------
		if ((MD->iBPNum != iBPnumT) || (MD->iBVNum != iBVnumT))            		{  	fprintf(stdout,"Patch or vertex number not matching: patchum1 %d  patchum2 %d       vertexnum1 %d    vertexnum2 %d\n", MD->iBPNum, iBPnumT, MD->iBVNum, iBVnumT);        			exit(10);        			}
        //-------------------------------------
		if (fread(TR->ivBG_V1_temp, sizeof(int), MD->iBPNum,fp1) != MD->iBPNum)	{	exit(10);	}
    	if (fread(TR->ivBG_V2_temp, sizeof(int), MD->iBPNum,fp1) != MD->iBPNum)	{	exit(10);	} 
    	if (fread(TR->ivBG_V3_temp, sizeof(int), MD->iBPNum,fp1) != MD->iBPNum)	{	exit(10);	} 
        //-------------------------------------
        fseek(fp1, 1L*sizeof(  int)*(long)MD->iBPNum,   SEEK_CUR); // this is fault section ID => but don't need it => skip over it
        fseek(fp1, 2L*sizeof(float)*(long)MD->iBVNum,   SEEK_CUR); // this would be local coordintates (from gridding) of the vertices..; not needed    
        //-------------------------------------
	    if (fread(VT->fvBG_PosE_temp,  sizeof(float),MD->iBVNum,fp1) != MD->iBVNum)	{	exit(10);	}
    	if (fread(VT->fvBG_PosN_temp,  sizeof(float),MD->iBVNum,fp1) != MD->iBVNum)	{	exit(10);	}
    	if (fread(VT->fvBG_PosZ_temp,  sizeof(float),MD->iBVNum,fp1) != MD->iBVNum)	{	exit(10);	}
    	//----------------------------------------
    	if (fread(TR->fvBG_CentE_temp, sizeof(float),MD->iBPNum,fp1) != MD->iBPNum)	{	exit(10);	}
    	if (fread(TR->fvBG_CentN_temp, sizeof(float),MD->iBPNum,fp1) != MD->iBPNum)	{	exit(10);	}    
    	if (fread(TR->fvBG_CentZ_temp, sizeof(float),MD->iBPNum,fp1) != MD->iBPNum)	{	exit(10);	}
		//-------------------------------------  
        fclose(fp1);     
    }
	//-----------------------------------------------------------------
    for (i = 0; i < (MD->iFVNum); i++)        //converting all km to meters, also shifting the indizes => Matlab starting with "1" while C starting with "0" => hence the -1 used here....    
    {   VT->fvFG_PosE_temp[i]     *= 1000.0;            VT->fvFG_PosN_temp[i]     *= 1000.0;         VT->fvFG_PosZ_temp[i]     *= 1000.0; 
    }
	for (i = 0; i < (MD->iBVNum); i++)            
    {   VT->fvBG_PosE_temp[i]     *= 1000.0;            VT->fvBG_PosN_temp[i]     *= 1000.0;         VT->fvBG_PosZ_temp[i]     *= 1000.0; 
    }
	for (i = 0; i < (MD->iFPNum); i++)            
    {	TR->ivFG_V1_temp[i]       -= 1;             	TR->ivFG_V2_temp[i]       -= 1;              TR->ivFG_V3_temp[i]       -= 1;                  
        TR->fvFG_CentE_temp[i]    *= 1000.0;         	TR->fvFG_CentN_temp[i]    *= 1000.0;         TR->fvFG_CentZ_temp[i]    *= 1000.0; 
        TR->ivFG_SegID_temp[i]    -= 1; 
	}
	for (i = 0; i < (MD->iBPNum); i++)            
    {	TR->ivBG_V1_temp[i]        -= 1;                TR->ivBG_V2_temp[i]       -= 1;              TR->ivBG_V3_temp[i]       -= 1;                  
        TR->fvBG_CentE_temp[i]  *= 1000.0;         		TR->fvBG_CentN_temp[i]    *= 1000.0;         TR->fvBG_CentZ_temp[i]    *= 1000.0; 
	}
	//-----------------------------------------------------------------

	
	
    for (i = 0; i < MD->ivF_OFFSET[MD->iRANK]; i++) //here the current static/dynamic friction coefficients are determines, also the current Dc value; further, the stressing rates are assigned
    {   
		iGlobPos                       = i + MD->ivF_START[MD->iRANK];
		fTemp0                         = cosf(SG->fvSlipRake[TR->ivFG_SegID_temp[iGlobPos]]) *SG->fvStrsRate[TR->ivFG_SegID_temp[iGlobPos]]; //get stressing rate in strike direction
		fTemp1                         = sinf(SG->fvSlipRake[TR->ivFG_SegID_temp[iGlobPos]]) *SG->fvStrsRate[TR->ivFG_SegID_temp[iGlobPos]]; //the stressing rate in dip direction
        gsl_vector_float_set(TR->fvFL_StrsRateStk, i, fTemp0); 
        gsl_vector_float_set(TR->fvFL_StrsRateDip, i, fTemp1);
        
		TR->fvFL_SlipRate_temp[i]  = SG->fvSlipRate[   TR->ivFG_SegID_temp[iGlobPos]]; //if a slip boundary condition was used, then it is put in here
        TR->fvFL_SlipRake_temp[i]  = SG->fvSlipRake[   TR->ivFG_SegID_temp[iGlobPos]]; //the strike and dip component respectively
        TR->fvFL_RefDcVal[i]       = SG->fvRefCritDist[TR->ivFG_SegID_temp[iGlobPos]]*MD->fLegLgth;

		TR->fvFL_RefDcVal_vari[i]  = SG->fvRefCritDist_vari[TR->ivFG_SegID_temp[iGlobPos]]/100.0; //the permitted random variation of Dc, and friction coefficients
		TR->fvFL_RefStaFric_vari[i]= SG->fvRefStaFr_vari[   TR->ivFG_SegID_temp[iGlobPos]]/100.0; //they are used here to determine the current valuesof Dc and friction coefficents;
		TR->fvFL_RefDynFric_vari[i]= SG->fvRefDynFr_vari[   TR->ivFG_SegID_temp[iGlobPos]]/100.0; //they are further used when those values are allowed to change between earthquakes
		TR->ivFL_FricLaw[i]        = SG->ivFricLawUSED[     TR->ivFG_SegID_temp[iGlobPos]];
		//------------------------
		if (TR->ivFL_FricLaw[i] == 1)		{	TR->fvFL_RefDcVal[i] = 0.0;		TR->fvFL_RefDcVal_vari[i]  = 0.0;			}
		//------------------------
        fTemp2                     = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
		TR->fvFL_CurDcVal[i]       = TR->fvFL_RefDcVal[i]  *(1.0 + TR->fvFL_RefDcVal_vari[i]*fTemp2);
        //------------------------
		fTemp2                     = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
        TR->fvFL_StaFric[i]        = TR->fvFL_RefStaFric[i]*(1.0 + TR->fvFL_RefStaFric_vari[i]*fTemp2);     
        TR->fvFL_CurFric[i]        = TR->fvFL_StaFric[i];
		//------------------------
		fTemp1                     = (TR->fvFL_RefStaFric[i] - TR->fvFL_RefDynFric[i])/TR->fvFL_RefStaFric[i]; //reference friction change as fraction of static coefficient
		fTemp2                     = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
        
		fTemp3                     = fTemp1*(1.0 + TR->fvFL_RefDynFric_vari[i]*fTemp2); //calculation of dyn friction coefficient looks a bit more complicated, this is because this value is defined relative to the static coefficent i.e., as a fraction change of that coefficient
        TR->fvFL_DynFric[i]        = TR->fvFL_StaFric[i]*(1.0 - fTemp3);
    }
    //-----------------------------------------------------------------

	if (MD->iUseVpVs == 1)           {   MD->fVpVsRatio = MD->fVp/MD->fVs;                 }			else                            {   MD->fVpVsRatio = 1.0;                			}              
	if (MD->iUseProp == 1)           {   MD->fDeltT     = 1.0*MD->fLegLgth/MD->fVp;        }			else                            {   MD->fDeltT     = FLT_MAX;      	     			}

    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void DefineMoreParas(struct MDstruct *MD, struct SGstruct *SG, struct TRstruct *TR, struct VTstruct *VT, int iPlot2Screen)
{	int 	i,				j,				iGlobPos;
	int    	iTemp0,			iTemp1;
	float	fTemp0,			fTemp1,			fTemp2;
	float	fP1[3],			fP2[3],			fP3[3];
	float	fP1P2[3],		fP1P3[3],		fSrcRcvVect[3];
	float	fvNrm[3],		fvStk[3],		fvDip[3];	
	float   feY[3],			feZ[3];
	//--------------------------------------------------------------------
	for (i = 0; i < MD->ivF_OFFSET[MD->iRANK];  i++)
    {   iGlobPos = i + MD->ivF_START[MD->iRANK];	
		fP1[0]   = VT->fvFG_PosE_temp[TR->ivFG_V1_temp[iGlobPos]];			fP1[1]   = VT->fvFG_PosN_temp[TR->ivFG_V1_temp[iGlobPos]];			fP1[2]   = VT->fvFG_PosZ_temp[TR->ivFG_V1_temp[iGlobPos]];
		fP2[0]   = VT->fvFG_PosE_temp[TR->ivFG_V2_temp[iGlobPos]];			fP2[1]   = VT->fvFG_PosN_temp[TR->ivFG_V2_temp[iGlobPos]];			fP2[2]   = VT->fvFG_PosZ_temp[TR->ivFG_V2_temp[iGlobPos]];
		fP3[0]   = VT->fvFG_PosE_temp[TR->ivFG_V3_temp[iGlobPos]];			fP3[1]   = VT->fvFG_PosN_temp[TR->ivFG_V3_temp[iGlobPos]];			fP3[2]   = VT->fvFG_PosZ_temp[TR->ivFG_V3_temp[iGlobPos]];
		fP1P2[0] = fP2[0] - fP1[0];                fP1P2[1] = fP2[1] - fP1[1];                fP1P2[2] = fP2[2] - fP1[2];
        fP1P3[0] = fP3[0] - fP1[0];                fP1P3[1] = fP3[1] - fP1[1];                fP1P3[2] = fP3[2] - fP1[2];
        fvNrm[0] = fP1P2[1]*fP1P3[2] - fP1P2[2]*fP1P3[1];
        fvNrm[1] = fP1P2[2]*fP1P3[0] - fP1P2[0]*fP1P3[2];
        fvNrm[2] = fP1P2[0]*fP1P3[1] - fP1P2[1]*fP1P3[0];
		//--------------------------------------------------------------------
		TR->fvFL_Area[i]       = 0.5*sqrtf( fvNrm[0]*fvNrm[0] +fvNrm[1]*fvNrm[1] +fvNrm[2]*fvNrm[2]); //is needed for moment/magnitue calculation
		TR->fvFL_RefNrmStrs[i] = -1.0e-6*(MD->fMedDense *MD->fg *fabs(TR->fvFG_CentZ_temp[iGlobPos])) + MD->fAddNrmStrs; // all in MPa now!!; again -1 to convert to compressive=negative 
		//--------------------------------------------------------------------
		//determine local coordinate system for currently selected (local) fault patch
		feY[0] = 0.0;                       feY[1] = 1.0;                               feY[2] = 0.0;
    	feZ[0] = 0.0;                       feZ[1] = 0.0;                               feZ[2] = 1.0;        
    	//-----------------------------------------  
    	fTemp0      = sqrtf(fvNrm[0]*fvNrm[0] +fvNrm[1]*fvNrm[1] +fvNrm[2]*fvNrm[2]);
    	fvNrm[0]    = fvNrm[0]/fTemp0; 	fvNrm[1] = fvNrm[1]/fTemp0;      		fvNrm[2] = fvNrm[2]/fTemp0;
    	if (fvNrm[2] < 0.0)     
		{   fvNrm[0] = -fvNrm[0];         	fvNrm[1] = -fvNrm[1];            			fvNrm[2] = -fvNrm[2];     }
    	//----------------------------------------- 
    	fvStk[0]    = feZ[1]*fvNrm[2] - feZ[2]*fvNrm[1];
    	fvStk[1]    = feZ[2]*fvNrm[0] - feZ[0]*fvNrm[2];
    	fvStk[2]    = feZ[0]*fvNrm[1] - feZ[1]*fvNrm[0];
    	// For horizontal elements ("Vnorm(3)" adjusts for Northward or Southward direction) 
    	fTemp0  = sqrtf(fvStk[0]*fvStk[0] +fvStk[1]*fvStk[1] +fvStk[2]*fvStk[2]);
    	if (fabs(fTemp0) < FLT_EPSILON)
		{   fvStk[0]= 0.0;                 fvStk[1] = feY[1]*fvNrm[2];        fvStk[2] = 0.0;             }
    	fTemp0       = sqrtf(fvStk[0]*fvStk[0] +fvStk[1]*fvStk[1] +fvStk[2]*fvStk[2]);
    	fvStk[0]    = fvStk[0]/fTemp0;      fvStk[1] = fvStk[1]/fTemp0;         fvStk[2] = fvStk[2]/fTemp0;
    	//----------------------------------------- 
    	fvDip[0]    = fvNrm[1]*fvStk[2] - fvNrm[2]*fvStk[1];
    	fvDip[1]    = fvNrm[2]*fvStk[0] - fvNrm[0]*fvStk[2];
    	fvDip[2]    = fvNrm[0]*fvStk[1] - fvNrm[1]*fvStk[0];
    	fTemp0      = sqrtf(fvDip[0]*fvDip[0] +fvDip[1]*fvDip[1] +fvDip[2]*fvDip[2]);
    	fvDip[0]    = fvDip[0]/fTemp0;      fvDip[1] = fvDip[1]/fTemp0;         fvDip[2] = fvDip[2]/fTemp0;
		//-------------------------------------------------------------------- 
		for (j = 0; j < MD->iFPNum; j++)
        {   
			fTemp0   = sqrtf( ( (TR->fvFG_CentE_temp[j]-TR->fvFG_CentE_temp[iGlobPos])*(TR->fvFG_CentE_temp[j]-TR->fvFG_CentE_temp[iGlobPos]) ) + ( (TR->fvFG_CentN_temp[j]-TR->fvFG_CentN_temp[iGlobPos])*(TR->fvFG_CentN_temp[j]-TR->fvFG_CentN_temp[iGlobPos]) ) + ( (TR->fvFG_CentZ_temp[j]-TR->fvFG_CentZ_temp[iGlobPos])*(TR->fvFG_CentZ_temp[j]-TR->fvFG_CentZ_temp[iGlobPos])) )/MD->fVp; // the distance between both 
            fTemp1   = fTemp0*MD->fVpVsRatio; //these are the travel times from current "source" to current "receiver"
 			//----------------------------------------- 
			iTemp0   = (int)(fTemp0/MD->fDeltT); //means it's rounding downm, cutting off whatever floating point contribution the travel times had
			iTemp1   = (int)(fTemp1/MD->fDeltT); 
            gsl_matrix_int_set(TR->imFGL_TTP, i, j, iTemp0);
			gsl_matrix_int_set(TR->imFGL_TTS, i, j, iTemp1);
            //----------------------------------------- 
            MD->iGlobTTmax   = (MD->iGlobTTmax > iTemp1) ? MD->iGlobTTmax : iTemp1;     
			//----------------------------------------- 
            fSrcRcvVect[0]   = TR->fvFG_CentE_temp[j] - TR->fvFG_CentE_temp[iGlobPos]; // this is vector from current patch (source, local) to receiver (global); this one points therefore from source to the receiver
            fSrcRcvVect[1]   = TR->fvFG_CentN_temp[j] - TR->fvFG_CentN_temp[iGlobPos]; //will be needed for vector projection => when using rupture propagation => which component of slip vector is in direct line towards receiver patch (mode II) and which component is perpendicular to that (mode III)
            fSrcRcvVect[2]   = TR->fvFG_CentZ_temp[j] - TR->fvFG_CentZ_temp[iGlobPos]; //this is therefore for rupture propagation and the transient signals that arrise when using different velocities for mode II and mode III
            
            if (i == iGlobPos)
            {  	gsl_matrix_float_set(TR->fmFGL_SrcRcvN, i, j, 0.0);
				gsl_matrix_float_set(TR->fmFGL_SrcRcvH, i, j, 0.0);
				gsl_matrix_float_set(TR->fmFGL_SrcRcvV, i, j, 0.0);
			}
            else
            {   fTemp0          = sqrtf(fSrcRcvVect[0]*fSrcRcvVect[0] +fSrcRcvVect[1]*fSrcRcvVect[1] +fSrcRcvVect[2]*fSrcRcvVect[2]);
                fSrcRcvVect[0] /= fTemp0;                fSrcRcvVect[1] /= fTemp0;                fSrcRcvVect[2] /= fTemp0; 
				fTemp0	        = fvNrm[0]*fSrcRcvVect[0] + fvNrm[1]*fSrcRcvVect[1] + fvNrm[2]*fSrcRcvVect[2]; //this rotates the vector from current source to receiver= fvNrm[0]*fSrcRcvVect[0] + fvNrm[1]*fSrcRcvVect[1] + fvNrm[2]*fSrcRcvVect[2];
				fTemp1          = fvStk[0]*fSrcRcvVect[0] + fvStk[1]*fSrcRcvVect[1] + fvStk[2]*fSrcRcvVect[2]; //into local coordinate system of the source patch(that is the idea...)
				fTemp2		    = fvDip[0]*fSrcRcvVect[0] + fvDip[1]*fSrcRcvVect[1] + fvDip[2]*fSrcRcvVect[2];         

                gsl_matrix_float_set(TR->fmFGL_SrcRcvN, i, j, fTemp0);
				gsl_matrix_float_set(TR->fmFGL_SrcRcvH, i, j, fTemp1);
				gsl_matrix_float_set(TR->fmFGL_SrcRcvV, i, j, fTemp2);        
    }   }	}

	if ((iPlot2Screen == 1) && (MD->iRANK == 0))
    {	fprintf(stdout,"Medium density:              %5.2f\n",MD->fMedDense);
		fprintf(stdout,"Shear modulus:               %5.2e\n",MD->fShearMod);
 		fprintf(stdout,"Poisson ratio:               %5.2f\n",MD->fPoisson);
 		fprintf(stdout,"Lambda:                      %5.2e\n",MD->fLambda);
		fprintf(stdout,"P-wave velocity:             %5.2f\n",MD->fVp);
		fprintf(stdout,"S-wave velocity:             %5.2f\n",MD->fVs);
		fprintf(stdout,"Added normal stress (MPa):   %5.2f\n",MD->fAddNrmStrs);
		fprintf(stdout,"Change friction between EQs: %d\n",MD->iChgBtwEQs);
		fprintf(stdout,"\nValues per Segment:\n");
		for (i = 0; i < MD->iFSegNum; i++)
		{	fprintf(stdout,"Segment# %d   SegmentSlipRate (m/yr):       %5.5f\n",i,   SG->fvSlipRate[i]);
			fprintf(stdout,"Segment# %d   SegmentSlipRake:              %5.2f\n",i,   SG->fvSlipRake[i]);
			fprintf(stdout,"Segment# %d   SegmentStressRate (MPa):      %5.2f\n",i, SG->fvStrsRate[i]);
			fprintf(stdout,"Segment# %d   Friction law used:            %d\n", i, SG->ivFricLawUSED[i]);
			fprintf(stdout,"Segment# %d   Reference static friction:    %4.4f\n", i, SG->fvRefStaFr[i]);
			fprintf(stdout,"Segment# %d   Static friction variation:    %4.4f\n", i, SG->fvRefStaFr_vari[i]);
			fprintf(stdout,"Segment# %d   Reference dynamic friction:   %4.4f\n", i, SG->fvRefDynFr[i]);
			fprintf(stdout,"Segment# %d   Dynamic friction variation:   %4.4f\n", i, SG->fvRefDynFr_vari[i]);
			fprintf(stdout,"Segment# %d   Critical slip distance (m):   %4.6f\n", i, SG->fvRefCritDist[i]*MD->fLegLgth);
			fprintf(stdout,"Segment# %d   Critical slip variation:      %4.4f\n\n", i, SG->fvRefCritDist_vari[i]);
	}	}

	return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
