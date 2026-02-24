/*SWT Interpretation Algorithm: Scale invariant filtering of bathymetry data*/
/*Last Altered:21/08/2008  -  John Hillier*/

/*COMMAND LINE*/
/*SWT.out
[(1)input.file distance,depth,longitude,latitude] 
[(2)file raw wavelet coefficients] 
[(3)file blocked at each scale] 
[(4)filtered output] 
[(5)combined selected across the scales]
[(6)Combined & blocked] 
*/


/*LIBRARIES*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*DEFINITIONS*/
#define MAXCOLUMNS 50000  /*Max number of locations in file*/
#define MAXLINE 1000  /*Max number characters in a line read from file*/


/*THE FUNCTION DEFINITIONS*/
int filter(float distance[], float depth[], float longitude[], float latitude[], float filtered[], char *argv[], int argc, int numberlines);
int formatarray(float another[]);
int formatarrayI(int another[]);
int getlineJH(char line[], FILE *fdata);
int Interpolate(float WavComb[], int numberlines, float WavFilt[], float depth[], float distance[], float WavCombScale[], int WavCombLeft[], int WavCombRight[]);
float quick_select(float arr[], int n);
int PostProcess(float WavComb[], int numberlines, float WavFilt[], float depth[], float distance[], float WavCombScale[], int WavCombLeft[], int WavCombRight[]);
int ReadData (char *argv[], float distance[], float depth[], float longitude[], float latitude[]);
int tempInterp(float distance[], float depth[], float tempFiltered[],int start,int end);
int wavCompute(float distance[], float depth[], float WavCoeff[], float scale, float wavWeighting[], int numberlines);
int wavInterpii(float WavCoeff[], int numberlines, float CoeffThreshold, float distance[], float scale);
int wavInterpComb(float WavComb[], int numberlines, float WavFilt[], float depth[], float distance[], float WavCombScale[], int WavCombLeft[], int WavCombRight[]);
int wavInterpCombii(char *argv[], float WavComb[], int numberlines, float distance[], float WavCombScale[], float ScaleInteract);
int wavelet(float distance[], float depth[], char *argv[], int numberlines);



/*MAIN*/
main(int argc, char *argv[])
{
	float distance[MAXCOLUMNS];
	float depth[MAXCOLUMNS];
	float longitude[MAXCOLUMNS];
	float latitude[MAXCOLUMNS];
	int numberlines;

/*Reading the data into arrays*/
	numberlines = (ReadData(argv, distance, depth, longitude, latitude));
	/*printf("numberlines %d\n", numberlines);*/
	printf("Read In The Data\n");

/*Wavelet - Filtering the data*/
	printf("Commencing Wavelet Filter\n");
	wavelet(distance, depth, argv, numberlines);
	printf("Done Wavelet Fiter\n");
}


/*FUNCTIONS - ALPHABETICAL*/

/*FORMATARRAY*/
/*formatarray: formats an array ready to be read into*/
int formatarray(float another[])

{
        int i;
	/*'<' is there as arrays start at [0]*/
	for (i = 0; i < MAXCOLUMNS; i++)
                another[i] = 0;
	return 0;
}

/*FORMATARRAYI*/
/*formatarrayI: formats an array ready to be read into*/
int formatarrayI(int another[])

{
        int i;
	/*'<' is there as arrays start at [0]*/
	for (i = 0; i < MAXCOLUMNS; i++)
                another[i] = 0;
	return 0;
}


/*GETLINEJH*/
/*getlineJH: reads a line into LINE[]*/
int getlineJH(char s[], FILE *fdata)

{
	int c, i;

/*printf("in getlineJH\n")*/;

	for (i=0; (c=getc(fdata))!=EOF && c!='\n'; ++i)
		s[i] = c;
	if (c == '\n')
	{
		s[i] = c;
		++i;
	}
	s[i] = '\0';

/*Is return the last thing that the loop wants to see?*/
/*printf("got line")*/
	return i;

}



/*INTERPOLATE*/
/*Interpolate: Takes the limits of the objects e.g. from WavCombLeft[] 
and creates a filtered value.*/
int Interpolate(float WavComb[], int numberlines, float WavFilt[], float depth[], float distance[], float WavCombScale[], int WavCombLeft[], int WavCombRight[])

{
  int e,f,g;		/*loop variable*/
  float gradient;		/*used in the linear interpolation of the space under objects*/
  

  /*initialise stuff*/
  gradient = 0;

  printf("\nDoing interpolation.......\n");

  /*Start of by re-setting WavFilt[] to depth[]*/
  /*avoids any carry through errors*/
  for (e = 0; e < numberlines; e++)
    {
      WavFilt[e] = depth[e];
      /*printf("WavFilt[%d] %f\n",e,WavFilt[e]);*/
    }

  /*In the objects, interpolate between the limits*/
  for (f = 0; f < numberlines; f++)
    {
      if (WavComb[f] > 0)
	{
	  /*Gradient for the line under the object*/
	  /*printf("Extrapolate under object centred on [%d]\n",f);*/
	  /*printf("Limits are (%f,%f) (%f,%f)\n",depth[WavCombLeft[f]],distance[WavCombLeft[f]],depth[WavCombRight[f]],distance[WavCombRight[f]]);*/
	  gradient = (depth[WavCombLeft[f]]-depth[WavCombRight[f]])/(distance[WavCombLeft[f]]-distance[WavCombRight[f]]);
	  /*printf("gradient is %f\n",gradient);*/
	  /*And do the interpolation*/
	  for (g = WavCombLeft[f]; g <= WavCombRight[f]; g++)
	    {
	      WavFilt[g] = depth[WavCombLeft[f]] + ((distance[g]-distance[WavCombLeft[f]])*gradient);
	    }
	}
    }

	return 0;
}





/*POSTPROCESS*/
/*PostProcess: Second generation of post processing routine*/
int PostProcess(float WavComb[], int numberlines, float WavFilt[], float depth[], float distance[], float WavCombScale[], int WavCombLeft[], int WavCombRight[])

{
  int a,b,c,d,e,f,h,i,j,k,l,m;		/*loop variables*/
  float gradient;		/*used in the linear interpolation of the space under objects*/
  float area, outline;          /*properties of the feature*/
  float extradist, extraarea;
  float ratio, ratiomax;        /*ratio: outline length / Area */
  float depthInterpL, depthInterpR;    /*interpolated depths used in calculations*/
  float maxratio, maxratioOld;         /*'success parameters' driving the routine*/
  int moveL, moveR;           /*number of data points to the left to try to shift the left hand edge*/
  int loopleft, loopright;        /*loop logical variable*/
  int mainloop, changes;          /*logic variable for the main loop and number of changes*/
  float objToDo[MAXCOLUMNS];      /*Array 1 = to do, 0 = done*/
  int loopobjects;                /*logic in the loop controlling the object that is being done*/
  float MaxCoeff;
  int MaxCoeffLocn;             /*location of max remaining coefficient*/
  int ObjectsCount;               /*How many object are there still to do*/
  int criteriaFulfilled;
  float widthFract, MoveMax, Move;
  int summitlocation[MAXCOLUMNS];   /*Array containing the summit locations of the features
				      as this is not necessarily where the max coeff is located*/
  float GradCentre, GradCleft, GradCright;    /*variables for computation of the slope-based restriction*/
  float height, heightmax, depthInterp; 
  float mingradient, RHSgradient, LHSgradient;   /*gradients of the extremities and their minimum allowed value*/
  float summitToEdge;
  int InsidePosnLHS, InsidePosnRHS;
  float outsideFract, TooFlatFract;           /*user-set geometric parameters for the slope based restriction*/
  float HeightDiffFract, ScaleMulLimit;       /*more user set parameters for real data*/
  int TooBig, TooDeep;                                 /*variable associated with ScaleMulLimit*/
  float difftoosmall;
  float diffRside, diffLside, diffratio;
  float LargestMinGradient;                    /*caps the 'mingradient' that represents the floodplain*/
  float LowestMaxMove;                         /*minimum on this movement restriction because of data density*/

  printf("\nIn Post-Processing .... \n");

  widthFract = 0.25;        /*[0.25]amount of width of object upto which it is allowed to search for
			      increases in area/outline*/
  outsideFract = 0.25;      /*[0.25]amount of object considered as outside for the slope decrease test*/
  TooFlatFract = 1/7.0;     /*[1/7]amount of flattening that is considered to be a flat area outside the object*/
                            /*To turn this off use -1 instead of 1/8.0*/
  HeightDiffFract = 1.25;   /*[1.25]maximum fraction that the height differences at the 2 sides can vary by*/
  difftoosmall = 0.1;       /*[200]height differences below this are too small to worry aobut on ship-track data*/

  ScaleMulLimit = 400.0;      /*[400 - deactivated]max times the originally estimated feature width that anything can get to*/
  LargestMinGradient = 7.0;  /*[7]Some things are really really spikey in the middle,
			       but I don't want this to prevent them spreading at the base*/
  LowestMaxMove = 5.0;      /*[5]as data are irregular, have got to be able to have a max move of a few data points*/

  /*initialise stuff*/
  gradient = 0;
  formatarray(objToDo);
  formatarrayI(summitlocation);

  /*Set up the array of objects to do*/
  for (b = 0; b < numberlines; b++)
    {
      if (WavComb[b] == 0)
	{
	  objToDo[b] = 0;
	}
      if (WavComb[b] > 0)
	{
	  objToDo[b] = 1;
	}
    }

  /*Adjust the limits of the objects so that it is correctly delimited*/
  loopobjects = 1;
  /*While there are some objects still to do,*/
  /*going from the largest first*/
  while (loopobjects == 1)
    {
      /*Find the largest remaining object*/
      MaxCoeff = 0;
      MaxCoeffLocn = 0;
      for (h = 0; h < numberlines; h++)
	{
	  if (objToDo[h] == 1 && WavComb[h] > MaxCoeff)
	    {
	      MaxCoeffLocn = h;
	      MaxCoeff = WavComb[h];
	    }
	}
      printf("MaxCoeffLocn %d MaxCoeff %f\n", MaxCoeffLocn, MaxCoeff);
      
      /*so the object I'm going to do is ......*/
      /*f is used for legacy reasons*/
      f = MaxCoeffLocn;

      /*So for this object  ........... */
      
      /*Work out what the initial parameters are*/
      /*Print some details*/
      printf("Limits are (%f,%f) (%f,%f)\n",depth[WavCombLeft[f]],distance[WavCombLeft[f]],depth[WavCombRight[f]],distance[WavCombRight[f]]);
      
      /*initialize some properties*/
      outline = 0;
      area = 0;
      
      /*Work out the properties*/
      /*gradient of line under feature*/
      gradient = (depth[WavCombLeft[f]]-depth[WavCombRight[f]])/(distance[WavCombLeft[f]]-distance[WavCombRight[f]]);
      printf("gradient is %f\n",gradient);
      /*The length of the outline*/
      for (i = WavCombLeft[f];i < WavCombRight[f]; i++)
	{
	  /*over the top - use pythag*/
	  extradist = sqrt(((depth[i]-depth[i+1])*(depth[i]-depth[i+1]))+((distance[i]-distance[i+1])*(distance[i]-distance[i+1])));
	  outline = outline + extradist;
	  /*printf("[%d] outline %f extradist %f\n",i, outline, extradist);*/
	}
      /*and the line underneath - use pythag*/
      outline = outline + sqrt(((depth[WavCombLeft[f]]-depth[WavCombRight[f]])*(depth[WavCombLeft[f]]-depth[WavCombRight[f]]))+((distance[WavCombLeft[f]]-distance[WavCombRight[f]])*(distance[WavCombLeft[f]]-distance[WavCombRight[f]])));
      printf("[%d] total outline (inc. base) %f \n",f, outline);
      /*The area of the feature*/
      /*There is no exception where the measured line crosses the extrapolated*/
      for (j = WavCombLeft[f];j < WavCombRight[f]; j++ )
	{
	  depthInterpL = depth[WavCombLeft[f]] + ((distance[j]-distance[WavCombLeft[f]])*gradient);
	  depthInterpR = depth[WavCombLeft[f]] + ((distance[j+1]-distance[WavCombLeft[f]])*gradient);
	  extraarea = 0.5*((depth[j]-depthInterpL) + (depth [j+1] - depthInterpR))*(distance[j+1]-distance[j]);
	  area = area + extraarea;
	  /*printf("[%d] area %f extraarea %f\n",j, area, extraarea);*/
	}
      printf("[%d] total area %f \n",f, area);
      ratio = area/outline;
      printf("area/outline %f\n",ratio);
      /*The slope of this intial part, which will remain at the centre of the object when edges expand*/
      /*Calculated from summit to the edges*/
      heightmax = 0;
      for (k = WavCombLeft[f]; k <= WavCombRight[f]; k++)
	{
	  depthInterp = depth[WavCombLeft[f]] + ((distance[k]-distance[WavCombLeft[f]])*gradient);
	  height = depth[k] - depthInterp;
	  if (height > heightmax)
	    {
	      heightmax = height;
	      summitlocation[f] = k;
	    }
	}
      /*fairly crude calculation of the gradients*/
      GradCright =  ((depth[summitlocation[f]])-(depth[WavCombRight[f]]))/((distance[WavCombRight[f]])-(distance[summitlocation[f]]));
      GradCleft= ((depth[summitlocation[f]])-(depth[WavCombLeft[f]]))/((distance[summitlocation[f]])-(distance[WavCombLeft[f]]));
      GradCentre = (GradCright+GradCleft)/2.0;
      printf("Summit Location: %f and height %f \n", distance[summitlocation[f]], depth[summitlocation[f]]);
      printf("Right edge location: %f and height %f \n",distance[WavCombRight[f]], depth[WavCombRight[f]]);
      printf("Left edge location: %f and height %f \n",distance[WavCombLeft[f]], depth[WavCombLeft[f]]);
      printf("GradCleft: %f GradCright: %f ",GradCleft, GradCright);
      printf("Average gradient of core of feature is %f\n", GradCentre);
      /*So the minimum gradient allowed at the edges is .....*/
      mingradient = GradCentre*TooFlatFract;
      printf("So minimum gradient allowed is %f\n", mingradient);
      /*If this minimum is really large, lower it*/
      if (mingradient > LargestMinGradient)
	{
	  mingradient = LargestMinGradient;
	  printf("Really steep middle part, so lowering abysaal plain cut-off gradient to %f\n", LargestMinGradient);
	}

      /*Main loop controlling the expanding of the object*/
      maxratio = ratio;
      mainloop = 1;
      while (mainloop == 1)
	{
	  /*reset the number of changes that have happened*/
	  changes = 0;
	  printf("\n\nNEW loop changes = %d\n", changes);
	  
	  
	  
	  /*Try Moving the left hand limit*/
	  /*Always move on first success*/
	  printf("Trying to move the left hand limit\n");
	  loopleft = 1;
	  moveL = 1;
	  while (loopleft == 1)
	    { 
	      /*initialize some properties*/
	      outline = 0;
	      area = 0;
	      
	      /*Print some details*/
	      printf("Limits are (%f,%f) (%f,%f)\n",depth[(WavCombLeft[f]-moveL)],distance[(WavCombLeft[f]-moveL)],depth[WavCombRight[f]],distance[WavCombRight[f]]);
	      
	      /*Work out the properties*/
	      /*gradient of line under feature*/
	      gradient = (depth[(WavCombLeft[f]-moveL)]-depth[WavCombRight[f]])/(distance[(WavCombLeft[f]-moveL)]-distance[WavCombRight[f]]);
	      printf("gradient is %f\n",gradient);
	      /*The length of the outline*/
	      for (i = (WavCombLeft[f]-moveL);i < WavCombRight[f]; i++)
		{
		  /*over the top - use pythag*/
		  extradist = sqrt(((depth[i]-depth[i+1])*(depth[i]-depth[i+1]))+((distance[i]-distance[i+1])*(distance[i]-distance[i+1])));
		  outline = outline + extradist;
		  /*printf("[%d] outline %f extradist %f\n",i, outline, extradist);*/
		}
	      /*and the line underneath - use pythag*/
	      outline = outline + sqrt(((depth[(WavCombLeft[f]-moveL)]-depth[WavCombRight[f]])*(depth[(WavCombLeft[f]-moveL)]-depth[WavCombRight[f]]))+((distance[(WavCombLeft[f]-moveL)]-distance[WavCombRight[f]])*(distance[(WavCombLeft[f]-moveL)]-distance[WavCombRight[f]])));
	      printf("[%d] total outline (inc. base) %f \n",f, outline);
	      /*The area of the feature*/
	      /*There is no exception where the measured line crosses the extrapolated*/
	      for (j = (WavCombLeft[f]-moveL);j < WavCombRight[f]; j++ )
		{
		  depthInterpL = depth[(WavCombLeft[f]-moveL)] + ((distance[j]-distance[(WavCombLeft[f]-moveL)])*gradient);
		  depthInterpR = depth[(WavCombLeft[f]-moveL)] + ((distance[j+1]-distance[(WavCombLeft[f]-moveL)])*gradient);
		  extraarea = 0.5*((depth[j]-depthInterpL) + (depth [j+1] - depthInterpR))*(distance[j+1]-distance[j]);
		  area = area + extraarea;
		  /*printf("[%d] area %f extraarea %f\n",j, area, extraarea);*/
		}
	      printf("[%d] total area %f \n",f, area);
	      ratio = area/outline;
	      printf("area/outline %f\n",ratio);
	      
	      /*print current state of movement*/
	      printf("moveL %d\n", moveL);
	      
	      /*Check that the gradient of the outside portion of the proposed move is not too flat*/
	      /*Make the proportion with respect to the summit*/
	      /*avoids going too far onto the abyssal plain*/
	      /*Left hand side*/
	      /*work out how far back towards the centre to take the outside part*/
	      LHSgradient = 0;
	      summitToEdge = (distance[summitlocation[f]]-distance[(WavCombLeft[f] - moveL)]);
	      for (l = (WavCombLeft[f] - moveL);(distance[l]-distance[(WavCombLeft[f] - moveL)])/summitToEdge < outsideFract;l++)
		{
		  InsidePosnLHS = l;
		}
	      /*make sure that InsidePosnLHS is moved inside by at least 1 data point*/
	      if (InsidePosnLHS == (WavCombLeft[f] - moveL))
		{
		  InsidePosnLHS++;
		}
	      LHSgradient = (depth[InsidePosnLHS]-depth[(WavCombLeft[f] - moveL)])/(distance[InsidePosnLHS]-distance[(WavCombLeft[f] - moveL)]);
	      printf("Summit location: %d\n", summitlocation[f]);
	      printf("WavCombLeft[%d] %d InsidePosnLHS: %d LHSgradient %f \n",f,WavCombLeft[f]-moveL,InsidePosnLHS, LHSgradient);
	      
	      /*Check that the object hasn't expanded too far*/
	      TooBig = 0;
	      printf("proposed left limit: %f\n", distance[WavCombLeft[f] - moveL]);
	      printf("summit %f km \n", distance[summitlocation[f]]);
	      printf("size of extension allowed %f km\n", ((WavCombScale[f]/4.0)*ScaleMulLimit) );
	      if (distance[WavCombLeft[f] - moveL] < distance[summitlocation[f]] - ((WavCombScale[f]/4.0)*ScaleMulLimit))
		{
		  TooBig = 1;
		  printf("TOO BIG to left\n");
		}
	      /*Check that one side hasn't gone down too far*/
	      TooDeep = 0;
	      diffRside = (depth[summitlocation[f]]-depth[WavCombRight[f]]);
	      diffLside = (depth[summitlocation[f]]-depth[WavCombLeft[f] - moveL]);
	      diffratio = diffLside/diffRside;
	      if (diffratio > HeightDiffFract)
		{
		  TooDeep = 1;
		  printf("TOO DEEP to left\n");
		}
	      /*correction so that data resolution doesn't unduly affect this*/
	      if(diffLside-diffRside < difftoosmall && diffLside-diffRside > 0)
		{
		  TooDeep = 0;
		  printf(".... but too small a difference to worry about\n");
		}

	      /*if the ratio is better.....*/
	      if (ratio > maxratio && LHSgradient > mingradient && TooBig == 0 && TooDeep == 0)
		{
		  /*change maxratio*/
		  maxratio = ratio;
		  /*and move the edge of the object*/
		  printf("WavCombLeft[%d] %d\n",f,WavCombLeft[f]);
		  WavCombLeft[f] = WavCombLeft[f] - moveL;
		  printf("WavCombLeft[%d] %d\n",f,WavCombLeft[f]);
		  /*and get out of the loopright loop on first success*/
		  loopleft = 0;
		  /*and note the change so that stay in main loop*/
		  changes++;
		}

	      /*if the ratio is worse, or gradient disallows a move ......*/
	      /*and it's possible, move the limit out another step*/
	      /*Criterion 1: Can't go beyond the limits of the data*/
	      /*Criterion 2: Can explore up to XX% of the current width of the feature*/
	      /*Criterion 3: Can't go beyond the limits of a larger coefficient object*/
	      criteriaFulfilled = 1;
	      /*Criterion 1*/
	      if (WavCombLeft[f] - (moveL +1) < 0)
		{
		  criteriaFulfilled = 0;
		  printf("FAIL 1: Tries to move outside data\n");
		}
	      /*Criterion 2*/
	      /*Version from Book16p26*/
	      /*MoveMax = WavCombScale[f]*widthFract;*/
	      MoveMax = (distance[WavCombRight[f]]-distance[WavCombLeft[f]])*widthFract;
	      /*Because this is real data allow 5km whatever*/
	      if (MoveMax < LowestMaxMove)
		{
		  MoveMax = LowestMaxMove;
		}
	      Move = distance[WavCombLeft[f]]-(distance[WavCombLeft[f]-(moveL+1)]);
	      printf("width: %f MoveMax %f Move %f \n", (distance[WavCombRight[f]]-distance[WavCombLeft[f]]), MoveMax, Move);
	      if (Move > MoveMax)
		{
		  criteriaFulfilled = 0;
		  printf("FAIL 2: Ties to move by more than %f of object width\n",widthFract);
		}
	      /*Criterion 3*/
	      /*if potenital move places it within another larger feature, then disallow*/
	      for (c = 0; c < numberlines; c++)
		{
		  /*for each other feature*/
		  if (WavComb[c] > 0 && c != f)
		    {
		      /*check if the position of the moving extremity*/
		      /*is between the limits of a larger object - in terms of its original coefficient*/
		      if ((WavCombLeft[f]-(moveL+1)) < WavCombRight[c] && (WavCombLeft[f]-(moveL+1)) > WavCombLeft[c] && WavComb[c] > WavComb[f])
			{
			  printf("FAIL 3: Tries to move inside other larger object\n");
			  criteriaFulfilled = 0;
			}
		    }
		}
	      
	      

	      /*And move the left limit, or get  out of the loop as appropriate*/
	      if (criteriaFulfilled == 1)
		if (ratio <= maxratio || LHSgradient < mingradient || TooDeep == 1)
		  {
		    {
		      /*Try moving the limit out another step*/
		      moveL = moveL + 1;
		    }
		  }
	      if (criteriaFulfilled == 0 || TooBig == 1)
		{
		  loopleft = 0;
		}	    
	    }
	  
	  
	  
	  /*Try Moving the righht hand limit*/
	  /*Always move on first success*/
	  printf("Trying to move the right hand limit\n");
	  loopright = 1;
	  moveR = 1;
	  while (loopright == 1)
	    { 
	      /*initialize some properties*/
	      outline = 0;
	      area = 0;
	      
	      /*Print some details*/
	      printf("Limits are (%f,%f) (%f,%f)\n",depth[WavCombLeft[f]],distance[WavCombLeft[f]],depth[(WavCombRight[f]+moveR)],distance[(WavCombRight[f]+moveR)]);
	      
	      /*Work out the properties*/
	      /*gradient of line under feature*/
	      gradient = (depth[WavCombLeft[f]]-depth[(WavCombRight[f]+moveR)])/(distance[WavCombLeft[f]]-distance[(WavCombRight[f]+moveR)]);
	      printf("gradient is %f\n",gradient);
	      /*The length of the outline*/
	      for (i = WavCombLeft[f];i < (WavCombRight[f]+moveR); i++)
		{
		  /*over the top - use pythag*/
		  extradist = sqrt(((depth[i]-depth[i+1])*(depth[i]-depth[i+1]))+((distance[i]-distance[i+1])*(distance[i]-distance[i+1])));
		  outline = outline + extradist;
		  /*printf("[%d] outline %f extradist %f\n",i, outline, extradist);*/
		}
	      /*and the line underneath - use pythag*/
	      outline = outline + sqrt(((depth[WavCombLeft[f]]-depth[(WavCombRight[f]+moveR)])*(depth[WavCombLeft[f]]-depth[(WavCombRight[f]+moveR)]))+((distance[WavCombLeft[f]]-distance[(WavCombRight[f]+moveR)])*(distance[WavCombLeft[f]]-distance[(WavCombRight[f]+moveR)])));
	      printf("[%d] total outline (inc. base) %f \n",f, outline);
	      /*The area of the feature*/
	      /*There is no exception where the measured line crosses the extrapolated*/
	      for (j = WavCombLeft[f];j < (WavCombRight[f]+moveR); j++ )
		{
		  depthInterpL = depth[WavCombLeft[f]] + ((distance[j]-distance[WavCombLeft[f]])*gradient);
		  depthInterpR = depth[WavCombLeft[f]] + ((distance[j+1]-distance[WavCombLeft[f]])*gradient);
		  extraarea = 0.5*((depth[j]-depthInterpL) + (depth [j+1] - depthInterpR))*(distance[j+1]-distance[j]);
		  area = area + extraarea;
		  /*printf("[%d] area %f extraarea %f\n",j, area, extraarea);*/
		}
	      printf("[%d] total area %f \n",f, area);
	      ratio = area/outline;
	      printf("area/outline %f\n",ratio);
	      
	      /*print current state of movement*/
	      printf("moveR %d\n", moveR);
	      
	      /*Check that the gradient of the outside portion of the proposed move is not too flat*/
	      /*Make the proportion with respect to the summit*/
	      /*avoids going too far onto the abyssal plain*/
	      /*Right hand side*/
	      /*work out how far back towards the centre to take the outside part*/
	      RHSgradient = 0;
	      summitToEdge = (distance[(WavCombRight[f] + moveR)] - distance[summitlocation[f]]);
	      for (m = (WavCombRight[f] + moveR);(distance[(WavCombRight[f] + moveR)] - distance[m])/summitToEdge < outsideFract;m--)
		{
		  InsidePosnRHS = m;
		}
	      /*make sure that InsidePosnRHS is moved inside by at least 1 data point*/
	      if (InsidePosnRHS == (WavCombRight[f] + moveR))
		{
		  InsidePosnRHS--;
		}
	      RHSgradient = (depth[InsidePosnRHS]-depth[(WavCombRight[f] + moveR)])/(distance[(WavCombRight[f] + moveR)]-distance[InsidePosnRHS]);
	      printf("Summit location: %d\n", summitlocation[f]);
	      printf("WavCombRight[%d] %d InsidePosnRHS: %d RHSgradient %f \n",f,WavCombRight[f]+moveR,InsidePosnRHS, RHSgradient);

	      /*Make sure that the object hasn't expanded too much*/
	      TooBig = 0;
	      printf("proposed right limit: %f\n", distance[WavCombRight[f] + moveR]);
	      printf("summit %f km \n", distance[summitlocation[f]]);
	      printf("size of extension allowed %f km\n", ((WavCombScale[f]/4.0)*ScaleMulLimit) );
	      if (distance[WavCombRight[f] + moveR] > distance[summitlocation[f]] + ((WavCombScale[f]/4.0)*ScaleMulLimit))
		{
		  TooBig = 1;
		  printf("TOO BIG to right\n");
		}
	      /*Check that one side hasn't gone down too far*/
	      TooDeep = 0;
	      diffRside = (depth[summitlocation[f]]-depth[WavCombRight[f] + moveR]);
	      diffLside = (depth[summitlocation[f]]-depth[WavCombLeft[f]]);
	      diffratio = diffRside/diffLside;
	      if (diffratio > HeightDiffFract)
		{
		  TooDeep = 1;
		  printf("TOO DEEP to right\n");
		}
	      /*correction so that data resolution doesn't unduly affect this*/
	      if(diffRside-diffLside < difftoosmall && diffRside-diffLside > 0)
		{
		  TooDeep = 0;
		  printf(".... but too small a difference to worry about\n");
		}




	      /*if the ratio is better.....*/
	      if (ratio > maxratio && RHSgradient > mingradient && TooBig == 0 && TooDeep ==0)
		{
		  /*change maxratio*/
		  maxratio = ratio;
		  /*and move the edge of the object*/
		  printf("WavCombRight[%d] %d\n",f,WavCombRight[f]);
		  WavCombRight[f] = WavCombRight[f] + moveR;
		  printf("WavCombRight[%d] %d\n",f,WavCombRight[f]);
		  /*and get out of the loopright loop - 1st success*/
		  loopright = 0;
		  /*and note the change so that stay in main loop*/
		  changes++;
		}
	      
	      /*if the ratio is worse......*/
	      /*and it's possible, move the limit out another step*/
	      /*Criterion 1: Can't go beyond the limits of the data*/
	      /*Criterion 2: Can explore up to XX% of the current width of the feature*/
	      /*Criterion 3: Can't go beyond the limits of a larger coefficient object*/
	      criteriaFulfilled = 1;
	      /*Criterion 1*/
	      if (WavCombRight[f] + (moveR +1) >= numberlines)
		{
		  criteriaFulfilled = 0;
		  printf("FAIL 1: Tries to move outside data\n");
		}
	      /*Criterion 2*/
	      /*Version from book16p26*/
	      /*MoveMax = WavCombScale[f]*widthFract;*/
	      MoveMax = (distance[WavCombRight[f]]-distance[WavCombLeft[f]])*widthFract;
	      /*Because this is real data allow 5km whatever*/
	      if (MoveMax < LowestMaxMove)
		{
		  MoveMax = LowestMaxMove;
		}
	      Move = (distance[WavCombRight[f]+(moveR + 1)])-distance[WavCombRight[f]];
	      printf("width: %f MoveMax %f Move %f\n", (distance[WavCombRight[f]]-distance[WavCombLeft[f]]), MoveMax, Move);
	      if (Move > MoveMax)
		{
		  criteriaFulfilled = 0;
		  printf("FAIL 2: Ties to move by more than %f of object width\n",widthFract);
		}
	      /*Criterion 3*/
	      /*if potenital move places it within another larger feature, then disallow*/
	      for (c = 0; c < numberlines; c++)
		{
		  /*for each other feature*/
		  if (WavComb[c] > 0 && c != f)
		    {
		      /*check if the position of the moving extremity*/
		      /*is between the limits of a larger object - in terms of its original coefficient*/
		      if ((WavCombRight[f]+(moveR+1)) < WavCombRight[c] && (WavCombRight[f]+(moveR+1)) > WavCombLeft[c] && WavComb[c] > WavComb[f])
			{
			  criteriaFulfilled = 0;
			  printf("FAIL 3: Tries to move inside other larger object\n");
			}
		    }
		}
	      


	      /*And move the left limit, or get  out of the loop as appropriate*/
	      if (criteriaFulfilled == 1)
		if (ratio <= maxratio || RHSgradient < mingradient || TooDeep == 1)
		  {
		    {
		      /*Try moving the limit out another step*/
		      moveR = moveR + 1;
		    }
		  }
	      if (criteriaFulfilled == 0 || TooBig == 1)
		{
		  loopright = 0;
		}
	    }


	  /*As a last resort try moving both sides at once by one step*/
	  /*This gets round some issues of very irregular data spacing*/
	  /*So, if there haven't been any changes so far on this loop.....*/
	    if (changes == 0)
	    {
	      /*Assume that I can try this as I assume that I can move any one individually*/
	      /*initialize some properties*/
	      outline = 0;
	      area = 0;
	      
	      /*Print some details*/
	      printf("Limits are (%f,%f) (%f,%f)\n",depth[WavCombLeft[f]-1],distance[WavCombLeft[f]-1],depth[(WavCombRight[f]+1)],distance[(WavCombRight[f]+1)]);
	      
	      /*Work out the properties*/
	      /*gradient of line under feature*/
	      gradient = (depth[WavCombLeft[f]-1]-depth[(WavCombRight[f]+1)])/(distance[WavCombLeft[f]-1]-distance[(WavCombRight[f]+1)]);
	      printf("gradient is %f\n",gradient);
	      /*The length of the outline*/
	      for (i = (WavCombLeft[f]-1);i < (WavCombRight[f]+1); i++)
		{
		  /*over the top - use pythag*/
		  extradist = sqrt(((depth[i]-depth[i+1])*(depth[i]-depth[i+1]))+((distance[i]-distance[i+1])*(distance[i]-distance[i+1])));
		  outline = outline + extradist;
		  /*printf("[%d] outline %f extradist %f\n",i, outline, extradist);*/
		}
	      /*and the line underneath - use pythag*/
	      outline = outline + sqrt(((depth[WavCombLeft[f]-1]-depth[(WavCombRight[f]+1)])*(depth[WavCombLeft[f]-1]-depth[(WavCombRight[f]+1)]))+((distance[WavCombLeft[f]-1]-distance[(WavCombRight[f]+1)])*(distance[WavCombLeft[f]-1]-distance[(WavCombRight[f]+1)])));
	      printf("[%d] total outline (inc. base) %f \n",f, outline);
	      /*The area of the feature*/
	      /*There is no exception where the measured line crosses the extrapolated*/
	      for (j = (WavCombLeft[f]-1);j < (WavCombRight[f]+1); j++ )
		{
		  depthInterpL = depth[WavCombLeft[f]-1] + ((distance[j]-distance[WavCombLeft[f]-1])*gradient);
		  depthInterpR = depth[WavCombLeft[f]-1] + ((distance[j+1]-distance[WavCombLeft[f]-1])*gradient);
		  extraarea = 0.5*((depth[j]-depthInterpL) + (depth [j+1] - depthInterpR))*(distance[j+1]-distance[j]);
		  area = area + extraarea;
		  /*printf("[%d] area %f extraarea %f\n",j, area, extraarea);*/
		}
	      printf("[%d] total area %f \n",f, area);
	      ratio = area/outline;
	      printf("area/outline %f\n",ratio);
	      
	      /*print current state of movement*/
	      printf("Trying to move BOTH at once\n");

	      /*Check that the gradient of the outside portion of the proposed move is not too flat*/
	      /*Make the proportion with respect to the summit*/
	      /*avoids going too far onto the abyssal plain*/
	      /*Left hand side*/
	      /*work out how far back towards the centre to take the outside part*/
	      LHSgradient = 0;
	      summitToEdge = (distance[summitlocation[f]]-distance[(WavCombLeft[f] - 1)]);
	      for (l = (WavCombLeft[f] - 1);(distance[l]-distance[(WavCombLeft[f] - 1)])/summitToEdge < outsideFract;l++)
		{
		  InsidePosnLHS = l;
		}
	      /*make sure that InsidePosnLHS is moved inside by at least 1 data point*/
	      if (InsidePosnLHS == (WavCombLeft[f] - 1))
		{
		  InsidePosnLHS++;
		}
	      LHSgradient = (depth[InsidePosnLHS]-depth[(WavCombLeft[f] - 1)])/(distance[InsidePosnLHS]-distance[(WavCombLeft[f] - 1)]);
	      printf("Summit location: %d\n", summitlocation[f]);
	      printf("WavCombLeft[%d] %d InsidePosnLHS: %d LHSgradient %f \n",f,WavCombLeft[f]-1,InsidePosnLHS, LHSgradient);

	      /*And the right hand side*/
	      /*work out how far back towards the centre to take the outside part*/
	      RHSgradient = 0;
	      summitToEdge = (distance[(WavCombRight[f] + 1)] - distance[summitlocation[f]]);
	      for (m = (WavCombRight[f] + 1);(distance[(WavCombRight[f] + 1)] - distance[m])/summitToEdge < outsideFract;m--)
		{
		  InsidePosnRHS = m;
		}
	      /*make sure that InsidePosnRHS is moved inside by at least 1 data point*/
	      if (InsidePosnRHS == (WavCombRight[f] + 1))
		{
		  InsidePosnRHS--;
		}
	      RHSgradient = (depth[InsidePosnRHS]-depth[(WavCombRight[f] + 1)])/(distance[(WavCombRight[f] + 1)]-distance[InsidePosnRHS]);
	      printf("Summit location: %d\n", summitlocation[f]);
	      printf("WavCombRight[%d] %d InsidePosnRHS: %d RHSgradient %f \n",f,WavCombRight[f]+1,InsidePosnRHS, RHSgradient);


	      /*Test to see if object has expanded too much*/
	      criteriaFulfilled = 1;
	      if (distance[WavCombRight[f] +1] > distance[summitlocation[f]] + ((WavCombScale[f]/4.0)*ScaleMulLimit))
		{
		  criteriaFulfilled = 0;
		  printf("FAIL 4: Tries to expand too much\n");
		}
	      if (distance[WavCombLeft[f] - 1] < distance[summitlocation[f]] - ((WavCombScale[f]/4.0)*ScaleMulLimit))
		{
		  criteriaFulfilled = 0;
		  printf("FAIL 4: Tries to expand too much\n");
		}
	      
	      /*Check that one side hasn't gone down too far*/
	      TooDeep = 0;
	      diffRside = (depth[summitlocation[f]]-depth[WavCombRight[f] + 1]);
	      diffLside = (depth[summitlocation[f]]-depth[WavCombLeft[f]] - 1);
	      diffratio = diffRside/diffLside;
	      if (diffratio > HeightDiffFract || (1/diffratio) > HeightDiffFract)
		{
		  TooDeep = 1;
		  printf("TOO DEEP to one side when both moved\n");
		}
	      /*correction so that data resolution doesn't unduly affect this*/
	      if(sqrt((diffRside-diffLside)*(diffRside-diffLside)) < difftoosmall)
		{
		  TooDeep = 0;
		  printf(".... but too small a difference to worry about\n");
		}

	      /*if the ratio is better, gradient doesn't get too shallow*/
	      /*and the object hasn't got too big*/
	      if (ratio > maxratio && RHSgradient > mingradient && LHSgradient > mingradient && criteriaFulfilled == 1 && TooDeep == 0)
		{
		  /*change maxratio*/
		  maxratio = ratio;
		  /*and move the edges of the object*/
		  printf("WavCombRight[%d] %d\n",f,WavCombRight[f]);
		  WavCombRight[f] = WavCombRight[f] + 1;
		  printf("WavCombRight[%d] %d\n",f,WavCombRight[f]);
		  printf("WavCombLeft[%d] %d\n",f,WavCombLeft[f]);
		  WavCombLeft[f] = WavCombLeft[f] -1;
		  printf("WavCombLeft[%d] %d\n",f,WavCombLeft[f]);
		  /*and note the change so that stay in main loop*/
		  changes++;
		}
	      if (ratio <= maxratio || RHSgradient < mingradient || LHSgradient < mingradient || criteriaFulfilled == 0 || TooDeep == 1)
		{
		  printf("Didn't work\n");
		}
	    }

	  


	  printf("Changes %d\n", changes);
	  mainloop = 0;
	  if (changes > 0)
	    {
	      mainloop = 1;
	    }
	  /*end of the main loop*/
	}

      /*At the end of the loopobjets loop*/
      /*alter the array storing objects still to do*/
      /*firstly to remove this object*/
      objToDo[f] = 0;
      /*... and secondly to remove any objects that are inside the new range of the bigger one*/
      for (d = 0; d < numberlines; d++)
	{
	  /*if there's an object*/
	  if (objToDo[d] == 1)
	    {
	      /*find out if the centre is inside (or even at edge of) the new range*/
	      if (d >= WavCombLeft[f] && d <= WavCombRight[f])
		{
		  /*and get rid of it*/
		  printf("lhs edge at [%d] rhs edge at [%d] so REMOVED object at [%d]\n",WavCombLeft[f],WavCombRight[f],d);
		  objToDo[d] = 0;
		  WavComb[d] = 0;
		  WavCombLeft[d] = 0;
		  WavCombRight[d] = 0;
		  WavCombScale[d] = 0;
		}
	      /*find out if initial guess at edge of small object now overlaps with new range of larger object*/
	      if (WavCombLeft[d] > WavCombLeft[f] && WavCombLeft[d] < WavCombRight[f] && objToDo[d] > 0)
		{
		  /*correct the initial range of the smaller object to avoid the clash*/
		  WavCombLeft[d] = WavCombRight[f];
		  printf("ADJUSTED lhs of object at [%d]\n",d);
		}
	      if (WavCombRight[d] > WavCombLeft[f] && WavCombRight[d] < WavCombRight[f] && objToDo[d] > 0)
		{
		  /*correct the initial range of the smaller object to avoid the clash*/
		  WavCombRight[d] = WavCombLeft[f];
		  printf("ADJUSTED rhs object at [%d]\n",d);
		}
	    }
	}

      /*Count the objects left to do - removals may have resulted in none being left*/
      ObjectsCount = 0;
      for (e = 0; e < numberlines; e++)
	{
	  if (objToDo[e] > 0)
	    {
	      ObjectsCount++; 
	    }
	}
      printf("%d objects remain\n\n\n", ObjectsCount);

      /*alter the loopobjects logic*/
      loopobjects = 0;
      if (ObjectsCount >= 1 && ObjectsCount >= 1)
	{
	  loopobjects = 1;
	}

    }

  return 0;
}



/*READ_DATA*/
/*ReadData: Does all appropriate to chuck arrays back to main*/
/*Doesn't use fscanf though*/
int ReadData (char *argv[], float distance[], float depth[], float longitude[], float latitude[])

{
	int len;
	int countlines;
	char line[MAXLINE];
	float column1, column2, column3, column4;

	FILE *fdata;
	fdata = fopen (argv [1], "r");


	formatarray(distance);
	formatarray(depth);
	formatarray(longitude);
	formatarray(latitude);


	countlines = 0;
	while ((len = getlineJH(line, fdata)) > 0)
		{
/*Reads 2 columns distance along track and depth*/
		sscanf(line, "%f %f", &column1, &column2);
/* x,z in meters*/
/*printf("In main while,%f %f\n", column1, column2);*/
		distance[countlines] = column1;
		depth[countlines] = column2;
		/*printf("%f %f %d\n",  distance[countlines], depth[countlines],countlines);*/
		countlines++;
		}

	/*printf("distance[0] %f\n", distance[0]);
	  printf("distance[MAXCOLUMNS -1] %f\n", distance[MAXCOLUMNS -1]);*/

	fclose (fdata);

	/*countlines is the number of lines of data e.g. 400 filling array from 0-399*/
	return countlines;
}


/*WAVELET*/
/*wavelet: Runs a wavelet across the time series to produce a scale, position
significance file, AND more importantly uses an estimate of the hierachy
to define the limits of objects an thus produce a regional residual filter.
Optional inclusion of shape checking of objects found at each scale before
their inclusion into the interpretation*/
/*NOTE: this is not true wavelet work, more a nice easy bastardisation that
I can understand and program*/
int wavelet(float distance[], float depth[], char *argv[], int numberlines)

{
	float WavFilt[MAXCOLUMNS];	/*The array to hold the output filtered points*/
	float WavComb[MAXCOLUMNS];	/*Holds the max approved coeffs from across the scales*/
	float WavCombScale[MAXCOLUMNS];	/*Holds the scale at which the coeffs above come from*/
	int WavCombLeft[MAXCOLUMNS];    /*Holds the array element of the left hand side extremity of each feature*/
	int WavCombRight[MAXCOLUMNS];   /*Holds the array element of the right hand side extremity of each feature*/

	float WavCoeff[MAXCOLUMNS];	/*Coeffs of wavelet transform at a given scale*/
	float WavWeighting[MAXCOLUMNS];	/*holds the weighting for the wavCompute
					changes every time it's called*/


	float scale;			/*scale currently being opperated at*/
	int i, j, k, l, m, n, p, q, r;			/*loop variable*/
	float CoeffThreshold;		/*Detection threshold for objects - same for all scales*/
	float maxscale;			/*max scale upon which to search*/
	float NumFeatures;               /*final output of number of features found*/
	float ScaleInteract;		 /*factor which governs the span of scales that are in direct competition*/

	FILE *frawCoeff;		/*Uninterpreted Coefficients file*/
	FILE *fwavObjCoeff;		/*Interpreted Coefficients file - objects as selected max coeff, all scales sequentially*/
	FILE *fwavFilt;			/*Filtered depths according to the wavelet method*/
	FILE *fwavObjCoeffComb;		/*Interpreted Coefficients file - max coeffs selected from across the scales*/
	FILE *fwavBlockedCoeff;		/*Interpreted Coefficients summed and blocked in*/
	FILE *fwavNum;	                /*count of number of features*/

/*Set things*/
	scale = 0.1;			/*The smallest scale - starting point*/
	CoeffThreshold = 0.0;		/*The thing that objects need to get above to
					  be objects*/
	maxscale = 3000;			/*[300] Maximum width of wavelet used*/
	ScaleInteract = 2;            /*Set a ScaleInteract factor which governs the span of scales that are in direct competition*/

/*Open the output files*/
	frawCoeff = fopen (argv [2], "w");
	fwavObjCoeff = fopen (argv [3], "w");
	fwavFilt = fopen (argv [4], "w");
	fwavObjCoeffComb = fopen (argv[5], "w");
	fwavBlockedCoeff = fopen (argv[6], "w");
	fwavNum = fopen (argv[7], "w");

/*Initialise the new arrays to be used for the multi-scale output*/
/*The other three are reinitialised at every scale in the first functions to call them*/
	formatarray(WavFilt);
	formatarray(WavComb);
	formatarray(WavCombScale);
	formatarrayI(WavCombLeft);
	formatarrayI(WavCombRight);

	while(scale <= maxscale)
	  {
	    printf("\nCurrent scale = %.1f km\n", scale);
	    
	    /*Calculate the wavelet coefficients for this scale - my INITIAL wavelet*/
	    wavCompute(distance, depth, WavCoeff, scale, WavWeighting, numberlines);
	    /*And send to file*/
	    /*printf("Printing The Raw Coefficients To File\n");*/
	    for(i = 0; i < numberlines; i++)
	      {
		/*printf("%d %f %f %f\n",i, distance[i], scale, WavCoeff[i]);*/
		fprintf(frawCoeff, "%f %f %f\n", distance[i], scale, WavCoeff[i]);
	      }

	    /*Figure out the positions of the objects for this scale*/
	    /*This reduces the information in WavCoeff[j] to the best estimate of the position of the objects*/
	    wavInterpii(WavCoeff, numberlines, CoeffThreshold, distance, scale);
	    /*And send to file*/
	    for(j = 0; j < numberlines; j++)
	      {
		if (WavCoeff[j] > 0)
		  {
		    fprintf(fwavObjCoeff, "%d %f %f %f\n", j, scale, WavCoeff[j], distance[j]);
		  }
	      }

	    /*A more efficient way of handling the scales*/
	    if (scale >= 0.1 && scale < 1)
	      scale = scale + 0.1;
	    if (scale >= 1 && scale < 25)
	      scale = scale + 1;
	    if (scale >= 25 && scale < 100)
	      scale = scale + 5;
	    if (scale >= 100 && scale < 500)
	      scale = scale + 10;
	    if (scale >= 500 && scale < 5000)
	      scale = scale + 100;
	    if (scale >= 5000 && scale < 50000)
	      scale = scale + 1000;
	    if (scale >= 50000 && scale < 500000)
	      scale = scale + 10000;
	    if (scale >= 500000)
	      scale = scale + 100000;
	  }
	/*Close the output file with all the raw coefficients in, and the one where only object locations remain*/
	fclose(frawCoeff);
	fclose(fwavObjCoeff);


	
	/*Combine the locations of objects across the scale*/
	printf("Combining information from all the scales........\n");
	wavInterpCombii(argv, WavComb, numberlines, distance, WavCombScale, ScaleInteract);

	/*The final product of wavInterpCombii - send to file*/
	for(p = 0; p < numberlines; p++)
	  {
	    if (WavComb[p] > 0)
	      {
		fprintf(fwavObjCoeffComb, "%f %f %f %f\n", distance[p], WavCombScale[p], WavComb[p], depth[p]);
	      }
	  }


/*Turn the coefficients into a regional by finding the objects and doing linear
 interpolation between their limits*/
 printf("\nCreating a regional filter from the objects at all scales\n");
 wavInterpComb(WavComb, numberlines, WavFilt, depth, distance, WavCombScale, WavCombRight, WavCombLeft);

/*Do the second stage of processing - finding the exact limits of the objects*/
 PostProcess(WavComb, numberlines, WavFilt, depth, distance, WavCombScale, WavCombRight, WavCombLeft);

 /*Re-interpolate underneath the objects*/
 Interpolate(WavComb, numberlines, WavFilt, depth, distance, WavCombScale, WavCombRight, WavCombLeft);



/*Got the blocked objects, now do some post processing to make their limits astetically acceptable*/
/*Firstly sort out any times when the cf pp 178 Book 2*/
     /*WavPostPro(distance, depth, WavFilt, WavComb, numberlines, CoeffThreshold, WavCombScale);*/
     /*WavPostProii(distance, depth, WavFilt, WavComb, numberlines, CoeffThreshold, WavCombScale, ScaleInteract);*/
     /*WavPostPro(distance, depth, WavFilt, WavComb, numberlines, CoeffThreshold, WavCombScale);*/

/*Do I need to use WavPostPro again to check that WavPostProii hasn't messed it up????????????*/
/*And send Filtered to file*/
		for(l = 0; l < numberlines; l++)
			fprintf(fwavFilt, "%f %f %f\n", distance[l], depth[l], WavFilt[l]);

/*Finally, how many objects have I got?*/
		NumFeatures = 0;
		for(r = 0; r < (numberlines -1); r++)
		  {
		    if (depth[r] == WavFilt[r] && depth[r+1] != WavFilt[r+1])
		      {
			NumFeatures++;
		      }
		  }
		fprintf(fwavNum, "%f\n", NumFeatures);
/*Close all files opened*/
	fclose(fwavFilt);
	fclose(fwavObjCoeffComb);
	fclose(fwavBlockedCoeff);
	fclose(fwavNum);

	return 0;
}


/*WAV_COMPUTE*/
/*wavCompute: For every point along a profile calculates a weighted average
according to my origional wavelet*/
int wavCompute(float distance[], float depth[], float WavCoeff[], float scale, float wavWeighting[], int numberlines)

{

	int i, j, k ,l;				/*loop variables*/
	int before, after;			/*Number position of the points at the
						start and end of the wavelet*/
	int toprint;                            /*printing logical variable - debugging*/
	float BegSum, MidSum, EndSum;		/*weighted averages for the three section of the wavelet*/
	float BegCount, MidCount, EndCount;	/*count the points for the weighted average*/
	float BegAv, MidAv, EndAv;              /*the weighted average*/
	float Coeff;				/*Weighted average at a scale or wavelet coeff of scale, position*/
	int bad;				/*To be triggered by REQUIREMENT 1 if it sets a coeff to zero
						so that further requirements don't change it back again
						0 = bad not active i.e. good, 1 = bad already ditched*/
	float diffleft, diffright;		/*REQUIREMENT 2: differences between the middle and
						right and left hand sides*/
	float diffFactor;			/*Factor for the difference between the two sides*/
	float ReqIIrecalc;			/* 0=don't need to recalculate coeff
						 1 = do need to recalculate coeff*/

/*Initialise the per cycle arrays first used here*/
	formatarray(wavWeighting);
	formatarray(WavCoeff);

/*Initialise any variables that obviously aren't in the code anywhere*/
	Coeff = 0;

	/*printf("wavCompute: Calcuating the wavelet coefficients\n");*/

	printf("numberlines %d\n", numberlines);
/*Use the i loop to circle through all the points along track in turn*/
	for(i = 0; i < numberlines; i++)
	  {

	    /*variable to print out a selected computation*/
	    toprint = 0;
	    if (i == 170)
	      toprint = 0;
	    /*change toprint to 1 in order to print things to help debugging*/

/*Initialize logical parameters for the requirements at the end of this loop*/
	    bad = 0;
	    ReqIIrecalc = 0;

/*Select array values of interest from within +- half scale*/
/*[before] is the first array element of interest, and [after] is the last*/
	    before = 0;
	    after = 0;
/*start at the center point and keep going back along track until you've gone
back a length of half the scale, or hit the begining of the data*/
	    for(k = i; distance[k] > (distance[i] - (scale/2.0)) && k >= 0; k--)
	      {
		before = k;
		if (toprint == 1)
		  {
		    printf("distance[k] %f distance[i] %f scale/2.0 %f\n", distance[k], distance[i], (scale/2.0));
		    printf("i %d before %d\n", i, before);
		  }
	      }
/*Same but forward*/
	    for(l = i; distance[l] < (distance[i] + (scale/2.0)) && l < numberlines; l++)
	      {
		after = l;
		if (toprint == 1)
		  {
		    printf("distance[l] %f distance[i] %f scale/2.0 %f\n", distance[l], distance[i], (scale/2.0));
		    printf("i %d after %d\n", i, after);
		  }
	      }
	/*printf("position(i) %d before %d after %d\n", i, before, after);*/


/*Work out the weighted average that is the coefficient*/
/*To cope with unequal data density this needs to be done in
the three sections independantly -  ungainly though it may seem*/
	    /*These variables are counts and sums for the 3 parts of the wavelet*/
		BegCount = 0;
		BegSum = 0;
		/*printf("Initial values: BegCount %f beginning %f\n", BegCount, beginning);*/
		MidCount = 0;
		MidSum = 0;
		/*printf("Initial values: MidCount %f middle %f\n", MidCount, middle);*/
		EndCount = 0;
		EndSum = 0;
		/*printf("Initial values: EndCount %f end %f\n", EndCount, end);*/
		for (j = before; j <= after; j++)
			{
/*In the first quarter the weighting is -1.0*/
/*Here computed as a mean, but a sum of wavWeighting[j]*depth[j] is the general way of doing it*/
			if (distance[j] > (distance[i] - (scale/2.0)))
				if (distance [j] < (distance[i] - (scale/4.0)))
					{
					wavWeighting[j] = -1.0;
					BegCount = BegCount + 1;
					BegSum = BegSum + depth[j];
					if (toprint == 1)
					  {
					    printf("j %d\n", j);
					    printf("Beg: wavWeighting[j] %f depth[j] %f\n", wavWeighting[j], depth[j]);
					    printf("BegCount %f beginning %f\n", BegCount, BegSum);
					  }
					}
/*In the middle half the weighting is 1.0*/
			if (distance[j] >= (distance[i] - (scale/4.0)))
				if (distance [j] <= (distance[i] + (scale/4.0)))
					{
					wavWeighting[j] = 1.0;
					MidCount = MidCount + 1;
					MidSum = MidSum + depth[j];
					if (toprint == 1)
					  {
					    printf("j %d\n", j);
					    printf("Mid: wavWeighting[j] %f depth[j] %f\n", wavWeighting[j], depth[j]);
					    printf("MidCount %f middle %f\n", MidCount, MidSum);
					  }
					}
/*In the third quarter the weighting is -1.0*/
			if (distance[j] > (distance[i] + (scale/4.0)))
				if (distance [j] < (distance[i] + (scale/2.0)))
					{
					wavWeighting[j] = -1.0;
					EndCount = EndCount + 1;
					EndSum = EndSum + depth[j];
					if (toprint == 1)
					  {
					    printf("j %d\n", j);
					    printf("End: wavWeighting[j] %f depth[j] %f\n", wavWeighting[j], depth[j]);
					    printf("endcount %f end %f\n", EndCount, EndSum);
					  }
					}
			}
		BegAv = BegSum/BegCount;
		MidAv = MidSum/MidCount;
		EndAv = EndSum/EndCount;
		if (toprint == 1)
		  {
		    printf("BegCount %f beginning %f av %f\n", BegCount, BegSum, BegAv);
		    printf("MidCount %f middle %f av %f\n", MidCount, MidSum, MidAv);
		    printf("EndCount %f end %f av %f\n", EndCount, EndSum, EndAv);
		  }

/*Actually work out the weighted average*/
/*To cope with unequal data density this needs to be done in
the three sections independantly*/
/*times 2.0 is so that the middle can cancel out both the beginning and the end*/
		Coeff = (-1*(BegSum/BegCount)) + ((2.0*MidSum)/MidCount) + (-1*(EndSum/EndCount));
		if (toprint == 1)
		  {
		    printf("Center point of wavelet(distance) %f wavelet coefficient %f\n", distance[i], Coeff);
		  }

/*Some spatial criteria so that only seamounts are found Book2 pp161-165*/
/*REQUIRMENT 1: That the edges must be lower than the middle*/
/*REMEMBER: depth is thought of as height with +ve up*/
		if(BegAv > MidAv)
			{
			
			  bad = 1;
			  Coeff = 0;
			  if (toprint == 1)
			    {
			      printf("SIDE TOO HIGH\n");
			    }
			}
		if(EndAv > MidAv)
			{
			  bad = 1;
			  Coeff = 0;
			  if (toprint == 1)
			    {
			      printf("SIDE TOO HIGH\n");
			    }
				
			}

/*REQUIREMENT 2: That if one side is much lower than the other one with respect to the
middle then a large value is still not produced*/

/*Set up as if the difference one side is > Factor times the difference the other side
then replace reduce the difference and recalculate the magnitude of the coefficient*/
		diffFactor = 1.0;		/*user set*/
		if(!bad && Coeff > 0)
		  {
/*Set up to get a positive difference if middle is highest*/
		    diffleft = sqrt((MidAv - BegAv)*(MidAv - BegAv));
		    diffright = sqrt((MidAv - EndAv)*(MidAv - EndAv));
		    if (toprint == 1)
		      {
			printf("SEEING IF NEED TO ADJUST SIZE OF DIFFERENCES\n");
			printf("BegAv %f MidAv %f EndAv %f\n",BegAv,MidAv,EndAv);
			printf("diffleft %f diffright %f \n",diffleft,diffright);
		      }
		    if (diffleft > (diffFactor*diffright))
		      {
			diffleft = diffright*diffFactor;
			ReqIIrecalc = 1;
			if (toprint == 1)
			  {
			    printf("LEFT adjusted\n");
			    printf("diffleft %f diffright %f \n",diffleft,diffright);
			  }
		      }
		    if (diffright > (diffFactor*diffleft))
		      {
			ReqIIrecalc = 1;
			diffright = diffleft*diffFactor;
			if (toprint == 1)
			  {
			    printf("RIGHT adjusted\n");
			    printf("diffleft %f diffright %f \n",diffleft,diffright);
			  }
		      }
		    /*If required recalculate the coefficient*/
		    /*but only if it hasn't been set to zero in REQUIREMENT 1*/
		    if (ReqIIrecalc == 1 && bad == 0)
		      {
			/*printf("Recalculating Coefficient in REQUIREMENT 2\n");*/
			/*NOTE: this calculation is equivalent to the one above*/
			Coeff = (diffright + diffleft);
			if (toprint == 1)
			  {
			    printf("Recalculating coeff\n");
			    printf("Coeff %f \n",Coeff);
			  }
		      }
		  }

/*And tuck the coefficient up nice and safe in an array*/
		WavCoeff[i] = Coeff;
		}

	return 0;
}







/*WAV_INTERPRET II*/
/*wavInterp: For one scale gets the array of raw coefficients and decides
where the objets are.  Changed to a slower but robust method*/
int wavInterpii(float WavCoeff[], int numberlines, float CoeffThreshold, float distance[], float scale)

{
  int j, k, l;	/*loop variables*/
  int pre;	/*place to start looking for conflicting variables*/
  
  
/*Initialise any variables that aren't obviously else where*/
  k = 0;
  l = 0;

/*Firstly set all Coeffs that are zero or below the threshold
to zero*/
  for (j = 0; j < numberlines; j++)
    {
      if (WavCoeff[j] < CoeffThreshold || WavCoeff[j] < 0)
	{
	  WavCoeff[j] = 0;
	}
    }
  
/*Then remove any conflicting coefficients*/
/*Set any smaller coefficient to zero*/
for (j = 0; j < numberlines; j++)
  if (WavCoeff[j] > 0)
    {
      /*find the point furthest with the lowest dist that could interact with WavCoeff[j]*/
      for (k = j; k >= 0 && distance[k] > (distance[j] - (scale/4.0)); k--)
	{
	  pre = k;
	}
      /*loop through all the other coefficients it could interract with*/
      for (l = pre; l < numberlines && distance[l] < (distance[j] + scale/4.0);l++)
	{
	  /*If WavCoeff[j] is smaller, set it to 0 and get out of the loop*/
	  if (WavCoeff[j] >=  WavCoeff[l] && j != l)
	    {
	      WavCoeff[l] = 0;
	    }
	  /*If WavCoeff[j] is larger, set the other one to 0*/
	  if (WavCoeff[j] < WavCoeff[l] && j != l)
	    {
	      WavCoeff[j] = 0;
	    }
	}
    }

	return 0;
}

/*WAV_INTERPRET_COMBINE*/
/*wavInterpComb: Finds the limits of the objects, and creates a filtered value.
My initial go at this will just be to find the non-zero patches of
WavComb*/
int wavInterpComb(float WavComb[], int numberlines, float WavFilt[], float depth[], float distance[], float WavCombScale[], int WavCombLeft[], int WavCombRight[])

{
  int a,b,c,d,e,f,g,k,h,i,j;		/*loop variable*/
  float gradient;		/*used in the linear interpolation of the space under objects*/
  float area, outline;          /*properties of the feature*/
  float extradist, extraarea;
  float ratio, ratiomax;        /*ratio: outline length / Area */
  float depthInterpL, depthInterpR;    /*interpolated depths used in calculations*/
  int countobjects;

  /*initialise stuff*/
  gradient = 0;
  countobjects = 0;

  printf("Finding the limits of the objects\n");
  /*Find the limits of the objects - store the location at the limit of the feature*/
  /*i.e. where height of the feature is zero*/
  /*NOTE: there can be no overlap at this stage*/
  for (k = 0; k < numberlines; k++)
    {
      /*printf("WavComb[%d] %f\n",k, WavComb[k]);*/
      if (WavComb[k] > 0)
	{
	  countobjects++;
	  /*Find the left hand side limit*/
	  for (a = k; a >=0 && distance[a] > distance[k] - (WavCombScale[k]/4.0); a--)
	    {
	      /*printf("Scale/4: %f\n",(WavCombScale[k]/4.0));*/
	      WavCombLeft[k] = a;
	      /*printf("Distance of left limit [%d]: %f\n",a,distance[a]);*/
	    }
	  /*If the lhs hasn't moved, move it by one as this must be OK if a coefficient can be calcuated*/
	  if (WavCombLeft[k] == k)
	    {
	      WavCombLeft[k]--;
	    }
	  /*Find the right hand side limit*/
	  for (b = k; b < numberlines && distance[b] < distance[k] + (WavCombScale[k]/4.0); b++)
	    {
	      /*printf("Scale/4: %f\n",(WavCombScale[k]/4.0));*/
	      WavCombRight[k] = b;
	      /*printf("Distance of right limit[%d]: %f\n",b, distance[b]);*/
	    }
	  /*If the rhs hasn't moved, move it by one as this must be OK if a coefficient can be calcuated*/
	  if (WavCombRight[k] == k)
	    {
	      WavCombRight[k]++;
	    }
	  /*printf("Left limit: %d, %f km\n", WavCombLeft[k], distance[WavCombLeft[k]]);*/
	  /*printf("Right limit: %d, %f km\n",WavCombRight[k], distance[WavCombRight[k]]);*/
	}
    }
  printf("No. OBJECTS is %d \n\n", countobjects);


  printf("Double-check that there is no overlap\n");
  /*For each of the non-zero coefficients indicating an object ........*/
  for (c = 0; c < numberlines; c++)
    {
      if (WavComb[c] > 0)
	{
	  /*go through all the other coefficients ....*/
	  for (d = 0;d < numberlines; d++)
	    {
	      /*And flag if any other left hand limit is between the limits for this object*/
	      if (WavCombLeft[d] >= WavCombLeft[c] && WavCombLeft[d] <= WavCombRight[c] && WavComb[d] > 0 && d != c)
		{
		  printf("WARNING - whilst centre points not in each other's ranges, ranges overlap\n");
		  printf("Left limit located at [%d] is inside limits of object at [%d]\n",d,c);
		  printf("%f (lhs of right object) %f (rhs of left object)\n", distance[WavCombLeft[d]], distance[WavCombRight[c]]);
		}
	      /*And flag if any other right hand limit is between the limits for this object*/
	      if (WavCombRight[d] >= WavCombLeft[c] && WavCombRight[d] <= WavCombRight[c] && WavComb[d] > 0 && d != c)
		{
		  printf("WARNING - whilst centre points not in each other's ranges, ranges overlap\n");
		  printf("Right limit located at [%d] is inside limits of object at [%d]\n",d,c);
		  printf("%f (rhs of left object) %f (lhs of right object)\n", distance[WavCombRight[d]], distance[WavCombLeft[c]]);
		}
	    }

	}
    }

 

  printf("\nDoing initial interpolation.......\n");

  /*Start of by setting WavFilt[] to depth[]*/
  for (e = 0; e < numberlines; e++)
    {
      WavFilt[e] = depth[e];
      /*printf("WavFilt[%d] %f\n",e,WavFilt[e]);*/
    }

  /*In the objects, interpolate between the limits*/
  for (f = 0; f < numberlines; f++)
    {
      if (WavComb[f] > 0)
	{
	  /*Gradient for the line under the object*/
	  /*printf("Extrapolate under object centred on [%d]\n",f);*/
	  /*printf("Limits are (%f,%f) (%f,%f)\n",depth[WavCombLeft[f]],distance[WavCombLeft[f]],depth[WavCombRight[f]],distance[WavCombRight[f]]);*/
	  gradient = (depth[WavCombLeft[f]]-depth[WavCombRight[f]])/(distance[WavCombLeft[f]]-distance[WavCombRight[f]]);
	  /*printf("gradient is %f\n",gradient);*/
	  /*And do the interpolation*/
	  for (g = WavCombLeft[f]; g <= WavCombRight[f]; g++)
	    {
	      WavFilt[g] = depth[WavCombLeft[f]] + ((distance[g]-distance[WavCombLeft[f]])*gradient);
	    }
	}
    }

  printf("Calculating properties of the object\n");
  /*In the objects, interpolate between the limits*/
  for (h = 0; h < numberlines; h++)
    {
      if (WavComb[h] > 0)
	{
	  outline = 0;
	  area = 0;
	  gradient = (depth[WavCombLeft[h]]-depth[WavCombRight[h]])/(distance[WavCombLeft[h]]-distance[WavCombRight[h]]);
	  /*The length of the outline*/
	  for (i = WavCombLeft[h];i < WavCombRight[h]; i++)
	    {
	      /*over the top - use pythag*/
	      extradist = sqrt(((depth[i]-depth[i+1])*(depth[i]-depth[i+1]))+((distance[i]-distance[i+1])*(distance[i]-distance[i+1])));
	      outline = outline + extradist;
	      /*printf("[%d] outline %f extradist %f\n",i, outline, extradist);*/
	    }
	  /*and the line underneath - use pythag*/
	  outline = outline + sqrt(((depth[WavCombLeft[h]]-depth[WavCombRight[h]])*(depth[WavCombLeft[h]]-depth[WavCombRight[h]]))+((distance[WavCombLeft[h]]-distance[WavCombRight[h]])*(distance[WavCombLeft[h]]-distance[WavCombRight[h]])));
	  printf("[%d] total outline (inc. base) %f \n",i, outline);
	  /*The area of the feature*/
	  /*There is no exception where the measured line crosses the extrapolated*/
	  for (j = WavCombLeft[h];j < WavCombRight[h]; j++ )
	    {
	      depthInterpL = depth[WavCombLeft[h]] + ((distance[j]-distance[WavCombLeft[h]])*gradient);
	      depthInterpR = depth[WavCombLeft[h]] + ((distance[j+1]-distance[WavCombLeft[h]])*gradient);
	      extraarea = 0.5*((depth[j]-depthInterpL) + (depth [j+1] - depthInterpR))*(distance[j+1]-distance[j]);
	      area = area + extraarea;
	      /*printf("[%d] area %f extraarea %f\n",j, area, extraarea);*/
	    }
	  printf("[%d] total area %f \n",j, area);
	  ratio = area/outline;
	  printf("area/outline %f\n",ratio);
	}
    }

	return 0;
}


/*WAV_INTERP_COMINE II*/
/*wavInterpCombii: combines information previously printed to file*/
int wavInterpCombii(char *argv[], float WavComb[], int numberlines, float distance[], float WavCombScale[], float ScaleInteract)

{
  int a,b,c,d,e,f,g, h;		/*loop variables*/
  int len;
  int linenum;
  char line[MAXLINE];
  float column1, column2, column3;
  FILE *fdata;		/*Interpreted Coefficients file - objects as selected max coeff, all scales sequentially*/
  int Arr1[MAXCOLUMNS];       /*arrays to read the file into - local*/
  float Arr2[MAXCOLUMNS];
  float Arr3[MAXCOLUMNS];
  float distdiff;           /*distance between 2 coefficients*/
  int PSFlogic, Rlogic;     /*whether or not the scales interract*/
  float PairScaleFactor;    /*difference in the scales between 2 coefficients*/
  float PairFlatnessFactor;
  float percent;
  FILE *fouttemp;
  FILE *fouttempB;
  float FlatnessFactor;      /*factor to determine what is too flat*/
  float FlatnessCoeffFactor; /*factor within which wavelet coeffs must be for pointiness to come into play*/
  float wide, narrow;

  /*As well as ScaleInteract, take out large features if excessively flat*/
  FlatnessFactor = 3.0;          /*smaller object has to be more than this many times pointier*/
  FlatnessCoeffFactor = 4.0;     /*if the smaller one has a coefficient within this multiple, it can take it out */

  /*Open the file for reading*/
  fdata = fopen (argv [3], "r");
  fouttemp = fopen ("look.temp", "w");
  fouttempB = fopen ("lookB.temp", "w");

  formatarrayI(Arr1);
  formatarray(Arr2);
  formatarray(Arr3);

  /*read the file into 3 arrays*/
  linenum = 0;
  while ((len = getlineJH(line, fdata)) > 0)
    {
      /*Reads 3 columns array location, scale, ceoefficient*/
      sscanf(line, "%f %f %f", &column1, &column2, &column3);
      Arr1[linenum] = column1;
      Arr2[linenum] = column2;
      Arr3[linenum] = column3;
      linenum++;
    }

  fclose (fdata);

  /*Now go through the arrays and sort out which points should be eliminated*/
  /*Go in reverse order through the array so that i'm taking the largest scale first*/
  for (a = linenum -1; a >= 0; a--)
    {
      /*printf("Coeff to check: %f %f %f\n", distance[Arr1[a]], Arr2[a], Arr3[a]);*/
      /*For each coefficient [a] remove it if it has a smaller coefficient than anything within its range*/
      for (b = 0; b < linenum && Arr3[a] != 0; b++)
	{
	  Rlogic = 0;
	  distdiff = sqrt((distance[Arr1[a]] - distance[Arr1[b]])*(distance[Arr1[a]] - distance[Arr1[b]]));
	  if (distdiff < Arr2[a]/4.0)
	    {
	      Rlogic = 1;
	      /*printf("Within range\n");*/
	    }
	  /*Now see if [a] should be eliminated by .....*/
	  /*1) Something thinner having a larger coefficient*/
	  if (Arr3[a] < Arr3[b] && Rlogic == 1 && a != b)
	    {
	      /*printf("%f %f %f vs %f %f %f\n", distance[Arr1[a]], Arr2[a], Arr3[a], distance[Arr1[b]], Arr2[b], Arr3[b]);*/
	      /*printf("coeff[a] %f is a smaller coeff than [b] %f \n\n", Arr3[a], Arr3[b]);*/
	      Arr3[a] = 0;
	    }
	  /*2) Something thinner being a decent size (e.g. 1/4 of coeff of larger scale) and substantially pointier*/
	  /*Flatnesses > 1 imply that the wider coefficient is flatter*/
	  wide = (Arr3[a]/Arr2[a]);
	  narrow = (Arr3[b]/Arr2[b]);
	  PairFlatnessFactor = narrow/wide;
	  /*The thinner feature must still have a decent size coefficient e.g. > a quarter of the wider one*/
	  if (Arr3[a] != 0 && Arr3[a] < ((Arr3[b])*FlatnessCoeffFactor) && Rlogic == 1 && a != b && PairFlatnessFactor > FlatnessFactor)
	    {
	      Arr3[a] = 0;
	      /*printf("[a] eliminated by small and pointy\n");
		printf("PairFlatnessFactor %f\n", PairFlatnessFactor);*/
	    }
	}



      /*If [a] is not to be eliminated, do some elimination*/
      if (Arr3[a] > 0)
	{
	  /*printf("NOT ELIMINATED %f %f %f\n", distance[Arr1[a]], Arr2[a], Arr3[a]);*/
	  for (c = 0; c < linenum; c++)
	    {
	      Rlogic = 0;
	      distdiff = sqrt((distance[Arr1[a]] - distance[Arr1[c]])*(distance[Arr1[a]] - distance[Arr1[c]]));
	      if (distdiff < Arr2[a]/4.0)
		{
		  Rlogic = 1;
		  /*printf("In Range: %f %f %f\n",distance[Arr1[c]], Arr2[c], Arr3[c]);*/
		}
	      /*Now, given the previous logic [c] should be eliminated if*/
	      /*in range with a smaller coefficient*/
	      /*If [c] should survive as pointy, it should have already eliminated [a]*/
	      /*Happens if [a] is nice and pointy and coeff[c] is smaller*/
	      if (Arr3[a] > Arr3[c] && a != c && Rlogic == 1 && Arr3[c] > 0)
		{
		  Arr3[c] = 0;
		  /*printf("And eliminate\n");*/
		}
	    }
	}
    }	

  /*See what this has done*/
  for (d = 0; d < linenum; d++)
    {
      if (Arr3[d] > 0)
	{
	  fprintf(fouttemp, "%f %f %f \n",distance[Arr1[d]],Arr2[d],Arr3[d]);
	}
    }
   
  
  


  /*Write the remaining data into WavComb[] and WavCombScale[]*/
  for (e = 0; e < linenum; e++)
    {
      if (Arr3[e] > 0)
	{
	  WavComb[Arr1[e]] = Arr3[e];
	  WavCombScale[Arr1[e]] = Arr2[e];
	}
    }


  fclose (fouttemp);
  fclose (fouttempB);

	return 0;
}





/*end of program*/






