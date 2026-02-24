/*Find_Smts: Interpreting bathymetry data in terms of frustums*/
/*05.07.04*/

/* CITATION */
/* If you use this code please cite the following paper */
/*Hillier, J. K. (2006) Seamount volcanism in space and time. Geophys. J. Int. 168(2) p877-889 doi:10.1111/j.1365-246-X.2006.03250.x. */

/*COMPILATION*/
/*Compile using    cc -V -lm Find_Smts_X.c -o Find_Smts_X.x*/

/*USAGE*/
/*Takes a line of bathymetry and it's filtered depth*/
/*and finds the objects and analyses them*/
/*input is file.dtfxy   d = dist t = topo f= filtered topo x = long y = lat*/
/*Find_Smts_X.x  input-file  -Options*/
/*Use input file with distance in km and depth (-ve values) in m*/
/*		-DO	Output File for display of frustums					*/
/*		-PO	Output File for the properties of frustums IN THE FOLLOWING ORDER  	*/
/*
			-No	- Number of the object
			-Ma	- Distance to shallowest point (km)
				- Distance to first and last points in object (km)
				- Long and Lat of the shallowest point (km)
			-Ad	- Dist (km) and depth (m) of four points of the
				  frustum: 0 = bottom left, 1 = top left
				  2 = top right, 3 = bottom right
			-Ar	- Basal radius and top radius of the frustum
				  (calculated as along track diameter/2.0)
			-As	- Slope (degrees) - av of 2 sides of frustum
			-Ah	- Depth to top and base of frustum (av of corners)
				  and pedestal height (difference between the two)
			-Md	- No. of data points in feature
			-Am	- Misfit of best fitting frustum
			-Al	- as Ad but longs and lats (decimal degrees)
			-Ac	- Long and lat of centre of frustum top
				  (av of points 1 and 2)
			-Av	- Measured Area (Trapezium approx from measured data)
				  Frustum Area (calculated from top/base radius and pedestal height)
				  Frustum Vol (calculated from top/base radius and pedestal height)
				  percent (Frustum Area/Measured Area)
			-TD	- Outputs 2 non-linearity parameters Book8pp135

												*/



/*LIBRARIES*/

#include <stdio.h>
#include <math.h>


/*DEFINITIONS*/
#define MAXOBJECTS 10000  /*Max number of seamount objects*/
#define MAXCOLUMNS 50000  /*Max number of locations in file*/
#define MAXTEMP 1000  /*Max number of locations to be in object temp array*/
#define PI 3.14159265358979323846
#define MVLOC1 3
#define MVLOC2 8
#define PRINTOBJECT 1

/*THE FUNCTION DEFINITIONS*/

int cleanfiles(int argc, char *argv[]);
char *CommandLine (int argc, char *argv[], char CommandPatterni, char CommandPatternii, int *found);
int DoMisfitFunction(int argc, char *argv[], int f, int g, float MoveLocnX[MVLOC1][MVLOC2], float MoveLocnY[MVLOC1][MVLOC2], int point, float PointCopyX[],float PointCopyY[], float FracDist[], float FracDepth[], int begin, int finish, int objectnumber, float MoveMisfit[MVLOC1][MVLOC2],float MotionUnit);
int findchunk(float distance[], float depth[], float filtered[], int start, int  numberlines);
int frustum(float distance[], float depth[], float FracDist[],  float FracDepth[], float filtered[], float longitude[], float latitude[], int begin, int finish, int argc, char *argv[], int objectnumber);
int object (int argc, char *argv[], int numberlines, float distance[], float depth[], float filtered[], float longitude[], float latitude[], float startDist[], float endDist[], int startEl[], int endEl[], int numberpoints[]);
int output(float distance[],float depth[],float longitude[], float latitude[], int begin,int finish,float PointX[],float  PointY[], float MinDepth, float MaxDepth, int argc, char *argv[], int objectnumber, int MinDepthLocn, float MISFIT, float FracDist[], float FracDepth[]);
float MisfitFunctionII(float PointX[], float PointY[], float FracDist[], float FracDepth[], int begin, int finish, int objectnumber, float ResampInc);
int normalise(int argc, char *argv[], float distance[], float depth[], float FracDist[], float FracDepth[], int begin, int finish, float MaxDepth, float MinDepth, int objectnumber);
int printarrayFFFFF(float a[], float b[], float c[], float d[], float e[]);
int ReadData (char *argv[], float distance[], float depth[], float filtered[], float longitude[], float latitude[]);
float resample(int argc, char *argv[], float FracDist[], float FracDepth[], float FracDistRS[], float FracDepthRS[], float NumRS,int begin, int finish);


/*MAIN*/
int main(int argc, char *argv[])
{
/*Data arrays*/
	float distance[MAXCOLUMNS] = {0.0};		/*Dist along track (km)*/
	float depth[MAXCOLUMNS] = {0.0};		/*Depth (m) -ve is down*/
	float filtered[MAXCOLUMNS] = {0.0};		/*Filtered depth (m) -ve is down*/
	float longitude[MAXCOLUMNS] = {0.0};		/*Longitude (decimal degrees)*/
	float latitude[MAXCOLUMNS] = {0.0};		/*Latitude (decimal degrees)*/
/*results arrays*/
	float startDist[MAXOBJECTS] = {0.0};		/*Distance along track to start of object
							i.e. last point where interp is >= depth that data*/
	float endDist[MAXOBJECTS] = {0.0};		/*Dist to end of object*/
	int startEl[MAXOBJECTS] = {0.0};		/*Array element containing start of Object*/
	int endEl[MAXOBJECTS] = {0.0};	   		/*Array element containing end of Object*/
	int numberpoints[MAXOBJECTS] = {0.0};		/*Number of data points in object (inc 1st and last)*/
	/*float questionA[MAXOBJECTS] = {0.0};*/
/*numbers required by functions*/
	int numberlines;				/*rows of input data*/
	int totalObjectnumber;

/*As the files are opened for append later, must remove them now*/
	cleanfiles(argc, argv);

/*Reading the data into arrays*/
	numberlines = (ReadData(argv, distance, depth, filtered, longitude, latitude));
/*NOTE: data in arrays elements 0 to numberlines - 1, where numberlies is the rows of data*/
	/*printarrayFFFFF(distance, depth, filtered, longitude, latitude);*/
	/*printf("Number of lines = %d\n", numberlines);*/

/*Finding the seamount objects and their properties*/
	totalObjectnumber = (object(argc, argv, numberlines, distance, depth, filtered, longitude, latitude, startDist, endDist, startEl, endEl, numberpoints));
	/*printf("After Object\n");*/
	
/*Output*/
	/*output(start, end, numberpoints, questionA, questionB, questionC, questionD, questionE, questionF, questionG, questionRA, questionRB, questionRC, questionRD,  argv, totalObjectnumber);*/

}


/*********************************************************************************/

/*FUNCTIONS - ALPHABETICAL*/

/*CLEANFILES*/
/*cleanfiles: because calling files for append later must remove any existing contents*/
int cleanfiles(int argc, char *argv[])

{
	int a, b;				/*loop variables*/
	int Display;				/*Variable to say if producing display output file*/
	char CommandPatterni;			/*variables for the command line grab*/
	char CommandPatternii;
	char *ReturnedPointer;			/*Pointer to string*/
	int foundArg;				/*Says if command line argument been sucessful*/

	/*Display output file*/
	/*Find out if needed and open*/
	/*-DO*/
	CommandPatterni = 'D';
	CommandPatternii = 'O';
	foundArg = 0;
	ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
	if (foundArg == 1)
		{
		printf("Output File %s Cleaned\n", ReturnedPointer);
		remove(ReturnedPointer);
		}

	/*Properties output file*/
	/*Find out if needed and open*/
	/*-PO*/
	CommandPatterni = 'P';
	CommandPatternii = 'O';
	foundArg = 0;
	ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
	if (foundArg == 1)
		{
		printf("Output Properties File %s Cleaned\n", ReturnedPointer);
		remove(ReturnedPointer);
		}

	return 0;
}

/*COMMANDLINE*/
/*Returns a pointer to the elemetnt of the array after the
one that has been found to have the pattern in that I want
to be matched.  The atof function back in the calling function then finds the
next numbers in the string thought to be suitable to be a double*/
char *CommandLine (int argc, char *argv[], char CommandPatterni, char CommandPatternii, int *found)
{
	char *retval;
/*pointer to a char*/
	int i;


	/*printf("Command Pattern %c%c\n", CommandPatterni, CommandPatternii);*/

	*found = 0; 				/*This function gets the pointer to found
						so need to specify that it's the value that we're
						dealing with here*/

	if (argc > 1)
		{
		for (i = 1; i < argc; i++)
			{
			/*printf("argv[being looked at] = %d\n", i);*/
			if (argv[i][0] == '-')
				{
				/*printf("a '-' has been found\n");*/
				if (argv[i][1] == CommandPatterni)
					{
					/*printf("%c has been found\n", CommandPatterni);*/
						if (argv[i][2] == CommandPatternii)
							{
							/*printf("%c has been found\n", CommandPatternii);*/
							retval = &argv[i][3];
							*found = 1;
							}
					}
/*because this is an array the name of the array is the pointer
to the array - Rule must break down when specifying elements of an
array because this will give the value of the element if the & is not
given*/
				}
			}
		}

	if (argc == 1)
		{
		retval = NULL;
/*NULL is (void *) 0long*/
		/*printf("retval is now NULL or %d\n", retval);*/
		}
	if (retval)
		{
		/*printf("retval now contains the value of a pointer = %d\n", retval);*/
		}

	return retval;
}

/*DOMISFITFUNCTION*/
/*DoMisfitFunction: Decides upon wether or not the move is valid
and therefore if it's worth calculating a misfit*/

int DoMisfitFunction(int argc, char *argv[], int f, int g, float MoveLocnX[MVLOC1][MVLOC2], float MoveLocnY[MVLOC1][MVLOC2], int point, float PointCopyX[],float PointCopyY[], float FracDist[], float FracDepth[], int begin, int finish,int objectnumber,float MoveMisfit[MVLOC1][MVLOC2], float MotionUnit)
{
	int todo;				/*Variable indicating if misfit is to be done
						if == 1 then do, if == 0 then don't do*/
	int t;
	float Gradient;				/*Gradient of the baseline from 1st data point*/
	float BaselineHeight;			/*Height of the baseline under point*/

	/*Default is to do*/
	todo = 1;

	/*printf("DoMisfitFunction: MemPosnLocnX = %f MemPosnLocnY = %f\n", MoveLocnX[f][g], MoveLocnY[f][g]);*/

	/*1) check that point is inside box and if not set todo to 0*/
	if(MoveLocnX[f][g] < 0.0 || MoveLocnX[f][g] > 1.0 || MoveLocnY[f][g] < 0.0 || MoveLocnY[f][g] > 1.0)
		{
		todo = 0;
		if (objectnumber == PRINTOBJECT)
				printf("OUTSIDEBOX\n");
		}
	/*printf("todo = %d\n", todo);*/



	/*2) Further restrict points to moving above the linear baseline
	between points 0 and 3*/
	/*Work out the height of the baseline for a given distance across the box*/
	BaselineHeight = FracDepth[begin] + (((MoveLocnX[f][g] - FracDist[begin])/(FracDist[finish] - FracDist[begin])) * (FracDepth[finish] - FracDepth[begin]));
	/*Now see if the point is below this*/
	if(MoveLocnY[f][g] < BaselineHeight)
		{
		todo = 0;
		if (objectnumber == PRINTOBJECT)
				printf("BaselineHeight = %f BELOW BASELINE\n", BaselineHeight);
		}




	/*printf("todo = %d\n", todo);*/

	/*3)Stop the points crossing over*/
	/*top points horizontal done*/
	/*top and bottom points not cross done*/
	if (point == 0)
		{
		/*Vertical with 1*/
		if (MoveLocnY[f][g] > PointCopyY[1])
			{
			todo = 0;
			if (objectnumber == PRINTOBJECT)
				printf("CROSS_OVER stopped\n");
			}
		}
	if (point == 1)
		{
		/*horizontal with 2*/
		if (MoveLocnX[f][g] > PointCopyX[2])
			{
			todo = 0;
			if (objectnumber == PRINTOBJECT)
				printf("CROSS_OVER stopped\n");
			}
		/*Vertical with 0*/
		if (MoveLocnY[f][g] < PointCopyY[0])
			{
			todo = 0;
			if (objectnumber == PRINTOBJECT)
				printf("CROSS_OVER stopped\n");
			}
		}
	if (point == 2)
		{
		/*horizontal with 1*/
		if (MoveLocnX[f][g] < PointCopyX[1])
			{
			todo = 0;
			if (objectnumber == PRINTOBJECT)
				printf("CROSS_OVER stopped\n");
			}
		/*Vertical with 3*/
		if (MoveLocnY[f][g] < PointCopyY[3])
			{
			todo = 0;
			if (objectnumber == PRINTOBJECT)
				printf("CROSS_OVER stopped\n");
			}
		}
	if (point == 3)
		{
		/*Vertical with 2*/
		if (MoveLocnY[f][g] > PointCopyY[2])
			{
			todo = 0;
			if (objectnumber == PRINTOBJECT)
				printf("CROSS_OVER stopped\n");
			}
		}



	/*printf("todo = %d\n", todo);*/

	/*4) Impose constraint that the top must be flat*/
	/*Can do an add in later to allow for a tilt in reality that's smaller
	than a fixed amount and smaller than either of the sides*/
	/*To do this link points 1 + 2 when they're moving up or down*/
	if (point == 1)
		{
		/*If upward move*/
		if (g == 0 || g == 1 || g == 2)
			{
			/*Adjust point 2 as well*/
			PointCopyY[2] = PointCopyY[2] + MotionUnit;
			if (objectnumber == PRINTOBJECT)
				{
				printf("Point1 moved Point2 adjusted up\n");
				}
			}
		/*If downward move*/
		if (g == 4 || g == 5 || g == 6)
			{
			/*Adjust point 2 as well*/
			PointCopyY[2] = PointCopyY[2] - MotionUnit;
			if (objectnumber == PRINTOBJECT)
				printf("Point1 moved Point2 adjusted down\n");
			}
		}
	if (point == 2)
		{
		/*If upward move*/
		if (g == 0 || g == 1 || g == 2)
			{
			/*Adjust point 1 as well*/
			PointCopyY[1] = PointCopyY[1] + MotionUnit;
			if (objectnumber == PRINTOBJECT)
				printf("Point2 moved Point1 adjusted up\n");
			}
		/*If downward move*/
		if (g == 4 || g == 5 || g == 6)
			{
			/*Adjust point 1 as well*/
			PointCopyY[1] = PointCopyY[1] - MotionUnit;
			if (objectnumber == PRINTOBJECT)
				printf("Point2 moved Point1 adjusted down\n");
			}
		}

	/*Summarise the whole affair*/
	if (objectnumber == PRINTOBJECT)
		{
		for (t = 0; t <=3; t++)
			{
			printf("PointCopyX[%d] = %f PointCopyY[%d] = %f\n",t, PointCopyX[t], t, PointCopyY[t]);
			}
		printf("todo = %d\n", todo);
		printf("DO? : MoveMisfit[%d][%d] = %f\n", f, g, MoveMisfit[f][g]);
		}


	return todo;
}


/*Findchunk*/
/*findchunk: determines the end of an object*/
/*returns next time difference is zero even if it's the*/
/*next point*/
int findchunk(float distance[], float depth[], float filtered[], int start, int  numberlines)

{
	int end, search;
	int test;
	int found;

	found = 0;
	end = 0;
/*the first data point to ask questions of should be the one*/
/*after the one the bug is sitting on*/
	search = start + 1;

	/*printf("In findchunk()\n");*/

/*found refers to having found then end of an object*/
/*this happens when the difference between real data and*/
/*filtered data is again zero*/
/*found is set at 0, so not 0 is 1 and the while loop runs*/
	while (!found)
		{
		if ((depth[search] - filtered[search]) <= 0)
			{
			end = search;
			found = 1;
			}

		if (search >= (numberlines - 1))
			{
/*instead of goind too far, it runs out of data*/
/*numberlines -1 is the last address with useful data in*/
                        end = numberlines;
/*do this so that when end is returned the function above knows that*/
/*it's time to break it's loop*/
                        found = 1;
                        }
		++search;
		}

	return end;
}

/*Frustum*/
/*frustum: Do an inversion to approximate the seamounts objects by frustums*/
int frustum(float distance[], float depth[], float FracDist[], float FracDepth[], float filtered[], float longitude[], float latitude[], int begin, int finish,int argc, char *argv[], int objectnumber)

{
	int b,c,d,point,e,f,g,h,i,j,k,l,m;	/*Loop variables*/
	int o,p,q,r,s,t,u,v,w,x;				/*Loop variables*/
	float MinDepth, MaxDepth;		/*Depth to top and bottom of object respectively*/
	int MinDepthLocn, MaxDepthLocn;		/*locations in array of the above*/
/*Location of the 4 points describing the frustum*/
/*points are 0= bottom left, 1= top left, 2 = top right, 3 = bottom right*/
/*x and y are standard cartesian , origin bottom left*/
	float PointX[4] = {-1.0};
	float PointY[4] = {-1.0};
	float PointCopyX[4] = {-1.0};
	float PointCopyY[4] = {-1.0};
/*To do with the misfit of the current configuration*/
	float MISFIT, testMISFIT;		/*Total Misfit*/
/*Amount to move each of the points*/
	float Move[4];
/*Moves around a point being tested*/
/*Three boxes around current best guess at 0 = 1/2 Move, 1 =  Move and 2=  2 x Move away*/
/*0 = top left, 1 =  top centre and so on*/
	float MoveMisfit[MVLOC1][MVLOC2];
	float MoveLocnX[MVLOC1][MVLOC2];
	float MoveLocnY[MVLOC1][MVLOC2];
	float Factor, Multiplier;		/*Factor by which expansion and contraction occurs
						and the multiplier it gives in the loops*/
	float MinMisfit;			/*Value of the minimum misfit*/
	int MinCont, MinNo;			/*array locations of the minimum misfit: MinCont = which contraction state
						MinNo = which number i.e. which direction*/
	int changedMisfit;			/*=0 as default, =1 if misfit better from moving a point*/
	int StillMoving;			/*indicates if all the points have their Move values reduced below MinMove*/
	float MinMove;				/*A level of movement which to give up on!*/
	int MaxItt;				/*Maximum number of itterations*/
	float MinMisfitStop;			/*Value of misfit to stop at because it's sufficiently good*/
/*Resample at equal sapcing along the lines - cf Book 8 pp53*/
	float FracDistRS[MAXTEMP];		/*Fractional Distance and depth of resampled values*/
	float FracDepthRS[MAXTEMP];
	float NumRS;				/*Number of resampled points to used*/
	int beginRS, finishRS;			/*Start and end places on the resampled arrays*/
	int DoMisfit;				/*variable that says wheter or not to bother calculating the misfit*/
	float MotionUnit;
	float ResampInc;			/*Increment for resampling*/



	if (objectnumber == PRINTOBJECT)
		printf("In Frustum\n");
/*Get the locations of the minimum and maximum depths*/

/*Convert depths and distances to fractional distances within object*/
/*Get min(shallow) and max(deep) depths - remember that -ve is still down*/
	MinDepth = -100000.0;
	MaxDepth = 100000.0;
	MinDepthLocn = -1;	/*cos don't have -ve array elements*/
	MaxDepthLocn = -1;
	for(b = begin; b <= finish; b++)
		{
		if(depth[b] < MaxDepth)
			{
			MaxDepth = depth[b];
			MaxDepthLocn = b;
			}
		if(depth[b] > MinDepth)
			{
			MinDepth = depth[b];
			MinDepthLocn = b;
			}
		}
		
	printf("Min depth %f max depth %f\n",MinDepth, MaxDepth);

/*Convert depths and distances to fractional distances within object*/
	normalise(argc, argv, distance, depth, FracDist, FracDepth, begin, finish, MaxDepth, MinDepth, objectnumber);


/*Resample the data so that the line not the data is fitted*/
	NumRS = 100.0;
	ResampInc = resample(argc, argv, FracDist, FracDepth, FracDistRS, FracDepthRS, NumRS, begin, finish);
	/*printf("ResampInc = %f\n", ResampInc);*/



/*Set the intial values of the 4 Points*/
/*Point 0*/
	/*Old Seed*/
	/*PointX[0] = 0.0;
	PointY[0] = 0.0;*/
	PointX[0] = FracDist[begin];
	PointY[0] = FracDepth[begin];
/*Point 1*/
	/*Old Seed*/
	/*PointX[1] = 0.0;*/
	/*seed new v1*/
	/*PointX[1] = 0.25;
	PointY[1] = 1.0;*/
	/*seed new v2*/
	PointX[1] = 0.0;
	PointY[1] = 0.75;
/*Point 2*/
	/*Old Seed*/
	/*PointX[2] = 1.0;*/
	/*seed new v1*/
	/*PointX[2] = 0.75;
	PointY[2] = 1.0;*/
	/*seed new v2*/
	PointX[2] = 1.0;
	PointY[2] = 0.75;
/*Point 3*/
	/*Old Seed*/
	/*PointX[3] = 1.0;
	PointY[3] = 0.0;*/
	PointX[3] = FracDist[finish];
	PointY[3] = FracDepth[finish];
/*if depth of begin or finish is > 0.75 then need to alter seed of [1] and [2]*/
	if (FracDepth[begin] >= 0.75 || FracDepth[finish] >= 0.75)
		{
		if (FracDepth[begin] >= FracDepth[finish])
			{
			PointY[1] = FracDepth[begin] + 0.01;
			PointY[2] = FracDepth[begin] + 0.01;
			}
		if (FracDepth[begin] < FracDepth[finish])
			{
			PointY[1] = FracDepth[finish] + 0.01;
			PointY[2] = FracDepth[finish] + 0.01;
			}
		/*The one proviso to this being if FracDepth[begin] == 1.0
		then SEED at top and move FracDepth[begin] down*/
		if (FracDepth[begin] > 0.99)
			{
			PointY[0] = 0.99;
			PointY[1] = 1.0;
			PointY[2] = 1.0;
			}
		if (FracDepth[finish] > 0.99)
			{
			PointY[1] = 1.0;
			PointY[2] = 1.0;
			PointY[2] = 0.99;
			}
		if (objectnumber == PRINTOBJECT)
			{
			printf("Forced to adjust the seed postition\n");
			printf("PointY[1] = %f PointY[2] = %f\n", PointY[1], PointY[2]);
			}
		}


/*Set an initial value of misfit - verification is being computed*/
	MISFIT = 1000000.0;

/*Work out the misfit of this set of points*/
	beginRS = 0;
	finishRS = NumRS - 1;
	/*if (objectnumber == PRINTOBJECT)
		printf("In Frustum before misfit: begin = %d finish = %d\n", begin, finish);*/
	MISFIT = MisfitFunctionII(PointX, PointY, FracDistRS, FracDepthRS, beginRS, finishRS, objectnumber, ResampInc);
	if (objectnumber == PRINTOBJECT)
		{
		printf("0 = (%f,%f)  1= (%f,%f) 2= (%f,%f) 3= (%f,%f)\n", PointX[0], PointY[0], PointX[1], PointY[1], PointX[2], PointY[2], PointX[3], PointY[3]);
		printf("Initial value of misfit = %f\n", MISFIT);
		}

/*Now I need to find the best fitting set of points*/
	/*Manipulate the points up to a set number of times*/
	/*Just once through to start with*/
	Move[0] = 0.1;
	Move[1] = 0.1;
	Move[2] = 0.1;
	Move[3] = 0.1;

/*do a number of itterations of the algorithm*/
/*Also pull out if MISFIT gets satisfactorily small*/
	MinMove = 0.001;
	MaxItt = 20;
	MinMisfitStop = 0.01;
	StillMoving = 1;
	for (d = 0; d <= MaxItt && MISFIT > MinMisfitStop && StillMoving == 1; d++)
		{
		if (objectnumber == PRINTOBJECT)
			printf("Itteration = %d\n", d);
		/*Loop through the 4 points*/
		for (point = 0; point <= 3; point++)
			{
			if (objectnumber == PRINTOBJECT)
				{
				printf("\nPoint %d\n", point);
				printf("Move[%d] = %f\n", point, Move[point]);
				for (t = 0; t <=3; t++)
					{
					printf("PointX[%d] = %f PointY[%d] = %f\n",t, PointX[t], t, PointY[t]);
					}
				testMISFIT = MisfitFunctionII(PointX, PointY, FracDistRS, FracDepthRS, beginRS, finishRS, objectnumber,ResampInc);
				printf("testMISFIT = %f\n\n", testMISFIT);
				}
			/*for each point remember to reset the values of MoveMisfit*/
			for(u = 0; u <= 2; u++)
				{
				for(v = 0; v <= 7; v++)
					{
					MoveMisfit[u][v] = -1.0;
					}
				}
			/*for each point need to reset the rest of the corners to
			the current best guess values*/
			for (c = 0; c <=3; c++)
				{
				PointCopyX[c] = PointX[c];
				PointCopyY[c] = PointY[c];
				}
			/*For a give point, work out the locations that it can move to*/
			/*cf book 8 pp 52*/
			/*Three loops Contraction [0], then Normal [1] then Expansion [2]*/
			Factor = 2.0;
			for(e = 0; e <= 2; e++)
				{
				if (e == 0)
					Multiplier = (1.0/Factor);
				if (e == 1)
					Multiplier = (1.0);
				if (e == 2)
					Multiplier = (1.0*Factor);
				/*X dirn*/
				MoveLocnX[e][0] = PointCopyX[point] - (Move[point]*Multiplier);
				MoveLocnX[e][1] = PointCopyX[point];
				MoveLocnX[e][2] = PointCopyX[point] + (Move[point]*Multiplier);
				MoveLocnX[e][3] = PointCopyX[point] + (Move[point]*Multiplier);
				MoveLocnX[e][4] = PointCopyX[point] + (Move[point]*Multiplier);
				MoveLocnX[e][5] = PointCopyX[point];
				MoveLocnX[e][6] = PointCopyX[point] - (Move[point]*Multiplier);
				MoveLocnX[e][7] = PointCopyX[point] - (Move[point]*Multiplier);
				/*Y dirn*/
				MoveLocnY[e][0] = PointCopyY[point] + (Move[point]*Multiplier);
				MoveLocnY[e][1] = PointCopyY[point] + (Move[point]*Multiplier);
				MoveLocnY[e][2] = PointCopyY[point] + (Move[point]*Multiplier);
				MoveLocnY[e][3] = PointCopyY[point];
				MoveLocnY[e][4] = PointCopyY[point] - (Move[point]*Multiplier);
				MoveLocnY[e][5] = PointCopyY[point] - (Move[point]*Multiplier);
				MoveLocnY[e][6] = PointCopyY[point] - (Move[point]*Multiplier);
				MoveLocnY[e][7] = PointCopyY[point];
				}
			/*Now run through the different points to get the misfit*/
			/*Work through contraction, normal, and expansion*/
			for(f =0; f <= 2; f++)
				{
				/*Work through points at this range away*/
				for(g = 0; g <= 7; g++)
					{
					/*In case any of my criteria have chaged locns of frustum points that
					are not the immediate concern here i.e. [point]*/
					/*So, anything I do to PointCopy, only applies for this move of this point*/
					for (s = 0; s <= 3; s++)
						{
						PointCopyX[s] = PointX[s];
						PointCopyY[s] = PointY[s];
						}

					/*A unit of vertical motion for this level of contraction is*/
					MotionUnit = (float)fabs(PointCopyY[point] - MoveLocnY[f][g]);
					if (objectnumber == PRINTOBJECT)
						{
						printf("MotionUnit = %f\n", MotionUnit);
						}

					/*Make these equate purely for the point of gettint stuff into MisfitFunctionII*/
					PointCopyX[point] = MoveLocnX[f][g];
					PointCopyY[point] = MoveLocnY[f][g];

					/*Decide if I actuall need to do the misfit*/
					/*Default is yes i.e DoMisfit == 1 until given a reason not to
					when DoMisfit == 0*/

					DoMisfit = DoMisfitFunction(argc, argv, f, g, MoveLocnX, MoveLocnY, point, PointCopyX, PointCopyY, FracDist, FracDepth, begin, finish, objectnumber, MoveMisfit, MotionUnit);

					if (objectnumber == PRINTOBJECT)
						{
						for (w = 0; w <=3; w++)
								{
								printf("PointCopyX[%d] = %f PointCopyY[%d] = %f\n",w, PointCopyX[w], w, PointCopyY[w]);
								}
						printf("MoveMisfit[%d][%d] = %f\n", f, g, MoveMisfit[f][g]);
						}

					/*Actually do the misfit*/
					if (DoMisfit)
						{
						MoveMisfit[f][g] = MisfitFunctionII(PointCopyX, PointCopyY, FracDistRS, FracDepthRS, beginRS, finishRS, objectnumber,ResampInc);
						}
					else
						{
						/*If not doing misfit because e.g. move is outside box then set
						the misfits to -1 (i.e. impossible value)*/
						MoveMisfit[f][g] = -1.0;
						}

					if (objectnumber == PRINTOBJECT)
						{
						for (x = 0; x <=3; x++)
								{
								printf("PointCopyX[%d] = %f PointCopyY[%d] = %f\n",x, PointCopyX[x], x, PointCopyY[x]);
								}
						printf("MoveMisfit[%d][%d] = %f\n\n", f, g, MoveMisfit[f][g]);
						}

					/*printf("MoveLocnX[%d][%d] = %f\tMoveLocnY[%d][%d] = %f\n", f, g, MoveLocnX[f][g], f, g, MoveLocnY[f][g]);
					printf("0 = (%f,%f)  1= (%f,%f) 2= (%f,%f) 3= (%f,%f)\n", PointCopyX[0], PointCopyY[0], PointCopyX[1], PointCopyY[1], PointCopyX[2], PointCopyY[2], PointCopyX[3], PointCopyY[3]);*/
					/*printf("MoveMisfit[%d][%d] = %f\n", f, g, MoveMisfit[f][g]);*/
					}
				}

			/*Print out the moves that a computation has been done for*/
			for(l = 0; l <= 2; l++)
				{
				for(m = 0; m <= 7; m++)
					{
					if(MoveMisfit[l][m] >= 0)
						{
						if (objectnumber == PRINTOBJECT)
							{
							printf("MoveMisfit[%d][%d] = %f\n", l, m, MoveMisfit[l][m]);
							}
						}
					}
				}


			/*Now select the minimum misfit (that's > 0 !) from the values*/
			/*Set min to a ridiculously high value*/
			MinMisfit = 1000000.0;
			for(j = 0; j <= 2; j++)
				{
				for(k = 0; k <= 7; k++)
					{
					if (MoveMisfit[j][k] < MinMisfit && MoveMisfit[j][k] >= 0)
						{
						MinMisfit = MoveMisfit[j][k];
						MinCont = j;
						MinNo = k;
						}
					}
				}
			if (objectnumber == PRINTOBJECT)
				printf("Minimum misfit = %f  cont = %d No. = %d\n", MinMisfit, MinCont, MinNo);

			/*Now Change the locations of the best guess points to include the minimum misfit found here*/
			changedMisfit = 0;
			if (MinMisfit < MISFIT)
				{
				if (objectnumber == PRINTOBJECT)
					printf("Point %d location changed from (%f,%f) to (%f,%f)\n", point, PointX[point], PointY[point], MoveLocnX[MinCont][MinNo], MoveLocnY[MinCont][MinNo]);
				PointX[point] = MoveLocnX[MinCont][MinNo];
				PointY[point] = MoveLocnY[MinCont][MinNo];
				/*Also adjust the second point in order to keep the top flat*/
				/*For now can assume that heights will be the same as going for pure flat*/
				/*Take the case of point == 1 first*/
				if (point == 1 || point == 2)
					{
					PointY[1] = MoveLocnY[MinCont][MinNo];
					PointY[2] = MoveLocnY[MinCont][MinNo];
					}
				if (objectnumber == PRINTOBJECT)
					printf("MISFIT reduced to %f\n",MinMisfit);
				MISFIT = MinMisfit;
				changedMisfit = 1;
				}
			else
				{
				if (objectnumber == PRINTOBJECT)
					printf("Existing point was at least equal best fit, no change\n");
				}

			/*Now set Move[] to be appropriate i.e search at any given point expands or contracts*/
			if (MinCont == 0 && changedMisfit == 1)
				{
				Move[point] =  Move[point]/Factor;
				if (objectnumber == PRINTOBJECT)
					printf("Move[%d] changed to %f\n", point, Move[point]);
				}
			if (MinCont == 2 && changedMisfit == 1)
				{
				Move[point] =  Move[point]*Factor;
				if (objectnumber == PRINTOBJECT)
					printf("Move[%d] changed to %f\n", point, Move[point]);
				}
			/*Contract if there has been no change*/
			if (changedMisfit == 0)
				{
				Move[point] =  Move[point]/Factor;
				if (objectnumber == PRINTOBJECT)
					printf("But, Move[%d] shrunk anyway to %f\n", point, Move[point]);
				}

			/*Report the reason for pulling out of the loop*/
			if (Move[0] < MinMove && Move[1] < MinMove && Move[2] < MinMove && Move[3] < MinMove)
				{
				StillMoving = 0;
				if (objectnumber == PRINTOBJECT)
					printf("Stopped because points stopped moving\n");
				}
			if (d == MaxItt)
				{
				if (objectnumber == PRINTOBJECT)
					printf("Stopped because of Itteration limit\n");
				}
			if (MISFIT <= MinMisfitStop)
				{
				if (objectnumber == PRINTOBJECT)
					printf("Stopped because fit is sufficiently good\n");
				}
			}
		}
		/*Output Various things to various files*/
		output(distance, depth, longitude, latitude, begin, finish, PointX, PointY, MinDepth, MaxDepth, argc, argv, objectnumber, MinDepthLocn, MISFIT, FracDist, FracDepth);


	return 0;
}

/*MISFITFUNCTION version ii*/
/*MisfitFunctionII: Works out a misfit for a set of coordinates given to it*/
/*in the new version last four inputs are followed by "RS" in the 'frustum' function*/
/*This version is a differnt fit algorithm based on an idea in Book8 pp56*/
/*NOTE: This is probably best seeded by points 1 & 2 being half way up on edges*/
float MisfitFunctionII(float PointX[], float PointY[], float FracDist[], float FracDepth[], int begin, int finish, int objectnumber, float ResampInc)

{
	int a,b,c,d,e;				/*loop variable*/
	float misfit[4];			/*Misfit for eash of the three lines left, centre, right*/
	float count[4];
	float Avmisfit[4];
	float maxHeight;			/*position of the highest point in the feature*/
	int maxHeightLocn;
	float Idist;				/*distance of line at height of data point*/
	float Idepth;				/*depth of line at height of data point*/
	float misfitTemp;			/*Misfit for each data point*/
	float SumMisfit;			/*Misfit for the section with the input frustum points*/
	float BaselineFracDist;			/*Baseline resampled at the same increments*/
	float BaselineFracDepth;		/*as the data line*/
	float BaseIncX, BaseIncY;		/*Increments in X and Y for baseline resampling*/
	float theta;				/*Angle of the baseline*/
	float CountBase;			/*number of points along baseline*/
	float betweenpercent, test, acrossA, acrossB;


	/*printf("0 = (%f,%f)  1= (%f,%f) 2= (%f,%f) 3= (%f,%f)\n", PointX[0], PointY[0], PointX[1], PointY[1], PointX[2], PointY[2], PointX[3], PointY[3]);*/

/*In this Line 1 is fit horizontally to all data before horizontal distance of maxheight that
are lower in height than point 1.  Line 3 is mirror. Line 2 is fit vertically to the remainder
or at the very least the point of maximum height*/

/*Set the initial values for the three sections i.e. 3 lines*/
/*Another element added to the array for fit of the baseline*/
	for (c = 0; c <= 3; c++)
		{
		misfit[c] = 0.0;
		count[c] = 0.0;
		Avmisfit[c] = 0.0;
		}

/*Do all (now first) 3 misfits at once - deciding point by point which on it contributes to*/
	/*Firstly need to know the location of the maxheight*/
	/*initialise values*/
	maxHeight = 0.0;
	maxHeightLocn = -1;
	for(a = begin; a <= finish; a++)
		{
		if(FracDepth[a] > maxHeight)
			{
			maxHeight = FracDepth[a];
			maxHeightLocn = a;
			}
		}
	/*The actual misfit loop*/
	if (objectnumber == PRINTOBJECT)
		{
		printf("begin = %d finish = %d\n", begin, finish);
		for(e = 0; e <= 3; e++)
			{
			printf("PointX[%d] = %f PointY[%d] = %f\n", e, PointX[e], e, PointY[e]);
			}
		}
	for (b = begin; b <= finish; b++)
		{
		if(b <= maxHeightLocn && FracDepth[b] < PointY[1])
			{
			/*contributes to the first line*/
			/*distance of line at height of data point*/
			/*Book 8 pp58*/
			Idist = PointX[0] + (((FracDepth[b]-PointY[0])/(PointY[1] - PointY[0]))*(PointX[1] - PointX[0]));
			misfitTemp = (float)fabs(Idist - FracDist[b]);
			count[0] = count[0] + 1.0;
			misfit[0] = misfit[0] + misfitTemp;
			/*if (objectnumber == PRINTOBJECT)
				printf("Point %d (%f,%f) counted to Line 1\n", b, FracDist[b], FracDepth[b]);*/
			}

		else if(b >= maxHeightLocn && FracDepth[b] < PointY[2])
			{
			/*contributes to the third line*/
			/*Book 8 pp58*/
			Idist = PointX[2] + (((FracDepth[b]-PointY[2])/(PointY[3] - PointY[2]))*(PointX[3] - PointX[2]));
			misfitTemp = (float)fabs(Idist - FracDist[b]);
			count[2] = count[2] + 1.0;
			misfit[2] = misfit[2] + misfitTemp;
			/*if (objectnumber == PRINTOBJECT)
				printf("Point %d (%f,%f) counted to Line 3\n", b, FracDist[b], FracDepth[b]);*/
			}

		else if((b <= maxHeightLocn && FracDepth[b] >= PointY[1]) || (b >= maxHeightLocn && FracDepth[b] >= PointY[2]))
			{
			/*printf("Point %d (%f,%f) not in Line 1 or Line 3\n", b, FracDist[b], FracDepth[b]);*/
			/*contributes to the second line*/
			/*Book 8 pp58*/
			Idepth = PointY[1] + (((FracDist[b] - PointX[1])/(PointX[2] - PointX[1]))*(PointY[2] - PointY[1]));
			misfitTemp = (float)fabs(Idepth - FracDepth[b]);
			count[1] = count[1] + 1.0;
			misfit[1] = misfit[1] + misfitTemp;
			/*if (objectnumber == PRINTOBJECT)
				printf("Point %d (%f,%f) counted to Line 2\n", b, FracDist[b], FracDepth[b]);*/
			}

		else
			{
			if (objectnumber == PRINTOBJECT)
				{
				printf("WARNING: Fitting logic has screwed up\n");
				printf("Point %d (%f,%f)\n", b, FracDist[b], FracDepth[b]);
				printf("Point 1  (%f,%f) and point 2 (%f,%f)\n", PointX[1], PointY[1], PointX[2], PointY[2]);
				}
			}
		}
	/*THINGS TO MAKE SURE THAT MISFITS ALWAYS HAVE A VALUE*/
	/*Make sure that line 2 always counts the maxheight point*/
	if (count[1] == 0)
		{
		/*Add the caluclation*/
		Idepth = PointY[1] + (((FracDist[maxHeightLocn] - PointX[1])/(PointX[2] - PointX[1]))*(PointY[2] - PointY[1]));
		misfitTemp = (float)fabs(Idepth - FracDepth[maxHeightLocn]);
		count[1] = count[1] + 1.0;
		misfit[1] = misfit[1] + misfitTemp;
		if (objectnumber == PRINTOBJECT)
			printf("Forcing Top line to include the Maxheight point\n");
		}
	/*Make sure that lines 1 and 3 always count the first and last point respectively*/
	if (count[0] == 0)
		{
		/*Add the caluclation*/
		Idist = PointX[0] + (((FracDepth[begin]-PointY[0])/(PointY[1] - PointY[0]))*(PointX[1] - PointX[0]));
		misfitTemp = (float)fabs(Idist - FracDist[begin]);
		count[0] = count[0] + 1.0;
		misfit[0] = misfit[0] + misfitTemp;
		if (objectnumber == PRINTOBJECT)
			printf("Forcing Left line to include the First point\n");
		}
	if (count[2] == 0)
		{
		/*Add the caluclation*/
		/*test = 0.0/2.7;
		printf("test = %f\n", test);*/
		acrossA = (FracDepth[finish] - PointY[2]);
		acrossB = (PointY[3] - PointY[2]);
		betweenpercent = (acrossA/acrossB);
		Idist = PointX[2] + (betweenpercent*(PointX[3] - PointX[2]));
		misfitTemp = (float)fabs(Idist - FracDist[finish]);
		count[2] = count[2] + 1.0;
		misfit[2] = misfit[2] + misfitTemp;
		if (objectnumber == PRINTOBJECT)
			{
			printf("acrossA = %f  acrossB = %f\n", acrossA, acrossB);
			printf("FracDepth[%d] = %f betweenpercent = %f \n", finish, FracDepth[finish], betweenpercent);
			printf("Forcing Right line to include the Last point\n");
			printf("Idist = %f misfitTemp = %f\n",Idist, misfitTemp);
			}
		}

/*Now do the extra misfit for the baseline*/
/*Resampinc is distance along line to resample at*/
	/*note for theta across is 1 unit*/
	theta = (float)atan(FracDepth[finish] - FracDepth[begin]);
	BaseIncX = ((float)cos(theta))*ResampInc;
	BaseIncY = ((float)sin(theta))*ResampInc;
	if (objectnumber == PRINTOBJECT)
		printf("ReampInc = %f, Xinc = %f, Yinc = %f\n", ResampInc, BaseIncX, BaseIncY);
	/*Set the first point*/
	BaselineFracDist = FracDist[begin];
	BaselineFracDepth = FracDepth[begin];
	for(count[3] = 0; BaselineFracDist <= 1.0; count[3]++)
		{
		/*The calculation*/
		/*see book 8 pp65*/
		Idepth = PointY[0] + (((BaselineFracDist - PointX[0])/(PointX[3] - PointX[0]))*(PointY[3] - PointY[0]));
		misfitTemp = (float)fabs(Idepth - BaselineFracDepth);
		count[3] = count[3] + 1.0;
		misfit[3] = misfit[3] + misfitTemp;
		/*Moving things on for the next calculation*/
		BaselineFracDist = BaselineFracDist + BaseIncX;
		BaselineFracDepth = BaselineFracDepth + BaseIncY;
		}

/*Get average misfits*/
	for (d = 0; d <= 3; d++)
		{
		Avmisfit[d] = misfit[d]/count[d];
		if (objectnumber == PRINTOBJECT)
			printf("count[%d] = %f, misfit[%d] =%f, Avmisfit[%d] = %f\n", d, count[d],d ,misfit[d],d, Avmisfit[d]);
		}

/*Sum the misfitts*/
	SumMisfit = Avmisfit[0] + Avmisfit[1] + Avmisfit[2] + Avmisfit[3];
		if (objectnumber == PRINTOBJECT)
			{
			printf("Avmisfit[0] = %f Avmisfit[1] = %f Avmisfit[2] = %f Avmisfit[3] = %f\n", Avmisfit[0],Avmisfit[1],Avmisfit[2],Avmisfit[3]);
			}

	/*Return the misfit for a given set of points*/
	return SumMisfit;
}

/*NORMALIZE*/
/*Normalize the distances and depths to fractions of the objects dimensions*/
int normalise(int argc, char *argv[], float distance[], float depth[], float FracDist[], float FracDepth[], int begin, int finish, float MaxDepth, float MinDepth, int objectnumber)

{
	int a;				/*Loop variables*/

	for(a = begin; a <= finish; a++)
		{
		FracDist[a] = (distance[a] - distance[begin])/(distance[finish] - distance[begin]);
		FracDepth[a] = (depth[a] - MaxDepth)/(MinDepth - MaxDepth);
		if (objectnumber == PRINTOBJECT)
			printf("begin %d finish %d, distance[%d] = %f depth[%d] = %f FracDist[%d] = %f FracDepth[%d] = %f\n", begin, finish, a, distance[a], a, depth[a], a, FracDist[a], a, FracDepth[a]);
		}

	return 0;
}


/*Object*/
/*object: anything with a point inbetween some zero's*/
/*is an object.  It has properties to be found out by questions*/
int object (int argc, char *argv[], int numberlines, float distance[], float depth[], float filtered[], float longitude[], float latitude[], float startDist[], float endDist[], int startEl[], int endEl[], int numberpoints[])

{
	int begin, finish;			/*variables to find beginning and end of object*/
						/*These are no. of array elements - start, end are distances*/
	int objectnumber;			/*number of the object*/
	int pointcount;				/*variable to find number of data points in object*/
	float FracDist[MAXCOLUMNS] = {0.0};	/*Fraction Dist along object*/
	float FracDepth[MAXCOLUMNS] = {0.0};	/*Fraction of height of object*/
						/*NOTE: these are set on the same array size as the input data
						for ease of linking back to actual values.  Points inside objects
						overwritten when frustum function is called*/

	begin = 0;			/*Start at element zero of data arrays*/
	objectnumber = 1;		/*Start with 1st object*/

	printf("In Object %d %d\n", begin, numberlines);

	while (begin < numberlines)
/*need to have < because numberlines is e.g. 4330 and arrays go to 4329*/
/*also when findchunk returns numberlines, it means it has*/
/*run out of data*/
	{
	/*printf("Before findchunk begin is %d finish is %d \n", begin, finish);*/
	finish = (findchunk(distance, depth, filtered, begin, numberlines));
	/*printf("After findchunk begin is %d finish is %d \n", begin, finish);*/
	 /*passes out an finish, begin stays the same*/

	if (finish < numberlines)
		{
		pointcount = ((finish - begin) + 1);		/* +1 is the fence-post factor*/
		if (pointcount > 2)
/*It's not an object if it's just 2 points*/
/*I don't want the line counted as an object*/
			{

/*Now the aim is to approximate it as a frustum*/
/*This is still OK with only 3 points, or should be, will be a good test non the less*/
/*NOTE: I'm doing this simplistically with no regional tilt correction to the data*/
			printf("Before Frustum\n");
			frustum(distance, depth, FracDist, FracDepth, filtered, longitude, latitude, begin, finish, argc, argv, objectnumber);

/*make sure all facts about object are in their relevant arrays*/
			startDist[objectnumber] = distance[begin];
			endDist[objectnumber] = distance[finish];
			startEl[objectnumber] = begin;
			endEl[objectnumber] = finish;
			numberpoints[objectnumber] = pointcount;
			if (objectnumber == PRINTOBJECT)
				{
				printf("objectnumber %d\nstart [%f] (begin = %d)\nend [%f] (finish = %d)\n numberpoints[%d]\n\n", objectnumber, startDist[objectnumber], startEl[objectnumber], endDist[objectnumber], endEl[objectnumber], numberpoints[objectnumber]);
				}
			objectnumber++;
			}
		}
	begin = finish;
/*needs to be outside the pointcount if as it needs to*/
/*increment even if they're not counted as an object*/
	}


	return (objectnumber - 1);
}

/*OUTPUT*/
/*output: Will put out various frustum data (and properties calculated from this) to files*/
/*I know that calling once for each output object is not the most efficient way to do this
but that shouldn't matter too much*/
int output(float distance[],float depth[],float longitude[], float latitude[], int begin,int finish,float PointX[],float  PointY[], float MinDepth, float MaxDepth,int argc, char *argv[], int objectnumber, int MinDepthLocn, float MISFIT, float FracDist[], float FracDepth[])

{
	int a, b, c, d, e, f,g,h;		/*loop variables*/
	int i,j;				/*loop variable*/
	float PointDist[4], PointDepth[4];	/*frustum conrners in the real world of d,t*/
	char CommandPatterni;			/*variables for the command line grab*/
	char CommandPatternii;
	char *ReturnedPointer;			/*Pointer to string*/
	int foundArg;				/*Says if command line argument been sucessful*/
	float AvRadiusTop, AvRadiusBottom;	/*Dist (km) from the centre to the edges of the top
						and bottom of the frustum*/
	float AvHeight;				/*Av height (m) between the top and bottom 2 corners of the frustum*/
	float AvSlope;				/*Av slope of the sides of the frustum*/
	float Tilt;				/*Tilt of the top of the frustum (degrees)*/
	float BaseRadius, TopRadius;		/*as it says*/
	float SlopeRHS, SlopeLHS, Slope;	/*Slopes of the sides of the approximated cone*/
	float up,across, upB, acrossB;
	float topDepth, bottomDepth, PedestalHeight;
	float PointLatitude[4], PointLongitude[4];	/*long and lat of the points*/
	float MeasuredArea;
	float FrustumArea;
	float TrapeziumAvDist, TrapBaseDepth, TrapTopDepth;	/*to work out measured area*/
	float ExtraArea;
	float percentbetween;
	float PointyConeHeight, WholeConeVol, TopBitVol, FrustumVol;
	float LatKm, LatRadians, LatRadius;			/*For the non-linearity thing*/
	float LongKm, maxheight, maxheightLong, maxheightLat;
	float dista, distb, distc;
	float intermediate, Angle;
	float DistAlong, DistAlongPart, ExtraTravelFrac;

	FILE *fdisplay;
	FILE *fproperties;

	/************************************************/
	/*Don't make any of the computations conditional*/
	/************************************************/

	/*Get the frustum parameters into the real world*/
		/*NOTE: MinDepth is depth at top of object i.e. shallowest point, depths -ve*/
		for(a = 0; a <= 3; a ++)
			{
			PointDist[a] = distance[begin] + ((distance[finish] - distance[begin])*PointX[a]);
			PointDepth[a] = MaxDepth + ((MinDepth-MaxDepth)*PointY[a]);
			}

	TopRadius = -1.0;
	BaseRadius = -1.0;
	BaseRadius = (PointDist[3] - PointDist[0])/2.0;
	TopRadius = (PointDist[2] - PointDist[1])/2.0;
	if(PointDist[2] == PointDist[1])
		{
		TopRadius = 0.0;
		}
	/*if (objectnumber == PRINTOBJECT)
		{
		for(h = 0; h <= 3; h++)
			{
			printf("PointDist[%d] = %f\n",h, PointDist[h]);
			}
		printf("BaseRadius = %f TopRadius= %f\n",BaseRadius, TopRadius);
		}*/

	/*Make sure that both of these are positive*/
	/*Note distance is in km, depth in m*/
	up = (PointDepth[1] - PointDepth[0]);
	across = (PointDist[1] - PointDist[0]);
	SlopeLHS = (float)atan(up/across);
	if (objectnumber == PRINTOBJECT)
		{
		printf("up = %f across = %f\n", up, across);
		printf("\nSlope LHS%.3f\t%.3f\t%.3f\n", up, across, SlopeLHS);
		}
	upB = (PointDepth[2] - PointDepth[3]);
	acrossB = ((PointDist[3] - PointDist[2]));
	SlopeRHS = (float)atan(upB/acrossB);
	/*And deal with the cases where points 0==1 and 2==3*/
	if (up < 0.001 || across < 0.001)
		{
		SlopeLHS = SlopeRHS;
		}
	if (upB < 0.001 || acrossB < 0.001)
		{
		SlopeRHS = SlopeLHS;
		}
	if (objectnumber == PRINTOBJECT)
		{
		printf("upB = %f acrossB = %f\n", upB, acrossB);
		printf("Slope RHS %.3f\t%.3f\t%.3f\n", upB, acrossB, SlopeRHS);
		}
	Slope = ((SlopeRHS + SlopeLHS)/2.0)*360.0/(2.0*PI);

	/*make it general so that I won't forget when adding tilt*/
	topDepth = (PointDepth[1] + PointDepth[2])/2.0;
	bottomDepth = (PointDepth[0] + PointDepth[3])/2.0;
	PedestalHeight = topDepth - bottomDepth;

	/*long lat of the points and the centre of the top*/
	/*take these from the data points either side*/
	for(e = 0; e <= 3; e++)
		{
		/*printf("begin = %d, Fracdist[%d] = %f, PointX[%d] = %f\n", begin, d, FracDist[d], e, PointX[e]);*/
		for(d = begin; d <= finish; d ++)
			{
			if (FracDist[d] <= PointX[e] && FracDist[d + 1] >= PointX[e])
				{
				PointLatitude[e] = 0.0;
				PointLongitude[e] = 0.0;
				/*printf("Working out point lats and longs\n");*/
				if (objectnumber == PRINTOBJECT)
					{
					printf("PointX[%d] = %f FracDist[%d] = %f FracDist[%d] = %f\n", e, PointX[e], d, FracDist[d],d +1, FracDist[(d+1)]);
					printf("long [%d] = %f lat[%d] = %f\n", d, longitude[d], d, latitude[d]);
					printf("long [%d] = %f lat[%d] = %f\n", d + 1, longitude[d+ 1], d+ 1, latitude[d+ 1]);
					}
				percentbetween = ((PointX[e] - FracDist[d])/(FracDist[(d+1)] - FracDist[d]));
				PointLatitude[e] = latitude[d] + (percentbetween*(latitude[(d+1)] - latitude[d]));
				PointLongitude[e] = longitude[d] + (percentbetween*(longitude[(d+1)] - longitude[d]));
				if (objectnumber == PRINTOBJECT)
					{
					printf("perceentbtween = %f\n", percentbetween);
					printf("Pointlong = %f PointLat = %f\n", PointLongitude[e], PointLatitude[e]);
					}
				}
			}
		}

	/*Compute the measured area (km^2)*/
	/*Set a few things to zero*/
	MeasuredArea = 0.0;
	TrapeziumAvDist = 0.0;
	percentbetween = 0.0;
	TrapBaseDepth = 0.0;
	TrapTopDepth = 0.0;
	ExtraArea = 0.0;
	FrustumArea = 0.0;
	PointyConeHeight = 0.0;
	WholeConeVol = 0.0;
	TopBitVol = 0.0;
	FrustumVol = 0.0;
	for (f = begin; f < finish; f++)
		{
		TrapeziumAvDist = (distance[f + 1] + distance[f])/2.0;
		percentbetween = ((TrapeziumAvDist - distance[begin])/(distance[finish] - distance[begin]));
		TrapBaseDepth = depth[begin] + (percentbetween*(depth[finish] - depth[begin]));
		TrapTopDepth = (depth[f] + depth[f + 1])/2.0;
		if (objectnumber == PRINTOBJECT)
			{
			printf("distance[f + 1] %f distance[f] %f percentbetween %f\n", distance[f + 1],distance[f], percentbetween);
			printf("depth[f + 1] %f depth[f] %f\n", depth[f + 1],depth[f]);
			printf("TrapeziumAvDist %f TrapBaseDepth %f TrapTopDepth %f\n", TrapeziumAvDist, TrapBaseDepth, TrapTopDepth);
			}
		ExtraArea = ((TrapTopDepth - TrapBaseDepth))*(distance[f + 1] - distance[f]);
		if (objectnumber == PRINTOBJECT)
			{
			printf("ExtraArea %f\n", ExtraArea);
			}
		MeasuredArea = MeasuredArea + ExtraArea;
		}
	if (objectnumber == PRINTOBJECT)
		{
		printf("MeasuredArea = %f\n", MeasuredArea);
		}
	/*Compute the Frustum area (km^2)*/
	FrustumArea = ((PedestalHeight)*TopRadius*2.0) + (0.5*(BaseRadius - TopRadius)*(PedestalHeight)*2.0);
	/*printf("PedestalHeight %f TopRadius %f BaseRadius %f \n", PedestalHeight, TopRadius, BaseRadius);
	printf("FrustumArea %f\n", FrustumArea);*/
	if (objectnumber == PRINTOBJECT)
		{
		printf("FrustumArea %f\n", FrustumArea);
		}
	/*Compute the Frustum volume (km^3)*/
	/*See Book 8 pp69*/
	if (TopRadius > 0.0)
		{
		PointyConeHeight = ((1.0/(1.0 - (TopRadius/BaseRadius)))*PedestalHeight);
		WholeConeVol = ((1.0/3.0)*PI*BaseRadius*BaseRadius*PointyConeHeight);
		TopBitVol = ((1.0/3.0)*PI*TopRadius*TopRadius*(PointyConeHeight - (PedestalHeight)));
		FrustumVol = WholeConeVol - TopBitVol;
		}
	if (TopRadius == 0.0)
		{
		FrustumVol = ((1.0/3.0)*PI*BaseRadius*BaseRadius*(PedestalHeight));
		}
	if (TopRadius == BaseRadius)
		{
		FrustumVol = PI*BaseRadius*BaseRadius*(PedestalHeight);
		}
	if (objectnumber == PRINTOBJECT)
		{
		printf("PedestalHeight %f\n ", PedestalHeight);
		printf("PointyConeHeight = %f WholeConeVol = %f TopBitVol = %f FrustumVol = %f\n", PointyConeHeight, WholeConeVol, TopBitVol, FrustumVol);
		}

	/*Calculate a non-linearity thing*/
	/*Firstly - the angle between begin-maxheight-finish*/
	/*Assume flat earth*/
	LatKm = (6370.0*2.0*PI)/360.0;
	LatRadians = (float)fabs(latitude[begin]*2.0*PI/360.0);
	LatRadius = 6370.0*(float)cos(LatRadians);
	LongKm = PI*LatRadius*2.0/360.0;
	if (objectnumber == PRINTOBJECT)
		{
		printf("latitude[begin] = %f, LatRadius %f, LongKm %f LatKm %f\n", latitude[begin], LatRadius, LongKm, LatKm);
		}
	/*Find the location of the max height*/
	maxheight = -100000.0;
	maxheightLong = -1.0;
	maxheightLat = -1.0;
	/*starts at begin +1 so that the first point can't also be
	the heighest and break things*/
	for(i = begin +1; i <= finish; i++)
		{
		if(depth[i] > maxheight)
			{
			maxheight = depth[i];
			maxheightLong = longitude[i];
			maxheightLat = latitude[i];
			}
		}
	if (objectnumber == PRINTOBJECT)
		{
		printf("maxheightLong = %f, maxheightLat %f\n", maxheightLong, maxheightLat);
		}
	/*Direction of start to maxheight*/
	/*Safest to use the cosine rule Book8pp135*/
	/*a*/
	across = ((longitude[begin]-longitude[finish])*LongKm);
	up = ((latitude[begin]-latitude[finish])*LatKm);
	dista = (float)sqrt((up*up)+(across*across));
	if (objectnumber == PRINTOBJECT)
		{
		printf("dist a ... (%f,%f) - (%f,%f)\n", longitude[begin], latitude[begin], longitude[finish], latitude[finish]);
		printf("up = %f across = %f dista = %f\n", up, across, dista);
		}
	/*b*/
	across = ((longitude[begin]-maxheightLong)*LongKm);
	up = ((latitude[begin]-maxheightLat)*LatKm);
	distb = (float)sqrt((up*up)+(across*across));
	if (objectnumber == PRINTOBJECT)
		{
		printf("dist b ... (%f,%f) - (%f,%f)\n", longitude[begin], latitude[begin], maxheightLong, maxheightLat);
		printf("up = %f across = %f distb = %f\n", up, across, distb);
		}
	/*c*/
	across = ((maxheightLong-longitude[finish])*LongKm);
	up = ((maxheightLat-latitude[finish])*LatKm);
	distc = (float)sqrt((up*up)+(across*across));
	if (objectnumber == PRINTOBJECT)
		{
		printf("dist c ... (%f,%f) - (%f,%f)\n", maxheightLong, maxheightLat, longitude[finish], latitude[finish]);
		printf("up = %f across = %f distc = %f\n", up, across, distc);
		}
	/*a sum*/
	intermediate = ((distb*distb)+(distc*distc)-(dista*dista))/(2.0*distb*distc);
	/*Sometimes numerical error can mean that b + c < a*/
	/*and if this is the case, it's a straight line i.e. intermediate should be -1 exactly*/
	/*this will then hopefully give me an output of PI radians*/
	if(dista > (distb + distc) || intermediate < -1.0)
		{
		intermediate = -1.0;
		}
	/*A*/
	Angle = ((float)acos(intermediate))*360.0/(2.0*PI);
	if (objectnumber == PRINTOBJECT)
		{
		printf("Intermediate %f Angle %f\n", intermediate, Angle);
		}
	/*Now do a difference in the distances*/
	/*dista from above is the direct distance*/
	DistAlong = 0.0;
	if (objectnumber == PRINTOBJECT)
		{
		printf("Working out 3D distance\n");
		}
	for(j = begin; j < finish; j++)
		{
		across = ((longitude[j]-longitude[j+1])*LongKm);
		up = ((latitude[j]-latitude[j+1])*LatKm);
		DistAlongPart = (float)sqrt((up*up)+(across*across));
		/*if (objectnumber == PRINTOBJECT)
			{
			printf("begin = %d (longitude[%d]-longitude[%d]) (%f-%f)\n",begin, j, j +1, longitude[j], longitude[j+1]);
			printf("begin = %d (latitude[%d]-latitude[%d])(%f-%f)\n",begin, j, j +1, latitude[j], latitude[j+1]);
			printf("DistAlongPart = %f\n", DistAlongPart);
			}*/
		DistAlong = DistAlong + DistAlongPart;
		}
	ExtraTravelFrac = DistAlong/dista;
	if (objectnumber == PRINTOBJECT)
		{
		printf("DistAlong %f / dista %f = ExtraTravelFrac %f\n", DistAlong, dista, ExtraTravelFrac);
		}

	/*Display Output file*/
	/*Find out if needed and open*/
	/*-DO*/
	CommandPatterni = 'D';
	CommandPatternii = 'O';
	foundArg = 0;
	ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
	if (foundArg == 1)
		{
		/*Start up the results file*/
		fdisplay = fopen (ReturnedPointer, "a");
		/*printf("Display Output File is %s\n", ReturnedPointer );*/

		/*Get the frustum parameters into the real world*/
		/*NOTE: MinDepth is depth at top of object i.e. shallowest point, depths -ve*/
		/*Output these*/
		for(b = 0; b <= 3; b ++)
			{
			if (b == 0)
				{
				fprintf(fdisplay, ">\n%f %f\n", PointDist[b], PointDepth[b]);
				if (objectnumber == PRINTOBJECT)
					printf(">\n%f %f\n", PointDist[b], PointDepth[b]);
				}
			else
				{
				fprintf(fdisplay, "%f %f\n", PointDist[b], PointDepth[b]);
				if (objectnumber == PRINTOBJECT)
					printf("%f %f\n", PointDist[b], PointDepth[b]);
				}
			}
		/*Close the file again*/
		/*Inside loop as can only be closed if has been opened*/
		fclose(fdisplay);
		}

	/*Properties Output file*/
	/*Find out if needed and open*/
	/*-PO*/
	CommandPatterni = 'P';
	CommandPatternii = 'O';
	foundArg = 0;
	ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
	if (foundArg == 1)
		{
		/*Start up the results file*/
		fproperties = fopen(ReturnedPointer, "a");
		/*printf("Properties Output File is %s\n", ReturnedPointer );*/
		/*Output contents*/
		/*Firstly the objectnumber*/
		CommandPatterni = 'N';
		CommandPatternii = 'o';
		foundArg = 0;
		ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
		if (foundArg == 1)
			{
			fprintf(fdisplay, "%d\t", objectnumber);
			if (objectnumber == PRINTOBJECT)
				printf("%d\t", objectnumber);
			}
		/*Properties of the measured object*/
		CommandPatterni = 'M';
		CommandPatternii = 'a';
		foundArg = 0;
		ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
		if (foundArg == 1)
			{
			/*Distance along track to the shallowest point*/
			fprintf(fdisplay, "%.3f\t", distance[MinDepthLocn]);
			if (objectnumber == PRINTOBJECT)
				printf("%.3f\t", distance[MinDepthLocn]);
			/*Distance object measured to start and stop*/
			fprintf(fdisplay, "%.3f\t%.3f\t", distance[begin],distance[finish]);
			if (objectnumber == PRINTOBJECT)
				printf("%.3f\t%.3f\t", distance[begin],distance[finish]);
			/*Long and lat of the Shallowest point*/
			fprintf(fdisplay, "%.3f\t%.3f\t", longitude[MinDepthLocn],latitude[MinDepthLocn]);
			if (objectnumber == PRINTOBJECT)
				printf("%.3f\t%.3f\t", longitude[MinDepthLocn],latitude[MinDepthLocn]);
			}

		/*Distance and depth of each of the points of the best-fit frustum*/
		CommandPatterni = 'A';
		CommandPatternii = 'd';
		foundArg = 0;
		ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
		if (foundArg == 1)
			{
			for(c = 0; c <= 3; c ++)
				{
				fprintf(fdisplay, "%.3f\t%.3f\t", PointDist[c], PointDepth[c]);
				if (objectnumber == PRINTOBJECT)
					printf("%.3f\t%.3f\t", PointDist[c], PointDepth[c]);
				}
			}
		/*Basal and top radii of the approximated cone*/
		CommandPatterni = 'A';
		CommandPatternii = 'r';
		foundArg = 0;
		ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
		if (foundArg == 1)
			{
			fprintf(fdisplay, "%.3f\t%.3f\t", BaseRadius, TopRadius);
			if (objectnumber == PRINTOBJECT)
				printf("%.3f\t%.3f\t", BaseRadius, TopRadius);
			}
		/*Average slope of the approximated cone*/
		CommandPatterni = 'A';
		CommandPatternii = 's';
		foundArg = 0;
		ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
		if (foundArg == 1)
			{
			fprintf(fdisplay, "%.3f\t", Slope);
			if (objectnumber == PRINTOBJECT)
				printf("%.3f\t", Slope);
			}
		/*Depth to top of frustum, to base and pedestal height*/
		CommandPatterni = 'A';
		CommandPatternii = 'h';
		foundArg = 0;
		ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
		if (foundArg == 1)
				{
				fprintf(fdisplay, "%.3f\t%.3f\t%.3f\t", topDepth, bottomDepth, PedestalHeight);
				if (objectnumber == PRINTOBJECT)
					{
					printf("%.3f\t%.3f\t%.3f\t", topDepth, bottomDepth, PedestalHeight);
					}
				}
		/*Number of data points in the feature*/
		CommandPatterni = 'M';
		CommandPatternii = 'd';
		foundArg = 0;
		ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
		if (foundArg == 1)
				{
				fprintf(fdisplay, "%d\t", finish - begin + 1);
				if (objectnumber == PRINTOBJECT)
					printf("%d\t", finish - begin + 1);
				}
		/*Misfit of best fit frustum*/
		CommandPatterni = 'A';
		CommandPatternii = 'm';
		foundArg = 0;
		ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
		if (foundArg == 1)
				{
				fprintf(fdisplay, "%.6f\t", MISFIT);
				if (objectnumber == PRINTOBJECT)
					printf("%.6f\t", MISFIT);
				}
		/*Longitude and latitude of each of the points of the best-fit frustum*/
		CommandPatterni = 'A';
		CommandPatternii = 'l';
		foundArg = 0;
		ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
		if (foundArg == 1)
			{
			for(f = 0; f <= 3; f ++)
				{
				fprintf(fdisplay, "%.3f\t%.3f\t", PointLongitude[f], PointLatitude[f]);
				if (objectnumber == PRINTOBJECT)
					printf("%.3f\t%.3f\t", PointLongitude[f], PointLatitude[f]);
				}
			}
		/*Long and lat of the centre of the flat top*/
		CommandPatterni = 'A';
		CommandPatternii = 'c';
		foundArg = 0;
		ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
		if (foundArg == 1)
				{
				fprintf(fdisplay, "%.3f\t%.3f\t", (PointLongitude[1] + PointLongitude[2])/2.0, (PointLatitude[1] + PointLatitude[2])/2.0);
				if (objectnumber == PRINTOBJECT)
					printf("%.3f\t%.3f\t", (PointLongitude[1] + PointLongitude[2])/2.0, (PointLatitude[1] + PointLatitude[2])/2.0);
				}
		/*Areas and volumes*/
		CommandPatterni = 'A';
		CommandPatternii = 'v';
		foundArg = 0;
		ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
		if (foundArg == 1)
				{
				/*fprintf(fdisplay, "%.3f\t%.3f\t%.3f\t", MeasuredArea, FrustumArea, FrustumVol);*/
				fprintf(fdisplay, "%.3f\t%.3f\t%.3f\t%.3f\t", MeasuredArea, FrustumArea, FrustumVol, (FrustumArea/MeasuredArea)*100.0);
				if (objectnumber == PRINTOBJECT)
					printf("%.3f\t%.3f\t%.3f\t", MeasuredArea, FrustumArea, FrustumVol);
				}
		/*Areas and volumes*/
		CommandPatterni = 'T';
		CommandPatternii = 'D';
		foundArg = 0;
		ReturnedPointer = (CommandLine(argc, argv, CommandPatterni, CommandPatternii, &foundArg));
		if (foundArg == 1)
				{
				/*fprintf(fdisplay, "%.3f\t%.3f\t", Angle, ExtraTravelFrac);*/
				fprintf(fdisplay, "%.3f\t%.3f\t", Angle, ExtraTravelFrac);
				if (objectnumber == PRINTOBJECT)
					printf("%.3f\t%.3f\t", Angle, ExtraTravelFrac);
				}
		/*and a newline to finish it off*/
		fprintf(fdisplay, "\n");
		/*Close the file again*/
		/*Inside loop as can only be closed if has been opened*/
		fclose(fproperties);
		}

	return 0;
}

/*printarrayFFFFF: prints the element of five float arrays and which element it
is.  Use for the INPUT, OUTPUT data*/
/*The definition*/
int printarrayFFFFF(float a[], float b[], float c[], float d[], float e[])
{
        int i;
        for (i = 0; i < MAXCOLUMNS; ++i)
                if (a[i])
                        printf("%d %f %f %f %f %f\n", i, a[i], b[i], c[i], d[i], e[i]);
        return 0;
}

/*ReadData*/
/*ReadData: Reads in $1 distance (km) $2 depth (m -ve down) $3 filtered depth (m -ve down)*/
/*                   $4 longitude (deg - decimal degrees i.e. 13.9E) $5 latitude (deg)*/
/*to arrays from file*/
int ReadData (char *argv[], float distance[], float depth[], float filtered[], float longitude[], float latitude[])

{

	float column1, column2, column3, column4, column5;
	int countlines;

	FILE *fdata;
	fdata = fopen (argv [1], "r");

	countlines = 0;
	while ((fscanf(fdata, "%f %f %f %f %f", &column1, &column2, &column3, &column4, &column5)) !=EOF)
		{
		distance[countlines] = column1;
		depth[countlines] = column2;
		filtered[countlines] = column3;
		longitude[countlines] = column4;
		latitude[countlines] = column5;
		countlines++;
		}

	fclose (fdata);

	return countlines;
}

/*RESAMPLE*/
/*Resample: resample the data at equal spacing along the line as it would be
drawn, this is so that the line and not the data points are fitted*/
/*Resample at equal sapcing along the lines - cf Book 8 pp53*/
float resample(int argc, char *argv[], float FracDist[], float FracDepth[], float FracDistRS[], float FracDepthRS[], float NumRS, int begin, int finish)

{
	int l,m,n;				/*loop variables*/
	float StepDist, LineDist;		/*Distance between points, and Distance along the line in object*/
	float ResampDist;			/*Distance (along line) at which any given point is resampled*/
	float ArrayPlace;			/*array element of the data arrays in interpolation*/
	float LineDistArray[MAXCOLUMNS];	/*the distance along line (in its seamount) for each data point*/
	float ResampFrac;			/*fraction of way between 2 data points that the resmapled
						point should be*/
	float ResampInc;



/*Work out the fractional distance and depth within the seamount of resampled points*/
/*Use the normalised shape*/
	/*Number of resampled points - includes fist and last points, so fencepost issues*/
	NumRS = 100.0;
	LineDist = 0.0;
	LineDistArray[begin] = 0;
	for (l = begin + 1; l <= finish; l++)
		{
		StepDist = (float)sqrt((FracDist[l-1] - FracDist[l])*(FracDist[l-1] - FracDist[l]) + (FracDepth[l-1] - FracDepth[l])*(FracDepth[l-1] - FracDepth[l]));
		/*printf("StepDist %f\n", StepDist);*/
		LineDist = LineDist + StepDist;
		LineDistArray[l] = LineDist;
		}
	/*to give NumRS - 1 points in addition to the first one*/
	ResampInc = LineDist/(NumRS - 1.0);
	/*printf("LineDist %f No. resamp points = %f Resamp Inc = %f\n", LineDist, NumRS, ResampInc);*/

	/*The actual resampling is pure linear interpolation*/
	/*Set the first and last points*/
	FracDistRS[0] = FracDist[begin];
	FracDepthRS[0] = FracDepth[begin];
	FracDistRS[(int)NumRS - 1] = FracDist[finish];
	FracDepthRS[(int)NumRS - 1] = FracDepth[finish];
	/*printf("FracDistRS[0] = %f, FracDepthRS[0] =%f, FracDistRS[%d] =%f, FracDepthRS[%d] = %f\n",FracDistRS[0], FracDepthRS[0], (int)NumRS - 1, FracDistRS[(int)NumRS - 1],(int)NumRS - 1, FracDepthRS[(int)NumRS - 1]);*/
	ArrayPlace = begin;
	for (m = 1; m <= (NumRS - 2); m++)
		{
		/*Determine the fractional Distance and Depth
		but must do this by finding which points it's between and how far*/
		ResampDist = m*ResampInc;
		/*Loop until just gone to point at beginning of interpolate pair*/
		n = ArrayPlace;
		while (LineDistArray[n] <= ResampDist)
			{
			n++;
			/*printf("ResampDist = %f,  LineDistArray[%d] = %f, LineDistArray[%d] = %f\n", ResampDist, n-1, LineDistArray[n-1], n, LineDistArray[n]);*/
			/*printf("n = %d\n", n);*/
			}
		/*printf("n = %d\n", n);*/
		/*The interpolation is now between array points n and n + 1*/
		/*The question is how far along*/
		ResampFrac = (ResampDist - LineDistArray[n-1])/(LineDistArray[n] - LineDistArray[n-1]);
		/*printf("Resamp between %f and %f - i.e. %f of the way between %d and %d\n", LineDistArray[n-1], LineDistArray[n], ResampFrac, n-1, n);*/
		/*So the resampled point is*/
		/*printf("Horizontal dist between %d and %d is %f, so %f percent of this is %f\n", n-1 , n, (FracDist[n]-FracDist[n-1]), 100*ResampFrac, (FracDist[n]-FracDist[n-1])*ResampFrac);*/
		FracDistRS[m] = FracDist[n-1] + (ResampFrac*(FracDist[n]-FracDist[n-1]));
		FracDepthRS[m] = FracDepth[n-1] + (ResampFrac*(FracDepth[n]-FracDepth[n-1]));
		/*printf("FracDistRS[%d] = %f, FracDepthRS[%d] = %f\n", m, FracDistRS[m],m, FracDepthRS[m]);*/
		/*The next resampled point may be in the same interval, but so as not to start the
		whole loop from the start set Array place back one*/
		ArrayPlace = n - 1;
		/*printf("\n");*/
		}
/*So now the resampled fractional arrays should be full from 0 to NumRS - 1*/
	/*printf("LineDist %f No. resamp points = %f Resamp Inc = %f\n", LineDist, NumRS, ResampInc);*/

	return ResampInc;
}
