#include <stdio.h>
#include <math.h>

/*
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
*/


//Specify a maximum crossing number of the tangles to be considered. This is used to intialize the size of the arrays used so that they are large enough to store all of the needed information.
#define MAXNN 15

//Global variables:
int gaussCrossingSigns[2*MAXNN];

//Does declaring an array of strings work in the same way?
//string ewingMillett[30][5];
int connectionsEM[2*MAXNN][2][5];

int numOfSubtangles;
int integerSubtangleConnectionsEM[2*MAXNN][2][5];
int integerSubtangleParametersEM[2*MAXNN][7];
int gaussIntegerSubtangleEM[2*MAXNN];
int barsGaussIntegerSubtangleEM[2];

int rationalSubtangleConnectionsEM[2*MAXNN][2][5];
int rationalSubtangleParametersEM[2*MAXNN][10];
int rationalTwistVectorsEM[2*MAXNN][2*MAXNN+1];
int gaussRationalSubtangleEM[2*MAXNN];
int barsGaussRationalSubtangleEM[2*MAXNN];

//Initialize an array to store information about which rational subtangles connect across endpoints after finding the rational planar diagram code.
int endpointConnectionsCornersRational[2][4];

//Initialize an array to store information about the canonical configurations of rational components, relative to a certain choice of the NW corner.
int rationalComponentsCanonicalConfiguration[MAXNN][4];

//int montesinosSubtangleConnectionsEM[2*MAXNN][2][5];
int montesinosSubtangleParametersEM[2*MAXNN][10];









/*
(NC, 7/22/18)
A quick recursive function to compute the greatest common denominator of two integers via the Euclidean algorithm.
Thank you to Wikipedia for the pseudocode!

THE STRANGE NAME OF THIS FUNCTION IS TO AVOID A FUTURE CONFLICT THAT MIGHT BE CAUSED BY MERGING THIS PROGRAM TO ONE THAT HAS THIS SAME FUNCTION ALREADY DEFINED

Function Inputs:
a: a nonegative integer.
b: another nonegative integer.

Function Outouts:
The greatest common divisor of a and b is returned as an integer output.

This function calls:
	nothing
This function is called by:
	calculatedExtendedFraction()
*/

int GCDem(int a, int b){
	if( b == 0 ){
		return a;
	} else {
		return GCDem(b, (a%b) );
	}
}



/*
(NC, 7/22/18)
This function will recover the corresponding rational number p/q for a rational tangle by simplifying the extended fraction defined by the twist vector.
This quantity is only meaningful for rational tangles, so this is checked as an input. The value of p and q uniquely classify the rational tangle.
These values will be stored as the first and second entries, respectively, of a global array rationalPQ[].
Note that these values should coincide with the invariants of the coloring matrix: a = p, c = q, and b = -p-q.
If p/q is negative, we will adopt the convention of making p positive and q negative.

THE STRANGE NAME OF THIS FUNCTION IS TO AVOID A FUTURE CONFLICT THAT MIGHT BE CAUSED BY MERGING THIS PROGRAM TO ONE THAT HAS THIS SAME FUNCTION ALREADY DEFINED

Necessary function inputs:
isRational: a boolean value to check whether the tangle under consideration is rational. If not, the function will exit.
twistVectLength: the length of the twsit vector array created in the buildTwistVector() function.
twistVector[]: the array of horizontal and vertical twists created the in buildTwistVector() function.

Function outputs (in the form of modified global variables):
rationalPQ[]: an array of two entries to store the values of p and q (the first and second entries, respectively).


This function calls:
	GCD()
This function is called by:
	findRationalTangles()

*/

void EMcalculateExtendedFractionPQ(bool isRational, int twistVectLength, int *twistVector, int *rationalPQ){
	
	if( isRational == true ){
		
		//Creat some local variables to store p and q. We will update the array rationalPQ[] at the end.
		int p;
		int q;
		
		//Create some local variables to store numerators, denominators, and GCDs if needed.
		int currentNumerator;
		int currentDenominator;
		int d;
		
		//In the special case that the twist vector has only one entry, define p and q accordingly.
		if( twistVectLength == 1 ){
			p = twistVector[0];
			q = 1;
			
		} else if( twistVectLength > 1 ){
			
			//Given a twist vector with two or more entries of the form (x1,x2,...,xn), consider the first part of the extended fraction: x1 + 1/x2 = (x1x2 + 1)/x2
			//Intialize the values of p and q by the values p = (x2x1 + 1)/d and q = x1/d, where d = GCD(x2x1+1,x1).
			currentNumerator = ( (twistVector[1]*twistVector[0]) + 1 );
			currentDenominator = twistVector[0];
			d = GCDem(currentNumerator,currentDenominator);
			
			p = currentNumerator/d;
			q = currentDenominator/d;
			
			//DEBUG
			//printf("\n p = %d, q = %d \n",p,q);
			
			//If the twist vector has more than two entries, the below while loop will move through all remaining portions of the extended fraction.
			//If the twist vector has only two entries, the while loop will not trigger.
			int k=2;
			while(k<twistVectLength){
				//At this point, p/q is the rationalization of the terms in the extended fraction considered so far.
				//The next term to consider has the following form: x + 1/(p/q) = x + q/p = (xp + q)/p.
				//Here, x denotes the current entry of twistVector[]. Use this formula to construct the next value of p and q.
				currentNumerator = ( (twistVector[k]*p) + q );
				currentDenominator = p;
				d = GCDem(currentNumerator,currentDenominator);
				
				p = currentNumerator/d;
				q = currentDenominator/d;
				
				//DEBUG
				//printf(" p = %d, q = %d \n",p,q);
				
				k++;
			}
			//By the end of this while loop, the extended fraction should be fully rationalized, with numerator p and denominator q.
			
		} else{
			//If the twist vector has a non-positive number of entries, then there is an error. Exit the function.
			return;
		}
		
		//Chcek the sign to follow the convention that p is always positive.
		if( p < 0 ){
			p = -1*p;
			q = -1*q;
		}
		
		//Update the global array with the values of p and q.
		rationalPQ[0] = p;
		rationalPQ[1] = q;
		
	} else {
		//The values of p and q are only meaningful for rational tangles. If the tangle is not rational, the function will exit.
		return;
	}
}


/*
(NC, 9/23/18)
Small function to print the information stored in the in the integerSubtangleParametersEM[][7] array.
This information is needed in various constructions, and this function is primarily intended for debugging purposes.
The rows in this array are indexed by the number of subtangles, and it is designed to store the following information for each subtangle:

[	0					, 	1				,	2								,	3			,	4			,	5				,	6				]
[ numOfCrossings/sign	,	twist axis type	,	horizontal/vertical twisting	,	fraction p	,	fraction q	,	internal parity	,	S2 direction	]

0: The absolute value is the number of crossings in the integer subtangle; these crossings all have the same orientation, denoted by the sign of this entry.
1: Twisting happens along either the (AD-BC) axis or along the (AC-BD) axis; indicate these by 1 and -1, respectively (this is equivalent to the "chainType" in the chaining function).
2: The twisting type, either horizontal (1) or vertical (-1), is determined internally for each subtangle on a case by case basis, following the convention that corner "A" denotes NW for the subtangle.
3: The numerator p of the corresponding rational fraction p/q, calculated depending on other parameters.
4: The denominator q of the corresponding rational fraction p/q, calculated depending on other parameters.
5: The internal parity (0, 1, or infinity(2)) of the given subtangle, following the convention that corner "A" denotes the NW. This IGNORES the orientation of the second strand.
6: Denotes the direction of the second strand in the subtangle, relative to the usual conventions for its parity: 1 if in the usual direction, or -1 if in the reverse direction.


Input:
numOfSubtanglesInput: the number of integer subtangles in the array (counting 1-crossing subtangles)
integerSubtangleParametersEM[][7]: described above

Outputs:
(only printf things)


This function calls:
	N/A
This function is called by:
	anywhere debugging is desired.

*/

void printIntegerSubtangleParametersEM(int numOfSubtanglesInput, int integerSubtangleParametersEM[][7]){
	
	printf("\n integerSubtangleParametersEM: \n");
	printf(" [ i ]\t [ 0\t 1\t 2\t 3\t 4\t 5\t 6\t]\n");
	printf(" sub- \t#cross \t twist \t horz/ \t frac \t frac \t sub- \t s2-");
	printf("\n tangle\t /sign \t axis \t vert \t p \t q \t parity\t direc \n\n");
	for(int i=0; i<numOfSubtanglesInput; i++){
		printf(" [ %d ]", i+1);
		printf("\t[");
		for(int j=0; j<7; j++){
			printf(" %d \t", integerSubtangleParametersEM[i][j]);
		}
		printf("]\n");
	}
	
}



/*
(NC, 10/8/18)
Small function to print the information stored in the in the rationalSubtangleParametersEM[][8] array.
This information is needed in various constructions, and this function is primarily intended for debugging purposes.
The rows in this array are indexed by the number of subtangles, and it is designed to store the following information for each subtangle:

[	0					, 	1								,	2					,	3			,	4			,	5				,	6				,	7						,	8		,	9				]
[ numOfCrossings		,	number of integer subgtangles	,	twist vector length	,	fraction p	,	fraction q	,	internal parity	,	S2 direction	,	subtangle "shape" sign	,	config	,	initial TV h/v	]

0: The number of crossings in the subtangle, unsigned
1: The number of integer sbutangles (counting simultaneous twisting?)
2: The number of nonzero entries in the rational twist vector
3: The numerator p of the corresponding rational fraction p/q, calculated depending on other parameters.
4: The denominator q of the corresponding rational fraction p/q, calculated depending on other parameters.
5: The internal parity (0, 1, or infinity(2)) of the given subtangle, following the convention that corner "A" denotes the NW. This IGNORES the orientation of the second strand.
6: Denotes the direction of the second strand in the subtangle, relative to the usual conventions for its parity: 1 if in the usual direction, or -1 if in the reverse direction.
7: The sign of the crossings in the corresponding rational tangle of this shape (which must have constant sign if irreducible), with the regular S2 direction.
8: The possible canonical configuration of the rational component (10 cases total--described elsewhere).
9: Whether the intial integer subtangle in the rational tangle is regarded as horizontl (1) or vertical (-1); this is used to construct the twist vector.


Input:
numOfSubtanglesInput: the number of integer subtangles in the array (counting 1-crossing subtangles)
rationalSubtangleParametersEM[][8]: described above

Outputs:
(only printf things)


This function calls:
	N/A
This function is called by:
	anywhere debugging is desired.

*/

void printRationalSubtangleParametersEM(int numOfSubtanglesInput, int rationalSubtangleParametersEM[][10]){
	
	printf("\n rationalSubtangleParametersEM: \n");
	printf(" [ i ]\t [ 0\t 1\t 2\t 3\t 4\t 5\t 6\t 7\t 8\t 9\t]\n");
	printf(" sub- \t#cross \t #int \t TwiVec\t frac \t frac \t sub- \t s2- \tsubtang\tcanon \tinitial");
	printf("\n tangle\t \t sub-t \t length\t p \t q \t parity\t direc \tsign \tconfig \tTV h/v\n\n");
	for(int i=0; i<numOfSubtanglesInput; i++){
		printf(" [ %d ]", i+1);
		printf("\t[");
		for(int j=0; j<10; j++){
			printf(" %d \t", rationalSubtangleParametersEM[i][j]);
		}
		printf("]\n");
	}
	
}



/*
(NC, 1/17/19)
Small function to print the information stored in the in the algebraicComponentsParameters[][10] array.
This information is needed in various constructions, and this function is primarily intended for debugging purposes.
The rows in this array are indexed by the number of subtangles, and it is designed to store the following information for each subtangle:

[	0				, 	1					,	2							,	3				,	4				,	5				,	6				,	7		,	8			,	9				]
[	joined stage	,	last stage index	,	number rational components	,	p / comp1 index	,	q / comp2 index	,	internal parity	,	S2 direction	,	config	,	J-sum/prod	, J-comp2 rotations	]

0: The stage at which this subtangle was joined from components, if it was joined. If the component is un-joined, this is 0.
1: The index of this subtangle in the preceding stage (if this subtangle was joined from two components last stage, the lower component index is the new subtangle index).
2: The number of rational components in the current "algebraic amalgam" subtangle.
3: If rational (one component), the numerator p of the fraction p/q; if algebraic, the index of the first joined component at the joined stage.
4: If rational (one component), the denominator q of the fraction p/q; if algebraic, the index of the second joined component at the joined stage.
5: The internal parity (0, 1, or infinity(2)) of the given subtangle, following the convention that corner "A" denotes the NW. This IGNORES the orientation of the second strand.
6: Denotes the direction of the second strand in the subtangle, relative to the usual conventions for its parity: 1 if in the usual direction, or -1 if in the reverse direction.
7: The possible canonical configuration of the subtangle (10 cases total--described elsewhere).
8: If a joined subtangle, whether the joining was a sum (1) or product (-1) (relative to the first component). This entry is 0 otherwise.
9: If a joined subtangle, the number of rotations applied to the second component before the joining.


Input:

numOfSubtanglesInput: the number of rational subtangle components in the original inp[ut array.
numOfStagesInput: the number of stages in the (joinings) used to build of the algebraic tangle.
algebraicComponentsParameters[][][10]: described above.

Outputs:
(only printf things)


This function calls:
	N/A
This function is called by:
	anywhere debugging is desired.

*/

void printAlgebraicComponentsParameters(int numOfSubtanglesInput, int numOfStagesInput, int algebraicComponentsParameters[MAXNN][MAXNN][10], bool lastStageOnly){
	
	//In some debugging cases, it might be helpful to see all stages at once, or only the most recent.
	int firstIndex;
	if( lastStageOnly ){
		firstIndex=numOfStagesInput;
	} else {
		firstIndex=0;
	}
	
	printf("\n algebraicComponentsParameters: \n");
	for(int i=firstIndex; i<=numOfStagesInput; i++){
		printf(" [ stage %d ]\n", i);
		printf(" [ i ]\t[ 0\t 1\t 2\t 3\t 4\t 5\t 6\t 7\t 8\t 9\t]\n");
		printf(" sub-\tjoin\tlast\t#ratio\tp or\tq or\tparity\tS2-\tcanon\tJ-sum/\tJ-comp2");
		printf("\n tangle\tstage\tstage\tcomps\tcomp1\tcomp2\t \tdirec\tconfig\tprod\trotates\n");
		for(int j=0; j<numOfSubtanglesInput-i; j++){
			printf(" [ %d ]", j+1);
			printf("\t[");
			for(int k=0; k<10; k++){
				printf(" %d \t", algebraicComponentsParameters[i][j][k]);
			}
			printf("]\n\n");
		}
	}
	
}




/*
(NC, 9/9/18)
Small function to print the entries in the connectionsEM arrays in an easy to follow format.
These store the information needed for the Ewing Millett code.
This is primarily for my own beneift with debugging.

Recall that the connectionsEM[][][] array is structured in the following way: connectionsEM[ crossing Gauss index ][2][5}

[ i ][ firt strand encountered in crossing over/under (+/-1) ][ Gauss index of the crossing connected to THIS INDEX corner of this crossing ]
	[ orientation of this crossing (+/-1)					][ corner (letter) of the crossing connected to THIS INDEX corner of this crossing ]

Inputs:
numOfCrossingsInput: the number of crossings of the currently considered tangle.
connectionsEM[][][]: the three dimensional array storing the information needed to reconstruct the Euling Millett array.

Outputs:
(only printf things)


This function calls:
	N/A
This function is called by:
	anywhere debugging is desired.

*/

void printArrayEM(int numOfCrossingsInput, int connectionsEMinput[][2][5]){
	
	printf("\n Ewing Millett Code Arrays: \n");
	for(int i=0; i<numOfCrossingsInput; i++){
		printf("\n [ %d ]", i+1);
		printf(" \t [ %d \t| ", connectionsEMinput[i][0][0]);
		for(int j=1; j<5; j++){
			printf(" %d ", connectionsEMinput[i][0][j]);
		}
		printf(" ]");
		printf("\n \t [ %d \t| ", connectionsEMinput[i][1][0]);
		for(int j=1; j<5; j++){
			printf(" %d ", connectionsEMinput[i][1][j]);
		}
		printf(" ]\n");
	}
	
}

/*
(NC, 9/9/18)
Small function to print the entries in the connectionsEM arrays in the standard written format of the Ewing Millett code.

Recall that the connectionsEM[][][] array is structured in the following way: connectionsEM[ crossing Gauss index ][2][5}

[ i ][ firt strand encountered in crossing over/under (+/-1) ][ Gauss index of the crossing connected to THIS INDEX corner of this crossing ]
	[ orientation of this crossing (+/-1)					][ corner (letter) of the crossing connected to THIS INDEX corner of this crossing ]

Inputs:
numOfCrossingsInput: the number of crossings of the currently considered tangle.
connectionsEM[][][]: the three dimensional array storing the information needed to reconstruct the Euling Millett array.

Outputs:
(only printf things)


This function calls:
	N/A
This function is called by:
	anywhere debugging is desired.

*/

void printStandardEM(int numOfCrossingsInput, int connectionsEM[][2][5]){
	
	printf("\n Standard Ewing Millett Code: \n");
	for(int i=0; i<numOfCrossingsInput; i++){
		//Print the orientation of this crossing.
		if( connectionsEM[i][1][0] > 0 ){
			printf("\t + ");
		} else {
			printf("\t - ");
		}
		for(int j=1; j<5; j++){
			//Print the Gauss index of the connected crossing at this corner
			printf("%d", connectionsEM[i][0][j]);
			//Print the letter corresponding this corner; four cases.
			int cornerLetter=connectionsEM[i][1][j];
			if( cornerLetter == 1 ){
				printf("a ");
			} else if ( cornerLetter == 2 ){
				printf("b ");
			} else if ( cornerLetter == 3 ){
				printf("c ");
			} else if ( cornerLetter == 4 ){
				printf("d ");
			} else {
				//If none of the above happens, there is a some bug somewhere. Flag this "Z" so it can identified at a glance.
				printf("Z ");
			}
		}
		printf("\n");
		
	}
	
}



/*
(NC, 9/16/18)
Small function to print the integerSubtangleEM arrays in a pretty format, simialr to the standard written EM code.
To avoid possible confusions with the labeling convention for corners of subtangles rather than crossings, capital letters will be used in this case.

Recall that the integerSubtangleConnectionsEM[][][] array is structured in the following way: integerSubtangleConnectionsEM[ crossing Gauss index ][2][5}

[ i ][ firt strand encountered in crossing over/under (+/-1) ][ GaussEM index of the crossing connected to THIS INDEX corner of this crossing ]
	[ number and sign of crossings in a subtangle chain		][ corner (letter) of the crossing connected to THIS INDEX corner of this crossing ]

Inputs:
numOfCrossingsInput: the number of crossings of the currently considered tangle.
integerSubtangleConnectionsEMinput[][][]: the three dimensional array storing the information needed to reconstruct the generalized EM code for integer subtangles.


This function calls:
	N/A
This function is called by:
	anywhere debugging is desired.

*/

void printPrettyIntegerSubtangleEM(int numOfSubtanglesInput, int integerSubtangleConnectionsEMinput[][2][5], int integerSubtangleParametersEM[][7]){
	
	printf("\n Integer Subtangle EM Code (pretty): \n");
	for(int i=0; i<numOfSubtanglesInput; i++){
		//Print fraction p/q corresponding to this subtangle.
		printf("\t[ %d / %d ]\t", integerSubtangleParametersEM[i][3],integerSubtangleParametersEM[i][4]);
		for(int j=1; j<5; j++){
			//Print the GaussEM index of the connected crossing at this corner
			printf("%d", integerSubtangleConnectionsEMinput[i][0][j]);
			//Print the letter corresponding this corner; four cases.
			int cornerLetter=integerSubtangleConnectionsEMinput[i][1][j];
			if( cornerLetter == 1 ){
				printf("A ");
			} else if ( cornerLetter == 2 ){
				printf("B ");
			} else if ( cornerLetter == 3 ){
				printf("C ");
			} else if ( cornerLetter == 4 ){
				printf("D ");
			} else {
				//If none of the above happens, there is a some bug somewhere. Flag this "z" so it can identified at a glance.
				printf("z ");
			}
		}
		
		//Print the internal parity of the subtangle and the whether or not the second subtangle strand is reversed relative to this parity.
		if( integerSubtangleParametersEM[i][5] == 2 ){
			printf("\t (Parity: infinity");
		} else {
			printf("\t (Parity: %d", integerSubtangleParametersEM[i][5]);
		}
		if( integerSubtangleParametersEM[i][6] < 0 ){
			printf("')\n");
		} else {
			printf(")\n");
		}
		
	}
	
}






/*
(NC, 10/15/18)
Small function to print the rationalSubtangleEM arrays in a pretty format, simialr to the standard written EM code.
To avoid possible confusions with the labeling convention for corners of subtangles rather than crossings, capital letters will be used in this case.

Recall that the rationalSubtangleConnectionsEM[][][] array is structured in the following way: integerSubtangleConnectionsEM[ crossing Gauss index ][2][5}

[ i ][ firt strand encountered in crossing over/under (+/-1) ][ GaussEM index of the crossing connected to THE INITIAL INDEX of this subtangle ]
	[ number and sign of crossings in a subtangle chain		][ corner (letter) of the subtangle connected to this subtangle ]

Inputs:
numOfCrossingsInput: the number of crossings of the currently considered tangle.
integerSubtangleConnectionsEMinput[][][]: the three dimensional array storing the information needed to reconstruct the generalized EM code for integer subtangles.


This function calls:
	N/A
This function is called by:
	anywhere debugging is desired.

*/

void printPrettyRationalSubtangleEM(int numOfSubtanglesInput, int rationalSubtangleConnectionsEMinput[][2][5], int rationalSubtangleParametersEM[][10], int twistVectorArray[][2*MAXNN+1]){
	
	printf("\n Rational Subtangle EM Code (pretty): \n");
	for(int i=0; i<numOfSubtanglesInput; i++){
		//Print fraction p/q corresponding to this subtangle.
		printf("\t[ %d / %d ]\t", rationalSubtangleParametersEM[i][3],rationalSubtangleParametersEM[i][4]);
		for(int j=1; j<5; j++){
			//Print the GaussEM index of the connected subtangle at this corner
			printf("%d", rationalSubtangleConnectionsEMinput[i][0][j]);
			//Print the letter corresponding this corner; four cases.
			int cornerLetter=rationalSubtangleConnectionsEMinput[i][1][j];
			if( cornerLetter == 1 ){
				printf("A ");
			} else if ( cornerLetter == 2 ){
				printf("B ");
			} else if ( cornerLetter == 3 ){
				printf("C ");
			} else if ( cornerLetter == 4 ){
				printf("D ");
			} else if ( cornerLetter == 0 ){
				printf("E ");
			} else {
				//If none of the above happens, there is a some bug somewhere. Flag this "z" so it can identified at a glance.
				printf("z ");
			}
		}
		printf("\t");
		
		//Not needed to for the notation, but could be helpful for debugging.
		//printf("Config: %d\t", rationalSubtangleParametersEM[i][8]);
		
		//Print the internal parity of the subtangle and the whether or not the second subtangle strand is reversed relative to this parity.
		if( rationalSubtangleParametersEM[i][5] == 2 ){
			printf("(P: infinity");
		} else {
			printf("(P: %d", rationalSubtangleParametersEM[i][5]);
		}
		if( rationalSubtangleParametersEM[i][6] < 0 ){
			printf("')\t");
		} else {
			printf(")\t");
		}
		if( rationalSubtangleParametersEM[i][5] != 2 ){
			printf("\t");
		}
		
		//Print the twist vector--this is technically redundant with the fraction, but is useful information at a glance for reconstructing the subtangle by hand.
		printf("TwiVec: [");
		for(int j=1; j<=twistVectorArray[i][0]; j++){
			printf(" %d", twistVectorArray[i][j]);
		}
		printf(" ]\n");
		
	}
	
}



/*
(NC, 1/20/19)
Small function to print the rationalSubtangleEM arrays with compass directions.
Note that, in order to use compass labels correctly, this should only be called with an appropriate array.

Recall that the rationalSubtangleConnectionsCompass[][][] array is structured in the following way:

[ i ][ canon config ][ GaussEM index of the crossing connected to THE INITIAL INDEX of this subtangle ]
	[ rotations NW	][ corner (letter) of the subtangle connected to this subtangle ]

Inputs:
numOfCrossingsInput: the number of crossings of the currently considered tangle.
rationalSubtangleConnectionsCompassInput[][][]: the three dimensional array storing the information needed to reconstruct the generalized EM code for rational subtangles.
rationalSubtangleParametersEM[][]


This function calls:
	N/A
This function is called by:
	anywhere debugging is desired.

*/

void printCompassRationalSubtangleEM(int numOfSubtanglesInput, int rationalSubtangleConnectionsCompassInput[][2][5], int rationalSubtangleParametersEM[][10], int twistVectorArray[][2*MAXNN+1]){
	
	printf("\n Rational Subtangle COMPASS Code: \n");
	for(int i=0; i<numOfSubtanglesInput; i++){
		//Print fraction p/q corresponding to this subtangle.
		printf("\t[ %d / %d ]\t", rationalSubtangleParametersEM[i][3],rationalSubtangleParametersEM[i][4]);
		for(int j=1; j<5; j++){
			//Print the GaussEM index of the connected subtangle at this corner
			printf("%d ", rationalSubtangleConnectionsCompassInput[i][0][j]);
			//Print the letter corresponding this corner; four cases.
			int cornerLetter=rationalSubtangleConnectionsCompassInput[i][1][j];
			if( cornerLetter == 1 ){
				printf("NW\t");
			} else if ( cornerLetter == 2 ){
				printf("NE\t");
			} else if ( cornerLetter == 3 ){
				printf("SE\t");
			} else if ( cornerLetter == 4 ){
				printf("SW\t");
			} else if( cornerLetter == 0 ){
				printf("EP\t");			
			} else {
				//If none of the above happens, there is a some bug somewhere. Flag this "z" so it can identified at a glance.
				printf("z\t");
			}
		}
		//printf("\t");
		
		//Not needed to for the notation, but could be helpful for debugging.
		//printf("Config: %d\t", rationalSubtangleParametersEM[i][8]);
		
		//Print the internal parity of the subtangle and the whether or not the second subtangle strand is reversed relative to this parity.
		if( rationalSubtangleParametersEM[i][5] == 2 ){
			printf("(P: infinity");
		} else {
			printf("(P: %d", rationalSubtangleParametersEM[i][5]);
		}
		if( rationalSubtangleParametersEM[i][6] < 0 ){
			printf("')\t");
		} else {
			printf(")\t");
		}
		if( rationalSubtangleParametersEM[i][5] != 2 ){
			printf("\t");
		}
		
		//Print the twist vector--this is technically redundant with the fraction, but is useful information at a glance for reconstructing the subtangle by hand.
		printf("TwiVec: [");
		for(int j=1; j<=twistVectorArray[i][0]; j++){
			printf(" %d", twistVectorArray[i][j]);
		}
		printf(" ]\n");
		
	}
	
}




/*
(NC, 1/12/19)
This function will use the labeling of the corners of a subtangle in the rational planar diagram code to determine the corresponding compass direction corner labeling.
There are three cases, depending on parity. In all cases, we assume that corner A is identified as NW.
This will also rotate the compass labels so that whatever original label is specified for newNWcorner is now in the NW position.

Inputs:
pairtyInput: the parity of the rational subtangle being considered, either 0, 1, or 2(infinity).
S2direcInput: the direction of the S2 strand of the input rational subtangle, either 1 (usual direction) or -1 (reverse direction).
newNWcorner: the original label of the corner which must now be rotated into the NW position, either 1, 2, 3, or 4 (corresponding to A, B, C, or D, respectively).
configurationInput: the canonical configuration of the input rational subtangle, as determined by the RPD code.

Outputs (in the form of modified array entries):
cornersOutput[]: an array to storing the old labels which now match the new compass positions. It will have the form [ * , NW , NE , SE , SW , CanonConfig , rotateNW ]. Entry 0 is left unmodified.


This function calls:
	N/A
This function is called by:
	detectIfMontesinos()

*/

void determineCompassCorners(int parityInput, int S2direcInput, int newNWcorner, int configurationInput, int cornersOutput[7]){
	
	//Initialize a variable to determine whether or not we need to change the NW corner.
	//The compass labels will be rotated so that whatever label matches newNWcorner is now in the NW position.
	//By default, set this to 0; if newNWcorner == 1 (A), then no rotation is necessary.
	int rotateNW=0;
	
	//In all case, NW = A by default.
	cornersOutput[1]=1;
	
	if( parityInput == 0 ){
		//Case 1: parity 0
		//NE = B
		cornersOutput[2]=2;
		
		if( newNWcorner == 2 ){
			rotateNW=1;
		}
		
		if( S2direcInput > 0 ){
			//If usual S2 direction, SE = C, SW = D
			cornersOutput[3]=3;
			cornersOutput[4]=4;
			
			if( newNWcorner == 3 ){
				rotateNW=2;
			} else if( newNWcorner == 4 ){
				rotateNW=3;
			}
			
		} else {
			//If reverse S2 direction, SE = D, SW = C
			cornersOutput[3]=4;
			cornersOutput[4]=3;
			
			if( newNWcorner == 3 ){
				rotateNW=3;
			} else if( newNWcorner == 4 ){
				rotateNW=2;
			}	
		}
		
	} else if( parityInput == 1 ){
		//Case 2: parity 1
		//SE = B
		cornersOutput[3]=2;
		
		if( newNWcorner == 2 ){
			rotateNW=2;
		}
		
		if( S2direcInput > 0 ){
			//If usual S2 direction, NE = C, SW = D
			cornersOutput[2]=3;
			cornersOutput[4]=4;
			
			if( newNWcorner == 3 ){
				rotateNW=1;
			} else if( newNWcorner == 4 ){
				rotateNW=3;
			}
			
		} else {
			//If reverse S2 direction, NE = D, SW = C
			cornersOutput[2]=4;
			cornersOutput[4]=3;
			
			if( newNWcorner == 3 ){
				rotateNW=3;
			} else if( newNWcorner == 4 ){
				rotateNW=1;
			}
		}
	} else if( parityInput == 2 ){
		//Case 3: parity infinity
		//SW = B
		cornersOutput[4]=2;
		
		if( newNWcorner == 2 ){
			rotateNW=3;
		}
		
		if( S2direcInput > 0 ){
			//If usual S2 direction, SE = C, NE = D
			cornersOutput[3]=3;
			cornersOutput[2]=4;
			
			if( newNWcorner == 3 ){
				rotateNW=2;
			} else if( newNWcorner == 4 ){
				rotateNW=1;
			}
			
		} else {
			//If reverse S2 direction, SE = D, NE = C
			cornersOutput[3]=4;
			cornersOutput[2]=3;
			
			if( newNWcorner == 3 ){
				rotateNW=1;
			} else if( newNWcorner == 4 ){
				rotateNW=2;
			}
		}
	} else {
		//Error case, shouldn't happen. Set all output corners to 0.
		cornersOutput[1]=0;
		cornersOutput[2]=0;
		cornersOutput[3]=0;
		cornersOutput[4]=0;
	}
	
	//Initialize variables for the old corners.
	int oldCorners[4];
	//Rotate the corner labels; repeat as needed.
	for(int i=0; i<rotateNW; i++){
		//Store the old corner information.
		oldCorners[0]=cornersOutput[1];
		oldCorners[1]=cornersOutput[2];
		oldCorners[2]=cornersOutput[3];
		oldCorners[3]=cornersOutput[4];
		//Update the rotated corner information.
		cornersOutput[1]=oldCorners[1];
		cornersOutput[2]=oldCorners[2];
		cornersOutput[3]=oldCorners[3];
		cornersOutput[4]=oldCorners[0];
	}
	
	//If there is rotating, there are a handful of cases for the new configuration after the rotations.
	//The configuration labels are essentially arbitrarily, so this is somewhat inelegant.
	//However, there are basically three cases, depending on the original configuration.
	if( configurationInput <= 0 ){
		//Case 1: -1 (integer subtangle, all configurations possible) or 0 (no canonical configurations), the configuration stays the same.
		cornersOutput[5] = configurationInput;
	} else if( (0<configurationInput) && (configurationInput<=4) ){
		//Case 2: 1, 2, 3, or 4; one configuration, using as NW either A, A*90deg, A*180deg, or A*270deg, respectively.
		//The new configuration is exactly rotated by rotateNW positions.
		cornersOutput[5] = ((configurationInput-rotateNW)%4);
		//Since the above returns the remainder upon division by 4, the configuration should actually be 4 when the above computation is 0.
		if( cornersOutput[5] == 0 ){
			cornersOutput[5]=4;
		}
	} else if( 4 < configurationInput ){
		//Case 3: 5, 6, 7, 8; two adjacent configurations, the first of which usese as NW either A, A*90deg, A*180deg, or A*270deg, respectively.
		cornersOutput[5] = 4+((configurationInput-rotateNW)%4);
		//By the same logic as above, the configuration should actually be 8 when the above computation is 4.
		if( cornersOutput[5] == 4 ){
			cornersOutput[5]=8;
		}
	}
	
	//Update entry [6] of the cornersOutput[] array to denote the number of rotations used.
	//This information is used to determine whether to use the fraction p/q or -q/p when giving the Montesinos construction.
	cornersOutput[6]=rotateNW;
	
	//DEBUG
	//printf("\n DEBUG: Line %d: rational component %d\n Old Configuration: %d , New Configuration: %d , Rotations: %d\n", __LINE__, cornersOutput[0]+1, configurationInput, cornersOutput[5], cornersOutput[6]);
	
}





/*
(NC, 9/9/18)
The following three functions are taken/adapted from trans.cpp. They are used to construct the Ewling Millet code of a tangle from its Gauss code.
Rather than delete something I don't fully understand, I am choosing to copy the parts I expect will work here to write a program that will do something similar.

*/

/*
ADAPTED FROM: findNextLetter in trans.cpp

This function calls:
	N/A
This function is called by:
	findEMconnections()
*/
int findNextLetter2(int posOrNeg, int crossingTypeInput)
{
	if (posOrNeg == 1)
	{
		return 3;
	}
	else
	{
		if (crossingTypeInput == 1)
		{
			return 2;
		}
		else
		{
			return 4;
		}
	}
}


/*
ADAPTED FROM: findCurrentLetter in trans.cpp

This function calls:
	N/A
This function is called by:
	findEMconnections()
*/
int findCurrentLetter2(int posOrNeg, int crossingTypeInput)
{
	if (posOrNeg == 1)
	{
		return 1;
	}
	else
	{
		if (crossingTypeInput == 1)
		{
			return 4;
		}
		else
		{
			return 2;
		}
	}
}


/*
(NC, 9/9/18)
This function constructs the three dimensional array connectionsEM[][2][5], which stores all of the information of the Ewing Millett code of the tangle.
It is built up from the Gauss code.

ADAPTED FROM: gaussToEwingMillet in trans.cpp

Inputs:
gaussInput[]: the usual Gauss code array of the tangle being considered.
crossingTypeInput[]: an array of length numOfCrossingsInput storing the oriented sign of each crossing, where the entries match the Gauss code index of the crossing.
numOfCrossingsInput: the number of crossings of the tangle being considered.

Outputs (in the form of modified global variables):
connectionsEM[numOfCrossingsInput][2][5]: a three dimensional array to store the information needed to construct the Euling Millett code of the tangle.
The general layout of this multi-dimensional array is as follows:

[ i ][ firt strand encountered in crossing over/under (+/-1) ][ Gauss index of the crossing connected to THIS INDEX corner of this crossing ]
	[ orientation of this crossing (+/-1)					][ corner (letter) of the crossing connected to THIS INDEX corner of this crossing ]

Above, i denotes the Gauss code index of the crossing.


This function calls:
	printArrayEM() (DEBUG)
	printStandardEM() (DEBUG)
This function is called by:
	buildGeneralizedEMCode()

*/

void findEMConnections(int* gaussInput, int* crossingTypeInput, int connectionsEM[][2][5], int numOfCrossingsInput){
	
	//Intialize an array to check the first time each crossing is encountered in the gaussInput[] array.
	//This is used to determine whether we first encounter the crossing along the overstrand or the understrand.
	//These are initially set to 0; after the crossing is encountered, the entry is set to 1.
	int checkCrossingFirstEncounter[numOfCrossingsInput] = {0};
	
	for (int i = 0; i < 2*numOfCrossingsInput; i++)
	{
		//The Gauss code index of this crossing is just the absolute value of this entry in the gaussInput[] array.
		int gaussIndex = abs(gaussInput[i]);
		
		//Check to see if the Gauss code indexing of this crossing has been encountered before.
		if( checkCrossingFirstEncounter[ gaussIndex - 1 ] == 0 ){
			//If not, set the entry tracking whether this crossing is first encountered along the overstrand or understrand (first column, first row)
			//If the gaussInput entry is positive, we encounter along the overstrand; otherwise, we encounter along the understrand.
			if( gaussInput[i] > 0 ){
				connectionsEM[ gaussIndex - 1 ][0][0] = 1;
			} else {
				connectionsEM[ gaussIndex - 1 ][0][0] = -1;
			}
			
			//Set the entry tracking this crossing's orientation (first colum, second row)
			if (crossingTypeInput[ gaussIndex - 1 ] == 1){
				connectionsEM[ gaussIndex - 1 ][1][0] = 1;
			} else {
				connectionsEM[ gaussIndex - 1 ][1][0] = -1;
			}
			
			//Update the check array to track that we have now encountered this crossing once already.
			checkCrossingFirstEncounter[ gaussIndex - 1 ] = 1;
		}
	}

	for (int i = 0; i < (2 * numOfCrossingsInput); i++)
	{
		int currentGauss = gaussInput[i];
		int nextGauss;
		if (i + 1 < (2 * numOfCrossingsInput)) {
			nextGauss = gaussInput[i + 1];
		} else {
			//If the current crossing happens to be the last one encountered on the srand, we loop back aground to the very first crossing.
			//Note that this applies for all three parites (0, 1, and infinity) by how we have chosen to take the closure.
			nextGauss = gaussInput[0];
		}
		
		//This determines the corner (letter in EM notation) of the current and the next crossing according the convention for the Euling Millet code.
		//The corner matching the outward pointing overstrand is always assigned letter a (corner 1); the remaining three corners are assigned subsequent letters moving clockwise from this corner.
		int currentLetter = findCurrentLetter2(currentGauss / abs(currentGauss), crossingTypeInput[abs(currentGauss) - 1]);
		int nextLetter = findNextLetter2(nextGauss / abs(nextGauss), crossingTypeInput[abs(nextGauss) - 1]);
		
		//Update the two corresponding entries of connectionsEM with this information on how the corners are connected.
		//Note that we update two entries at a time because connections always come in pairs.
		
		//CURRENT CROSSING:
		//FIRST ROW: Gauss index of the connected crossing (nextGauss), connected to the current corner (currentLetter)
		connectionsEM[ abs(currentGauss) - 1 ][0][ currentLetter ] = abs(nextGauss);
		//SECOND ROW: corner (nextLetter) of the connected crossing (nextGauss), connected to the current corner (currentLetter)
		connectionsEM[ abs(currentGauss) - 1 ][1][ currentLetter ] = nextLetter;
		
		//NEXT CROSSING :
		//FIRST ROW: Gauss index of the connected crossing (currentGauss), connected to the next corner (nextLetter)
		connectionsEM[ abs(nextGauss) - 1 ][0][ nextLetter ] = abs(currentGauss);
		//SECOND ROW: corner (currentLetter) of the connected crossing (currentGauss), connected to the next corner (nextLetter)
		connectionsEM[ abs(nextGauss) - 1 ][1][ nextLetter ] = currentLetter;
		
	}
	
	//DEBUG:
	//printArrayEM(numOfCrossingsInput, connectionsEM);
	//printStandardEM(numOfCrossingsInput, connectionsEM);
	
}


/*
(NC, 9/9/18)
This function creates an array with numOfCrossingsInput number of entries to store the sign (1 or -1) of each crossing in the gaussInput array.
The index is set to match the "Gauss code index" of each crossing, so that crossing i in the Gauss code has sign gaussCrossingSigns[i-1].
Note that this information is equivalent to what is stored in orientedSignsGaussInput[], but truncated and indexed differently.

Inputs:
gaussInput[]: the usual array storing the Gauss code the tangle.
oirentedSignGaussInput[]: the array storing the information of the crossing signs for the Gauss code of the input tangle.
numOfCrossingsInput: the number of crossings in the input tangle

Outputs (in the form of modified globabl variables):
gaussCrossingSigns[]: array of length numOfCrossingsInput to store the crossing signs, as described above.


This function calls:
	N/A
This function is called by:
	buildGeneralizedEMCode()

*/

void findGaussCrossingSigns(int *gaussInput, int *orientedSignGaussInput, int numOfCrossingsInput, int *gaussCrossingSigns){
	
	//Each crossing in the Gauss code array occurs twice in the array: once with an even index and once with an odd index.
	//The loop below merely iterates through the entries of even index, and checks the orientedSignGaussInput array to determine sign.
	//The index of gaussCrossingSigns is hence the Gauss code index of the crossing (-1 to start with 0), with sign determined by the corresponding orientedSignGaussInput entry.
	for(int i=0; i<numOfCrossingsInput; i++){
		gaussCrossingSigns[ abs(gaussInput[2*i]) - 1 ] = orientedSignGaussInput[2*i];
	}
	
	//DEBUG
	/*
	printf("\n gaussCrossingSigns: \n [");
	for(int i=0; i<numOfCrossingsInput; i++){
		printf(" %d ", gaussCrossingSigns[i]);
	}
	printf("]\n");
	*/
	
	
}



/*
(NC, 9/9/18)
This function reindexes the corners of the crossings according to a convention more suitable for subtangles.
In the usual Euling Millett notation, the exit of the overstrand is assigned the letter "a", and the remaining corners are assigned "b", "c", "d" by moving clockwise from "a".
This convention doesn't make sense in general for subtangles, so I am adopting a new one based on the entrance/exit of the two strands in the subtangle.
I will use the following convention: "a" = S1 enter; "b" = S1 exit; "c" = S2 enter; "d" = S2 exit.
This is just a permutation of the corners of each crossing, now regarded as a 1-crossing subtangle.
There are four possible cases, depending on whether the tangle is first entered along the overstrand or the understrand, and the sign of the crossing.
This information is stored in connectionsEM[i][0][0] and connectionsEM[i][1][0], respectively.
The NEW corners that the OLD corners get mapped to are specified via integers a, b, c, and d for each of these four cases (that is, the permuation).

Inputs:
numOfCrossingsInput: the number of crossings of the tangle being considered.
connectionsEM[][2][5]: the three dimensional array storing the information needed to construct the Euling Millett code of the tangle.

Outputs (in the form of modified global variables):
integerSubtangleConnectionsEM[][2][5]: the three dimensional array with the labeling reindexed to match the description above.

This function calls:
	printArrayEM() (DEBUG)
This function is called by:
	buildIntegerSubtangleEMCode()

*/

void reindexEMforSubtangles(int numOfCrossingsInput, int connectionsEM[][2][5], int integerSubtangleConnectionsEM[][2][5]){
	
	//Specify integers to store NEW corners.
	int a;
	int b;
	int c;
	int d;
	//These are set to 1, 2, 3, or 4 (representing a, b, c, or d, respectively) to describe the NEW corner mapped to.
	
	for(int i=0; i<numOfCrossingsInput; i++){
		
		//CASE 1: enter along overstrand, positive crossing
		if( (connectionsEM[i][0][0] == 1) && (connectionsEM[i][1][0] == 1) ){
			a=2;
			b=3;
			c=1;
			d=4;
		//CASE 2: enter along understrand, negative crossing
		} else if( (connectionsEM[i][0][0] == 1) && (connectionsEM[i][1][0] == -1) ){
			a=2;
			b=4;
			c=1;
			d=3;
		//CASE 3: enter along understrand, positive crossing	
		} else if( (connectionsEM[i][0][0] == -1) && (connectionsEM[i][1][0] == 1) ){
			a=4;
			b=1;
			c=3;
			d=2;
		//CASE 4: enter along understrand, positive crossing	
		} else if( (connectionsEM[i][0][0] == -1) && (connectionsEM[i][1][0] == -1) ){
			a=4;
			b=2;
			c=3;
			d=1;
		} else {
			//If some reason none of these things happened, which should not happen unless there is a bug, just set the permutation to be "identity".
			a=1;
			b=2;
			c=3;
			d=4;
		}
		
		//The information regarding whether we enter the strand on the over/understrand and the crossing orienation are the same in this case.	
		integerSubtangleConnectionsEM[i][0][0] = connectionsEM[i][0][0];
		integerSubtangleConnectionsEM[i][1][0] = connectionsEM[i][1][0];
		
		//Row 0 of this array stores the Gauss index of the connected crossings; these indices don't change, but their positions in the row do.
		//The NEW corners are set to equal the OLD corners, as specified by the chosen permutation.
		integerSubtangleConnectionsEM[i][0][a] = connectionsEM[i][0][1];
		integerSubtangleConnectionsEM[i][0][b] = connectionsEM[i][0][2];
		integerSubtangleConnectionsEM[i][0][c] = connectionsEM[i][0][3];
		integerSubtangleConnectionsEM[i][0][d] = connectionsEM[i][0][4];
		
		//Row 1 of this array stores the corner of the connected crossings; the positions AND the labels change.
		//The NEW corners are set equal to the OLD corners in the same way as above.
		integerSubtangleConnectionsEM[i][1][a] = connectionsEM[i][1][1];
		integerSubtangleConnectionsEM[i][1][b] = connectionsEM[i][1][2];
		integerSubtangleConnectionsEM[i][1][c] = connectionsEM[i][1][3];
		integerSubtangleConnectionsEM[i][1][d] = connectionsEM[i][1][4];
		
		//Unfortunately, accounting for the label changing in Row 1 is not as straightforward since each connected crossing can have a different permutation.
		//To work around this, search through ALL Row1 entries for every crossing to see if they match this crossing index; this crossing must show up in exactly 4 connections.
		//Whenever we find a connection, we will need relabel it according to the case found here.
		//HOWEVER, not all of the corresponding entries have even been defined yet! We can't do this relabeling until the entire array is defined.
		//This doesn't seem efficient, but it should make the re-labeling work, while the above takes care of the rearranging.
			
	}
	
	
	//A final for loop to relabel the connected corners in Row1, as described above.
	//This seems SUPER inefficient.
	for(int i=0; i<numOfCrossingsInput; i++){
		//CASE 1: enter along overstrand, positive crossing
		if( (connectionsEM[i][0][0] == 1) && (connectionsEM[i][1][0] == 1) ){
			a=2;
			b=3;
			c=1;
			d=4;
		//CASE 2: enter along understrand, negative crossing
		} else if( (connectionsEM[i][0][0] == 1) && (connectionsEM[i][1][0] == -1) ){
			a=2;
			b=4;
			c=1;
			d=3;
		//CASE 3: enter along understrand, positive crossing	
		} else if( (connectionsEM[i][0][0] == -1) && (connectionsEM[i][1][0] == 1) ){
			a=4;
			b=1;
			c=3;
			d=2;
		//CASE 4: enter along understrand, positive crossing	
		} else if( (connectionsEM[i][0][0] == -1) && (connectionsEM[i][1][0] == -1) ){
			a=4;
			b=2;
			c=3;
			d=1;
		} else {
			//If some reason none of these things happened, which should not happen unless there is a bug, just set the permutation to be "identity".
			a=1;
			b=2;
			c=3;
			d=4;
		}
		
		for(int j=0; j<numOfCrossingsInput; j++){
			for(int k=1; k<5; k++){
				//Check Row0 first to see if the Gauss code index of the crossing matches i.
				if(integerSubtangleConnectionsEM[j][0][k] == (i+1)){
					//If it does, re-label the Row1 entry that specifies the corner connection.
					//This should happen exactly four times, once per corner.
					if( integerSubtangleConnectionsEM[j][1][k] == 1 ){
						integerSubtangleConnectionsEM[j][1][k] = a;
					} else if( integerSubtangleConnectionsEM[j][1][k] == 2 ){
						integerSubtangleConnectionsEM[j][1][k] = b;
					} else if( integerSubtangleConnectionsEM[j][1][k] == 3 ){
						integerSubtangleConnectionsEM[j][1][k] = c;
					} else if( integerSubtangleConnectionsEM[j][1][k] == 4 ){
						integerSubtangleConnectionsEM[j][1][k] = d;
					}
				}
			}
		}
	}
	
	
	//DEBUG
	/*
	printf("Orignal EM Code:");
	printArrayEM(numOfCrossingsInput,connectionsEM);
	printf("New EM Code, with corners labeled by: a=S1enter, b=S1exit, c=S2enter, d=S2exit");
	printArrayEM(numOfCrossingsInput,integerSubtangleConnectionsEM);
	*/
	
}



/*
(NC, 9/9/18)
Function to chain crossings together while building the integer subtangle variant of the Euling Millet code.
(Work in progress--rewrite description later--mention where this is called?

Inputs:

Outputs:


This function calls:
	N/A
This function is called by:
	buildIntegerSubtanlgeEMCode()

*/

void chainCrossings(int &numOfSubtanglesInput, int *gaussInput, int integerSubtangleConnectionsEM[][2][5], int *barsInput, int integerSubtangleParametersEM[][7], int *gaussIntegerSubtangleEM, int *barsGaussIntegerSubtangleEM){
	
	//Setup a local array to store integerSubtangleConnectionsEM[][2][5].
	//The array will be refined and stored in the local array, and finally replace the original array at the end.
	//Create a new variable to track the number of subtangles obtained by refining the EM code, initialized as zero.
	//Create a new variable to track the first "bar" placement (that is, the number of subtangles passed through before exiting the tangle the first time).
	int integerSubtangleConnectionsEMlocal[numOfSubtanglesInput][2][5];
	int subtanglesRefined=0;
	
	
	//In order to avoid incorrectly chaining across endpoints, we need to define a bunch arrays to track this information.
	//This is a lot. I'm not happy about it but we must account for this in the general case.
	bool furtherChainingPossible;
	
	//Initialize an array ot specify the subtangles connected at the endpoint arcs and the connected corners of these subtangles.
	//These are calculated from the gaussInput array and information about the barsInput.
	//The first column always refers to the first enter endpoint; column 2 refers to the first exit endpoint; colum 3 to the second enter endpoint; and column 4 to the second exit endpoint.
	int endpointConnectionsCorners[2][4];
	
	//The first crossing listed in the gaussInput array is always connected to the first enter endpoint arc, and always at corner A.
	//Likewise, the last crossing in gaussInput array is always conneccted to the second exit endpoint arc, and always at corner D.
	endpointConnectionsCorners[0][0]=abs(gaussInput[0]);
	endpointConnectionsCorners[1][0]=1;
	endpointConnectionsCorners[0][3]=abs(gaussInput[(2*numOfSubtanglesInput)-1]);
	endpointConnectionsCorners[1][3]=4;
	
	//The remaining two endpoints depend on whether this is the first or second time the adjacent subtangle has been encountered.
	int cornerExitIndex1=abs(gaussInput[barsInput[0]-1]);
	int cornerEnterIndex2=abs(gaussInput[barsInput[0]]);
	endpointConnectionsCorners[0][1]=cornerExitIndex1;
	endpointConnectionsCorners[0][2]=cornerEnterIndex2;
	
	//DEBUG:
	//printf("\n cornerExitIndex1 = %d \n cornerEnterIndex2 = %d \n", cornerExitIndex1,cornerEnterIndex2);
	//printf("\n endpointConnectionsCorners:\n\t S1-enter: %d \n\t S1-exit: %d \n\t S2-enter: %d \n\t s2-exit: %d \n", endpointConnectionsCorners[0][0],endpointConnectionsCorners[0][1],endpointConnectionsCorners[0][2],endpointConnectionsCorners[0][3]);
	
	int checkIndex;
	
	//Count the number of times cornerExitIndex1 shows up.
	checkIndex=0;
	for(int i=0; i<barsInput[0]; i++){
		if( abs(gaussInput[i]) == cornerExitIndex1 ){
			checkIndex++;
		}
	}
	//If it only showed up once, the first exit endpoint connects at corner B; if twice, it connects at corner D.
	if( checkIndex < 2 ){
		endpointConnectionsCorners[1][1]=2;
	} else {
		endpointConnectionsCorners[1][1]=4;
	}
	
	//Repeat this count and check for the number of times cornerEnterIndex2 shows up.
	checkIndex=0;
	for(int i=0; i<=barsInput[0]; i++){
		if( abs(gaussInput[i]) == cornerEnterIndex2 ){
			checkIndex++;
		}
	}
	//If it only showed up once, the second enter endpoint connects at corner A; if twice, it connects at corner C.
	if( checkIndex < 2 ){
		endpointConnectionsCorners[1][2]=1;
	} else {
		endpointConnectionsCorners[1][2]=3;
	}
	
	//DEBUG
	/*
	printf("\n endpointConnectionsCorners: Line 987\n\t [");
	for(int i=0; i<4; i++){
		printf(" %d ", endpointConnectionsCorners[0][i]);
	}
	printf("] \n\t [");
	for(int i=0; i<4; i++){
		printf(" %d ", endpointConnectionsCorners[1][i]);
	}
	printf("]\n");
	*/
		
	
	//While refining, we only want to look for "forward" chaining.
	//There are two types of forward chaining connections: b=c and b=d.
	//(This chaining corresponds to Twist Axis (AD-BC) and (AC-BD), respectively)
	//If we start a chain with one connection, we want to continue the chain with only this connection.
	//Track the chaining type with an integer chainType: 0 denotes no chain type, 1 denotes b=c, and -1 denotes b=d.
	//As more crossings are added to the chain, increment the magnitude of chainType.
	//Reset it to 0 for each new unchained crossing looked at.
	int chainType=0;
	
	//If (internal) vertical twisting happens while chaining crossings, the direction of S2 alternates based on the number of crossings chained.
	//This happens in four of the eight cases. Initialize a boolean variable to track this.
	bool alternatingS2direction;
	
	//For indexing purposes, each chain will share the un-chained indexing of the first crossing.
	//After constructing the refinement, we will go through and reindex by chains (that is, by integer subtangles).
	
	//Variable to store the index of the first crossing in a chain, used as the integer subtangle index; to be updated later.
	int currentIntegerSubtangleIndex;
	
	for(int i=0; i<numOfSubtanglesInput; i++){
		
		//Each iteration, increment the number of refined subtangles, and specify the current integer subtangle index.
		subtanglesRefined++;
		currentIntegerSubtangleIndex=i;
		
		for(int j=0; j<5; j++){
			integerSubtangleConnectionsEMlocal[subtanglesRefined-1][0][j]=integerSubtangleConnectionsEM[i][0][j];
			integerSubtangleConnectionsEMlocal[subtanglesRefined-1][1][j]=integerSubtangleConnectionsEM[i][1][j];
		}
		
		//Determine the crossing sign and track this in the integerSubtangleParametersEM[] array.
		integerSubtangleParametersEM[subtanglesRefined-1][0]=integerSubtangleConnectionsEM[i][1][0];
		
		//We need to include checks here to avoid a few special cases of unwated chaining being "detected" at endpoints.
		//If a crossing is chained to itself, or a crossing is incorrectly chained across an endpoint, we want to exclude this.
		//Note that the former case is accounted for in the second case, which can be elminated by checking against the barsInput array.
		//Thanks to forward chaining, it is sufficient to check only that the current crossing does not precede an endpoint.
		//To check this, we compare the current crossing index i+1 with the last crossing index (numOfSubtangles) and the crossing index before end of the first strand (gaussInput[barsInput[0]-1]).
		//UPDATE-- Above analysis is incorrect; check for endpoint chaining elsewhere.
		if( (integerSubtangleConnectionsEM[i][0][2] == integerSubtangleConnectionsEM[i][0][3]) && (i+1<numOfSubtangles) ){
			chainType=1;
		} else if ( (integerSubtangleConnectionsEM[i][0][2] == integerSubtangleConnectionsEM[i][0][4]) && (i+1<numOfSubtangles) ){
			chainType=-1;
		} else {
			//If neither of the above happens, no forward chaining is possible for this crossing.
			chainType=0;
		}
		
		//This is unfotunate, but the above check is not general enough. It turns out that checking against the endpoint connections takes a lot of work.
		//Check against the entries in the endpointConnectionsCorners array.
		if( chainType == 1 ){
			//In this case, B-C twisting occurs; make sure it doesn't happen over the endpoint. If it does, set the chainType to 0.
			for(int j=0; j<4; j++){
				if( (endpointConnectionsCorners[0][j]-1) == i ){
					if( (endpointConnectionsCorners[1][j]==2) || (endpointConnectionsCorners[1][j]==3) ){
						//If crossing i is adjacent to an endpoint, and if B-C twisting occurs, and if the endpoint corner is B or C, then this twisting is over an endpoint and not allowed.
						//Set the chainType to 0.
						chainType=0;
					}
				}
			}
		} else if( chainType == -1 ){
			//In this case, B-D twisting occurs; make sure it doesn't happen over the endpoint. If it does, set the chainType to 0.
			for(int j=0; j<4; j++){
				if( (endpointConnectionsCorners[0][j]-1) == i ){
					if( (endpointConnectionsCorners[1][j]==2) || (endpointConnectionsCorners[1][j]==4) ){
						//If crossing i is adjacent to an endpoint, and if B-D twisting occurs, and if the endpoint corner is B or D, then this twisting is over an endpoint and not allowed.
						//Set the chainType to 0.
						chainType=0;
					}
				}
			}
		}
		
		
		//It is also possible that and endpoint crossing gets chained into an integer subtangle. In this case, we need to update the endpointConnectionsCorners array to reflect this.
		//This only matters if there is chaining, and only if chaining an endpoint crossing. The updated index is the first crossing in the chain (currentIntegerSubtangleIndex in this case).
		//Since both chain types involve twisting corner B, the chained crossing index is the index of the crossing connected to this corner.
		//For integer subtangle chaining, the corner connected to the endpoint does NOT change (unless I'm mistaken, but future me can deal with that headache).
		//ERROR TO FIX--updating  the endpoint crossing indices here will throw off the later indexing, update endpoint labels while udpating all of the crossing indices?
		/*
		if( chainType != 0 ){
			for(int j=0; j<4; j++){
				if( (endpointConnectionsCorners[0][j]) == integerSubtangleConnectionsEM[i][0][2] ){
					endpointConnectionsCorners[0][j] = (currentIntegerSubtangleIndex+1);
					printf("\n corner index update at 1088");	
				}
			}
		}
		*/
		
		//DEBUG
		/*
		printf("\n endpointConnectionsCorners: Line 1088 Check \n\t [");
		for(int i=0; i<4; i++){
			printf(" %d ", endpointConnectionsCorners[0][i]);
		}
		printf("] \n\t [");
		for(int i=0; i<4; i++){
			printf(" %d ", endpointConnectionsCorners[1][i]);
		}
		printf("]\n");
		*/
		
		
		//Update integerSubtangleParametersEM[][1] to store the chainType, which denotes the twist axis.
		//Axis (AD-BC) corresponds to chain type 1; axis (AC-BD) corresponds to chain type -1; no chaining sets the axis entry to 0.
		integerSubtangleParametersEM[subtanglesRefined-1][1]=chainType;
		
		//Reset the boolean tracking whether the direction of S2 alternates at each iteration.
		alternatingS2direction=false;
		
		//Next, fill in as much of the integerSubtangleParmetersEM[][] array as possible. There are up to eight cases.
		//First, split into cases depending on crossing sign:
		if( integerSubtangleParametersEM[subtanglesRefined-1][0] > 0 ){
			//Poisitive Crossing:
			
			//Second, split into more cases depending on whether the first strands enters along the over or understrand (tracked in integerSubtangleConnections[i][0][0]):
			if( integerSubtangleConnectionsEM[i][0][0] > 0 ){
				//Enter along the overstrand:
				
				//S2 has its direction reversed:
				integerSubtangleParametersEM[subtanglesRefined-1][6]=-1;
				
				//Third, split into cases depending on the twist axis:
				if( chainType > 0 ){
					//Twist Axis: (AD-BC), corresponding to chainType 1:
					//Internal vertical (-1) twisting in this case:
					integerSubtangleParametersEM[subtanglesRefined-1][2]=-1;
					//In this case, the direction of S2 alternates:
					alternatingS2direction=true;
				} else if( chainType < 0 ){
					//Twist Axis: (AC-BD), corresponding to chainType -1:
					//Internal horizontal (1) twisting in this case:
					integerSubtangleParametersEM[subtanglesRefined-1][2]=1;				
				} else {
					//In this case, no chaining happens; this single crossing constitutes its own "chain".
					//It can later be joined along either "axis", so set the twist axis entry to be 0.
					integerSubtangleParametersEM[subtanglesRefined-1][2]=0;
				}
				
			} else {
				//Enter along the understrand:
				
				//S2 has its direction preserved:
				integerSubtangleParametersEM[subtanglesRefined-1][6]=1;
				
				//Third, split into cases depending on the twist axis:
				if( chainType > 0 ){
					//Twist Axis: (AD-BC), corresponding to chainType 1:
					//Internal horizontal (1) twisting in this case:
					integerSubtangleParametersEM[subtanglesRefined-1][2]=1;
				} else if( chainType < 0 ){
					//Twist Axis: (AC-BD), corresponding to chainType -1:
					//Internal vertical (-1) twisting in this case:
					integerSubtangleParametersEM[subtanglesRefined-1][2]=-1;
					//In this case, the direction of S2 alternates:
					alternatingS2direction=true;
				} else {
					//In this case, no chaining happens; this single crossing constitutes its own "chain".
					//It can later be joined along either "axis", so set the twist axis entry to be 0.
					integerSubtangleParametersEM[subtanglesRefined-1][2]=0;
				}
				
			}
		} else {
			//Negative Crossing:
			
			//Second, split into more cases depending on whether the first strands enters along the over or understrand (tracked in integerSubtangleConnections[i][0][0]):
			if( integerSubtangleConnectionsEM[i][0][0] > 0 ){
				//Enter along the overstrand:
				
				//S2 has its direction preserved:
				integerSubtangleParametersEM[subtanglesRefined-1][6]=1;
				
				//Third, split into cases depending on the twist axis:
				if( chainType > 0 ){
					//Twist Axis: (AD-BC), corresponding to chainType 1:
					//Internal horizontal (1) twisting in this case:
					integerSubtangleParametersEM[subtanglesRefined-1][2]=1;
				} else if( chainType < 0 ){
					//Twist Axis: (AC-BD), corresponding to chainType -1:
					//Internal vertical (-1) twisting in this case:
					integerSubtangleParametersEM[subtanglesRefined-1][2]=-1;
					//In this case, the direction of S2 alternates:
					alternatingS2direction=true;
				} else {
					//In this case, no chaining happens; this single crossing constitutes its own "chain".
					//It can later be joined along either "axis", so set the twist axis entry to be 0.
					integerSubtangleParametersEM[subtanglesRefined-1][2]=0;
				}
				
			} else {
				//Enter along the understrand:
				
				//S2 has its direction reversed:
				integerSubtangleParametersEM[subtanglesRefined-1][6]=-1;
				
				//Third, split into cases depending on the twist axis:
				if( chainType > 0 ){
					//Twist Axis: (AD-BC), corresponding to chainType 1:
					//Internal vertical (-1) twisting in this case:
					integerSubtangleParametersEM[subtanglesRefined-1][2]=-1;
					//In this case, the direction of S2 alternates:
					alternatingS2direction=true;
				} else if( chainType < 0 ){
					//Twist Axis: (AC-BD), corresponding to chainType -1:
					//Internal horizontal (1) twisting in this case:
					integerSubtangleParametersEM[subtanglesRefined-1][2]=1;
				} else {
					//In this case, no chaining happens; this single crossing constitutes its own "chain".
					//It can later be joined along either "axis", so set the twist axis entry to be 0.
					integerSubtangleParametersEM[subtanglesRefined-1][2]=0;
				}
			}
		}
		
		//At this point, all of the entries in integerSubtangleParametersEM have been updated except for p, q, and the internal parity.
		//These depend on the number of chained crossings and whether this is even or odd, and depends on the twist axis.
		//Wait to update these values until after chaining.
		
		//If chaining is possible, continue to check for possible crossings to add to the chain.
		if( chainType == 1 ){
			
			currentIntegerSubtangleIndex=i;
			
			//If the chaining condition is still satisfied, keep iterating.
			while( chainType!=0 ){
				
				//Move the EM code index to the next crossing:
				i++;
				
				//Combine the next crossing with the current chain.
				//Forward chaining with: b=c
				//Replace the "b" corner of the current chain with that of the next crossing:
				integerSubtangleConnectionsEMlocal[subtanglesRefined-1][0][2]=integerSubtangleConnectionsEM[i][0][2];
				integerSubtangleConnectionsEMlocal[subtanglesRefined-1][1][2]=integerSubtangleConnectionsEM[i][1][2];
				//Likewise, replace the "c" corner of the current chain with that of the next crossing:
				integerSubtangleConnectionsEMlocal[subtanglesRefined-1][0][3]=integerSubtangleConnectionsEM[i][0][3];
				integerSubtangleConnectionsEMlocal[subtanglesRefined-1][1][3]=integerSubtangleConnectionsEM[i][1][3];
				
				//If this chaining is vertical, the direction of S2 alternates each time a new crossing is added.
				if( alternatingS2direction ){
					integerSubtangleParametersEM[subtanglesRefined-1][6] = ( -1 * integerSubtangleParametersEM[subtanglesRefined-1][6] );
				}
				
				//Now check if the next crossing can be added to the chain with the same chain type (b=c).
				//Also check that we haven't run out of possible crossings to chain.
				//ALSO check that this would not cause the subtangle to chain across endpoints.
				if( (integerSubtangleConnectionsEM[i][0][2] == integerSubtangleConnectionsEM[i][0][3]) && (i+1<numOfSubtangles) ){
					
					//Set a boolean to indicate that further chaining is possible, but check that this doesn't chain over endpoints.
					furtherChainingPossible=true;
						
					//Again, the above check is not general enough to rule out endpoint twisting.
					//Check against the entries in the endpointConnectionsCorners array.
					if( chainType > 0 ){
						//In this case, B-C twisting occurs; make sure it doesn't happen over the endpoint.
						for(int j=0; j<4; j++){
							if( (endpointConnectionsCorners[0][j]-1) == i ){
								if( (endpointConnectionsCorners[1][j]==2) || (endpointConnectionsCorners[1][j]==3) ){
									//If crossing i is adjacent to an endpoint, and if B-C twisting occurs, and if the endpoint corner is B or C, then this twisting is over an endpoint and not allowed.
									//Set the endpoint chaining boolean to false.
									furtherChainingPossible=false;
								}
							}
						}
					}
					
					//It is also possible that an endpoint crossing gets chained into an integer subtangle. In this case, we need to update the endpointConnectionsCorners array to reflect this.
					//This only matters if there is chaining, and only if chaining an endpoint crossing. The updated index is the first crossing in the chain (i+1 in this case). ERROR--NEED TO DEBUG
					//Since both chain types involve twisting corner B, the chained crossing index is the index of the crossing connected to this corner.
					//For integer subtangle chaining, the corner connected to the endpoint does NOT change (unless I'm mistaken, but future me can deal with that headache).
					if( furtherChainingPossible==true ){
						/*
						for(int j=0; j<4; j++){
							if( (endpointConnectionsCorners[0][j]) == integerSubtangleConnectionsEM[i][0][2] ){
								endpointConnectionsCorners[0][j] = (currentIntegerSubtangleIndex+1);
								printf("\n corner index update at 1281");
							}
						}
						*/
						chainType++;
						
						//DEBUG
						/*
						printf("\n endpointConnectionsCorners: Line 1281 Check \n\t [");
						for(int i=0; i<4; i++){
							printf(" %d ", endpointConnectionsCorners[0][i]);
						}
						printf("] \n\t [");
						for(int i=0; i<4; i++){
							printf(" %d ", endpointConnectionsCorners[1][i]);
						}
						printf("]\n");
						*/
					}
					
				} else {
					furtherChainingPossible=false;
				}
				
				if( furtherChainingPossible == false ){
					//If no more crossings can be chained, we would still like to track the total number of chained crossings.
					//This will be tracked in the array integerSubtangleParametersEM[], along with the sign of the chained crossings.
					//Note that chained crossings must share the same sign, else they would be eliminated by an earlier R2 check.
					//THIS WILL NOT BE CORRECT IF USED WITH A TANGLE THAT CAN BE SIMPLIFIED IN THIS WAY--this assumes the tangle is already simplified in this regard.
					//Scale the array entry by the number of chained crossings.
					integerSubtangleParametersEM[subtanglesRefined-1][0] = ( abs(chainType) + 1 ) * integerSubtangleParametersEM[subtanglesRefined-1][0];
					
					//Then set the chainType to zero to terminate the while loop,
					chainType=0;
				}		
			}
			
		} else if( chainType == -1 ){
			
			//If the chaining condition is still satisfied, keep iterating.
			while( chainType!=0 ){
				
				//Move the EM code index to the next crossing:
				i++;
				
				//Combine the next crossing with the current chain.
				//Forward chaining with: b=d
				//Replace the "b" corner of the current chain with that of the next crossing:
				integerSubtangleConnectionsEMlocal[subtanglesRefined-1][0][2]=integerSubtangleConnectionsEM[i][0][2];
				integerSubtangleConnectionsEMlocal[subtanglesRefined-1][1][2]=integerSubtangleConnectionsEM[i][1][2];
				//Likewise, replace the "d" corner of the current chain with that of the next crossing:
				integerSubtangleConnectionsEMlocal[subtanglesRefined-1][0][4]=integerSubtangleConnectionsEM[i][0][4];
				integerSubtangleConnectionsEMlocal[subtanglesRefined-1][1][4]=integerSubtangleConnectionsEM[i][1][4];
				
				//If this chaining is vertical, the direction of S2 alternates each time a new crossing is added.
				if( alternatingS2direction ){
					integerSubtangleParametersEM[subtanglesRefined-1][6] = ( -1 * integerSubtangleParametersEM[subtanglesRefined-1][6] );
				}
				
				//Now check if the next crossing can be added to the chain with the same chain type (b=d).
				//Also check that we haven't run out of possible crossings to chain.
				//ALSO check that this would not cause the subtangle to chain across endpoints.
				if( (integerSubtangleConnectionsEM[i][0][2] == integerSubtangleConnectionsEM[i][0][4]) && (i+1<numOfSubtangles) ){
					
					//Set a boolean to indicate that further chaining is possible, but check that this doesn't chain over endpoints.
					furtherChainingPossible=true;
						
					//Again, the above check is not general enough to rule out endpoint twisting.
					//Check against the entries in the endpointConnectionsCorners array.
					if( chainType < 0 ){
						//In this case, B-D twisting occurs; make sure it doesn't happen over the endpoint.
						for(int j=0; j<4; j++){
							if( (endpointConnectionsCorners[0][j]-1) == i ){
								if( (endpointConnectionsCorners[1][j]==2) || (endpointConnectionsCorners[1][j]==4) ){
									//If crossing i is adjacent to an endpoint, and if B-D twisting occurs, and if the endpoint corner is B or D, then this twisting is over an endpoint and not allowed.
									//Set further chaining boolean to false.
									furtherChainingPossible=false;
								}
							}
						}
					}
					
					//It is also possible that and endpoint crossing gets chained into an integer subtangle. In this case, we need to update the endpointConnectionsCorners array to reflect this.
					//This only matters if there is chaining, and only if chaining an endpoint crossing. The updated index is the first crossing in the chain (i+1 in this case).
					//Since both chain types involve twisting corner B, the chained crossing index is the index of the crossing connected to this corner.
					//For integer subtangle chaining, the corner connected to the endpoint does NOT change (unless I'm mistaken, but future me can deal with that headache).
					//ERROR TO FIX
					if( furtherChainingPossible==true ){
						/*
						for(int j=0; j<4; j++){
							if( (endpointConnectionsCorners[0][j]-1) == integerSubtangleConnectionsEM[i][0][2] ){
								endpointConnectionsCorners[0][j] = (i+1);
							}
						}
						*/
						chainType--;
					}
					
					
				
				} else {
					furtherChainingPossible=false;
				}
				
				if( furtherChainingPossible == false ){
					//If no more crossings can be chained, we would still like to track the total number of chained crossings.
					//This will be tracked in the array integerSubtangleParametersEM[], along with the sign of the chained crossings.
					//Note that chained crossings must share the same sign, else they would be eliminated by an earlier R2 check.
					//THIS WILL NOT BE CORRECT IF USED WITH A TANGLE THAT CAN BE SIMPLIFIED IN THIS WAY--this assumes the tangle is already simplified in this regard.
					//Scale the array entry by the number of chained crossings.
					integerSubtangleParametersEM[subtanglesRefined-1][0] = ( abs(chainType) + 1 ) * integerSubtangleParametersEM[subtanglesRefined-1][0];
					
					//Then set the chainType to zero to terminate the while loop.
					chainType=0;
				}		
			}
		}
		
		
		//Returning to the integerSubtangleParametersEM[][] array, we specify the values of p and q for the corresponding rational fraction p/q.
		//If the integer subtangle has n crossings, this will either be 1/n or n/1, with the sign determined by the crossing sign AND whether or not S2 is reversed.
		//This depends on whether the internal twisting is considered horizontal (n/1) or vertical (1/n).
		if( integerSubtangleParametersEM[subtanglesRefined-1][2] > 0 ){
			//This case corresponds to horizontal twisting; the numerator p is set equal to the number of crossings, while the denominator q is set to 1.
			integerSubtangleParametersEM[subtanglesRefined-1][3] = abs(integerSubtangleParametersEM[subtanglesRefined-1][0]);
			integerSubtangleParametersEM[subtanglesRefined-1][4] = 1;
		} else {
			//This case corresponds to vertical twisting (or a 1-crossing chain); the numerator p is set to 1, while the denominator q is set equal to the number of crossings.
			integerSubtangleParametersEM[subtanglesRefined-1][3] = 1;
			integerSubtangleParametersEM[subtanglesRefined-1][4] = abs(integerSubtangleParametersEM[subtanglesRefined-1][0]);
		}
		//The overall sign of the fraction is tracked in q, by convention. To match the usual parity convention, this depends both on the internal sign of the crossings AND whether or not S2 is reversed.
		//If S2 is not reversed, the sign of q matches the sign of the crossings; otherwise, the sign of q is opposite the sign of the crossings.
		//Effectively, this means that q matches the sign of the product of integerSubtangleParametersEM[] entries 0 and 6.
		if( ( integerSubtangleParametersEM[subtanglesRefined-1][0] * integerSubtangleParametersEM[subtanglesRefined-1][6] ) < 0 ){
			integerSubtangleParametersEM[subtanglesRefined-1][4] = -1*(integerSubtangleParametersEM[subtanglesRefined-1][4]);
		}
		//SPECIAL CASE: 
		//When there are an even number of vertical twists, the sign must be changed one more time to obtain the correct sign of the shape of that subtangle (due to how products with parity 1 tangles change signs).
		//Note that this will only happen when q is divisible by 2 (as q is defined above, integerSubtangleParametersEM[4]).
		if( (integerSubtangleParametersEM[subtanglesRefined-1][4] % 2) == 0 ){
			integerSubtangleParametersEM[subtanglesRefined-1][4] = -1*(integerSubtangleParametersEM[subtanglesRefined-1][4]);
		}
		
		//Update the internal parity of the subtangle. If the number of crossings is odd, the parity is always 1.
		//If the number of crossings is even, the parity is either 0 or infinity(2) depending on the twist type (horizontal or vertical).
		if ( ( integerSubtangleParametersEM[subtanglesRefined-1][0] % 2) == 0 ){
			//If there are an even number of crossings, the subtangle has parity 0 if it twists horizontally and parity infinity(2) if it twists vertically.
			if( integerSubtangleParametersEM[subtanglesRefined-1][2] == 1 ){
				//Horzitonal twist type: partiy 0
				integerSubtangleParametersEM[subtanglesRefined-1][5]=0;
			} else if( integerSubtangleParametersEM[subtanglesRefined-1][2] == -1 ){
				//Vertical twist type: parity infinity(2)
				integerSubtangleParametersEM[subtanglesRefined-1][5]=2;
			}
		} else {
			//If the number of crossing is odd, the subtangle parity is 1.
			integerSubtangleParametersEM[subtanglesRefined-1][5]=1;
		}
			
	}
	
	//DEBUG
	/*
	printf("\n After chaining crossings; integer subtangles, but not reindexed:");
	printArrayEM(subtanglesRefined,integerSubtangleConnectionsEMlocal);
	*/
			
	//At this point, the chaining is completed, but the chains are indexed by the Gauss code index of the first crossing in the chain.
	//Create an array to track these indices.
	int oldIntegerSubtangleIndex[subtanglesRefined];
	oldIntegerSubtangleIndex[0]=1;
	for (int i=1; i<subtanglesRefined; i++){
		//To determine the next old index, shift over the previous old index by the number of crossings chained at that step.
		oldIntegerSubtangleIndex[i] = ( oldIntegerSubtangleIndex[i-1] + abs(integerSubtangleParametersEM[i-1][0]) );
	}
	//DEBUG
	/*
	printf("\n oldIntegerSubtangleIndex:\n [");
	for(int i=0; i<subtanglesRefined; i++){
		printf(" %d ", oldIntegerSubtangleIndex[i]);
	}
	printf("] \n");
	printf("\n Old Bars: [ %d %d ]\n", barsInput[0],barsInput[1]);
	*/
	
	//We would like to reindex this to be increasing on the chains, rather than the by the old Guass index of the first crossing in a chain.
	//After chaining, the only surviving (old) indices of a chain reference the first and last crossing of the "links" of the chain.
	//Array entries which reference the old index of the last link of a chain needs to be set to match the first link.
	int currentOldIndex;
	int nextOldIndex;
	for(int i=0; i<subtanglesRefined; i++){
		//In principle, we will "collapse" non-matching old indices along each chain to the index of the first crossing.
		//The the currentOldIndex is the first index of the current chain.
		currentOldIndex=oldIntegerSubtangleIndex[i];
		//The nextOldIndex is one less than the first index of the next chain (equivalently, this is the old Gauss index of the last crossing in the chain, which connected chains will reference).
		//As an artifact of the chaining construction, non-matching indices in a chain can only take this value
		if( (i+1) < subtanglesRefined ){
			nextOldIndex = (oldIntegerSubtangleIndex[i+1] - 1 );
		} else {
			//In the special case where we consider the last chain, the nextOldIndex is the original numOfSubtanglesInput
			nextOldIndex = numOfSubtanglesInput;
		}
		//Next, check to see if any references are made to nextOldIndex; if they are, set this reference to currentOldIndex.
		//Effectively, this will make it so that every reference to a chain shares the same index, in this case the original index of the first crossing in the chain.
		for(int j=0; j<subtanglesRefined; j++){	
			for(int k=1; k<5; k++){
				if( integerSubtangleConnectionsEMlocal[j][0][k] == nextOldIndex ){	
					integerSubtangleConnectionsEMlocal[j][0][k] = currentOldIndex;
				}	
			}	
		}
		//Repeat this check for references to endpoint corners; if references are made to nextOldIndex, update these to currentOldIndex.
		for(int j=0; j<4; j++){
			if( endpointConnectionsCorners[0][j] == nextOldIndex ){
				endpointConnectionsCorners[0][j] = currentOldIndex;
			}
		}
		//DEBUG
		/*
		printf("\n endpointConnectionsCorners: Line 1497\n\t [");
		for(int i=0; i<4; i++){
			printf(" %d ", endpointConnectionsCorners[0][i]);
		}
		printf("] \n\t [");
		for(int i=0; i<4; i++){
			printf(" %d ", endpointConnectionsCorners[1][i]);
		}
		printf("]\n");
		*/
		
	}
	//DEBUG
	/*
	printf("\n After matching chain indices, but before re-labeling the indexing:");
	printArrayEM(subtanglesRefined,integerSubtangleConnectionsEMlocal);
	*/
	
	//After the indices of each chain match, we re-label so that indexing is increasing on the chains.
	for(int i=0; i<subtanglesRefined; i++){
		for( int j=1; j<5; j++){
			for( int k=0; k<subtanglesRefined; k++){
				if( integerSubtangleConnectionsEMlocal[i][0][j] == oldIntegerSubtangleIndex[k] ){
					integerSubtangleConnectionsEMlocal[i][0][j] = (k+1);
				}
			}
		}
	}
	//Also update the references to the endpoint corners.
	for(int j=0; j<4; j++){
		for(int k=0; k<subtanglesRefined; k++){
			if( endpointConnectionsCorners[0][j] == oldIntegerSubtangleIndex[k] ){
				endpointConnectionsCorners[0][j] = (k+1);
			}
		}
	}
	//DEBUG
	/*
	printf("\n endpointConnectionsCorners: Line 1535\n\t [");
	for(int i=0; i<4; i++){
		printf(" %d ", endpointConnectionsCorners[0][i]);
	}
	printf("] \n\t [");
	for(int i=0; i<4; i++){
		printf(" %d ", endpointConnectionsCorners[1][i]);
	}
	printf("]\n");
	*/
	
	//Originally, entry [i][1][0] in the EM array referred to the sign of the crossing encountered.
	//Modify this entry in the integer subtangle EM array so that its magnitude is the number of crossings in the subtangle, and it's sign denotes the sign of the crossings.
	//This information is currently stored in the integerSubtangleParametersEM array.
	for(int i=0; i<subtanglesRefined; i++){
		integerSubtangleConnectionsEMlocal[i][1][0]=integerSubtangleParametersEM[i][0];
	}
	
	//DEBUG
	/*
	printf("\n integerSubtangleParametersEM: [");
	for(int i=0; i<subtanglesRefined; i++){
		printf(" %d ", integerSubtangleParametersEM[i]);
	}
	printf("]\n");
	printf("\n After fully re-indexing and re-labeling the EM code on integer subtangles (chains),\n and tracking the number and sign of the crossings per subtangle:\n");
	printArrayEM(subtanglesRefined,integerSubtangleConnectionsEMlocal);
	*/
	
	
	//It would helpful to have an analogue of the Guass code for integer subtangles so that one can see how the subtangles are encountered at a glance.
	//Note that this loses a lot of information about the actual subtangles, and over/under-strand information does not make sense in this context.
	//Initialize a (poorly named) check array to track when we encounter each subtangle; entries are intially set to 0, but changed to 1 after encountering the corresponding subtangle.
	//This is used to determine wh NOT NEEDED? DELETE?
	int gaussIntegerSubtangleEMcheckArray[subtanglesRefined] = {0};
	
	gaussIntegerSubtangleEM[0]=1;
	
	for(int i=1; i<2*subtanglesRefined; i++){
		//The next subtangle encountered is determined by the exit strand of the previous subtangle encountered, of which there are two possibilities.
		if( gaussIntegerSubtangleEMcheckArray[ gaussIntegerSubtangleEM[i-1] - 1 ] == 0 ){
			//If the previous subtangle has only been encountered once so far, the next subtangle is connected to its S1-exit (row 0, column 2)
			gaussIntegerSubtangleEM[i] = integerSubtangleConnectionsEMlocal[ gaussIntegerSubtangleEM[i-1]-1 ][0][2];
			//Update the previous subtangle to indicate it has been encountered.
			gaussIntegerSubtangleEMcheckArray[ gaussIntegerSubtangleEM[i-1] - 1 ]++;
		} else {
			//If the previous subtangle has been encountered twice now, the next subtangle is connected to its S2-exit (row 0, column 4).
			gaussIntegerSubtangleEM[i] = integerSubtangleConnectionsEMlocal[ gaussIntegerSubtangleEM[i-1]-1 ][0][4];
		}
	}
	
	//New bar placement, looking at integer subtangles rather than crossings:
	int newFirstBarIndex=0;
	for(int i=0; i<barsInput[0]; i++){
		newFirstBarIndex++;
		//This super convoluted thing below is attempting to track "skipped" crossings while determining after which subtangle the first tangle strand exits.
		//Yeah it's a mess, but it works! Probably. (confirm with some wierd looking tangle)
		i = i + ( abs( integerSubtangleParametersEM[gaussIntegerSubtangleEM[ newFirstBarIndex-1 ]-1][0] ) - 1 );
	}
	
	barsGaussIntegerSubtangleEM[0]=newFirstBarIndex;
	barsGaussIntegerSubtangleEM[1]=2*subtanglesRefined;
	
	//DEBUG
	/*
	printf("\n Gauss code analog for integer subtangles:\n [");
	for(int i=0; i<2*subtanglesRefined; i++){
		printf(" %d ", gaussIntegerSubtangleEM[i]);
	}
	printf("]\n");
	printf("\n New Bars: [ %d %d ]\n", barsGaussIntegerSubtangleEM[0],barsGaussIntegerSubtangleEM[1]);
	*/
	
	//Finally, overwrite the input integerSubtangleConnectionsEM[][2][5] array with the local array storing the chained crossings.
	numOfSubtanglesInput=subtanglesRefined;
	for(int i=0; i<subtanglesRefined; i++){
		for(int j=0; j<5; j++){
			integerSubtangleConnectionsEM[i][0][j]=integerSubtangleConnectionsEMlocal[i][0][j];
			integerSubtangleConnectionsEM[i][1][j]=integerSubtangleConnectionsEMlocal[i][1][j];
		}
	}
	
}



/*
(NC, 9/9/18)
Prototype crossing-chaining function--used to build up integer subtangles.
This is a work in progress, rewrite the description after I figure out how it's supposed to work.

Inputs:

Outputs (in the form of modified global variables):
numOfSubtangles
integerSubtangleConnectionsEM[][2][5]
integerSubtangleParametersEM[][7]
gaussIntegerSubtangleEM[]
barsGaussIntegerSubtangleEM[2]


This function calls:
	reindexEMforSubtangles()
	chainCrossings()
	printArrayEM() (DEBUG)
	printIntegerSubtangleParametersEM() (DEBUG)
This function is called by:
	buildGeneralizedEMCode()

*/

void buildIntegerSubtanlgeEMCode(int numOfCrossingsInput, int *gaussInput, int connectionsEM[][2][5], int *barsInput, int &numOfSubtangles, int integerSubtangleConnectionsEM[][2][5], int integerSubtangleParametersEM[][7], int *gaussIntegerSubtangleEM, int *barsGaussIntegerSubtangleEM){
	
	//Initialize the subtangle EM code to be the same as the usual EM code.
	//As the code is "refined" by chaining crossings when possible, this subtangle EM code will be updated.
	//This way, the original EM code of the tangle still exists if needed.
	numOfSubtangles=numOfCrossingsInput;
	reindexEMforSubtangles(numOfCrossingsInput,connectionsEM,integerSubtangleConnectionsEM);
	
	//The majority of this construction involves chaining crossings, if possible. This is checked in the following function.
	chainCrossings(numOfSubtangles,gaussInput,integerSubtangleConnectionsEM,barsInput,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	//DEBUG:
	printIntegerSubtangleParametersEM(numOfSubtangles,integerSubtangleParametersEM);
	
}




/*
(NC, 9/26/18)
Function to detect twisting along the specified corners of the current rational amalgam.
(WRITE A DETAILED DESCRIPTION AFTER IT WORKS)

If twisting does occur, return the integer number of twists that happen (SIGNED BY SLOPE?)
If twisting does not occur, return 0.


Inputs:

Outputs:


This function calls:

This function is called by:
	buildGeneralizedEMCode()

*/
int detectCornerTwisting(int subtangleIndex, int cornerIndex1, int cornerIndex2, int rationalSubtangleConnectionsEM[][2][5], int *integerSubtangleJoins, int integerSubtangleParameters[][7], int *gaussSubtangleInput, int *barsInput){
	
	//Run through a series of checks to determine whether there is twisting at the two input corners.
	//If any check fails, the function will skip the remaining checks and return false.
	//If all chceks succeed, the function will return true and end.
	
	//Initialize a variable to represent the index of the integer subtangle that is possibly being joined, and variables for the corners possibly joined at.
	//Note that possibleJoinedIndex starts at 1, not 0; remember to subtract 1 to obtain the corresponding array entry.
	int possibleJoinedIndex = rationalSubtangleConnectionsEM[subtangleIndex][0][cornerIndex1];
	int possibleJoinedCorner1 = rationalSubtangleConnectionsEM[subtangleIndex][1][cornerIndex1];
	int possibleJoinedCorner2 = rationalSubtangleConnectionsEM[subtangleIndex][1][cornerIndex2];
	
	//The number and type of twists, if they occur, is specified by the values of p and q in the rational fraction for the possible joined integer subtangle.
	//For an integer subtangle, either p is 1 or q is 1, and the sign of the twists is tracked in q by convention. Hence, the signed number of twists is exactly the product p*q.
	//Initialize this value here (it always exists by construction); if all of the checks succeed, the function returns this number.
	int pq = (integerSubtangleParameters[possibleJoinedIndex-1][3] * integerSubtangleParameters[possibleJoinedIndex-1][4]);
	
	
	//CHECK BELOW IS NOW OBSELTED BY A BETTER CHECK ELSEWHERE, REMOVE IF CERTAIN?
	bool endpointJoin=false;
	/*
	//What follows is an unfortunately long and convoluted bit of code needed to eliminate a special case where twisting is incorrectly "detected" across tangle endpoints.
	//There are a number of very specific ways this could happen, but we don't want exclude legal twisting, so the careful checking requires a fair bit of work.
	
	//It is possible that twisting may be incorrectly identified across endpoints. This is an unusual case that we want to rule out.
	//To check against this, intialize variables for the index of the subtangle adjacent to each endpoint, based on where the two strands enter/exit.
	//This information is extracted from the gaussEM array with integer subtangles and the corresponding bar information.
	int endpoint1SubtangleIndex = gaussSubtangleInput[0];
	int endpoint2SubtangleIndex = gaussSubtangleInput[barsInput[0]-1];
	int endpoint3SubtangleIndex = gaussSubtangleInput[barsInput[0]];
	int endpoint4SubtangleIndex = gaussSubtangleInput[barsInput[1]-1];
	//However, it's not enough just to know that subtangles are connected across an endpoint; it is possible that they are ALSO legally joined in a dtectable twist.
	//In this case, we also need to compare against the joined corners. The first and fourth endpoint subtangles always join at corners A and D, repsectively.
	//The corners of the second and third endpoint subtangles depend on whether this is the first or second encounter of the subtangles. This information must be inferred from the gaussEM array.
	int endpoint1SubtangleCorner = 1;
	int endpoint4SubtangleCorner = 4;
	int endpoint2SubtangleCorner;
	int endpoint3SubtangleCorner;
	
	//Count how many times the endpoint2SubtangleIndex subtangle is encountered while moving along the first strand.
	int trackSubtangleEncounter=0;
	for(int i=0; i<barsInput[0]; i++){
		if( gaussSubtangleInput[i] == endpoint2SubtangleIndex ){
			trackSubtangleEncounter++;
		}
	}
	if( trackSubtangleEncounter == 1 ){
		//If this subtangle was encountered only once along the first strand, then corner B joins the endpoint.
		endpoint2SubtangleCorner=2;
	} else {
		//If the subtangle was encountered twice along the first strand, then corner D joins the endpoint.
		endpoint2SubtangleCorner=4;
	}
	//Now repeat this check to count how many times the endpoint3SubtangleIndex subtangle is encountered while moving along the first strand.
	//In this case, this subtangle is the first encountered on the second strand; hence, it could be encountered once or not at all on the first strand.
	trackSubtangleEncounter=0;
	for(int i=0; i<barsInput[0]; i++){
		if( gaussSubtangleInput[i] == endpoint3SubtangleIndex ){
			trackSubtangleEncounter++;
		}
	}
	if( trackSubtangleEncounter == 0 ){
		//If this subtangle was not encountered along the first strand, then corner A joins the endpoint.
		endpoint3SubtangleCorner=1;
	} else {
		//If the subtangle was encountered once along the first strand, then corner C joins the endpoint.
		endpoint3SubtangleCorner=3;
	}
	
	//Depending on which is the current subtangle index and which is the possible joined index, there are a handful of cases that a twist across an endpoint might be identified.
	//Intialize a boolean as false to track this; set it to true if any of these happens; this must be false when checked later in this function.
	bool endpointJoin=false;
	if( ((subtangleIndex+1)==endpoint1SubtangleIndex) && (possibleJoinedIndex==endpoint4SubtangleIndex) ){
		//endpoint1 connects endpoint4:
		if( ((cornerIndex1==endpoint1SubtangleCorner) && (possibleJoinedCorner1==endpoint4SubtangleCorner)) || ((cornerIndex2==endpoint1SubtangleCorner) && (possibleJoinedCorner2==endpoint4SubtangleCorner)) ){
			//Check if either of the possible joined connections goes across the endpoint.
			endpointJoin=true;
		}
	} else if( ((subtangleIndex+1)==endpoint4SubtangleIndex) && (possibleJoinedIndex==endpoint1SubtangleIndex) ){
		//endpoint4 connects endpoint1:
		if( ((cornerIndex1==endpoint4SubtangleCorner) && (possibleJoinedCorner1==endpoint1SubtangleCorner)) || ((cornerIndex2==endpoint4SubtangleCorner) && (possibleJoinedCorner2==endpoint1SubtangleCorner)) ){
			//Check if either of the possible joined connections goes across the endpoint.
			endpointJoin=true;
		}
	} else if( ((subtangleIndex+1)==endpoint2SubtangleIndex) && (possibleJoinedIndex==endpoint3SubtangleIndex) ){
		//endpoint2 connects endpoint3:
		if( ((cornerIndex1==endpoint2SubtangleCorner) && (possibleJoinedCorner1==endpoint3SubtangleCorner)) || ((cornerIndex2==endpoint2SubtangleCorner) && (possibleJoinedCorner2==endpoint3SubtangleCorner)) ){
			//Check if either of the possible joined connections goes across the endpoint.
			endpointJoin=true;
		}
	} else if( ((subtangleIndex+1)==endpoint3SubtangleIndex) && (possibleJoinedIndex==endpoint2SubtangleIndex) ){
		//endpoint3 connects endpoint2:
		if( ((cornerIndex1==endpoint3SubtangleCorner) && (possibleJoinedCorner1==endpoint2SubtangleCorner)) || ((cornerIndex2==endpoint3SubtangleCorner) && (possibleJoinedCorner2==endpoint2SubtangleCorner)) ){
			//Check if either of the possible joined connections goes across the endpoint.
			endpointJoin=true;
		}
	}
	
	//DEBUG
	//printf("\n\t Endpoint Subtangle Indices: %d , %d , %d , %d \n",endpoint1SubtangleIndex,endpoint2SubtangleIndex,endpoint3SubtangleIndex,endpoint4SubtangleIndex);
	//printf("\n DEBUG FLAG: corner1 = %d , corner2 = %d \n", cornerIndex1,cornerIndex2);
	*/
	
	
	
	//CHECK 1: Check if the two input corner indices connect to the same subtangle.
	if( rationalSubtangleConnectionsEM[subtangleIndex][0][cornerIndex1] == rationalSubtangleConnectionsEM[subtangleIndex][0][cornerIndex2] ){
		//CHECK 2: Check that these corners do not join the current tangle to itself, nor across the endpoints.
		if( ((possibleJoinedIndex-1) != subtangleIndex) && (endpointJoin == false) ){
			//CHECK 3: Check if the connected subtangle is an integer subtangle that hasn't already been joined to some other rational amalgam.
			if( integerSubtangleJoins[possibleJoinedIndex-1] == 1 ){
				//CHECK 4: Check if this proposed integer subtangle has been joined along it's twist axis.
				//There are two cases here, depending the on the twist axis type of the connected integer subtangle.
				//Note that a type 0 twist axis could be vertical or horizontal (this corresponds to a one crossing subtangle).
				if( integerSubtangleParameters[possibleJoinedIndex-1][1] > -1 ){
					//Twist Axis: (AD-BC). The possible joined corners have different indices, but these must be either (AD) or (BC) for twisting to happen.
					//This holds if and only if the sum of these is 5 (A+D=1+4=5, or B+C=2+3=5).
					if( (possibleJoinedCorner1+possibleJoinedCorner2) == 5 ){
						//Twisting happens, and the number and type of twists is pq:
						return abs(pq);
						//SIGNED?
					}
				} 
				//The previous check is NOT mutually exclusive with this one in the case of 1-crossing subtangles, even though the actually twist axis (if it exists) is unique in general.
				if( integerSubtangleParameters[possibleJoinedIndex-1][1] < 1 ){
					//Twist Axis: (AC-BD). The possible joined corners have different indices, but these must be either (AC) or (BD) for twisting to happen.
					//This holds if and only if the sum of these is even (A+C=1+3=4, or B+D=2+4=6).
					if( ( (possibleJoinedCorner1+possibleJoinedCorner2) % 2 ) == 0 ){
						//Twisting happens, and the number and type of twists is pq:
						return abs(pq);
						//SIGNED?
					}
				}
			}
		}
	}
	
	return 0;
	
}



/*
(NC, 10/8/18)
This is a small function used to determine which endpoint (A,B,C, or D) is connected to an input endpoint.
It takes an integer as input (1, 2, 3, or 4, respectively), and returns the corresponding integer of the connected endpoint.

Inputs:
inputEndpoint: an integer corresponding to the endpoint A, B, C, or D

Outputs:
connectedEnpoint: the integer corresponding to the connected endpoint

This function calls:
	N/A
This function is called by:
	joinIntegerToRational()

*/

int determineConnectedEndpoint(int inputEndpoint){
	
	if( inputEndpoint == 1 ){
		//A -> B
		return 2;
	} else if( inputEndpoint == 2 ){
		//B -> A
		return 1;
	} else if( inputEndpoint == 3 ){
		//C -> D
		return 4;
	} else if( inputEndpoint == 4 ){
		//D -> C
		return 3;
	} else {
		//Error case, this will only happen if the function input is incorrect.
		return 0;
	}
	
}



/*
(NC, 10/9/18)
This is a small function used to update the index references in the rationalSubtangleConnectionsEM array.
The (initial) index of a rational amalgam is the integer subtangle index of the first integer subtangle in the amalgam.
After joining an integer subtangle to the amalgam, existing references to the index of the joined subtangle must be modified to match the index of the amalgam.
Later, after all joinings have been completed, the labeling of the rational subtangles will ned to be update.


Inputs:
numOfSubtanglesInput: the number of subtangles in the array; this is what we iterate through.
initialIndex: the integer subtangle index of the initial integer subtangle in the rational amalgam.
joinedIndex: the integer subtangle index of the integer subtangle being joined to them amalgam (to be overwritten).
rationalSubtangleConnectionsEM[][2][5]: the multi-array of connection information for the rational subtangle generalized EM code.

Outputs (in the form of modified global variables):
rationalSubtangleConnectionsEM[][2][5]: any reference to a connection to joinedIndex is now made to initialIndex.


This function calls:
	N/A
This function is called by:
	joinIntegerToRational()

*/

void updateIntegerJoinToRationalIndex(int numOfSubtanglesInput, int initialIndex, int joinedIndex, int rationalSubtangleConnectionsEM[][2][5]){
	
	//For each subtangle in the rationalSubtanlgeConnectionsEM array, look at all of the connected subtangle indices (row 0, column entries 1 through 4).
	//If any of these matches the index of the joined subtangle, replaced by the index of the initial subtangle.
	//Note that subtangle index begins at 0 in the array, but the "displayed" index begins at 1; hence, we add 1 to the index in the check to account for this.
	for(int i=0; i<numOfSubtanglesInput; i++){
		for(int j=1; j<5; j++){
			if( rationalSubtangleConnectionsEM[i][0][j] == (joinedIndex+1) ){
				rationalSubtangleConnectionsEM[i][0][j] = initialIndex+1;
			}
		}
	}
	
}



/*
(NC, 10/8/18)
Prototype of a function to join an integer subtangle into a rational amalgam.
This function presupposes that possible twists have already been correctly detected, and it just worries about the joining.
This is a work in progress, I will rewrite the description after I figure out how to make it work.


Inputs:
rightTwist
botTwist
leftTwist
topTwist
&cornerNW
&cornerNE
&cornerSE
&cornerSW
&twistCount
rationalTwistVectorArray[]
integerSubtangleJoins[]
initialIndex
rationalSubtangleConnectionsEM[][2][5]
rationalSubtangleParametersEM[][8]
integerSubtangleParametersEM[][7]
&rationalSubtanglesRefined
numOfIntegerSubtangles
endpointConnectionsCorners[2][4]
gaussRationalSubtangleEM[]
barsGaussRationalSubtangleEM[2]


This function calls:
	determineConnectedEndpoint()
	updateIntegerJoinToRationalIndex()
This function is called by:
	buildRationalSubtangleEMCode()

*/
void joinIntegerToRational(int rightTwist, int botTwist, int leftTwist, int topTwist, int &cornerNW, int &cornerNE, int &cornerSE, int &cornerSW, int &twistCount, int rationalTwistVectorArray[][31], int *integerSubtangleJoins, int initialIndex, int rationalSubtangleConnectionsEM[][2][5], int rationalSubtangleParametersEM[][10], int integerSubtangleParametersEM[][7], int &rationalSubtanglesRefined, int numOfIntegerSubtangles, int endpointConnectionsCorners[2][4], int *gaussRationalSubtangleEM, int *barsGaussRationalSubtangleEM){
	
	//Initialize a variable to represent the index of the possible joined integer subtangle. The value of this variable will depend on the case.
	int joinedIntegerSubtangle1Index;
	int joinedIntegerSubtangle2Index;
	int joinedSubtangle1Corner1;
	int joinedSubtangle1Corner2;
	int joinedSubtangle2Corner1;
	int joinedSubtangle2Corner2;
	int localTwistCount;
	
	//Initialize a variable to store the current number of refined rational subtangles.
	//This is used when removing entries from the rational Gauss EM array.
	int oldRationalSubtanglesRefined=rationalSubtanglesRefined;
	
	//Ok, my apologies to future generations of mathematicians who may need to decipher this, but the labeling I'm about it introduce is terrible.
	//In essence, I need variables to track the index AND corresponding corner of the "subtangle connected to the other end of the newly connected strand of the joined integer subtangle".
	//In the rationalSubtangleConnectionsEM array, previous references to the newly joined subtangle must be updated to reflecte the initial index of the rational amalgam.
	//The corners of the joined subtangle might not be the same as those of the amalgam; so we need to update not only the index, but the corner information as well.
	//These variables track this information and are used to update the connection information when the joinings are made in the rationalSubtangleConnectionsEM array.
	//There are eight such variables in total, two (an index and a corner) for each possible joined corner.
	int joinedSubtangle1Corner1ConnectedCrossIndex;
	int joinedSubtangle1Corner1ConnectedCrossIndexCorner;
	int joinedSubtangle1Corner2ConnectedCrossIndex;
	int joinedSubtangle1Corner2ConnectedCrossIndexCorner;
	int joinedSubtangle2Corner1ConnectedCrossIndex;
	int joinedSubtangle2Corner1ConnectedCrossIndexCorner;
	int joinedSubtangle2Corner2ConnectedCrossIndex;
	int joinedSubtangle2Corner2ConnectedCrossIndexCorner;
	
	//Initialize variables to represent twisted corners, which start as the values of the input corners.
	//If twisting happens along any of these corners, the index may need to be updated.
	//At the end of the function, replace the original corners with these updated ones.
	int newCornerNW = cornerNW;
	int newCornerNE = cornerNE;
	int newCornerSE = cornerSE;
	int newCornerSW = cornerSW;
	
	if( (rightTwist>0) && (leftTwist>0) ){
		//CASE 1: two subtangles joined RIGHT AND LEFT
		localTwistCount = rightTwist+leftTwist;
		rationalSubtangleParametersEM[initialIndex][0] = rationalSubtangleParametersEM[initialIndex][0]+localTwistCount;
		rationalSubtangleParametersEM[initialIndex][1] = rationalSubtangleParametersEM[initialIndex][1]+2;
		integerSubtangleJoins[initialIndex] = integerSubtangleJoins[initialIndex]+2;
		rationalSubtanglesRefined = rationalSubtanglesRefined-2;
		
		//Determine the index of the two joined subtangles, and determine the corners of these tangles.
		//Update the corresponding entries in the integerSubtangleJoins array to indicate these have been joined.
		joinedIntegerSubtangle1Index = rationalSubtangleConnectionsEM[initialIndex][0][cornerNE]-1;
		joinedIntegerSubtangle2Index = rationalSubtangleConnectionsEM[initialIndex][0][cornerNW]-1;
		integerSubtangleJoins[joinedIntegerSubtangle1Index]--;
		integerSubtangleJoins[joinedIntegerSubtangle2Index]--;
		
		//Right twist corners:
		joinedSubtangle1Corner1 = rationalSubtangleConnectionsEM[initialIndex][1][cornerNE];
		joinedSubtangle1Corner2 = rationalSubtangleConnectionsEM[initialIndex][1][cornerSE];
		//Left twist corners:
		joinedSubtangle2Corner1 = rationalSubtangleConnectionsEM[initialIndex][1][cornerNW];
		joinedSubtangle2Corner2 = rationalSubtangleConnectionsEM[initialIndex][1][cornerSW];
		
		
		//Connected-connected subtangle information (I'm sorry):
		//RIGHT:
		joinedSubtangle1Corner1ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][0][determineConnectedEndpoint(joinedSubtangle1Corner1)]-1;
		joinedSubtangle1Corner1ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][1][determineConnectedEndpoint(joinedSubtangle1Corner1)];
		joinedSubtangle1Corner2ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][0][determineConnectedEndpoint(joinedSubtangle1Corner2)]-1;
		joinedSubtangle1Corner2ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][1][determineConnectedEndpoint(joinedSubtangle1Corner2)];
		//LEFT:
		joinedSubtangle2Corner1ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle2Index][0][determineConnectedEndpoint(joinedSubtangle2Corner1)]-1;
		joinedSubtangle2Corner1ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle2Index][1][determineConnectedEndpoint(joinedSubtangle2Corner1)];
		joinedSubtangle2Corner2ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle2Index][0][determineConnectedEndpoint(joinedSubtangle2Corner2)]-1;
		joinedSubtangle2Corner2ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle2Index][1][determineConnectedEndpoint(joinedSubtangle2Corner2)];
		
		//New connections after joining: overwrite the twisted corners with the new connections
		//RIGHT: join to NE and SE (remember to +1 since the arrays start at 0 internally)
		rationalSubtangleConnectionsEM[initialIndex][0][cornerNE] = joinedSubtangle1Corner1ConnectedCrossIndex+1;
		rationalSubtangleConnectionsEM[initialIndex][0][cornerSE] = joinedSubtangle1Corner2ConnectedCrossIndex+1;
		//LEFT: join to NW and SW (remember to +1 since the arrays start at 0 internally)
		rationalSubtangleConnectionsEM[initialIndex][0][cornerNW] = joinedSubtangle2Corner1ConnectedCrossIndex+1;
		rationalSubtangleConnectionsEM[initialIndex][0][cornerSW] = joinedSubtangle2Corner2ConnectedCrossIndex+1;
		
		//Update the connected corner information for the rational amalgam.
		//RIGHT:
		rationalSubtangleConnectionsEM[initialIndex][1][cornerNE] = joinedSubtangle1Corner1ConnectedCrossIndexCorner;
		rationalSubtangleConnectionsEM[initialIndex][1][cornerSE] = joinedSubtangle1Corner2ConnectedCrossIndexCorner;
		//LEFT:
		rationalSubtangleConnectionsEM[initialIndex][1][cornerNW] = joinedSubtangle2Corner1ConnectedCrossIndexCorner;
		rationalSubtangleConnectionsEM[initialIndex][1][cornerSW] = joinedSubtangle2Corner2ConnectedCrossIndexCorner;
		
		//Set the corner connection information of the connected-connected subtangle to match that of the rational amalgam.
		//RIGHT:
		rationalSubtangleConnectionsEM[joinedSubtangle1Corner1ConnectedCrossIndex][1][joinedSubtangle1Corner1ConnectedCrossIndexCorner] = cornerNE;
		rationalSubtangleConnectionsEM[joinedSubtangle1Corner2ConnectedCrossIndex][1][joinedSubtangle1Corner2ConnectedCrossIndexCorner] = cornerSE;
		//LEFT:
		rationalSubtangleConnectionsEM[joinedSubtangle2Corner1ConnectedCrossIndex][1][joinedSubtangle2Corner1ConnectedCrossIndexCorner] = cornerNW;
		rationalSubtangleConnectionsEM[joinedSubtangle2Corner2ConnectedCrossIndex][1][joinedSubtangle2Corner2ConnectedCrossIndexCorner] = cornerSW;
		
		
		//Update references to joined subtangle index:
		//RIGHT:
		updateIntegerJoinToRationalIndex(numOfIntegerSubtangles,initialIndex,joinedIntegerSubtangle1Index,rationalSubtangleConnectionsEM);
		//LEFT:
		updateIntegerJoinToRationalIndex(numOfIntegerSubtangles,initialIndex,joinedIntegerSubtangle2Index,rationalSubtangleConnectionsEM);
		
		
		//After joining, determine the NEW NE/SE/SW/NW corner information by cases:
		//RIGHT TWIST: an EVEN number of twists leaves everything fixed; we only worry about odd twists.
		//The WEST corners stay fixed; the EAST corners stay fixed if the number of twists is even but swap if its odd.
		if( (rightTwist%2) != 0 ){
			newCornerNE=cornerSE;
			newCornerSE=cornerNE;
			//For ODD twists, the parity and S2 direction of the tangle also might change after joining.
			if( rationalSubtangleParametersEM[initialIndex][5] == 0 ){
				//CASE R.1: pre-joined rational tangle had partiy 0: new parity is 1 and S2 stays the same
				rationalSubtangleParametersEM[initialIndex][5] = 1;
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 1 ){
				//CASE R.2: pre-joined rational tangle had partiy 1: new parity is 0 and S2 stays the same
				rationalSubtangleParametersEM[initialIndex][5] = 0;
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 2 ){
				//CASE R.3: pre-joined rational tangle had parity infinity: the parity stays the same and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} 
		}
		//LEFT TWIST: an EVEN number of twists leaves everything fixed; we only worry about odd twists.
		//The EAST corners stay fixed; the WEST corners stay fixed if the number of twists is even but swap if its odd.
		if( (leftTwist%2) != 0 ){
			newCornerNW=cornerSW;
			newCornerSW=cornerNW;
			//For ODD twists, the parity and S2 direction of the tangle also might change after joining.
			if( rationalSubtangleParametersEM[initialIndex][5] == 0 ){
				//CASE L.1: pre-joined rational tangle had partiy 0: new parity is 1 and S2 stays the same
				rationalSubtangleParametersEM[initialIndex][5] = 1;
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 1 ){
				//CASE L.2: pre-joined rational tangle had partiy 1: new parity is 0 and S2 stays the same
				rationalSubtangleParametersEM[initialIndex][5] = 0;
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 2 ){
				//CASE L.3: pre-joined rational tangle had parity infinity: the parity stays the same and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			}
			
			//THE NW CORNER IS SWAPPED WITH THE SW:
			//We "rotate" the entire tangle 90 degrees clockwise so force the NW corner to match corner A at the end of the check.
		}	
		
		
	} else if( (botTwist>0) && (topTwist>0) ){
		//CASE 2: two subtangles joined BOTTOM AND TOP
		localTwistCount = botTwist+topTwist;
		rationalSubtangleParametersEM[initialIndex][0] = rationalSubtangleParametersEM[initialIndex][0]+localTwistCount;
		rationalSubtangleParametersEM[initialIndex][1] = rationalSubtangleParametersEM[initialIndex][1]+2;
		integerSubtangleJoins[initialIndex] = integerSubtangleJoins[initialIndex]+2;
		rationalSubtanglesRefined = rationalSubtanglesRefined-2;
		
		//Determine the index of the two joined subtangles, and determine the corners of these tangles.
		//Update the corresponding entries in the integerSubtangleJoins array to indicate these have been joined.
		joinedIntegerSubtangle1Index = rationalSubtangleConnectionsEM[initialIndex][0][cornerSE]-1;
		joinedIntegerSubtangle2Index = rationalSubtangleConnectionsEM[initialIndex][0][cornerNE]-1;
		integerSubtangleJoins[joinedIntegerSubtangle1Index]--;
		integerSubtangleJoins[joinedIntegerSubtangle2Index]--;
		
		//Bottom twist corners:
		joinedSubtangle1Corner1 = rationalSubtangleConnectionsEM[initialIndex][1][cornerSW];
		joinedSubtangle1Corner2 = rationalSubtangleConnectionsEM[initialIndex][1][cornerSE];
		//Top twist corners:
		joinedSubtangle2Corner1 = rationalSubtangleConnectionsEM[initialIndex][1][cornerNW];
		joinedSubtangle2Corner2 = rationalSubtangleConnectionsEM[initialIndex][1][cornerNE];
		
		
		//Connected-connected subtangle information (I'm sorry):
		//BOTTOM:
		joinedSubtangle1Corner1ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][0][determineConnectedEndpoint(joinedSubtangle1Corner1)]-1;
		joinedSubtangle1Corner1ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][1][determineConnectedEndpoint(joinedSubtangle1Corner1)];
		joinedSubtangle1Corner2ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][0][determineConnectedEndpoint(joinedSubtangle1Corner2)]-1;
		joinedSubtangle1Corner2ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][1][determineConnectedEndpoint(joinedSubtangle1Corner2)];
		//TOP:
		joinedSubtangle2Corner1ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle2Index][0][determineConnectedEndpoint(joinedSubtangle2Corner1)]-1;
		joinedSubtangle2Corner1ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle2Index][1][determineConnectedEndpoint(joinedSubtangle2Corner1)];
		joinedSubtangle2Corner2ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle2Index][0][determineConnectedEndpoint(joinedSubtangle2Corner2)]-1;
		joinedSubtangle2Corner2ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle2Index][1][determineConnectedEndpoint(joinedSubtangle2Corner2)];
		
		//New connections after joining: overwrite the twisted corners with the new connections
		//BOTTOM: join to SW and SE (remember to +1 since the arrays start at 0 internally)
		rationalSubtangleConnectionsEM[initialIndex][0][cornerSW] = joinedSubtangle1Corner1ConnectedCrossIndex+1;
		rationalSubtangleConnectionsEM[initialIndex][0][cornerSE] = joinedSubtangle1Corner2ConnectedCrossIndex+1;
		//TOP: join to NW and NE (remember to +1 since the arrays start at 0 internally)
		rationalSubtangleConnectionsEM[initialIndex][0][cornerNW] = joinedSubtangle2Corner1ConnectedCrossIndex+1;
		rationalSubtangleConnectionsEM[initialIndex][0][cornerNE] = joinedSubtangle2Corner2ConnectedCrossIndex+1;
		
		//Update the connected corner information for the rational amalgam.
		//BOTTOM:
		rationalSubtangleConnectionsEM[initialIndex][1][cornerSW] = joinedSubtangle1Corner1ConnectedCrossIndexCorner;
		rationalSubtangleConnectionsEM[initialIndex][1][cornerSE] = joinedSubtangle1Corner2ConnectedCrossIndexCorner;
		//TOP:
		rationalSubtangleConnectionsEM[initialIndex][1][cornerNW] = joinedSubtangle2Corner1ConnectedCrossIndexCorner;
		rationalSubtangleConnectionsEM[initialIndex][1][cornerNE] = joinedSubtangle2Corner2ConnectedCrossIndexCorner;
		
		//Set the corner connection information of the connected-connected subtangle to match that of the rational amalgam.
		//BOTTOM:
		rationalSubtangleConnectionsEM[joinedSubtangle1Corner1ConnectedCrossIndex][1][joinedSubtangle1Corner1ConnectedCrossIndexCorner] = cornerSW;
		rationalSubtangleConnectionsEM[joinedSubtangle1Corner2ConnectedCrossIndex][1][joinedSubtangle1Corner2ConnectedCrossIndexCorner] = cornerSE;
		//TOP:
		rationalSubtangleConnectionsEM[joinedSubtangle2Corner1ConnectedCrossIndex][1][joinedSubtangle2Corner1ConnectedCrossIndexCorner] = cornerNW;
		rationalSubtangleConnectionsEM[joinedSubtangle2Corner2ConnectedCrossIndex][1][joinedSubtangle2Corner2ConnectedCrossIndexCorner] = cornerNE;
		
		
		//Update references to joined subtangle index:
		//BOTTOM:
		updateIntegerJoinToRationalIndex(numOfIntegerSubtangles,initialIndex,joinedIntegerSubtangle1Index,rationalSubtangleConnectionsEM);
		//TOP:
		updateIntegerJoinToRationalIndex(numOfIntegerSubtangles,initialIndex,joinedIntegerSubtangle2Index,rationalSubtangleConnectionsEM);
		
		
		//After joining, determine the NEW NE/SE/SW/NW corner information by cases:
		//BOTTOM TWIST: an EVEN number of twists leaves everything fixed; we only worry about odd twists.
		//The NORTH corners stay fixed; the SOUTH corners stay fixed if the number of twists is even but swap if its odd.
		if( (botTwist%2) != 0 ){
			newCornerSW=cornerSE;
			newCornerSE=cornerSW;
			//For ODD twists, the parity and S2 direction of the tangle also might change after joining.
			if( rationalSubtangleParametersEM[initialIndex][5] == 0 ){
				//CASE B.1: pre-joined rational tangle had partiy 0: the parity stays the same and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 1 ){
				//CASE B.2: pre-joined rational tangle had partiy 1: new parity is infinity and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][5] = 2;
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 2 ){
				//CASE B.3: pre-joined rational tangle had parity infinity: new parity is 1 and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][5] = 1;
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} 
		}
		//TOP TWIST: an EVEN number of twists leaves everything fixed; we only worry about odd twists.
		//The SOUTH corners stay fixed; the NORTH corners stay fixed if the number of twists is even but swap if its odd.
		if( (topTwist%2) != 0 ){
			newCornerNW=cornerNE;
			newCornerNE=cornerNW;
			//For ODD twists, the parity and S2 direction of the tangle also might change after joining.
			if( rationalSubtangleParametersEM[initialIndex][5] == 0 ){
				//CASE T.1: pre-joined rational tangle had partiy 0: the parity stays the same and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 1 ){
				//CASE T.2: pre-joined rational tangle had partiy 1: new parity is infinity and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][5] = 2;
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 2 ){
				//CASE T.3: pre-joined rational tangle had parity infinity: new parity is 1 and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][5] = 1;
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} 
			
			//THE NW CORNER IS SWAPPED WITH THE NE:
			//We "rotate" the entire tangle 90 degrees counterclockwise so force the NW corner to match corner A at the end of the check.
		}	
		
		
	} else if( rightTwist>0 ){
		//CASE 3: one subtangle joined RIGHT
		localTwistCount = rightTwist;
		rationalSubtangleParametersEM[initialIndex][0] = rationalSubtangleParametersEM[initialIndex][0]+localTwistCount;
		rationalSubtangleParametersEM[initialIndex][1]++;
		integerSubtangleJoins[initialIndex]++;
		rationalSubtanglesRefined--;
		
		//Determine the index of the joined tangle and update the corresponding integerSubtangleJoins array entry.
		joinedIntegerSubtangle1Index = rationalSubtangleConnectionsEM[initialIndex][0][cornerNE]-1;
		integerSubtangleJoins[joinedIntegerSubtangle1Index]--;
		//Right twist corners:
		joinedSubtangle1Corner1 = rationalSubtangleConnectionsEM[initialIndex][1][cornerNE];
		joinedSubtangle1Corner2 = rationalSubtangleConnectionsEM[initialIndex][1][cornerSE];
		
		//Connected-connected subtangle information (I'm sorry):
		joinedSubtangle1Corner1ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][0][determineConnectedEndpoint(joinedSubtangle1Corner1)]-1;
		joinedSubtangle1Corner1ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][1][determineConnectedEndpoint(joinedSubtangle1Corner1)];
		joinedSubtangle1Corner2ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][0][determineConnectedEndpoint(joinedSubtangle1Corner2)]-1;
		joinedSubtangle1Corner2ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][1][determineConnectedEndpoint(joinedSubtangle1Corner2)];
		
		//DEBUG:
		//printf("\n joinedSubtangle1: %d+1 \n Corner1: %d , Corner2: %d \n", joinedIntegerSubtangle1Index,joinedSubtangle1Corner1,joinedSubtangle1Corner2);
		//printf("\n Connected-connected 1-1: %d+1 , corner: %d \n Connected-connected 1-2: %d+1 , corner: %d \n", joinedSubtangle1Corner1ConnectedCrossIndex,joinedSubtangle1Corner1ConnectedCrossIndexCorner,joinedSubtangle1Corner2ConnectedCrossIndex,joinedSubtangle1Corner2ConnectedCrossIndexCorner);
		
		//New connections after joining: overwrite the twisted corners with the new connections
		//RIGHT: join to NE and SE (remember to +1 since the arrays start at 0 internally)
		rationalSubtangleConnectionsEM[initialIndex][0][cornerNE] = joinedSubtangle1Corner1ConnectedCrossIndex+1;
		rationalSubtangleConnectionsEM[initialIndex][0][cornerSE] = joinedSubtangle1Corner2ConnectedCrossIndex+1;
		//Update the connected corner information for the rational amalgam.
		rationalSubtangleConnectionsEM[initialIndex][1][cornerNE] = joinedSubtangle1Corner1ConnectedCrossIndexCorner;
		rationalSubtangleConnectionsEM[initialIndex][1][cornerSE] = joinedSubtangle1Corner2ConnectedCrossIndexCorner;
		
		//Set the corner connection information of the connected-connected subtangle to match that of the rational amalgam.
		rationalSubtangleConnectionsEM[joinedSubtangle1Corner1ConnectedCrossIndex][1][joinedSubtangle1Corner1ConnectedCrossIndexCorner] = cornerNE;
		rationalSubtangleConnectionsEM[joinedSubtangle1Corner2ConnectedCrossIndex][1][joinedSubtangle1Corner2ConnectedCrossIndexCorner] = cornerSE;
		
		//Update references to joined subtangle index:
		updateIntegerJoinToRationalIndex(numOfIntegerSubtangles,initialIndex,joinedIntegerSubtangle1Index,rationalSubtangleConnectionsEM);
		
		//After joining, determine the NEW NE/SE/SW/NW corner information by cases:
		//RIGHT TWIST: an EVEN number of twists leaves everything fixed; we only worry about odd twists.
		//The WEST corners stay fixed; the EAST corners stay fixed if the number of twists is even but swap if its odd.
		if( (rightTwist%2) != 0 ){
			newCornerNE=cornerSE;
			newCornerSE=cornerNE;
			//For ODD twists, the parity and S2 direction of the tangle also might change after joining.
			if( rationalSubtangleParametersEM[initialIndex][5] == 0 ){
				//CASE R.1: pre-joined rational tangle had partiy 0: new parity is 1 and S2 stays the same
				rationalSubtangleParametersEM[initialIndex][5] = 1;
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 1 ){
				//CASE R.2: pre-joined rational tangle had partiy 1: new parity is 0 and S2 stays the same
				rationalSubtangleParametersEM[initialIndex][5] = 0;
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 2 ){
				//CASE R.3: pre-joined rational tangle had parity infinity: the parity stays the same and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} 
		}
			
		
	} else if( botTwist>0 ){
		//CASE 4: one subtangle joined BOTTOM
		localTwistCount = botTwist;
		rationalSubtangleParametersEM[initialIndex][0] = rationalSubtangleParametersEM[initialIndex][0]+localTwistCount;
		rationalSubtangleParametersEM[initialIndex][1]++;
		integerSubtangleJoins[initialIndex]++;
		rationalSubtanglesRefined--;
		
		//Determine the index of the joined tangle and update the corresponding integerSubtangleJoins array entry.
		joinedIntegerSubtangle1Index = rationalSubtangleConnectionsEM[initialIndex][0][cornerSE]-1;
		integerSubtangleJoins[joinedIntegerSubtangle1Index]--;
		//Bottom twist corners:
		joinedSubtangle1Corner1 = rationalSubtangleConnectionsEM[initialIndex][1][cornerSW];
		joinedSubtangle1Corner2 = rationalSubtangleConnectionsEM[initialIndex][1][cornerSE];
		
		//Connected-connected subtangle information (I'm sorry):
		joinedSubtangle1Corner1ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][0][determineConnectedEndpoint(joinedSubtangle1Corner1)]-1;
		joinedSubtangle1Corner1ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][1][determineConnectedEndpoint(joinedSubtangle1Corner1)];
		joinedSubtangle1Corner2ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][0][determineConnectedEndpoint(joinedSubtangle1Corner2)]-1;
		joinedSubtangle1Corner2ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][1][determineConnectedEndpoint(joinedSubtangle1Corner2)];
		
		//New connections after joining: overwrite the twisted corners with the new connections
		//BOTTOM: join to SW and SE (remember to +1 since the arrays start at 0 internally)
		rationalSubtangleConnectionsEM[initialIndex][0][cornerSW] = joinedSubtangle1Corner1ConnectedCrossIndex+1;
		rationalSubtangleConnectionsEM[initialIndex][0][cornerSE] = joinedSubtangle1Corner2ConnectedCrossIndex+1;
		//Update the connected corner information for the rational amalgam.
		rationalSubtangleConnectionsEM[initialIndex][1][cornerSW] = joinedSubtangle1Corner1ConnectedCrossIndexCorner;
		rationalSubtangleConnectionsEM[initialIndex][1][cornerSE] = joinedSubtangle1Corner2ConnectedCrossIndexCorner;
		
		//Set the corner connection information of the connected-connected subtangle to match that of the rational amalgam.
		rationalSubtangleConnectionsEM[joinedSubtangle1Corner1ConnectedCrossIndex][1][joinedSubtangle1Corner1ConnectedCrossIndexCorner] = cornerSW;
		rationalSubtangleConnectionsEM[joinedSubtangle1Corner2ConnectedCrossIndex][1][joinedSubtangle1Corner2ConnectedCrossIndexCorner] = cornerSE;
		
		//Update references to joined subtangle index:
		updateIntegerJoinToRationalIndex(numOfIntegerSubtangles,initialIndex,joinedIntegerSubtangle1Index,rationalSubtangleConnectionsEM);
		
		//After joining, determine the NEW NE/SE/SW/NW corner information by cases:
		//BOTTOM TWIST: an EVEN number of twists leaves everything fixed; we only worry about odd twists.
		//The NORTH corners stay fixed; the SOUTH corners stay fixed if the number of twists is even but swap if its odd.
		if( (botTwist%2) != 0 ){
			newCornerSW=cornerSE;
			newCornerSE=cornerSW;
			//For ODD twists, the parity and S2 direction of the tangle also might change after joining.
			if( rationalSubtangleParametersEM[initialIndex][5] == 0 ){
				//CASE B.1: pre-joined rational tangle had partiy 0: the parity stays the same and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 1 ){
				//CASE B.2: pre-joined rational tangle had partiy 1: new parity is infinity and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][5] = 2;
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 2 ){
				//CASE B.3: pre-joined rational tangle had parity infinity: new parity is 1 and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][5] = 1;
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} 
		}
		
		
	} else if( leftTwist>0 ){
		//CASE 5: one subtangle joined LEFT
		localTwistCount = leftTwist;
		rationalSubtangleParametersEM[initialIndex][0] = rationalSubtangleParametersEM[initialIndex][0]+localTwistCount;
		rationalSubtangleParametersEM[initialIndex][1]++;
		integerSubtangleJoins[initialIndex]++;
		rationalSubtanglesRefined--;
		
		//Determine the index of the joined tangle and update the corresponding integerSubtangleJoins array entry.
		joinedIntegerSubtangle1Index = rationalSubtangleConnectionsEM[initialIndex][0][cornerNW]-1;
		integerSubtangleJoins[joinedIntegerSubtangle1Index]--;
		//Left twist corners:
		joinedSubtangle1Corner1 = rationalSubtangleConnectionsEM[initialIndex][1][cornerNW];
		joinedSubtangle1Corner2 = rationalSubtangleConnectionsEM[initialIndex][1][cornerSW];
		
		//Connected-connected subtangle information (I'm sorry):
		joinedSubtangle1Corner1ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][0][determineConnectedEndpoint(joinedSubtangle1Corner1)]-1;
		joinedSubtangle1Corner1ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][1][determineConnectedEndpoint(joinedSubtangle1Corner1)];
		joinedSubtangle1Corner2ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][0][determineConnectedEndpoint(joinedSubtangle1Corner2)]-1;
		joinedSubtangle1Corner2ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][1][determineConnectedEndpoint(joinedSubtangle1Corner2)];
		
		//New connections after joining: overwrite the twisted corners with the new connections
		//LEFT: join to NW and SW (remember to +1 since the arrays start at 0 internally)
		rationalSubtangleConnectionsEM[initialIndex][0][cornerNW] = joinedSubtangle1Corner1ConnectedCrossIndex+1;
		rationalSubtangleConnectionsEM[initialIndex][0][cornerSW] = joinedSubtangle1Corner2ConnectedCrossIndex+1;
		//Update the connected corner information for the rational amalgam.
		rationalSubtangleConnectionsEM[initialIndex][1][cornerNW] = joinedSubtangle1Corner1ConnectedCrossIndexCorner;
		rationalSubtangleConnectionsEM[initialIndex][1][cornerSW] = joinedSubtangle1Corner2ConnectedCrossIndexCorner;
		
		//Set the corner connection information of the connected-connected subtangle to match that of the rational amalgam.
		rationalSubtangleConnectionsEM[joinedSubtangle1Corner1ConnectedCrossIndex][1][joinedSubtangle1Corner1ConnectedCrossIndexCorner] = cornerNW;
		rationalSubtangleConnectionsEM[joinedSubtangle1Corner2ConnectedCrossIndex][1][joinedSubtangle1Corner2ConnectedCrossIndexCorner] = cornerSW;
		
		//Update references to joined subtangle index:
		updateIntegerJoinToRationalIndex(numOfIntegerSubtangles,initialIndex,joinedIntegerSubtangle1Index,rationalSubtangleConnectionsEM);
		
		//After joining, determine the NEW NE/SE/SW/NW corner information by cases:
		//LEFT TWIST: an EVEN number of twists leaves everything fixed; we only worry about odd twists.
		//The EAST corners stay fixed; the WEST corners stay fixed if the number of twists is even but swap if its odd.
		if( (leftTwist%2) != 0 ){
			newCornerNW=cornerSW;
			newCornerSW=cornerNW;
			//For ODD twists, the parity and S2 direction of the tangle also might change after joining.
			if( rationalSubtangleParametersEM[initialIndex][5] == 0 ){
				//CASE L.1: pre-joined rational tangle had partiy 0: new parity is 1 and S2 stays the same
				rationalSubtangleParametersEM[initialIndex][5] = 1;
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 1 ){
				//CASE L.2: pre-joined rational tangle had partiy 1: new parity is 0 and S2 stays the same
				rationalSubtangleParametersEM[initialIndex][5] = 0;
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 2 ){
				//CASE L.3: pre-joined rational tangle had parity infinity: the parity stays the same and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			}
			
			//THE NW CORNER IS SWAPPED WITH THE SW:
			//We "rotate" the entire tangle 90 degrees clockwise so force the NW corner to match corner A at the end of the check.
		}	
		
		
	} else if( topTwist>0 ){
		//CASE 6: one subtangle joined TOP
		localTwistCount = topTwist;
		rationalSubtangleParametersEM[initialIndex][0] = rationalSubtangleParametersEM[initialIndex][0]+localTwistCount;
		rationalSubtangleParametersEM[initialIndex][1]++;
		integerSubtangleJoins[initialIndex]++;
		rationalSubtanglesRefined--;
		
		//Determine the index of the joined tangle and update the corresponding integerSubtangleJoins array entry.
		joinedIntegerSubtangle1Index = rationalSubtangleConnectionsEM[initialIndex][0][cornerNE]-1;
		integerSubtangleJoins[joinedIntegerSubtangle1Index]--;
		//Top twist corners:
		joinedSubtangle1Corner1 = rationalSubtangleConnectionsEM[initialIndex][1][cornerNW];
		joinedSubtangle1Corner2 = rationalSubtangleConnectionsEM[initialIndex][1][cornerNE];
		
		//Connected-connected subtangle information (I'm sorry):
		joinedSubtangle1Corner1ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][0][determineConnectedEndpoint(joinedSubtangle1Corner1)]-1;
		joinedSubtangle1Corner1ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][1][determineConnectedEndpoint(joinedSubtangle1Corner1)];
		joinedSubtangle1Corner2ConnectedCrossIndex = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][0][determineConnectedEndpoint(joinedSubtangle1Corner2)]-1;
		joinedSubtangle1Corner2ConnectedCrossIndexCorner = rationalSubtangleConnectionsEM[joinedIntegerSubtangle1Index][1][determineConnectedEndpoint(joinedSubtangle1Corner2)];
		
		//New connections after joining: overwrite the twisted corners with the new connections
		//TOP: join to NW and NE (remember to +1 since the arrays start at 0 internally)
		rationalSubtangleConnectionsEM[initialIndex][0][cornerNW] = joinedSubtangle1Corner1ConnectedCrossIndex+1;
		rationalSubtangleConnectionsEM[initialIndex][0][cornerNE] = joinedSubtangle1Corner2ConnectedCrossIndex+1;
		//Update the connected corner information for the rational amalgam.
		rationalSubtangleConnectionsEM[initialIndex][1][cornerNW] = joinedSubtangle1Corner1ConnectedCrossIndexCorner;
		rationalSubtangleConnectionsEM[initialIndex][1][cornerNE] = joinedSubtangle1Corner2ConnectedCrossIndexCorner;
		
		//Set the corner connection information of the connected-connected subtangle to match that of the rational amalgam.
		rationalSubtangleConnectionsEM[joinedSubtangle1Corner1ConnectedCrossIndex][1][joinedSubtangle1Corner1ConnectedCrossIndexCorner] = cornerNW;
		rationalSubtangleConnectionsEM[joinedSubtangle1Corner2ConnectedCrossIndex][1][joinedSubtangle1Corner2ConnectedCrossIndexCorner] = cornerNE;
		
		//Update references to joined subtangle index:
		updateIntegerJoinToRationalIndex(numOfIntegerSubtangles,initialIndex,joinedIntegerSubtangle1Index,rationalSubtangleConnectionsEM);
		
		//After joining, determine the NEW NE/SE/SW/NW corner information by cases:
		//TOP TWIST: an EVEN number of twists leaves everything fixed; we only worry about odd twists.
		//The SOUTH corners stay fixed; the NORTH corners stay fixed if the number of twists is even but swap if its odd.
		if( (topTwist%2) != 0 ){
			newCornerNW=cornerNE;
			newCornerNE=cornerNW;
			//For ODD twists, the parity and S2 direction of the tangle also might change after joining.
			if( rationalSubtangleParametersEM[initialIndex][5] == 0 ){
				//CASE T.1: pre-joined rational tangle had partiy 0: the parity stays the same and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 1 ){
				//CASE T.2: pre-joined rational tangle had partiy 1: new parity is infinity and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][5] = 2;
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} else if( rationalSubtangleParametersEM[initialIndex][5] == 2 ){
				//CASE T.3: pre-joined rational tangle had parity infinity: new parity is 1 and S2 is reversed
				rationalSubtangleParametersEM[initialIndex][5] = 1;
				rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
			} 
			
			//THE NW CORNER IS SWAPPED WITH THE NE:
			//We "rotate" the entire tangle 90 degrees counterclockwise so force the NW corner to match corner A at the end of the check.
		}	
		
		
	} else {
		//If none of these happen, this tangle has no twisting.
		localTwistCount=0;
	}
	
	
	if( localTwistCount > 0 ){
		//Increment the length of the current rational twist vector, and add a new entry to the vector. This entry is localTwistCount.
		rationalTwistVectorArray[initialIndex][0]++;
		rationalTwistVectorArray[initialIndex][rationalTwistVectorArray[initialIndex][0]] = localTwistCount;
		
		//Also update the "twist vector length" entry in the rationalSubtangleParemetersEM array.
		rationalSubtangleParametersEM[initialIndex][2]++;
		
	}
	
	
	//It's possible that an endpoint subtangle got joined to the rational amalgam.
	//If so, we need to update the corresponding entries in the endpointConnectionsCorners array.
	//The NW and SW endpoints will always match an A or D corner, respectively.
	//The NE and SE could be B/D and A/C, and this CAN CHANGE after joining an integer endpoint tangle to a rational amalgam.
	//Unfortunately, we need to check each of the six cases above separately for how to update the endpoint corners.
	if( localTwistCount > 0 ){
		//If any twisting happened:
		if( (rightTwist>0) && (leftTwist>0) ){
			//CASE 1: right and left twisting
			//RIGHT: check joinedIntegerSubtangle1
			for(int i=0; i<4; i++){
				//Check if the joined integer subtangle1 was any endpoint tangle.
				if( (endpointConnectionsCorners[0][i]-1) == joinedIntegerSubtangle1Index ){
					//If yes, we replace the reference to this joined subtangle with the initial index (+1, to match the indexing)
					endpointConnectionsCorners[0][i] = (initialIndex+1);
					//Next, check to see if the endpoint corner was opposite one of the joined corners (it must be for one of them).
					if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle1Corner1) ){
						//If the other side of the FIRST joined subtangle1 corner is an endpoint, the new endpoint corner is cornerNE.
						endpointConnectionsCorners[1][i] = cornerNE;
					} else if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle1Corner2) ){
						//If the other side of the SECOND joined subtangle1 corner is an endpoint, the new endpoint corner is cornerSE.
						endpointConnectionsCorners[1][i] = cornerSE;
					}
				}
			}
			//LEFT: check joinedIntegerSubtangle2
			for(int i=0; i<4; i++){
				//Check if the joined integer subtangle2 was any endpoint tangle.
				if( (endpointConnectionsCorners[0][i]-1) == joinedIntegerSubtangle2Index ){
					//If yes, we replace the reference to this joined subtangle with the initial index (+1, to match the indexing)
					endpointConnectionsCorners[0][i] = (initialIndex+1);
					//Next, check to see if the endpoint corner was opposite one of the joined corners (it must be for one of them).
					if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle2Corner1) ){
						//If the other side of the FIRST joined subtangle1 corner is an endpoint, the new endpoint corner is cornerNW.
						endpointConnectionsCorners[1][i] = cornerNW;
					} else if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle2Corner2) ){
						//If the other side of the SECOND joined subtangle1 corner is an endpoint, the new endpoint corner is cornerSW.
						endpointConnectionsCorners[1][i] = cornerSW;
					}
				}
			}
			
		} else if( (botTwist>0) && (topTwist>0) ){
			//CASE 2: bottom and top twisting
			//BOTTOM: check joinedIntegerSubtangle1
			for(int i=0; i<4; i++){
				//Check if the joined integer subtangle1 was any endpoint tangle.
				if( (endpointConnectionsCorners[0][i]-1) == joinedIntegerSubtangle1Index ){
					//If yes, we replace the reference to this joined subtangle with the initial index (+1, to match the indexing)
					endpointConnectionsCorners[0][i] = (initialIndex+1);
					//Next, check to see if the endpoint corner was opposite one of the joined corners (it must be for one of them).
					if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle1Corner1) ){
						//If the other side of the FIRST joined subtangle1 corner is an endpoint, the new endpoint corner is cornerSW.
						endpointConnectionsCorners[1][i] = cornerSW;
					} else if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle1Corner2) ){
						//If the other side of the SECOND joined subtangle1 corner is an endpoint, the new endpoint corner is cornerSE.
						endpointConnectionsCorners[1][i] = cornerSE;
					}
				}
			}
			//TOP: check joinedIntegerSubtangle2
			for(int i=0; i<4; i++){
				//Check if the joined integer subtangle2 was any endpoint tangle.
				if( (endpointConnectionsCorners[0][i]-1) == joinedIntegerSubtangle2Index ){
					//If yes, we replace the reference to this joined subtangle with the initial index (+1, to match the indexing)
					endpointConnectionsCorners[0][i] = (initialIndex+1);
					//Next, check to see if the endpoint corner was opposite one of the joined corners (it must be for one of them).
					if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle2Corner1) ){
						//If the other side of the FIRST joined subtangle1 corner is an endpoint, the new endpoint corner is cornerNW.
						endpointConnectionsCorners[1][i] = cornerNW;
					} else if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle2Corner2) ){
						//If the other side of the SECOND joined subtangle1 corner is an endpoint, the new endpoint corner is cornerNE.
						endpointConnectionsCorners[1][i] = cornerNE;
					}
				}
			}
			
		} else if( rightTwist>0 ){
			//CASE 3: only right twisting
			for(int i=0; i<4; i++){
				//Check if the joined integer subtangle1 was any endpoint tangle.
				if( (endpointConnectionsCorners[0][i]-1) == joinedIntegerSubtangle1Index ){
					//If yes, we replace the reference to this joined subtangle with the initial index (+1, to match the indexing)
					endpointConnectionsCorners[0][i] = (initialIndex+1);
					//Next, check to see if the endpoint corner was opposite one of the joined corners (it must be for one of them).
					if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle1Corner1) ){
						//If the other side of the FIRST joined subtangle1 corner is an endpoint, the new endpoint corner is cornerNE.
						endpointConnectionsCorners[1][i] = cornerNE;
					} else if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle1Corner2) ){
						//If the other side of the SECOND joined subtangle1 corner is an endpoint, the new endpoint corner is cornerSE.
						endpointConnectionsCorners[1][i] = cornerSE;
					}
				}
			}
			
		} else if( botTwist>0 ){
			//CASE 4: only bottom twisting
			for(int i=0; i<4; i++){
				//Check if the joined integer subtangle1 was any endpoint tangle.
				if( (endpointConnectionsCorners[0][i]-1) == joinedIntegerSubtangle1Index ){
					//If yes, we replace the reference to this joined subtangle with the initial index (+1, to match the indexing)
					endpointConnectionsCorners[0][i] = (initialIndex+1);
					//Next, check to see if the endpoint corner was opposite one of the joined corners (it must be for one of them).
					if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle1Corner1) ){
						//If the other side of the FIRST joined subtangle1 corner is an endpoint, the new endpoint corner is cornerSW.
						endpointConnectionsCorners[1][i] = cornerSW;
					} else if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle1Corner2) ){
						//If the other side of the SECOND joined subtangle1 corner is an endpoint, the new endpoint corner is cornerSE.
						endpointConnectionsCorners[1][i] = cornerSE;
					}
				}
			}
			
		} else if( leftTwist>0 ){
			//CASE 5: only left twisting
			for(int i=0; i<4; i++){
				//Check if the joined integer subtangle1 was any endpoint tangle.
				if( (endpointConnectionsCorners[0][i]-1) == joinedIntegerSubtangle1Index ){
					//If yes, we replace the reference to this joined subtangle with the initial index (+1, to match the indexing)
					endpointConnectionsCorners[0][i] = (initialIndex+1);
					//Next, check to see if the endpoint corner was opposite one of the joined corners (it must be for one of them).
					if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle1Corner1) ){
						//If the other side of the FIRST joined subtangle1 corner is an endpoint, the new endpoint corner is cornerNW.
						endpointConnectionsCorners[1][i] = cornerNW;
					} else if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle1Corner2) ){
						//If the other side of the SECOND joined subtangle1 corner is an endpoint, the new endpoint corner is cornerSW.
						endpointConnectionsCorners[1][i] = cornerSW;
					}
				}
			}
			
		} else if( topTwist>0 ){
			//CASE 6: only top twisting
			for(int i=0; i<4; i++){
				//Check if the joined integer subtangle1 was any endpoint tangle.
				if( (endpointConnectionsCorners[0][i]-1) == joinedIntegerSubtangle1Index ){
					//If yes, we replace the reference to this joined subtangle with the initial index (+1, to match the indexing)
					endpointConnectionsCorners[0][i] = (initialIndex+1);
					//Next, check to see if the endpoint corner was opposite one of the joined corners (it must be for one of them).
					if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle1Corner1) ){
						//If the other side of the FIRST joined subtangle1 corner is an endpoint, the new endpoint corner is cornerNW.
						endpointConnectionsCorners[1][i] = cornerNW;
					} else if( endpointConnectionsCorners[1][i] == determineConnectedEndpoint(joinedSubtangle1Corner2) ){
						//If the other side of the SECOND joined subtangle1 corner is an endpoint, the new endpoint corner is cornerNE.
						endpointConnectionsCorners[1][i] = cornerNE;
					}
				}
			}
		}
	}
	
	
	//Update the rational Gauss code information to delete references to the joined tangles.
	//If there is any subtangle joining, the rational Gauss EM array must be updated.
	int removedIndexCount=0;
	if( localTwistCount > 0 ){
		//If multiple twisting happens at the same time, then two Gauss array entries must be removed; otherwise, only one needs to be.
		if( ((rightTwist>0) && (leftTwist>0)) || ((botTwist>0) && (topTwist>0)) ){
			for(int i=0; i<(2*oldRationalSubtanglesRefined)-removedIndexCount; i++){
				if( gaussRationalSubtangleEM[i] == (joinedIntegerSubtangle1Index+1) ){
					//If the joined subtangle1 index matches a rational Gauss code entry, remove it and shift all remaining entries down one position.
					removedIndexCount++;
					//printf("\n DEBUG FLAG: i = %d , removedIndexCount = %d \n", i,removedIndexCount);
					for(int j=i; j<(2*oldRationalSubtanglesRefined)-removedIndexCount; j++){
						gaussRationalSubtangleEM[j]=gaussRationalSubtangleEM[j+1];
					}
					//Also decrease the values of the bars correspondingly. The last bar is always decreased, but the first is only decreased if the removed subtangle precedes it.
					if( i < barsGaussRationalSubtangleEM[0] ){
						barsGaussRationalSubtangleEM[0]--;
					}
					barsGaussRationalSubtangleEM[1]--;
					
					//Decrement i to account for the shifting.
					i--;
				}
			}
			//Repeat the same check for the joined integer subtangle2.
			for(int i=0; i<(2*oldRationalSubtanglesRefined)-removedIndexCount; i++){
				if( gaussRationalSubtangleEM[i] == (joinedIntegerSubtangle2Index+1) ){
					//If the joined subtangle1 index matches a rational Gauss code entry, remove it and shift all remaining entries down one position.
					removedIndexCount++;
					for(int j=i; j<(2*oldRationalSubtanglesRefined)-removedIndexCount; j++){
						gaussRationalSubtangleEM[j]=gaussRationalSubtangleEM[j+1];
					}
					//Also decrease the values of the bars correspondingly. The last bar is always decreased, but the first is only decreased if the removed subtangle precedes it.
					if( i < barsGaussRationalSubtangleEM[0] ){
						barsGaussRationalSubtangleEM[0]--;
					}
					barsGaussRationalSubtangleEM[1]--;
					
					//Decrement i to account for the shifting.
					i--;
				}
			}
		} else {
			for(int i=0; i<(2*oldRationalSubtanglesRefined)-removedIndexCount; i++){
				if( gaussRationalSubtangleEM[i] == (joinedIntegerSubtangle1Index+1) ){
					//If the joined subtangle1 index matches a rational Gauss code entry, remove it and shift all remaining entries down one position.
					removedIndexCount++;
					for(int j=i; j<(2*oldRationalSubtanglesRefined)-removedIndexCount; j++){
						gaussRationalSubtangleEM[j]=gaussRationalSubtangleEM[j+1];
					}
					//Also decrease the values of the bars correspondingly. The last bar is always decreased, but the first is only decreased if the removed subtangle precedes it.
					if( i < barsGaussRationalSubtangleEM[0] ){
						barsGaussRationalSubtangleEM[0]--;
					}
					barsGaussRationalSubtangleEM[1]--;
					
					//printf("\n DEBUG FLAG: i = %d , removedIndexCount = %d \n",i,removedIndexCount);
					//Decrement i to account for the shifting.
					i--;	
				}
			}
		}
	}
	
	//DEBUG:
	/*
	printf("\n DEBUG FLAG: joinedIntegerSubtangle1Index: %d ",joinedIntegerSubtangle1Index);
	printf("\n gaussRationalSubtangleEM: [");
	for(int i=0; i<2*numOfIntegerSubtangles; i++){
		printf(" %d ", gaussRationalSubtangleEM[i]);
	}
	printf("]\n");
	*/
	
	
	//Twisting shuffles which of the endpoints match NW/NE/SE/SW, but the indetifications are determined by case considerations on the number and type of twists.
	//Overwite the original corner information with the new. Note that this does not change rationSubtangleConnectionsEM, it just tracks how we identify the corners.
	//If the NW corner was changed, we ALSO rotate the tangle to force corner A to be in the NW position.
	//This can happen if there are an ODD number of left twists OR top twists; note that these two cannot happen simultaneously.
	if( (leftTwist%2) != 0 ){
		//Odd number of left twists: rotate tangle 90 degrees clockwise.
		cornerNW=newCornerSW;
		cornerNE=newCornerNW;
		cornerSE=newCornerNE;
		cornerSW=newCornerSE;
		
		//One consequence of a 90 degree rotation is that the twist sign of the rational amalgam changes.
		//A second consequence of rotating is that we switch the type of initial integer tangle (horizontal or vertical).
		//These are tracked in rationalSubtangleParametersEM[initialIndex][] array entries 7 and 9, respectively.
		rationalSubtangleParametersEM[initialIndex][7] = -1*rationalSubtangleParametersEM[initialIndex][7];
		rationalSubtangleParametersEM[initialIndex][9] = -1*rationalSubtangleParametersEM[initialIndex][9];
		
		//Rotating the tangle can also change the parity and the direction of S2.
		if( rationalSubtangleParametersEM[initialIndex][5] == 0 ){
			//Case 1: parity 0 changes to infinity; S2 stays the same
			rationalSubtangleParametersEM[initialIndex][5] = 2;
		} else if( rationalSubtangleParametersEM[initialIndex][5] == 1 ){
			//Case 2: parity 1 stays the same; S2 is reversed
			rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
		} else if( rationalSubtangleParametersEM[initialIndex][5] == 2 ){
			//Case 3: parity infinity chages to 0; S2 stays the same
			rationalSubtangleParametersEM[initialIndex][5] = 0;
		} else {
			//Error Case: this shouldn't happen, so print large numbers to make this error identifiable at a glance if it does.
			rationalSubtangleParametersEM[initialIndex][5]=100;
			rationalSubtangleParametersEM[initialIndex][6]=100;
		}
		
	} else if( (topTwist%2) != 0 ){
		//Odd number of top twists: rotate the tangle 90 degrees counterclockwise.
		cornerNW=newCornerNE;
		cornerSW=newCornerNW;
		cornerSE=newCornerSW;
		cornerNE=newCornerSE;
		
		//One consequence of a 90 degree rotation is that the twist sign of the rational amalgam changes.
		//A second consequence of rotating is that we switch the type of initial integer tangle (horizontal or vertical).
		//These are tracked in rationalSubtangleParametersEM[initialIndex][] array entries 7 and 9, respectively.
		rationalSubtangleParametersEM[initialIndex][7] = -1*rationalSubtangleParametersEM[initialIndex][7];
		rationalSubtangleParametersEM[initialIndex][9] = -1*rationalSubtangleParametersEM[initialIndex][9];
		
		//Rotating the tangle can also change the parity and the direction of S2.
		if( rationalSubtangleParametersEM[initialIndex][5] == 0 ){
			//Case 1: parity 0 changes to infinity; S2 stays the same
			rationalSubtangleParametersEM[initialIndex][5] = 2;
		} else if( rationalSubtangleParametersEM[initialIndex][5] == 1 ){
			//Case 2: parity 1 stays the same; S2 is reversed
			rationalSubtangleParametersEM[initialIndex][6] = -1*rationalSubtangleParametersEM[initialIndex][6];
		} else if( rationalSubtangleParametersEM[initialIndex][5] == 2 ){
			//Case 3: parity infinity chages to 0; S2 stays the same
			rationalSubtangleParametersEM[initialIndex][5] = 0;
		} else {
			//Error Case: this shouldn't happen, so print large numbers to make this error identifiable at a glance if it does.
			rationalSubtangleParametersEM[initialIndex][5]=100;
			rationalSubtangleParametersEM[initialIndex][6]=100;
		}
		
	} else {
		//If corner A was not changed, we leave the tangle as is and set these to be the corners.
		cornerNW=newCornerNW;
		cornerNE=newCornerNE;
		cornerSE=newCornerSE;
		cornerSW=newCornerSW;
	}
	
	
	//Update the number of crossings in the current rational amalgam? Is this something worth tracking?
	twistCount = (twistCount+localTwistCount);
	
			
}





/*
(NC, 9/23/18)
Prototype of a function to build up the generalized EM code for rational subtangles.
This is again a work in progress, and I will rewrite the description after I figure out how to make it work.

Inputs:
numOfSubtangles
integerSubtangleConnectionsEM[][2][5]
integerSubtangleParametersEM[][7]
gaussIntegerSubtangleEM[]
barsGaussIntegerSubtangleEM[2]
rationalTwistVectorArray[][31]

Outputs (in the form of modified global variables):
numOfSubtangles
rationalSubtangleConnectionsEM[][2][5]
rationalSubtangleParametersEM[][7]
rationalTwistVectorsEM[][]
gaussRationalSubtangleEM[]
barsGaussRationalSubtangleEM[2]


This function calls:
	detectCornerTwisting()
	joinIntegerToRational()
This function is called by:
	buildGeneralizedEMCode()

*/

void buildRationalSubtangleEMCode(int &numOfSubtangles, int integerSubtangleConnectionsEM[][2][5], int integerSubtangleParametersEM[][7], int *gaussIntegerSubtangleEM, int *barsGaussIntegerSubtangleEM, int rationalSubtangleConnectionsEM[][2][5], int rationalSubtangleParametersEM[][10], int rationalTwistVectorsEM[][31], int *gaussRationalSubtangleEM, int *barsGaussRationalSubtangleEM, int endpointConnectionsCornersRational[2][4], int rationalComponentsCanonicalConfiguration[][4]){
	
	//Initialize the gaussRationalSubtangleEM array and bars to be the same as the gaussIntegerSubtangleEM array and bars.
	//As subtangles are joined, the rational array will be refined accordingly to removed joined indices.
	for(int i=0; i<2*numOfSubtangles; i++){
		gaussRationalSubtangleEM[i]=gaussIntegerSubtangleEM[i];
	}
	barsGaussRationalSubtangleEM[0]=barsGaussIntegerSubtangleEM[0];
	barsGaussRationalSubtangleEM[1]=barsGaussIntegerSubtangleEM[1];
	
	//Initialize an array ot specify the subtangles connected at the endpoint arcs and the connected corners of these subtangles.
	//These are calculated from the gaussIntegerSubtangleEM array and information about the bars.
	//The first column always refers to the first enter endpoint (NW by defualt); column 2 refers to the first exit endpoint; colum 3 to the second enter endpoint; and column 4 to the second exit endpoint.
	int endpointConnectionsCorners[2][4];
	
	//The first subtangle listed in the GaussEM array is always connected to the first enter endpoint arc, and always at corner A.
	//Likewise, the last subtangle in GaussEM array is always conneccted to the second exit endpoint arc, and always at corner D.
	endpointConnectionsCorners[0][0]=gaussIntegerSubtangleEM[0];
	endpointConnectionsCorners[1][0]=1;
	endpointConnectionsCorners[0][3]=gaussIntegerSubtangleEM[(2*numOfSubtangles)-1];
	endpointConnectionsCorners[1][3]=4;
	
	//The remaining two endpoints depend on whether this is the first or second time the adjacent subtangle has been encountered.
	int cornerExitIndex1=gaussIntegerSubtangleEM[barsGaussIntegerSubtangleEM[0]-1];
	int cornerEnterIndex2=gaussIntegerSubtangleEM[barsGaussIntegerSubtangleEM[0]];
	endpointConnectionsCorners[0][1]=cornerExitIndex1;
	endpointConnectionsCorners[0][2]=cornerEnterIndex2;
	
	//DEBUG:
	//printf("\n cornerExitIndex1 = %d \n cornerEnterIndex2 = %d \n", cornerExitIndex1,cornerEnterIndex2);
	//printf("\n endpointConnectionsCorners:\n\t S1-enter: %d \n\t S1-exit: %d \n\t S2-enter: %d \n\t s2-exit: %d \n", endpointConnectionsCorners[0][0],endpointConnectionsCorners[0][1],endpointConnectionsCorners[0][2],endpointConnectionsCorners[0][3]);
	
	int checkIndex;
	
	//Count the number of times cornerExitIndex1 shows up.
	checkIndex=0;
	for(int i=0; i<barsGaussIntegerSubtangleEM[0]; i++){
		if( gaussIntegerSubtangleEM[i] == cornerExitIndex1 ){
			checkIndex++;
		}
	}
	//If it only showed up once, the first exit endpoint connects at corner B; if twice, it connects at corner D.
	if( checkIndex < 2 ){
		endpointConnectionsCorners[1][1]=2;
	} else {
		endpointConnectionsCorners[1][1]=4;
	}
	
	//Repeat this count and check for the number of times cornerEnterIndex2 shows up.
	checkIndex=0;
	for(int i=0; i<=barsGaussIntegerSubtangleEM[0]; i++){
		if( gaussIntegerSubtangleEM[i] == cornerEnterIndex2 ){
			checkIndex++;
		}
	}
	//If it only showed up once, the second enter endpoint connects at corner A; if twice, it connects at corner C.
	if( checkIndex < 2 ){
		endpointConnectionsCorners[1][2]=1;
	} else {
		endpointConnectionsCorners[1][2]=3;
	}
	
	/*
	//DEBUG
	printf("\n endpointConnectionsEcorners: \n\t [");
	for(int i=0; i<4; i++){
		printf(" %d ", endpointConnectionsCorners[0][i]);
	}
	printf("] \n\t [");
	for(int i=0; i<4; i++){
		printf(" %d ", endpointConnectionsCorners[1][i]);
	}
	printf("]\n");
	*/
	
	//Initialize the rationalSubtangleConnectionsEM array as the integerSubtangleConnectionsEM array.
	//We continue to refine this array as more integer subtangles are joined into rational amalgams.
	for(int i=0; i<numOfSubtangles; i++){
		for(int j=0; j<5; j++){
			rationalSubtangleConnectionsEM[i][0][j] = integerSubtangleConnectionsEM[i][0][j];
			rationalSubtangleConnectionsEM[i][1][j] = integerSubtangleConnectionsEM[i][1][j];
		}
	}
	
	//Intialize a variable to track the number of refined rational tangles.
	//This starts out equal to the number of integer subtangles, but every time an integer subtangle gets joined into a rational amalgam, we reduce this count by 1.
	int rationalSubtanglesRefined=numOfSubtangles;
	
	//Define variables to track to the corners of the current rational amalgam. This is updated as more integer subtangles are joined.
	//These indicate the indicate the corresponding index of this corner under the generalized EM labeling (A, B, C, or D, which are represented as 1,2,4, or 4, respectivel).
	//We will follow the convention that corner A is always set to be NW.
	int cornerNW;
	int cornerNE;
	int cornerSE;
	int cornerSW;
	
	//Initialize variables for the types of twisting which could occur; update these later.
	int topTwist=0;
	int botTwist=0;
	int leftTwist=0;
	int rightTwist=0;
	
	//Define an array of length numOfSubtangles to track whether or not each integer subtangle has been joined into a rational amalgam.
	//Initialize each entry in this array as 1. If an integer tangle is joined to a rational amalgam, set the corresponding integerSubtangle entry to 0.
	//Only integer subtangles of entry 1 can freely be joined to rational subtangles; a 0 indicates it part of some other amalgam.
	//Set the entry corresponding to the initial integer subtangle in a rational amalgam to be the total number of joined integer subtangles.
	int integerSubtangleJoins[numOfSubtangles];
	for(int i=0; i<numOfSubtangles; i++){
		integerSubtangleJoins[i]=1;
	}
	
	//Initialize an array of twist vectors for each rational amalgam. These are built up while building the amalgam, and are used to recover a corresponding rational number.
	//There are at most numOfSubtangles twist vectors (rows), if each integer subtangle is a rational subtangle that cannot be built up any further.
	//Each twist vector has length at most numOfSubtangles, if there is a single twist vector consisting of all the integer subtangles (such as a rational tangle in canonical form).
	//Adopt the convention that the very first entry of every row denotes the length (number of nonzero entries) in the twist vector following it.
	//The sign of this very first entry will denote sign* of the crossings in the rational tangle as a whole (since a rational tangle of non-constant sign could be simplified).
	//This sign* is with respect to the usual orientation of the second strand of the subtangle, which may be reversed relative to the original tangle.
	int rationalTwistVectorArray[numOfSubtangles][2*MAXNN+1];
	int localTwistCount;
	int localInitialRationalIndex;
	
	//Set initial entries for the rationalSubtanleParametersEM array; these will be updated as amalgams are refined.
	for(int i=0; i<numOfSubtangles; i++){
		//0: crossing number and type:
		rationalSubtangleParametersEM[i][0]=abs(integerSubtangleParametersEM[i][0]);
		//1: number of integer subtangles, by default 1 at the beginning.
		rationalSubtangleParametersEM[i][1]=1;
		//2: twist vector length, by default 1 at the beginning.
		rationalSubtangleParametersEM[i][2]=1;
		//3: p; 4: q
		rationalSubtangleParametersEM[i][3]=integerSubtangleParametersEM[i][3];
		rationalSubtangleParametersEM[i][4]=integerSubtangleParametersEM[i][4];
		//5: internal parity:
		rationalSubtangleParametersEM[i][5]=integerSubtangleParametersEM[i][5];
		//6: S2 direction:
		rationalSubtangleParametersEM[i][6]=integerSubtangleParametersEM[i][6];
		//7: overall subtangle sign, as determined by shape: this is exactly the sign of q (integerSubtangleParametersEM[i][4]), regarding the integer tangle as a rational tangle.
		if( integerSubtangleParametersEM[i][4] < 0 ){
			rationalSubtangleParametersEM[i][7] = -1;
		} else {
			rationalSubtangleParametersEM[i][7] = 1;
		}
		//printf("\n\t DEBUG sign of rational subtangle shape: %d", rationalSubtangleParametersEM[i][7]);
		//8: The number of possible canonical configurations of the tangle? Still figuring out the best way to look at this info.
		rationalSubtangleParametersEM[i][8]=-1;
		
		//9: The type of twisting used in the initial integer subtangle, either horizontal (1) or vertical (-1). This is used later to build the twist vector.
		//By default, this is the same as whether the integer subtangle is horizontal or vertical, but it could be updated while building the rational amalgam.
		rationalSubtangleParametersEM[i][9]=integerSubtangleParametersEM[i][2];
		//One crossing subtangles will be regarded as horizontal.
		if( rationalSubtangleParametersEM[i][9] == 0 ){
			rationalSubtangleParametersEM[i][9]=1;
		}
	}
	
	/*
	//UPDATE (NC, 1/14/19) -- This entire block has been reworked below.
	//Initialze the twist vectors for each integer subtangle as though they were a rational subtangle; these will be updated as subtangles are joined.
	//Single crossing integer subtangles are by default set to a single length twist vector, with entry 1.
	//If an integer subtangle is internally regarded as horizontal twisting, it has initial length 1, with first entry the number of twists.
	//If an integer subtangle is internally regarded as vertical twisting, it has initial length 2, with first entry 1 and second entry the number of twists minus 1.
	//The above fact works thanks to certain structural properties of rational tanlges, specifically that the "first" twist can be regarded as either vertical or horizontal.
	//This convention establishes that the first twist is always regarded as horizontal in the twist vector.
	//The overall sign of the twist vector entries will be specified later, but it will depend on internal crossing sign of the initial integer subtangle and whether S2 is reversed or not.
	for(int i=0; i<numOfSubtangles;i++){
		if( integerSubtangleParametersEM[i][2] == 0 ){
			//Case 1: 1 crossing integer subtangle
			rationalTwistVectorArray[i][0]=1;
			rationalTwistVectorArray[i][1]=1;
		} else if( integerSubtangleParametersEM[i][2] == 1 ){
			//Case 2: horizontal integer subtangle
			rationalTwistVectorArray[i][0]=1;
			rationalTwistVectorArray[i][1]=abs(integerSubtangleParametersEM[i][0]);
		} else if( integerSubtangleParametersEM[i][2] == -1 ){
			//Case 3: vertical integer subtangle
			rationalTwistVectorArray[i][0]=2;
			rationalTwistVectorArray[i][1]=1;
			rationalTwistVectorArray[i][2]=abs(integerSubtangleParametersEM[i][0])-1;
			//Also increment this length in the corresponding entry of the rationalSubtangleParametersEM array.
			rationalSubtangleParametersEM[i][2]++;
		} else {
			//Error Case: this shouldn't happen, but having a garbage value here will explode the program when iterating up to the length, so set the length equal to 0.
			rationalTwistVectorArray[i][0]=0;
			rationalSubtangleParametersEM[i][2]=0;
		}
	}
	*/
	
	//The twist type (horizontal or vertical) of the initial integer tangle is now being tracked in rationalSubtangleParametersEM[i][9]; this potentially changes while building up the subtangle.
	//By default, we will ignore modifying the twist vector to distinguish these until after the entire rational subtnangle is finished; at the start, all twist vectors have length 1.
	for(int i=0; i<numOfSubtangles; i++){
		rationalTwistVectorArray[i][0]=1;
		rationalTwistVectorArray[i][1]=abs(integerSubtangleParametersEM[i][0]);
	}
	
	
	/*
	//DEBUG
	printRationalSubtangleParametersEM(rationalSubtanglesRefined,rationalSubtangleParametersEM);
	printArrayEM(rationalSubtanglesRefined,rationalSubtangleConnectionsEM);
	
	//DEBUG:
	printf("\n Twist Vectors: \n");
	for(int i=0; i<numOfSubtangles; i++){
		printf(" [ %d ]\t[", i+1);
		for(int j=1; j<=rationalTwistVectorArray[i][0]; j++){
			printf(" %d ", rationalTwistVectorArray[i][j]);
		}
		printf("]\n");
	}
	*/
	
	//Initialize a storage array; this is used to help match the canonical configurations with a rotated NW label in some special cases.
	int oldRationalComponentsCanonicalConfiguration[5];
	
	for(int i=0; i<numOfSubtangles; i++){
		
		//This array interates through the original list of integer subtangles and checks to see if each is the intial tangle is some rational amalgam.
		//It is possible that the current integer tangle might already have been joined to some other rational amalgam in a previous step, in which case we don't need to consider it separately.
		//Check to see if this is true before looking for further joinings; if it has already been joined, nothing else happens and we move to the next iteration of the for loop.
		if ( integerSubtangleJoins[i] == 1 ){
			
			//By defualt, any choice of NW yields a canonical configuration for an integer tangle.
			//Entry rationalComponentsCanonicalConfiguration[i][0] denotes the number of integer components, which is 1 in this case.
			for(int j=0; j<5; j++){
				rationalComponentsCanonicalConfiguration[i][j] = 1;
			}
			
			//We look to see if this integer subtangle is the initial subtangle in some rational amalgam.
			//Reset the local twist count; use this to track how many integer subtangles get joined.
			localInitialRationalIndex=i;
			localTwistCount=0;
			
			//For everything: NW = A (1)
			cornerNW=1;
			
			//The index of the remaining corners requires case considerations.
			//These depend on the (internal) parity of the current subtangle (integer or rational amalgam) and whether or not S2 is reversed.
			if( rationalSubtangleParametersEM[i][5] == 0 ){
				//Parity 0 case: NE = B (2)
				cornerNE=2;
				if( rationalSubtangleParametersEM[i][6] > 0 ){
					//If S2 is not reversed, then SE = C (3) and SW = D (4).
					cornerSE=3;
					cornerSW=4;
				} else {
					//If S2 is reversed, then SE = D (4) and SW = C (3)
					cornerSE=4;
					cornerSW=3;
				}
			} else if( rationalSubtangleParametersEM[i][5] == 1 ){
				//Parity 1 case: SE = B (2)
				cornerSE=2;
				if( rationalSubtangleParametersEM[i][6] > 0 ){
					//If S2 is not reversed, then NE = C (3) and SW = D (4).
					cornerNE=3;
					cornerSW=4;
				} else {
					//If S2 is reversed, then NE = D (4) and SW = C (3)
					cornerNE=4;
					cornerSW=3;
				}
			} else if( rationalSubtangleParametersEM[i][5] == 2 ){
				//Parity infinity (2) case: SW = B (2)
				cornerSW=2;
				if( rationalSubtangleParametersEM[i][6] > 0 ){
					//If S2 is not reversed, then SE = C (3) and NW = D (4).
					cornerSE=3;
					cornerNE=4;
				} else {
					//If S2 is reversed, then SE = D (4) and NW = C (3)
					cornerSE=4;
					cornerNE=3;
				}
			}
			//Tracking the corners in this way allows for a more elegant check while looking for twisting.
			
			//DEBUG
			//printf("\n DEBUG FLAG: parity = %d \n", rationalSubtangleParametersEM[i][5]);
			
			//If any twisting happens, it will be detected by the detectCronerTwisting() function, which returns the number of twists.
			rightTwist=detectCornerTwisting(i,cornerNE,cornerSE,rationalSubtangleConnectionsEM,integerSubtangleJoins,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
			botTwist=detectCornerTwisting(i,cornerSW,cornerSE,rationalSubtangleConnectionsEM,integerSubtangleJoins,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
			leftTwist=detectCornerTwisting(i,cornerNW,cornerSW,rationalSubtangleConnectionsEM,integerSubtangleJoins,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
			topTwist=detectCornerTwisting(i,cornerNW,cornerNE,rationalSubtangleConnectionsEM,integerSubtangleJoins,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
			
			//DEBUG
			//printf("\n DEBUG TWIST CHECKS--BEFORE CHECKING ENDPOINTS: Subtangle %d \n rightTwist = %d \n botTwist = %d \n leftTwist = %d \n topTwist = %d \n",i+1,rightTwist,botTwist,leftTwist,topTwist);
			
			//In some cases, subtangles which are joined across endpoint arcs in the tangle closure might satisfy twist check conditions.
			//We do not want to classify these as twists; we compare against the information in the endpointConnectionsCorners array.
			//Unfortunately, this check requires a number of case considerations and feels generally inelegant, which is weird since this seems like an uncommon fringe case.
			//If we find that any such conditions are satisfied for one of the four types of twisting, we reset that twist to 0.
			for(int j=0; j<4; j++){
				if( (endpointConnectionsCorners[0][j]-1) == localInitialRationalIndex ){
					//Right Twists:
					if( (endpointConnectionsCorners[1][j]==cornerNE) || (endpointConnectionsCorners[1][j]==cornerSE) ){
						rightTwist=0;
					}
					//Bottom Twists:
					if( (endpointConnectionsCorners[1][j]==cornerSW) || (endpointConnectionsCorners[1][j]==cornerSE) ){
						botTwist=0;
					}
					//Left Twists:
					if( (endpointConnectionsCorners[1][j]==cornerNW) || (endpointConnectionsCorners[1][j]==cornerSW) ){
						leftTwist=0;
					}
					//Top Twists:
					if( (endpointConnectionsCorners[1][j]==cornerNW) || (endpointConnectionsCorners[1][j]==cornerNE) ){
						topTwist=0;
					}
				}
			}
			//If the current subtangle showed up as one of the endpoint tangles, any twisting involving the corner at that endpoint is set to 0.
			//Note that, if an endpoint subtangle gets joined to a rational amalgam, the index of that tangle in the endpointConnectionsCorners array will need to be updated.
			
			//DEBUG
			//printf("\n DEBUG TWIST CHECKS--AFTER CHECKING ENDPOINTS: Subtangle %d \n rightTwist = %d \n botTwist = %d \n leftTwist = %d \n topTwist = %d \n",i+1,rightTwist,botTwist,leftTwist,topTwist);
			
			//If there is any tiwsting, update the configuration of the component as needed.
			if( (rightTwist+botTwist+leftTwist+topTwist) > 0 ){
				//If there are both left and right twists, no configuaration is canonical. Set every configuration to 0, but indicate that 2 new integer components are added.
				if( (rightTwist>0) && (leftTwist>0) ){
					rationalComponentsCanonicalConfiguration[i][0]+=2;
					for(int j=1; j<5; j++){
						rationalComponentsCanonicalConfiguration[i][j] = 0;
					}
				}
				//Do the same as above if there are both top and bot twists.
				if( (topTwist>0) && (botTwist>0) ){
					rationalComponentsCanonicalConfiguration[i][0]+=2;
					for(int j=1; j<5; j++){
						rationalComponentsCanonicalConfiguration[i][j] = 0;
					}
				}
				
				//If there is twisting, but we avoid the two cases above, then there is twisting in one and only one direction.
				if( rightTwist > 0 ){
					rationalComponentsCanonicalConfiguration[i][0]++;
					//If right twisting, invaldiate canonical configurations 2 and 3;
					rationalComponentsCanonicalConfiguration[i][2] = 0;
					rationalComponentsCanonicalConfiguration[i][3] = 0;
				} else if( botTwist > 0 ){
					rationalComponentsCanonicalConfiguration[i][0]++;
					//If bot twisting, invaldiate canonical configurations 3 and 4;
					rationalComponentsCanonicalConfiguration[i][3] = 0;
					rationalComponentsCanonicalConfiguration[i][4] = 0;
				} else if( leftTwist > 0 ){
					rationalComponentsCanonicalConfiguration[i][0]++;
					//If left twisting, invaldiate canonical configurations 1 and 4;
					rationalComponentsCanonicalConfiguration[i][1] = 0;
					rationalComponentsCanonicalConfiguration[i][4] = 0;
					//If the number of left twists is odd, corner A moves to the existing SW position, and then the labels rotate 90 degrees so that this is regarded as NW again.
					if( (leftTwist%2) == 1 ){
						for(int j=1; j<5; j++){
							oldRationalComponentsCanonicalConfiguration[j] = rationalComponentsCanonicalConfiguration[i][j];
						}
						rationalComponentsCanonicalConfiguration[i][1] = oldRationalComponentsCanonicalConfiguration[4];
						rationalComponentsCanonicalConfiguration[i][2] = oldRationalComponentsCanonicalConfiguration[1];
						rationalComponentsCanonicalConfiguration[i][3] = oldRationalComponentsCanonicalConfiguration[2];
						rationalComponentsCanonicalConfiguration[i][4] = oldRationalComponentsCanonicalConfiguration[3];
					}
				} else if( topTwist > 0 ){
					rationalComponentsCanonicalConfiguration[i][0]++;
					//If top twisting, invaldiate canonical configurations 1 and 2;
					rationalComponentsCanonicalConfiguration[i][1] = 0;
					rationalComponentsCanonicalConfiguration[i][2] = 0;
					//If the number of top twists is odd, corner A moves to the existing NE position, and then the labels rotate 90 degrees so that this is regarded as NW again.
					if( (topTwist%2) == 1 ){
						for(int j=1; j<5; j++){
							oldRationalComponentsCanonicalConfiguration[j] = rationalComponentsCanonicalConfiguration[i][j];
						}
						rationalComponentsCanonicalConfiguration[i][1] = oldRationalComponentsCanonicalConfiguration[2];
						rationalComponentsCanonicalConfiguration[i][2] = oldRationalComponentsCanonicalConfiguration[3];
						rationalComponentsCanonicalConfiguration[i][3] = oldRationalComponentsCanonicalConfiguration[4];
						rationalComponentsCanonicalConfiguration[i][4] = oldRationalComponentsCanonicalConfiguration[1];
					}
				}
				
				//DEBUG
				//printf("\n DEBUG FLAG: line 3163:\n rationalComponentsCanonicalConfiguration: [%d] [ %d | %d %d %d %d ]", i+1, rationalComponentsCanonicalConfiguration[i][0], rationalComponentsCanonicalConfiguration[i][1], rationalComponentsCanonicalConfiguration[i][2], rationalComponentsCanonicalConfiguration[i][3], rationalComponentsCanonicalConfiguration[i][4]);
			}
			
			
			joinIntegerToRational(rightTwist,botTwist,leftTwist,topTwist,cornerNW,cornerNE,cornerSE,cornerSW,localTwistCount,rationalTwistVectorArray,integerSubtangleJoins,localInitialRationalIndex,rationalSubtangleConnectionsEM,rationalSubtangleParametersEM,integerSubtangleParametersEM,rationalSubtanglesRefined,numOfSubtangles,endpointConnectionsCorners,gaussRationalSubtangleEM,barsGaussRationalSubtangleEM);
			
			//DEBUG
			//printf("\n DEBUG-- rationalSubtangleConnectionsEM array after FIRST joining:");
			//printArrayEM(numOfSubtangles,rationalSubtangleConnectionsEM);
			
			//Continue to search for possible integer subtangles that can be joined to the growing rational amalgam.
			while( localTwistCount > 0 ){
				
				//Reset the local twist count variable every iteration. If more twists are added, this gets updated again.
				localTwistCount=0;
				
				//Now that connections in the rationalSubtangleConnectionsEM array have been updated into a larger rational amalgam, we want to check if any further integer subtangles can be joined.
				//The NW/NE/SE/SW corner information is updated while joining things in the joinIntegerToRational() function.
				//If any twisting happens, it will be detected by the detectCronerTwisting() function, which returns the number of twists.
				rightTwist=detectCornerTwisting(i,cornerNE,cornerSE,rationalSubtangleConnectionsEM,integerSubtangleJoins,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
				botTwist=detectCornerTwisting(i,cornerSW,cornerSE,rationalSubtangleConnectionsEM,integerSubtangleJoins,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
				leftTwist=detectCornerTwisting(i,cornerNW,cornerSW,rationalSubtangleConnectionsEM,integerSubtangleJoins,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
				topTwist=detectCornerTwisting(i,cornerNW,cornerNE,rationalSubtangleConnectionsEM,integerSubtangleJoins,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
				
				//DEBUG
				/*
				printf("\n Check for MORE twisting--BEFORE CHECKING ENDPOINTS: Subtangle %d \n rightTwist = %d \n botTwist = %d \n leftTwist = %d \n topTwist = %d \n",i+1,rightTwist,botTwist,leftTwist,topTwist);
				printf("\n\t CORNERS: NW = %d , NE = %d , SE = %d , SW = %d \n", cornerNW, cornerNE, cornerSE, cornerSW);
				printArrayEM(numOfSubtangles,rationalSubtangleConnectionsEM);
				*/
				
				//The endpointConnectionsCorners array should have been updated after the previous joining.
				//Check again against the entries in this array to rule out the possibility of endpoint twisting.
				for(int j=0; j<4; j++){
					if( (endpointConnectionsCorners[0][j]-1) == localInitialRationalIndex ){
						//Right Twists:
						if( (endpointConnectionsCorners[1][j]==cornerNE) || (endpointConnectionsCorners[1][j]==cornerSE) ){
							rightTwist=0;
						}
						//Bottom Twists:
						if( (endpointConnectionsCorners[1][j]==cornerSW) || (endpointConnectionsCorners[1][j]==cornerSE) ){
							botTwist=0;
						}
						//Left Twists:
						if( (endpointConnectionsCorners[1][j]==cornerNW) || (endpointConnectionsCorners[1][j]==cornerSW) ){
							leftTwist=0;
						}
						//Top Twists:
						if( (endpointConnectionsCorners[1][j]==cornerNW) || (endpointConnectionsCorners[1][j]==cornerNE) ){
							topTwist=0;
						}
					}
				}
				
				//DEBUG
				//printf("\n Check for MORE twisting--AFTER CHECKING ENDPOINTS: Subtangle %d \n rightTwist = %d \n botTwist = %d \n leftTwist = %d \n topTwist = %d \n",i+1,rightTwist,botTwist,leftTwist,topTwist);
				
				//If there is any tiwsting, update the configuration of the component as needed.
				if( (rightTwist+botTwist+leftTwist+topTwist) > 0 ){
					//If there are both left and right twists, no configuaration is canonical. Set every configuration to 0, but indicate that 2 new integer components are added.
					if( (rightTwist>0) && (leftTwist>0) ){
						rationalComponentsCanonicalConfiguration[i][0]+=2;
						for(int j=1; j<5; j++){
							rationalComponentsCanonicalConfiguration[i][j] = 0;
						}
					}
					//Do the same as above if there are both top and bot twists.
					if( (topTwist>0) && (botTwist>0) ){
						rationalComponentsCanonicalConfiguration[i][0]+=2;
						for(int j=1; j<5; j++){
							rationalComponentsCanonicalConfiguration[i][j] = 0;
						}
					}
					
					if( rightTwist > 0 ){
						rationalComponentsCanonicalConfiguration[i][0]++;
						//If right twisting, invaldiate canonical configurations 2 and 3;
						rationalComponentsCanonicalConfiguration[i][2] = 0;
						rationalComponentsCanonicalConfiguration[i][3] = 0;
					} else if( botTwist > 0 ){
						rationalComponentsCanonicalConfiguration[i][0]++;
						//If bot twisting, invaldiate canonical configurations 3 and 4;
						rationalComponentsCanonicalConfiguration[i][3] = 0;
						rationalComponentsCanonicalConfiguration[i][4] = 0;
					} else if( leftTwist > 0 ){
						rationalComponentsCanonicalConfiguration[i][0]++;
						//If left twisting, invaldiate canonical configurations 1 and 4;
						rationalComponentsCanonicalConfiguration[i][1] = 0;
						rationalComponentsCanonicalConfiguration[i][4] = 0;
						//If the number of left twists is odd, corner A moves to the existing SW position, and then the labels rotate 90 degrees so that this is regarded as NW again.
						if( (leftTwist%2) == 1 ){
							for(int j=1; j<5; j++){
								oldRationalComponentsCanonicalConfiguration[j] = rationalComponentsCanonicalConfiguration[i][j];
							}
							rationalComponentsCanonicalConfiguration[i][1] = oldRationalComponentsCanonicalConfiguration[4];
							rationalComponentsCanonicalConfiguration[i][2] = oldRationalComponentsCanonicalConfiguration[1];
							rationalComponentsCanonicalConfiguration[i][3] = oldRationalComponentsCanonicalConfiguration[2];
							rationalComponentsCanonicalConfiguration[i][4] = oldRationalComponentsCanonicalConfiguration[3];
						}
					} else if( topTwist > 0 ){
						rationalComponentsCanonicalConfiguration[i][0]++;
						//If top twisting, invaldiate canonical configurations 1 and 2;
						rationalComponentsCanonicalConfiguration[i][1] = 0;
						rationalComponentsCanonicalConfiguration[i][2] = 0;
						//If the number of top twists is odd, corner A moves to the existing NE position, and then the labels rotate 90 degrees so that this is regarded as NW again.
						if( (topTwist%2) == 1 ){
							for(int j=1; j<5; j++){
								oldRationalComponentsCanonicalConfiguration[j] = rationalComponentsCanonicalConfiguration[i][j];
							}
							rationalComponentsCanonicalConfiguration[i][1] = oldRationalComponentsCanonicalConfiguration[2];
							rationalComponentsCanonicalConfiguration[i][2] = oldRationalComponentsCanonicalConfiguration[3];
							rationalComponentsCanonicalConfiguration[i][3] = oldRationalComponentsCanonicalConfiguration[4];
							rationalComponentsCanonicalConfiguration[i][4] = oldRationalComponentsCanonicalConfiguration[1];
						}
					}
					//DEBUG
					//printf("\n DEBUG FLAG: line 3279\n rationalComponentsCanonicalConfiguration: [%d] [ %d | %d %d %d %d ]", i+1, rationalComponentsCanonicalConfiguration[i][0], rationalComponentsCanonicalConfiguration[i][1], rationalComponentsCanonicalConfiguration[i][2], rationalComponentsCanonicalConfiguration[i][3], rationalComponentsCanonicalConfiguration[i][4]);					
				}
				
				
				//Attempt to iteratively join possible integer subtangles to the amalgam until no more twists are detected.
				joinIntegerToRational(rightTwist,botTwist,leftTwist,topTwist,cornerNW,cornerNE,cornerSE,cornerSW,localTwistCount,rationalTwistVectorArray,integerSubtangleJoins,localInitialRationalIndex,rationalSubtangleConnectionsEM,rationalSubtangleParametersEM,integerSubtangleParametersEM,rationalSubtanglesRefined,numOfSubtangles,endpointConnectionsCorners,gaussRationalSubtangleEM,barsGaussRationalSubtangleEM);
				
			}
			
			//There are 10 possible cases denoting which configurations of a given rational component are canonical.
			if( (rationalComponentsCanonicalConfiguration[i][1]==1) && (rationalComponentsCanonicalConfiguration[i][2]==1) && (rationalComponentsCanonicalConfiguration[i][3]==1) && (rationalComponentsCanonicalConfiguration[i][4]==1) ){
				//CASE -1: [1,1,1,1] (4 canonical configurations, integer tangle)
				rationalSubtangleParametersEM[i][8] = -1;
			} else if( (rationalComponentsCanonicalConfiguration[i][1]==0) && (rationalComponentsCanonicalConfiguration[i][2]==0) && (rationalComponentsCanonicalConfiguration[i][3]==0) && (rationalComponentsCanonicalConfiguration[i][4]==0) ){
				//CASE 0: [0,0,0,0] (no canonical configurations)
				rationalSubtangleParametersEM[i][8] = 0;
			} else if( (rationalComponentsCanonicalConfiguration[i][1]==1) && (rationalComponentsCanonicalConfiguration[i][2]==0) && (rationalComponentsCanonicalConfiguration[i][3]==0) && (rationalComponentsCanonicalConfiguration[i][4]==0) ){
				//CASE 1: [1,0,0,0] (only 1 canonical configuration, with NW=A)
				rationalSubtangleParametersEM[i][8] = 1;
			} else if( (rationalComponentsCanonicalConfiguration[i][1]==0) && (rationalComponentsCanonicalConfiguration[i][2]==1) && (rationalComponentsCanonicalConfiguration[i][3]==0) && (rationalComponentsCanonicalConfiguration[i][4]==0) ){
				//CASE 2: [0,1,0,0] (only 1 canonical configuration, with NW=A*90deg)
				rationalSubtangleParametersEM[i][8] = 2;
			} else if( (rationalComponentsCanonicalConfiguration[i][1]==0) && (rationalComponentsCanonicalConfiguration[i][2]==0) && (rationalComponentsCanonicalConfiguration[i][3]==1) && (rationalComponentsCanonicalConfiguration[i][4]==0) ){
				//CASE 3: [0,0,1,0] (only 1 canonical configuration, with NW=A*180deg)
				rationalSubtangleParametersEM[i][8] = 3;
			} else if( (rationalComponentsCanonicalConfiguration[i][1]==0) && (rationalComponentsCanonicalConfiguration[i][2]==0) && (rationalComponentsCanonicalConfiguration[i][3]==0) && (rationalComponentsCanonicalConfiguration[i][4]==1) ){
				//CASE 4: [0,0,0,1] (only 1 canonical configuration, with NW=A*270deg)
				rationalSubtangleParametersEM[i][8] = 4;
			} else if( (rationalComponentsCanonicalConfiguration[i][1]==1) && (rationalComponentsCanonicalConfiguration[i][2]==1) && (rationalComponentsCanonicalConfiguration[i][3]==0) && (rationalComponentsCanonicalConfiguration[i][4]==0) ){
				//CASE 5: [1,1,0,0] (2 canonical configurations, with NW=A or A*90deg)
				rationalSubtangleParametersEM[i][8] = 5;
			} else if( (rationalComponentsCanonicalConfiguration[i][1]==0) && (rationalComponentsCanonicalConfiguration[i][2]==1) && (rationalComponentsCanonicalConfiguration[i][3]==1) && (rationalComponentsCanonicalConfiguration[i][4]==0) ){
				//CASE 6: [0,1,1,0] (2 canonical configurations, with NW=A*90deg or A*180deg)
				rationalSubtangleParametersEM[i][8] = 6;
			} else if( (rationalComponentsCanonicalConfiguration[i][1]==0) && (rationalComponentsCanonicalConfiguration[i][2]==0) && (rationalComponentsCanonicalConfiguration[i][3]==1) && (rationalComponentsCanonicalConfiguration[i][4]==1) ){
				//CASE 7: [0,0,1,1] (2 canonical configurations, with NW=A*180deg or A*270deg)
				rationalSubtangleParametersEM[i][8] = 7;
			} else if( (rationalComponentsCanonicalConfiguration[i][1]==1) && (rationalComponentsCanonicalConfiguration[i][2]==0) && (rationalComponentsCanonicalConfiguration[i][3]==0) && (rationalComponentsCanonicalConfiguration[i][4]==1) ){
				//CASE 8: [1,0,0,1] (2 canonical configurations, with NW=A*270deg or A)
				rationalSubtangleParametersEM[i][8] = 8;
			}
		}
		
	}
	
	//DEBUG:
	//TO BE RE-INDEXED AT SOME POINT TO ACCOUNT FOR JOINED SUBTANGLES?
	//printf("\n DEBUG FLAG: line %d: rational subtangle parameters BEFORE re-index", __LINE__);
	//printRationalSubtangleParametersEM(numOfSubtangles,rationalSubtangleParametersEM);
	//printArrayEM(numOfSubtangles,rationalSubtangleConnectionsEM);
	
	//After the joining is completed, we must finalize the twist vector.
	//If the initial integer subtangle in a rational amalgam is regarded as vertical instead of horizontal, we must modify the twist vector slightly to start on a horizontal twist.
	for(int i=0; i<numOfSubtangles; i++){
		if( rationalSubtangleParametersEM[i][9] == -1 ){
			//Increase the length of the twist vector by 1.
			rationalTwistVectorArray[i][0]++;
			rationalSubtangleParametersEM[i][2]++;
			//Working backwards, shift all twist vector entries one position over to the right.
			for(int j=rationalTwistVectorArray[i][0]; j>1; j--){
				rationalTwistVectorArray[i][j]=rationalTwistVectorArray[i][j-1];
			}
			//Lastly, set the very first twist vector entry to 1, and reduce the second twist vector entry by 1.
			rationalTwistVectorArray[i][1]=1;
			rationalTwistVectorArray[i][2]--;
		}
	}
	
	//Next, the sign must be specified on the entries in the twist vectors.
	//We assume that each twist vector entry must have the same sign (else, the corresponding rational tangle could be reduced).
	//Hence, the sign of each entry will match the twist sign* of the crossings in the initial integer subtangle in the amalgam.
	//In this case, we mean the sign* of the shape of the tangle (that is, the twists), assuming S2 were not reversed, which is possible.
	//This information is tracked in entry 7 of the rationalSubtangleParametersEM array.
	//Note that this was originally determined by the sign of q from the initial integer subtanlge (integerSubtangleParametersEM[i][4]), but this sign can change while joining to the rational amalgam.
	for(int i=0; i<numOfSubtangles; i++){
		if( rationalSubtangleParametersEM[i][7] < 0 ){
			for(int j=1; j<=rationalTwistVectorArray[i][0]; j++){
				rationalTwistVectorArray[i][j] = -1*rationalTwistVectorArray[i][j];
			}
		}
	}
	//The twist vector convention requires that the last entry denote a horizontal twist.
	//Since an earlier choice forces the first entry in the twist vector to denote a horziontal twist as well, permissible twist vectors must have an odd number of entries.
	//If one of the twist vectors obtained has an even number of entries, the last entry must be vertical (since the twists alternate between horizontal and vertical).
	//In this case, we may append a single entry of 0 to the end of the twist vector to denote a horizontal twist, which forces the twist vector to satisfy the above requirements.
	for(int i=0; i<numOfSubtangles; i++){
		if( (rationalTwistVectorArray[i][0]%2) == 0 ){
			//The 0 entrey of the rationalTwistVectorArray[i] is the length of the twist vector; if this length is even, we increase it by 1 and append an entry of 0.
			rationalTwistVectorArray[i][0]++;
			rationalSubtangleParametersEM[i][2]++;
			rationalTwistVectorArray[i][rationalTwistVectorArray[i][0]]=0;
		}
	} 
	
	/*
	//DEBUG
	printf("\n Twist Vectors: \n");
	for(int i=0; i<numOfSubtangles; i++){
		printf(" [ %d ]\t[", i+1);
		for(int j=1; j<=rationalTwistVectorArray[i][0]; j++){
			printf(" %d ", rationalTwistVectorArray[i][j]);
		}
		printf("]\n");
	}
	
	//DEBUG:
	//WHEN RE-LABELING INDEXING: collapse the index of anything with integerSubtangleJoins entry == 0. Do this at the VERY end of the joining.
	printf("\n Integer Subtangle Joins: [");
	for(int i=0; i<numOfSubtangles; i++){
		printf(" %d ", integerSubtangleJoins[i]);
	}
	printf("]\n");
	*/
	
	//Now use the twist vectors obtained earlier to compute the corresponding rational number p/q for each of the rational subtangles.
	//These will be stored in the rationalSubtangleParametersEM array.
	//AFTER SETTING UP THE RE-LABELING INDEX THING, CHANGE THE UPPER BOUND ON ITERATION TO rationalSubtanglesRefined RATHER THAN numOFSubtangles
	int rationalPQ[2];
	int localTwistVector[30];
	for(int i=0; i<numOfSubtangles; i++){
		//Setup a local twist vector to be consistent with the required notation of the function that calculates the extended fraction.
		for(int j=1; j<=rationalTwistVectorArray[i][0]; j++){
			localTwistVector[j-1]=rationalTwistVectorArray[i][j];
		}
		EMcalculateExtendedFractionPQ(true,rationalTwistVectorArray[i][0],localTwistVector,rationalPQ);
		//Set the corresponding entries of the rationalSubtangleParametersEM array eqaul to the p and q found by the function.
		rationalSubtangleParametersEM[i][3]=rationalPQ[0];
		rationalSubtangleParametersEM[i][4]=rationalPQ[1];
	}
	
	//DEBUG
	//printRationalSubtangleParametersEM(numOfSubtangles,rationalSubtangleParametersEM);
	
	//DEBUG:
	/*
	printf("\n\n DEBUG line 3396: rationalSubtangleConnectionsEM array BEFORE removing dead:\n");
	printArrayEM(numOfSubtangles,rationalSubtangleConnectionsEM);
	printRationalSubtangleParametersEM(numOfSubtangles,rationalSubtangleParametersEM);
	*/
	
	
	//DEBUG
	/*
	printf("\n Rational GaussEM Code: [");
	for(int i=0; i<2*rationalSubtanglesRefined; i++){
		printf(" %d ",gaussRationalSubtangleEM[i]);
	}
	printf("]");
	*/
	
	//Now go through and reindex the rationalSubtangleConnectionsEM array and rationalSubtangleParametersEM array to eliminate joined tangles.
	int joinedIntegerSubtangles=0;
	for(int i=0; i<numOfSubtangles; i++){
		if( integerSubtangleJoins[i] == 0 ){
			//In this case, the subtangle is fully joined elsewhere and we can safely eliminate it.
			joinedIntegerSubtangles++;
		} else {
			//If this tangle is not eliminated, move its entries in the array up by the number of joined subtangles.
			//This way, only the joined subtangle portions of the array are overwitten; everything else just shifts up.
			for(int j=0; j<5; j++){
				rationalSubtangleConnectionsEM[i-joinedIntegerSubtangles][0][j] = rationalSubtangleConnectionsEM[i][0][j];
				rationalSubtangleConnectionsEM[i-joinedIntegerSubtangles][1][j] = rationalSubtangleConnectionsEM[i][1][j];
			}
			//Do the same thing for the rationalSubtangleParametersEM array.
			for(int j=0; j<10; j++){
				rationalSubtangleParametersEM[i-joinedIntegerSubtangles][j] = rationalSubtangleParametersEM[i][j];
			}
			//Do it for the twist vectors too, while we're at it.
			for(int j=0; j<=rationalTwistVectorArray[i][0]; j++){
				rationalTwistVectorArray[i-joinedIntegerSubtangles][j] = rationalTwistVectorArray[i][j];
			}
		}
	}
	
	//DEBUG:
	/*
	printf("\n\n DEBUG line 3437: COLLAPSED rationalSubtangleConnectionsEM array AFTER removing dead labels:\n");
	printArrayEM(rationalSubtanglesRefined,rationalSubtangleConnectionsEM);
	printRationalSubtangleParametersEM(rationalSubtanglesRefined,rationalSubtangleParametersEM);
	*/
	/*
	//DEBUG
	printf("\n Twist Vectors: \n");
	for(int i=0; i<rationalSubtanglesRefined; i++){
		printf(" [ %d ]\t[", i+1);
		for(int j=1; j<=rationalTwistVectorArray[i][0]; j++){
			printf(" %d ", rationalTwistVectorArray[i][j]);
		}
		printf("]\n");
	}
	*/
	//DEBUG
	//printPrettyRationalSubtangleEM(rationalSubtanglesRefined,rationalSubtangleConnectionsEM,rationalSubtangleParametersEM,rationalTwistVectorArray);
	
	
	//The subtangles in the rational Gauss EM code are currently indexed by the integer subtangle index of the corresponding initial subtangle.
	//We would like to reindex these to be increasing on subtangles, without skipping indices.
	//Initialize a variable to track the number of skipped indices.
	int skippedIndex = 0;
	int referenceIndex[rationalSubtanglesRefined] = {0};
	int failsafe=0;
	for(int i=1; i<=rationalSubtanglesRefined; i++){
		//Use a while-loop to to continue searching until all necessary indices are referenced.
		while( (referenceIndex[i-1] == 0) && (failsafe<(2*numOfSubtangles)) ){
			//For each subtangle index, iterate through all entries of the gaussRationalSubtangleEM[] array and update references to callapse missing indices.
			for(int j=0; j<(2*rationalSubtanglesRefined); j++){
				//Each time skipped indices are identified, shift the entry magnitudes down by the number of skips.
				if( gaussRationalSubtangleEM[j] == (i+skippedIndex) ){
					referenceIndex[i-1]++;
					gaussRationalSubtangleEM[j] = i;
				}
			}
			//If there were references identified, then we would also like to update these in the rationalSubtangleConnectionsEM[][] array.
			//We would also like to update these in the endpointConnectionsCorners[][] array.
			if( referenceIndex[i-1] > 0 ){
				//Update rationalSubtangleConnectionsEM[][]:
				for(int j=0; j<rationalSubtanglesRefined; j++){
					for(int k=1; k<5; k++){
						if( rationalSubtangleConnectionsEM[j][0][k] == (i+skippedIndex) ){
							rationalSubtangleConnectionsEM[j][0][k] = i;
						}
					}
				}
				//Update endpointConnectionsCorners[][]:
				for(int j=0; j<4; j++){
					if( endpointConnectionsCorners[0][j] == (i+skippedIndex) ){
						endpointConnectionsCorners[0][j] = i;
					}
				}
			} else {
				//However, if no references were identified, iterate the number of skips. Eventually, this process must terminate.
				skippedIndex++;
			}
			failsafe++;
		}
		//DEBUG
		//printf("\n LOOP %d , failsafe = %d ", i, failsafe);
	}
	//DEBUG
	/*
	printf("\n Rational GaussEM Code: [");
	for(int i=0; i<2*rationalSubtanglesRefined; i++){
		printf(" %d ",gaussRationalSubtangleEM[i]);
	}
	printf("]");
	*/
	
	//After all of the rational subtangles have been determined and the indexing on everything updated, set numOfSubtangles equal to rationalSubtanglesRefined, the number of rational subtangles.
	//Also set the rationalTwistVectorsEM array to the locally defined rationalTwistVectorArray.
	numOfSubtangles=rationalSubtanglesRefined;
	for(int i=0; i<rationalSubtanglesRefined; i++){
		for(int j=0; j<=rationalTwistVectorArray[i][0]; j++){
			rationalTwistVectorsEM[i][j] = rationalTwistVectorArray[i][j];
		}
	}
	//DEBUG
	//printPrettyRationalSubtangleEM(numOfSubtangles,rationalSubtangleConnectionsEM,rationalSubtangleParametersEM,rationalTwistVectorsEM);
	
	//Populate the entries in the endpointConnectionsCornersRational[][] array with the endpoint connections information.
	for(int i=0; i<4; i++){
		endpointConnectionsCornersRational[0][i]=endpointConnectionsCorners[0][i];
		endpointConnectionsCornersRational[1][i]=endpointConnectionsCorners[1][i];
	}
	
}






/*
(NC, 1/21/19)
This function will convert the subtangle connections array from the enter/exit labeling (A/B/C/D) to the compass labeling (NW/NE/SE/SW).
We assume that the endpoint connections in the input array are already identified.
This was written for use converting the labeling on rational components, but this could be modified for use elsewhere.

Inputs:
numOfSubtanglesInput: the number of component subtangles in the connections arrays.
connectionsEnterExitInput[][2][5]: the array storing the EM connections info using the enter/exit labeling; this includes identified endpoints.
rationalSubtangleParametersEM[][10]: the array of parameters for the rational subtangles used in the connectionsEnterExitInput[][] arrays, used to to determine the new configurations and rotations.
connectionsCompass[][2][5]: an array to store the modified connections labeling, using compass direction labels.

Outputs (in the form of a modified input array):
connectionsCompass[][]


This function calls:
	determineCompassCorners()
This function is called by
	buildAlgebraicSubtangleEMCode()
*/

void convertConnectionsToCompass(int numOfSubtanglesInput, int connectionsEnterExitInput[][2][5], int connectionsCompass[][2][5], int rationalSubtangleParametersEM[][10]){
	
	//Arrays to store the corner information of the component subtangles.
	//Note that these will bee by default use the internal choice of corner A as NW.
	//Entry [i][0] denotes the Gauss index of the ith subtangle in the sum; the four middle entries denotes the corner (A,B,C,D) matches NW/NE/SE/SW.
	//Entry [i][5] denotes the canonical configuration of the component rational tangle, with respect to the choice of NW endpoint.
	//Entry [i][6] denotes the number of clockwise rotations on labels were needed for the NW corner of the subtangle to be regarded in the preferred position.
	int subtangleCorners[numOfSubtanglesInput][7];
	
	//Initialize a storage array for modifying the compass labels of the corners.
	int cornerStorageArray[7];
	
	//Fill in the initial compass regarding the internal choice of A as NW in each case.
	//Note that these directions might be modified if there is joining.
	for(int i=0; i<numOfSubtanglesInput; i++){
		
		//Use the determineCompassCorners() function to find the compass labels of the ith rational subtangle.
		//Store these in the cornerStorageArray[], and then use these to update the subtangleCorners[0][] array,
		cornerStorageArray[i]=i;
		determineCompassCorners(rationalSubtangleParametersEM[i][5],rationalSubtangleParametersEM[i][6],1,rationalSubtangleParametersEM[i][8],cornerStorageArray);
		
		subtangleCorners[i][1]=cornerStorageArray[1];
		subtangleCorners[i][2]=cornerStorageArray[2];
		subtangleCorners[i][3]=cornerStorageArray[3];
		subtangleCorners[i][4]=cornerStorageArray[4];
		
		subtangleCorners[i][5]=cornerStorageArray[5];
		subtangleCorners[i][6]=cornerStorageArray[6];
		
		//DEBUG:
		//printf("\n Subtangle Corners:\n Rational Component %d: NW = %d , NE = %d , SE = %d , SW = %d \n", i, subtangleCorners[i][1], subtangleCorners[i][2], subtangleCorners[i][3], subtangleCorners[i][4]);
		
	}
	
	//DEBUG
	/*
	printf("\n subtangleCorners:\n");
	for(int i=0; i<numOfSubtanglesInput; i++){
		printf(" [%d]\t [", i+1);
		for(int j=1; j<5; j++){
			printf(" %d", subtangleCorners[i][j]);
		}
		printf(" ]\n");
	}
	*/
	
	//Initialize an array to track the connections in terms of the compass corners, and fill these in.
	int rationalSubtangleConnectionsCompass[numOfSubtangles][2][5];
	for(int i=0; i<numOfSubtanglesInput; i++){
		connectionsCompass[i][0][0]=subtangleCorners[i][5];
		connectionsCompass[i][1][0]=subtangleCorners[i][6];
		for(int j=1; j<5; j++){
			//Display connections based on compass directions now, so the index of the subtangle connected at the [ * , NW , NE , SE , SW ] corners.
			connectionsCompass[i][0][j] = connectionsEnterExitInput[i][0][subtangleCorners[i][j]];
			//Update the connected corner to reflect these new labels; I'm sorry that this is unreadable.
			if( connectionsEnterExitInput[i][1][subtangleCorners[i][j]] == 0 ){
				connectionsCompass[i][1][j]=0;
			} else {
				for(int k=1; k<5; k++){
					if( subtangleCorners[connectionsCompass[i][0][j]-1][k]==connectionsEnterExitInput[i][1][subtangleCorners[i][j]] ){
						connectionsCompass[i][1][j]=k;
					}
				}
			}
		}
	}
	
}




/*
(NC, 1/21/19)
This function will rotate the labeling of a single, specified component in the connections array with the compass labeling.
It will also update the labels of components connected to this rotated subtangle.


Inputs:
numOfSubtanglesInput: the number of subtangles, used to determine the size of the connections array.
connectionsCompass[][2][5]: the array of connections information, using the compass labels.
rotatedComponentIndex: the index of the subtangle component being rotated.
numOfRotations: the number of rotations applied to the component being rotated.

Outputs (in the form of modified input arrays):
connectionsCompass[][2][5]: any references made to the rotated component are now updated to account for the rotation.


This function calls:
	N/A
This function is called by:
	buildAlgebraicSubtangleEMCode()

*/

void rotateComponentCompassConnections(int numOfSubtanglesInput, int connectionsCompass[][2][5], int rotatedComponentIndex, int numOfRotations){
	
	//Initialze some variables to store information needed for the rotation.
	int shift;
	int rotationStorage[2][5];
	
	//Record the number of rotations of this component in connectionsCompass[rotatedComponentIndex][1][0],
	connectionsCompass[rotatedComponentIndex][1][0]=numOfRotations;
	
	for(int i=1; i<5; i++){
		rotationStorage[0][i]=connectionsCompass[rotatedComponentIndex][0][i];
		rotationStorage[1][i]=connectionsCompass[rotatedComponentIndex][1][i];	
	}
	for(int i=1; i<5; i++){
		//The general formula for the permutation is x -> ( x+3 - R mod 4 ) + 1, where x is the old index and R is comp2rotations.
		//This is ugly becuase of how the modular arithmetic works, and since the labeling was originally chosen to be 1-4 rather than 0-3.
		shift = (((i+3)-numOfRotations)%4)+1;
		//DEBUG
		//printf("\n BOOP: i = %d , shift = %d", i, shift);
		connectionsCompass[rotatedComponentIndex][0][shift] = rotationStorage[0][i];
		connectionsCompass[rotatedComponentIndex][1][shift] = rotationStorage[1][i];
	}
	//Next, go through the connections arrays and update any connections to roteted component to reflect the new corner labels.
	for(int i=0; i<numOfSubtanglesInput; i++){
		for(int j=1; j<5; j++){
			//Check for connections to the rotated tangle, excluding endpoint connections and self-connections of the rotated tangle.
			if( (connectionsCompass[i][0][j]==rotatedComponentIndex+1) && (connectionsCompass[i][1][j]!=0) && (i!=rotatedComponentIndex) ){
				//Determine the shift based on the old label of the corner being updated, using the same formula given above.
				shift = (((connectionsCompass[i][1][j]+3)-numOfRotations)%4)+1;
				connectionsCompass[i][1][j] = shift;
			}
		}
	}
	
	//CONSIDER LATER WHETHER I NEED THIS
	//Update the information ragarding the canonical configuration of the rotated component; the configuration is specified relative to the choice of NW.
	//There are 4ish cases, depending on how many canonical configurations there are.
	int newConfig;
	if( connectionsCompass[rotatedComponentIndex][0][0] <= 0 ){
		//Cases 1 and 2: any configuration (-1) or no configuration (0); the configuration does not change.
		newConfig=connectionsCompass[rotatedComponentIndex][0][0];
	} else if( (connectionsCompass[rotatedComponentIndex][0][0]>0) && (connectionsCompass[rotatedComponentIndex][0][0]<=4) ){
		//Case 3: one configuration (1,2,3,4); rotate the configuration mod 4 (accounting for starting at 1 rather than 0)
		//Formula: x -> ((x+3)-R mod 4) + 1.
		newConfig = (((connectionsCompass[rotatedComponentIndex][0][0]+3)-numOfRotations)%4)+1;
	} else if( connectionsCompass[rotatedComponentIndex][0][0] > 4 ){
		//Case 4: two adjacent configurations (5,6,7,8); rotate the configuration mod 4 (accounting for starting at 1 rather than 0), and shift up.
		//Formula: x -> ((x+3)-R mod 4) + 5.
		newConfig = (((connectionsCompass[rotatedComponentIndex][0][0]+3)-numOfRotations)%4)+5;
	} else {
		//This shouldn't happen, but if for some reason the previous configuration was not specified, we will say there is no canonical configuration.
		newConfig=0;
	}
	connectionsCompass[rotatedComponentIndex][0][0]=newConfig;
}




/*
(NC, 1/21/19)
This function determines the new canonical configuration of a subtangle from two joined components.
The canonical configuration can be visualized as a 4-tuple, [ NW , NE , SE , SW ], where each entry is 1 if the tangle can be viewed as canonical after rotating to place that corner in the NW, and 0 if not.
There are ten possible configurations (building on certain properties of rational tangles), each of which is assigned the following label.

[ 1 , 1 , 1 , 1 ] = -1
[ 0 , 0 , 0 , 0 ] = 0
[ 1 , 0 , 0 , 0 ] = 1
[ 0 , 1 , 0 , 0 ] = 2
[ 0 , 0 , 1 , 0 ] = 3
[ 0 , 0 , 0 , 1 ] = 4
[ 1 , 1 , 0 , 0 ] = 5
[ 0 , 1 , 1 , 0 ] = 6
[ 0 , 0 , 1 , 1 ] = 7
[ 1 , 0 , 0 , 1 ] = 8

For compactness elsewhere, the canonical configuration is usually track by its integer label, but this makes modifying the configuration really ugly.
This function converts back to the 4-tuple form of the configuration of the two components.
When joining two components, the new configuration is really the intersection of the old configurations.
Equivalently, the only corners which are canonical after joining are those which were canonical for both components before joining.
We use this fact to compute a new 4-tuple; it has a 1 in any position where both components had a 1, and a 0 otherwise.
We then convert this 4-tuple back to one of the 10 labels above, and assign this as the configuration of the joined component.

This conversion makes the function a little bit unwieldy, but hopefully the logic is clear.
However, the nice thing is, this is independent of the type of joining (horizontal or vertical) and independent of the order of the components, assuming things have been properly rotated ahead of time.


Function Inputs:
firstConfig: the original configuration of the first joind component.

Function Outputs:
This function returns an integer value (one of the ten choices shown above) corresponding to the new configuration of the joined subtangle.

This function calls:
	N/A
This function is called by:
	joinAlgebraicSubtangles()
*/

int determineJoinedConfiguration(int firstConfig, int secondConfig){
	
	//Initialize arrays for the 4-tuples corresponding to the two input configurations.
	int inputConfigs[2];
	int tupleConfigs[2][4];
	int newTupleConfig[4];
	
	inputConfigs[0]=firstConfig;
	inputConfigs[1]=secondConfig;
	
	//Construct the corresponding 4-tuple for each configuration, following the conventions introduced above.
	//It's a bunch of ugly cases, I'm sorry.
	for(int i=0; i<2; i++){
		if( inputConfigs[i] == -1 ){
			tupleConfigs[i][0]=1;
			tupleConfigs[i][1]=1;
			tupleConfigs[i][2]=1;
			tupleConfigs[i][3]=1;
		} else if( inputConfigs[i] == 0 ){
			tupleConfigs[i][0]=0;
			tupleConfigs[i][1]=0;
			tupleConfigs[i][2]=0;
			tupleConfigs[i][3]=0;
		} else if( inputConfigs[i] == 1 ){
			tupleConfigs[i][0]=1;
			tupleConfigs[i][1]=0;
			tupleConfigs[i][2]=0;
			tupleConfigs[i][3]=0;
		} else if( inputConfigs[i] == 2 ){
			tupleConfigs[i][0]=0;
			tupleConfigs[i][1]=1;
			tupleConfigs[i][2]=0;
			tupleConfigs[i][3]=0;
		} else if( inputConfigs[i] == 3 ){
			tupleConfigs[i][0]=0;
			tupleConfigs[i][1]=0;
			tupleConfigs[i][2]=1;
			tupleConfigs[i][3]=0;
		} else if( inputConfigs[i] == 4 ){
			tupleConfigs[i][0]=0;
			tupleConfigs[i][1]=0;
			tupleConfigs[i][2]=0;
			tupleConfigs[i][3]=1;
		} else if( inputConfigs[i] == 5 ){
			tupleConfigs[i][0]=1;
			tupleConfigs[i][1]=1;
			tupleConfigs[i][2]=0;
			tupleConfigs[i][3]=0;
		} else if( inputConfigs[i] == 6 ){
			tupleConfigs[i][0]=0;
			tupleConfigs[i][1]=1;
			tupleConfigs[i][2]=1;
			tupleConfigs[i][3]=0;
		} else if( inputConfigs[i] == 7 ){
			tupleConfigs[i][0]=0;
			tupleConfigs[i][1]=0;
			tupleConfigs[i][2]=1;
			tupleConfigs[i][3]=1;
		} else if( inputConfigs[i] == 8 ){
			tupleConfigs[i][0]=0;
			tupleConfigs[i][1]=1;
			tupleConfigs[i][2]=1;
			tupleConfigs[i][3]=0;
		} else {
			//There shouldn't be any other possibilities, but if one of the configurations is not accounted for above, just return 0 (no canonical configuration).
			return 0;
		}
	}
	
	//Use the two tuple configurations to determine the new configuration tuple after joining.
	//For each entry in the new tuple, if the two joined tuples both have an entry of 1 in this position, then so does the new tuple; otherwise, this entry is 0.
	for(int i=0; i<4; i++){
		if( (tupleConfigs[0][i] == 1) && (tupleConfigs[1][i] == 1) ){
			newTupleConfig[i]=1;
		} else {
			newTupleConfig[i]=0;
		}
	}
	
	//Finally, use the newTupleConfig[] to determine the corresponding integer label. The function returns this value.
	//Again, this is a bunch of ugly cases, but the logic is hopefully clear.
	if( (newTupleConfig[0]==1) && (newTupleConfig[1]==1) && (newTupleConfig[2]==1) && (newTupleConfig[3]==1) ){
		return -1;
	} else if( (newTupleConfig[0]==0) && (newTupleConfig[1]==0) && (newTupleConfig[2]==0) && (newTupleConfig[3]==0) ){
		return 0;
	} else if( (newTupleConfig[0]==1) && (newTupleConfig[1]==0) && (newTupleConfig[2]==0) && (newTupleConfig[3]==0) ){
		return 1;
	} else if( (newTupleConfig[0]==0) && (newTupleConfig[1]==1) && (newTupleConfig[2]==0) && (newTupleConfig[3]==0) ){
		return 2;
	} else if( (newTupleConfig[0]==0) && (newTupleConfig[1]==0) && (newTupleConfig[2]==1) && (newTupleConfig[3]==0) ){
		return 3;
	} else if( (newTupleConfig[0]==0) && (newTupleConfig[1]==0) && (newTupleConfig[2]==0) && (newTupleConfig[3]==1) ){
		return 4;
	} else if( (newTupleConfig[0]==1) && (newTupleConfig[1]==1) && (newTupleConfig[2]==0) && (newTupleConfig[3]==0) ){
		return 5;
	} else if( (newTupleConfig[0]==0) && (newTupleConfig[1]==1) && (newTupleConfig[2]==1) && (newTupleConfig[3]==0) ){
		return 6;
	} else if( (newTupleConfig[0]==0) && (newTupleConfig[1]==0) && (newTupleConfig[2]==1) && (newTupleConfig[3]==1) ){
		return 7;
	} else if( (newTupleConfig[0]==1) && (newTupleConfig[1]==0) && (newTupleConfig[2]==0) && (newTupleConfig[3]==1) ){
		return 8;
	} else {
		//Nothing else should be possible, but if so, return 0 (no configurations possible)
		return 0;
	}
	
	
}



/*
(NC, 1/21/19)
This function determines the parity of the subtangle obtained after joining the two components.
This depends on the parity of the components and the type of join (vertical or horizontal).
This function ALSO accounts for the rotations of the second component, if any, which is not updated when the rotated labels are updated.


Inputs:
joinType: the type of joining, horizontal (+1) or vertical (-1).
firstParityInput: the input parity of the first component.
secondParityInput: the input parity of the second component.
secondComponentRotations: the number of rotations of the second component.

Outputs:
This function returns an integer value (0, 1, or 2(infinitiy) corresponding to the parity of the joined subtangle.


This function calls:
	N/A
This function is called by:
	joinAlgebraicSubtangles()	
*/

int determineJoinedParity(int joinType, int firstParityInput, int secondParityInput, int secondComponentRotations){
	
	//Initialize the parities of the two components.
	int firstParity = firstParityInput;
	int secondParity;
	//There are a few cases for the parity, depending on the number of rotations.
	if( secondParityInput == 1 ){
		//If the second subtangle has parity 1, then it still has parity 1 regardless of the number of rotations.
		secondParity = 1;
	} else {
		//If the second subtangle has parity 0 or infinity and the number of rotations is even, the parity stays the same.
		//If the number of rotations is odd, then parity 0 becomes parity infinity and vice versa.
		if( (secondComponentRotations%2) == 0 ){
			secondParity = secondParityInput;
		} else {
			if( secondParityInput == 0 ){
				secondParity = 2;
			} else if( secondParityInput == 2 ){
				secondParity = 0;
			}
		}
	}
	
	//Initialize the joined parity.
	int parityJoin;
	
	//For each of the two types of joining, there are (3 choose 2)=9 cases for the joined parity depending on the parity of each of the joined components.
	//This is a bunch of cases, but fortunately listing them all is straightforward.
	if( joinType == 1 ){
		//Horizontal join:
		if( (firstParity==0) && (secondParity==0) ){
			//Case 1: 0+0=0
			parityJoin=0;
		} else if( (firstParity==0) && (secondParity==1) ){
			//Case 2: 0+1=1
			parityJoin=1;
		} else if( (firstParity==0) && (secondParity==2) ){
			//Case 3: 0+infty=infty
			parityJoin=2;
		} else if( (firstParity==1) && (secondParity==0) ){
			//Case 4: 1+0=1
			parityJoin=1;
		} else if( (firstParity==1) && (secondParity==1) ){
			//Case 5: 1+1=0
			parityJoin=0;
		} else if( (firstParity==1) && (secondParity==2) ){
			//Case 6: 1+infty=infty
			parityJoin=2;
		} else if( (firstParity==2) && (secondParity==0) ){
			//Case 7: infty+0=infty
			parityJoin=2;
		} else if( (firstParity==2) && (secondParity==1) ){
			//Case 8: infty+1=infty
			parityJoin=2;
		} else if( (firstParity==2) && (secondParity==2) ){
			//Case 9: infty+infty= NOT VALID
			parityJoin=100;
		}
	} else if( joinType == -1 ){
		//Vertical join:
		if( (firstParity==0) && (secondParity==0) ){
			//Case 1: 0*0= NOT VALID
			parityJoin=100;
		} else if( (firstParity==0) && (secondParity==1) ){
			//Case 2: 0*1=0
			parityJoin=0;
		} else if( (firstParity==0) && (secondParity==2) ){
			//Case 3: 0*infty=0
			parityJoin=0;
		} else if( (firstParity==1) && (secondParity==0) ){
			//Case 4: 1*0=0
			parityJoin=0;
		} else if( (firstParity==1) && (secondParity==1) ){
			//Case 5: 1*1=infty
			parityJoin=2;
		} else if( (firstParity==1) && (secondParity==2) ){
			//Case 6: 1*infty=1
			parityJoin=1;
		} else if( (firstParity==2) && (secondParity==0) ){
			//Case 7: infty*0=0
			parityJoin=0;
		} else if( (firstParity==2) && (secondParity==1) ){
			//Case 8: infty*1=1
			parityJoin=1;
		} else if( (firstParity==2) && (secondParity==2) ){
			//Case 9: infty*infty=infty
			parityJoin=2;
		}
	}
	
	//Return the whatever the parityJoin found above is.
	return parityJoin;
	
}




/*
(NC, 1/21/19)
This function joins two subtangles into a single algebraic component, and updates the corresponding connection information in the compass connections array.
The new subtangle uses the lower index of the two joined subtangles (this is always the first of the two components).
References to the second joined component are now replaced with this new index.


Inputs:
numOfSubtanglesInput: the number of subtangles in connections array (before joining); used to determine the size of the array.
joinType: the type of joining, sum (+1) or product (-1); the first component is always on the left/top, respectively.
firstJoinIndex: the index of the first of the two components being joined.
secondJoinIndex: the index of the second of the two components being joined.
compassConnections[][][]: the array storing the subtangle connection infromation, using the compass labeling.
stageNumber: the current total number of joins accounted for, including this one.
algebraicComponentsParameters[][][]: the parameters of the algebraic subtangles, as determined by stage.

Outputs:


This function calls:
	determineJoinedConfiguration()
	determineJoinedParity()
This function is called by:
	buildAlgebraicSubtangleEMCode()


*/

void joinAlgebraicSubtangles(int numOfSubtanglesInput, int joinType, int firstJoinIndex, int secondJoinIndex, int secondComponentRotations, int compassConnections[][2][5], int stageNumber, int algebraicComponentsParameters[MAXNN][MAXNN][10]){
	
	//First, join second component to the first in the connections arrays.
	if( joinType == 1 ){
		//Case 1: horizontal joining (sum).
		//Replace the NE and SE connection corners of the first component with the information in the NE and SE connection corners of the second component.
		compassConnections[firstJoinIndex][0][2]=compassConnections[secondJoinIndex][0][2];
		compassConnections[firstJoinIndex][0][3]=compassConnections[secondJoinIndex][0][3];
		compassConnections[firstJoinIndex][1][2]=compassConnections[secondJoinIndex][1][2];
		compassConnections[firstJoinIndex][1][3]=compassConnections[secondJoinIndex][1][3];
	} else if( joinType == -1 ){
		//Case 2: vertical joining (product).
		//Replace the SW and SE connection corners of the first component with the information in the SW and SE connection corners of the second component.
		compassConnections[firstJoinIndex][0][4]=compassConnections[secondJoinIndex][0][4];
		compassConnections[firstJoinIndex][0][3]=compassConnections[secondJoinIndex][0][3];
		compassConnections[firstJoinIndex][1][4]=compassConnections[secondJoinIndex][1][4];
		compassConnections[firstJoinIndex][1][3]=compassConnections[secondJoinIndex][1][3];
	}
	
	//Next, determine the configuration of this new subtangle based on the configurations of the two components.
	//This is determined case by case in the determineJoinedConfiguration() function; see the description there for the full explanation.
	int newConfig;
	newConfig = determineJoinedConfiguration(compassConnections[firstJoinIndex][0][0],compassConnections[secondJoinIndex][0][0]);
	compassConnections[firstJoinIndex][0][0]=newConfig;
	
	//Now go through the algebraicComponentsParameters array for this stage number and update the information about the joining.
	//The firstJoinIndex row now refers to the new subtangle obtained after joining; all other rows are initially the same as the preceding stage.
	//These will be collapsed afterward to remove the secondJoinIndex subtangle row.
	//WORK IN PROGRESS
	for(int i=0; i<numOfSubtanglesInput; i++){
		if( i == firstJoinIndex ){
			//0: Stage joined, in this case the current stageNumber.
			algebraicComponentsParameters[stageNumber][i][0]=stageNumber;
			//1: The index of this component at the preceding stage; in this case firstJoinIndex.
			algebraicComponentsParameters[stageNumber][i][1]=firstJoinIndex;
			//2: The number of rational components, in this case 1 more than the preceding number.
			algebraicComponentsParameters[stageNumber][i][2]=algebraicComponentsParameters[stageNumber-1][firstJoinIndex][2]+1;
			//3: The index of the first joined component (since this subtangle is algebraic).
			algebraicComponentsParameters[stageNumber][i][3]=firstJoinIndex;
			//4: The index of the second joined component (since this subtangle is algebraic).
			algebraicComponentsParameters[stageNumber][i][4]=secondJoinIndex;
			//5: The internal parity of the subtangle, which is determined by a handful of case considerations accounted for in the determineJoinedParity() function.
			algebraicComponentsParameters[stageNumber][i][5]=determineJoinedParity(joinType,algebraicComponentsParameters[stageNumber-1][firstJoinIndex][5],algebraicComponentsParameters[stageNumber-1][secondJoinIndex][5],secondComponentRotations);
			//There is a bug here, sort it out!
			printf("\n BOOP: firstParity = %d , secondParity = %d , newParity = %d\n", algebraicComponentsParameters[stageNumber-1][firstJoinIndex][5], algebraicComponentsParameters[stageNumber-1][secondJoinIndex][5], algebraicComponentsParameters[stageNumber][i][5]);
			//6: The direction of the second strand in the subtangle, which depends on the joined components, the rotations of the second component, and on the join type.
			//DEAL WITH THIS LATER
			algebraicComponentsParameters[stageNumber][i][6]=100;
			//7: The canonical configuration of the current componet (10 possibilities, relative to the current choice of A=NW); newConfig as determined above.
			algebraicComponentsParameters[stageNumber][i][7]=newConfig;
			//8: The join type, horizontal (+1) or vertical (-1); this is stored in joinType in this function.
			algebraicComponentsParameters[stageNumber][i][8]=joinType;
			//9: The number of rotations of the corners of the second component in the joining; this is stored in secondComponentRotations in this function.
			algebraicComponentsParameters[stageNumber][i][0]=secondComponentRotations;
		} else {
			for(int j=0; j<10; j++){
				algebraicComponentsParameters[stageNumber][i][j]=algebraicComponentsParameters[stageNumber-1][i][j];
			}
		}
	}
	//DEBUG
	printf("\n DEBUG FLAG: Line %d, initial algebraicComponentsParameters at stage %d:", __LINE__, stageNumber);
	printAlgebraicComponentsParameters(numOfSubtanglesInput+stageNumber-1,stageNumber,algebraicComponentsParameters,true);
	
	//Next, go through all connections array entries and replace any references to the second component with the first component.
	for(int i=0; i<numOfSubtangles; i++){
		for(int j=1; j<5; j++){
			if( compassConnections[i][0][j] == (secondJoinIndex+1) ){
				compassConnections[i][0][j] = (firstJoinIndex+1);
			}
		}
	}
	
	//Next, collapse the connections array to eliminate the row for the second component.
	for(int i=(secondJoinIndex+1); i<numOfSubtangles; i++){
		for(int j=0; j<5; j++){
			compassConnections[i-1][0][j] = compassConnections[i][0][j];
			compassConnections[i-1][1][j] = compassConnections[i][1][j];
		}
	}	
	
	//Next, decrease the array index of any subtangles with index higher than the removed second component.
	//Effectively, this means going through and decreasing any reference to these updated indices by one.
	for(int i=(secondJoinIndex+1); i<=numOfSubtangles; i++){
		for(int j=0; j<numOfSubtangles-1; j++){
			for(int k=1; k<5; k++){
				if( compassConnections[j][0][k] == i ){
					compassConnections[j][0][k]--;
				}
			}
		}
	}
	
	

	
}





/*
(NC, 1/17/19)
This function is ued to build up the a version of the planar diagram code in terms of algebraic subtangles.


Inputs:

Outputs:


This function calls:
	printAlgebraicComponentsParameters()
	printArrayEM()
	printCompassRationalSubtangleEM()
	determineCompassCorners()
	rotateComponentCompassConnections()
This function is called by:
	buildGeneralizedEMCode()

*/
void buildAlgebraicSubtangleEMCode(int numOfSubtangles, int rationalSubtangleConnectionsEM[][2][5], int rationalSubtangleParametersEM[][10], int rationalTwistVectorsEM[][2*MAXNN+1], int *gaussRationalSubtangleEM, int *barsGaussRationalSubtangleEM, int endpointConnectionsCornersRational[2][4], int rationalComponentsCanonicalConfiguration[][4]){
	
	//Initialize an array to track the parameters of the algebraic components as they are joined in stages, and a varible to track the number of parameters.
	int algebraicComponentsParameters[MAXNN][MAXNN][10];
	int algebraicSubtanglesRefined=numOfSubtangles;
	int numOfStages=0;
	
	//Fill in the initial parameter information, regarding each rational subtangle as an algebraic component.
	for(int i=0; i<numOfSubtangles; i++){
		//0: Stage joined, by default 0 for unjoined rational subtangles.
		algebraicComponentsParameters[0][i][0]=0;
		//1: The index of this component at the preceding stage; in this case i.
		algebraicComponentsParameters[0][i][1]=i;
		//2: The number of rational components, by default 1.
		algebraicComponentsParameters[0][i][2]=1;
		//3: The numerator p from p/q (since this component is rational).
		algebraicComponentsParameters[0][i][3]=rationalSubtangleParametersEM[i][3];
		//4: The denominator q from p/q (since this component is rational).
		algebraicComponentsParameters[0][i][4]=rationalSubtangleParametersEM[i][4];
		//5: The internal parity of the subtangle, relative to corner A=NW.
		algebraicComponentsParameters[0][i][5]=rationalSubtangleParametersEM[i][5];
		//6: The direction of the second strand in the subtangle.
		algebraicComponentsParameters[0][i][6]=rationalSubtangleParametersEM[i][6];
		//7: The canonical configuration of the current componet (10 possibilities, relative to the current choice of A=NW).
		algebraicComponentsParameters[0][i][7]=rationalSubtangleParametersEM[i][8];
		//8: If a subtangle constructed from joined components, information denoting a tangle sum or product. In this case, 0 by defualt (no joins).
		algebraicComponentsParameters[0][i][8]=0;
		//9: If a subtangle constructed from joined components, information denoting the number of rotations of the corners of the second component. In this case, 0 by defualt (no joins).
		algebraicComponentsParameters[0][i][0]=0;
	}
	
	//DEBUG
	printf("\n DEBUG FLAG: Line %d, initial algebraicComponentsParameters at stage %d:", __LINE__, numOfStages);
	printAlgebraicComponentsParameters(numOfSubtangles,numOfStages,algebraicComponentsParameters,true);
	
	//Initialize a boolean in the special case where a local knot is detected, and another boolean to track when an indentified tangle has a non-canonical form.
	//If something strange happens where this tangle fails to be identified as Montesinos, flag this will a boolean as well.
	bool localKnotDetected = false;
	bool nonCanonical = false;
	bool somethingStrangeHappened = false;
	
	//Initialize some local arrays to store the rationalSubtangleConnectionsEM data.
	int rationalSubtangleConnectionsEMlocal[numOfSubtangles][2][5];
	for(int i=0; i<numOfSubtangles; i++){
		for(int j=0; j<5; j++){
			rationalSubtangleConnectionsEMlocal[i][0][j] = rationalSubtangleConnectionsEM[i][0][j];
			rationalSubtangleConnectionsEMlocal[i][1][j] = rationalSubtangleConnectionsEM[i][1][j];
		}
	}
	//Replace each endpoint corner with a 0; this way, it will be possible to avoid chaining across corners.
	rationalSubtangleConnectionsEMlocal[ endpointConnectionsCornersRational[0][0] - 1 ][1][ endpointConnectionsCornersRational[1][0] ] = 0;
	rationalSubtangleConnectionsEMlocal[ endpointConnectionsCornersRational[0][1] - 1 ][1][ endpointConnectionsCornersRational[1][1] ] = 0;
	rationalSubtangleConnectionsEMlocal[ endpointConnectionsCornersRational[0][2] - 1 ][1][ endpointConnectionsCornersRational[1][2] ] = 0;
	rationalSubtangleConnectionsEMlocal[ endpointConnectionsCornersRational[0][3] - 1 ][1][ endpointConnectionsCornersRational[1][3] ] = 0;
		
	//DEBUG
	//printArrayEM(numOfSubtangles,rationalSubtangleConnectionsEMlocal);
	
	//Initialize an array to track the connections in terms of the compass corners, and fill these in.
	int rationalSubtangleConnectionsCompass[numOfSubtangles][2][5];
	convertConnectionsToCompass(numOfSubtangles,rationalSubtangleConnectionsEMlocal,rationalSubtangleConnectionsCompass,rationalSubtangleParametersEM);
	
	//DEBUG
	printf("\n COMPASS CONNECTIONS: ");
	printArrayEM(numOfSubtangles,rationalSubtangleConnectionsCompass);
	printCompassRationalSubtangleEM(numOfSubtangles,rationalSubtangleConnectionsCompass,rationalSubtangleParametersEM,rationalTwistVectorsEM);
	
	
	//Initialize a boolean to track possible joines.
	int moreJoinsPossible = true;
	
	//Set up an array to store information about possible joining; specifically, track the type of join, the index of the joined subtangles, and new NW corner of the second.
	//The array will have the form [ join type , comp1 index , comp2 index , comp2 newNW ], where join type is =/-1 to denote horizontal/vertical joining relative to the first subtangle.
	int joinInfo[4];
	
	//While joining is possible, iterate through subtangle until a join is found.
	//If a join is found. update the arrays to reflect this and increment the stage.
	//If a join is not found, end the search.
	//The search continues until no joins are found, or this is only one component remaining (in which case the tangle is algebraic).
	while( moreJoinsPossible ){
		
		//Note that we only search for forward joining, so we do not check the very last subtangle for joins.
		//Also, multiple joins coulb be possible at the same stage, so break the loop after finding a join; these would be accounted for in later stages.
		for(int i=0; i<algebraicSubtanglesRefined-1; i++){
			//START AT i=0; that was just for testing something else
			
			if( (rationalSubtangleConnectionsCompass[i][0][2]==rationalSubtangleConnectionsCompass[i][0][3])
				&& (rationalSubtangleConnectionsCompass[i][1][2]!=0) && (rationalSubtangleConnectionsCompass[i][1][3]!=0)
				&& ( (rationalSubtangleConnectionsCompass[i][0][4]!=rationalSubtangleConnectionsCompass[i][0][3]) || (rationalSubtangleConnectionsCompass[i][1][4]==0) ) ){
				//If ( NE == SE ), and they are not joined across endpoints, and ( ( SW != SE ), except if it is also joind across an endpoint ), then horizontal joining.
				
				joinInfo[0] = 1;
				joinInfo[1] = i;
				joinInfo[2] = rationalSubtangleConnectionsCompass[i][0][2]-1;
				joinInfo[3] = rationalSubtangleConnectionsCompass[i][1][2];
				
				moreJoinsPossible = true;
				break;
				
			} else if( (rationalSubtangleConnectionsCompass[i][0][4]==rationalSubtangleConnectionsCompass[i][0][3])
						&& (rationalSubtangleConnectionsCompass[i][1][4]!=0) && (rationalSubtangleConnectionsCompass[i][1][3]!=0)
						&& ( (rationalSubtangleConnectionsCompass[i][0][2]!=rationalSubtangleConnectionsCompass[i][0][3]) || (rationalSubtangleConnectionsCompass[i][1][2]==0) ) ){
						//If ( SW == SE ), and they are not joined across endpoints, and ( ( NE != SE ), except if it is also joind across an endpoint ), then vertical joining.
				
				joinInfo[0] = -1;
				joinInfo[1] = i;
				joinInfo[2] = rationalSubtangleConnectionsCompass[i][0][4]-1;
				joinInfo[3] = rationalSubtangleConnectionsCompass[i][1][4];
				
				moreJoinsPossible = true;
				break;
				
			} else if( (rationalSubtangleConnectionsCompass[i][0][2]==rationalSubtangleConnectionsCompass[i][0][3]) && (rationalSubtangleConnectionsCompass[i][0][4]==rationalSubtangleConnectionsCompass[i][0][3])
						&& (rationalSubtangleConnectionsCompass[i][1][2]!=0) && (rationalSubtangleConnectionsCompass[i][1][3]!=0) && (rationalSubtangleConnectionsCompass[i][1][4]!=0) ){
						//If ( NE == SE ) AND ( SW == SE ), and none of these is joined across an endpoint, then this is a special case of a local knot; we do not consider this a possible join.
				localKnotDetected=true;
				moreJoinsPossible = false;
				joinInfo[0]=0;
				
			} else {
				//If the above conditions are not satisfied, then joining is not possible.
				moreJoinsPossible = false;
				joinInfo[0]=0;
			}
			
		}
		
		//If all components have been joined, and there is only one rational component left, then there is no need to check for further joins.
		if( algebraicSubtanglesRefined == 1 ){
			moreJoinsPossible = false;
			joinInfo[0]=0;
		}
		
		//If joining is detected, update the rationalSubtangleConnectionsCompass[][][] and the algebraicComponentsParameters[][][] arrays to reflect this joining, and move to the next stage.
		if( joinInfo[0] != 0 ){
			//First, we possibly need to rotate the second componnet before we may properly join it.
			//If rotating, we must also update any connections to this component in the rationalSubtangleConnectionsCompass[][][] array.
			//Do this via the rotateCompassConnections() function.
			rotateComponentCompassConnections(numOfSubtangles,rationalSubtangleConnectionsCompass,joinInfo[2],joinInfo[3]-1);

			//DEBUG:
			//printf("\n DEBUG FLAG: Line %d (after rotations)\n Join Number %d, join type %d detected between components %d and %d:", __LINE__, KillLoopTest+1, joinInfo[0], joinInfo[1]+1, joinInfo[2]+1);
			//printArrayEM(algebraicSubtanglesRefined,rationalSubtangleConnectionsCompass);
			
			
			//Now we may proceed with connections.
			//The join type is accounted for in the joining function.
			numOfStages++;
			joinAlgebraicSubtangles(algebraicSubtanglesRefined,joinInfo[0],joinInfo[1],joinInfo[2],joinInfo[3],rationalSubtangleConnectionsCompass,numOfStages,algebraicComponentsParameters);	
			algebraicSubtanglesRefined--;
				
			//DEBUG
			printf("\n DEBUG FLAG: Line %d \n Stage %d, after type %d joining of components %d and %d:", __LINE__, numOfStages, joinInfo[0], joinInfo[1]+1, joinInfo[2]+1);
			printArrayEM(algebraicSubtanglesRefined,rationalSubtangleConnectionsCompass);
				
		} else {
			//No joining is possible.
			moreJoinsPossible = false;
		}
		
		
	}
	
}



/*
(NC, 9/9/18)
This is a small function to be called when constructing the generalized Euling Millet Code of the tangle.
All of the heavy lifting in the constructions is handled by other functions, but things are consolidated here.
Some DEBUG printf statements are also called here as needed.
This function is currently a work in progress (NC, 9/9/18)

First: call findGaussCrossingSigns() to obtain an array of the crossing sign orienations indexed by the Gauss code index of each crossing (needed to construct the EM code).
Second: call findEMConnections() to construct regular EM code of a tangle from the Gauss code.
Third: call buildIntegerSubtangleEMCode() to construct a generalized vairant of the EM code, using integer subtangles rather than crossings.
Fourth: rational subtangles?


Inputs:

Outputs:


This function calls:
	findGaussCrossingSigns()
	findEMConnections()
	buildIntegerSubtangleEMCode()
	printArrayEM() (DEBUG)
	printStandardEM() (DEBUG)
This function is called by:
	main()

*/

void buildGeneralizedEMCode(int *gaussInput, int *orientedSignGaussInput, int *barsInput, int numOfCrossingsInput, int *gaussCrossingSigns, int connectionsEM[][2][5], int &numOfSubtangles, int integerSubtangleConnectionsEM[][2][5], int integerSubtangleParametersEM[][7], int *gaussIntegerSubtangleEM, int *barsGaussIntegerSubtangleEM){
		
	//Construct the usual Ewing Millett code of the tangle.
	findGaussCrossingSigns(gaussInput,orientedSignGaussInput,numOfCrossingsInput,gaussCrossingSigns);
	findEMConnections(gaussInput,gaussCrossingSigns,connectionsEM,numOfCrossingsInput);
	
	//DEBUG:
	printf("\n Input Gauss Code: [");
	for(int i=0; i<2*numOfCrossingsInput; i++){
		printf(" %d ", gaussInput[i]);
	}
	printf("]");
	printf("\n barsInput: [ %d %d ]\n",barsInput[0],barsInput[1]);
	printf("\n STANDARD EWING MILLETT CODE:");
	printArrayEM(numOfCrossingsInput, connectionsEM);
	printStandardEM(numOfCrossingsInput, connectionsEM);
	
	
	//DEBUG
	printf("\n\n\n INTEGER SUBTANGLE GENERALIZED EM CODE:\n");
	
	//Construct the generalized version of the EM code for integer subtangles.
	buildIntegerSubtanlgeEMCode(numOfCrossingsInput,gaussInput,connectionsEM,barsInput,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	//DEBUG
	//printf("\n INTEGER SUBTANGLE GENERALIZED EM CODE:");
	printArrayEM(numOfSubtangles,integerSubtangleConnectionsEM);
	printPrettyIntegerSubtangleEM(numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM);
	printf("\n Integer GaussEM Code: [");
	for(int i=0; i<2*numOfSubtangles; i++){
		printf(" %d ",gaussIntegerSubtangleEM[i]);
	}
	printf("]");
	printf("\n barsGaussIntegerSubtangleEM: [ %d %d ] \n\n",barsGaussIntegerSubtangleEM[0],barsGaussIntegerSubtangleEM[1]);
	
	
	//UPDATE THE INPUTS OF THIS FUNCTION TO INCLUDE THE RATIONAL VARIABLES, EVEN THOUGH RIGHT NOW THIS IS AVOIDED THANKS TO MAKING THESE VARIABLES GLOBAL, THAT COULD CHANGE IN THE FUTURE
	//Construct the generalized version of the EM code for rational subtangles.
	buildRationalSubtangleEMCode(numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM,rationalSubtangleConnectionsEM,rationalSubtangleParametersEM,rationalTwistVectorsEM,gaussRationalSubtangleEM,barsGaussRationalSubtangleEM,endpointConnectionsCornersRational,rationalComponentsCanonicalConfiguration);
	
	//DEBUG:
	printf("\n\n\n RATIONAL SUBTANGLE GENERALIZED EM CODE:\n");
	printRationalSubtangleParametersEM(numOfSubtangles,rationalSubtangleParametersEM);
	printArrayEM(numOfSubtangles,rationalSubtangleConnectionsEM);
	printPrettyRationalSubtangleEM(numOfSubtangles,rationalSubtangleConnectionsEM,rationalSubtangleParametersEM,rationalTwistVectorsEM);
	printf("\n Rational GaussEM Code: [");
	for(int i=0; i<2*numOfSubtangles; i++){
		printf(" %d ",gaussRationalSubtangleEM[i]);
	}
	printf("]");
	printf("\n barsGaussRationalSubtangleEM: [ %d %d ] \n\n",barsGaussRationalSubtangleEM[0],barsGaussRationalSubtangleEM[1]);
	
	
	
	//DEBUG
	printf("\n\n\n ALGEBRAIC SUBTANGLE GENERALIZED EM CODE: \n");
	
	buildAlgebraicSubtangleEMCode(numOfSubtangles,rationalSubtangleConnectionsEM,rationalSubtangleParametersEM,rationalTwistVectorsEM,gaussRationalSubtangleEM,barsGaussRationalSubtangleEM,endpointConnectionsCornersRational,rationalComponentsCanonicalConfiguration);

	
}



/*
(NC, 1/8/19)
This function is used to determine whether an input tangle is Montesinos or not based on its Rational Planar Diagram code.

Inputs:

Outputs:


This function calls:
	determineCompassCorners()
This function is called by:
	main()

*/

void detectIfMontesinos(int numOfSubtangles, int rationalSubtangleConnectionsEM[][2][5], int rationalSubtangleParametersEM[][10], int rationalTwistVectorsEM[][2*MAXNN+1], int *gaussRationalSubtangleEM, int *barsGaussRationalSubtangleEM, int endpointConnectionsCornersRational[2][4], int rationalComponentsCanonicalConfiguration[][4]){
	
	//Initialize an integer value to determine if the input tangle is Montesinos or not.
	//Use 1 to indicate a horizontal Montesinos tangle, -1 to indicate a vertical Montesinos tangle, and 0 to denote non-Montesinos tangle
	//Initialize a boolean in the special case where a local knot is detected, and another boolean to track when an indentified Montesinos tangle has a non-canonical form.
	//If something strange happens where this tangle fails to be identified as Montesinos, flag this will a boolean as well.
	int possibleMontesinos;
	bool localKnotDetected = false;
	bool nonCanonicalMontesinos = false;
	bool somethingStrangeHappened = false;
	
	//Endpoint connections information.
	//DEBUG
	/*
	printf("\n Endpoint Connections Information: \n Subtangle: \t [");
	for(int i=0; i<4; i++){
		printf(" %d ", endpointConnectionsCornersRational[0][i]);
	}
	printf("] \n Corner: \t [");
	for(int i=0; i<4; i++){
		printf(" %d ", endpointConnectionsCornersRational[1][i]);
	}
	printf("]\n");
	*/
	
	//Initialize some local arrays to store the rationalSubtangleConnectionsEM data.
	int rationalSubtangleConnectionsEMlocal[numOfSubtangles][2][5];
	for(int i=0; i<numOfSubtangles; i++){
		for(int j=0; j<5; j++){
			rationalSubtangleConnectionsEMlocal[i][0][j] = rationalSubtangleConnectionsEM[i][0][j];
			rationalSubtangleConnectionsEMlocal[i][1][j] = rationalSubtangleConnectionsEM[i][1][j];
		}
	}
	//Replace each endpoint corner with a 0; this way, it will be possible to avoid chaining across corners.
	rationalSubtangleConnectionsEMlocal[ endpointConnectionsCornersRational[0][0] - 1 ][1][ endpointConnectionsCornersRational[1][0] ] = 0;
	rationalSubtangleConnectionsEMlocal[ endpointConnectionsCornersRational[0][1] - 1 ][1][ endpointConnectionsCornersRational[1][1] ] = 0;
	rationalSubtangleConnectionsEMlocal[ endpointConnectionsCornersRational[0][2] - 1 ][1][ endpointConnectionsCornersRational[1][2] ] = 0;
	rationalSubtangleConnectionsEMlocal[ endpointConnectionsCornersRational[0][3] - 1 ][1][ endpointConnectionsCornersRational[1][3] ] = 0;
	
	//DEBUG
	//printArrayEM(numOfSubtangles,rationalSubtangleConnectionsEMlocal);
	
	
	//Arrays to store the corner information of the rational components.
	//Note that these will be determined relative to the natural description of a Montesinos tangle, which may differ from the internal description of the rational subtangles.
	//Entry [i][0] denotes the Gauss index of the ith subtangle in the sum; the four middle entries denotes the corner (A,B,C,D) matches NW/NE/SE/SW.
	//Entry [i][5] denotes the canonical configuration of the component rational tangle, with respect to the choice of NW endpoint.
	//Entry [i][6] denotes the number of clockwise rotations on labels were needed for the NW corner of the subtangle to be regarded in the usual Montesinos position.
	int subtangleCorners[numOfSubtangles][7];
	
	//By default, the first component in the possible Montesinos chain always has NW = A and has consistent orientation.
	subtangleCorners[0][0]=0;
	
	//Initialize a storage array for modifying the compass labels of the corners.
	int cornerStorageArray[7];
	
	//DEBUG
	//printf("\n Component 1: Guass index %d \n", subtangleCorners[0][0]+1);
	
	//If there is only one subtangle, the tangle is rational by default, and we will not consider it Montesinos in this case.
	if( numOfSubtangles < 2 ){
		possibleMontesinos = 0;
		
		//We will still compute the compass directions to describe the single rational component.
		//Use the determineCompassCorners() function to find the compass labels of the first rational subtangle.
		//Store these in the cornerStorageArray[], and then use these to update the subtangleCorners[0][] array,
		cornerStorageArray[0]=subtangleCorners[0][0];
		determineCompassCorners(rationalSubtangleParametersEM[0][5],rationalSubtangleParametersEM[0][6],1,rationalSubtangleParametersEM[0][8],cornerStorageArray);
		
		subtangleCorners[0][1]=cornerStorageArray[1];
		subtangleCorners[0][2]=cornerStorageArray[2];
		subtangleCorners[0][3]=cornerStorageArray[3];
		subtangleCorners[0][4]=cornerStorageArray[4];
		
		subtangleCorners[0][5]=cornerStorageArray[5];
		subtangleCorners[0][6]=cornerStorageArray[6];
			
	} else {
		
		//Use the determineCompassCorners() function to find the compass labels of the first rational subtangle.
		//Store these in the cornerStorageArray[], and then use these to update the subtangleCorners[0][] array,
		cornerStorageArray[0]=subtangleCorners[0][0];
		determineCompassCorners(rationalSubtangleParametersEM[0][5],rationalSubtangleParametersEM[0][6],1,rationalSubtangleParametersEM[0][8],cornerStorageArray);
		
		subtangleCorners[0][1]=cornerStorageArray[1];
		subtangleCorners[0][2]=cornerStorageArray[2];
		subtangleCorners[0][3]=cornerStorageArray[3];
		subtangleCorners[0][4]=cornerStorageArray[4];
		
		subtangleCorners[0][5]=cornerStorageArray[5];
		subtangleCorners[0][6]=cornerStorageArray[6];
		
		//DEBUG:
		//printf("\n Subtangle Corners:\n Rational Component 1: NW = %d , NE = %d , SE = %d , SW = %d \n", subtangleCorners[0][1], subtangleCorners[0][2], subtangleCorners[0][3], subtangleCorners[0][4]);
		
		//The tangle can only be Montesinos if each rational component shares exactly two (non-endpoint) connections with the next component, excluding possible endpoint connections.
		//If three connections are shared, and none of these cross over endpoints, then there is a local knot, in which case this is not Montesinos either.
		if( (rationalSubtangleConnectionsEMlocal[0][0][subtangleCorners[0][3]]==rationalSubtangleConnectionsEMlocal[0][0][subtangleCorners[0][2]])
			&& (rationalSubtangleConnectionsEMlocal[0][0][subtangleCorners[0][3]]==rationalSubtangleConnectionsEMlocal[0][0][subtangleCorners[0][4]]) ){
			//if ( SE == NE ) && ( SE == SW ), and there is not joining across endpoints, then this tangle has a local knot, and it is not possible to be Montesinos.
			if( (rationalSubtangleConnectionsEMlocal[0][1][subtangleCorners[0][2]]!=0)
				&& (rationalSubtangleConnectionsEMlocal[0][1][subtangleCorners[0][3]]!=0)
				&& (rationalSubtangleConnectionsEMlocal[0][1][subtangleCorners[0][4]]!=0) ){
				localKnotDetected = true;
				possibleMontesinos = 0;
				//If there is not a local knot, then ( SE == NE ) && ( SE == SW ) can only happen if at least one of these three connections involves joining across endpoints.
				//However, this could still be a Montesinos tangle as long as the endpoint joining only happens with the connection not involved with the Montesinos chain.
			} else if( (rationalSubtangleConnectionsEMlocal[0][1][subtangleCorners[0][2]]!=0)
						&& (rationalSubtangleConnectionsEMlocal[0][1][subtangleCorners[0][3]]!=0)
						&& (rationalSubtangleConnectionsEMlocal[0][1][subtangleCorners[0][4]]==0) ){
				//If both NE and SE connections are not across an endpoint, but the SW connection is, then horiozntal Montesinos is possible.
				possibleMontesinos = 1;	
			} else if( (rationalSubtangleConnectionsEMlocal[0][1][subtangleCorners[0][2]]==0)
						&& (rationalSubtangleConnectionsEMlocal[0][1][subtangleCorners[0][3]]!=0)
						&& (rationalSubtangleConnectionsEMlocal[0][1][subtangleCorners[0][4]]!=0)){
				//If both the SE and SW connections are not across an endpoint, but the NE connection is, then vertical Montesinos is possible.
				possibleMontesinos = -1;
			} else {
				//There should not be any other possibilities, but if something else does happen, this tangle probably isn't Montesinos, but for a strange resaon. Flag this just in case.
				possibleMontesinos = 0;
				somethingStrangeHappened = true;
			}
			//DEBUG:
			//printf("\n Debug Flag: Line %d: possibleMontesinos = %d \n", __LINE__, possibleMontesinos);
		} else if( rationalSubtangleConnectionsEMlocal[0][0][subtangleCorners[0][3]]==rationalSubtangleConnectionsEMlocal[0][0][subtangleCorners[0][2]] ){
			//if SE == NE, then there is possible horizontal joining.
			possibleMontesinos = 1;;
		} else if( rationalSubtangleConnectionsEMlocal[0][0][subtangleCorners[0][3]]==rationalSubtangleConnectionsEMlocal[0][0][subtangleCorners[0][4]] ){
			//if SE == SW, then there is possible vertical joining.
			possibleMontesinos = -1;
		} else {
			//Otherwise, Montesinos is not possible.
			possibleMontesinos = 0;
		}
	}
	
	//Intialize a variable to track the Guass code index of the next joined subtangle.
	//Initialize another variable to track the possible kind of joining of the next subtangle.
	int nextJoinedSubtangle;
	int nextPossibleMontesinos;
	
	//Initialze a variable to track the the corner label for the NW position of the next joined subtangle.
	//This is defined if possible Montesinos joining is detected.
	int nextSubtangleCornerNW;
	
	for(int i=1; i<numOfSubtangles; i++){
		if( possibleMontesinos != 0 ){
			//The index of the next subtangle in the Montesinos chain is the subtangle connected to the SE corner of the preceding subtangle.
			//Note that a -1 is needed to match the Gauss code index starting at 0.
			nextJoinedSubtangle = rationalSubtangleConnectionsEMlocal[subtangleCorners[i-1][0]][0][subtangleCorners[i-1][3]]-1;
			subtangleCorners[i][0] = nextJoinedSubtangle;
			cornerStorageArray[0] = i;
			
			//DEBUG
			//printf("\n Component %d: Gauss index %d \n", i+1, nextJoinedSubtangle+1);
			
			if( possibleMontesinos == 1 ){
				//If possible horizontal Montesinos:
				
				//The NW corner of the next component in the Montesinos chain is the corner connected to the NE corner of the preceding component in the chain.
				nextSubtangleCornerNW = rationalSubtangleConnectionsEMlocal[subtangleCorners[i-1][0]][1][subtangleCorners[i-1][2]];
				
				//Use this information and the determineCompassCorners() function to find the rest of the corners.
				determineCompassCorners(rationalSubtangleParametersEM[subtangleCorners[i][0]][5],rationalSubtangleParametersEM[subtangleCorners[i][0]][6],nextSubtangleCornerNW,rationalSubtangleParametersEM[subtangleCorners[i][0]][8],cornerStorageArray);
		
				subtangleCorners[i][1]=cornerStorageArray[1];
				subtangleCorners[i][2]=cornerStorageArray[2];
				subtangleCorners[i][3]=cornerStorageArray[3];
				subtangleCorners[i][4]=cornerStorageArray[4];
				
				subtangleCorners[i][5]=cornerStorageArray[5];
				subtangleCorners[i][6]=cornerStorageArray[6];
				
			} else if( possibleMontesinos == -1 ){
				//If possible vertical Montesinos:
				
				//The NW corner of the next component in the Montesinos chain is the corner connected to the SW corner of the preceding component in the chain.
				nextSubtangleCornerNW = rationalSubtangleConnectionsEMlocal[subtangleCorners[i-1][0]][1][subtangleCorners[i-1][4]];
				
				//Use this information and the determineCompassCorners() function to find the rest of the corners.
				determineCompassCorners(rationalSubtangleParametersEM[subtangleCorners[i][0]][5],rationalSubtangleParametersEM[subtangleCorners[i][0]][6],nextSubtangleCornerNW,rationalSubtangleParametersEM[subtangleCorners[i][0]][8],cornerStorageArray);
		
				subtangleCorners[i][1]=cornerStorageArray[1];
				subtangleCorners[i][2]=cornerStorageArray[2];
				subtangleCorners[i][3]=cornerStorageArray[3];
				subtangleCorners[i][4]=cornerStorageArray[4];
				
				subtangleCorners[i][5]=cornerStorageArray[5];
				subtangleCorners[i][6]=cornerStorageArray[6];
				
			}
			//DEBUG:
			//printf("\n Subtangle Corners:\n Rational Component %d: NW = %d , NE = %d , SE = %d , SW = %d \n", i+1, subtangleCorners[i][1], subtangleCorners[i][2], subtangleCorners[i][3], subtangleCorners[i][4]);
			
			
			//Next, we check if there is possible chaining of rational subtangles into a Motnesinos tangle.
			//We do not check for this if all subtangles have been accounted for already.
			if( (i+1) < numOfSubtangles ){
				//The tangle can only be Montesinos if each rational component shares exactly two connections with next component, excluding possible endpoint connections.
				//If three connections are shared, there is a local knot, in which case this is not Montesinos either.
				if( (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][0][subtangleCorners[i][3]]==rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][0][subtangleCorners[i][2]]) 
					&& (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][0][subtangleCorners[i][3]]==rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][0][subtangleCorners[i][4]]) ){
					//if ( SE == NE ) && ( SE == SW ), and we are not joining across an endpoint, then this is a local knot, and Montesinos is not possible.
					if( (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][1][subtangleCorners[i][2]]!=0)
						&& (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][1][subtangleCorners[i][3]]!=0)
						&& (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][1][subtangleCorners[i][4]]!=0) ){
						//If any of the three corners join across an endpoint, this entry is 0; check against this.
						localKnotDetected = true;
						nextPossibleMontesinos = 0;
						//If there is not a local knot, then ( SE == NE ) && ( SE == SW ) can only happen if at least one of these three connections involves joining across endpoints.
						//However, this could still be a Montesinos tangle as long as the endpoint joining only happens with the connection not involved with the Montesinos chain.
					} else if( (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][1][subtangleCorners[i][2]]!=0)
								&& (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][1][subtangleCorners[i][3]]!=0)
								&& (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][1][subtangleCorners[i][4]]==0) ){
						//If both NE and SE connections are not across an endpoint, but the SW connection is, then horiozntal Montesinos is possible.
						nextPossibleMontesinos = 1;
					} else if( (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][1][subtangleCorners[i][2]]==0)
								&& (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][1][subtangleCorners[i][3]]!=0)
								&& (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][1][subtangleCorners[i][4]]!=0) ){
						//If both the SE and SW connections are not across an endpoint, but the NE connection is, then vertical Montesinos is possible.
						possibleMontesinos = -1;
					} else {
						//There should not be any other possibilities, but if something else does happen, this tangle probably isn't Montesinos, but for a strange resaon. Flag this just in case.
						nextPossibleMontesinos = 0;
						somethingStrangeHappened = true;
					}
				} else if( rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][0][subtangleCorners[i][3]]==rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][0][subtangleCorners[i][2]] ){
					//if SE == NE, and endpoint joining does not occur, then there is possible horizontal joining.
					//Check to make sure joining does not happen across an endpoint.
					if( (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][1][subtangleCorners[i][3]]!=0) && (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][1][subtangleCorners[i][2]]!=0) ){
						nextPossibleMontesinos = 1;
					} else {
						nextPossibleMontesinos = 0;
						//printf("\n UN-check");
					}
					//printf("\n CHECK \n");
				} else if( rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][0][subtangleCorners[i][3]]==rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][0][subtangleCorners[i][4]] ){
					//if SE == SW, and endpoint joining does not occur, then there is possible vertical joining.
					if( (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][1][subtangleCorners[i][3]]!=0) && (rationalSubtangleConnectionsEMlocal[subtangleCorners[i][0]][1][subtangleCorners[i][4]]!=0) ){
						nextPossibleMontesinos = -1;
					} else {
						nextPossibleMontesinos = 0;
					}
				} else {
					//Otherwise, Montesinos is not possible.
					nextPossibleMontesinos = 0;
				}
				
				//If the nextPossibleMontesinos is 0, then the current tangle is not Montesisnos and we stop searching.
				//If the nextPossibleMontesinos is nonzero but does not agree with  possibleMontesinos, then flag this as non-Montesinos, but this would be strange since both types should not be detected at the same time.
				//If the nextPossibleMontesinos agrees with possibleMontesinos (1 or -1 for horizontal/vertical), then we leave possibleMontesinos as is and keep checking for joins to the Montesinos chain.
				if( nextPossibleMontesinos == 0 ){
					possibleMontesinos = 0;
				} else if( (nextPossibleMontesinos!=0) && (nextPossibleMontesinos!=possibleMontesinos) ){
					possibleMontesinos = 0;
					somethingStrangeHappened = true;
				}
					
			}
			
		}
		
	}
	
	//At this point, it has been determined whether or not the input tangle is Montesinos.
	//Next, if it is Montesinos, we check whether the tangle is in canononical form or not based on the configurations of the components; this information is stored in subtangleCorners[i][5].
	//To be canonical, each component must have a configuration of either -1, 1, 5, or 8 when viewed with its orientation matching the Montesinos chain.
	//Recall that these are the labels for the configurations in which the current NW corner matches a rational tangle in canonical form.
	for(int i=0; i<numOfSubtangles; i++){
		if( (subtangleCorners[i][5] != -1) && (subtangleCorners[i][5] != 1) && (subtangleCorners[i][5] != 5) && (subtangleCorners[i][5] != 8) ){
			nonCanonicalMontesinos = true;
		}
	}
	
	
	//Ininitialize an array to track the fraction p/q for each component of the possible Montesinos tangle.
	//Initialize an array to track the twist sign of each component.
	//Both of these 
	//The convention for a canonical montesinos tangle is that each rational component in the sum has a fraction satisfying 0 < p/q < 1.
	//Also initialize a boolean to track if this condition is satisfied; the tangle cannot be in canonical form if this holds.
	int componentFraction[numOfSubtangles][2];
	bool noncanonFraction = false;
	
	//If the tangle is in fact Montesinos, determine the fractions p/q of each component.
	if( possibleMontesinos !=0 ){
		for(int i=0; i<numOfSubtangles; i++){
			//The corresponding fraction will depend on whether or not the internal orientation of the component matches the orientation of the Montesinos chain.
			if( (subtangleCorners[i][6]%2) == 0 ){
				//If an even number of 90 degree clockwise rotations where needed to identify the NW corner of a rational component in the Montesinos construction, use fraction p/q.
				componentFraction[i][0] = rationalSubtangleParametersEM[subtangleCorners[i][0]][3];
				componentFraction[i][1] = rationalSubtangleParametersEM[subtangleCorners[i][0]][4];
			} else {
				//If an odd number of 90 degree clockwise rotations where needed to identify the NW corner of a rational component in the Montesinos construction, use fraction -q/p.
				//Note that we need to check against signs to remain consistent with the convetion that only the denominator is negative.
				if( rationalSubtangleParametersEM[subtangleCorners[i][0]][4] > 0 ){
					componentFraction[i][0] = rationalSubtangleParametersEM[subtangleCorners[i][0]][4];
					componentFraction[i][1] = -1*rationalSubtangleParametersEM[subtangleCorners[i][0]][3];
				} else {
					componentFraction[i][0] = -1*rationalSubtangleParametersEM[subtangleCorners[i][0]][4];
					componentFraction[i][1] = rationalSubtangleParametersEM[subtangleCorners[i][0]][3];
				}
			}
			
			//If the fraction of any component does not satisfy 1 < p/q < 0, then this fraction is non-canonical, and we set the boolean noncanonFraction to true.
			//In this case, we check whether q < 0 or p > q.
			if( (componentFraction[i][1]<0) || (componentFraction[i][0]>componentFraction[i][1]) ){
				noncanonFraction = true;
				//DEBUG
				//printf("\n BOOP %d", i);
			}
		}
	}
	
	//Special case for rational tangles, which only happens if there is a single component.
	if( numOfSubtangles == 1 ){
		componentFraction[0][0] = rationalSubtangleParametersEM[0][3];
		componentFraction[0][1] = rationalSubtangleParametersEM[0][4];
	}
	
	//If it is Montesinos, we would also like to describe the construction.
	if( possibleMontesinos == 1 ){
		printf("\n\n The current tangle is horizontal Montesinos,");
		if( nonCanonicalMontesinos || noncanonFraction ){
			printf(" but it is NOT canonical.\n");
			if( nonCanonicalMontesinos == true ){
				printf(" At least one component is not in a canonical configuration in the sum.\n");
			}
			if( noncanonFraction == true ){
				printf(" At least one component has a fraction which does not satisfy 0 < p/q < 1.\n");
			}
		} else {
			printf(" and it's also in CANONICAL form!\n");
		}
		printf("\n construction: \t");
		for(int i=0; i<numOfSubtangles; i++){
			printf("(%d/%d)", componentFraction[i][0], componentFraction[i][1]);
			if( (i+1) < numOfSubtangles ){
				printf("+");
			} else {
				printf("\n");
			}
		}
	} else if( possibleMontesinos == -1 ){
		printf("\n\n Current tangle is vertical Montesinos,");
		if( nonCanonicalMontesinos || noncanonFraction ){
			printf(" but it is NOT canonical.\n");
			if( nonCanonicalMontesinos == true ){
				printf(" At least one component is not in a canonical configuration in the sum.\n");
			}
			if( noncanonFraction == true ){
				printf(" At least one component has a fraction which does not satisfy 0 < p/q < 1.\n");
			}
		} else {
			printf(" and it's also in CANONICAL form!\n");
		}
		printf("\n construction: \t");
		for(int i=0; i<numOfSubtangles; i++){
			printf("(%d/%d)", componentFraction[i][0], componentFraction[i][1]);
			if( (i+1) < numOfSubtangles ){
				printf("*");
			} else {
				printf("\n");
			}
		}
	} else {
		if( numOfSubtangles == 1 ){
			printf("\n\n The current tangle is rational,");
			if( nonCanonicalMontesinos == true ){
				printf("\n but it isn't in canonical form.\n");
				printf("\n Equivalent to:\t (%d/%d)\n", componentFraction[0][0], componentFraction[0][1]);
			} else {
				printf("\n and it is in canononical form.\n");
				printf("\n Fraction:\t (%d/%d)\n", componentFraction[0][0], componentFraction[0][1]);
			}
		} else {
			printf("\n\n Current tangle is NOT Montesinos.\n");
			if( localKnotDetected ){
				printf(" A local knot was detected in this tangle.\n");
			}
			if( somethingStrangeHappened ){
				printf(" Something strange happened, look into why this tangle wasn't Montesinos.\n");
			}
		}
	
	}
	
}







int main(){
	
	//List of test inputs.
	
	int numOfCrossings1 = 3;
	int gauss1[2*numOfCrossings1] = {1,-2,3,-1,2,-3};
	int bars1[2] = {3,6};
	int orientedSignGauss1[2*numOfCrossings1] = {-1,-1,-1,-1,-1,-1};
	int parity1=1;
	
	//buildGeneralizedEMCode(gauss1,orientedSignGauss1,bars1,numOfCrossings1,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	int numOfCrossings2 = 4;
	int gauss2[2*numOfCrossings2] = {-1,2,3,-4,-2,1,4,-3};
	int bars2[2]={2,8};
	int orientedSignGauss2[2*numOfCrossings2] = {1,1,1,1,1,1,1,1};
	
	//buildGeneralizedEMCode(gauss2,orientedSignGauss2,bars2,numOfCrossings2,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	//THE REALLY GOOD ONE FROM MY OLD NOTEBOOK
	int numOfCrossings3 = 8;
	int gauss3[2*numOfCrossings3] = {1,-2,3,-1,-4,-5,6,4,7,-6,5,-7,-8,-3,2,8};
	int bars3[2]={8,16};
	int orientedSignGauss3[2*numOfCrossings3] = {-1,1,1,-1,1,1,1,1,-1,1,1,-1,1,1,1,1};
	int parity3=2;
	
	//buildGeneralizedEMCode(gauss3,orientedSignGauss3,bars3,numOfCrossings3,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//algebraic-COMPARE WITH BELOW
	int numOfCrossings4 = 6;
	int gauss4[2*numOfCrossings4] = {1,-2,3,4,-5,6,-4,5,-6,-1,2,-3};
	int bars4[2]={6,12};
	int orientedSignGauss4[2*numOfCrossings4] = {-1,-1,-1,1,1,1,1,1,1,-1,-1,-1};
	int parity4=0;
	
	//buildGeneralizedEMCode(gauss4,orientedSignGauss4,bars4,numOfCrossings4,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//rational-COMPARE WITH ABOVE
	int numOfCrossings5 = 6;
	int gauss5[2*numOfCrossings5] = {1,-2,3,4,-5,6,-6,5,-4,-1,2,-3};
	int bars5[2]={6,12};
	int orientedSignGauss5[2*numOfCrossings5] = {-1,-1,-1,1,1,1,1,1,1,-1,-1,-1};
	int parity5=0;
	
	//buildGeneralizedEMCode(gauss5,orientedSignGauss5,bars5,numOfCrossings5,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	int numOfCrossings6 = 4;
	int gauss6[2*numOfCrossings6] = {1,-2,3,-4,2,-1,4,-3};
	int bars6[2]={2,8};
	int orientedSignGauss6[2*numOfCrossings6] = {-1,-1,1,1,-1,-1,1,1};
	int parity6=0;
	
	//buildGeneralizedEMCode(gauss6,orientedSignGauss6,bars6,numOfCrossings6,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//Weird example from 10/7/18 in notes--good illustration that it works?
	int numOfCrossings7 = 8;
	int gauss7[2*numOfCrossings7] = {-1,2,-3,4,5,-6,6,-5,-7,3,-4,7,-8,1,-2,8};
	int bars7[2]={6,16};
	int orientedSignGauss7[2*numOfCrossings7] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	
	//buildGeneralizedEMCode(gauss7,orientedSignGauss7,bars7,numOfCrossings7,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//5 crossings rational tangle (canonical)
	int numOfCrossings8 = 5;
	int gauss8[2*numOfCrossings8] = {1,-2,3,-4,4,-5,2,-1,5,-3};
	int bars8[2]={4,10};
	int orientedSignGauss8[2*numOfCrossings8] = {1,1,-1,-1,-1,-1,1,1,-1,-1};
	
	//buildGeneralizedEMCode(gauss8,orientedSignGauss8,bars8,numOfCrossings8,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//5 crossing rational tangle (non-canonical)
	int numOfCrossings9 = 5;
	int gauss9[2*numOfCrossings9] = {1,-2,3,-4,2,-5,4,-3,5,-1};
	int bars9[2]={6,10};
	int orientedSignGauss9[2*numOfCrossings9] = {1,1,-1,-1,1,1,-1,-1,1,1};
	
	//buildGeneralizedEMCode(gauss9,orientedSignGauss9,bars9,numOfCrossings9,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//4 crossing rational tangle (non-canonical)
	int numOfCrossings10 = 4;
	int gauss10[2*numOfCrossings10] = {1,-2,3,-4,2,-1,4,-3};
	int bars10[2]={5,8};
	int orientedSignGauss10[2*numOfCrossings10] = {-1,-1,1,1,-1,-1,1,1};
	
	//buildGeneralizedEMCode(gauss10,orientedSignGauss10,bars10,numOfCrossings10,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	
	
	//6 crossing Montesinos tangle for comp
	int numOfCrossings11 = 6;
	int gauss11[2*numOfCrossings11] = {-1,2,-3,4,-5,6,-4,5,-6,1,-2,3};
	int bars11[2]={6,12};
	int orientedSignGauss11[2*numOfCrossings11] = {1,1,1,1,1,1,1,1,1,1,1,1};
	
	//buildGeneralizedEMCode(gauss11,orientedSignGauss11,bars11,numOfCrossings11,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//6 crossing Rational tangle for comp
	int numOfCrossings12 = 6;
	int gauss12[2*numOfCrossings12] = {-1,2,-3,4,-5,6,-6,5,-4,1,-2,3};
	int bars12[2]={6,12};
	int orientedSignGauss12[2*numOfCrossings12] = {1,1,1,1,1,1,1,1,1,1,1,1};
	
	//buildGeneralizedEMCode(gauss12,orientedSignGauss12,bars12,numOfCrossings12,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//12 crossing generalized Montesinos tangle for comp
	int numOfCrossings13 = 12;
	int gauss13[2*numOfCrossings13] = {1,-2,3,-4,5,-6,7,-8,9,-9,8,-7,10,-11,12,-3,2,-1,6,-5,4,-12,11,-10};
	int bars13[2]={9,24};
	int orientedSignGauss13[2*numOfCrossings13] = {1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1};
	
	//buildGeneralizedEMCode(gauss13,orientedSignGauss13,bars13,numOfCrossings13,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//10 crossing Algebraic tangle for comp
	//KNOWN BUG-THIRD SUBTANLGE HAS SIGN ERROR? ----- FIXED, there was an issue with how rational_q was defined in the special case of an even number of vertial twists
	int numOfCrossings14 = 10;
	int gauss14[2*numOfCrossings14] = {-1,2,-3,4,-5,1,-2,3,-6,7,-8,9,-10,6,-7,8,-9,10,-4,5};
	int bars14[2]={10,20};
	int orientedSignGauss14[2*numOfCrossings14] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	
	//buildGeneralizedEMCode(gauss14,orientedSignGauss14,bars14,numOfCrossings14,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//11 crossing Non-Algebraic tangle for comp
	//KNOWN BUG-CHAINING ACROSS ENDPOINTS ------ FIXED, there were multiple issues with labeling corner crossings in the integerEM function
	int numOfCrossings15 = 11;
	int gauss15[2*numOfCrossings15] = {-1,2,-3,4,-5,6,-7,8,-2,1,-9,10,-11,7,-6,5,-4,3,-8,11,-10,9};
	int bars15[2]={15,22};
	int orientedSignGauss15[2*numOfCrossings15] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	
	//buildGeneralizedEMCode(gauss15,orientedSignGauss15,bars15,numOfCrossings15,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//9 crossing Montesinos tangle
	int numOfCrossings16 = 9;
	int gauss16[2*numOfCrossings16] = {-1,2,-3,4,-5,6,-7,5,-4,3,-6,7,-8,9,-2,1,-9,8};
	int bars16[2]={7,18};
	int orientedSignGauss16[2*numOfCrossings16] = {1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,1,1,-1,-1};
	
	//buildGeneralizedEMCode(gauss16,orientedSignGauss16,bars16,numOfCrossings16,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//8 crossing non-canonical Montesinos tangle
	int numOfCrossings17 = 8;
	int gauss17[2*numOfCrossings17] = {-1,2,-3,4,-2,1,-5,6,-7,8,-6,5,-8,7,-4,3};
	int bars17[2]={8,16};
	int orientedSignGauss17[2*numOfCrossings17] = {-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1};
	
	//buildGeneralizedEMCode(gauss17,orientedSignGauss17,bars17,numOfCrossings17,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//8 crossing canonical Montesinos tangle (compare with above--rational components are equivalent)
	int numOfCrossings18 = 8;
	int gauss18[2*numOfCrossings18] = {-1,2,-3,4,-5,6,-4,3,-6,5,-7,8,-2,1,-8,7};
	int bars18[2]={4,16};
	int orientedSignGauss18[2*numOfCrossings18] = {1,1,1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1};
	
	//buildGeneralizedEMCode(gauss18,orientedSignGauss18,bars18,numOfCrossings18,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//8 crossing Montesinos tangle with a rotated rational component (compare with above two)
	int numOfCrossings19 = 8;
	int gauss19[2*numOfCrossings19] = {-1,2,-3,4,-5,6,-2,1,-6,5,-7,8,-4,3,-8,7};
	int bars19[2]={10,16};
	int orientedSignGauss19[2*numOfCrossings19] = {1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1};
	
	//buildGeneralizedEMCode(gauss19,orientedSignGauss19,bars19,numOfCrossings19,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//10 crossing Montesinos tangle with a rotated rational component--non-canonical (compare with above three)
	int numOfCrossings20 = 10;
	int gauss20[2*numOfCrossings20] = {-1,2,-3,4,-5,6,-7,8,-2,1,-8,7,-9,10,-6,5,-10,9,-4,3};
	int bars20[2]={12,20};
	int orientedSignGauss20[2*numOfCrossings20] = {1,1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1};
	
	//buildGeneralizedEMCode(gauss20,orientedSignGauss20,bars20,numOfCrossings20,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//10 crossing generalized Montesinos tangle--remove exterior twists before classifying?
	
	int numOfCrossings21 = 10;
	int gauss21[2*numOfCrossings21] = {-1,2,-3,4,-5,6,-7,8,-6,5,-9,3,-4,9,-10,1,-2,10,-8,7};
	int bars21[2]={6,20};
	int orientedSignGauss21[2*numOfCrossings21] = {-1,-1,-1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	
	//buildGeneralizedEMCode(gauss21,orientedSignGauss21,bars21,numOfCrossings21,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
		
	
	//13 crossing horizontal Montesinos tangle -- experiment with the function that classifies this.
	int numOfCrossings22 = 13;
	int gauss22[2*numOfCrossings22] = {-1,2,-3,4,-5,6,-2,1,-6,5,-7,8,-9,10,-11,3,-4,12,-13,11,-10,13,-12,7,-8,9};
	int bars22[2]={10,26};
	int orientedSignGauss22[2*numOfCrossings22] = {1,1,1,1,-1,-1,1,1,-1,-1,1,1,1,1,1,1,1,-1,-1,1,1,-1,-1,1,1,1};
	
	buildGeneralizedEMCode(gauss22,orientedSignGauss22,bars22,numOfCrossings22,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//10 crossing canononical Montesinos tangle
	int numOfCrossings23 = 10;
	int gauss23[2*numOfCrossings23] = {-1,2,-3,4,-5,6,-7,8,-9,10,-6,5,-4,1,-2,3,-8,9,-10,7};
	int bars23[2]={4,20};
	int orientedSignGauss23[2*numOfCrossings23] = {-1,-1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,1};
	
	//buildGeneralizedEMCode(gauss23,orientedSignGauss23,bars23,numOfCrossings23,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	
	
	
	//detectIfMontesinos(numOfSubtangles,rationalSubtangleConnectionsEM,rationalSubtangleParametersEM,rationalTwistVectorsEM,gaussRationalSubtangleEM,barsGaussRationalSubtangleEM,endpointConnectionsCornersRational,rationalComponentsCanonicalConfiguration);
	
}












