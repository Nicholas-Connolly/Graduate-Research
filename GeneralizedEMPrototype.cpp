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



//Global variables:
int gaussCrossingSigns[30];

//Does declaring an array of strings work in the same way?
//string ewingMillett[30][5];
int connectionsEM[30][2][5];

int numOfSubtangles;
int integerSubtangleConnectionsEM[30][2][5];
int integerSubtangleParametersEM[30][7];
int gaussIntegerSubtangleEM[30];
int barsGaussIntegerSubtangleEM[2];

int rationalSubtangleConnectionsEM[30][2][5];
int rationalSubtangleParametersEM[30][8];
int rationalTwistVectorsEM[30][30+1];
int gaussRationalSubtangleEM[30];
int barsGaussRationalSubtangleEM[30];









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

[	0					, 	1								,	2					,	3			,	4			,	5				,	6				,	7						]
[ numOfCrossings		,	number of integer subgtangles	,	twist vector length	,	fraction p	,	fraction q	,	internal parity	,	S2 direction	,	subtangle "shape" sign	]

0: The number of crossings in the subtangle, unsigned
1: The number of integer sbutangles (counting simultaneous twisting?)
2: The number of nonzero entries in the rational twist vector
3: The numerator p of the corresponding rational fraction p/q, calculated depending on other parameters.
4: The denominator q of the corresponding rational fraction p/q, calculated depending on other parameters.
5: The internal parity (0, 1, or infinity(2)) of the given subtangle, following the convention that corner "A" denotes the NW. This IGNORES the orientation of the second strand.
6: Denotes the direction of the second strand in the subtangle, relative to the usual conventions for its parity: 1 if in the usual direction, or -1 if in the reverse direction.
7: The sign of the crossings in the corresponding rational tangle of this shape (which must have constant sign if irreducible), with the regular S2 direction.


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

void printRationalSubtangleParametersEM(int numOfSubtanglesInput, int rationalSubtangleParametersEM[][8]){
	
	printf("\n rationalSubtangleParametersEM: \n");
	printf(" [ i ]\t [ 0\t 1\t 2\t 3\t 4\t 5\t 6\t 7\t]\n");
	printf(" sub- \t#cross \t #int \t TwiVec\t frac \t frac \t sub- \t s2- \t subtang");
	printf("\n tangle\t \t sub-t \t length\t p \t q \t parity\t direc \t sign \n\n");
	for(int i=0; i<numOfSubtanglesInput; i++){
		printf(" [ %d ]", i+1);
		printf("\t[");
		for(int j=0; j<8; j++){
			printf(" %d \t", rationalSubtangleParametersEM[i][j]);
		}
		printf("]\n");
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

void printPrettyRationalSubtangleEM(int numOfSubtanglesInput, int rationalSubtangleConnectionsEMinput[][2][5], int rationalSubtangleParametersEM[][8], int twistVectorArray[][31]){
	
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
			} else {
				//If none of the above happens, there is a some bug somewhere. Flag this "z" so it can identified at a glance.
				printf("z ");
			}
		}
		
		//Print the internal parity of the subtangle and the whether or not the second subtangle strand is reversed relative to this parity.
		if( rationalSubtangleParametersEM[i][5] == 2 ){
			printf("\t (Parity: infinity");
		} else {
			printf("\t (Parity: %d", rationalSubtangleParametersEM[i][5]);
		}
		if( rationalSubtangleParametersEM[i][6] < 0 ){
			printf("')\t");
		} else {
			printf(")\t");
		}
		
		//Print the twist vector--this is technically redundant with the fraction, but is useful information at a glance for reconstructing the subtangle by hand.
		printf("Twist Vector: [");
		for(int j=1; j<=twistVectorArray[i][0]; j++){
			printf(" %d ", twistVectorArray[i][j]);
		}
		printf("]\n");
		
	}
	
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
	
	for(int i=0; i<numOfSubtanglesInput; i++){
		
		//Each iteration, increment the number of refined subtangles.
		subtanglesRefined++;
		
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
		if( (integerSubtangleConnectionsEM[i][0][2] == integerSubtangleConnectionsEM[i][0][3]) && (i+1<numOfSubtangles) && ((i+1) != abs(gaussInput[barsInput[0]-1]) ) ){
			chainType=1;
		} else if ( (integerSubtangleConnectionsEM[i][0][2] == integerSubtangleConnectionsEM[i][0][4]) && (i+1<numOfSubtangles) && ((i+1) != abs(gaussInput[barsInput[0]-1]) ) ){
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
		//This only matters if there is chaining, and only if chaining an endpoint crossing. The updated index is the first crossing in the chain (i+1 in this case).
		//Since both chain types involve twisting corner B, the chained crossing index is the index of the crossing connected to this corner.
		//For integer subtangle chaining, the corner connected to the endpoint does NOT change (unless I'm mistaken, but future me can deal with that headache).
		if( chainType != 0 ){
			for(int j=0; j<4; j++){
				if( (endpointConnectionsCorners[0][j]-1) == integerSubtangleConnectionsEM[i][0][2] ){
					endpointConnectionsCorners[0][j] = (i+1);
				}
			}
		}
		
		
		
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
				if( (integerSubtangleConnectionsEM[i][0][2] == integerSubtangleConnectionsEM[i][0][3]) && (i+1<numOfSubtangles) && ((i+1) != abs(gaussInput[barsInput[0]-1])) ){
					
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
									//Set the endpoint chaining boolean to true.
									furtherChainingPossible=false;
								}
							}
						}
					}
					
					//It is also possible that an endpoint crossing gets chained into an integer subtangle. In this case, we need to update the endpointConnectionsCorners array to reflect this.
					//This only matters if there is chaining, and only if chaining an endpoint crossing. The updated index is the first crossing in the chain (i+1 in this case).
					//Since both chain types involve twisting corner B, the chained crossing index is the index of the crossing connected to this corner.
					//For integer subtangle chaining, the corner connected to the endpoint does NOT change (unless I'm mistaken, but future me can deal with that headache).
					if( furtherChainingPossible==true ){
						for(int j=0; j<4; j++){
							if( (endpointConnectionsCorners[0][j]-1) == integerSubtangleConnectionsEM[i][0][2] ){
								endpointConnectionsCorners[0][j] = (i+1);
							}
						}
						chainType++;
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
				if( (integerSubtangleConnectionsEM[i][0][2] == integerSubtangleConnectionsEM[i][0][4]) && (i+1<numOfSubtangles) && ((i+1) != abs(gaussInput[barsInput[0]-1])) ){
					
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
					if( furtherChainingPossible==true ){
						for(int j=0; j<4; j++){
							if( (endpointConnectionsCorners[0][j]-1) == integerSubtangleConnectionsEM[i][0][2] ){
								endpointConnectionsCorners[0][j] = (i+1);
							}
						}
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
		//The nextOldIndex is one less than the first index of the next chain.
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
void joinIntegerToRational(int rightTwist, int botTwist, int leftTwist, int topTwist, int &cornerNW, int &cornerNE, int &cornerSE, int &cornerSW, int &twistCount, int rationalTwistVectorArray[][31], int *integerSubtangleJoins, int initialIndex, int rationalSubtangleConnectionsEM[][2][5], int rationalSubtangleParametersEM[][8], int integerSubtangleParametersEM[][7], int &rationalSubtanglesRefined, int numOfIntegerSubtangles, int endpointConnectionsCorners[2][4], int *gaussRationalSubtangleEM, int *barsGaussRationalSubtangleEM){
	
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
	//In the rationalSubtangleConnectionsEM array, previous references to the newly joined subtangle must be updated to reflected the initial index of the rational amalgam.
	//The corners of the joined subtangle might not be the same as those of the amalgram; so we need to update not only the index, but the corner information as well.
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

void buildRationalSubtangleEMCode(int &numOfSubtangles, int integerSubtangleConnectionsEM[][2][5], int integerSubtangleParametersEM[][7], int *gaussIntegerSubtangleEM, int *barsGaussIntegerSubtangleEM, int rationalSubtangleConnectionsEM[][2][5], int rationalSubtangleParametersEM[][8], int rationalTwistVectorsEM[][31], int *gaussRationalSubtangleEM, int *barsGaussRationalSubtangleEM){
	
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
	int rationalTwistVectorArray[numOfSubtangles][31];
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
		//7: overall subtangle sign, as determined by shape: the sign of the crossing type, flipped depending on if S2 is reveresed.
		rationalSubtangleParametersEM[i][7] = ( integerSubtangleParametersEM[i][6] * (integerSubtangleParametersEM[i][0]/abs(integerSubtangleParametersEM[i][0])) );
		//printf("\n\t DEBUG sign of rational subtangle shape: %d", rationalSubtangleParametersEM[i][7]);
	}
	
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
	
	for(int i=0; i<numOfSubtangles; i++){
		
		//This array interates through the original list of integer subtangles and checks to see if each is the intial tangle is some rational amalgam.
		//It is possible that the current integer tangle might already have been joined to some other rational amalgam in a previous step, in which case we don't need to consider it separately.
		//Check to see if this is true before looking for further joinings; if it has already been joined, nothing else happens and we move to the next iteration of the for loop.
		if ( integerSubtangleJoins[i] == 1 ){
			
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
				
				
				//Attempt to iteratively join possible integer subtangles to the amalgam until no more twists are detected.
				joinIntegerToRational(rightTwist,botTwist,leftTwist,topTwist,cornerNW,cornerNE,cornerSE,cornerSW,localTwistCount,rationalTwistVectorArray,integerSubtangleJoins,localInitialRationalIndex,rationalSubtangleConnectionsEM,rationalSubtangleParametersEM,integerSubtangleParametersEM,rationalSubtanglesRefined,numOfSubtangles,endpointConnectionsCorners,gaussRationalSubtangleEM,barsGaussRationalSubtangleEM);
				
			}
					
		}
		
	}
	
	//DEBUG:
	//TO BE RE-INDEXED AT SOME POINT TO ACCOUNT FOR JOINED SUBTANGLES?
	//printRationalSubtangleParametersEM(numOfSubtangles,rationalSubtangleParametersEM);
	//printArrayEM(numOfSubtangles,rationalSubtangleConnectionsEM);
	
	//After the joining is completed, the sign must be specified on the entries in the twist vectors.
	//We assume that each twist vector entry must have the same sign (else, the corresponding rational tangle could be reduced).
	//Hence, the sign of each entry will match the sign* of the crossings in the initial integer subtangle in the amalgam.
	//In this case, we mean the sign* of the shape of the tangle, assumeing S2 were not reversed, which is possible.
	//To account for this case, the sign of the entries will be determined by the sign of the product of the sign of the crossings in the initial integer tange and the S2 direction variable.
	//This information is already tracked in entry 7 of the rationalSubtangleParametersEM array.
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
	
	
	//Now go through and reindex the rationalSubtangleConnectionsEM array and rationalSubtangleParametersEM array to eliminate joined tangles and collapse the labeling.
	int joinedIntegerSubtangles=0;
	for(int i=0; i<numOfSubtangles; i++){
		if( integerSubtangleJoins[i] == 0 ){
			//In this case, the subtangle is fully joind elsewhere and we can safely eliminate it.
			joinedIntegerSubtangles++;
			//Lower the index of anything higher than eliminated index by 1.
			for(int j=0; j<numOfSubtangles; j++){
				for(int k=1; k<5; k++){
					if( rationalSubtangleConnectionsEM[j][0][k] >= (i+1) ){
						rationalSubtangleConnectionsEM[j][0][k]--;
					}
				}
			}
		} else {
			//If this tangle is not eliminated, move its entries in the array up by the number of joined subtangles.
			//This way, only the joined subtangle portions of the array are overwitten; everything else just shifts up.
			for(int j=0; j<5; j++){
				rationalSubtangleConnectionsEM[i-joinedIntegerSubtangles][0][j] = rationalSubtangleConnectionsEM[i][0][j];
				rationalSubtangleConnectionsEM[i-joinedIntegerSubtangles][1][j] = rationalSubtangleConnectionsEM[i][1][j];
			}
			//Do the same thing for the rationalSubtangleParametersEM array.
			for(int j=0; j<8; j++){
				rationalSubtangleParametersEM[i-joinedIntegerSubtangles][j] = rationalSubtangleParametersEM[i][j];
			}
			//Do it for the twist vectors too, while we're at it.
			for(int j=0; j<=rationalTwistVectorArray[i][0]; j++){
				rationalTwistVectorArray[i-joinedIntegerSubtangles][j] = rationalTwistVectorArray[i][j];
			}
		}
	}
	
	/*
	//DEBUG:
	printf("\n\n DEBUG: COLLAPSED rationalSubtangleConnectionsEM array after re-labeling:\n");
	printArrayEM(rationalSubtanglesRefined,rationalSubtangleConnectionsEM);
	printRationalSubtangleParametersEM(rationalSubtanglesRefined,rationalSubtangleParametersEM);
	//DEBUG
	printf("\n Twist Vectors: \n");
	for(int i=0; i<rationalSubtanglesRefined; i++){
		printf(" [ %d ]\t[", i+1);
		for(int j=1; j<=rationalTwistVectorArray[i][0]; j++){
			printf(" %d ", rationalTwistVectorArray[i][j]);
		}
		printf("]\n");
	}
	
	//DEBUG
	printPrettyRationalSubtangleEM(rationalSubtanglesRefined,rationalSubtangleConnectionsEM,rationalSubtangleParametersEM,rationalTwistVectorArray);
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
	
	//The subtangles in the rational Gauss EM code are currently indexed by the integer subtangle index of the corresponding initial subtangle.
	//We would like to reindex these to be increasing on subtangles.
	
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
	buildRationalSubtangleEMCode(numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM,rationalSubtangleConnectionsEM,rationalSubtangleParametersEM,rationalTwistVectorsEM,gaussRationalSubtangleEM,barsGaussRationalSubtangleEM);
	
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

	
}










int main(){
	
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
	
	buildGeneralizedEMCode(gauss9,orientedSignGauss9,bars9,numOfCrossings9,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
	//4 crossing rational tangle (non-canonical)
	int numOfCrossings10 = 4;
	int gauss10[2*numOfCrossings10] = {1,-2,3,-4,2,-1,4,-3};
	int bars10[2]={5,8};
	int orientedSignGauss10[2*numOfCrossings10] = {-1,-1,1,1,-1,-1,1,1};
	
	//buildGeneralizedEMCode(gauss10,orientedSignGauss10,bars10,numOfCrossings10,gaussCrossingSigns,connectionsEM,numOfSubtangles,integerSubtangleConnectionsEM,integerSubtangleParametersEM,gaussIntegerSubtangleEM,barsGaussIntegerSubtangleEM);
	
	
}












