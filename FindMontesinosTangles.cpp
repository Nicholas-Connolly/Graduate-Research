#include <stdio.h>
#include <math.h>


//The maximum crossing number of the tangles being considered.
//The dimensions of the globally defined Gauss code arrays are defined relative to this; the arrays need to be large enough to hold the Gauss code information for a tangle with this many crossings.
//Note that the number of possible compositions grows exponetially with this bound, this could potentially lead to memory leaks or something?
#define MAXN 11
//Define the composition list bound (CLB) to be 2^(MAXN-1). This needed to specify the correct array sizes.
//This must be set manually, along with MAXN, I can't find a way to define it implicitly from MAXN just using preprocessor commands.
#define CLB 1024

//By default, we will assume a maximum crossing number of MAXN for a Montesinos tangle; arrays are then 2*MAXN, large enough to hold the corresponding Gauss code.
int gaussSUM[2*MAXN];
int barsSUM[2];
int orientedSignGaussSUM[2*MAXN];
int numOfCrossingsSUM;
int paritySUM;
int gaussPROD[2*MAXN];
int barsPROD[2];
int orientedSignGaussPROD[2*MAXN];
int numOfCrossingsPROD;
int parityPROD;

//Global varaibles use to store a locally generated component rational tangle. Overwritten as needed.
int twistCrossings;
int gaussTwist[2*MAXN];
int barsTwist[2];
int orientedSignGaussTwist[2*MAXN];
int parityTwist;

int rationalPQ[2];

//Three dimensional array to store the following information:
//First Index: the integer being partitioned.
//Second Two Indices: a matrix with rows corresponding to a distinct partition of the first index (for a maximum integer n, this matrix needs to be at least 2^(n-1) by n).
//Note that the number of nonzero entries for each composition varies by row.
int compositionList[MAXN+1][CLB][MAXN+1];
//The first index refers to the integer being partitioned. The second index is the number of terms in the corresponding row of the matrix of compositions.
int rationalTwistVectLength[MAXN+1][CLB];
//An array to store the positive and negative versions of a given twist vector. This is used as local storage and overwritten as needed.
int sameSignTwistVector[2][MAXN+1];

//Arrays to store the Gauss code information for the locally generated component rational tangles.
int compCross[MAXN];
int compGauss[MAXN][2*MAXN];
int compGaussSign[MAXN][2*MAXN];
int compBars[MAXN][2];
int compParity[MAXN];
int compRationalPQ[MAXN][2];

//Global arrays to store the Gauss code information of a Montesinos tangle as it's built up. This is overwritten each iteration.
int montCross;
int montGauss[2*MAXN];
int montGaussSign[2*MAXN];
int montBars[2];
int montParity;

//Global arrays to store the Gauss code information of a generalized Montesinos tangle, built up iteratively and overwritten.
int crossGM;
int gaussGM[2*MAXN];
int orientedSignGaussGM[2*MAXN];
int barsGM[2];
int parityGM;

//Variable to track the number of Montesinos and generalized Montesinos tangles found.
//Note that regular Montesinos tangles are included when counting all of the generalized Montesinos tangles, so this variable counts everything.
int totalMontesinosTanglesFound=0;
int totalGMTanglesFound=0;












/*
(NC, 7/22/18)
A quick recursive function to compute the greatest common denominator of two integers via the Euclidean algorithm.
Thank you to Wikipedia for the pseudocode!

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

int GCD(int a, int b){
	if( b == 0 ){
		return a;
	} else {
		return GCD(b, (a%b) );
	}
}



/*
(NC, 7/22/18)
This function will recover the corresponding rational number p/q for a rational tangle by simplifying the extended fraction defined by the twist vector.
This quantity is only meaningful for rational tangles, so this is checked as an input. The value of p and q uniquely classify the rational tangle.
These values will be stored as the first and second entries, respectively, of a global array rationalPQ[].
Note that these values should coincide with the invariants of the coloring matrix: a = p, c = q, and b = -p-q.
If p/q is negative, we will adopt the convention of making p positive and q negative.

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

void calculateExtendedFractionPQ(bool isRational, int twistVectLength, int *twistVector, int *rationalPQ){
	
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
			d = GCD(currentNumerator,currentDenominator);
			
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
				d = GCD(currentNumerator,currentDenominator);
				
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
(NC, 7/26/18)
Small function to shift the magnitude of an index by the specified amount without needing to check for the sign of the index explicitly.
*/

int shiftIndexMagnitude(int oldIndex, int shiftAmount){
	
	int newIndex;
	
	if( oldIndex >= 0 ){
		newIndex = (oldIndex + shiftAmount);
	} else if( oldIndex < 0 ){
		newIndex = (oldIndex - shiftAmount); 
	}
	
	return newIndex;
}


/*
(NC, 7/29/18)
This function constructs the Gauss code of a product of two tangles using the Gauss code of the left and right tangles. There are nine possible cases depending on the parities of the two tangles being multiplied.
The Gauss code of both the input and output tangles is stored across three arrays: gauss[], bars[], and orientedSignGauss[]. The number of crossings and tangle parity must also be accounted for.
The function is only intended for use with 2-string tangles. The two tangles being multiplied are destinguished by TOP and BOT, and the result is referred to as PROD.

Function inputs:
numOfCrossingsTOP: the number of crossings of the TOP tangle.
numOfCrossingsBOT: the number of crossings of the BOT tangle.
gaussTOP[]: an array storing the (ordered) Gauss code of the TOP tangle. Each entry contains the index of the crossing encountered. Positive entries denote an over-crossing; negative entries denote an under-crossing.
gaussBOT[]: an array storing the (ordered) Gauss code the BOT tangle, with the same convention described above.
orientedSignGaussTOP[]: an array of the same length as gaussTOP[] containing information on the orientation of each crossing encountered (either 1 or -1). Each entry matches the corresponding entry in gaussTOP[].
orientedSignGaussBOT[]: an array of the same length as gaussBOT[] containing the information on the orientation of each crossing encountered, with the same convention as described above.
barsTOP[]: an array with 2 entries to denote where each of the two strings exits the TOP tangle. The first/second entry denotes the last crossing encountered along the first/second string (regardless of crossing index).
barsBOT[]: an array with 2 entries to denote where each of the two strings exits the BOT tangle, with the same convention described above.
parityTOP: the parity of the TOP tangle, either 0, 1, or 2 (infinity).
parityBOT: the parity of the BOT tangle, either 0, 1, or 2 (infinity).

Function outputs (in the form of modified global variables):
numOfCrossingsPROD: the number of crossings of the PROD tangle.
gaussPROD[]: an array storing the (ordered) Gauss code the PROD tangle, with the same convention described above.
orientedSignGaussPROD[]: an array of the same length as gaussPROD[] containing the information on the orientation of each crossing encountered, with the same convention as described above.
barsPROD[]: an array with 2 entries to denote where each of the two strings exits the PROD tangle, with the same convention described above.
parityPROD: the parity of the PROD tangle, either 0, 1, or 2 (infinity).

This function calls:
	nothing
This function is called by:
	buildTwistTangle()

*/

void tangleProduct(int numOfCrossingsTOP, int *gaussTOP, int *barsTOP, int *orientedSignGaussTOP, int parityTOP, int numOfCrossingsBOT, int *gaussBOT, int *barsBOT, int *orientedSignGaussBOT, int parityBOT, int &numOfCrossingsPROD, int *gaussPROD, int *barsPROD, int *orientedSignGaussPROD, int &parityPROD){
	
	numOfCrossingsPROD = numOfCrossingsTOP + numOfCrossingsBOT;
	
	//Local arrays to store and manipulate the input information.
	int gaussLocalPROD[numOfCrossingsPROD];
	int orientedSignGaussLocalPROD[numOfCrossingsPROD];
	int barsCombined[4];
	
	//Initialize arrays for each of the four strands.
	int firstStrandTOPlength = barsTOP[0];
	int secondStrandTOPlength = barsTOP[1] - barsTOP[0];
	int firstStrandBOTlength = barsBOT[0];
	int secondStrandBOTlength = barsBOT[1] - barsBOT[0];
	int strand1[firstStrandTOPlength];
	int strand2[secondStrandTOPlength];
	int strand3[firstStrandBOTlength];
	int strand4[secondStrandBOTlength];
	int orientedSignStrand1[firstStrandTOPlength];
	int orientedSignStrand2[secondStrandTOPlength];
	int orientedSignStrand3[firstStrandBOTlength];
	int orientedSignStrand4[secondStrandBOTlength];
	
	//Initilize variables to represent the highest crossing index encounterd for each of the four strands.
	//Set the default values to be 0, this allows for strands with no crossings to be considered.
	//These are checked when updating indices of later strands.
	int strand1MaxIndex=0;
	int strand2MaxIndex=0;
	int strand3MaxIndex=0;
	int strand4MaxIndex=0;
	
	//Strand 1: first strand TOP
	for(int i=0;i<firstStrandTOPlength;i++){
		strand1[i]=gaussTOP[i];
		orientedSignStrand1[i]=orientedSignGaussTOP[i];
		//Update the strand1MaxIndex if needed.
		if( abs(strand1[i]) > strand1MaxIndex ){
			strand1MaxIndex = abs(strand1[i]);
		}
	}
	//Strand 2: second strand TOP
	for(int i=0;i<secondStrandTOPlength;i++){
		strand2[i]=gaussTOP[i+firstStrandTOPlength];
		orientedSignStrand2[i]=orientedSignGaussTOP[i+firstStrandTOPlength];
		//Update the strand2MaxIndex if needed.
		if( abs(strand2[i]) > strand2MaxIndex ){
			strand2MaxIndex = abs(strand2[i]);
		}
	}
	//Strand 3: first strand BOT
	for(int i=0;i<firstStrandBOTlength;i++){
		strand3[i]=gaussBOT[i];
		orientedSignStrand3[i]=orientedSignGaussBOT[i];
		//Update the strand3MaxIndex if needed.
		if( abs(strand3[i]) > strand3MaxIndex ){
			strand3MaxIndex = abs(strand3[i]);
		}
	}
	//Strand 4: second strand BOT
	for(int i=0;i<secondStrandBOTlength;i++){
		strand4[i]=gaussBOT[i+firstStrandBOTlength];
		orientedSignStrand4[i]=orientedSignGaussBOT[i+firstStrandBOTlength];
		//Update the strand4MaxIndex if needed.
		if( abs(strand4[i]) > strand4MaxIndex ){
			strand4MaxIndex = abs(strand4[i]);
		}
	}
	
	//Re-define strand lengths for more compact and intuitive notation:
	int strand1Length = firstStrandTOPlength;
	int strand2Length = secondStrandTOPlength;
	int strand3Length = firstStrandBOTlength;
	int strand4Length = secondStrandBOTlength;
	
	
	//DEBUG:
	/*
	printf("\n Strand 1: [");
	for(int i=0; i<strand1Length; i++){
		printf(" %d ", strand1[i]);
	}
	printf("]");
	printf("\n Strand 2: [");
	for(int i=0; i<strand2Length; i++){
		printf(" %d ", strand2[i]);
	}
	printf("]");
	printf("\n Strand 3: [");
	for(int i=0; i<strand3Length; i++){
		printf(" %d ", strand3[i]);
	}
	printf("]");
	printf("\n Strand 4: [");
	for(int i=0; i<strand4Length; i++){
		printf(" %d ", strand4[i]);
	}
	printf("]\n\n");
	
	printf(" Strand 1 Max: %d \n", strand1MaxIndex);
	printf(" Strand 2 Max: %d \n", strand2MaxIndex);
	printf(" Strand 3 Max: %d \n", strand3MaxIndex);
	printf(" Strand 4 Max: %d \n", strand4MaxIndex);
	*/
	
	
	//Each of the nine cases below depends on the parity of the constituent tangles.
	//The resulting product of tangles is constructed by connecting the four strands defined above in certain order, depending on the case.
	//As the strings are connected, the crossing indices need to be tracked and updated.
	//The array barsCombined[] will be used to track the last crossing on each of the four strands, in the order in which they are combined.
	//Each case will be prefaced by the certain order of the four strands.
	//By convetion, Strand 1 (firstStrandTOP) is always the first of the four strands in the new tangle.
	
	barsCombined[0]=barsTOP[0];
	for(int i=0; i<strand1Length; i++){
		gaussPROD[i]=strand1[i];
		orientedSignGaussPROD[i]=orientedSignStrand1[i];
	}

	
	//CASE 1:
	if( (parityTOP == 0) && (parityBOT == 0) ){
		//If two parity 0 tangles are multiplied, the resulting structure is not a 2-string tangle since it has some sort of closed-loop thing in the middle.
		//I don't know whether or not this structure is of interest, but this possibility will not be considered allowed for the purpose of this function.
		printf("\n ERROR: The product of two parity 0 tangles is no longer a 2-string tangle.\n");
		return;
	
	//CASE 2:	
	} else if( (parityTOP == 0) && (parityBOT == 1) ){
		
		//Strand Order: S1-S3r-S2r-S4
		barsCombined[1]=barsTOP[0]+barsBOT[0];
		barsCombined[2]=barsTOP[1]+barsBOT[0];
		barsCombined[3]=barsTOP[1]+barsBOT[1];
		
		//Index Shifts:
		//All crossings on S3 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by strand3MaxIndex.
		//Previously encountered crossings on S4 by strand1MaxIndex; previously un-encountered crossings on S4 by strand2MaxIndex.
		
		//Orientation Reversals:
		//Strands S3 and S2 are reversed. The order the strands are encountered in is not changed.
		//Any TOP crossing on both S1 and S2 has its orientation reversed; likewise, any BOT crossings on both S3 and S4 has its orientation reversed.
		
		
		//REVERSAL CONSTRUCTION:
		//Initilize the reversals of S2 and S3.
		int strand2reverse[strand2Length];
		int strand3reverse[strand3Length];
		int orientedSignStrand2reverse[strand2Length];
		int orientedSignStrand3reverse[strand3Length];
		//List the crossings in the opposite order on each strand:
		for(int i=strand2Length; i>0; i--){
			strand2reverse[i-1]=strand2[strand2Length-i];
			orientedSignStrand2reverse[i-1]=orientedSignStrand2[strand2Length-i];
		}
		for(int i=strand3Length; i>0; i--){
			strand3reverse[i-1]=strand3[strand3Length-i];
			orientedSignStrand3reverse[i-1]=orientedSignStrand3[strand3Length-i];
		}
		
		
		//DEBUG:
		/*
		printf("\n\n BEFORE RE-LABEL");
		printf("\n strand2reverse: [");
		for(int i=0;i<strand2Length;i++){
			printf(" %d ", strand2reverse[i]);
		}
		printf("]");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]");
		printf("\n strand4: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4[i]);
		}
		printf("]\n");
		*/
		
		//Since the crossings on each of S2 and S3 are now encountered in the opposite order, the indices in the reverse arrays must be re-labeled accordingly.
		//Arrays to specfiy if an entry of the reverse arrays has been modified or not; used to avoid undoing the updated labeling.
		int modifyCheckS2r[strand2Length];
		int modifyCheckS3r[strand3Length];
		int modifyCheckS4[strand4Length];
		//Initialize both arrays as only entires of 0. Each time a label is updated, change entry to 1.
		for(int i=0; i<strand2Length; i++){
			modifyCheckS2r[i]=0;
		}
		for(int i=0; i<strand3Length; i++){
			modifyCheckS3r[i]=0;
		}
		for(int i=0; i<strand4Length; i++){
			modifyCheckS4[i]=0;
		}
		//The crossing indices must be modified so that they are increasing while traveling in reverse order.
		//Iterate through each entry in the strand and relabel the crossings as encountered.
		//Check if any other entries in the strand match this one; update those entries correspondingly.
		//Check against the modifyCheck array to avoid modifying an entry that has already been accounted for.
		
		//First relabel the crossing indices on S2r. Any crossing which is on both S1 and S2 has already been accounted for.
		//Thus, we need only consider those crossing indices which are greater than strand1MaxIndex.
		int S2crossCount=strand1MaxIndex;
		int storeCross;
		for(int i=0; i<strand2Length; i++){
			//We only need to consider those crossings which are BOTH not shared with S1 and which have not already been modified.
			if( (modifyCheckS2r[i] == 0) && ( abs(strand2reverse[i]) > strand1MaxIndex ) ){
				S2crossCount++;
				storeCross=strand2reverse[i];
				if(storeCross>0){
					strand2reverse[i]=S2crossCount;
				} else {
					strand2reverse[i]=(-1*S2crossCount);
				}
				modifyCheckS2r[i]=1;
				//After relabeling this entry, check for corresponding entries later along S2r.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand2Length; j++){
					if ( (abs(strand2reverse[j]) == abs(storeCross)) && (modifyCheckS2r[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand2reverse[j]=(-1*strand2reverse[i]);
						modifyCheckS2r[j]=1;
					}
				}
			}
		}
		//At this point, all crossings indices on S2r not shared with S1 are relabled to be increasing.
		
		//NOTE: the order to run these loops is important. In this case, we travel along S3r, then S4.
		//The same crossing index could show up either twice on one strand, or once on both strands; both possibilities must be accounted for.
		//Iterate first along S3r, whenever an entry is relabeled, check for corresponding entries both on S3r AND on S4 and update these as well.
		//Since S4 is not reversed, only the S4 crossings shared with S3r need to be updated. The remaining S4 crossings stay as is.
		int reverseCrossCount=0;
		for(int i=0; i<strand3Length; i++){
			//If modifyCheckS3r[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS3r[i] == 0){
				reverseCrossCount++;
				storeCross=strand3reverse[i];
				if(storeCross>0){
					strand3reverse[i]=reverseCrossCount;
				} else {
					strand3reverse[i]=(-1*reverseCrossCount);
				}
				modifyCheckS3r[i]=1;
				//After relabeling this entry, check for corresponding entries later along S3r.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand3Length; j++){
					if ( (abs(strand3reverse[j]) == abs(storeCross)) && (modifyCheckS3r[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand3reverse[j]=(-1*strand3reverse[i]);
						modifyCheckS3r[j]=1;
					}
				}
				//Next, check for later corresponding entries along S4 with the same conventions as above.
				for(int j=0; j<strand4Length; j++){
					if( (abs(strand4[j]) == abs(storeCross)) && (modifyCheckS4[j] == 0) ){
						strand4[j]=(-1*strand3reverse[i]);
						modifyCheckS4[j]=1;
					}
				}
			}
		}
		//The above for-loop accounts for all crossings along S3r, including those which also show up on S4. We do not need to modify any other crossings on S4.
		//At this point, all labels of the strands in the BOT tangle are updated.
		
		//DEBUG:
		/*
		printf("\n AFTER RE-LABEL");
		printf("\n strand2reverse: [");
		for(int i=0;i<strand2Length;i++){
			printf(" %d ", strand2reverse[i]);
		}
		printf("]");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]");
		printf("\n strand4: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4[i]);
		}
		printf("]\n");
		*/
		
		
		//Since S2 was reversed but S1 was not, any crossing which was on both S1 and S2 has its orientation reversed. This must be accounted for in both orientedSign arrays.
		//Is is sufficient to use a nested for-loop to search for any crossing indices occuring on both strands, and update the corresponding orientedSign entries.
		for(int i=0; i<strand1Length; i++){
			for(int j=0; j<strand2Length; j++){
				//Each crossing index in the TOP tangle occurs twice; if this occurence happens to be on the two different strands, then the orientation is reversed accordingly.
				if( abs(strand1[i]) == abs(strand2reverse[j]) ){
					orientedSignStrand1[i] = (-1*orientedSignStrand1[i]);
					orientedSignStrand2reverse[j] = (-1*orientedSignStrand2reverse[j]);
				}
			}
		}
		//Likewise, since S3 was reversed but S4 was not, we must again account for this in the orientedSign arrays in the same way as above.
		for(int i=0; i<strand3Length; i++){
			for(int j=0; j<strand4Length; j++){
				//Each crossing index in the BOT tangle occurs twice; if this occurence happens to be on the two different strands, then the orientation is reversed accordingly.
				if( abs(strand3reverse[i]) == abs(strand4[j]) ){
					orientedSignStrand3reverse[i] = (-1*orientedSignStrand3reverse[i]);
					orientedSignStrand4[j] = (-1*orientedSignStrand4[j]);
				}
			}
		}
		//If the orientation has been changed for any crossings on S1, we need to update the orientedSignGaussSUM[] array since this was earlier set equal to orientedSignStrand1.
		for(int i=0; i<strand1Length; i++){
			orientedSignGaussPROD[i]=orientedSignStrand1[i];
		}
		
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussPROD[i] = shiftIndexMagnitude(strand3reverse[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussPROD[i]=orientedSignStrand3reverse[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2reverse[i-barsCombined[1]]) <= strand1MaxIndex ){
				gaussPROD[i] = strand2reverse[i-barsCombined[1]];
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand2reverse[i-barsCombined[1]],strand3MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand2reverse[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note that crossings already encountered on S3r are shifted to match.
			if( abs(strand4[i-barsCombined[2]]) <= strand3MaxIndex ){
				gaussPROD[i] = shiftIndexMagnitude(strand4[i-barsCombined[2]],strand1MaxIndex);
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand4[i-barsCombined[2]],strand2MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand4[i-barsCombined[2]];
		}
		
		parityPROD=0;
		barsPROD[0]=barsCombined[0];
		barsPROD[1]=barsCombined[3];
			
	
	//CASE 3:	
	} else if( (parityTOP == 0) && (parityBOT == 2) ){
		
		//Strand Order: S1-S4-S2-S3
		barsCombined[1]=barsTOP[0]+(barsBOT[1]-barsBOT[0]);
		barsCombined[2]=barsTOP[1]+(barsBOT[1]-barsBOT[0]);
		barsCombined[3]=barsTOP[1]+barsBOT[1];
		
		//Index Shifts:
		//All crossings on S4 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by RELABELEDstrand4MaxIndex.
		//Previously encountered crossings on S3 by strand1MaxIndex; previously un-encountered crossings on S4 by strand2MaxIndex.
		
		//Orientation Reversals:
		//The strands S3 and S4 are now encountered in reverse order, even though the orientation along these strands is not changed.
		//Effectively, this means that crossing indices must be relabeled on both strands to reflect this, but crossing orientations do not change.
		
		
		//DEBUG:
		/*
		printf("\n\n BEFORE RE-LABEL");
		printf("\n strand3: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3[i]);
		}
		printf("]\n");
		printf(" strand4: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4[i]);
		}
		printf("]\n");
		*/
		
		//Since strands S3 and S4 are now travelled in the opposite order, the crossing labels must be updated to be increasing along S4 and then along S3.
		//Arrays to specfiy if an entry of S3 or S4 has been modified or not; used to avoid undoing the updated labeling.
		int modifyCheckS3[strand3Length];
		int modifyCheckS4[strand4Length];
		//Initialize the arrays with only entires of 0. Each time a label is updated, change entry to 1.
		for(int i=0; i<strand3Length; i++){
			modifyCheckS3[i]=0;
		}
		for(int i=0; i<strand4Length; i++){
			modifyCheckS4[i]=0;
		}
		//NOTE: the order to run these loops is important. In this case, we travel along S4, then S3.
		//The same crossing index could show up either twice on one strand, or once on both strands; both possibilities must be accounted for.
		//Iterate first along S4, whenever an entry is re-labeled, check for corresponding entries both on S4 AND on S3 and update these as well.
		//After this, iterate along S3 and update any remaining crossing labels that are not shared with S4.
		int crossCount=0;
		int storeCross;
		for(int i=0; i<strand4Length; i++){
			//If modifyCheckS4[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS4[i] == 0){
				crossCount++;
				storeCross=strand4[i];
				if(storeCross>0){
					strand4[i]=crossCount;
				} else {
					strand4[i]=(-1*crossCount);
				}
				modifyCheckS4[i]=1;
				//After relabeling this entry, check for corresponding entries later along S4.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand4Length; j++){
					if ( (abs(strand4[j]) == abs(storeCross)) && (modifyCheckS4[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand4[j]=(-1*strand4[i]);
						modifyCheckS4[j]=1;
					}
				}
				//Next, check for later corresponding entries along S3 with the same conventions as above.
				for(int j=0; j<strand3Length; j++){
					if( (abs(strand3[j]) == abs(storeCross)) && (modifyCheckS3[j] == 0) ){
						strand3[j]=(-1*strand4[i]);
						modifyCheckS3[j]=1;
					}
				}
			}
		}
		//The above for-loop accounts for all crossings along S4, including those which also show up on S3.
		//Now we move along S3 and update any remaining crossings, with the same conventions as above.
		for(int i=0; i<strand3Length; i++){
			//If modifyCheckS3r[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS3[i] == 0){
				crossCount++;
				storeCross=strand3[i];
				if(storeCross>0){
					strand3[i]=crossCount;
				} else {
					strand3[i]=(-1*crossCount);
				}
				modifyCheckS3[i]=1;
				for(int j=0; j<strand3Length; j++){
					if ( (abs(strand3[j]) == abs(storeCross)) && (modifyCheckS3[j] == 0) ){
						strand3[j]=(-1*strand3[i]);
						modifyCheckS3[j]=1;
					}
				}
			}
		}
		//At this point, all labels of the crossings along the strands in the BOT tangle are updated to travel first along S4 and then along S3.
		//In order to correctly shift the crossing index when building the combined tangle, we also need to record the maximum crossing index encountered on the relabeled S4.
		int RELABELEDstrand4MaxIndex=0;
		for(int i=0; i<strand4Length; i++){
			if( abs(strand4[i]) > RELABELEDstrand4MaxIndex ){
				RELABELEDstrand4MaxIndex = abs(strand4[i]);
			}
		}
		
		//DEBUG:
		/*
		printf("\n AFTER RE-LABEL");
		printf("\n strand3: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3[i]);
		}
		printf("]");
		printf("\n strand4: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4[i]);
		}
		printf("]\n");
		*/
		
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussPROD[i] = shiftIndexMagnitude(strand4[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussPROD[i]=orientedSignStrand4[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2[i-barsCombined[1]]) <= strand1MaxIndex ){
				gaussPROD[i] = strand2[i-barsCombined[1]];
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand2[i-barsCombined[1]],RELABELEDstrand4MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand2[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note crossings already encountered on strand4 are shifted to match.
			if( abs(strand3[i-barsCombined[2]]) <= RELABELEDstrand4MaxIndex ){
				gaussPROD[i] = shiftIndexMagnitude(strand3[i-barsCombined[2]],strand1MaxIndex);
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand3[i-barsCombined[2]],strand2MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand3[i-barsCombined[2]];
		}
		
		parityPROD=0;
		barsPROD[0]=barsCombined[0];
		barsPROD[1]=barsCombined[3];
		
		
	//CASE 4:	
	} else if( (parityTOP == 1) && (parityBOT == 0) ){
		
		//Strand Order: S1-S3r-S2r-S4
		barsCombined[1]=barsTOP[0]+barsBOT[0];
		barsCombined[2]=barsTOP[1]+barsBOT[0];
		barsCombined[3]=barsTOP[1]+barsBOT[1];
		
		//Index Shifts:
		//All crossings on S3 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by strand3MaxIndex.
		//Previously encountered crossings on S4 by strand1MaxIndex; previously un-encountered crossings on S4 by strand2MaxIndex.
		
		//Orientation Reversals:
		//Strands S3 and S2 are reversed. The order the strands are encountered in is not changed.
		//Any TOP crossing on both S1 and S2 has its orientation reversed; likewise, any BOT crossings on both S3 and S4 has its orientation reversed.
		
		
		//REVERSAL CONSTRUCTION:
		//Initilize the reversals of S2 and S3.
		int strand2reverse[strand2Length];
		int strand3reverse[strand3Length];
		int orientedSignStrand2reverse[strand2Length];
		int orientedSignStrand3reverse[strand3Length];
		//List the crossings in the opposite order on each strand:
		for(int i=strand2Length; i>0; i--){
			strand2reverse[i-1]=strand2[strand2Length-i];
			orientedSignStrand2reverse[i-1]=orientedSignStrand2[strand2Length-i];
		}
		for(int i=strand3Length; i>0; i--){
			strand3reverse[i-1]=strand3[strand3Length-i];
			orientedSignStrand3reverse[i-1]=orientedSignStrand3[strand3Length-i];
		}
		
		
		//DEBUG:
		/*
		printf("\n\n BEFORE RE-LABEL");
		printf("\n strand2reverse: [");
		for(int i=0;i<strand2Length;i++){
			printf(" %d ", strand2reverse[i]);
		}
		printf("]");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]");
		printf("\n strand4: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4[i]);
		}
		printf("]\n");
		*/
		
		
		//Since the crossings on each of S2 and S3 are now encountered in the opposite order, the indices in the reverse arrays must be re-labeled accordingly.
		//Arrays to specfiy if an entry of the reverse arrays has been modified or not; used to avoid undoing the updated labeling.
		int modifyCheckS2r[strand2Length];
		int modifyCheckS3r[strand3Length];
		int modifyCheckS4[strand4Length];
		//Initialize both arrays as only entires of 0. Each time a label is updated, change entry to 1.
		for(int i=0; i<strand2Length; i++){
			modifyCheckS2r[i]=0;
		}
		for(int i=0; i<strand3Length; i++){
			modifyCheckS3r[i]=0;
		}
		for(int i=0; i<strand4Length; i++){
			modifyCheckS4[i]=0;
		}
		//The crossing indices must be modified so that they are increasing while traveling in reverse order.
		//Iterate through each entry in the strand and relabel the crossings as encountered.
		//Check if any other entries in the strand match this one; update those entries correspondingly.
		//Check against the modifyCheck array to avoid modifying an entry that has already been accounted for.
		
		//First relabel the crossing indices on S2r. Any crossing which is on both S1 and S2 has already been accounted for.
		//Thus, we need only consider those crossing indices which are greater than strand1MaxIndex.
		int S2crossCount=strand1MaxIndex;
		int storeCross;
		for(int i=0; i<strand2Length; i++){
			//We only need to consider those crossings which are BOTH not shared with S1 and which have not already been modified.
			if( (modifyCheckS2r[i] == 0) && ( abs(strand2reverse[i]) > strand1MaxIndex ) ){
				S2crossCount++;
				storeCross=strand2reverse[i];
				if(storeCross>0){
					strand2reverse[i]=S2crossCount;
				} else {
					strand2reverse[i]=(-1*S2crossCount);
				}
				modifyCheckS2r[i]=1;
				//After relabeling this entry, check for corresponding entries later along S2r.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand2Length; j++){
					if ( (abs(strand2reverse[j]) == abs(storeCross)) && (modifyCheckS2r[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand2reverse[j]=(-1*strand2reverse[i]);
						modifyCheckS2r[j]=1;
					}
				}
			}
		}
		//At this point, all crossings indices on S2r not shared with S1 are relabled to be increasing.
		
		//NOTE: the order to run these loops is important. In this case, we travel along S3r, then S4.
		//The same crossing index could show up either twice on one strand, or once on both strands; both possibilities must be accounted for.
		//Iterate first along S3r, whenever an entry is relabeled, check for corresponding entries both on S3r AND on S4 and update these as well.
		//Since S4 is not reversed, only the S4 crossings shared with S3r need to be updated. The remaining S4 crossings stay as is.
		int reverseCrossCount=0;
		for(int i=0; i<strand3Length; i++){
			//If modifyCheckS3r[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS3r[i] == 0){
				reverseCrossCount++;
				storeCross=strand3reverse[i];
				if(storeCross>0){
					strand3reverse[i]=reverseCrossCount;
				} else {
					strand3reverse[i]=(-1*reverseCrossCount);
				}
				modifyCheckS3r[i]=1;
				//After relabeling this entry, check for corresponding entries later along S3r.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand3Length; j++){
					if ( (abs(strand3reverse[j]) == abs(storeCross)) && (modifyCheckS3r[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand3reverse[j]=(-1*strand3reverse[i]);
						modifyCheckS3r[j]=1;
					}
				}
				//Next, check for later corresponding entries along S4 with the same conventions as above.
				for(int j=0; j<strand4Length; j++){
					if( (abs(strand4[j]) == abs(storeCross)) && (modifyCheckS4[j] == 0) ){
						strand4[j]=(-1*strand3reverse[i]);
						modifyCheckS4[j]=1;
					}
				}
			}
		}
		//The above for-loop accounts for all crossings along S3r, including those which also show up on S4. We do not need to modify any other crossings on S4.
		//At this point, all labels of the strands in the BOT tangle are updated.
		
		//DEBUG:
		/*
		printf("\n AFTER RE-LABEL");
		printf("\n strand2reverse: [");
		for(int i=0;i<strand2Length;i++){
			printf(" %d ", strand2reverse[i]);
		}
		printf("]");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]");
		printf("\n strand4: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4[i]);
		}
		printf("]\n");
		*/
		
		
		//Since S2 was reversed but S1 was not, any crossing which was on both S1 and S2 has its orientation reversed. This must be accounted for in both orientedSign arrays.
		//Is is sufficient to use a nested for-loop to search for any crossing indices occuring on both strands, and update the corresponding orientedSign entries.
		for(int i=0; i<strand1Length; i++){
			for(int j=0; j<strand2Length; j++){
				//Each crossing index in the TOP tangle occurs twice; if this occurence happens to be on the two different strands, then the orientation is reversed accordingly.
				if( abs(strand1[i]) == abs(strand2reverse[j]) ){
					orientedSignStrand1[i] = (-1*orientedSignStrand1[i]);
					orientedSignStrand2reverse[j] = (-1*orientedSignStrand2reverse[j]);
				}
			}
		}
		//Likewise, since S3 was reversed but S4 was not, we must again account for this in the orientedSign arrays in the same way as above.
		for(int i=0; i<strand3Length; i++){
			for(int j=0; j<strand4Length; j++){
				//Each crossing index in the BOT tangle occurs twice; if this occurence happens to be on the two different strands, then the orientation is reversed accordingly.
				if( abs(strand3reverse[i]) == abs(strand4[j]) ){
					orientedSignStrand3reverse[i] = (-1*orientedSignStrand3reverse[i]);
					orientedSignStrand4[j] = (-1*orientedSignStrand4[j]);
				}
			}
		}
		//If the orientation has been changed for any crossings on S1, we need to update the orientedSignGaussPROD[] array since this was earlier set equal to orientedSignStrand1.
		for(int i=0; i<strand1Length; i++){
			orientedSignGaussPROD[i]=orientedSignStrand1[i];
		}
		
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussPROD[i] = shiftIndexMagnitude(strand3reverse[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussPROD[i]=orientedSignStrand3reverse[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2reverse[i-barsCombined[1]]) <= strand1MaxIndex ){
				gaussPROD[i] = strand2reverse[i-barsCombined[1]];
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand2reverse[i-barsCombined[1]],strand3MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand2reverse[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note that crossings already encountered on S3r are shifted to match.
			if( abs(strand4[i-barsCombined[2]]) <= strand3MaxIndex ){
				gaussPROD[i] = shiftIndexMagnitude(strand4[i-barsCombined[2]],strand1MaxIndex);
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand4[i-barsCombined[2]],strand2MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand4[i-barsCombined[2]];
		}
		
		parityPROD=0;
		barsPROD[0]=barsCombined[2];
		barsPROD[1]=barsCombined[3];
		
		
	//CASE 5:
	} else if( (parityTOP == 1) && (parityBOT == 1) ){
		
		//Strand Order: S1-S4-S3r-S2r
		barsCombined[1]=barsTOP[0]+(barsBOT[1]-barsBOT[0]);
		barsCombined[2]=barsTOP[0]+barsBOT[1];
		barsCombined[3]=barsTOP[1]+barsBOT[1];
		
		//Index Shifts:
		//All crossings on S3 and S4 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by strand4MaxIndex.
		
		//Orientation Reversals:
		//Strands S3 and S2 are reversed, and S4 is now encountered before S3.
		//Any TOP crossing on both S1 and S2 has its orientation reversed; likewise, any BOT crossing on both S3 and S4 has its orientation reversed.
		
		
		//REVERSAL CONSTRUCTION:
		//Initilize the reversals of S2 and S3.
		int strand2reverse[strand2Length];
		int strand3reverse[strand3Length];
		int orientedSignStrand2reverse[strand2Length];
		int orientedSignStrand3reverse[strand3Length];
		//List the crossings in the opposite order on each strand:
		for(int i=strand2Length; i>0; i--){
			strand2reverse[i-1]=strand2[strand2Length-i];
			orientedSignStrand2reverse[i-1]=orientedSignStrand2[strand2Length-i];
		}
		for(int i=strand3Length; i>0; i--){
			strand3reverse[i-1]=strand3[strand3Length-i];
			orientedSignStrand3reverse[i-1]=orientedSignStrand3[strand3Length-i];
		}
		
		
		//DEBUG:
		/*
		printf("\n\n BEFORE RE-LABEL");
		printf("\n strand2reverse: [");
		for(int i=0;i<strand2Length;i++){
			printf(" %d ", strand2reverse[i]);
		}
		printf("]");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]");
		printf("\n strand4: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4[i]);
		}
		printf("]\n");
		*/
		
		
		//Since the crossings on each of S2 and S3 are now encountered in the opposite order, the indices in the reverse arrays must be relabeled accordingly.
		//Arrays to specfiy if an entry of the reverse arrays has been modified or not; used to avoid undoing the updated labeling.
		int modifyCheckS2r[strand2Length];
		int modifyCheckS3r[strand3Length];
		int modifyCheckS4[strand4Length];
		//Initialize both arrays as only entires of 0. Each time a label is updated, change entry to 1.
		for(int i=0; i<strand2Length; i++){
			modifyCheckS2r[i]=0;
		}
		for(int i=0; i<strand3Length; i++){
			modifyCheckS3r[i]=0;
		}
		for(int i=0; i<strand4Length; i++){
			modifyCheckS4[i]=0;
		}
		//The crossing indices must be modified so that they are increasing while traveling in reverse order.
		//Iterate through each entry in the strand and relabel the crossings as encountered.
		//Check if any other entries in the strand match this one; update those entries correspondingly.
		//Check against the modifyCheck array to avoid modifying an entry that has already been accounted for.
		
		//First relabel the crossing indices on S2r. Any crossing which is on both S1 and S2 has already been accounted for.
		//Thus, we need only consider those crossing indices which are greater than strand1MaxIndex.
		int S2crossCount=strand1MaxIndex;
		int storeCross;
		for(int i=0; i<strand2Length; i++){
			//We only need to consider those crossings which are BOTH not shared with S1 and which have not already been modified.
			if( (modifyCheckS2r[i] == 0) && ( abs(strand2reverse[i]) > strand1MaxIndex ) ){
				S2crossCount++;
				storeCross=strand2reverse[i];
				if(storeCross>0){
					strand2reverse[i]=S2crossCount;
				} else {
					strand2reverse[i]=(-1*S2crossCount);
				}
				modifyCheckS2r[i]=1;
				//After relabeling this entry, check for corresponding entries later along S2r.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand2Length; j++){
					if ( (abs(strand2reverse[j]) == abs(storeCross)) && (modifyCheckS2r[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand2reverse[j]=(-1*strand2reverse[i]);
						modifyCheckS2r[j]=1;
					}
				}
			}
		}
		//At this point, all crossings indices on S2r not shared with S1 are relabled to be increasing.
		
		//NOTE: the order to run these loops is important. In this case, we travel along S4, then S3r.
		//The same crossing index could show up either twice on one strand, or once on both strands; both possibilities must be accounted for.
		//Iterate first along S4, whenever an entry is relabeled, check for corresponding entries both on S4 AND on S3r and update these as well.
		//Next, iterate along S3r and update the index of any remaining crossings not previously encountered on S4.
		int reverseCrossCount=0;
		for(int i=0; i<strand4Length; i++){
			//If modifyCheckS4[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS4[i] == 0){
				reverseCrossCount++;
				storeCross=strand4[i];
				if(storeCross>0){
					strand4[i]=reverseCrossCount;
				} else {
					strand4[i]=(-1*reverseCrossCount);
				}
				modifyCheckS4[i]=1;
				//After relabeling this entry, check for corresponding entries later along S4.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand4Length; j++){
					if ( (abs(strand4[j]) == abs(storeCross)) && (modifyCheckS4[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand4[j]=(-1*strand4[i]);
						modifyCheckS4[j]=1;
					}
				}
				//Next, check for later corresponding entries along S3r with the same conventions as above.
				for(int j=0; j<strand3Length; j++){
					if( (abs(strand3reverse[j]) == abs(storeCross)) && (modifyCheckS3r[j] == 0) ){
						strand3reverse[j]=(-1*strand4[i]);
						modifyCheckS3r[j]=1;
					}
				}
			}
		}
		//The above for-loop accounts for all crossings along S4, including those which also show up on S3r.
		//Next, we iterate through S3r and update the label of any crossings not encountered on S4.
		for(int i=0; i<strand3Length; i++){
			//If modifyCheckS3r[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS3r[i] == 0){
				reverseCrossCount++;
				storeCross=strand3reverse[i];
				if(storeCross>0){
					strand3reverse[i]=reverseCrossCount;
				} else {
					strand3reverse[i]=(-1*reverseCrossCount);
				}
				modifyCheckS3r[i]=1;
				//After relabeling this entry, check for corresponding entries later along S3r.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand3Length; j++){
					if ( (abs(strand3reverse[j]) == abs(storeCross)) && (modifyCheckS3r[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand3reverse[j]=(-1*strand3reverse[i]);
						modifyCheckS3r[j]=1;
					}
				}
			}
		}
		//At this point, all labels of the strands in the BOT tangle are updated.
		
		//DEBUG:
		/*
		printf("\n AFTER RE-LABEL");
		printf("\n strand2reverse: [");
		for(int i=0;i<strand2Length;i++){
			printf(" %d ", strand2reverse[i]);
		}
		printf("]");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]");
		printf("\n strand4: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4[i]);
		}
		printf("]\n");
		*/
		
		
		//Since S2 was reversed but S1 was not, any crossing which was on both S1 and S2 has its orientation reversed. This must be accounted for in both orientedSign arrays.
		//It is sufficient to use a nested for-loop to search for any crossing indices occuring on both strands, and update the corresponding orientedSign entries.
		for(int i=0; i<strand1Length; i++){
			for(int j=0; j<strand2Length; j++){
				//Each crossing index in the TOP tangle occurs twice; if this occurence happens to be on the two different strands, then the orientation is reversed accordingly.
				if( abs(strand1[i]) == abs(strand2reverse[j]) ){
					orientedSignStrand1[i] = (-1*orientedSignStrand1[i]);
					orientedSignStrand2reverse[j] = (-1*orientedSignStrand2reverse[j]);
				}
			}
		}
		//Likewise, since S3 was reversed but S4 was not, we must again account for this in the orientedSign arrays in the same way as above.
		for(int i=0; i<strand3Length; i++){
			for(int j=0; j<strand4Length; j++){
				//Each crossing index in the BOT tangle occurs twice; if this occurence happens to be on the two different strands, then the orientation is reversed accordingly.
				if( abs(strand3reverse[i]) == abs(strand4[j]) ){
					orientedSignStrand3reverse[i] = (-1*orientedSignStrand3reverse[i]);
					orientedSignStrand4[j] = (-1*orientedSignStrand4[j]);
				}
			}
		}
		//If the orientation has been changed for any crossings on S1, we need to update the orientedSignGaussPROD[] array since this was earlier set equal to orientedSignStrand1.
		for(int i=0; i<strand1Length; i++){
			orientedSignGaussPROD[i]=orientedSignStrand1[i];
		}
		
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussPROD[i] = shiftIndexMagnitude(strand4[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussPROD[i]=orientedSignStrand4[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			gaussPROD[i] = shiftIndexMagnitude(strand3reverse[i-barsCombined[1]],strand1MaxIndex);
			orientedSignGaussPROD[i]=orientedSignStrand3reverse[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2reverse[i-barsCombined[2]]) <= strand1MaxIndex ){
				gaussPROD[i] = strand2reverse[i-barsCombined[2]];
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand2reverse[i-barsCombined[2]],strand4MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand2reverse[i-barsCombined[2]];
		}
		
		parityPROD=2;
		barsPROD[0]=barsCombined[1];
		barsPROD[1]=barsCombined[3];
		
		
	//CASE 6:	
	} else if( (parityTOP == 1) && (parityBOT == 2) ){

		//Strand Order: S1-S4r-S2-S3
		barsCombined[1]=barsTOP[0]+(barsBOT[1]-barsBOT[0]);
		barsCombined[2]=barsTOP[1]+(barsBOT[1]-barsBOT[0]);
		barsCombined[3]=barsTOP[1]+barsBOT[1];
		
		//Index Shifts:
		//All crossings S4 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by RELABELEDstrand4MaxIndex.
		//Previously encountered crossigns on S3 by strand1MaxIndex; previously un-encountered crossings on S3 by strand2maxIndex.
		
		//Orientation Reversals:
		//Only strand S4 is reversed, and S4 is now encountered before S3.
		//Any BOT crossing on both S3 and S4 has its orientation reversed.
		
		
		//REVERSAL CONSTRUCTION:
		//Initilize the reversal of S4.
		int strand4reverse[strand4Length];
		int orientedSignStrand4reverse[strand4Length];
		//List the crossings in the opposite order on each strand:
		for(int i=strand4Length; i>0; i--){
			strand4reverse[i-1]=strand4[strand4Length-i];
			orientedSignStrand4reverse[i-1]=orientedSignStrand4[strand4Length-i];
		}
		
		
		//DEBUG:
		/*
		printf("\n\n BEFORE RE-LABEL");
		printf("\n strand3: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3[i]);
		}
		printf("]");
		printf("\n strand4reverse: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4reverse[i]);
		}
		printf("]\n");
		*/
		
		
		//Since the crossings on S4 are now encountered in the opposite order, the indices in the reverse array must be relabeled accordingly.
		//Arrays to specfiy if an entry of the reverse arrays has been modified or not; used to avoid undoing the updated labeling.
		int modifyCheckS3[strand3Length];
		int modifyCheckS4r[strand4Length];
		//Initialize both arrays as only entires of 0. Each time a label is updated, change entry to 1.
		for(int i=0; i<strand3Length; i++){
			modifyCheckS3[i]=0;
		}
		for(int i=0; i<strand4Length; i++){
			modifyCheckS4r[i]=0;
		}
		//The crossing indices must be modified so that they are increasing while traveling in reverse order.
		//Iterate through each entry in the strand and relabel the crossings as encountered.
		//Check if any other entries in the strand match this one; update those entries correspondingly.
		//Check against the modifyCheck array to avoid modifying an entry that has already been accounted for.
		
		//NOTE: the order to run these loops is important. In this case, we travel along S4r, then S3.
		//The same crossing index could show up either twice on one strand, or once on both strands; both possibilities must be accounted for.
		//Iterate first along S4r, whenever an entry is relabeled, check for corresponding entries both on S4r AND on S3 and update these as well.
		//Next, iterate along S3 and update the index of any remaining crossings not previously encountered on S4r.
		int reverseCrossCount=0;
		int storeCross;
		for(int i=0; i<strand4Length; i++){
			//If modifyCheckS4r[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS4r[i] == 0){
				reverseCrossCount++;
				storeCross=strand4reverse[i];
				if(storeCross>0){
					strand4reverse[i]=reverseCrossCount;
				} else {
					strand4reverse[i]=(-1*reverseCrossCount);
				}
				modifyCheckS4r[i]=1;
				//After relabeling this entry, check for corresponding entries later along S4r.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand4Length; j++){
					if ( (abs(strand4reverse[j]) == abs(storeCross)) && (modifyCheckS4r[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand4reverse[j]=(-1*strand4reverse[i]);
						modifyCheckS4r[j]=1;
					}
				}
				//Next, check for later corresponding entries along S3 with the same conventions as above.
				for(int j=0; j<strand3Length; j++){
					if( (abs(strand3[j]) == abs(storeCross)) && (modifyCheckS3[j] == 0) ){
						strand3[j]=(-1*strand4reverse[i]);
						modifyCheckS3[j]=1;
					}
				}
			}
		}
		//The above for-loop accounts for all crossings along S4r, including those which also show up on S3.
		//Next, we iterate through S3 and update the label of any crossings not encountered on S4r.
		for(int i=0; i<strand3Length; i++){
			//If modifyCheckS3[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS3[i] == 0){
				reverseCrossCount++;
				storeCross=strand3[i];
				if(storeCross>0){
					strand3[i]=reverseCrossCount;
				} else {
					strand3[i]=(-1*reverseCrossCount);
				}
				modifyCheckS3[i]=1;
				//After relabeling this entry, check for corresponding entries later along S3.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand3Length; j++){
					if ( (abs(strand3[j]) == abs(storeCross)) && (modifyCheckS3[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand3[j]=(-1*strand3[i]);
						modifyCheckS3[j]=1;
					}
				}
			}
		}
		//At this point, all labels of the crossings along the strands in the BOT tangle are updated to travel first along S4r and then along S3.
		//In order to correctly shift the crossing index when building the combined tangle, we also need to record the maximum crossing index encountered on the relabeled S4r.
		int RELABELEDstrand4MaxIndex=0;
		for(int i=0; i<strand4Length; i++){
			if( abs(strand4reverse[i]) > RELABELEDstrand4MaxIndex ){
				RELABELEDstrand4MaxIndex = abs(strand4reverse[i]);
			}
		}
		
		//DEBUG:
		/*
		printf("\n BEFORE RE-LABEL");
		printf("\n strand3: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3[i]);
		}
		printf("]");
		printf("\n strand4reverse: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4reverse[i]);
		}
		printf("]\n");
		*/
		
		
		//Since S4 was reversed but S3 was not, any crossing which was on both S3 and S4 has its orientation reversed. This must be accounted for in both orientedSign arrays.
		//It is sufficient to use a nested for-loop to search for any crossing indices occuring on both strands, and update the corresponding orientedSign entries.
		for(int i=0; i<strand3Length; i++){
			for(int j=0; j<strand4Length; j++){
				//Each crossing index in the BOT tangle occurs twice; if this occurence happens to be on the two different strands, then the orientation is reversed accordingly.
				if( abs(strand3[i]) == abs(strand4reverse[j]) ){
					orientedSignStrand3[i] = (-1*orientedSignStrand3[i]);
					orientedSignStrand4reverse[j] = (-1*orientedSignStrand4reverse[j]);
				}
			}
		}
		
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussPROD[i] = shiftIndexMagnitude(strand4reverse[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussPROD[i]=orientedSignStrand4reverse[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2[i-barsCombined[1]]) <= strand1MaxIndex ){
				gaussPROD[i] = strand2[i-barsCombined[1]];
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand2[i-barsCombined[1]],RELABELEDstrand4MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand2[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note that we shift the index of any crossing previously encountered on S4r by strand1MaxIndex to be consistent.
			//The other crossings are shifted by strand2MaxIndex.
			if( abs(strand3[i-barsCombined[2]]) <= RELABELEDstrand4MaxIndex ){
				gaussPROD[i] = shiftIndexMagnitude(strand3[i-barsCombined[2]],strand1MaxIndex);
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand3[i-barsCombined[2]],strand2MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand3[i-barsCombined[2]];
		}
		
		parityPROD=1;
		barsPROD[0]=barsCombined[1];
		barsPROD[1]=barsCombined[3];
		
		
	//CASE 7:	
	} else if( (parityTOP == 2) && (parityBOT == 0) ){
		
		//Strand Order: S1-S3-S2-S4
		barsCombined[1]=barsTOP[0]+barsBOT[0];
		barsCombined[2]=barsTOP[1]+barsBOT[0];
		barsCombined[3]=barsTOP[1]+barsBOT[1];
		
		//Index Shifts:
		//All crossings on S3 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by strand3MaxIndex.
		//Previously encountered crossings on S4 by strand1MaxIndex; previously un-encountered crossings on S4 by strand2MaxIndex.
		
		//Orientation Reversals:
		//None.
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussPROD[i] = shiftIndexMagnitude(strand3[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussPROD[i]=orientedSignStrand3[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2[i-barsCombined[1]]) <= strand1MaxIndex ){
				gaussPROD[i] = strand2[i-barsCombined[1]];
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand2[i-barsCombined[1]],strand3MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand2[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note that we shift any crossings on S4 previously encountered on S3 by strand1MaxIndex to be consistent.
			//Other crossings on S4 not previously encountered on S3 are shifted by strand2MaxIndex instead.
			if( abs(strand4[i-barsCombined[2]]) <= strand3MaxIndex ){
				gaussPROD[i] = shiftIndexMagnitude(strand4[i-barsCombined[2]],strand1MaxIndex);
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand4[i-barsCombined[2]],strand2MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand4[i-barsCombined[2]];
		}
		
		parityPROD=0;
		barsPROD[0]=barsCombined[2];
		barsPROD[1]=barsCombined[3];
		
		
	//CASE 8:	
	} else if( (parityTOP == 2) && (parityBOT == 1) ){
		
		//Strand Order: S1-S3-S2r-S4
		barsCombined[1]=barsTOP[0]+barsBOT[0];
		barsCombined[2]=barsTOP[1]+barsBOT[0];
		barsCombined[3]=barsTOP[1]+barsBOT[1];
		
		//Index Shifts:
		//All crossings on S3 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by strand3MaxIndex.
		//Previously encountered crossings on S4 by strand1MaxIndex; previously un-encountered crossings on S4 by strand2MaxIndex.
		
		//Orientation Reversals:
		//Only strand S2 has its orientation reversed. Any TOP crossing occuring on both S1 and S2 has its orientation reversed.
		
		
		//REVERSAL CONSTRUCTION:
		//Initilize the reversal of S2.
		int strand2reverse[strand2Length];
		int orientedSignStrand2reverse[strand2Length];
		//List the crossings in the opposite order on the strand:
		for(int i=strand2Length; i>0; i--){
			strand2reverse[i-1]=strand2[strand2Length-i];
			orientedSignStrand2reverse[i-1]=orientedSignStrand2[strand2Length-i];
		}
		
		
		//DEBUG:
		/*
		printf("\n\n BEFORE RE-LABEL");
		printf("\n strand2reverse: [");
		for(int i=0;i<strand2Length;i++){
			printf(" %d ", strand2reverse[i]);
		}
		printf("]\n");
		*/
		
		
		//Since the crossings on S2 are now encountered in the opposite order, the indices in the reverse array must be relabeled accordingly.
		//Arrays to specfiy if an entry of the reverse arrays has been modified or not; used to avoid undoing the updated labeling.
		int modifyCheckS2r[strand2Length];
		//Initialize both arrays as only entires of 0. Each time a label is updated, change entry to 1.
		for(int i=0; i<strand2Length; i++){
			modifyCheckS2r[i]=0;
		}
		//The crossing indices must be modified so that they are increasing while traveling in reverse order.
		//Iterate through each entry in the strand and relabel the crossings as encountered.
		//Check if any other entries in the strand match this one; update those entries correspondingly.
		//Check against the modifyCheck array to avoid modifying an entry that has already been accounted for.
		
		//First relabel the crossing indices on S2r. Any crossing which is on both S1 and S2 has already been accounted for.
		//Thus, we need only consider those crossing indices which are greater than strand1MaxIndex.
		int S2crossCount=strand1MaxIndex;
		int storeCross;
		for(int i=0; i<strand2Length; i++){
			//We only need to consider those crossings which are BOTH not shared with S1 and which have not already been modified.
			if( (modifyCheckS2r[i] == 0) && ( abs(strand2reverse[i]) > strand1MaxIndex ) ){
				S2crossCount++;
				storeCross=strand2reverse[i];
				if(storeCross>0){
					strand2reverse[i]=S2crossCount;
				} else {
					strand2reverse[i]=(-1*S2crossCount);
				}
				modifyCheckS2r[i]=1;
				//After relabeling this entry, check for corresponding entries later along S2r.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand2Length; j++){
					if ( (abs(strand2reverse[j]) == abs(storeCross)) && (modifyCheckS2r[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand2reverse[j]=(-1*strand2reverse[i]);
						modifyCheckS2r[j]=1;
					}
				}
			}
		}
		//At this point, all crossings indices on S2r not shared with S1 are relabled to be increasing.
		
		
		//DEBUG:
		/*
		printf("\n AFTER RE-LABEL");
		printf("\n strand2reverse: [");
		for(int i=0;i<strand2Length;i++){
			printf(" %d ", strand2reverse[i]);
		}
		printf("]\n");
		*/
		
		
		//Since S2 was reversed but S1 was not, any crossing which was on both S1 and S2 has its orientation reversed. This must be accounted for in both orientedSign arrays.
		//It is sufficient to use a nested for-loop to search for any crossing indices occuring on both strands, and update the corresponding orientedSign entries.
		for(int i=0; i<strand1Length; i++){
			for(int j=0; j<strand2Length; j++){
				//Each crossing index in the TOP tangle occurs twice; if this occurence happens to be on the two different strands, then the orientation is reversed accordingly.
				if( abs(strand1[i]) == abs(strand2reverse[j]) ){
					orientedSignStrand1[i] = (-1*orientedSignStrand1[i]);
					orientedSignStrand2reverse[j] = (-1*orientedSignStrand2reverse[j]);
				}
			}
		}
		//If the orientation has been changed for any crossings on S1, we need to update the orientedSignGaussPROD[] array since this was earlier set equal to orientedSignStrand1.
		for(int i=0; i<strand1Length; i++){
			orientedSignGaussPROD[i]=orientedSignStrand1[i];
		}
		
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussPROD[i] = shiftIndexMagnitude(strand3[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussPROD[i]=orientedSignStrand3[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2reverse[i-barsCombined[1]]) <= strand1MaxIndex ){
				gaussPROD[i] = strand2reverse[i-barsCombined[1]];
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand2reverse[i-barsCombined[1]],strand3MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand2reverse[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note that we shift any crossings on S4 previously encountered on S3 by strand1MaxIndex to be consistent.
			//Other crossings on S4 not previously encountered on S3 are shifted by strand2MaxIndex instead.
			if( abs(strand4[i-barsCombined[2]]) <= strand3MaxIndex ){
				gaussPROD[i] = shiftIndexMagnitude(strand4[i-barsCombined[2]],strand1MaxIndex);
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand4[i-barsCombined[2]],strand2MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand4[i-barsCombined[2]];
		}
		
		parityPROD=1;
		barsPROD[0]=barsCombined[1];
		barsPROD[1]=barsCombined[3];
		
		
	//CASE 9:	
	} else if( (parityTOP == 2) && (parityBOT == 2) ){
		
		//Strand Order: S1-S3-S4-S2
		barsCombined[1]=barsTOP[0]+barsBOT[0];
		barsCombined[2]=barsTOP[0]+barsBOT[1];
		barsCombined[3]=barsTOP[1]+barsBOT[1];
		
		//Index Shifts:
		//All crossings on S3 and S4 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by strand4MaxIndex.
		
		//Orientation Reversals:
		//None.
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussPROD[i] = shiftIndexMagnitude(strand3[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussPROD[i]=orientedSignStrand3[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			gaussPROD[i] = shiftIndexMagnitude(strand4[i-barsCombined[1]],strand1MaxIndex);
			orientedSignGaussPROD[i]=orientedSignStrand4[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2[i-barsCombined[2]]) <= strand1MaxIndex ){
				gaussPROD[i] = strand2[i-barsCombined[2]];
			} else {
				gaussPROD[i] = shiftIndexMagnitude(strand2[i-barsCombined[2]],strand4MaxIndex);
			}
			orientedSignGaussPROD[i]=orientedSignStrand2[i-barsCombined[2]];
		}
		
		parityPROD=2;
		barsPROD[0]=barsCombined[1];
		barsPROD[1]=barsCombined[3];
		
		
	//ERROR CASE:	
	} else {
		//This is only possible if either of the two input tangles has parity other than 0, 1, or 2 (infinity). This is not allowed, so the function stops here.
		printf("\n ERROR: It appears that one (or both) of the input tangles has parity other than 0, 1, or infinity.\n");
		return;
	}
	
	
	//DEBUG:
	/*
	printf("\n barsCombined: [");
	for(int i=0;i<4;i++){
		printf(" %d ", barsCombined[i]);
	}
	printf("]\n");
	printf("\n gaussPROD: [");
	for(int i=0; i<(2*numOfCrossingsPROD); i++){
		printf(" %d ", gaussPROD[i]);
	}
	printf("]");
	printf("\n orientedSignGaussPROD: [");
	for(int i=0; i<(2*numOfCrossingsPROD); i++){
		printf(" %d ", orientedSignGaussPROD[i]);
	}
	printf("]");
	printf("\n barsPROD: [ %d %d ]", barsPROD[0],barsPROD[1]);
	printf("\n parityPROD: %d \n", parityPROD);
	*/
	
	
}




/*
(NC, 7/26/18)
This function constructs the Gauss code of a sum of two tangles using the Gauss code of the left and right tangles. There are nine possible cases depending on the parities of the two tangles being added.
The Gauss code of both the input and output tangles is stored across three arrays: gauss[], bars[], and orientedSignGauss[]. The number of crossings and tangle parity must also be accounted for.
The function is only intended for use with 2-string tangles. The two tangles being added are destinguished by LEFT and RIGHT, and the result is referred to as SUM.

Function inputs:
numOfCrossingsLEFT: the number of crossings of the LEFT tangle.
numOfCrossingsRIGHT: the number of crossings of the RIGHT tangle.
gaussLEFT[]: an array storing the (ordered) Gauss code of the LEFT tangle. Each entry contains the index of the crossing encountered. Positive entries denote an over-crossing; negative entries denote an under-crossing.
gaussRIGHT[]: an array storing the (ordered) Gauss code the RIGHT tangle, with the same convention described above.
orientedSignGaussLEFT[]: an array of the same length as gaussLEFT[] containing information on the orientation of each crossing encountered (either 1 or -1). Each entry matches the corresponding entry in gaussLEFT[].
orientedSignGaussRIGHT[]: an array of the same length as gaussRIGHT[] containing the information on the orientation of each crossing encountered, with the same convention as described above.
barsLEFT[]: an array with 2 entries to denote where each of the two strings exits the LEFT tangle. The first/second entry denotes the last crossing encountered along the first/second string (regardless of crossing index).
barsRIGHT[]: an array with 2 entries to denote where each of the two strings exits the RIGHT tangle, with the same convention described above.
parityLEFT: the parity of the LEFT tangle, either 0, 1, or 2 (infinity).
parityRIGHT: the parity of the RIGHT tangle, either 0, 1, or 2 (infinity).

Function outputs (in the form of modified global variables):
numOfCrossingsSUM: the number of crossings of the SUM tangle.
gaussSUM[]: an array storing the (ordered) Gauss code the SUM tangle, with the same convention described above.
orientedSignGaussSUM[]: an array of the same length as gaussSUM[] containing the information on the orientation of each crossing encountered, with the same convention as described above.
barsSUM[]: an array with 2 entries to denote where each of the two strings exits the SUM tangle, with the same convention described above.
paritySUM: the parity of the SUM tangle, either 0, 1, or 2 (infinity).


This function calls:
	nothing
This function is called by:
	buildTwistTangle()


*/

void tangleSum(int numOfCrossingsLEFT, int *gaussLEFT, int *barsLEFT, int *orientedSignGaussLEFT, int parityLEFT, int numOfCrossingsRIGHT, int *gaussRIGHT, int *barsRIGHT, int *orientedSignGaussRIGHT, int parityRIGHT, int &numOfCrossingsSUM, int *gaussSUM, int *barsSUM, int *orientedSignGaussSUM, int &paritySUM){
	
	numOfCrossingsSUM = numOfCrossingsLEFT + numOfCrossingsRIGHT;
	
	//Local arrays to store and manipulate the input information.
	int gaussLocalSUM[numOfCrossingsSUM];
	int orientedSignGaussLocalSUM[numOfCrossingsSUM];
	int barsCombined[4];
	
	//Initialize arrays for each of the four strands.
	int firstStrandLEFTlength = barsLEFT[0];
	int secondStrandLEFTlength = barsLEFT[1] - barsLEFT[0];
	int firstStrandRIGHTlength = barsRIGHT[0];
	int secondStrandRIGHTlength = barsRIGHT[1] - barsRIGHT[0];
	int strand1[firstStrandLEFTlength];
	int strand2[secondStrandLEFTlength];
	int strand3[firstStrandRIGHTlength];
	int strand4[secondStrandRIGHTlength];
	int orientedSignStrand1[firstStrandLEFTlength];
	int orientedSignStrand2[secondStrandLEFTlength];
	int orientedSignStrand3[firstStrandRIGHTlength];
	int orientedSignStrand4[secondStrandRIGHTlength];
	
	//Initilize variables to represent the highest crossing index encounterd for each of the four strands.
	//Set the default values to be 0, this allows for strands with no crossings to be considered.
	//These are checked when updating indices of later strands.
	int strand1MaxIndex=0;
	int strand2MaxIndex=0;
	int strand3MaxIndex=0;
	int strand4MaxIndex=0;
	
	//Strand 1: first strand LEFT
	for(int i=0;i<firstStrandLEFTlength;i++){
		strand1[i]=gaussLEFT[i];
		orientedSignStrand1[i]=orientedSignGaussLEFT[i];
		//Update the strand1MaxIndex if needed.
		if( abs(strand1[i]) > strand1MaxIndex ){
			strand1MaxIndex = abs(strand1[i]);
		}
	}
	//Strand 2: second strand LEFT
	for(int i=0;i<secondStrandLEFTlength;i++){
		strand2[i]=gaussLEFT[i+firstStrandLEFTlength];
		orientedSignStrand2[i]=orientedSignGaussLEFT[i+firstStrandLEFTlength];
		//Update the strand2MaxIndex if needed.
		if( abs(strand2[i]) > strand2MaxIndex ){
			strand2MaxIndex = abs(strand2[i]);
		}
	}
	//Strand 3: first strand RIGHT
	for(int i=0;i<firstStrandRIGHTlength;i++){
		strand3[i]=gaussRIGHT[i];
		orientedSignStrand3[i]=orientedSignGaussRIGHT[i];
		//Update the strand3MaxIndex if needed.
		if( abs(strand3[i]) > strand3MaxIndex ){
			strand3MaxIndex = abs(strand3[i]);
		}
	}
	//Strand 4: second strand RIGHT
	for(int i=0;i<secondStrandRIGHTlength;i++){
		strand4[i]=gaussRIGHT[i+firstStrandRIGHTlength];
		orientedSignStrand4[i]=orientedSignGaussRIGHT[i+firstStrandRIGHTlength];
		//Update the strand4MaxIndex if needed.
		if( abs(strand4[i]) > strand4MaxIndex ){
			strand4MaxIndex = abs(strand4[i]);
		}
	}
	
	//Re-define strand lengths for more compact and intuitive notation:
	int strand1Length = firstStrandLEFTlength;
	int strand2Length = secondStrandLEFTlength;
	int strand3Length = firstStrandRIGHTlength;
	int strand4Length = secondStrandRIGHTlength;
	
	
	//DEBUG:
	/*
	printf("\n Strand 1: [");
	for(int i=0; i<strand1Length; i++){
		printf(" %d ", strand1[i]);
	}
	printf("]");
	printf("\n Strand 2: [");
	for(int i=0; i<strand2Length; i++){
		printf(" %d ", strand2[i]);
	}
	printf("]");
	printf("\n Strand 3: [");
	for(int i=0; i<strand3Length; i++){
		printf(" %d ", strand3[i]);
	}
	printf("]");
	printf("\n Strand 4: [");
	for(int i=0; i<strand4Length; i++){
		printf(" %d ", strand4[i]);
	}
	printf("]\n\n");
	
	printf(" Strand 1 Max: %d \n", strand1MaxIndex);
	printf(" Strand 2 Max: %d \n", strand2MaxIndex);
	printf(" Strand 3 Max: %d \n", strand3MaxIndex);
	printf(" Strand 4 Max: %d \n", strand4MaxIndex);
	*/
	
	
	//Each of the nine cases below depends on the parity of the constituent tangles.
	//The resulting sum of tangles is constructed by connecting the four strands defined above in certain order, depending on the case.
	//As the strings are connected, the crossing indices need to be tracked and updated.
	//The array barsCombined[] will be used to track last crossing on each of the four strands, in the order in which they are combined.
	//Each case will be prefaced by the certain order of the four strands.
	//By convetion, Strand 1 (firstStrandLEFT) is always the first of the four strands in the new tangle.
	
	barsCombined[0]=barsLEFT[0];
	for(int i=0; i<strand1Length; i++){
		gaussSUM[i]=strand1[i];
		orientedSignGaussSUM[i]=orientedSignStrand1[i];
	}
	
	//CASE 1:
	if( (parityLEFT == 0) && (parityRIGHT == 0) ){
		
		//Strand Order: S1-S3-S4-S2
		barsCombined[1]=barsLEFT[0]+barsRIGHT[0];
		barsCombined[2]=barsLEFT[0]+barsRIGHT[1];
		barsCombined[3]=barsLEFT[1]+barsRIGHT[1];
		
		//Index Shifts:
		//All crossings on S3 and S4 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by strand4MaxIndex.
		
		//Orientation Reversals:
		//None.
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussSUM[i] = shiftIndexMagnitude(strand3[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussSUM[i]=orientedSignStrand3[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			gaussSUM[i] = shiftIndexMagnitude(strand4[i-barsCombined[1]],strand1MaxIndex);
			orientedSignGaussSUM[i]=orientedSignStrand4[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2[i-barsCombined[2]]) <= strand1MaxIndex ){
				gaussSUM[i] = strand2[i-barsCombined[2]];
			} else {
				gaussSUM[i] = shiftIndexMagnitude(strand2[i-barsCombined[2]],strand4MaxIndex);
			}
			orientedSignGaussSUM[i]=orientedSignStrand2[i-barsCombined[2]];
		}
		
		paritySUM=0;
		barsSUM[0]=barsCombined[1];
		barsSUM[1]=barsCombined[3];
	
	//CASE 2:	
	} else if( (parityLEFT == 0) && (parityRIGHT == 1) ){
		//NOTE: the construction for case2 is the same as for case 1, except for paritySUM.
		
		//Strand Order: S1-S3-S4-S2
		barsCombined[1]=barsLEFT[0]+barsRIGHT[0];
		barsCombined[2]=barsLEFT[0]+barsRIGHT[1];
		barsCombined[3]=barsLEFT[1]+barsRIGHT[1];
		
		//Index Shifts:
		//All crossings on S3 and S4 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by strand4MaxIndex.
		
		//Orientation Reversals:
		//None.
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussSUM[i] = shiftIndexMagnitude(strand3[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussSUM[i]=orientedSignStrand3[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			gaussSUM[i] = shiftIndexMagnitude(strand4[i-barsCombined[1]],strand1MaxIndex);
			orientedSignGaussSUM[i]=orientedSignStrand4[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2[i-barsCombined[2]]) <= strand1MaxIndex ){
				gaussSUM[i] = strand2[i-barsCombined[2]];
			} else {
				gaussSUM[i] = shiftIndexMagnitude(strand2[i-barsCombined[2]],strand4MaxIndex);
			}
			orientedSignGaussSUM[i]=orientedSignStrand2[i-barsCombined[2]];
		}
		
		paritySUM=1;
		barsSUM[0]=barsCombined[1];
		barsSUM[1]=barsCombined[3];
	
	//CASE 3:	
	} else if( (parityLEFT == 0) && (parityRIGHT == 2) ){
		
		//Strand Order: S1-S3-S2-S4
		barsCombined[1]=barsLEFT[0]+barsRIGHT[0];
		barsCombined[2]=barsLEFT[1]+barsRIGHT[0];
		barsCombined[3]=barsLEFT[1]+barsRIGHT[1];
		
		//Index Shifts:
		//All crossings on S3 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by strand3MaxIndex.
		//Previously encountered crossings on S4 by strand1MaxIndex; previously un-encountered crossings on S4 by strand2MaxIndex.
		
		//Orientation Reversals:
		//None.
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussSUM[i] = shiftIndexMagnitude(strand3[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussSUM[i]=orientedSignStrand3[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2[i-barsCombined[1]]) <= strand1MaxIndex ){
				gaussSUM[i] = strand2[i-barsCombined[1]];
			} else {
				gaussSUM[i] = shiftIndexMagnitude(strand2[i-barsCombined[1]],strand3MaxIndex);
			}
			orientedSignGaussSUM[i]=orientedSignStrand2[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note crossings already encountered on strand3 are shifted to match.
			if( abs(strand4[i-barsCombined[2]]) <= strand3MaxIndex ){
				gaussSUM[i] = shiftIndexMagnitude(strand4[i-barsCombined[2]],strand1MaxIndex);
			} else {
				gaussSUM[i] = shiftIndexMagnitude(strand4[i-barsCombined[2]],strand2MaxIndex);
			}
			orientedSignGaussSUM[i]=orientedSignStrand4[i-barsCombined[2]];
		}
		
		paritySUM=2;
		barsSUM[0]=barsCombined[2];
		barsSUM[1]=barsCombined[3];
		
	//CASE 4:	
	} else if( (parityLEFT == 1) && (parityRIGHT == 0) ){
		
		//Strand Order: S1-S4r-S3r-S2
		barsCombined[1]=barsLEFT[0]+(barsRIGHT[1]-barsRIGHT[0]);
		barsCombined[2]=barsLEFT[0]+barsRIGHT[1];
		barsCombined[3]=barsLEFT[1]+barsRIGHT[1];
		
		//Index Shifts:
		//All crossings on S3 and S4 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by strand4MaxIndex.
		
		//Orientation Reversals:
		//The RIGHT tangle is totally reversed: S3 and S4 both have orientations reversed and are encountered in reverse order.
		//Except for updating the indices, orientedSigns for S3 and S4 do not change since all crossing arcs are reversed.
		
		//REVERSAL CONSTRUCTION:
		//Initilize the reversals of S3 and S4:
		int strand3reverse[strand3Length];
		int strand4reverse[strand4Length];
		int orientedSignStrand3reverse[strand3Length];
		int orientedSignStrand4reverse[strand4Length];
		//List the crossings in the opposite order on each strand:
		for(int i=strand3Length; i>0; i--){
			strand3reverse[i-1]=strand3[strand3Length-i];
			orientedSignStrand3reverse[i-1]=orientedSignStrand3[strand3Length-i];
		}
		for(int i=strand4Length; i>0; i--){
			strand4reverse[i-1]=strand4[strand4Length-i];
			orientedSignStrand4reverse[i-1]=orientedSignStrand4[strand4Length-i];
		}
		
		
		//DEBUG:
		/*
		printf("\n\n BEFORE RE-LABEL");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]");
		printf("\n strand4reverse: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4reverse[i]);
		}
		printf("]\n");
		*/
		
		//Since the crossings on the RIGHT tangle are now encountered in the opposite order, the indices in the reverse arrays must be re-labeled accordingly.
		//Arrays to to specfiy if an entry of strand3reverse or strand4reverse has been modified or not; used to avoid undoing the updated labeling.
		int modifyCheckS3r[strand3Length];
		int modifyCheckS4r[strand4Length];
		//Initialize both arrays as only entires of 0. Each time a label is updated, change entry to 1.
		for(int i=0; i<strand3Length; i++){
			modifyCheckS3r[i]=0;
		}
		for(int i=0; i<strand4Length; i++){
			modifyCheckS4r[i]=0;
		}
		//The crossing indices must be modified so that they are increasing while traveling in reverse order.
		//Iterate through each entry in the strand and relabel the crossings as encountered.
		//Check if any other entries in the strand match this one; update those entries correspondingly.
		//Check against the modifyCheck array to avoid modifying an entry that has already been accounted for.
		
		//NOTE: the order to run these loops is important. In this case, we travel along S4r, then S3r.
		//The same crossing index could show up either twice on one strand, or once on both strands; both possibilities must be accounted for.
		//Iterate first alond S4r, whenever an entry is re-labeled, check for corresponding entries both on S4r AND on S3r and update these as well.
		//After finishing S4r, then move to S3r. Any crossing on both strands should already be accounted for; we now need only check along S3r.
		int reverseCrossCount=0;
		int storeCross;
		for(int i=0; i<strand4Length; i++){
			//If modifyCheckS4r[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS4r[i] == 0){
				reverseCrossCount++;
				storeCross=strand4reverse[i];
				if(storeCross>0){
					strand4reverse[i]=reverseCrossCount;
				} else {
					strand4reverse[i]=(-1*reverseCrossCount);
				}
				modifyCheckS4r[i]=1;
				//After relabeling this entry, check for corresponding entries later along S4r.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand4Length; j++){
					if ( (abs(strand4reverse[j]) == abs(storeCross)) && (modifyCheckS4r[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand4reverse[j]=(-1*strand4reverse[i]);
						modifyCheckS4r[j]=1;
					}
				}
				//Next, check for later corresponding entries along S3r with the same conventions as above.
				for(int j=0; j<strand3Length; j++){
					if( (abs(strand3reverse[j]) == abs(storeCross)) && (modifyCheckS3r[j] == 0) ){
						strand3reverse[j]=(-1*strand4reverse[i]);
						modifyCheckS3r[j]=1;
					}
				}
			}
		}
		//The above for-loop accounts for all crossings along S4r, including those which also show up on S3r.
		//Now we move along S3r and update any remaining crossings, with the same conventions as above.
		for(int i=0; i<strand3Length; i++){
			//If modifyCheckS3r[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS3r[i] == 0){
				reverseCrossCount++;
				storeCross=strand3reverse[i];
				if(storeCross>0){
					strand3reverse[i]=reverseCrossCount;
				} else {
					strand3reverse[i]=(-1*reverseCrossCount);
				}
				modifyCheckS3r[i]=1;
				for(int j=0; j<strand3Length; j++){
					if ( (abs(strand3reverse[j]) == abs(storeCross)) && (modifyCheckS3r[j] == 0) ){
						strand3reverse[j]=(-1*strand3reverse[i]);
						modifyCheckS3r[j]=1;
					}
				}
			}
		}
		//At this point, all labels of the strands in the RIGHT tangle are updated to travel in the reverse direction, with increasing crossing index.
		
		//DEBUG:
		/*
		printf("\n AFTER RE-LABEL");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]");
		printf("\n strand4reverse: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4reverse[i]);
		}
		printf("]\n");
		*/
		
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussSUM[i] = shiftIndexMagnitude(strand4reverse[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussSUM[i]=orientedSignStrand4reverse[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			gaussSUM[i] = shiftIndexMagnitude(strand3reverse[i-barsCombined[1]],strand1MaxIndex);
			orientedSignGaussSUM[i]=orientedSignStrand3reverse[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2[i-barsCombined[2]]) <= strand1MaxIndex ){
				gaussSUM[i] = strand2[i-barsCombined[2]];
			} else {
				gaussSUM[i] = shiftIndexMagnitude(strand2[i-barsCombined[2]],strand4MaxIndex);
			}
			orientedSignGaussSUM[i]=orientedSignStrand2[i-barsCombined[2]];
		}
		
		paritySUM=1;
		barsSUM[0]=barsCombined[1];
		barsSUM[1]=barsCombined[3];
		
	//CASE 5:
	} else if( (parityLEFT == 1) && (parityRIGHT == 1) ){
		
		//Strand Order: S1-S4r-S3r-S2
		barsCombined[1]=barsLEFT[0]+(barsRIGHT[1]-barsRIGHT[0]);
		barsCombined[2]=barsLEFT[0]+barsRIGHT[1];
		barsCombined[3]=barsLEFT[1]+barsRIGHT[1];
		
		//Index Shifts:
		//All crossings on S3 and S4 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by strand4MaxIndex.
		
		//Orientation Reversals:
		//The RIGHT tangle is totally reversed: S3 and S4 both have orientations reversed and are encountered in reverse order.
		//Except for updating the indices, orientedSigns for S3 and S4 do not change since all crossing arcs are reversed.
		
		//REVERSAL CONSTRUCTION:
		//Initilize the reversals of S3 and S4:
		int strand3reverse[strand3Length];
		int strand4reverse[strand4Length];
		int orientedSignStrand3reverse[strand3Length];
		int orientedSignStrand4reverse[strand4Length];
		//List the crossings in the opposite order on each strand:
		for(int i=strand3Length; i>0; i--){
			strand3reverse[i-1]=strand3[strand3Length-i];
			orientedSignStrand3reverse[i-1]=orientedSignStrand3[strand3Length-i];
		}
		for(int i=strand4Length; i>0; i--){
			strand4reverse[i-1]=strand4[strand4Length-i];
			orientedSignStrand4reverse[i-1]=orientedSignStrand4[strand4Length-i];
		}
		
		//DEBUG:
		/*
		printf("\n\n BEFORE RE-LABEL");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]");
		printf("\n strand4reverse: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4reverse[i]);
		}
		printf("]\n");
		*/
		
		//Since the crossings on the RIGHT tangle are now encountered in the opposite order, the indices in the reverse arrays must be re-labeled accordingly.
		//Arrays to to specfiy if an entry of strand3reverse or strand4reverse has been modified or not; used to avoid undoing the updated labeling.
		int modifyCheckS3r[strand3Length];
		int modifyCheckS4r[strand4Length];
		//Initialize both arrays as only entires of 0. Each time a label is updated, change entry to 1.
		for(int i=0; i<strand3Length; i++){
			modifyCheckS3r[i]=0;
		}
		for(int i=0; i<strand4Length; i++){
			modifyCheckS4r[i]=0;
		}
		//The crossing indices must be modified so that they are increasing while traveling in reverse order.
		//Iterate through each entry in the strand and relabel the crossings as encountered.
		//Check if any other entries in the strand match this one; update those entries correspondingly.
		//Check against the modifyCheck array to avoid modifying an entry that has already been accounted for.
		
		//NOTE: the order to run these loops is important. In this case, we travel along S4r, then S3r.
		//The same crossing index could show up either twice on one strand, or once on both strands; both possibilities must be accounted for.
		//Iterate first alond S4r, whenever an entry is re-labeled, check for corresponding entries both on S4r AND on S3r and update these as well.
		//After finishing S4r, then move to S3r. Any crossing on both strands should already be accounted for; we now need only check along S3r.
		int reverseCrossCount=0;
		int storeCross;
		for(int i=0; i<strand4Length; i++){
			//If modifyCheckS4r[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS4r[i] == 0){
				reverseCrossCount++;
				storeCross=strand4reverse[i];
				if(storeCross>0){
					strand4reverse[i]=reverseCrossCount;
				} else {
					strand4reverse[i]=(-1*reverseCrossCount);
				}
				modifyCheckS4r[i]=1;
				//After relabeling this entry, check for corresponding entries later along S4r.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand4Length; j++){
					if ( (abs(strand4reverse[j]) == abs(storeCross)) && (modifyCheckS4r[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand4reverse[j]=(-1*strand4reverse[i]);
						modifyCheckS4r[j]=1;
					}
				}
				//Next, check for later corresponding entries along S3r with the same conventions as above.
				for(int j=0; j<strand3Length; j++){
					if( (abs(strand3reverse[j]) == abs(storeCross)) && (modifyCheckS3r[j] == 0) ){
						strand3reverse[j]=(-1*strand4reverse[i]);
						modifyCheckS3r[j]=1;
					}
				}
			}
		}
		//The above for-loop accounts for all crossings along S4r, including those which also show up on S3r.
		//Now we move along S3r and update any remaining crossings, with the same conventions as above.
		for(int i=0; i<strand3Length; i++){
			//If modifyCheckS3r[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS3r[i] == 0){
				reverseCrossCount++;
				storeCross=strand3reverse[i];
				if(storeCross>0){
					strand3reverse[i]=reverseCrossCount;
				} else {
					strand3reverse[i]=(-1*reverseCrossCount);
				}
				modifyCheckS3r[i]=1;
				for(int j=0; j<strand3Length; j++){
					if ( (abs(strand3reverse[j]) == abs(storeCross)) && (modifyCheckS3r[j] == 0) ){
						strand3reverse[j]=(-1*strand3reverse[i]);
						modifyCheckS3r[j]=1;
					}
				}
			}
		}
		//At this point, all labels of the strands in the RIGHT tangle are updated to travel in the reverse direction, with increasing crossing index.
		
		//DEBUG:
		/*
		printf("\n AFTER RE-LABEL");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]");
		printf("\n strand4reverse: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4reverse[i]);
		}
		printf("]\n");
		*/
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussSUM[i] = shiftIndexMagnitude(strand4reverse[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussSUM[i]=orientedSignStrand4reverse[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			gaussSUM[i] = shiftIndexMagnitude(strand3reverse[i-barsCombined[1]],strand1MaxIndex);
			orientedSignGaussSUM[i]=orientedSignStrand3reverse[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2[i-barsCombined[2]]) <= strand1MaxIndex ){
				gaussSUM[i] = strand2[i-barsCombined[2]];
			} else {
				gaussSUM[i] = shiftIndexMagnitude(strand2[i-barsCombined[2]],strand4MaxIndex);
			}
			orientedSignGaussSUM[i]=orientedSignStrand2[i-barsCombined[2]];
		}
		
		paritySUM=0;
		barsSUM[0]=barsCombined[1];
		barsSUM[1]=barsCombined[3];
		
	//CASE 6:	
	} else if( (parityLEFT == 1) && (parityRIGHT == 2) ){
		
		//Strand Order: S1-S3r-S2-S4
		barsCombined[1]=barsLEFT[0]+barsRIGHT[0];
		barsCombined[2]=barsLEFT[1]+barsRIGHT[0];
		barsCombined[3]=barsLEFT[1]+barsRIGHT[1];
		
		//Index Shifts:
		//All crossings on S3 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by strand3MaxIndex.
		//Previously encountered crossings on S4 by strand1MaxIndex; previously un-encountered crossings on S4 by strand2MaxIndex.
		
		//Orientation Reversals:
		//Only S3 is reversed. Any crossings shared by S3 and S4 have their orientations reversed.
		
		//REVERSAL CONSTRUCTION:
		//Initilize the reversal of S3:
		int strand3reverse[strand3Length];
		int orientedSignStrand3reverse[strand3Length];
		//List the crossings in the opposite order on S3:
		for(int i=strand3Length; i>0; i--){
			strand3reverse[i-1]=strand3[strand3Length-i];
			orientedSignStrand3reverse[i-1]=orientedSignStrand3[strand3Length-i];
		}
		
		
		//DEBUG:
		/*
		printf("\n\n BEFORE RE-LABEL");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]\n");
		printf(" strand4: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4[i]);
		}
		printf("]\n");
		*/
		
		//Since the crossings along S3 tangle are now encountered in the opposite order, the indices in the reverse array must be re-labeled accordingly.
		//Array to specfiy if an entry of strand3reverse has been modified or not; used to avoid undoing the updated labeling.
		//Since some crossings are possibly shared between S3 and S4, we need a check array for S4 as well.
		int modifyCheckS3r[strand3Length];
		int modifyCheckS4[strand4Length];
		//Initialize the array with only entires of 0. Each time a label is updated, change entry to 1.
		for(int i=0; i<strand3Length; i++){
			modifyCheckS3r[i]=0;
		}
		for(int i=0; i<strand4Length; i++){
			modifyCheckS4[i]=0;
		}
		//The crossing indices must be modified so that they are increasing while traveling in reverse order.
		//Iterate through each entry in the strand3 and relabel the crossings as encountered.
		//Check if any other entries in the strand match this one; update those entries correspondingly.
		//Check against the modifyCheck array to avoid modifying an entry that has already been accounted for.
		
		//NOTE: the order to run these loops is important. In this case, we travel along S3r, then S4.
		//The same crossing index could show up either twice on one strand, or once on both strands; both possibilities must be accounted for.
		//Iterate first alond S3r, whenever an entry is re-labeled, check for corresponding entries both on S3r AND on S4 and update these as well.
		//The entries in S4 do not need to be modified except for crossings shared with S3.
		int reverseCrossCount=0;
		int storeCross;
		for(int i=0; i<strand3Length; i++){
			//If modifyCheckS3r[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS3r[i] == 0){
				reverseCrossCount++;
				storeCross=strand3reverse[i];
				if(storeCross>0){
					strand3reverse[i]=reverseCrossCount;
				} else {
					strand3reverse[i]=(-1*reverseCrossCount);
				}
				modifyCheckS3r[i]=1;
				//After relabeling this entry, check for corresponding entries later along S3r.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand3Length; j++){
					if ( (abs(strand3reverse[j]) == abs(storeCross)) && (modifyCheckS3r[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand3reverse[j]=(-1*strand3reverse[i]);
						modifyCheckS3r[j]=1;
					}
				}
				//Next, check for later corresponding entries along S4 with the same conventions as above.
				for(int j=0; j<strand4Length; j++){
					if( (abs(strand4[j]) == abs(storeCross)) && (modifyCheckS4[j] == 0) ){
						strand4[j]=(-1*strand3reverse[i]);
						modifyCheckS4[j]=1;
					}
				}
			}
		}
		//The above for-loop accounts for all crossings along S3r, including those which also show up on S4.
		//The remaining crossings on S4, which do not overlap with S3, should not be modified.
		
		//DEBUG:
		/*
		printf("\n AFTER RE-LABEL");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]");
		printf("\n strand4: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4[i]);
		}
		printf("]\n");
		*/
		
		//Since S3 was reversed but S4 was not, any crossing which was on both S3 and S4 has its orientation reversed. This must be accounted for in both orientedSign arrays.
		//Is is sufficient to use a nested for-loop to search for any crossing indices occuring on both strands, and update the corresponding orientedSign enties.
		for(int i=0; i<strand3Length; i++){
			for(int j=0; j<strand4Length; j++){
				//Each crossing index in the RIGHT tangle occurs twice; if this occurence happens to be on the two different strands, then the orientation is reversed accordingly.
				if( abs(strand3reverse[i]) == abs(strand4[j]) ){
					orientedSignStrand3reverse[i] = (-1*orientedSignStrand3reverse[i]);
					orientedSignStrand4[j] = (-1*orientedSignStrand4[j]);
				}
			}
		}
		
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussSUM[i] = shiftIndexMagnitude(strand3reverse[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussSUM[i]=orientedSignStrand3reverse[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2[i-barsCombined[1]]) <= strand1MaxIndex ){
				gaussSUM[i] = strand2[i-barsCombined[1]];
			} else {
				gaussSUM[i] = shiftIndexMagnitude(strand2[i-barsCombined[1]],strand3MaxIndex);
			}
			orientedSignGaussSUM[i]=orientedSignStrand2[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note crossings already encountered on strand3 are shifted to match.
			if( abs(strand4[i-barsCombined[2]]) <= strand3MaxIndex ){
				gaussSUM[i] = shiftIndexMagnitude(strand4[i-barsCombined[2]],strand1MaxIndex);
			} else {
				gaussSUM[i] = shiftIndexMagnitude(strand4[i-barsCombined[2]],strand2MaxIndex);
			}
			orientedSignGaussSUM[i]=orientedSignStrand4[i-barsCombined[2]];
		}
		
		paritySUM=2;
		barsSUM[0]=barsCombined[2];
		barsSUM[1]=barsCombined[3];
		
	//CASE 7:	
	} else if( (parityLEFT == 2) && (parityRIGHT == 0) ){
		
		//Strand Order: S1-S4-S2-S3
		barsCombined[1]=barsLEFT[0]+(barsRIGHT[1]-barsRIGHT[0]);
		barsCombined[2]=barsLEFT[1]+(barsRIGHT[1]-barsRIGHT[0]);
		barsCombined[3]=barsLEFT[1]+barsRIGHT[1];
		
		//Index Shifts:
		//All crossings on S4 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by RELABELEDstrand4MaxIndex.
		//Previously encountered crossings on S3 by strand1MaxIndex; previously un-encountered crossings on S4 by strand2MaxIndex.
		
		//Orientation Reversals:
		//The strands S3 and S4 are now encountered in reverse order, even though the orientation along these strands is not changed.
		//Effectively, this means that crossing indices must be relabeled on both strands to reflect this, but crossing orientations do not change.
		
		
		//DEBUG:
		/*
		printf("\n\n BEFORE RE-LABEL");
		printf("\n strand3: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3[i]);
		}
		printf("]\n");
		printf(" strand4: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4[i]);
		}
		printf("]\n");
		*/
		
		//Since strands S3 and S4 are now travelled in the opposite order, the crossing labels must be updated to be increasing along S4 and then along S3.
		//Arrays to specfiy if an entry of S3 or S4 has been modified or not; used to avoid undoing the updated labeling.
		int modifyCheckS3[strand3Length];
		int modifyCheckS4[strand4Length];
		//Initialize the arrays with only entires of 0. Each time a label is updated, change entry to 1.
		for(int i=0; i<strand3Length; i++){
			modifyCheckS3[i]=0;
		}
		for(int i=0; i<strand4Length; i++){
			modifyCheckS4[i]=0;
		}
		//NOTE: the order to run these loops is important. In this case, we travel along S4, then S3.
		//The same crossing index could show up either twice on one strand, or once on both strands; both possibilities must be accounted for.
		//Iterate first along S4, whenever an entry is re-labeled, check for corresponding entries both on S4 AND on S3 and update these as well.
		//After this, iterate along S3 and update any remaining crossing labels that are not shared with S4.
		int crossCount=0;
		int storeCross;
		for(int i=0; i<strand4Length; i++){
			//If modifyCheckS4[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS4[i] == 0){
				crossCount++;
				storeCross=strand4[i];
				if(storeCross>0){
					strand4[i]=crossCount;
				} else {
					strand4[i]=(-1*crossCount);
				}
				modifyCheckS4[i]=1;
				//After relabeling this entry, check for corresponding entries later along S4.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand4Length; j++){
					if ( (abs(strand4[j]) == abs(storeCross)) && (modifyCheckS4[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand4[j]=(-1*strand4[i]);
						modifyCheckS4[j]=1;
					}
				}
				//Next, check for later corresponding entries along S3 with the same conventions as above.
				for(int j=0; j<strand3Length; j++){
					if( (abs(strand3[j]) == abs(storeCross)) && (modifyCheckS3[j] == 0) ){
						strand3[j]=(-1*strand4[i]);
						modifyCheckS3[j]=1;
					}
				}
			}
		}
		//The above for-loop accounts for all crossings along S4, including those which also show up on S3.
		//Now we move along S3 and update any remaining crossings, with the same conventions as above.
		for(int i=0; i<strand3Length; i++){
			//If modifyCheckS3r[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS3[i] == 0){
				crossCount++;
				storeCross=strand3[i];
				if(storeCross>0){
					strand3[i]=crossCount;
				} else {
					strand3[i]=(-1*crossCount);
				}
				modifyCheckS3[i]=1;
				for(int j=0; j<strand3Length; j++){
					if ( (abs(strand3[j]) == abs(storeCross)) && (modifyCheckS3[j] == 0) ){
						strand3[j]=(-1*strand3[i]);
						modifyCheckS3[j]=1;
					}
				}
			}
		}
		//At this point, all labels of the crossings along the strands in the RIGHT tangle are updated to travel first along S4 and then along S3.
		//In order to correctly shift the crossing index when building the combined tangle, we also need to record the maximum crossing index encountered on the relabeled S4.
		int RELABELEDstrand4MaxIndex=0;
		for(int i=0; i<strand4Length; i++){
			if( abs(strand4[i]) > RELABELEDstrand4MaxIndex ){
				RELABELEDstrand4MaxIndex = abs(strand4[i]);
			}
		}
		
		//DEBUG:
		/*
		printf("\n AFTER RE-LABEL");
		printf("\n strand3: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3[i]);
		}
		printf("]");
		printf("\n strand4: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4[i]);
		}
		printf("]\n");
		*/
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussSUM[i] = shiftIndexMagnitude(strand4[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussSUM[i]=orientedSignStrand4[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2[i-barsCombined[1]]) <= strand1MaxIndex ){
				gaussSUM[i] = strand2[i-barsCombined[1]];
			} else {
				gaussSUM[i] = shiftIndexMagnitude(strand2[i-barsCombined[1]],RELABELEDstrand4MaxIndex);
			}
			orientedSignGaussSUM[i]=orientedSignStrand2[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note crossings already encountered on strand4 are shifted to match.
			if( abs(strand3[i-barsCombined[2]]) <= RELABELEDstrand4MaxIndex ){
				gaussSUM[i] = shiftIndexMagnitude(strand3[i-barsCombined[2]],strand1MaxIndex);
			} else {
				gaussSUM[i] = shiftIndexMagnitude(strand3[i-barsCombined[2]],strand2MaxIndex);
			}
			orientedSignGaussSUM[i]=orientedSignStrand3[i-barsCombined[2]];
		}
		
		paritySUM=2;
		barsSUM[0]=barsCombined[0];
		barsSUM[1]=barsCombined[3];
		
	//CASE 8:	
	} else if( (parityLEFT == 2) && (parityRIGHT == 1) ){
		
		//Strand Order: S1-S3r-S2r-S4r
		barsCombined[1]=barsLEFT[0]+barsRIGHT[0];
		barsCombined[2]=barsLEFT[1]+barsRIGHT[0];
		barsCombined[3]=barsLEFT[1]+barsRIGHT[1];
		
		//Index Shifts:
		//All crossings on S3 by strand1MaxIndex.
		//Previously un-encountered crossings on S2 by strand3MaxIndex.
		//Previously encountered crossings on S4 by strand1MaxIndex; previously un-encountered crossings on S4 by strand2MaxIndex.
		
		//Orientation Reversals:
		//Strands S2, S3, and S4 are now all traveled in reverse. The order the strands are encountered in is not changed.
		//Any LEFT crossing on both S1 and S2 has its orientation reversed; since both S3 and S4 are reversed, no RIGHT crossing has its orientation reversed.
		
		
		
		//REVERSAL CONSTRUCTION:
		//Initilize the reversals of S2, S3, and S4:
		int strand2reverse[strand2Length];
		int strand3reverse[strand3Length];
		int strand4reverse[strand4Length];
		int orientedSignStrand2reverse[strand2Length];
		int orientedSignStrand3reverse[strand3Length];
		int orientedSignStrand4reverse[strand4Length];
		//List the crossings in the opposite order on each strand:
		for(int i=strand2Length; i>0; i--){
			strand2reverse[i-1]=strand2[strand2Length-i];
			orientedSignStrand2reverse[i-1]=orientedSignStrand2[strand2Length-i];
		}
		for(int i=strand3Length; i>0; i--){
			strand3reverse[i-1]=strand3[strand3Length-i];
			orientedSignStrand3reverse[i-1]=orientedSignStrand3[strand3Length-i];
		}
		for(int i=strand4Length; i>0; i--){
			strand4reverse[i-1]=strand4[strand4Length-i];
			orientedSignStrand4reverse[i-1]=orientedSignStrand4[strand4Length-i];
		}
		
		
		//DEBUG:
		/*
		printf("\n\n BEFORE RE-LABEL");
		printf("\n strand2reverse: [");
		for(int i=0;i<strand2Length;i++){
			printf(" %d ", strand2reverse[i]);
		}
		printf("]");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]");
		printf("\n strand4reverse: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4reverse[i]);
		}
		printf("]\n");
		*/
		
		
		//Since the crossings on each of S2, S3, and S4 are now encountered in the opposite order, the indices in the reverse arrays must be re-labeled accordingly.
		//Arrays to specfiy if an entry of the reverse arrays has been modified or not; used to avoid undoing the updated labeling.
		int modifyCheckS2r[strand2Length];
		int modifyCheckS3r[strand3Length];
		int modifyCheckS4r[strand4Length];
		//Initialize both arrays as only entires of 0. Each time a label is updated, change entry to 1.
		for(int i=0; i<strand2Length; i++){
			modifyCheckS2r[i]=0;
		}
		for(int i=0; i<strand3Length; i++){
			modifyCheckS3r[i]=0;
		}
		for(int i=0; i<strand4Length; i++){
			modifyCheckS4r[i]=0;
		}
		//The crossing indices must be modified so that they are increasing while traveling in reverse order.
		//Iterate through each entry in the strand and relabel the crossings as encountered.
		//Check if any other entries in the strand match this one; update those entries correspondingly.
		//Check against the modifyCheck array to avoid modifying an entry that has already been accounted for.
		
		//First relabel the crossing indices on S2r. Any crossing which is on both S1 and S2 has already been accounted for.
		//Thus, we need only consider those crossing indices which are greater than strand1MaxIndex.
		int S2crossCount=strand1MaxIndex;
		int storeCross;
		for(int i=0; i<strand2Length; i++){
			//We only need to consider those crossings which are BOTH not shared with S1 and which have not already been modified.
			if( (modifyCheckS2r[i] == 0) && ( abs(strand2reverse[i]) > strand1MaxIndex ) ){
				S2crossCount++;
				storeCross=strand2reverse[i];
				if(storeCross>0){
					strand2reverse[i]=S2crossCount;
				} else {
					strand2reverse[i]=(-1*S2crossCount);
				}
				modifyCheckS2r[i]=1;
				//After relabeling this entry, check for corresponding entries later along S2r.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand2Length; j++){
					if ( (abs(strand2reverse[j]) == abs(storeCross)) && (modifyCheckS2r[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand2reverse[j]=(-1*strand2reverse[i]);
						modifyCheckS2r[j]=1;
					}
				}
			}
		}
		//At this point, all crossings indices on S2r not shared with S1 are relabled to be increasing.
		
		//NOTE: the order to run these loops is important. In this case, we travel along S3r, then S4r.
		//The same crossing index could show up either twice on one strand, or once on both strands; both possibilities must be accounted for.
		//Iterate first along S3r, whenever an entry is relabeled, check for corresponding entries both on S3r AND on S4r and update these as well.
		//After finishing S3r, then move to S4r. Any crossing on both strands should already be accounted for; we now need only check along S4r.
		int reverseCrossCount=0;
		for(int i=0; i<strand3Length; i++){
			//If modifyCheckS3r[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS3r[i] == 0){
				reverseCrossCount++;
				storeCross=strand3reverse[i];
				if(storeCross>0){
					strand3reverse[i]=reverseCrossCount;
				} else {
					strand3reverse[i]=(-1*reverseCrossCount);
				}
				modifyCheckS3r[i]=1;
				//After relabeling this entry, check for corresponding entries later along S3r.
				//Check against the modifyCheck arrays to avoid accidentally overwriting entries that have already been updated.
				for(int j=0; j<strand3Length; j++){
					if ( (abs(strand3reverse[j]) == abs(storeCross)) && (modifyCheckS3r[j] == 0) ){
						//If a matching crossing index is found, update it to match the earlier entry.
						//Note that these entries must have the same magnitude but opposite signs (the under/over-strand information), which is why this entry is defined as below.
						strand3reverse[j]=(-1*strand3reverse[i]);
						modifyCheckS3r[j]=1;
					}
				}
				//Next, check for later corresponding entries along S4r with the same conventions as above.
				for(int j=0; j<strand4Length; j++){
					if( (abs(strand4reverse[j]) == abs(storeCross)) && (modifyCheckS4r[j] == 0) ){
						strand4reverse[j]=(-1*strand3reverse[i]);
						modifyCheckS4r[j]=1;
					}
				}
			}
		}
		//The above for-loop accounts for all crossings along S3r, including those which also show up on S4r.
		//Now we move along S4r and update any remaining crossings, with the same conventions as above.
		for(int i=0; i<strand4Length; i++){
			//If modifyCheckS4r[i] is nonzero, this entry has already been updated and we may skip it.
			if(modifyCheckS4r[i] == 0){
				reverseCrossCount++;
				storeCross=strand4reverse[i];
				if(storeCross>0){
					strand4reverse[i]=reverseCrossCount;
				} else {
					strand4reverse[i]=(-1*reverseCrossCount);
				}
				modifyCheckS4r[i]=1;
				for(int j=0; j<strand4Length; j++){
					if ( (abs(strand4reverse[j]) == abs(storeCross)) && (modifyCheckS4r[j] == 0) ){
						strand4reverse[j]=(-1*strand4reverse[i]);
						modifyCheckS4r[j]=1;
					}
				}
			}
		}
		//At this point, all labels of the strands in the RIGHT tangle are updated to travel in the reverse direction, with increasing crossing index.
		
		//DEBUG:
		/*
		printf("\n AFTER RE-LABEL");
		printf("\n strand2reverse: [");
		for(int i=0;i<strand2Length;i++){
			printf(" %d ", strand2reverse[i]);
		}
		printf("]");
		printf("\n strand3reverse: [");
		for(int i=0;i<strand3Length;i++){
			printf(" %d ", strand3reverse[i]);
		}
		printf("]");
		printf("\n strand4reverse: [");
		for(int i=0;i<strand4Length;i++){
			printf(" %d ", strand4reverse[i]);
		}
		printf("]\n");
		*/
		
		
		//Since S2 was reversed but S1 was not, any crossing which was on both S1 and S2 has its orientation reversed. This must be accounted for in both orientedSign arrays.
		//Is is sufficient to use a nested for-loop to search for any crossing indices occuring on both strands, and update the corresponding orientedSign entries.
		for(int i=0; i<strand1Length; i++){
			for(int j=0; j<strand2Length; j++){
				//Each crossing index in the LEFT tangle occurs twice; if this occurence happens to be on the two different strands, then the orientation is reversed accordingly.
				if( abs(strand1[i]) == abs(strand2reverse[j]) ){
					orientedSignStrand1[i] = (-1*orientedSignStrand1[i]);
					orientedSignStrand2reverse[j] = (-1*orientedSignStrand2reverse[j]);
				}
			}
		}
		//If the orientation has been changed for any crossings on S1, we need to update the orientedSignGaussSUM[] array since this was earlier set equal to orientedSignStrand1.
		for(int i=0; i<strand1Length; i++){
			orientedSignGaussSUM[i]=orientedSignStrand1[i];
		}
		
		
		for(int i=barsCombined[0]; i < barsCombined[1]; i++){
			gaussSUM[i] = shiftIndexMagnitude(strand3reverse[i-barsCombined[0]],strand1MaxIndex);
			orientedSignGaussSUM[i]=orientedSignStrand3reverse[i-barsCombined[0]];
		}
		for(int i=barsCombined[1]; i < barsCombined[2]; i++){
			//Note that we don't want to shift the crossing index of any crossings that were already encountered on strand1.
			if( abs(strand2reverse[i-barsCombined[1]]) <= strand1MaxIndex ){
				gaussSUM[i] = strand2reverse[i-barsCombined[1]];
			} else {
				gaussSUM[i] = shiftIndexMagnitude(strand2reverse[i-barsCombined[1]],strand3MaxIndex);
			}
			orientedSignGaussSUM[i]=orientedSignStrand2reverse[i-barsCombined[1]];
		}
		for(int i=barsCombined[2]; i < barsCombined[3]; i++){
			//Note that crossings already encountered on S3r are shifted to match.
			if( abs(strand4reverse[i-barsCombined[2]]) <= strand3MaxIndex ){
				gaussSUM[i] = shiftIndexMagnitude(strand4reverse[i-barsCombined[2]],strand1MaxIndex);
			} else {
				gaussSUM[i] = shiftIndexMagnitude(strand4reverse[i-barsCombined[2]],strand2MaxIndex);
			}
			orientedSignGaussSUM[i]=orientedSignStrand4reverse[i-barsCombined[2]];
		}
		
		paritySUM=2;
		barsSUM[0]=barsCombined[0];
		barsSUM[1]=barsCombined[3];
		
	//CASE 9:	
	} else if( (parityLEFT == 2) && (parityRIGHT == 2) ){
		//If two parity infinity tangles are added, the resulting structure is not a 2-string tangle since it has some sort of closed-loop thing in the middle.
		//I don't know whether or not this structure is of interest, but this possibility will not be considered allowed for the purpose of this function.
		printf("\n ERROR: The sum of two parity infinity tangles is no longer a 2-string tangle.\n");
		return;
		
	//ERROR CASE:	
	} else {
		//This is only possible if either of the two input tangles has parity other than 0, 1, or 2 (infinity). This is not allowed, so the function stops here.
		printf("\n ERROR: It appears that one (or both) of the input tangles has parity other than 0, 1, or infinity.\n");
		return;
	}
	
	
	//DEBUG:
	/*
	printf("\n barsCombined: [");
	for(int i=0;i<4;i++){
		printf(" %d ", barsCombined[i]);
	}
	printf("]\n");
	printf("\n gaussSUM: [");
	for(int i=0; i<(2*numOfCrossingsSUM); i++){
		printf(" %d ", gaussSUM[i]);
	}
	printf("]");
	printf("\n orientedSignGaussSUM: [");
	for(int i=0; i<(2*numOfCrossingsSUM); i++){
		printf(" %d ", orientedSignGaussSUM[i]);
	}
	printf("]");
	printf("\n barsSUM: [ %d %d ]", barsSUM[0],barsSUM[1]);
	printf("\n paritySUM: %d \n", paritySUM);
	*/
	
	
}




/*
(NC, 7/30/18)
This is a function designed to print the Gauss code of a 2-string tangle in the standard form.
It is intended to work with 2-string tangles of parity 0, 1, or 2 (infinity).
Note that the Guass code doesn't have to be ordered for this function to work, but that is generally a helpful convention.

Function Inputs:
numOfCrossingsInput: the number of crossings of the tangle being considered.
gaussInput[]: an array of length 2*numOfCrossingsInput storing indices of the tangle's crossing. A positive index indicates an overcrossing; a negative index indicates an undercrossing.
barsInput[]: an array of length 2 where each entry indicates the last crossing encountered (regardless of index) before the end of the first or second strand, respectively.
orientedSignGaussInput[]: an array of length 2*numOfCrossingsInput that stores the orientation of each crossing: 1 indicates a positive orientation; -1 indicates a negative orientation.
parityInput: the parity of the input tangle. This is allowed to be 0, 1, or 2 (infinity).

Functoin Outputs:
This functions prints out the ordered Gauss code of the input tangle.
It can also output the number of crossings and the parity, if desired, but this functionality is currently commented out.


This function calls:
	nothing
This function is called by:
	findRationalTangles()
	
*/

void printTangleGaussCode(int numOfCrossingsInput, int *gaussInput, int *barsInput, int *orientedSignGaussInput, int parityInput){
	
	//In each case, the Gauss code constructed is the Gauss code of the oriented link obtained by taking either the numerator or denominator closure of the tangle.
	//It lists the outer circle first, oriented clockwise and beginning with the crossing at NW endpoint, and then travels along the strand starting from this same NW endpoint crossing.
	
	printf("\n Gauss Code: ");
	
	/*
	//Toggle this if the parity and crossing number are desired.
	if( (parityInput == 0) || (parityInput == 1)){
		printf("parity %d tangle with %d crossings", parityInput,numOfCrossingsInput);
	} else if( parityInput == 2 ){
		printf("parity infinity tangle with %d crossings",numOfCrossingsInput);
	} else {
		printf("unknown parity tangle with %d crossings",numOfCrossingsInput);
	}
	*/
	
	if( parityInput == 0 ){
		
		//Use the denominator closure for a parity 0 tangle.
		printf("\n a1-a2+a3-a4+|b1-)");
		
		if ( numOfCrossingsInput == 0 ){
			printf("(b2+b3-)");
		}
	
		for (int i = 0; i < (2*numOfCrossingsInput); i++){
			
			if (i == barsInput[0]){
				printf("(b2+b3-)");
			}
			
			if (gaussInput[i] > 0){
				printf("a");
			}else{
				printf("b");
			}

			printf("%d", abs(gaussInput[i])+4);
		
			if (orientedSignGaussInput[i] < 0){
				printf("-");
			}else{
				printf("+");
			}
		}
		
		printf("(b4+");
		
	} else if( parityInput == 1 ){
		
		//Use the denominator closure for a parity 1 tangle.
		printf("\n a1-a2-a3+a4+|b1-)");
		
		if ( numOfCrossingsInput == 0 ){
			printf("(b3+b2-)");
		}
	
		for (int i = 0; i < (2*numOfCrossingsInput); i++){
			
			if (i == barsInput[0]){
				printf("(b3+b2-)");
			}
			
			if (gaussInput[i] > 0){
				printf("a");
			}else{
				printf("b");
			}

			printf("%d", abs(gaussInput[i])+4);
		
			if (orientedSignGaussInput[i] < 0){
				printf("-");
			}else{
				printf("+");
			}
		}
		
		printf("(b4+");
		
	} else if( parityInput == 2 ){
		
		//Use the numerator closure for a parity infinity tangle.
		printf("\n a1-a2+a3-a4+|b1-)");
		
		if ( numOfCrossingsInput == 0){
			printf("(b4+b3-)");
		}
	
		for (int i = 0; i < (2*numOfCrossingsInput); i++){
			
			if (i == barsInput[0]){
				printf("(b4+b3-)");
			}
			
			if (gaussInput[i] > 0){
				printf("a");
			}else{
				printf("b");
			}

			printf("%d", abs(gaussInput[i])+4);
		
			if (orientedSignGaussInput[i] < 0){
				printf("-");
			}else{
				printf("+");
			}
		}
		
		printf("(b2+");
		
	} else{
		printf("ERROR - input parity is something other than 0, 1, or 2 (infinity)");
	}
	
	printf("\n");
	
}





/*
(NC, 7/30/18)
This function is designed to construct a rational tangle from a twist vector.
In principle, any rational tangle can be built from the [0] or [infinity] tangle via a sequence of horizontal and vertical twists.
Adding a twist to a tangle is effectively the same thing as adding or multiplying by a parity 1 tangle with one crossing.
There are two such parity 1 tangles, denoted [1] and [-1], which are differentiated by the "slope" of the overstrand.
The distinction is as follows: if the first strand is the over-strand, this is the [1] tangle; if the first strand is the under-strand, this is the [-1].
These two tangles form the building blocks of all rational tangles. Hence, I will refer to them as the "atomic" tangles, mostly because that sounds really cool.

It is worth noting that, while the primary purpose of this function is to construct rational tangles, it can also be used with other tangle families.
This function merely constructs the tangle defined by a twist vector, but this new tangle can be added to tangles other than the [0] tangle.
This could be useful if one wanted to show that two non-rational tangles are in the same family; one could show that adding some twists to one tangle results in the other.

Function Inputs:
twistVectLength: the length of the twistVector[] array.
twistVector: an array of length twistVectLength storing information about the twists to be added.

This function will also take the global arrays and variables storing information about the SUM and PROD tangles as input. This amounts to five inputs each.
While this isn't elegant, these will be used and overwritten a number of times as local storage arrays while constructing the twist tangle.

The twistVector may be regarded as a tuple, (x1,x2,...,xn). Each entry in the vector refers to an integer number of horizontal or vertical twists.
By convention, the last entry xn is always assumed to refer to a horizontal twist. The last entry (but no other) is allowed to be 0 to satisfy this convention.
Otherwise, the entries in the vector alternate between horizontal and vertical twists. The sign of the entry refers to the "slope" of the twists (either the [1] or [-1] tangle).
Note that this convention means the first entry in a vector of odd length is a horizontal twist; the first entry in a vector of even length is a vertical twist.

Function Outputs (in the form of modified global variables):
twistCrossings: the number of crossings of the constructed tangle.
gaussTwist[]: an array storing the (ordered) Gauss code of the constructed tangle. Each entry contains the index of the crossing encountered. Positive entries denote an over-crossing; negative entries denote an under-crossing.
barsTwist[]: an array with 2 entries to denote where each of the two strings exits the constructed tangle. The first/second entry denotes the last crossing encountered along the first/second string (regardless of crossing index).
orientedSignGaussTwist[]: an array of the same length as gaussTwist[] containing information on the orientation of each crossing encountered (either 1 or -1). Each entry matches the corresponding entry in gaussTwist[].
parityTwist: the parity of the constructed tangle, either 0, 1, or 2 (infinity).


This function calls:
	tangleSum
	tangleProduct
This function is called by:
	findRationalTangles()

*/

void buildTwistTangle(int twistVectLength, int *twistVector, int &twistCrossings, int *gaussTwist, int *barsTwist, int *orientedSignGaussTwist, int &parityTwist, int &numOfCrossingsSUM, int *gaussSUM, int *barsSUM, int *orientedSignGaussSUM, int &paritySUM, int &numOfCrossingsPROD, int *gaussPROD, int *barsPROD, int *orientedSignGaussPROD, int &parityPROD){
	
	//Initialize the Gauss code for the two atomic tangles:
	//Positive Slope: [1]
	int atomicPositiveGauss[2] = {1,-1};
	int atomicPositiveBars[2] = {1,2};
	int atomicPositiveOrientedSign[2] = {-1,-1};
	//Negative Slope: [-1]
	int atomicNegativeGauss[2] = {-1,1};
	int atomicNegativeBars[2] = {1,2};
	int atomicNegativeOrientedSign[2] = {1,1};
	//Note that both atomic tangles have parity 1 and only a single crossing, so there is no special reason to define variables for the crossing and parity of these.
	
	//The twist tangle to be constructed will have a number of crossings equal to the sum of the magnitudes of the twistVector entries. Compute this here.
	twistCrossings=0;
	for(int i=0; i<twistVectLength; i++){
		twistCrossings = (twistCrossings + abs(twistVector[i]));
	}
	
	//By convention, the entries in the twist vector alternate detween denoting horizontal and vertical twists, with the last entry always denoting a horizontal twist.
	//As was observed in the preamble, the first twist is horizontal if the twistVector has odd length, and vertical if the twistVector has even length.
	//Initialze a variable to denote the twist to perform: positive denotes horizontal, negative denotes vertical.
	int nextTwist;
	if( (twistVectLength%2) == 1 ){
		//horizontal twist:
		nextTwist=1;
	} else if( (twistVectLength%2) == 0 ){
		//vertical twist:
		nextTwist=-1;
	}
	
	//Initialize local arrays for building the twist tangle. At the end of the function, the output global variables will be set equal to these.
	//If the first twist is horizontal, the default tangle will be set to the [0] tangle.
	//If the first twist is vertical, the default tangle will be set to the [infinity] tangle.
	//Note that neither of these has any crossings, we only need this so that there is a starting tangle to perform twists on.
	//Since there are no crossings, the gauss[] and orientedSign[] arrays are empty, but we still specify crossings, bars[], and parity.
	int twistCrossingsTemp=0;
	int gaussTwistTemp[2*twistCrossings];
	int barsTwistTemp[2] = {0,0};
	int orientedSignGaussTwistTemp[2*twistCrossings];
	int parityTwistTemp;
	if(nextTwist>0){
		parityTwistTemp=0;
	} else if(nextTwist<0){
		parityTwistTemp=2;
	}
	
	
	for(int i=0; i<twistVectLength; i++){
		
		if(nextTwist>0){
			//If nextTwist>0, we use horizontal twists (SUM function).
			for(int j=0; j<abs(twistVector[i]); j++){
				if(twistVector[i]>0){
					//If twistVector[j]>0, we use [1] (positive atomic tangle)
					tangleSum(twistCrossingsTemp,gaussTwistTemp,barsTwistTemp,orientedSignGaussTwistTemp,parityTwistTemp,1,atomicPositiveGauss,atomicPositiveBars,atomicPositiveOrientedSign,1,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
					
				} else if(twistVector[i]<0){
					//If twistVecotr[j]<0, we use [-1] (negative atomic tangle)
					tangleSum(twistCrossingsTemp,gaussTwistTemp,barsTwistTemp,orientedSignGaussTwistTemp,parityTwistTemp,1,atomicNegativeGauss,atomicNegativeBars,atomicNegativeOrientedSign,1,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
					
				} else {
					//The only other possibility is if twistVector[j]=0, in which case there is no twist so nothing happens.
					//With the allowed convention, a twistVector entry of 0 is only allowed if it's the last entry.
				}
				
				//After performing the twist, set the temporary arrays equal to the SUM arrays, then repeat as needed.
				twistCrossingsTemp=numOfCrossingsSUM;
				for(int k=0; k<(2*numOfCrossingsSUM); k++){
					gaussTwistTemp[k]=gaussSUM[k];
					orientedSignGaussTwistTemp[k]=orientedSignGaussSUM[k];
				}
				barsTwistTemp[0]=barsSUM[0];
				barsTwistTemp[1]=barsSUM[1];
				parityTwistTemp=paritySUM;
			}
				
		} else if(nextTwist<0){
			//If nextTwist<0, we use vertical twists (PROD function).
			for(int j=0; j<abs(twistVector[i]); j++){
				if(twistVector[i]>0){
					//If twistVector[j]>0, we use [1] (positive atomic tangle)
					tangleProduct(twistCrossingsTemp,gaussTwistTemp,barsTwistTemp,orientedSignGaussTwistTemp,parityTwistTemp,1,atomicPositiveGauss,atomicPositiveBars,atomicPositiveOrientedSign,1,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
					
				} else if(twistVector[i]<0){
					//If twistVecotr[j]<0, we use [-1] (negative atomic tangle)
					tangleProduct(twistCrossingsTemp,gaussTwistTemp,barsTwistTemp,orientedSignGaussTwistTemp,parityTwistTemp,1,atomicNegativeGauss,atomicNegativeBars,atomicNegativeOrientedSign,1,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
					
				} else {
					//The only other possibility is if twistVector[j]=0, in which case there is no twist so nothing happens.
					//With the allowed convention, a twistVector entry of 0 is only allowed if it's the last entry.
				}
				
				//After performing the twist, set the temporary arrays equal to the SUM arrays, then repeat as needed.
				twistCrossingsTemp=numOfCrossingsPROD;
				for(int k=0; k<(2*numOfCrossingsPROD); k++){
					gaussTwistTemp[k]=gaussPROD[k];
					orientedSignGaussTwistTemp[k]=orientedSignGaussPROD[k];
				}
				barsTwistTemp[0]=barsPROD[0];
				barsTwistTemp[1]=barsPROD[1];
				parityTwistTemp=parityPROD;
			}
		}
		
		nextTwist = (-1*nextTwist);
	}
	
	//After the conclusion of these nested loops, the constructed twist tangle is stored in the temporary arrays. Update the corresponding global variables with these values.
	for(int i=0; i<(2*twistCrossings); i++){
		gaussTwist[i]=gaussTwistTemp[i];
		orientedSignGaussTwist[i]=orientedSignGaussTwistTemp[i];
	}
	barsTwist[0]=barsTwistTemp[0];
	barsTwist[1]=barsTwistTemp[1];
	parityTwist=parityTwistTemp;
	
	
	//DEBUG
	/*
	printf("\n\n TWIST TANGLE:\n");
	printf("\n twistVector: [");
	for(int i=0; i<twistVectLength; i++){
		printf(" %d ", twistVector[i]);
	}
	printf("]");
	printTangleGaussCode(twistCrossings,gaussTwist,barsTwist,orientedSignGaussTwist,parityTwist);
	*/
	
	
	
}



/*
(NC, 8/1/18)
This is a very simple function to compute a^b. This is only inteded to work with non-negative integers.

Function Inputs:
a: a non-negative integer.
b: another non-negative integer.

Function Outputs:
This function returns the integer a^b.


This function calls:
	none
This function is called by:
	findRationalTangles()
	miscellaneous functions?

*/

int ApowerB(int a, int b){
	int count=1;
	for(int i=0; i<b; i++){
		count = a*count;
	}
	return count;
}




/*
(NC, 7/31/18)
The total number of compositions of an integer n is 2^(n-1). This is a short function to calculate this value with a for loop and then return the number of compositions.

Function Inputs:
n: an integer.

Function Outputs:
2^n is returned as an integer.


This function calls:
	none
This function is called by:
	twistVectorSigns()
	findRationalTangles()
	
*/

int countCompositions(int n){
	int count=1;
	for(int i=0; i<(n-1); i++){
		count = 2*count;
	}
	return count;
}


/*
(NC, 7/31/18)
For a given integer n, this function finds all 2^(n-1) distinct compositions of n (recall that a composition of an integer is an ordered partition of that integer).
It does this recursively by filling in the entries of a matrix. Each row of the matrix is a distinct composition. The number of nonzero entries in a row varies.

Function Inputs:
n: the integer being considered.
nFixed: a storage index used to specify the correct first index of the array of matrices while updateding recursively. This should be n when calling this function elsewhere.
rowIndex: the default row index. This is needed for the function to update the matrix recursively. This should be 0 when calling this function elsewhere.
colIndex: the default column index. This is also needed for the function to update the matrix recursively. This should be 0 when calling this function elsewhere.

Function Ouputs (in the form of modified global variables):
compositionList[][][]: a triple-indexed array (really an array of matrices) to store the compositions. The first index is the integer being considered, in this case n.
The second two indices define a matrix, each row of which is a distinct composition of n. The dimensions will need to be changed accordingly for n>9.


This function calls:
	itself (recursive)
This function is called by:
	findPositiveTwistVectors()
	findRationalTangles()

*/

void findCompositions(int n, int nFixed, int compositionList[MAXN+1][CLB][MAXN+1], int rowIndex, int colIndex){
	
	int localRowIndex=rowIndex;
	int localColIndex=colIndex;
	
	for(int k=0; k<n; k++){
		
		int j;
		j=countCompositions(k);
		
		for(int i=0; i<j; i++){
			compositionList[nFixed][localRowIndex+i][localColIndex] = (n-k);
		}
		
		findCompositions(k,nFixed,compositionList,localRowIndex,localColIndex+1);
		
		localRowIndex = (localRowIndex+j);
		
	}
	
	
}


/*
(NC, 7/31/18)
A function to build up all positive rational twist vectors with up to maxN crossings. It uses the compositions generated by the findCompositions() function to do this.
Later functions consider how to add signs to these twist vectors. This is where the compositionList[][][] and rationalTwistVectLength[][] arrays are actually populated.

Function Inputs:
maxN: the maximum integer whose compositions are considered.

Function Outputs 9in the form of modified global variables):
compositionList[][][]: an array of matrices, where the rows of any fixed matrix store individual compositions.
rationalTwistVectLength[][]: an array storing the length of each individual compositions.


This function calls:
	findCompositions()
	countCompositions()
This function is called by:
	findRationalTangles()

*/

void findPositiveTwistVectors(int compositionList[MAXN+1][CLB][MAXN+1], int rationalTwistVectLength[MAXN+1][CLB], int maxN){
	
	int countNonzero;
	int localCompositionCount;
	
	for(int n=0; n<=maxN; n++){
		//Populate the three dimensional arrays used to store the integer compositions.
		findCompositions(n,n,compositionList,0,0);
		localCompositionCount=countCompositions(n);
		
		//For a given composition of n, count the number of nonzero entries and store this as the length of the corresponding rational twist vector.
		//The above description is subject to the convention that rational twist vectors have odd length; if the length is even, increase it by 1.
		for(int i=0; i<localCompositionCount; i++){
			countNonzero=0;
			for(int j=0; j<n; j++){
				if( compositionList[n][i][j] != 0 ){
					countNonzero++;
				}
			}
			if( (countNonzero%2) == 1 ){
				//Odd length vector:
				rationalTwistVectLength[n][i] = countNonzero;
			} else if( (countNonzero%2) == 0 ){
				//Even length vector; increment by 1:
				rationalTwistVectLength[n][i] = countNonzero+1;
			}
		}
		
		//DEBUG:
		/*
		printf("\n compositionList: n = %d", n);
		for(int i=0; i<localCompositionCount; i++){
			//calculateExtendedFractionPQ(isRational,rationalTwistVectLength[n][i],compositionList[n][i],rationalPQ);
			printf("\n [");
			for(int j=0; j<rationalTwistVectLength[n][i]; j++){
				printf(" %d ", compositionList[n][i][j]);
			}
			printf("]");
		}
		*/
	}	
}


/*
(NC, 8/1/18)
This recursive function is desgined to find all possible ways of changing the signs on the entries of a vector.
In general, for a vector of length m, there are 2^m possible ways to give signs to the entries.
This function stores all of the possibilities in a matrix, where each row corresponds to a distinct signed vector.

Function Inputs:
vectorLength: the length of the vector being considered.
vector[]: an array of length vectorLength; the vector we are considering.
startEntry: used to specify the index of the vector entry to modify; this is used by the recursive part of the function. When calling the function elsewhere, this should be 0.
signVectorCount: a global or static variable used to index the signed vectors obtained in the program, updated recursively. When calling the function elsewhere, this should be 0.

Function Outputs (in the form of modified global variables):
signVectorMatrix: a double indexed array to store the signed vectors. Each row is a different signed vector.
By default, the dimensions are set to be [2^9][9], so this will work with vectors of length at most 9.
If one wants to use this function with larger vectors, increase these dimensions accordingly.


This function calls:
	itself (recursive)
This function is called by:
	twistVectorSigns()

*/

void vectorSignChanges(int vectorLength, int *vector, int startEntry, int signVectorMatrix[2*CLB][MAXN+1], int &signVectorCount){
	
	//Set the current row of the vectorMatrix to be the vector plugged into this function.
	for(int i=startEntry; i<vectorLength; i++){
		signVectorMatrix[signVectorCount][i] = vector[i];
	}
	//Increment signVectorCount, which is the index of the current signVectorMatrix row.
	//As a global variable, this is updated each time the function is called. It is used to index the different possible sign changes.
	signVectorCount++;
	
	for(int j=startEntry; j<vectorLength; j++){
		
		//Initialize the next row of the vectorMatrix as the vector plugged into this function; then change the sign of a single entry.
		for(int i=0; i<vectorLength; i++){
			signVectorMatrix[signVectorCount][i] = vector[i];
		}
		signVectorMatrix[signVectorCount][j] = (-1*signVectorMatrix[signVectorCount][j]);
		
		//Run the function again, shifting the start term to be j+1. The "vector" used is now signVectorMatrix[signVectorCount].
		vectorSignChanges(vectorLength,signVectorMatrix[signVectorCount],j+1,signVectorMatrix,signVectorCount);
		
	}
	
}



/*
(NC, 8/1/18)
This function takes a potential twist vector and finds all possible combinations of signs on the terms of this vector.
In general, for a vector of length m, there are 2^m possible ways to sign the entries of this vector. This function finds all of these and stores them in global array.
Note that every time this function is called, this global array is overwritten, so this must be accounted for when calling this function a number of times in a loop.

(NC, 8/3/18)
ALL TWIST VECTORS WITH BOTH POSITIVE AND NEGATIVE TERMS GENERATE REDUNDANT RATIONAL TANGLES.
While these technically define rational tangles, such tangles can be reduced to a simpler rational tangle via Reidemeister moves.
Hence, we need only consider twist vectors of all positive terms or all negative terms for the purpose of this program.
This functionality might have applications elsewhere, so I won't delete it, but it is no longer called when constructing all rational tangles.


Function Inputs:
n: in the context of this function, this is the integer whose compositions we consider; this is also the number of crossings of the corresponding rational tangle.
twistVectorIndex: the compositions are stored in an array of matrices; for fixed matrix, this index denotes the row of the matrix (where each row is one possible compositions).
compositionList: this is the array of matrices storing the compositions; in the context of this function, the first index is fixed as n and the second index is fixed as twistVectorIndex.
twistVectLength: this is the length of the composition stored in the compositionList[n][twistVectorIndex].

Function Outputs (in the form of modified global variables):
twistSigns[][]: this a 2-dimensional array used to store the possible signed twist vectors; it is needed for input into the vectorSignChanges() function.


This function calls:
	vectorSignChanges()
This function is called by:
	none

*/

void twistVectorSigns(int n, int twistVectorIndex, int compositionList[MAXN+1][CLB][MAXN+1], int twistVectLength, int twistSigns[2*CLB][MAXN+1]){
	
	//Recall that, by convention, rational twist vectors have odd length, and sometimes this means the last entry will be zero.
	//Since zero cannot have its sign changed, we reduce the number of entries checked by 1 if there is a zero entry.
	//We need only check that last entry since this is the only entry allowed to be zero.
	int numNonzeroTerms;
	if( compositionList[n][twistVectorIndex][twistVectLength-1] == 0 ){
		numNonzeroTerms = twistVectLength-1;
	} else {
		numNonzeroTerms=twistVectLength;
	}
	
	
	//DEBUG
	/*
	printf("\n LENGTH: %d \n NON-ZERO: %d \n", twistVectLength, numNonzeroTerms);
	printf(" ORIGINAL TWIST VECTOR: [");
	for(int i=0; i<twistVectLength; i++){
		printf(" %d ", compositionList[n][twistVectorIndex][i]);
	}
	printf("]\n");
	*/
	
	
	//Construct the signed twist vectors and store them in the twistSigns[][] matrix.
	int signedTwistVectorCount=0;
	vectorSignChanges(numNonzeroTerms,compositionList[n][twistVectorIndex],0,twistSigns,signedTwistVectorCount);
	
	
	int numberOfSignedTwistVectors;
	numberOfSignedTwistVectors=ApowerB(2,numNonzeroTerms);
	
	
	//Note that, if the input twistVector had a last entry of 0, we need to go back and add this zero entry to all of the signed twist vectors.
	if( numNonzeroTerms<twistVectLength ){
		for(int i=0; i<numberOfSignedTwistVectors; i++){
			twistSigns[i][twistVectLength-1]=0;
		}
	}
	
	
	//DEBUG
	/*
	printf("\n Composition %d (of n = %d):\n Total of %d  of potential signed twist vectors.", twistVectorIndex, n, numberOfSignedTwistVectors);
	for(int i=0; i<numberOfSignedTwistVectors; i++){
		printf("\n Twist Vector %d: [",i);
		for(int j=0; j<twistVectLength; j++){
			printf(" %d ", twistSigns[i][j]);
		}
		printf("]");
	}
	*/
		
}



/*
(NC, 8/3/18)
This function takes a positive twist vector and generates the corresponding negative twist vector.
Both of these are stored as the first and second rows, respectively, of the two row matrix sameSignTwistVector[][].

Function Inputs:
n: in the context of this function, this is the integer whose compositions we consider; this is also the number of crossings of the corresponding rational tangle.
twistVectorIndex: the compositions are stored in an array of matrices; for fixed matrix, this index denotes the row of the matrix (where each row is one possible compositions).
compositionList: this is the array of matrices storing the compositions; in the context of this function, the first index is fixed as n and the second index is fixed as twistVectorIndex.
twistVectLength: this is the length of the composition stored in the compositionList[n][twistVectorIndex].

Function Outputs (in the form of modified global variables):
sameSignTwistVector[][]: this a 2-dimensional array used to store the positive and negative twist vectors.


This function calls:
	none
This function is called by UPDATE
	findRationalTangles()

*/

void twistVectorPositiveNegative(int n, int twistVectorIndex, int compositionList[MAXN+1][CLB][MAXN+1], int twistVectLength, int sameSignTwistVector[2][MAXN+1]){
	
	//Store the input (positive) twist vector sameSignTwistVector[0].
	//Create a negative twist vector by swapping the signs and store this as sameSignTwistVector[1].
	for(int i=0; i<twistVectLength; i++){
		sameSignTwistVector[0][i] = compositionList[n][twistVectorIndex][i];
		sameSignTwistVector[1][i] = (-1*compositionList[n][twistVectorIndex][i]);
	}
	//Store the negative twist vector as sameSignTwistVector
	
}





/*
(NC, 12/18/18)
Description pending


This function calls:
	nothing
	
This function is called by:
	findGeneralizedMontesinosTangles()


*/

void writeGMTangleToFileSQL(int indexA, int indexB, int indexC, int indexD, int coreCrossings, int crossingsInput, int *gaussInput, int *barsInput, int *orientedSignGaussInput, int parityInput, int numOfComponents, int compRationalPQ[][2], int *inputTwistVector, int inputTwistVectorLength){
	
	int endpointTwists=(crossingsInput-coreCrossings);
	
	FILE *tangleFile3;
	tangleFile3 = fopen("MontesinosTangleGeneralizedListSQL.txt", "a");
	
	//gm_tangle_number, program_index, crossing_number
	fprintf(tangleFile3, "INSERT INTO generalized_montesinos_tangle_list(gm_tangle_number,program_index,crossing_number,tangle_parity,core_subtangle,twist_vector,gauss_code,core_crossing_number,endpoint_twist_number) VALUES(%d,'%d.%d.%d.%d',%d,", totalGMTanglesFound, indexA, indexB, indexC, indexD, crossingsInput);
	
	//tangle_parity
	if( (parityInput==0) || (parityInput==1) ){
		fprintf(tangleFile3, "'%d',", parityInput);
	} else if( parityInput==2 ) {
		fprintf(tangleFile3, "'infinity',");
	} else {
		fprintf(tangleFile3, "'unknown',");
	}
	
	//core_subtangle
	fprintf(tangleFile3, "'(%d/%d)", compRationalPQ[0][0], compRationalPQ[0][1]);
	for(int i=1; i<numOfComponents; i++){
		fprintf(tangleFile3, "+(%d/%d)", compRationalPQ[i][0], compRationalPQ[i][1]);
	}
	fprintf(tangleFile3, "',");
	
	//twist_vector
	fprintf(tangleFile3, "'[");
	if( inputTwistVectorLength == 0 ){
		fprintf(tangleFile3, " 0");
	} else {
		for(int i=0; i<inputTwistVectorLength; i++){
			fprintf(tangleFile3, " %d", inputTwistVector[i]);
		}
	}
	fprintf(tangleFile3, " ]',");
	
	//gauss_code
	if( parityInput == 0 ){
		
		//Use the denominator closure for a parity 0 tangle.
		fprintf(tangleFile3, "'a1-a2+a3-a4+|b1-)");
		
		if ( crossingsInput == 0 ){
			fprintf(tangleFile3, "(b2+b3-)");
		}
	
		for (int i = 0; i < (2*crossingsInput); i++){
			
			if (i == barsInput[0]){
				fprintf(tangleFile3, "(b2+b3-)");
			}
			
			if (gaussInput[i] > 0){
				fprintf(tangleFile3, "a");
			}else{
				fprintf(tangleFile3, "b");
			}

			fprintf(tangleFile3, "%d", abs(gaussInput[i])+4);
		
			if (orientedSignGaussInput[i] < 0){
				fprintf(tangleFile3, "-");
			}else{
				fprintf(tangleFile3, "+");
			}
		}
		
		fprintf(tangleFile3, "(b4+',");
		
	} else if( parityInput == 1 ){
		
		//Use the denominator closure for a parity 1 tangle.
		fprintf(tangleFile3, "'a1-a2-a3+a4+|b1-)");
		
		if ( crossingsInput == 0 ){
			fprintf(tangleFile3, "(b3+b2-)");
		}
	
		for (int i = 0; i < (2*crossingsInput); i++){
			
			if (i == barsInput[0]){
				fprintf(tangleFile3, "(b3+b2-)");
			}
			
			if (gaussInput[i] > 0){
				fprintf(tangleFile3, "a");
			}else{
				fprintf(tangleFile3, "b");
			}

			fprintf(tangleFile3, "%d", abs(gaussInput[i])+4);
		
			if (orientedSignGaussInput[i] < 0){
				fprintf(tangleFile3, "-");
			}else{
				fprintf(tangleFile3, "+");
			}
		}
		
		fprintf(tangleFile3, "(b4+',");
		
	} else if( parityInput == 2 ){
		
		//Use the numerator closure for a parity infinity tangle.
		fprintf(tangleFile3, "'a1-a2+a3-a4+|b1-)");
		
		if ( crossingsInput == 0){
			fprintf(tangleFile3, "(b4+b3-)");
		}
	
		for (int i = 0; i < (2*crossingsInput); i++){
			
			if (i == barsInput[0]){
				fprintf(tangleFile3, "(b4+b3-)");
			}
			
			if (gaussInput[i] > 0){
				fprintf(tangleFile3, "a");
			}else{
				fprintf(tangleFile3, "b");
			}

			fprintf(tangleFile3, "%d", abs(gaussInput[i])+4);
		
			if (orientedSignGaussInput[i] < 0){
				fprintf(tangleFile3, "-");
			}else{
				fprintf(tangleFile3, "+");
			}
		}
		
		fprintf(tangleFile3, "(b2+',");
		
	} else{
		fprintf(tangleFile3, "NULL,");
	}
	
	//core_crossing_number, endpoint_twists
	fprintf(tangleFile3, "%d,%d);\n", coreCrossings, endpointTwists);
	
	fclose(tangleFile3);
	
}




/*
(NC, 12/18/18)
This function records the rational number (p/q) and the order of each of the rational components used in the construction of the current Montesinos tangle in the form a SQL command for later input into a database.
This is primarily to allow for a comparison between the list of Montesinos tangles and the list of rational tangles.

This database table has the following columns:
component_index
montesinos_tangle_number
component_number
rational_p
rational_q


Function Inputs:
numOfComponents: the number of rational components used in the construction of the current Montesinos tangle.
compRationalPQ[][2]: an array storing the numerator p and denominator q of each of the rational components.

Function Outputs:
This function records the information in the MontesinosRationalComponentsSQL.txt, located in the same directory as this program.


This function calls:
	nothing
This function is called by:
	buildMontesinosFromComponents()

*/

void writeRationalComponentsToFileSQL(int numOfComponents, int compRationalPQ[][2]){
	
	FILE *tangleFile2;
	tangleFile2 = fopen("MontesinosRationalComponentsSQL.txt", "a");
	
	//Insert a row for each rational component.
	for(int i=0; i<numOfComponents; i++){
		fprintf(tangleFile2, "INSERT INTO montesinos_rational_components(component_index,montesinos_tangle_number,component_number,rational_p,rational_q) VALUES('%d.%d',%d,%d,%d,%d);\n", totalMontesinosTanglesFound, i+1, totalMontesinosTanglesFound, i, compRationalPQ[i][0], compRationalPQ[i][1]);
	}
	
	fclose(tangleFile2);
	
}




/*
(NC, 12/18/18)
This function records the information associated to an input 2-string Montesinos tangle in the form of a SQL command for later input into a database.

This database table has the following columns:
montesinos_tangle_number
crossing_number
tangle_parity
rational_components_number
rational_components_composition
construction
twist_vector_construction
gauss_code


Function Inputs:
"Input" arrays: the arrays storing the Gauss code informtion for the Montesinos tangle.
numOfComponents: the number of rational components which make up the Montesinos tangle.
composition[]: the integer composition describing the break down of the Montesinos tangle into its rational components.
compTwistVectors[][]: arrays storing the twist vectors of the rational components.
compTwistVectorsLength[]: an array storing the lengths of the twists vectors of the rational components.
compRationalPQ[][2]: array storing the numerator and denominator of the fraction corresponding to each rational component.

Function Outputs:
This function records the information in the MontesinosTangleListSQL.txt, located in the same directory as this program.


This function calls:
	nothing
This function is called by:
	buildMontesinosFromComponents()

*/

void writeMontesinosTangleToFileSQL(int crossingsInput, int *gaussInput, int *barsInput, int *orientedSignGaussInput, int parityInput, int *composition, int numOfComponents, int compTwistVectors[][MAXN-2], int *compTwistVectorsLength, int compRationalPQ[][2]){
	
	FILE *tangleFile1;
	tangleFile1 = fopen("MontesinosTangleListSQL.txt", "a");
	
	//montesinos_tangle_number, crossing_number
	fprintf(tangleFile1, "INSERT INTO montesinos_tangle_list(montesinos_tangle_number,crossing_number,tangle_parity,rational_components_number,rational_components_composition,construction,twist_vector_construction,gauss_code) VALUES(%d,%d,", totalMontesinosTanglesFound, crossingsInput);
	
	//tangle_parity
	if( (parityInput==0) || (parityInput==1) ){
		fprintf(tangleFile1, "'%d',", parityInput);
	} else if( parityInput==2 ) {
		fprintf(tangleFile1, "'infinity',");
	} else {
		fprintf(tangleFile1, "'unknown',");
	}
	
	//rational_components_number
	fprintf(tangleFile1, "%d,", numOfComponents);
	
	//rational_components_composition
	fprintf(tangleFile1, "'[");
	for(int i=0; i<numOfComponents; i++){
		fprintf(tangleFile1, " %d", composition[i]);
	}
	fprintf(tangleFile1, " ]',");
	
	//construction
	fprintf(tangleFile1, "'(%d/%d)", compRationalPQ[0][0], compRationalPQ[0][1]);
	for(int i=1; i<numOfComponents; i++){
		fprintf(tangleFile1, "+(%d/%d)", compRationalPQ[i][0], compRationalPQ[i][1]);
	}
	fprintf(tangleFile1, "',");
	
	//twist_vector_construction
	fprintf(tangleFile1, "'[");
	for(int i=0; i<compTwistVectorsLength[0]; i++){
		fprintf(tangleFile1, " %d", compTwistVectors[0][i]);
	}
	fprintf(tangleFile1, " ]");
	for(int j=1; j<numOfComponents; j++){
		fprintf(tangleFile1, "+[");
		for(int i=0; i<compTwistVectorsLength[j]; i++){
			fprintf(tangleFile1, " %d", compTwistVectors[j][i]);
		}
		fprintf(tangleFile1, " ]',");
	}
	
	//gauss_code
	if( parityInput == 0 ){
		
		//Use the denominator closure for a parity 0 tangle.
		fprintf(tangleFile1, "'a1-a2+a3-a4+|b1-)");
		
		if ( crossingsInput == 0 ){
			fprintf(tangleFile1, "(b2+b3-)");
		}
	
		for (int i = 0; i < (2*crossingsInput); i++){
			
			if (i == barsInput[0]){
				fprintf(tangleFile1, "(b2+b3-)");
			}
			
			if (gaussInput[i] > 0){
				fprintf(tangleFile1, "a");
			}else{
				fprintf(tangleFile1, "b");
			}

			fprintf(tangleFile1, "%d", abs(gaussInput[i])+4);
		
			if (orientedSignGaussInput[i] < 0){
				fprintf(tangleFile1, "-");
			}else{
				fprintf(tangleFile1, "+");
			}
		}
		
		fprintf(tangleFile1, "(b4+');\n");
		
	} else if( parityInput == 1 ){
		
		//Use the denominator closure for a parity 1 tangle.
		fprintf(tangleFile1, "'a1-a2-a3+a4+|b1-)");
		
		if ( crossingsInput == 0 ){
			fprintf(tangleFile1, "(b3+b2-)");
		}
	
		for (int i = 0; i < (2*crossingsInput); i++){
			
			if (i == barsInput[0]){
				fprintf(tangleFile1, "(b3+b2-)");
			}
			
			if (gaussInput[i] > 0){
				fprintf(tangleFile1, "a");
			}else{
				fprintf(tangleFile1, "b");
			}

			fprintf(tangleFile1, "%d", abs(gaussInput[i])+4);
		
			if (orientedSignGaussInput[i] < 0){
				fprintf(tangleFile1, "-");
			}else{
				fprintf(tangleFile1, "+");
			}
		}
		
		fprintf(tangleFile1, "(b4+');\n");
		
	} else if( parityInput == 2 ){
		
		//Use the numerator closure for a parity infinity tangle.
		fprintf(tangleFile1, "'a1-a2+a3-a4+|b1-)");
		
		if ( crossingsInput == 0){
			fprintf(tangleFile1, "(b4+b3-)");
		}
	
		for (int i = 0; i < (2*crossingsInput); i++){
			
			if (i == barsInput[0]){
				fprintf(tangleFile1, "(b4+b3-)");
			}
			
			if (gaussInput[i] > 0){
				fprintf(tangleFile1, "a");
			}else{
				fprintf(tangleFile1, "b");
			}

			fprintf(tangleFile1, "%d", abs(gaussInput[i])+4);
		
			if (orientedSignGaussInput[i] < 0){
				fprintf(tangleFile1, "-");
			}else{
				fprintf(tangleFile1, "+");
			}
		}
		
		fprintf(tangleFile1, "(b2+');\n");
		
	} else{
		fprintf(tangleFile1, "NULL);\n");
	}
	
	fclose(tangleFile1);
	
}





/*
(NC, 12/14/18)
This function constructs a generalized Montesinos tangle from a Montesinos tangle to use as the core and a twist vector to describe the endpoint twisting.
Horizontal and vertical twisting is equivalent to adding to or multiplying by a one crossing tangle.
There are two such one crossing tangles (the "atomic" tangles defined below): positive and negative, depending on the crossing sign.
These are used with the tangleSum() and tangleProduct() functions to together with the input Montesinos tangle to construct the generalized Montesinos tangle.
The type and ammount of twisting is all described in the twist vector.


Function Inputs:
twistVectorInput[]: an array storing the twist vector information.
twistVectorInputLength: the length of the twistVectorInput[] array.
"mont" arrays: the Gauss code information for the Montesinos tangle representative.
"GM" arrays: the Gauss code information for the generalized Montesinos tangle to be constructed.

Function Outputs (in the form of modified global variables):
The "GM" arrays are updated to describe the desired generalized Montesinos tangle.


This function calls:
	tangleSum()
	tangleProduct()

This function is called by:
	findGeneralizedMontesinosTangles()

*/

void twistMontesinosTangle(int *twistVectorInput, int twistVectorInputLength, int montCross, int montParity, int *montBars, int *montGauss, int *montGaussSign, int &crossGM, int *gaussGM, int *orientedSignGaussGM, int *barsGM, int &parityGM){
	
	//Initialize the Gauss code for the two atomic tangles:
	//Positive Slope: [1]
	int atomicPositiveGauss[2] = {1,-1};
	int atomicPositiveBars[2] = {1,2};
	int atomicPositiveOrientedSign[2] = {-1,-1};
	//Negative Slope: [-1]
	int atomicNegativeGauss[2] = {-1,1};
	int atomicNegativeBars[2] = {1,2};
	int atomicNegativeOrientedSign[2] = {1,1};
	//Note that both atomic tangles have parity 1 and only a single crossing, so there is no special reason to define variables for the crossing and parity of these.
	
	//Initialize the GM Gauss code arrays as the mont Gauss code arrays.
	crossGM=montCross;
	parityGM=montParity;
	barsGM[0]=montBars[0];
	barsGM[1]=montBars[1];
	for(int i=0; i<2*montCross; i++){
		gaussGM[i]=montGauss[i];
		orientedSignGaussGM[i]=montGaussSign[i];
	}
	
	//An integer to track the number of twists to be performed at each step, and a second integer to track the sign of the twists.
	int twistNum;
	int twistSign;
	
	//Iterate through each entry of the twist vector, performing the indicated twists for each entry.
	for(int t=0; t<twistVectorInputLength; t++){
		
		twistNum=abs(twistVectorInput[t]);
		
		if( twistVectorInput[t] < 0 ){
			twistSign=-1;
		} else {
			twistSign=1;
		}
		
		//By convention, the twist vectors will always have odd length, beginning and ending with a horizontal twist.
		//Since t begins at zero, we use horizontal twisting when t is even and vertical twisting when t is odd.
		if( (t%2) == 0 ){
			//Horizontal Twisting:
			for(int i=0; i<twistNum; i++){
				if( twistSign > 0 ){
					//Positive Twist:
					tangleSum(crossGM,gaussGM,barsGM,orientedSignGaussGM,parityGM,1,atomicPositiveGauss,atomicPositiveBars,atomicPositiveOrientedSign,1,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
				} else if( twistSign < 0 ){
					//Negative Twist:
					tangleSum(crossGM,gaussGM,barsGM,orientedSignGaussGM,parityGM,1,atomicNegativeGauss,atomicNegativeBars,atomicNegativeOrientedSign,1,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
				}
				
				//Set the GM Gauss code arrays equal to the SUM Gauss code arrays before continuing.
				crossGM=numOfCrossingsSUM;
				parityGM=paritySUM;
				barsGM[0]=barsSUM[0];
				barsGM[1]=barsSUM[1];
				for(int j=0; j<2*crossGM; j++){
					gaussGM[j]=gaussSUM[j];
					orientedSignGaussGM[j]=orientedSignGaussSUM[j];
				}
				//The for loop repeats until all twists have be incorporated.
				
			}
			
		} else {
			//Vertical Twisting: (identical to the above case, except using tangleProd and the PROD arrays.
			for(int i=0; i<twistNum; i++){
				if( twistSign > 0 ){
					//Positive Twist:
					tangleProduct(crossGM,gaussGM,barsGM,orientedSignGaussGM,parityGM,1,atomicPositiveGauss,atomicPositiveBars,atomicPositiveOrientedSign,1,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
				} else if( twistSign < 0 ){
					//Negative Twist:
					tangleProduct(crossGM,gaussGM,barsGM,orientedSignGaussGM,parityGM,1,atomicNegativeGauss,atomicNegativeBars,atomicNegativeOrientedSign,1,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
				}
				
				//Set the GM Gauss code arrays equal to the PROD Gauss code arrays before continuing.
				crossGM=numOfCrossingsPROD;
				parityGM=parityPROD;
				barsGM[0]=barsPROD[0];
				barsGM[1]=barsPROD[1];
				for(int j=0; j<2*crossGM; j++){
					gaussGM[j]=gaussPROD[j];
					orientedSignGaussGM[j]=orientedSignGaussPROD[j];
				}
				//The for loop repeats until all twists have be incorporated.
			}
			
		}
		//At the conclusion of this loop, all twists from entry t of the twist vector should be included.
	}
	//At the conclusion of this loop, all twists in the twist vector should be included.
	
}






/*
(NC, 12/14/18)
This function generates all Montesinos tangles (up to maxN crossing number) in the family with the representative Montesinos tangle defined in the "mont" Gauss code arrays.
It does this by fixing possible twist vectors, and then using these to determine how to twist the endpoints of the representative Montesinos tangle.
The actual twisting is computed by calling the twistMontesinosTangle() function.

Twist vectors follow the usual convention that the last entry indicates a horizontal twist and the number of entries is odd. Hence, the first twist vector entry also denotes a horizontal twist.
Each composition obtained from the compositionList[][][] arrays can be shifted one position right, making the first entry 0 and adding or removing a last entry of zero.
In this way, the first nonzero twisting is vertical, but the twist vector still satisfies the above convention.
Twist vector entries must have constant sign; negative twisting can happen, but only when the first nozero twist is vertical, else the tangle can be reduced.


Function Inputs:
maxN: the maximum crossing number of the generalizd Montesinos tangles to be constructed.
compositionList[][][]: the arrays storing possible integer compositions, to be used as twist vectors.
rationalTwistVectLength[][]: arrays storing the length of the corresponding rational twist vectors stored in the compositionList array.
numOfComponents: the number of components comprising representative Montesinos tangle.
compRationalPQ[][2]: the rational numerator and denominator of the fraction corresponding to the rational tangle components of the representative Montesinos tangle.
"mont" arrays: arrays to store the Gauss code information of the representative Montesinos tangle.
"GM" arrays: arrays to store the Gauss code information of the generalized Montesinos tangle.

Function Outputs:
A bunch of printf statements to report the constructed tangle, including:
program index, crossing number, parity, construction, Gauss code.

The generalized Montesinos tangles are printed with a program index of the form A.B.C.D, where:
A denotes the index of the representative Montesinos tangle.
B denotes the number of additional twists.
C denotes the composition number.
D denotes a shifted twist vector (1) or a negative shifted twist vector (2); otherwise D is 0 by default.


This function calls:
	ApowerB()
	twistMontesinosTangle()

This function is called by:
	buildMontesinosFromComponents()

*/

void findGeneralizedMontesinosTangles(int maxN, int compositionList[MAXN+1][CLB][MAXN+1], int rationalTwistVectLength[MAXN+1][CLB], int totalMontesinosTanglesFound, int numOfComponents, int compRationalPQ[][2], int montCross, int montParity, int *montBars, int *montGauss, int *montGaussSign, int &crossGM, int *gaussGM, int *orientedSignGaussGM, int *barsGM, int &parityGM){
	
	//The maximum number of endpoint twists which can be used is the difference of maxN with crossing number of core Montesinos tangle.
	int maxTwists=(maxN-montCross);
	//Initialize an an array to store local twist vectors.
	//Note that the length is "+2" to account for possible padding by zeros, which is necessary to generate some configurations while adhering to the rule that the first and last entries denote horizontal twists.
	int localTwistVector[maxTwists+2]={0};
	int localTwistVectorLength;
	
	int shiftTwistVector[maxTwists+2]={0};
	int shiftTwistVectorLength;
	
	int compositionMax;
	
	
	//While this is technically unessary, report the original Montesions tangle as a generalized Montesinos tanlge with a twist vector of [ 0 ].
	//This way, the table of generalized Montesinos tanlges will include the Montesinos tangles.
	printf("\n Generalized Montesinos Tangle %d.0.0.0:\n", totalMontesinosTanglesFound);
	printf("\n Crossings: \t %d", montCross);
	if( montParity == 2 ){
		printf("\n Parity: \t infinity");
	} else {
		printf("\n Parity: \t %d", montParity);
	}
	printf("\n Construction: \t (%d/%d)", compRationalPQ[0][0], compRationalPQ[0][1]);
	for(int j=1; j<numOfComponents; j++){
		printf("+(%d/%d)", compRationalPQ[j][0], compRationalPQ[j][1]);
	}
	printf(" o [ 0 ]\n");
	printTangleGaussCode(montCross,montGauss,montBars,montGaussSign,montParity);
	printf("\n ---------- \n");
	
	//Also record this information in the form of SQL input commands.
	//Record the tangle information as a SQL input command.
	writeGMTangleToFileSQL(totalMontesinosTanglesFound,0,0,0,montCross,montCross,montGauss,montBars,montGaussSign,montParity,numOfComponents,compRationalPQ,localTwistVector,0);
	
	
	//Iterate through the possible number of twists that can be included.
	for(int t=1; t<=maxTwists; t++){
		
		//For each choice of i, the number of compositions of t is 2^(t-1); this gives specifies the the number of times next loop iterates.
		compositionMax=ApowerB(2,(t-1));
		for(int i=0; i<compositionMax; i++){
			
			for(int j=0; j<rationalTwistVectLength[t][i]; j++){
				localTwistVector[j]=compositionList[t][i][j];
			}
			localTwistVectorLength=rationalTwistVectLength[t][i];
			
			//Construct the generalized Montesinos tangle by performing the indicated twists on the current Montesinos tangle.
			twistMontesinosTangle(localTwistVector,localTwistVectorLength,montCross,montParity,montBars,montGauss,montGaussSign,crossGM,gaussGM,orientedSignGaussGM,barsGM,parityGM);
			totalGMTanglesFound++;
			
			//DEBUG
			/*
			printf("\n Generalized Montesinos tangle:");
			printf("\n Crossing Number: %d", crossGM);
			printf("\n Parity: %d", parityGM);
			printf("\n Bars: [ %d %d ]", barsGM[0], barsGM[1]);
			printf("\n Gauss Array: [");
			for(int k=0; k<2*crossGM; k++){
				printf(" %d", gaussGM[k]);
			}
			printf(" ]");
			printf("\n Gauss OSigns: [");
			for(int k=0; k<2*crossGM; k++){
				printf(" %d", orientedSignGaussGM[k]);
			}
			printf(" ]\n");
			*/
			
			//Record the tangle information as a SQL input command.
			writeGMTangleToFileSQL(totalMontesinosTanglesFound,t,i,0,montCross,crossGM,gaussGM,barsGM,orientedSignGaussGM,parityGM,numOfComponents,compRationalPQ,localTwistVector,localTwistVectorLength);
			
			//Print out the tangle information.
			printf("\n Generalized Montesinos Tangle %d.%d.%d.0:\n", totalMontesinosTanglesFound, t, i);
			printf("\n Crossings: \t %d", crossGM);
			if( parityGM == 2 ){
				printf("\n Parity: \t infinity");
			} else {
				printf("\n Parity: \t %d", parityGM);
			}
			printf("\n Construction: \t (%d/%d)", compRationalPQ[0][0], compRationalPQ[0][1]);
			for(int j=1; j<numOfComponents; j++){
				printf("+(%d/%d)", compRationalPQ[j][0], compRationalPQ[j][1]);
			}
			printf(" o [");
			for(int j=0; j<localTwistVectorLength; j++){
				printf(" %d", localTwistVector[j]);
			}
			printf(" ]\n");
			printTangleGaussCode(crossGM,gaussGM,barsGM,orientedSignGaussGM,parityGM);
			printf("\n ---------- \n");
			
			
			//The entries in each twist vector can be shifted to so that the first nonzero twist is vertical.
			if( localTwistVector[localTwistVectorLength-1] == 0 ){
				//If the last entry is 0, shift all entries to right by one position, and set the first entry to 0.
				shiftTwistVectorLength=localTwistVectorLength;
				shiftTwistVector[0]=0;
				for(int j=1; j<shiftTwistVectorLength; j++){
					shiftTwistVector[j]=localTwistVector[j-1];
				}
			} else {
				//If the last entry is not zero, then pad the first and last entries by a zero.
				shiftTwistVectorLength=localTwistVectorLength+2;
				shiftTwistVector[0]=0;
				shiftTwistVector[shiftTwistVectorLength-1]=0;
				for(int j=1; j<shiftTwistVectorLength-1; j++){
					shiftTwistVector[j]=localTwistVector[j-1];
				}
			}
			
			twistMontesinosTangle(shiftTwistVector,shiftTwistVectorLength,montCross,montParity,montBars,montGauss,montGaussSign,crossGM,gaussGM,orientedSignGaussGM,barsGM,parityGM);
			totalGMTanglesFound++;
			
			//Record the tangle information as a SQL input command.
			writeGMTangleToFileSQL(totalMontesinosTanglesFound,t,i,1,montCross,crossGM,gaussGM,barsGM,orientedSignGaussGM,parityGM,numOfComponents,compRationalPQ,shiftTwistVector,shiftTwistVectorLength);
			
			//Print out the tangle information.
			printf("\n Generalized Montesinos Tangle %d.%d.%d.1:\n", totalMontesinosTanglesFound, t, i);
			printf("\n Crossings: \t %d", crossGM);
			if( parityGM == 2 ){
				printf("\n Parity: \t infinity");
			} else {
				printf("\n Parity: \t %d", parityGM);
			}
			printf("\n Construction: \t (%d/%d)", compRationalPQ[0][0], compRationalPQ[0][1]);
			for(int j=1; j<numOfComponents; j++){
				printf("+(%d/%d)", compRationalPQ[j][0], compRationalPQ[j][1]);
			}
			printf(" o [");
			for(int j=0; j<shiftTwistVectorLength; j++){
				printf(" %d", shiftTwistVector[j]);
			}
			printf(" ]\n");
			printTangleGaussCode(crossGM,gaussGM,barsGM,orientedSignGaussGM,parityGM);
			printf("\n ---------- \n");
			
			
			//TWist vcetors must have constant sign; negative twists can also be used, but only when the first nonzero twist is vertical. Else, the tangle is reducible.
			for(int j=0; j<shiftTwistVectorLength; j++){
				shiftTwistVector[j] = -1*shiftTwistVector[j];
			}
			
			twistMontesinosTangle(shiftTwistVector,shiftTwistVectorLength,montCross,montParity,montBars,montGauss,montGaussSign,crossGM,gaussGM,orientedSignGaussGM,barsGM,parityGM);
			totalGMTanglesFound++;
			
			//Record the tangle information as a SQL input command.
			writeGMTangleToFileSQL(totalMontesinosTanglesFound,t,i,2,montCross,crossGM,gaussGM,barsGM,orientedSignGaussGM,parityGM,numOfComponents,compRationalPQ,shiftTwistVector,shiftTwistVectorLength);
			
			//Print out the tangle information.
			printf("\n Generalized Montesinos Tangle %d.%d.%d.2:\n", totalMontesinosTanglesFound, t, i);
			printf("\n Crossings: \t %d", crossGM);
			if( parityGM == 2 ){
				printf("\n Parity: \t infinity");
			} else {
				printf("\n Parity: \t %d", parityGM);
			}
			printf("\n Construction: \t (%d/%d)", compRationalPQ[0][0], compRationalPQ[0][1]);
			for(int j=1; j<numOfComponents; j++){
				printf("+(%d/%d)", compRationalPQ[j][0], compRationalPQ[j][1]);
			}
			printf(" o [");
			for(int j=0; j<shiftTwistVectorLength; j++){
				printf(" %d", shiftTwistVector[j]);
			}
			printf(" ]\n");
			printTangleGaussCode(crossGM,gaussGM,barsGM,orientedSignGaussGM,parityGM);
			printf("\n ---------- \n");
			
		}
		
	}
	
}




/*
(NC, 12/9/18)
After fixing a composition for the Montesinos tangle, and fixing compositions for each of the rational components, this function is called to build all of the subtangles and put them together into a Montesinos tangle.
The "mont" Gauss code arrays are used to store Gauss code information of the Montesinos tangle as it is built up one component at a time.
Each rational subtangle component is built from its input twist vector using buildTwistTangle(); the Gauss code information for this component is stored in the "twist" Gauss code arrays.
This component is the summed with the existing Montesinos tangle in the "mont" arrays using the tangleSum() function,
This process repeats, iterating through each of the rational components. The final Montesinos tangle Gauss code information is stored across the "mont" arrays.

Note that some combinations of components might not yield a valid tangle if the construction would sum two parity infinity tangles.
Before tangleSum() is called, the parity of the tangles being summed is checked first; if they are both parity infinity, the construction is not possible and so the function terminates here.
If the construction is valid, it prints out the construction information and writes it to a file for later input into a MySQL database.

Next, the function moves on to generating all generalized Montesinos tangles which have this Montesinos tangle as a core representative.


Function Inputs:
compTwistVectors[][]: the multi-array storing the twist vectors corresponding to each of the rational components.
compTwistVectorsLength[]: the array storing the length (number of nonzero terms) of each twist vector in the composition.
montComposition[]: the array storing the composition of the Montesinos tangle currently being constructed.
numOfComponents: the number of rational components used in constructing this Montesinos tangle.


This function calls:
	buildTwistTangle()
	calculateExtendedFractionPQ()
	tangleSum()
	printTangleGaussCode()
	writeMontesinosTangleToFileSQL()
	findGeneralizedMontesinosTangles()

This function is called by:
	findRationalComponents

*/

void buildMontesinosFromComponents(int compTwistVectors[][MAXN-2], int *compTwistVectorsLength, int *montComposition, int numOfComponents, int maxN){
	
	//Initialize the Montesinos tangle Gauss code arrays as the zero tangle.
	montCross=0;
	montBars[0]=0;
	montBars[1]=0;
	montParity=0;
	
	//Build up the remaining rational components one at a time from their twist vectors and iteratively sum these to build up the full Montesinos tangle.
	for(int i=0; i<numOfComponents; i++){
		buildTwistTangle(compTwistVectorsLength[i],compTwistVectors[i],twistCrossings,gaussTwist,barsTwist,orientedSignGaussTwist,parityTwist,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
		calculateExtendedFractionPQ(true,compTwistVectorsLength[i],compTwistVectors[i],compRationalPQ[i]);
		
		//Sum the component constructed above with the existing Montesinos tangle.
		//Note that the rational component obtained above is stored across the "Twist" tangle arrays.
		//We must first check to make sure that we are not summing two parity infinity tangles. If we are, this is not a valid construction of a Montesinos tangle and we must eliminate it.
		//We check that at least one of the parities of the two tangles being summed is NOT equal to 2 (which we recall is used to denote infinity).
		if( (montParity!=2) || (parityTwist!=2) ){
			tangleSum(montCross,montGauss,montBars,montGaussSign,montParity,twistCrossings,gaussTwist,barsTwist,orientedSignGaussTwist,parityTwist,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
		} else {
			//If the current construction would cause two parity infinity tangles to be summed, the construction is invalid since that would not result in a 2-string tangle.
			//In this case, exit the loop and move on to considering the next possible construction.
			return;
		}
		
		//Assuming the construction is still valid, we update the Montesinos tangles arrays.
		//The tangle sum obtained above is stored in the "SUM" tangle arrays.
		//Overwrite the "mont" tangle arrays with this information.
		montCross=numOfCrossingsSUM;
		montParity=paritySUM;
		montBars[0]=barsSUM[0];
		montBars[1]=barsSUM[1];
		for(int j=0; j<(2*numOfCrossingsSUM); j++){
			montGauss[j]=gaussSUM[j];
			montGaussSign[j]=orientedSignGaussSUM[j];
		}
			
	}
	
	totalMontesinosTanglesFound++;
	totalGMTanglesFound++;
	
	//Print all of the information regarding the constructed Montesinos tangle.
	printf("\n ********** \n ");
	printf("\n Montesinos Tangle %d:\n", totalMontesinosTanglesFound);
	printf("\n Crossings: \t %d", montCross);
	if(montParity==2){
		printf("\n Parity: \t infinity");
	} else {
		printf("\n Parity: \t %d", montParity);
	}
	printf("\n Components: \t %d", numOfComponents);
	printf("\n Composition: \t [");
	for(int j=0; j<numOfComponents; j++){
		printf(" %d", montComposition[j]);
	}
	printf(" ]\n");
	printf("\n Construction: \t (%d/%d)", compRationalPQ[0][0], compRationalPQ[0][1]);
	for(int j=1; j<numOfComponents; j++){
		printf("+(%d/%d)", compRationalPQ[j][0], compRationalPQ[j][1]);
	}
	printf("\n\n TwistVectors: \t [");
	for(int k=0; k<compTwistVectorsLength[0]; k++){
		printf(" %d", compTwistVectors[0][k]);
	}
	printf(" ]");
	for(int j=1; j<numOfComponents; j++){
		printf("+[");
		for(int k=0; k<compTwistVectorsLength[j]; k++){
			printf(" %d", compTwistVectors[j][k]);
		}
		printf(" ]");
	}
	printf("\n");
	printTangleGaussCode(montCross,montGauss,montBars,montGaussSign,montParity);
	
	printf("\n ---------- \n");
	
	
	//Write the information for this tangle to a file in the form a SQL input command for later use importing into a database.
	writeMontesinosTangleToFileSQL(montCross,montGauss,montBars,montGaussSign,montParity,montComposition,numOfComponents,compTwistVectors,compTwistVectorsLength,compRationalPQ);
	
	//Write the information pertaining to the component rational tangles in the form of a SQL input command.
	writeRationalComponentsToFileSQL(numOfComponents,compRationalPQ);
	
	
	//Next, use the Montesinos tangle found above as the core to generate possible generalized Montesinos tangles.
	findGeneralizedMontesinosTangles(maxN,compositionList,rationalTwistVectLength,totalMontesinosTanglesFound,numOfComponents,compRationalPQ,montCross,montParity,montBars,montGauss,montGaussSign,crossGM,gaussGM,orientedSignGaussGM,barsGM,parityGM);
	
		
}




/*
(NC, 12/9/18)
For a given Montesinos composition, this function searches for all rational tangle components which could be used in this composition.
There are many possible integer compositions which can be used as twist vectors for each components, and the number of components can vary.
Thus, this function searches for one component at a time and then calls itself to move on to the next component, until it reaches the last component.
A rational component is specified by a choice of twist vector, which is used to explicitly construct the component.
After choosing the last component, this function calls buildMontesinosFromComponents() to actually construct the Gauss code for the Montesinos tangle with these rational components.

Function Inputs:
componentIndex: the index the rational component being considered in the current composition.
componentIndexMax: the maximum component index, used to stop the recursive part of this function after all components are fixed.
compTwistVectors[][]: a multi-array used to store the twist vector corresponding to each component, used to construct the rational subtangle later.
compTwistVectorsLength[]: an array storing the length (number of nonzero entries) in each of the twist vectors stored in the above multi-array.
montComposition[]: an array storing the composition of the Montesinos tangle currently being considered.
compositionList[][][]: the multi-array storing all possible compositions.
componentCompositionMaxCount[]: an array storing the number of possible compositions for each component; this is used to specify an upper bound for each of the for-loops which iterate through compositions.


This function calls:
	itslef (recursive)
	buildMontesinosFromComponents()

This function is called by:
	findMontesinosTangles()

*/

void findRationalComponents(int componentIndex, int componentIndexMax, int compTwistVectors[][MAXN-2], int *compTwistVectorsLength, int *montComposition, int compositionList[MAXN+1][CLB][MAXN+1], int *componentCompositionMaxCount, int maxN){
	
	//We consider compositions corresponding to a rational tangle with this many crossings.
	int localRationalCrossingNumber=montComposition[componentIndex];
	
	if( componentIndex < componentIndexMax-1 ){
		//In this case, we only use positive twist vectors and we call the function again to construct the next component.
		
		//Iterate through all possible compositions of the appropriate integer.
		for(int j=0; j<componentCompositionMaxCount[componentIndex]; j++){
			
			//Set the corresponding entries of the compTwistVectors[componentIndex][] array to match the current composition.
			compTwistVectorsLength[componentIndex]=rationalTwistVectLength[localRationalCrossingNumber][j];
			for(int i=0; i<compTwistVectorsLength[componentIndex]; i++){
				compTwistVectors[componentIndex][i]=compositionList[localRationalCrossingNumber][j][i];			
			}
			
			//The component rational tangles are subject to the convention that the last entry of the corresponding twist vector is 0.
			//If this condition is met, the function calls itself again to generate the next rational component. Otherwise, it skips this composition.
			if( compTwistVectors[componentIndex][compTwistVectorsLength[componentIndex]-1] == 0 ){
				findRationalComponents(componentIndex+1,componentIndexMax,compTwistVectors,compTwistVectorsLength,montComposition,compositionList,componentCompositionMaxCount,maxN);
			}
		}
		
		
	} else if( componentIndex == componentIndexMax-1 ){
		//This happens when constructing the last component. In this special case, twist vectors can be positive or negative.
		//This is also where actual construction of the Montesinos tangle is done.
		
		//Iterate through all possible compositions of the appropriate integer.
		for(int j=0; j<componentCompositionMaxCount[componentIndex]; j++){
			
			//Set the corresponding entries of the comptTwistVectors[componentIndex][] array to match the current composition.
			compTwistVectorsLength[componentIndex]=rationalTwistVectLength[localRationalCrossingNumber][j];
			for(int i=0; i<compTwistVectorsLength[componentIndex]; i++){
				compTwistVectors[componentIndex][i]=compositionList[localRationalCrossingNumber][j][i];			
			}
			
			//Check to see if the last twist vector entry is 0.
			//Technically, we don't need this condition on the last rational component, but it will allow us a unique representative of the core subtangle.
			if( compTwistVectors[componentIndex][compTwistVectorsLength[componentIndex]-1] == 0 ){
				
				//DEBUG
				/*
				printf("\n Montesinos Tangle Decomposition:\n montComposition: [");
				for(int l=0; l<componentIndexMax; l++){
					printf(" %d ", montComposition[l]);
				}
				printf("]");
				printf("\n Rational Components:");
				for(int l=0; l<componentIndexMax; l++){
					printf("\n Component %d: \t [", l);
					for(int r=0; r<compTwistVectorsLength[l]; r++){
						printf(" %d ", compTwistVectors[l][r]);
					}
					printf("]");
				}
				printf("\n");
				*/
				
				//Construct the Montesinos tangle corresponding to this choice of rational components.
				buildMontesinosFromComponents(compTwistVectors,compTwistVectorsLength,montComposition,componentIndexMax,maxN);				
				
				//Change the sign of the last twist vector and construct the other Montesinos tangle corresponding to thiese rational components.
				//(UPDATE 12/13/18) NOTE NEEDED-OPPOSITE SIGNS CAN BE ACCOUNTED FOR IN THE ENDPOINT TWISTING
				/*
				for(int i=0; i<compTwistVectorsLength[componentIndex]; i++){
					compTwistVectors[componentIndex][i]=(-1*compTwistVectors[componentIndex][i]);
				}
				buildMontesinosFromComponents(compTwistVectors,compTwistVectorsLength,montComposition,componentIndexMax);
				*/
				
			}
			
		}
		
		
	} else {
		//Error case, should never happen.
		return;
	}
	
}



/*
(NC, 12/9/18)
For each crossing number up to the specified maximum, this function generates all possible compositions of this number which could correspond to a Montesinos tangle.
The number of terms in a composition is the number of possible rational components; each of these terms is the crossing number of a rational subtangle.
At minumum, there most be at least two components (otherwise the Montesinos tangle is rational), and each component must have at least two crossings.
This function calls findPositiveTwistVectors() to generate all possible compositions. It then iterates through each composition and checks if the above conditions are satisfied.
If the composition is good, the function calls findRationalComponents() to iterate through all possible rational tangles which could be used as components of the Montesinos tangle.

Function Inputs:
maxN: the maximum crossing number.
minN: the minimum crossing number (5 at the smallest for a non-rational Montesinos tangle, but can be increased to generate partial lists).
compositionList[][][]: a multi-array used to store all possible compositions.
rationalTwistVectLength[][]: a multi-array used to store the length (number of non-zero terms) of a composition used as a twist vector.
sameSignTwistVector[][]: a storage array for changing the signs on the terms in a twist vector; needed for other functions.
totalRationalTanglesFound: an artifact from previous versions of this function that is no longer used?
rationalPQ: another artifact
compRationalPQ: an array used to store the rational number corresponding to the rational tangle components.

Global arrays for manipulating the Gauss code of tangles across multiple functions: twist, SUM, PROD, comp, mont
For each of thes, information is stored on: crossing number, parity, gauss, oriented sign gauss, bars.


This function calls:
	findPositiveTwistVectors()
	findRationalComponents()

This function is called by:
	main()

*/

void findMontesinosTangles(int maxN, int minN, int compositionList[MAXN+1][CLB][MAXN+1], int rationalTwistVectLength[MAXN+1][CLB], int sameSignTwistVector[2][MAXN+1], int &totalRationalTanglesFound, int *rationalPQ, int &twistCrossings, int *gaussTwist, int *barsTwist, int *orientedSignGaussTwist, int &parityTwist, int &numOfCrossingsSUM, int *gaussSUM, int *barsSUM, int *orientedSignGaussSUM, int &paritySUM, int &numOfCrossingsPROD, int *gaussPROD, int *barsPROD, int *orientedSignGaussPROD, int &parityPROD, int *compCross, int compGauss[][2*MAXN], int compGaussSign[][2*MAXN], int compBars[][2], int *compParity, int compRationalPQ[][2], int montCross, int *montGauss, int *montGaussSign, int *montBars, int montParity){
	
	//printf("\n Hello again, World!");
	
	int unsignedCompositionNumber=0;
	for(int i=0; i<=maxN; i++){
		unsignedCompositionNumber=unsignedCompositionNumber+countCompositions(i);
	}
	
	//printf("\n *************** \n\n Total number of unsigned compositions of integers up to %d: %d \n", maxN, unsignedCompositionNumber);
	
	//Construct the array of matrices storing unsigned compositions of integers up to maxN.
	findPositiveTwistVectors(compositionList,rationalTwistVectLength,maxN);
	
	//Since each rational component must have at least two crossings, the floor of maxN/2 is an upper bound for the number of possible rational components.
	//This is used to place a bound on the size of the multi-arrays needed to store the components.
	int maxComponents = floor(maxN/2);
	
	//Local arrays to store twist vectors. Used to generate the components.
	int compTwistVectors[maxComponents][MAXN-2];
	int compTwistVectorsLength[maxComponents];
	
	int compositionMax;
	int compNumLocal;
	bool goodComposition;
	
	//Iterate first through the crossing number of the final Montesinos tangle, set to be maxN by default.
	for(int n=minN; n<=maxN; n++){
		
		printf("\n ********** \n List of Montesinos tangles with %d crossings \n ********** \n", n);
		
		//The total number of compositions for each choice of n is 2^(n-1).
		compositionMax=ApowerB(2,(n-1));
		
		//Iterate through the possible compositions, stored in compositionList[n][][].
		//NOTE: The we start with index i=1 to exclude considering the possibility of a Monetesions tangle with one component, which is always listed as the [n][0][] twist vector.
		for(int i=1; i<compositionMax; i++){
			
			//Determine the number of nonzero entries in the composition.
			if( compositionList[n][i][rationalTwistVectLength[n][i]-1] == 0 ){
				compNumLocal = (rationalTwistVectLength[n][i] - 1);
			} else {
				compNumLocal = rationalTwistVectLength[n][i];
			}
			
			//Set a boolen variable "goodComposition" to determine if this composition specifies a corresponding Montesinos tangle.
			//The composition is good as long as each entry is at least 2.
			goodComposition=true;
			for(int j=0; j<compNumLocal; j++){
				if( compositionList[n][i][j] < 2 ){
					goodComposition=false;
				}
			}
			
			
			//DEBUG
			/*
			if( goodComposition ){
				printf("\n Composition: compNumLocal = %d , [", compNumLocal);
				for(int j=0; j<compNumLocal; j++){
					printf(" %d ", compositionList[n][i][j]);
				}
				printf("]\n");
			}
			*/
			
			int montComposition[compNumLocal];
			for(int j=0; j<compNumLocal; j++){
				montComposition[j]=compositionList[n][i][j];
			}
			
			//If this is a good composition, use this to determine possible rational components.
			if( goodComposition ){
				
				//Determine the number of possible rational tangles for each component.
				int componentCompositionMaxCount[compNumLocal];
				for(int j=0; j<compNumLocal; j++){
					componentCompositionMaxCount[j]=ApowerB(2,(montComposition[j]-1));
				}
				
				//DEBUG
				/*
				printf("\n montComposition: [");
				for(int l=0; l<compNumLocal; l++){
					printf(" %d ", montComposition[l]);
				}
				printf("]");
				printf("\n componentCompositionMaxCount: [");
				for(int l=0; l<compNumLocal; l++){
					printf(" %d ", componentCompositionMaxCount[l]);
				}
				printf("]\n");
				*/
				
				//Call recursive function to move through all possible combinations of rational components.
				//The actual construction of the Montesinos tangles is performed in the called functions.
				findRationalComponents(0,compNumLocal,compTwistVectors,compTwistVectorsLength,montComposition,compositionList,componentCompositionMaxCount,maxN);
				
			}
			
		}
		
	}
	
}






int main(){
	
	/*
	The defualt maximum number of crossings considered (for a given component rational tangle) is 9.
	This can be increased if desired, but if doing so, one must be sure to update dimensions of the global arrays where they are defined in the preamble AND where they are defined in the function inputs.
	In general, if an integer n is chosen to be the maximum number of crossings, the global arrays should have at least the following dimensions to store the necessary data:
	
	int compositionList[n+1][2^(n-1)][n+1];
	int rationalTwistVectLength[n][2^(n-1)];
	int twistSigns[2^n][n+1];
	int sameSignTwistVector[2][n+1];
	
	The default maximum number of component rational tangles is 9.
	The Gauss code for the component rational tangles used are constructed from twist vectors as needed. These Gauss codes are stored across several multi-arrays.
	This number can be increased if desired, but the corresponding dimensions on arrays and functions must be updated.
	In general, if an integer k is chosen to be the maximum number of components, the global arrays should have at least the following dimensions to store the neecessary data:
	
	int compCross[k]
	int compGauss[k][2*n]
	int compGaussSign[k][2*n]
	int compBars[k][2]
	int compParity[k]
	*/
	
	//Maximum and minimum crossing numbers of the rational components. Can be altered to restrict the range of crossing numbers of Montesinos tangles considered.
	//Note that the smallest non-rational Montesinos tangle has 5 crossings. The value of maxN cannot exceed MAXN, else the array dimensions will be too small.
	int maxN=6;
	//Maximum number of component rational tangles.
	int minN=5;
	
	//printf("\n Hello, World!");
	
	findMontesinosTangles(maxN,minN,compositionList,rationalTwistVectLength,sameSignTwistVector,totalMontesinosTanglesFound,rationalPQ,twistCrossings,gaussTwist,barsTwist,orientedSignGaussTwist,parityTwist,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD,compCross,compGauss,compGaussSign,compBars,compParity,compRationalPQ,montCross,montGauss,montGaussSign,montBars,montParity);
	
	printf("\n ********** \n Total tangles found: %d\n", totalGMTanglesFound);
	
	
}



