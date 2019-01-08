#include <stdio.h>
#include <math.h>

int gaussSUM[30];
int barsSUM[30];
int orientedSignGaussSUM[30];
int numOfCrossingsSUM;
int paritySUM;
int gaussPROD[30];
int barsPROD[30];
int orientedSignGaussPROD[30];
int numOfCrossingsPROD;
int parityPROD;

int twistCrossings;
int gaussTwist[30];
int barsTwist[30];
int orientedSignGaussTwist[30];
int parityTwist;



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
		printf("\n\n BEFORE RE-LABEL");
		printf("\n strand2reverse: [");
		for(int i=0;i<strand2Length;i++){
			printf(" %d ", strand2reverse[i]);
		}
		printf("]\n");
		
		
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
		printf("\n AFTER RE-LABEL");
		printf("\n strand2reverse: [");
		for(int i=0;i<strand2Length;i++){
			printf(" %d ", strand2reverse[i]);
		}
		printf("]\n");
		
		
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
*/

void printTangleGaussCode(int numOfCrossingsInput, int *gaussInput, int *barsInput, int *orientedSignGaussInput, int parityInput){
	
	//In each case, the Gauss code constructed is the Gauss code of the oriented link obtained by taking either the numerator or denominator closure of the tangle.
	//It lists the outer circle first, oriented clockwise and beginning with the crossing at NW endpoint, and then travels along the strand starting from this same NW endpoint crossing.
	
	printf("\n Gauss Code: ");
	if( (parityInput == 0) || (parityInput == 1)){
		printf("parity %d tangle with %d crossings\n", parityInput,numOfCrossingsInput);
	} else if( parityInput == 2 ){
		printf("parity infinity tangle with %d crossings\n",numOfCrossingsInput);
	} else {
		printf("unknown parity tangle with %d crossings\n",numOfCrossingsInput);
	}
	
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
	printf("\n\n TWIST TANGLE:\n");
	printf("\n twistVector: [");
	for(int i=0; i<twistVectLength; i++){
		printf(" %d ", twistVector[i]);
	}
	printf("]");
	printTangleGaussCode(twistCrossings,gaussTwist,barsTwist,orientedSignGaussTwist,parityTwist);
	
	
	
}






int main(){
	
	//Test inputs
	
	int numOfCrossings1=3;
	int gauss1[6] = {-1,2,-3,1,-2,3};
	int bars1[2] = {2,6};
	int orientedSignGauss1[6] = {-1,-1,-1,-1,-1,-1};
	int parity1=0;
	
	//printTangleGaussCode(numOfCrossings1,gauss1,bars1,orientedSignGauss1,parity1);
	
	
	int numOfCrossings2=3;
	int gauss2[6] = {1,-2,3,-1,2,-3};
	int bars2[2] = {2,6};
	int orientedSignGauss2[6] = {1,1,1,1,1,1};
	int parity2=0;
	
	
	int numOfCrossings3=3;
	int gauss3[6] = {-1,2,-3,3,-2,1};
	int bars3[2] = {3,6};
	int orientedSignGauss3[6] = {1,1,1,1,1,1};
	int parity3=1;
	
	int twistVectLength1=1;
	int twistVector1[1] = {-3};
	
	//printTangleGaussCode(numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3);
	
	//buildTwistTangle(twistVectLength1,twistVector1,twistCrossings,gaussTwist,barsTwist,orientedSignGaussTwist,parityTwist,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
	
	
	
	int numOfCrossings4=3;
	int gauss4[6] = {-1,2,-3,1,-2,3};
	int bars4[2] = {2,6};
	int orientedSignGauss4[6] = {1,1,1,1,1,1};
	int parity4=2;
	
	
	int numOfCrossings5=3;
	int gauss5[6] = {-1,2,3,1,-2,-3};
	int bars5[2] = {2,6};
	int orientedSignGauss5[6] = {-1,-1,1,-1,-1,1};
	int parity5=0;
	
	
	int numOfCrossings6=2;
	int gauss6[4] = {-1,2,-2,1};
	int bars6[2] = {2,4};
	int orientedSignGauss6[4] = {1,1,1,1};
	int parity6=0;
	
	
	int numOfCrossings7=3;
	int gauss7[6] = {1,-2,-3,-1,2,3};
	int bars7[2] = {2,6};
	int orientedSignGauss7[6] = {-1,-1,1,-1,-1,1};
	int parity7=2;
	
	
	int numOfCrossings8=3;
	int gauss8[6] = {1,-2,3,-3,2,-1};
	int bars8[2] = {3,6};
	int orientedSignGauss8[6] = {-1,-1,-1,-1,-1,-1};
	int parity8=1;
	
	
	
	//Too many crossings, this one is tedious to check by hand
	int numOfCrossings9=6;
	int gauss9[12] = {-1,2,-3,4,-5,1,-6,3,-4,6,-2,5};
	int bars9[2] = {9,12};
	int orientedSignGauss9[12] = {-1,-1,1,1,-1,-1,1,1,1,1,-1,-1};
	int parity9=1;
	
	//printTangleGaussCode(numOfCrossings9,gauss9,bars9,orientedSignGauss9,parity9);
	
	//This is the [0] tangle; the SUM and PROD functions still work with tangles of zero crossings.
	int numOfCrossings0=0;
	int gauss0[0];
	int bars0[2] = {0,0};
	int orientedSignGauss0[0];
	int parity0=0;
	
	//printTangleGaussCode(numOfCrossings0,gauss0,bars0,orientedSignGauss0,parity0);
	
	
	
	//Below is a list of test combinations of tangles of various parities using either the sum or product function.
	//These are used as sample cases for debugging and each has been verified by hand.
	
	//tangles: 0+0; parities: 0,0
	//tangleSum(numOfCrossings0,gauss0,bars0,orientedSignGauss0,parity0,numOfCrossings0,gauss0,bars0,orientedSignGauss0,parity0,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 1+0; parities: 0,0
	//tangleSum(numOfCrossings1,gauss1,bars1,orientedSignGauss1,parity1,numOfCrossings0,gauss0,bars0,orientedSignGauss0,parity0,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 1+2; parities: 0,0
	//tangleSum(numOfCrossings1,gauss1,bars1,orientedSignGauss1,parity1,numOfCrossings2,gauss2,bars2,orientedSignGauss2,parity2,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 1+3; parities: 0,1
	//tangleSum(numOfCrossings1,gauss1,bars1,orientedSignGauss1,parity1,numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 1+4; parities: 0,2
	//tangleSum(numOfCrossings1,gauss1,bars1,orientedSignGauss1,parity1,numOfCrossings4,gauss4,bars4,orientedSignGauss4,parity4,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 3+1; parities: 1,0
	//tangleSum(numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossings1,gauss1,bars1,orientedSignGauss1,parity1,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 3+2; parities: 1,0
	//tangleSum(numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossings2,gauss2,bars2,orientedSignGauss2,parity2,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 3+5; parities: 1,0
	//tangleSum(numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossings5,gauss5,bars5,orientedSignGauss5,parity5,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 3+3; parities: 1,1
	//tangleSum(numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 3+6; parities: 1,0
	//tangleSum(numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossings6,gauss6,bars6,orientedSignGauss6,parity6,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 3+4; parities: 1,2
	//tangleSum(numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossings4,gauss4,bars4,orientedSignGauss4,parity4,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 3+7; parities: 1,2
	//tangleSum(numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossings7,gauss7,bars7,orientedSignGauss7,parity7,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 4+1; parities: 2,0
	//tangleSum(numOfCrossings4,gauss4,bars4,orientedSignGauss4,parity4,numOfCrossings1,gauss1,bars1,orientedSignGauss1,parity1,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 3+8; parities: 1,1
	//tangleSum(numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossings8,gauss8,bars8,orientedSignGauss8,parity8,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 4+3; parities: 2,1
	//tangleSum(numOfCrossings4,gauss4,bars4,orientedSignGauss4,parity4,numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 4+7; parities: 2,2
	//tangleSum(numOfCrossings4,gauss4,bars4,orientedSignGauss4,parity4,numOfCrossings7,gauss7,bars7,orientedSignGauss7,parity7,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM);
	
	//tangles: 1*3; parities: 0,1
	//tangleProduct(numOfCrossings1,gauss1,bars1,orientedSignGauss1,parity1,numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
	
	//tangles: 1*4; parities: 0,2
	//tangleProduct(numOfCrossings1,gauss1,bars1,orientedSignGauss1,parity1,numOfCrossings4,gauss4,bars4,orientedSignGauss4,parity4,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
	
	//tangles: 3*1; parities: 1,0
	//tangleProduct(numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossings1,gauss1,bars1,orientedSignGauss1,parity1,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
	
	//tangles: 9*1; parities: 1,0
	//tangleProduct(numOfCrossings9,gauss9,bars9,orientedSignGauss9,parity9,numOfCrossings1,gauss1,bars1,orientedSignGauss1,parity1,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
	
	//tangles: 3*3; parities: 1,1
	//tangleProduct(numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
	
	//tangles: 3*4; parities: 1,2
	//tangleProduct(numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossings4,gauss4,bars4,orientedSignGauss4,parity4,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
	
	//tangles: 3*9; parities: 1,1
	//tangleProduct(numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossings9,gauss9,bars9,orientedSignGauss9,parity9,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
	
	//tangles: 4*1; parities: 2,0
	//tangleProduct(numOfCrossings4,gauss4,bars4,orientedSignGauss4,parity4,numOfCrossings1,gauss1,bars1,orientedSignGauss1,parity1,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
	
	//tangles: 4*3; parities: 2,1
	//tangleProduct(numOfCrossings4,gauss4,bars4,orientedSignGauss4,parity4,numOfCrossings3,gauss3,bars3,orientedSignGauss3,parity3,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
	
	//tangles: 4*4; parities: 2,2
	//tangleProduct(numOfCrossings4,gauss4,bars4,orientedSignGauss4,parity4,numOfCrossings4,gauss4,bars4,orientedSignGauss4,parity4,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
	
	
	
	//Twist Vector Tests:
	
	int twistVectorA[3] = {2,-1,0};
	
	buildTwistTangle(3,twistVectorA,twistCrossings,gaussTwist,barsTwist,orientedSignGaussTwist,parityTwist,numOfCrossingsSUM,gaussSUM,barsSUM,orientedSignGaussSUM,paritySUM,numOfCrossingsPROD,gaussPROD,barsPROD,orientedSignGaussPROD,parityPROD);
	
	
	
	/*
	
	int testNumber1;
	int testNumber2;
	
	testNumber1 = shiftIndexMagnitude(5,3);
	testNumber2 = shiftIndexMagnitude(-5,3);
	
	printf("\n\n Test 1: %d , Test 2: %d ", testNumber1, testNumber2);
	
	*/
	
}


















