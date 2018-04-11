//  Created by alastair.maxwell@glasgow.ac.uk on 22/03/2018.
//  Copyright © 2018 University of Glasgow. All rights reserved.
//  HDSP -- simple DSP module for Huntington Disease amplicon data
//  Implementation of python functions from ScaleHD in C++

// Compilation information
// OS X 10.11.6 G++ Apple LLVM 7.3.0 (clang-703.0.31)
// Seqan header library version 2.4.0, ZLIB 1.2.11 
// $ g++-7 -std=c++17 -lz -D SEQAN_HAS_ZLIB -fconcepts -I /usr/local/include

#include <seqan/bam_io.h>                                                                                               // Seqan library (requires ZLIB in /usr/local/include)
#include <sys/stat.h>
#include <algorithm>                                                                                                    // Stat library for determining file existence
#include <iostream>                                                                                                     // Other obvious imports and namespace
#include <iterator>
#include <vector>
#include <difflib.h>																									// Difflib library (include)
using namespace seqan;

bool fileExists(const std::string & file){                                                                             // Function to determine if a file exists
    struct stat buf;
    return (stat(file.c_str(), &buf) == 0);
}

class individual_contig{                                                                                               // Class for storing an individual contig
private:                                                                                                                // Vector of all reads for that contig within private member
    std::string contigReferenceLabel;                                                                        			// Reference label from assembly (e.g. 17_1_1_7_2)
    std::vector<std::string> alignedReadsPresent;                                                           			// Vector of all reads found for this reference contig
    unsigned long numberOfContigReads;
public:
    void setReferenceLabel(std::string inputReferenceLabel){                                                     		// Set reference label function
        contigReferenceLabel = inputReferenceLabel;
        }
    void appendToReadVector(std::string alignedReadAppendTarget){                                                 		// Append an aligned read to this contig's vector
        alignedReadsPresent.push_back(alignedReadAppendTarget);
        }
    void updateReadNumber(unsigned long update){
        numberOfContigReads = update;
    }
    std::string getReferenceLabel(){                                                                        			// Return reference label
        return contigReferenceLabel;
        }
    std::vector<std::string> getAlignedReadsVector(){                                                        			// Return alignment vector
        return alignedReadsPresent;
        }
    unsigned long getNumberOfContigReads(){
        return numberOfContigReads;
    }
};

double similarityCalculator(std::string inSplit, std::string inMask){
	double differenceRatio = difflib::SequenceMatcher(inSplit, inMask).ratio();
	return differenceRatio;
}

std::vector<int> getRepeatTract(std::vector<std::string> inputTriplets, std::string currentMask){

	// Score the entire read against the current mask
	std::vector<std::pair<std::string, double>> currTract;
	for(std::string currentSplit : inputTriplets){
		double currentScore = similarityCalculator(currentSplit, currentMask);
		std::pair<std::string, double> comparisonResults (currentSplit, currentScore);
		currTract.push_back(comparisonResults);
	}

	//Anchors.. find beginning of repeat tract
	int regionStart = 0; int regionEnd = 0;
	for(int i = 0; i < currTract.size(); i++){
		try{
			if(currTract[i].second == 1.0){
				if(regionStart == 0){
					if(currTract[i+1].second == 1.0 and currTract[i+2].second == 1.0){
						regionStart = i;
					}
				}
			}
			if(currTract[i].second == 1.0){
				regionEnd = i+1;
			}
		}
		catch (Exception const & e){
			;
		}
	}

	// Try to fill a vector with a range(), start to end of this tract
	int targetFirstPassSize = (regionEnd-regionStart);
	std::vector<int> firstPassRange(targetFirstPassSize);
	std::iota(firstPassRange.begin(), firstPassRange.end(), regionStart);

	// Create iterators for current triplet, previous, next and next+1
	// Create our sub slice vector, to be filled with data from our iterators
	std::vector<std::pair<std::string,double>>::const_iterator prevItr, currItr, nextItr, followingItr;
	std::vector<std::pair<std::string, double>> subData;
	currItr = currTract.begin();
	while(currItr != currTract.end()){

		// Current iteration into a pair, append to our subslice vector
		std::pair<std::string, double> currSub = {currItr->first, currItr->second};
		subData.push_back(currSub);

		// Check our current iterator is != the beginning of our data structure
		// If so, we can decrement safely to get our previous iterator
		std::pair<std::string, double> prevSub;
		if(currItr != currTract.begin()){
			prevItr = currItr; --prevItr;
			prevSub = std::make_pair(prevItr->first, prevItr->second);
			subData.push_back(prevSub);
		}

		// Since we are in a while loop, we know currItr is != vector.end()
		// Check if incremented, iterator would not go out of bounds
		nextItr = currItr;
		std::pair<std::string, double> nextSub;
		if(++nextItr != currTract.end()){
			nextSub = std::make_pair(nextItr->first, nextItr->second);
			subData.push_back(nextSub);

			// Since we are in a while loop, we know nextItr is != vector.end()
			// Check if incremented, iterator would not go out of bounds
			followingItr = nextItr;
			std::pair<std::string, double> followingSub;
			if(++followingItr != currTract.end()){
				followingSub = std::make_pair(followingItr->first, followingItr->second);
				subData.push_back(followingSub);
			}
		}
		++currItr;
	}

	//Loop over rough range, remove items from subData which are not good matches for current mask
	for(int j = regionStart; j < regionEnd; j++){
		// If the current triplet score is not 1.0, check for downstream poor scores
		// i.e. we think the end of repeat region is here but was missed
		if (currTract[j].second != 1.0){
			int subScore = 0;
			for (auto item : subData){
				if (item.second == 1.0){
					subScore += 1;
				}
			}
			// Remove the entries where matches were false positives
			if (subScore != 3){
				for(int tripletIndex : firstPassRange){
					if(tripletIndex == j){
						int pos = std::find(firstPassRange.begin(), firstPassRange.end(), j) - firstPassRange.begin();
						firstPassRange.erase(firstPassRange.begin() + pos);
					}
				}
			}
		}
	}

	// Some downstream matches may still exist so..
	// remove anything outside of >1 streak in pass
	int diff = 0; int flaggedIDX = 0;
	std::vector<int> lengthOfFPRange(firstPassRange.size());
	std::iota(lengthOfFPRange.begin(), lengthOfFPRange.end(), 0);
	for(int k = 0; k < firstPassRange.size(); k++){
		try{
			diff = std::abs(firstPassRange[k+1]-firstPassRange[k]);
		}
		catch (Exception const & e){
			;
		}
		if(diff > 1 and flaggedIDX == 0){
			flaggedIDX = firstPassRange[k]+1;
		}
		for(auto index : firstPassRange){
			if(flaggedIDX != 0 and index > flaggedIDX){
				for(int tripletIndex : firstPassRange){
					if(tripletIndex == index){
						int pos = std::find(firstPassRange.begin(), firstPassRange.end(), index) - firstPassRange.begin();
						firstPassRange.erase(firstPassRange.begin() + pos);
					}
				}
			}
		}
	}

	return firstPassRange;
}

std::vector<int> getCCTRepeatTract(std::vector<std::string> inputTriplets, std::string currentMask, int anchor){

	// Get all triplets after the end of CCG tract (anchor)
	std::vector<std::pair<int, std::string>> postAnchor;
	for(int i = 0; i < inputTriplets.size(); i++){
		if(i > anchor){
			std::pair<int, std::string> postAnchorPair = std::make_pair(i, inputTriplets[i]);
			postAnchor.push_back(postAnchorPair);
		}
	}

	// If similarity matches the current mask, add index to tract
	std::vector<int> cctTract;
	for(auto item : postAnchor){
		if(similarityCalculator(currentMask, item.second) == 1.0){
			cctTract.push_back(item.first);
		}
	}

	// Loop over our CCT values.. remove wrongly added CCT outside tract
	// while not at the end of the vector, check current index and following index
	// calculate difference, flag if > 1 (CCT gapped regions never observed)
	int idxDiff = 0; int flaggedIDX = 0;
	std::vector<int>::const_iterator currItr, nextItr;
	currItr = cctTract.begin();
	while(currItr != cctTract.end()){
		int currIndex = currItr - cctTract.begin();
		nextItr = currItr;
		if(++nextItr != cctTract.end()){
			int nextIndex = nextItr - cctTract.begin();
			idxDiff = std::abs(*nextItr - *currItr);
			if(idxDiff > 1 and flaggedIDX == 0){
				flaggedIDX = *nextItr;
			}
		}
		std::advance(currItr,1);
	}

	// Delete the offending indices from our vector
	for(int i=0; i < cctTract.size(); i++) {
		if (flaggedIDX !=0 and cctTract[i] >= flaggedIDX){
			cctTract.erase(cctTract.begin() + i);
		}
	}

	return cctTract;
}

bool rotationCheck(std::string string1, std::string string2){

	//String sizes
	int size1 = string1.length();
	int size2 = string2.length();

	// If not the same, then invalid intv
	if(size1 != size2){
		return false;
	}

	// Double our known typical intervening to search for our matches
	std::string temp = string1 + string1;

	// If our input string is within double(string1) then a rotation of
	// the typical intervening sequence exists, return true
	if(temp.find(string2) != std::string::npos){
		return true;
	}else{
		return false;
	}
}

int scanReferenceReads(std::vector<individual_contig*> assemblyTargets){

	// Flags for ScaleHD
	bool fatalReadAllele = false;

	//
	// Check read count per ref
	if (assemblyTargets[0]->getNumberOfContigReads() < 200) {
		fatalReadAllele = true;
	}
	if (assemblyTargets[1]->getNumberOfContigReads() < 100) {
		fatalReadAllele = true;
	}
	if (assemblyTargets[2]->getNumberOfContigReads() < 50) {
		if ((unsigned) (assemblyTargets[2]->getNumberOfContigReads() - 45) <= (55 - 45)) { ;
		} else {
			fatalReadAllele = true;
		}
	}

	//
	// Iterate over top 3 aligned references in this assembly
	// Fetch the reads aligned to the current reference
	for (individual_contig *investigateContig : assemblyTargets) {
		std::cerr << "scanReferenceReads on : " << investigateContig->getReferenceLabel() << std::endl;

		// Counts of atypical/typical reads and other information
		int typicalCount = 0;
		int atypicalCount = 0;
		std::vector<std::string> intvPopulation;
		std::vector<std::string> fpFlanks;
		std::vector<std::string> tpFlanks;
		std::vector<int> refCAG;
		std::vector<int> refCCG;
		std::vector<int> refCCT;

		// For every sequenceRead in this contig, get the aligned sequence
		// Split into triplet sliding window list, remove any triplets that are < 3
		for (std::string currentRead : investigateContig->getAlignedReadsVector()) {

			std::cerr << "\n>>new read..." << std::endl;

			// Vector of all substrings of this seqread
			std::vector<std::string> sequenceWindows;
			for (unsigned i = 0; i < currentRead.length(); i += 3) {
				std::string tripletSlice = currentRead.substr(i, 3);
				if (tripletSlice.length() == 3) {
					sequenceWindows.push_back(currentRead.substr(i, 3));
				}
			}

			// Get repeat regions for CAG and CCG; similarity mask scores for current window
			// Any regions that are (idx > end_of_region) are truncated
			std::vector<std::string> cagMasks = {"CAG", "AGC", "GCA"};
			std::vector<std::string> ccgMasks = {"CCG", "CGC", "GCC"};
			std::vector<std::string> cctMasks = {"CCT", "CTC", "TCC"};

			// CAG/CCG masking -- sort vectors of tuples by repeat tract length
			std::vector<std::pair<std::string, std::vector<int>>> cagTracts;
			std::vector<std::pair<std::string, std::vector<int>>> ccgTracts;
			std::vector<std::pair<std::string, std::vector<int>>> cctTracts;
			for (std::string mask : cagMasks){
				std::vector<int> cagResults = getRepeatTract(sequenceWindows, mask);
				std::pair<std::string, std::vector<int>> cagToAppend = {mask, cagResults};
				cagTracts.push_back(cagToAppend);
			}
			for (std::string mask : ccgMasks){
				std::vector<int> ccgResults = getRepeatTract(sequenceWindows, mask);
				std::pair<std::string, std::vector<int>> ccgToAppend = {mask, ccgResults};
				ccgTracts.push_back(ccgToAppend);
			}

			// Sort our tract vectors
			std::sort(cagTracts.begin(), cagTracts.end(), [](auto &left, auto &right){
				return left.second.size() > right.second.size();
			});
			std::sort(ccgTracts.begin(), ccgTracts.end(), [](auto &left, auto &right){
				return left.second.size() > right.second.size();
			});

			// Select most common
			std::vector<int> cagTract = cagTracts.front().second;
			std::vector<int> ccgTract = ccgTracts.front().second;

			// CCT Masking/Intervening calculation
			// Sort CCT / select most common again
			std::string intvString; std::string fivePrimeFlankStr; std::string threePrimeFlankStr;
			for (std::string mask : cctMasks){
				int anchor = ccgTract.back();
				std::vector<int> indvCCTResult = getCCTRepeatTract(sequenceWindows, mask, anchor);
				std::pair<std::string, std::vector<int>> indvCCTPair;
				indvCCTPair = std::make_pair(mask, indvCCTResult);
				cctTracts.push_back(indvCCTPair);

			}
			std::sort(cctTracts.begin(), cctTracts.end(), [](auto &left, auto &right){
				return left.second.size() > right.second.size();
			});
			std::vector<int> cctTract = cctTracts.front().second;

			// Create range for intervening sequence (i.e. range between CAG.back and CCG.front)
			// Skip if malloc thrown (CAG[-1] idx > CCG[0] idx from bad quality read)
			int cagSuffix; int ccgPrefix; int intvSize;
			try{
				cagSuffix = cagTract.back() + 1;
				ccgPrefix = ccgTract.front();
				intvSize = ccgPrefix - cagSuffix;
			}catch(Exception const & e){;}
			std::vector<int> interveningRange(std::abs(intvSize));
			std::iota(interveningRange.begin(), interveningRange.end(), cagSuffix);


			// Create range for FP flank
			int cagPrefix;
			try{
				cagPrefix = cagTract.front() - 1;
			}catch(Exception const & e){;}
			std::vector<int> fivePrimeRange(abs(cagPrefix));
			std::iota(fivePrimeRange.begin(), fivePrimeRange.end(), 0);

			//Create a range for TP flank
			int cctSuffix; int readEnd; int threePrimeSize;
			try{
				if(cctTract.size() != 0){
					cctSuffix = cctTract.back() + 1;
					readEnd = sequenceWindows.size();
					threePrimeSize = readEnd - cctSuffix;
				}
			}catch(Exception const & e){;}
			std::vector<int> threePrimeRange(std::abs(threePrimeSize));
			std::iota(threePrimeRange.begin(), threePrimeRange.end(), cctSuffix);

			// Append lengths to contig-run
			refCAG.push_back(cagTract.size());
			refCCG.push_back(ccgTract.size());
			refCCT.push_back(cctTract.size());

			// Count FP, TP and intervening occurrences
			for(int i = 0; i < sequenceWindows.size(); i++){
				if(std::find(fivePrimeRange.begin(), fivePrimeRange.end(), i) != fivePrimeRange.end()){
					fivePrimeFlankStr += sequenceWindows[i];
				}
				if(std::find(threePrimeRange.begin(), threePrimeRange.end(), i) != threePrimeRange.end()){
					threePrimeFlankStr += sequenceWindows[i];
				}
				if(std::find(interveningRange.begin(), interveningRange.end(), i) != interveningRange.end()){
					intvString += sequenceWindows[i];
				}
			}

			// Check INTV structure
			std::string typicalIntv = "CAACAGCCGCCA";
			bool substringPresence = rotationCheck(typicalIntv, intvString);
			if(substringPresence){
				intvString = "CAACAGCCGCCA";
			}
			if (intvString != "CAACAGCCGCCA"){
				atypicalCount += 1;
			}else{
				typicalCount +=1;
			}

			// Append strings to relevant tract
			fpFlanks.push_back(fivePrimeFlankStr);
			tpFlanks.push_back(threePrimeFlankStr);
			intvPopulation.push_back(intvString);





		}
	}
	return 0;
}

int main(int argc, char const ** argv){                                                                                 // Main function
	const char* inputBamFilePath = "/Users/alastairm/Desktop/hdsptest/test_assembly.bam";                               // Path to test file.. TODO replace with stdin
	seqan::BamFileIn inputBamFileObject;

	if(!fileExists(toCString(inputBamFilePath))){                                                                       // Check that the file exists..
		std::cerr << "Input file: " << inputBamFilePath << " does not exist." << std::endl;                             // File doesn't exist, quit program
		return 1;
	}
	else{
		std::cerr << "Input file: " << inputBamFilePath << " exists OK!" << std::endl;                                  // File exists, check if openable by Seqan
		if (!seqan::open(inputBamFileObject, toCString(inputBamFilePath))){
			std::cerr << "Error: Could not open " << inputBamFilePath << std::endl;
			std::cerr << "Do you have ZLIB in your include folder? " << std::endl;                                      // ZLIB required to be present for C++ pre-processor
		}                                                                                                               // for G++, "-D SEQAN_HAS_ZLIB" passed in CLI
		else{
			std::cerr << "File: " << inputBamFilePath << " opened OK!\n" << std::endl;
		}
	}

	std::vector<individual_contig*> assemblyContigObjects;                                                              // Vector of all contig record objects we will create/fill
	std::vector<seqan::CharString> observedContigs;                                                                     // Vector for adhering to only novel contigs

	try{                                                                                                                // Try scrape data from our sequence assembly
		seqan::BamHeader alignmentHeader;                                                                               // Create header object
		seqan::readHeader(alignmentHeader, inputBamFileObject);                                                         // Read header to preserve starting indice
		while(!atEnd(inputBamFileObject)){                                                                              // Loop over open sequence assembly data
			seqan::BamAlignmentRecord alignmentRecord;                                                                  // Create object for current alignment record (sequence read)
			seqan::readRecord(alignmentRecord, inputBamFileObject);                                                     // Fill BamAlignmentRecord with data via readRecord()
			seqan::CharString initContigName = getContigName(alignmentRecord, inputBamFileObject);                          // Get contig name for current alignment record
			seqan::CharString initSequenceRead = alignmentRecord.seq;                                                       // Get DNA sequence for current alignment record
			//Convert to string
			char * intermediateContigName = toCString(initContigName);
			std::string contigName = std::string(intermediateContigName);
			char * intermediateSequenceRead = toCString(initSequenceRead);
			std::string sequenceRead = std::string(intermediateSequenceRead);

			if(std::find(observedContigs.begin(), observedContigs.end(), contigName) != observedContigs.end()){         // if contigName is observed, loop to find contig object
				for(individual_contig* currentContig : assemblyContigObjects){                                          // then when the correct contig is found, we append that object's seq-read vector
					if(contigName == currentContig->getReferenceLabel()){
						currentContig->appendToReadVector(sequenceRead);                                                // Read appending to object's vector
						currentContig->updateReadNumber(currentContig->getAlignedReadsVector().size());
					}
				}
			}
			else{                                                                                                       // The current contigName was not found in observedContigs
				individual_contig* newObservation = new individual_contig;                                              // Create new object ON THE HEAP, fill with data
				newObservation->setReferenceLabel(contigName);                                                          // Contig name
				newObservation->appendToReadVector(sequenceRead);                                                       // Add read to object's sequenceRead vector
				assemblyContigObjects.push_back(newObservation);                                                        // Append our job-wide vector with this new object
				observedContigs.push_back(contigName);                                                                  // Update our vector of observed contigs
				newObservation->updateReadNumber(newObservation->getAlignedReadsVector().size());
			}
		}
	}
	catch(Exception const & e){                                                                                         // We caught an exception from any of the above
		std::cout << "Caught exception: " << e.what() << std::endl;                                                     // Inform user
		return 1;                                                                                                       // Quit status 1 (error)
	}
	std::sort(assemblyContigObjects.begin(), assemblyContigObjects.end(), [](const auto& lhs, const auto& rhs){         // Sort our vector of contig objects by number of reads
		return lhs->getNumberOfContigReads() > rhs->getNumberOfContigReads();
	});
	assemblyContigObjects.resize(3);                                                                                    // Only keep the first 5 contigs

	// Here we try to replicate ScaleHD python functions
	// C++ assemblyContigObjects = Py assembly_targets
	scanReferenceReads(assemblyContigObjects);

	std::cerr << "Finished! Exiting." << std::endl;                                                                     // Done processing on this sequence assembly!
	return 0;                                                                                                           // Quit status 0 (success)
}



