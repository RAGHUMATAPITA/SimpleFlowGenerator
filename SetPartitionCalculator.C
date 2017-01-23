//-------------------------------------------------------
// Code to compute all of the partitions of a given set.
//
// J. Orjuela-Koop
//
// Jan. 2017
//-------------------------------------------------------

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>

#include <sys/time.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TMath.h"
#include "TComplex.h"

using namespace std;

//-------------------------------------
// Structures and Classes
//-------------------------------------

struct ComplexAverage {
	int numElements = 0;

	TComplex value = TComplex(0, 0);

	void addElement(TComplex newElement)
	{
		value = value + newElement;
		numElements = numElements + 1;
	}

	TComplex getValue()
	{
		return value;
	}

	float getAverageReal()
	{
		return value.Re() / numElements;
	}

	int getElements()
	{
		return numElements;
	}
};

struct PSet {
	vector<int> elements;

	void addElement(int newElement)
	{
		//First, check that the new element is not already in the set
		for (int i = 0; i < elements.size(); i++)
		{
			if (newElement == elements[i]) return;
		}

		elements.push_back(newElement);
	}

	int setSize()
	{
		return elements.size();
	}

	int getElement(int i)
	{
		return elements[i];
	}

	int getSetIndex()
	{
		int setIndex = 0;

		for (int i = 0; i < elements.size(); i++)
		{
			setIndex += pow(10, elements.size() - i - 1) * elements[i];
		}

		return setIndex;
	}

	void printSet()
	{
		cout << " [ ";
		for (int i = 0; i < elements.size(); i++)
		{
			cout << elements[i] << " ";
		}
		cout << "] ";
	}
};

struct SetList {
	vector<PSet> list;

	void addSet(PSet s)
	{
		list.push_back(s);
	}

	int listSize()
	{
		return list.size();
	}

	PSet getSet(int i)
	{
		return list[i];
	}

	void clear()
	{
		list.clear();
	}

	void printList()
	{
		cout << "{";

		for (int i = 0; i < list.size(); i++)
		{
			PSet s = list[i];
			s.printSet();
		}

		cout << "}" << endl;
	}
};

//-------------------------------------
// Variables
//-------------------------------------

//Number of particles per event
const int MULT = 6;

//Set verbosity on or off
bool verbosity = false;

//TProfile to take the average of the cumulants over many events
TProfile *avgCumulant;

TProfile *f1234;
TProfile *f13;
TProfile *f24;
TProfile *f14;
TProfile *f23;

//Vector to store the partitions of a given set, where each partition is a SetList object
vector<SetList> setPartitions;

//Vector to store all subsets of a given set
vector<PSet> setSubsets;

//Map to store the histograms for the distributions of the products of all possible particle combinations (subsets)
std::map <int, TProfile*> mapParticleComboHisto;

//Map to store the complex averages of the products of all possible particle combinations (subsets)
std::map <int, ComplexAverage> mapParticleComboComplex;

//Number of particles to use in cumulant analysis
int kOrder;

//Number of recursive steps to compute all the partitions of the set
int recursiveStep = 0;

//-------------------------------------
// Functions
//-------------------------------------

/*
 *
 */
void initializeHistograms()
{
	for (int i = 0; i < setSubsets.size(); i++)
	{
		//Get the integer corresponding to the set elements
		//E.g., "123" = 100*1 + 10*2 + 1*3
		PSet s = setSubsets[i];
		int setIndex = s.getSetIndex();

		//Create new histogram and add to map
		mapParticleComboHisto[setIndex] = new TProfile(Form("%i", setIndex), Form("%i", setIndex), 1, 0, 100);
	}

	for (int i = 0; i < setSubsets.size(); i++)
	{
		PSet s = setSubsets[i];
		int setIndex = s.getSetIndex();

		ComplexAverage ca;
		mapParticleComboComplex[setIndex] = ca;
	}

	//Initialize TProfile to take the average cumulant over many events
	avgCumulant = new TProfile("avgCumulant", "avgCumulant", 1, -1000, 1000);
	f1234 = new TProfile("f1234", "f1234", 1, -1000, 1000);
	f13 = new TProfile("f13", "f13", 1, -1000, 1000);
	f24 = new TProfile("f24", "f24", 1, -1000, 1000);
	f14 = new TProfile("f14", "f14", 1, -1000, 1000);
	f23 = new TProfile("f23", "f23", 1, -1000, 1000);
}

/*
 * Print all partitions found for a given set
 */
void printSetPartitions()
{
	cout << endl << "--> Found " << setPartitions.size() << " Partitions of the Set { ";

	for (int i = 0; i < kOrder; i++)
	{
		if (i == kOrder - 1)
		{
			cout << i + 1 << " }";
		}
		else
		{
			cout << i + 1 << ", ";
		}
	}

	cout << endl << endl;

	for (int i = 0; i < setPartitions.size(); i++)
	{
		SetList partition = setPartitions[i];

		partition.printList();
	}
}


/*
 * Print all subsets found for a given set
 */
void printSetSubsets()
{
	cout << endl << "--> Found " << setSubsets.size() << " Subsets of the Set { ";

	for (int i = 0; i < kOrder; i++)
	{
		if (i == kOrder - 1)
		{
			cout << i + 1 << " }";
		}
		else
		{
			cout << i + 1 << ", ";
		}
	}

	cout << endl << endl;

	for (int i = 0; i < setSubsets.size(); i++)
	{
		PSet subset = setSubsets[i];

		subset.printSet();
	}

	cout << endl;
}


/*
 * Determine the equality of two given sets.
 * Sets must be the same size and contain the same elements
 * regardless of their order.
 */
bool setEquality(PSet s1, PSet s2)
{
	//First check that the sets are the same size
	if (s1.setSize() != s2.setSize()) return false;

	for (int i = 0; i < s1.setSize(); i++)
	{
		int found = 0;

		for (int j = 0; j < s2.setSize(); j++)
		{
			int element1 = s1.getElement(i);
			int element2 = s2.getElement(j);

			//If an element from set1 is found in set2, remove it from set2
			if (element1 == element2) found++;
		}

		if (found != 1) return false;
	}

	return true;
}


/*
 * Determine the equality of two given lists.
 * List must be the same size and contain the same sets
 * regardless of their order.
 */
bool listEquality(SetList l1, SetList l2)
{
	//If the list sizes are different, return false
	if (l1.listSize() != l2.listSize()) return false;

	for (int i = 0; i < l1.listSize(); i++)
	{
		int found = 0;

		for (int j = 0; j < l2.listSize(); j++)
		{
			PSet s1 = l1.getSet(i);
			PSet s2 = l2.getSet(j);

			if (setEquality(s1, s2)) found++;
		}

		if (found != 1) return false;
	}

	return true;
}


/*
 * Carry out the union of two given sets,
 * returning a new set
 */
PSet setUnion(PSet s1, PSet s2)
{
	PSet sOut;

	for (int i = 0; i < s1.setSize(); i++)
	{
		sOut.addElement(s1.getElement(i));
	}

	for (int i = 0; i < s2.setSize(); i++)
	{
		sOut.addElement(s2.getElement(i));
	}

	return sOut;
}


/*
 * Insert a new partition into the list of partitions
 * only if it has not been previously encountered
 */
void insertNewPartition(SetList newPartition)
{
	int found = 0;

	for (int i = 0; i < setPartitions.size(); i++)
	{
		SetList oldPartition = setPartitions[i];

		if (listEquality(newPartition, oldPartition)) found++;
	}

	if (found > 0) return;

	setPartitions.push_back(newPartition);
}


/*
 * Insert a new set into the vector of subsets
 * only if it hasn't been encountered before
 */
void insertNewSubset(PSet subset)
{
	for (int i = 0; i < setSubsets.size(); i++)
	{
		PSet oldSet = setSubsets[i];

		if (setEquality(subset, oldSet)) return;
	}

	setSubsets.push_back(subset);
}


/*
 * Recursive algorithm that actually computes the partitons
 * TODO: Document
 */
void findPartitions(SetList sList)
{
	recursiveStep++;

	//Add the given list as starting point
	insertNewPartition(sList);

	//Create new empty list to store the result of merging
	SetList mergedList;

	//Condition to break out of recursion
	//if (sList.listSize() == 1) return;

	for (int i = 0; i < sList.listSize(); i++)
	{
		PSet s1 = sList.getSet(i);
		insertNewSubset(s1);

		for (int j = 0; j < sList.listSize(); j++)
		{
			if (i == j) continue;

			PSet sAux = setUnion(s1, sList.getSet(j));
			mergedList.addSet(sAux);
			insertNewSubset(sAux);

			for (int k = 0; k < sList.listSize(); k++)
			{
				if (j == k || i == k) continue;

				mergedList.addSet(sList.getSet(k));
			}

			//Save newly found partition
			insertNewPartition(mergedList);

			//Recursive call
			findPartitions(mergedList);

			mergedList.clear();
		}
	}
}


/*
 * Implementation of the following formula for the cumulant using the partitions
 *
 * k(x_1, x_2, ..., x_k) = \sum_{\pi} (|\pi|-1)! (-1)^{|\pi|-1} \prod_{B \in \pi} \langle \prod_{i \in B} x_i \rangle
 */
void computeKCumulant()
{
	float cumulant = 0;

	//Total number of partitions of the set
	int numPartitions = setPartitions.size();

	for (int p = 0; p < numPartitions; p++)
	{
		//Get the partition (i.e., list)
		SetList partition = setPartitions[p];

		//Number of parts (i.e., sets) in the partion
		int numParts = partition.listSize();

		float factor1 = TMath::Factorial(TMath::Abs(numParts) - 1);
		float factor2 = TMath::Power(-1.0, TMath::Abs(numParts) - 1);

		if (verbosity)
		{
			cout << "--> Partition " << p << ": ";
			partition.printList();
			cout << endl;
			cout << "     Factor1 = " << factor1 << endl;
			cout << "     Factor2 = " << factor2 << endl;
		}

		//Now loop over all blocks (i.e., sets) in the current partition
		float factor3 = 1.0;

		for (int j = 0; j < numParts; j++)
		{
			//Get the block and its index
			PSet block     = partition.getSet(j);
			int blockIndex = block.getSetIndex();

			factor3 = factor3 * mapParticleComboHisto[blockIndex]->GetMean();
		}

		if (verbosity) cout << "     Factor3 = " << factor3 << endl << endl;


		cumulant = cumulant + factor1 * factor2 * factor3;
	}

	if (verbosity) cout << "CUMULANT = " << cumulant << endl;
}


/*
 * Implementation of the following formula for the cumulant using the partitions
 *
 * k(x_1, x_2, ..., x_k) = \sum_{\pi} (|\pi|-1)! (-1)^{|\pi|-1} \prod_{B \in \pi} \langle \prod_{i \in B} x_i \rangle
 */
void computeComplexKCumulant()
{
	float cumulant = 0;

	//Total number of partitions of the set
	int numPartitions = setPartitions.size();

	for (int p = 0; p < numPartitions; p++)
	{
		//Get the partition (i.e., list)
		SetList partition = setPartitions[p];

		//Number of parts (i.e., sets) in the partion
		int numParts = partition.listSize();

		float factor1 = TMath::Factorial(TMath::Abs(numParts) - 1);
		float factor2 = TMath::Power(-1.0, TMath::Abs(numParts) - 1);

		if (verbosity)
		{
			cout << "--> Partition " << p << ": ";
			partition.printList();
			cout << endl;
			cout << "     Factor1 = " << factor1 << endl;
			cout << "     Factor2 = " << factor2 << endl;
		}

		//Now loop over all blocks (i.e., sets) in the current partition
		float factor3 = 1.0;

		for (int j = 0; j < numParts; j++)
		{
			//Get the block and its index
			PSet block     = partition.getSet(j);
			int blockIndex = block.getSetIndex();

			factor3 = factor3 * mapParticleComboComplex[blockIndex].getAverageReal();
		}

		if (verbosity) cout << "     Factor3 = " << factor3 << endl << endl;


		cumulant = cumulant + factor1 * factor2 * factor3;
	}

	if (verbosity) cout << "CUMULANT = " << cumulant << endl;

	avgCumulant->Fill(cumulant, 1);
}


/*
 *
 */
void runEvent2PC(std::vector<float> phi)
{
	TComplex phipair[2];

	for (int i = 0; i < MULT; i++)
	{
		phipair[0] = TComplex(TMath::Cos(2 * phi[i]), TMath::Sin(2 * phi[i]));

		for (int j = 0; j < MULT; j++)
		{
			phipair[1] = TComplex(TMath::Cos(-2 * phi[j]), TMath::Sin(-2 * phi[j]));

			if (i == j) continue;

			for (int k = 0; k < setSubsets.size(); k++)
			{
				PSet s = setSubsets[k];
				int index = s.getSetIndex();
				TComplex content = TComplex(1, 0);

				//Get the constituent digits of the index
				int num = index;
				while (num > 0)
				{
					content = content * phipair[num % 10 - 1];
					num = num / 10.0;
				}

				mapParticleComboComplex[index].addElement(content);
			}
		}
	}

	computeComplexKCumulant();
}

/*
 *
 */
void runEvent4PC(std::vector<float> phi)
{
	TComplex phigroup[4];

	for (int i = 0; i < MULT; i++)
	{
		phigroup[0] = TComplex(TMath::Cos(4 * phi[i]), TMath::Sin(4 * phi[i]));

		for (int j = 0; j < MULT; j++)
		{
			phigroup[1] = TComplex(TMath::Cos(4 * phi[j]), TMath::Sin(4 * phi[j]));

			for (int k = 0; k < MULT; k++)
			{
				phigroup[2] = TComplex(TMath::Cos(-4 * phi[k]), TMath::Sin(-4 * phi[k]));

				for (int l = 0; l < MULT; l++)
				{
					phigroup[3] = TComplex(TMath::Cos(-4 * phi[l]), TMath::Sin(-4 * phi[l]));

					if (i == j || i == k || i == l || j == k || j == l || k == l) continue;

					//cout << "PARTICLE PHIS  = " << phigroup[0] << "  " << phigroup[1] << "  " << phigroup[2] << "  " << phigroup[3] << endl;

					for (int m = 0; m < setSubsets.size(); m++)
					{
						PSet s = setSubsets[m];
						int index = s.getSetIndex();
						TComplex content = TComplex(1, 0);

						//Get the constituent digits of the index
						int num = index;
						while (num > 0)
						{
							content = content * phigroup[num % 10 - 1];
							num = num / 10.0;
						}

						mapParticleComboComplex[index].addElement(content);

						if (index == 13)
						{
							f13->Fill(content.Re(), 1);
						}

						if (index == 24)
						{
							f24->Fill(content.Re(), 1);
						}

						if (index == 14)
						{
							f14->Fill(content.Re(), 1);
						}

						if (index == 23)
						{
							f23->Fill(content.Re(), 1);
						}

						if (index == 1234)
						{
							f1234->Fill(content.Re(), 1);
						}
					}
				}
			}
		}
	}

	computeComplexKCumulant();
}


/*
 *
 */
void runEvent6PC(std::vector<float> phi)
{
	TComplex phigroup[6];

	for (int i = 0; i < MULT; i++)
	{
		phigroup[0] = TComplex(TMath::Cos(6 * phi[i]), TMath::Sin(6 * phi[i]));

		for (int j = 0; j < MULT; j++)
		{
			phigroup[1] = TComplex(TMath::Cos(6 * phi[j]), TMath::Sin(6 * phi[j]));

			for (int k = 0; k < MULT; k++)
			{
				phigroup[2] = TComplex(TMath::Cos(6 * phi[k]), TMath::Sin(6 * phi[k]));

				for (int l = 0; l < MULT; l++)
				{
					phigroup[3] = TComplex(TMath::Cos(-6 * phi[l]), TMath::Sin(-6 * phi[l]));

					for (int m = 0; m < MULT; m++)
					{
						phigroup[4] = TComplex(TMath::Cos(-6 * phi[m]), TMath::Sin(-6 * phi[m]));

						for (int n = 0; n < MULT; n++)
						{
							phigroup[5] = TComplex(TMath::Cos(-6 * phi[n]), TMath::Sin(-6 * phi[n]));

							if (i == j || i == k || i == l || i == m || i == n || j == k || j == l || j == m || j == n || k == l || k == m || k == n || l == m || l == n || m == n) continue;

							for (int t = 0; t < setSubsets.size(); t++)
							{
								PSet s = setSubsets[t];
								int index = s.getSetIndex();
								TComplex content = TComplex(1, 0);

								//Get the constituent digits of the index
								int num = index;
								while (num > 0)
								{
									content = content * phigroup[num % 10 - 1];
									num = num / 10.0;
								}

								mapParticleComboComplex[index].addElement(content);
							}
						}
					}
				}
			}
		}
	}

	computeComplexKCumulant();
}


/*
 * Load the particles from synthetic v2 generator
 */
void loadSyntheticParticles()
{
	//Read in the tree from the particle generator and extract the phi branch
	//There is an array of phi values for every event (i.e., every entry in the tree)
	TFile file("simpletree_10k.root");
	TTreeReader reader("simpletree", &file);
	TTreeReaderArray<float> d_phi(reader, "phi");

	//Run the cumulant calculation on an event-by-event basis
	int nEvent = 0;
	std::vector<float> eventPhiValues;
	while (reader.Next())
	{
		//if (nEvent > 0) break;

		if (nEvent % 100 == 0) cout << "----> Processing event " << nEvent << endl;

		for (int i = 0; i < MULT; i++)
		{
			eventPhiValues.push_back(d_phi[i]);
		}

		if (kOrder == 2)
		{
			runEvent2PC(eventPhiValues);
		}
		else if (kOrder == 4)
		{
			runEvent4PC(eventPhiValues);
		}
		else if(kOrder == 6)
		{
			runEvent6PC(eventPhiValues);
		}

		eventPhiValues.clear();
		nEvent++;
	}
}


void loadParticles()
{
	float phi1[5] = {1, 2, 3, 4, 5};

	float phipair[2];

	for (int i = 0; i < 5; i++)
	{
		phipair[0] = phi1[i];

		for (int j = 0; j < 5; j++)
		{
			phipair[1] = phi1[j];

			if (i == j) continue;

			for (int k = 0; k < setSubsets.size(); k++)
			{
				PSet s = setSubsets[k];
				int index = s.getSetIndex();
				float content = 1;

				//Get the constituent digits of the index
				int num = index;
				while (num > 0)
				{
					content = content * phipair[num % 10 - 1];
					num = num / 10;
				}

				mapParticleComboHisto[index]->Fill(content, 1);
			}
		}
	}
}


/*
 * Load the partitions for {1,2,3,4,5,6} from file instead
 * of computing them recursively every time the program is run.
 * This saves a lot of time!
 */
void load6PCPartitions()
{

}

/*
 * Entry point for the library
 * The order of k-particle cumulants is given as parameter
 */
void SetPartitionCalculator(int k)
{
	kOrder = k;

	//Initialize a set for the number of desired cumulants
	vector<PSet> sets;

	for (int i = 0; i < kOrder; i++)
	{
		PSet p;
		p.addElement(i + 1);
		sets.push_back(p);
	}

	SetList l1;

	for (int i = 0; i < kOrder; i++)
	{
		l1.addSet(sets[i]);
	}

	//Find partitions and subsets for the case of k-particle cumulants
	cout << "Finding partitions for " << k << " elements" << endl << endl;
	findPartitions(l1);

	//Once found, define histograms to store the distribution of the product of particles
	//for all possible combination of particles.
	//The possible combinations are defined by the subsets themselves
	cout << "Initializing complex variables for bookkeeping" << endl << endl;
	initializeHistograms();

	//Add particles for cumulant calculation
	loadSyntheticParticles();

	//Compute cumulant
	//cout << "Sqrt(< k_" << kOrder << " >) = " << pow(-1 * avgCumulant->GetMean(), 0.25) << endl << endl;
	cout << "Sqrt(< k_" << kOrder << " >) = " << pow(avgCumulant->GetMean(), 0.5) << endl << endl;
	if(kOrder == 6)
	{
		cout << "k_6 = " << avgCumulant->GetMean() << endl;
	}

	//Print partitions and subsets, for diagnostic purposes
	printSetPartitions();
	printSetSubsets();
}