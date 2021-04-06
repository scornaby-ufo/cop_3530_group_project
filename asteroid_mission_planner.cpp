//COP 3530
//Project 3: Asteroid Mission Planner
#include<iostream>
#include<vector>
#include <random>
using namespace std;

const int SCIENCE_MISSION = 1;
const int MINING_MISSION = 2;
const int RANDOM_MISSION = 3;

const int QUICKSORT_ALGORITHM = 1;
const int AVL_TREE_ALGORITHM = 2;
const int HEAP_ALGORITHM = 3;

class Asteroid
{
	//Add useful entries from data file here, like name, size, etc...
	public:
		int id;
	
};

class AsteroidMissionRating
{
	public:
		Asteroid asteroid;
		double missionRating;
};

void runQuickSort(int missionType, int resultCount){
	cout << "Quicksort not implemented yet" << endl;
}

void runTreeAlgorithm(int missionType, int resultCount){
	cout << "Tree not implemented yet" << endl;
}

void runHeapAlgorithm(vector<Asteroid>& asteroids, int missionType, int resultCount){
	cout << "Heap not implemented yet" << endl;
	
	default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0,1.0);
	
	AsteroidMissionRating asteroidHeap[resultCount];
	int heapSize = 0;
	
	for(Asteroid asteroid: asteroids){
		AsteroidMissionRating newRatedAsteroid;
		newRatedAsteroid.asteroid = asteroid;
		if(missionType == RANDOM_MISSION){
			newRatedAsteroid.missionRating = distribution(generator);
		}
		
		//Sort into heap
		
		//If the heap is not yet full, insert and bubble up
		if(heapSize < resultCount){
			asteroidHeap[heapSize]=newRatedAsteroid;

			int currentHeapPoint = heapSize;
			while(true){
				//We have reached the top of the heap, time to exit the loop
				if(currentHeapPoint == 0){
					heapSize++;
					break;
				}
				int heapComparisonPoint = currentHeapPoint % 2 == 1 ? (currentHeapPoint - 1)/2 : (currentHeapPoint - 2) /2;
				if(asteroidHeap[heapComparisonPoint].missionRating > newRatedAsteroid.missionRating){
					AsteroidMissionRating temp = asteroidHeap[heapComparisonPoint];
					asteroidHeap[heapComparisonPoint] = newRatedAsteroid;
					asteroidHeap[currentHeapPoint] = temp;
					currentHeapPoint = heapComparisonPoint;
				}
				else{
					heapSize++;
					break;
				}
			}
		}
		else{ //If the heap is full, see if the new value is better than the current min
			if(newRatedAsteroid.missionRating > asteroidHeap[0].missionRating){
				//Swap the value into the "min" slot and bubble it down as needed
				asteroidHeap[0]=newRatedAsteroid;
				int heapComparisonPoint = 0;
				while(heapComparisonPoint * 2 + 2 <= heapSize){
					AsteroidMissionRating leftChild = asteroidHeap[heapComparisonPoint * 2 + 1];
					AsteroidMissionRating rightChild = asteroidHeap[heapComparisonPoint * 2 + 2];
					if(newRatedAsteroid.missionRating < leftChild.missionRating && newRatedAsteroid.missionRating < rightChild.missionRating){
						break;
					}
					else if(leftChild.missionRating < rightChild.missionRating){
						asteroidHeap[heapComparisonPoint] = leftChild;
						asteroidHeap[heapComparisonPoint * 2 + 1] = newRatedAsteroid;
						heapComparisonPoint = heapComparisonPoint * 2 + 1;
					}
					else{
						asteroidHeap[heapComparisonPoint] = rightChild;
						asteroidHeap[heapComparisonPoint * 2 + 2] = newRatedAsteroid;
						heapComparisonPoint = heapComparisonPoint * 2 + 2;
					}
				}
			}
		}
		
		
	}
	
	
	//Now extract and reverse the min heap to get the max N results
	AsteroidMissionRating sortedResults[resultCount];
	int resultsSortedSoFar = 0;
	
	
	while(heapSize > 0){
				//Get minimum match from front of heap, swap in last value, resort heap
				sortedResults[resultCount - 1 - resultsSortedSoFar] = asteroidHeap[0];
				resultsSortedSoFar++;
				
				AsteroidMissionRating newRoot = asteroidHeap[heapSize];
				asteroidHeap[0]=newRoot;
				heapSize--;
				int heapComparisonPoint = 0;
				while(heapComparisonPoint * 2 + 2 <= heapSize){
					AsteroidMissionRating leftChild = asteroidHeap[heapComparisonPoint * 2 + 1];
					AsteroidMissionRating rightChild = asteroidHeap[heapComparisonPoint * 2 + 2];
					if(newRoot.missionRating < leftChild.missionRating && newRoot.missionRating < rightChild.missionRating){
						break;
					}
					else if(leftChild.missionRating < rightChild.missionRating){
						asteroidHeap[heapComparisonPoint] = leftChild;
						asteroidHeap[heapComparisonPoint * 2 + 1] = newRoot;
						heapComparisonPoint = heapComparisonPoint * 2 + 1;
					}
					else{
						asteroidHeap[heapComparisonPoint] = rightChild;
						asteroidHeap[heapComparisonPoint * 2 + 2] = newRoot;
						heapComparisonPoint = heapComparisonPoint * 2 + 2;
					}
				}
	}
	
	cout << "Sorted results" << endl;
	for(int i = 0; i < resultCount; i++){
		cout << sortedResults[i].asteroid.id << ": " << sortedResults[i].missionRating << endl;
	}
	

	
}

int main(){
	cout << "Loading asteroid data..." << endl;
	
	//Asteroid import code goes here
	vector<Asteroid> asteroids;
	for(int i = 0; i < 500000; i++)
	{
		Asteroid sampleAsteroid;
		sampleAsteroid.id = i;
		asteroids.push_back(sampleAsteroid);
	}
	
	
	cout << "Choose your mission criteria" << endl;
	cout << "\t1) Science Mission: Asteroids with near neighbors" << endl;
	cout << "\t2) Mining Mission: Large, metal rich easy to land on asteroids" << endl;
	cout << "\t3) Random Mission: Assign asteroids a random mission value" << endl;
	
	int missionTypeChoice;
	cin >> missionTypeChoice;
	
	cout << "How many results do you want? [1 to 600,000]" << endl;
	
	int resultCount;
	cin >> resultCount;
	
	cout << "Choose algorithm for sorting results" << endl;
	cout << "\t1)Quicksort" << endl;
	cout << "\t2)AVL Tree" << endl;
	cout << "\t3)Heap" << endl;
	
	int algorithmTypeChoice;
	cin >> algorithmTypeChoice;
	
	if(algorithmTypeChoice == QUICKSORT_ALGORITHM){
		runQuickSort(missionTypeChoice, resultCount);
	}
	else if(algorithmTypeChoice == AVL_TREE_ALGORITHM){
		runTreeAlgorithm(missionTypeChoice, resultCount);
	}
	else if(algorithmTypeChoice == HEAP_ALGORITHM){
		runHeapAlgorithm(asteroids, missionTypeChoice, resultCount);
	}
	
}