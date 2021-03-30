//COP 3530
//Project 3: Asteroid Mission Planner
#include<iostream>
using namespace std;

const int SCIENCE_MISSION = 1;
const int MINING_MISSION = 2;
const int RANDOM_MISSION = 3;

const int QUICKSORT_ALGORITHM = 1;
const int AVL_TREE_ALGORITHM = 2;
const int HEAP_ALGORITHM = 3;

void runQuickSort(int missionType, int resultCount){
	cout << "Quicksort not implemented yet" << endl;
}

void runTreeAlgorithm(int missionType, int resultCount){
	cout << "Tree not implemented yet" << endl;
}

void runHeapAlgorithm(int missionType, int resultCount){
	cout << "Heap not implemented yet" << endl;
}

int main(){
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
		runHeapAlgorithm(missionTypeChoice, resultCount);
	}
	
}