//COP 3530
//Project 3: Asteroid Mission Planner
#include<iostream>
#include<fstream>
#include<vector>
#include <random>
#include <chrono>
#include<string>
#include<cmath>
#include<iomanip>
using namespace std;
using namespace std::chrono;

#define PI 3.14159265
#define DAYS_PER_YEAR 365.25636;

const int SCIENCE_MISSION = 1;
const int MINING_MISSION = 2;
const int RANDOM_MISSION = 3;

const int QUICKSORT_ALGORITHM = 1;
const int AVL_TREE_ALGORITHM = 2;
const int HEAP_ALGORITHM = 3;

const int APRIL_2021_MJD = 59320; //Number of days since November 17, 1858. Used as a base time in astronomy for historical purposes
const double MAX_NEIGHBOR_DISTANCE = .01; //How many astronomical units apart can two asteroids be while still being considered "near neighbors" for purpose of a mission?
const double MAX_NEIGHBOR_DISTANCE_SQUARED = pow(MAX_NEIGHBOR_DISTANCE, 2); //Comparing distances can be made slightly faster in some cases by comparing squares and thus avoiding multiple square roots

class Asteroid
{
	//Add useful entries from data file here, like name, size, etc...
	public:
		int id;
		string name; //Official name of the asteroid
		double a; //semi major axis of the asteroid's orbit
		double e; //eccentricity of the asteroid's orbit
		double i; //inclination of the asteroid's orbit
		double om; //longitude of ascending node
		double w; //argument of perihelion
		double ma; //mean anomaly of the asteroid's orbit
		double magnitude; //Observed visual size of the asteroid, known more often than actual size
		int epoch_mjd; //Modified julian date for the reference time the asteroids orbital data
		double x; //x coordinate in sun centered coordinate system
		double y; //y coordinate in sun centered coordinate system
		double z; //z coordinate in sun centered coordinate system
		double missionRating;

		void calculatePosition() {
			double M = currentM(ma, a, epoch_mjd - APRIL_2021_MJD);
			double E = eccentricAnomaly(M, e);
			double v = approximateTrueAnomaly(M, e);
			double r = helioDistance(a, e, v);
			double theta = v + w;
			double omRad = om * PI / 180;
			double thetaRad = theta * PI / 180;
			double iRad = i * PI / 180;

			x = r * (cos(omRad) * cos(thetaRad) - sin(omRad) * sin(thetaRad) * cos(iRad));
			y = r * (sin(omRad) * cos(thetaRad) + cos(omRad) * sin(thetaRad) * cos(iRad));
			z = r * sin(thetaRad) * sin(iRad);

		}

	private:
		//Calculate the current anomaly of an asteroid based on how many days it has been since it's original reference
		double currentM(double Mref, double a, int days) {
			return Mref + 2 * PI / sqrt(pow(a, 3)) * (double)days/ DAYS_PER_YEAR;
		}
		//Use newton's method to calculate a very good approximation of the eccentric anomaly of the asteroid's orbit
		double eccentricAnomaly(double M, double e) {
			double E = M;
			double newE = (M - e * (E * cos(E * PI / 180) - sin(E * PI / 180))) / (1 - e * cos(E * PI / 180));
			int count = 0;
			while (count < 1000 && abs(newE - E) > .0001) {
				E = newE;
				newE = (M - e * (E * cos(E * PI / 180) - sin(E * PI / 180))) / (1 - e * cos(E * PI / 180));
				count++;
				if (count == 1000) {
					cout << "WARNING: Non converging eccentric anomaly. Asteroid orbit prediction for " << name << " may be of low quality" << endl;
				}
			}

			return newE;
		}

		//Use a fourier transform found on Wikipedia to find a close enough approximation of the orbit's true anomaly
		double approximateTrueAnomaly(double M, double e) {
			double radM = M * PI / 180;
			return M + (2 * e - 1 / 4 * pow(e, 3)) * sin(radM)
				+ 5 / 4 * pow(e, 2) * sin(2 * radM)
				+ 13 / 12 * pow(e, 3) * sin(3 * M);
		}

		//Use the semi-major axis, the eccentricity and the calculated true anomaly of the asteroid to find the distance from the sun
		double helioDistance(double a, double e, double v) {
			return a * (1 - pow(e, 2)) / (1 + e * cos(v * PI / 180));
		}

};

class AsteroidNeighborTreeNode {
	public:
		Asteroid asteroid;
		AsteroidNeighborTreeNode* leftChild;
		AsteroidNeighborTreeNode* rightChild;
};

// structure data type called 'AVLT' used to group the asteroid mission rating, name, and left and right children in AVL tree 
struct AVLT {
    double missionRating = -1;
    string name = "";
    struct AVLT* left;
    struct AVLT* right;
};

// declare AVLT structure 
AVLT * tree;

// AVL tree class with methods that enable asteroids to be sorted by their mission rating 
class AVLTree {
public:
    // Constructor
    AVLTree() {
        tree = NULL;
    }

    // Getters
    int getHeight(AVLT * );
    int getBalanceFactor(AVLT * );
	
    // Balance the tree after a node(asteroid) insertion 
    AVLT * balanceTree(AVLT * );
 
    // Insert Node
    void insertNode(AVLT*, Asteroid *);
    
    // Find In-Order Successor
    AVLT * inOrderSuccessor(AVLT * );
    
    // Tree Traversal - visits nodes with asteroid mission ratings from smallest to largest
    void inOrderTraversal(AVLT *, vector<pair<double, string>>& values);

    // Tree Rotations
    AVLT * rightRightRotation(AVLT * );
    AVLT * leftLeftRotation(AVLT * );
    AVLT * leftRightRotation(AVLT * );
    AVLT * rightLeftRotation(AVLT * );
};

AsteroidNeighborTreeNode* addAsteroidNeighborTreeNode(AsteroidNeighborTreeNode* root, Asteroid asteroid) {
	if (root == nullptr) {
		AsteroidNeighborTreeNode* newNode = new AsteroidNeighborTreeNode();
		newNode->asteroid = asteroid;
		return newNode;
	}
	else {
		if (asteroid.x <= root->asteroid.x) {
			root->leftChild = addAsteroidNeighborTreeNode(root->leftChild, asteroid);
		}
		else {
			root->rightChild = addAsteroidNeighborTreeNode(root->rightChild, asteroid);
		}
		return root;
	}

}

int countNearNeighbors(AsteroidNeighborTreeNode* node, Asteroid asteroid) {
	if (node == nullptr)
		return 0;

	//Assume asteroid is too far away for eithr of its children nodes to be worth exploring
	int leftChildCount = 0;
	int rightChildCount = 0;

	//If the asteroid is within MAX_NEIGHBOR_DISTANCE to the left then it's left childe might also be a nieghbor. Recurse!
	if (node->asteroid.x >= asteroid.x - MAX_NEIGHBOR_DISTANCE) {
		leftChildCount = countNearNeighbors(node->leftChild, asteroid);
	}
	
	//If the asteroid is within MAX_NEIGHBOR_DISTANCE to the right then it's right child might also be a neighbor. Recurse!
	if (node->asteroid.x <= asteroid.x + MAX_NEIGHBOR_DISTANCE) {
		rightChildCount = countNearNeighbors(node->rightChild, asteroid);
	}

	if (node->asteroid.id == asteroid.id) {
		return leftChildCount + rightChildCount;
	}
	else {
		if (pow(node->asteroid.x - asteroid.x, 2) + pow(node->asteroid.y - asteroid.y, 2) + pow(node->asteroid.z - asteroid.z, 2) <= MAX_NEIGHBOR_DISTANCE_SQUARED) {
			return leftChildCount + rightChildCount + 1;
		}
		else {
			return leftChildCount + rightChildCount;
		}
	}
}

void rateAsteroids(vector<Asteroid>& asteroids, int missionType) {

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	uniform_real_distribution<double> distribution(0.0, 1.0);

	if (missionType == RANDOM_MISSION) {
		for (int i = 0; i < asteroids.size(); i++) {
			asteroids[i].missionRating = distribution(generator);
		}
		return;
	}
	//Calculate how many other asteroids are near neighbors, allowing a single mission to take pictures of multiple objects
	if (missionType == SCIENCE_MISSION) {

		system_clock::time_point startTime = system_clock::now();
		cout << " Calculating near neighbors for asteroids..." << endl;
		
		//Comparing every asteroid to every asteroid to find neighbors is N^2 run time and not realistic for a million asteroids
		//Sort asteroids into a tree first so we can more efficiently look at only asteroids that are within a reasonable X distance of ourselvevs
		AsteroidNeighborTreeNode* neighborTree = nullptr;
		for (Asteroid asteroid : asteroids) {
			neighborTree = addAsteroidNeighborTreeNode(neighborTree, asteroid);
		}
		//Now find neighbor count
		cout << "Done creating near neighbor tree" << endl;
		cout << "Use tree to calculate near neighbors" << endl;
		for (int i = 0; i < asteroids.size(); i++) {
			asteroids[i].missionRating = (double)countNearNeighbors(neighborTree, asteroids[i]);
			if ((i+1) % 10000 == 0) {
				cout << "Near neighbors calculated for " << i + 1 << " asteroids..." << endl;
			}
		}
		cout << "All near neighbor counts calculated" << endl;
		system_clock::time_point endTime = system_clock::now();

		cout << "Calculating near neighbors took: " << duration_cast<milliseconds>(endTime - startTime).count() << " milliseconds" << endl;
		return;
	}
	//Detect large magnitude asteroids (easy to land on) divided by how far away it is from the earths average raidus of 1 au
	if (missionType == MINING_MISSION) {
		for (int i = 0; i < asteroids.size(); i++) {
			asteroids[i].missionRating = asteroids[i].magnitude/(sqrt(pow(asteroids[i].x,2)+pow(asteroids[i].y,2)+pow(asteroids[i].z,2))-1);
		}
		return;
	}
}

void runQuickSort(vector<Asteroid>& asteroids, int l, int r) {
	if (r - l == 1) {
		if (asteroids[l].missionRating > asteroids[r].missionRating) {
			swap(asteroids[l], asteroids[r]);
		}
		return;
	}
	int p = (r + l) / 2;
	int begin = l;
	int end = r;

	swap(asteroids[p], asteroids[r]);
	p = r;

	while (l <= r) {
		while ((asteroids[l].missionRating < asteroids[p].missionRating)) { l++; }
		while ((r >= l) && (asteroids[r].missionRating >= asteroids[p].missionRating)) { r--; }
		if (r > l) {
			swap(asteroids[l], asteroids[r]);
		}
	}
	swap(asteroids[p], asteroids[l]);

	//	p = l;
	if ((l - begin) > 1) { runQuickSort(asteroids, begin, l - 1); }
	if ((end - l) > 1) { runQuickSort(asteroids, l + 1, end); }
	return;
}

void runTreeAlgorithm(int missionType, int resultCount){
	cout << "Tree not implemented yet" << endl;
}

void runHeapAlgorithm(vector<Asteroid>& asteroids, int resultCount){
	
	Asteroid* asteroidHeap = new Asteroid[resultCount];
	int heapSize = 0;
	
	for(Asteroid& asteroid: asteroids){
		//Sort into heap		
		//If the heap is not yet full, insert and bubble up
		if(heapSize < resultCount){
			asteroidHeap[heapSize]=asteroid;

			int currentHeapPoint = heapSize;
			while(true){
				//We have reached the top of the heap, time to exit the loop
				if(currentHeapPoint == 0){
					heapSize++;
					break;
				}
				int heapComparisonPoint = currentHeapPoint % 2 == 1 ? (currentHeapPoint - 1)/2 : (currentHeapPoint - 2) /2;
				if(asteroidHeap[heapComparisonPoint].missionRating > asteroid.missionRating){
					Asteroid temp = asteroidHeap[heapComparisonPoint];
					asteroidHeap[heapComparisonPoint] = asteroid;
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
			if(asteroid.missionRating > asteroidHeap[0].missionRating){
				//Swap the value into the "min" slot and bubble it down as needed
				asteroidHeap[0]=asteroid;
				int heapComparisonPoint = 0;
				while(heapComparisonPoint * 2 + 2 <= heapSize){
					Asteroid leftChild = asteroidHeap[heapComparisonPoint * 2 + 1];
					if (heapComparisonPoint * 2 + 2 == heapSize) {
						if (asteroid.missionRating < leftChild.missionRating)
							break;
						else {
							asteroidHeap[heapComparisonPoint] = leftChild;
							asteroidHeap[heapComparisonPoint * 2 + 1] = asteroid;
							heapComparisonPoint = heapComparisonPoint * 2 + 1;
						}
					}
					else {
						Asteroid rightChild = asteroidHeap[heapComparisonPoint * 2 + 2];
						if (asteroid.missionRating < leftChild.missionRating && asteroid.missionRating < rightChild.missionRating) {
							break;
						}
						else if (leftChild.missionRating < rightChild.missionRating) {
							asteroidHeap[heapComparisonPoint] = leftChild;
							asteroidHeap[heapComparisonPoint * 2 + 1] = asteroid;
							heapComparisonPoint = heapComparisonPoint * 2 + 1;
						}
						else {
							asteroidHeap[heapComparisonPoint] = rightChild;
							asteroidHeap[heapComparisonPoint * 2 + 2] = asteroid;
							heapComparisonPoint = heapComparisonPoint * 2 + 2;
						}
					}
				}
			}
		}
		
		
	}
	
	//Now extract and reverse the min heap to get the max N results
	Asteroid* sortedResults = new Asteroid[resultCount];
	int resultsSortedSoFar = 0;
	
	
	while(heapSize > 0){
				//Get minimum match from front of heap, swap in last value, resort heap
				sortedResults[resultCount - 1 - resultsSortedSoFar] = asteroidHeap[0];
				resultsSortedSoFar++;
				
				Asteroid newRoot = asteroidHeap[heapSize-1];
				asteroidHeap[0]=newRoot;
				heapSize--;
				int heapComparisonPoint = 0;
				while(heapComparisonPoint * 2 + 2 <= heapSize){
					Asteroid leftChild = asteroidHeap[heapComparisonPoint * 2 + 1];
					Asteroid rightChild = asteroidHeap[heapComparisonPoint * 2 + 2];
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
	cout << "Asteroid Official Name : Calculated Mission Value" << endl;
	for(int i = 0; i < resultCount; i++){
		cout << sortedResults[i].name << " : " << sortedResults[i].missionRating << endl;
	}
	

	
}

int main(){
	cout << "Loading asteroid data..." << endl;
	
	//Asteroid import code goes here
	vector<Asteroid> asteroids;

	ifstream asteroidReader("asteroids_all.csv");
	if (!asteroidReader.is_open()) {
		cout << "Error opening data file" << endl;
		return 0;
	}

	string line;
	getline(asteroidReader, line); //Read and discard the header line
	int count = 0;
	cout << "Loading asteroid data and computing orbits..." << endl;
	while (getline(asteroidReader, line))
	{
		string fields[32];
		//First field is the only quoted, non-numeric field so treat it differently
		int nextPos = line.find("\",");
		fields[0] = line.substr(1, nextPos -1);


		//Everything else is a number, so we can predictably loop through them
		int prevPos = nextPos + 2;
		int fieldCount = 1;
		while (true) {
			nextPos = line.find(',', prevPos);
			if (nextPos > 0) {
				fields[fieldCount] = line.substr(prevPos, nextPos - prevPos);
			}
			else {
				fields[fieldCount] = line.substr(prevPos);
				break;
			}
			
			fieldCount++;
			prevPos = nextPos + 1;
		}

		//Now create an asteroid from this
		Asteroid asteroid;
		asteroid.id = count + 2; //Track line number from original file
		asteroid.name = fields[0];
		asteroid.a = stod(fields[1]);
		asteroid.e = stod(fields[2]);
		asteroid.i = stod(fields[3]);
		asteroid.om = stod(fields[4]);
		asteroid.w = stod(fields[5]);
		asteroid.ma = stod(fields[31]);
		asteroid.epoch_mjd = stoi(fields[25]);
		asteroid.magnitude = fields[14] != "" ? stod(fields[14]) : 0.0;

		asteroid.calculatePosition();

		asteroids.push_back(asteroid);

		count++;
		if (count % 10000 == 0) {
			cout << count << " asteroids loaded..." << endl;
		}
		if (count == 100000) {
			break;
		}
	}
	asteroidReader.close();

	cout << "Asteroids loaded and orbital locations computed!" << endl;

	/*for(int i = 0; i < 500000; i++)
	{
		Asteroid sampleAsteroid;
		sampleAsteroid.id = i;
		asteroids.push_back(sampleAsteroid);
	}*/
	
	
	cout << "Choose your mission criteria" << endl;
	cout << "\t1) Science Mission: Asteroids with near neighbors" << endl;
	cout << "\t2) Mining Mission: Large, metal rich easy to land on asteroids" << endl;
	cout << "\t3) Random Mission: Assign asteroids a random mission value" << endl;
	
	int missionTypeChoice;
	cin >> missionTypeChoice;
	rateAsteroids(asteroids, missionTypeChoice);


	cout << "How many results do you want? [1 to 600,000]" << endl;
	
	int resultCount;
	cin >> resultCount;
	
	cout << "Choose algorithm for sorting results" << endl;
	cout << "\t1)Quicksort" << endl;
	cout << "\t2)AVL Tree" << endl;
	cout << "\t3)Heap" << endl;
	
	int algorithmTypeChoice;
	cin >> algorithmTypeChoice;
	
	system_clock::time_point startTime = system_clock::now();
	
	if (algorithmTypeChoice == QUICKSORT_ALGORITHM) {
		runQuickSort(asteroids, 0, asteroids.size() - 1);
		cout << "Sorted results" << endl;
		int lastAsteroid = asteroids.size() - 1;
		for (int i = 0; i < resultCount; i++) {
			cout  << asteroids[lastAsteroid - i].name << " : " << asteroids[lastAsteroid - i].missionRating << endl;
		}

	}
	else if(algorithmTypeChoice == AVL_TREE_ALGORITHM){
		runTreeAlgorithm(missionTypeChoice, resultCount);
	}
	else if(algorithmTypeChoice == HEAP_ALGORITHM){
		runHeapAlgorithm(asteroids, resultCount);
	}
	
	system_clock::time_point endTime = system_clock::now();
	
	cout << "Algorithm took: " << duration_cast<milliseconds>(endTime - startTime).count() << " milliseconds" << endl;
	return 0;
	
}
