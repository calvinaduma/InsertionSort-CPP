/*
 * Name: Calvin Aduma
 * Date Submitted:
 * Lab Section: 002
 * Assignment Name: Finding the CLosest Pair of Points
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

struct point
{
    double x;
    double y;
};

/*Implement the following function
  Reads in a file specified by the parameter
  Format of file: #of points on first line
                  remaining (#of points) lines: x-value and y-value of point
                  one point per line
                  x-value and y-value are double-precision values and
                  bounded by the unit square 0.0 <= x,y < 1.0
  Should use "spatial hashing" where the #of cells scales with the #of points
  and find the distance between the closest pair of points which will be
  returned as a double type value
*/
double closestPair(string filename);
double distanceFunc(double,double,double,double);
int factorial(int);
void insertionSort(double[],int);

int main()
{
    double min;
    string filename;
    cout << "File with list of points within unit square: ";
    cin >> filename;
    min = closestPair(filename);
    cout << setprecision(16);
    cout << "Distance between closest pair of points: " << min << endl;
    return 0;
}

double closestPair(string filename){
    int numOfPoints, b, i=0;
    double p1, p2;
    struct point coordinate;
    string input;
    ifstream afile;
    afile.open(filename);
    afile >> numOfPoints;
    b = ceil(sqrt(numOfPoints)); // sqrt is taken from the number of points to create a perfect bxb grid of points
    double interval = 1.0/b;
    vector<point> plane[b][b];
    //point distanceGrid[b][b];
    
    // read file into plane
    int increase = 0;
    string xInput, yInput;
    while (afile >> coordinate.x >> coordinate.y){
      //  afile >> xInput >> yInput;
        //p1 = stod(xInput);
        //p2 = stod(yInput);
        //coordinate = {p1,p2};
        int xcoordinate = floor(coordinate.x/interval);
        int ycoordinate = floor(coordinate.y/interval);
        plane[xcoordinate][ycoordinate].push_back(coordinate);
        //increase++;
    }
    //int arraySize = (factorial(numOfPoints)) / (factorial(2)*factorial(numOfPoints - 2)); // this gets the total number of distances that can be calculated from n number of points. Combinatorics formula for nCr.
    //double distanceArray[arraySize];
    //double finalDistanceArray[arraySize]; // over-estimates number of indices.
    int finalIterator=0;
    double minDistance = 2.0;
    // get distances in grid
    for (int row=0; row<b-1; row++){
        for (int column =0; column<b-1; column++){
            int size = plane[row][column].size();  // size of vector = # of points in grid
            vector<point> tempPointGrid;  // array to hold each individual points in vector
            int i;
            double curDistance;
            tempPointGrid.clear();
            i=0;
            // gets each point stored in each vector of the array
            for (auto pointV:plane[row][column]){
                tempPointGrid.push_back(pointV);
                i++;
            }
            for (auto pointV:plane[row][column+1]){
                tempPointGrid.push_back(pointV);
                i++;
            }
            for (auto pointV:plane[row+1][column]){
                tempPointGrid.push_back(pointV);
                i++;
            }
            for (auto pointV:plane[row+1][column+1]){
                tempPointGrid.push_back(pointV);
                i++;
            }
            // gets the distance betweeen each pairs of points in the same grid
            //int iteratorSize = factorial(size); // number of distances per number of points in grid
            //int iterator = 0;
            for (int iteratorOne=0; iteratorOne<i - 1; iteratorOne++){ // minus 1 from overall size so that it will not exceed bounds
                for (int iteratorTwo=iteratorOne + 1; iteratorTwo<i; iteratorTwo++){
                   curDistance = distanceFunc(tempPointGrid[iteratorOne].x,tempPointGrid[iteratorTwo].x,tempPointGrid[iteratorOne].y,tempPointGrid[iteratorTwo].y);
                   if (curDistance < minDistance)
                   {
                       minDistance = curDistance;
                   }
                   //iterator++;
                }
            }
            // sorts distance from greatest to least.
            //insertionSort(distanceArray,iteratorSize);
            
            // puts smalles distance of grid into finalDistanceArray
            //finalDistanceArray[finalIterator] = distanceArray[0]; // smallest distance per grid // c
            //finalIterator++;
        }
    }
    // unnecessary effort to divide the plane into 8 groups and get smallest when a single array can store all distances and record smallest
    //int finalDistanceArraySize = sizeof(finalDistanceArray)/sizeof(finalDistanceArray[0]);
    //insertionSort(finalDistanceArray,finalDistanceArraySize);
    afile.close();
    //return finalDistanceArray[finalDistanceArraySize-1];
    return minDistance;
}

double distanceFunc(double x1, double x2, double y1, double y2){
    double d;
    double data = (pow((x1-x2),2) + pow((y1-y2),2));
    d = sqrt(data);
    return d;
}

int factorial(int n){
    if(n == 0) { return 1; }
    else { return n*factorial(n-1); }
}

void insertionSort(double arr[], int arraySize){
    for (int i=1; i<arraySize; i++){
        double key = arr[i];
        int j = i-1;
        while(key<arr[j] && j>=0){
            arr[j+1] = arr[j];
            j--;
        }
        arr[j+1] = key;
    }
}