#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <numeric>
#include <math.h>
#include <algorithm>

using namespace std;

//readFile is responsible for obtaining the starting index the input file which is stored in startingIndex
double readFile(string filename, vector<double>& inputData)
{
	ifstream f;
	string line;
	double valid;
	double startingIndex = 0;

	f.open(filename);
	if (!f)  // If file not found, exit for now but should prompt the user to re-type the filename, prolly another function 
	{
		cerr << "File not found. Exiting program."; // @JG @ Kobe Maybe add something that will loop back to open file
		exit(1);
	}

	// Read first line
	getline(f, line);
	stringstream ss;

	// assign the first line to ss
	ss << line;

	// Get first entity (separated by spaces of first line) and assign it to startIndex
	ss >> startingIndex; // If invalid ung first value, gawa nalang tayo ng 'Invalid file content'
	if ((ss >> valid)) // If there is still another element ie input 1 input 2 (e.g. -3 5)
	{
		inputData.push_back(valid);
	}
	else // If none, then the first value is already the first datapoint and startingIndex = 0
	{
		inputData.push_back(startingIndex);
		startingIndex = 0;
	}

	while (getline(f, line)) // Read succeeding lines until EOF is reached
	{
		stringstream ss;
		ss << line;
		if (ss >> valid)
		{
			inputData.push_back(valid);
		}
	}

	f.close();

	return startingIndex;
}

bool contains(vector<double> vec, double &elem)
{
    bool result = false;
    if(find(vec.begin(), vec.end(), elem) != vec.end() )
    {
      result = true;
    }
    return result;
}

//crossCorrelation responsible for creating the array of possible crosscorrelation values given two signals
vector <double> crossCorrelation(vector<double> x, vector<double> y, int startX, int startY, int shift, double z, int &outputStartIndex)
{
  vector <double> Z;
  int endX = startX + x.size() - 1;
  int endY = startY + y.size() - 1;
  
  int ccStartIndex = startX - endY;
  int ccEndIndex = endX - startY;

  int duration = ccEndIndex - ccStartIndex + 1;

  for(shift = ccStartIndex; shift <= ccEndIndex; shift++) 
  {
    z = 0;
    for (int i = 0; i < x.size(); i++)
    {
      if(startX < startY || startX == startY)
      {
        if((contains(y, y[i - (startY - startX) - shift]) == 1))
          z += x[i] * y[i - (startY - startX) - shift];
        else
          continue; 
      }
      else if(startX > startY)
      {
        if((contains(y, y[i + (startX - startY) - shift]) == 1))
          z += x[i] * y[i + (startX - startY) - shift];
        else
          continue;
      }
    }
    Z.push_back(z);
  }

  outputStartIndex = ccStartIndex;

  return Z;
}

// showData prints out all the elements in a vector
// Got this from https://stackoverflow.com/questions/10750057/how-do-i-print-out-the-contents-of-a-vector
// range-based for loop

void showData(vector<double> data)
{
	for (double dataPoint : data)
	{
		cout << dataPoint << endl;
	} 
}

// average is responsible for taking the averages of the elements in a given vector. This is used to refine the input signal vector
// Got the average accumulate here: https://stackoverflow.com/questions/28574346/find-average-of-input-to-vector-c

double average(vector<double> data)
{
	double ave;

	ave = accumulate(data.begin(), data.end(), 0.0)/data.size();
	cout << "The average is " << ave << endl;

	return ave;
}

// subtractAverages subtracts the average from the ith vector element.
vector<double> subtractAverages(vector<double> inputData)
{
	vector<double> newData;
	double ave;
	ave = average(inputData);

	for (int i = 0; i < inputData.size(); i++)
	{
		newData.push_back(inputData[i] - ave);
	}

	return newData;
}

// autocorrelation takes the crosscorrelation of the input signal vector to itself while letting l = 0
double autocorrelation(vector<double> inputData)
{
	double sum = 0;

	for (int i = 0; i < inputData.size(); i++)
	{
		sum += inputData[i] * inputData[i];
	}

	return sum;
}

// normalizationCoefficients returns the normalization coefficient needed for normalizing the resulting crosscorrelation; reducing them into the values that range from -1 to 1
double normalizationCoefficients(vector<double> x, vector<double> y)
{
	return sqrt(autocorrelation(x)*autocorrelation(y));
}

// normalizedCrosscorrelation outputs the resulting vector containing all the possible normalized crosscorrelation values for a given pair of signal vectors
vector <double> normalizedCrosscorrelation(double normCoef, vector <double> crosscorreVec)
{
  vector <double> normVec;

  for (double dataPoint : crosscorreVec)
	{
    normVec.push_back(dataPoint / normCoef);
	} 

  return normVec;
}

// outputFile writes the values from the normalizedCrosscorrelation vector to a given filename in .txt format
void outputFile(string filename, int startIndex, vector<double> outputData)
{
  ofstream f;
  cout << "Program is saving your output to: " << filename << endl;
  filename = filename + ".txt";

  f.open(filename);

  if (startIndex != 0) // To copy input format of index data_1
  {
    f << startIndex << ' '; 
  }

  for (double output : outputData)
  {
    f << output << endl;
  }

  f.close();

}


// main function
int main()
{
  vector <double> xRaw;
  vector <double> yRaw;

  string xFilename;
	cout << "Please input filename for x: " << endl;
	cin >> xFilename;

	string yFilename;
	cout << "Please input filename for y: " << endl;
	cin >> yFilename;

  string outputFilename;
  cout << "Please input your desired output filename for the normalized cross-correlated signals (No need to include .txt): " << endl;
  cin >> outputFilename;

  int startX = readFile(xFilename, xRaw);
  int startY = readFile(yFilename, yRaw);

  vector <double> xSubtractedAverages = subtractAverages(xRaw); 
  vector <double> ySubtractedAverages = subtractAverages(yRaw);

  int shift;
  int outputStartIndex;
  double z;

  double autocorrelationX = autocorrelation(xSubtractedAverages);
  double autocorrelationY = autocorrelation(ySubtractedAverages);
  double normCoef = normalizationCoefficients(xSubtractedAverages, ySubtractedAverages);

  cout << "X Autocorrelation: " << autocorrelationX << endl;
  cout << "Y Autocorrelation: " << autocorrelationY << endl;
  cout << "Normalization Coefficient: " << normCoef << endl;

  vector <double> result = crossCorrelation(xSubtractedAverages, ySubtractedAverages, startX, startY, shift, z, outputStartIndex);

  vector <double> normVec = normalizedCrosscorrelation(normCoef, result);

  outputFile(outputFilename, outputStartIndex, normVec);
}