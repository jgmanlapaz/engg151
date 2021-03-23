// This Project on the Discrete Fourier Transform
// was made by Kent March Kobe Bismark, Gerson Cruz
// and your favorite student, Juan Glicerio Manlapaz :) 

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;

// This function returns a boolean whether the input 
// is a valid input without extra characters
bool isDouble(string input, double& number)
{
  stringstream ss(input);

  ss >> number;

  // If the stringstream is indeed a double with no extra characters
  // and is at the end of the line
  if (!ss.fail() && ss.eof())
  {
    return true;
  }
  else
  {
    return false;
  }
}

// readFile is responsible for obtaining the starting 
// index of the input file which is stored in startingIndex
// readFile also creates a vector containing all the inputs
double readFile(string filename, vector<double>& inputData, 
bool& error)
{
  ifstream f;
  string line;
  string extraTest;
  double valid;
  double startingIndex = 0;
  error = false;
  
  f.open(filename);
  // If file not found, prompt the user to 
  // re-type the filename
  if (!f)  
  {
    cerr << "File not found. Please input a valid filename.\n";
    error = true;
    return 0;
  }
  
  // Read first line
  getline(f, line);
  stringstream ss;
  
  // Assign the first line to ss
  ss << line;
  
  // Get first entity (separated by spaces of first line)
  // and assign it to startIndex
  // Code below does the checking too if first entity is invalid
  ss >> extraTest;
  if (isDouble(extraTest, valid))
  {
    startingIndex = valid;
  }

  else // Prompts the user to check the file and re-input.
  { 
    cout << "\nInvalid input. Please check your input file.\n";
    cout << "Here is the invalid line: " << line << endl; 
    // Here is the invalid line; 
    // outputs empty space if file is empty.
    cout << "Check your file and input the filename again.\n";
    error = true;
    return 0;
  }

  // If there is still another element 
  // i.e. input 1 input 2 (e.g. -3 5)
  if ((ss >> valid)) 
  {
    inputData.push_back(valid);
  }
  // If none, then the first value is already 
  // the first datapoint and startingIndex = 0
  else 
  {
    inputData.push_back(startingIndex);
    startingIndex = 0;
  }
  
  // Read succeeding lines until EOF is reached
  while (getline(f, line)) 
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

// Converts analog frequencies to digital frequencies 
// using the sampling frequency
double digitalFrequency(double analogFreq, double samplingFreq)
{
  return 2 * M_PI * (analogFreq/samplingFreq);
}

// Gets the SUMPRODUCT of the cosine/real part
double autocorr_cos(vector<double> signalData, double digitalFreq)
{
  double cosVec = 0;
  for (int i = 0; i < signalData.size(); i++)
  {
    cosVec = cosVec + cos(digitalFreq * i) * cos(digitalFreq * i);
  }

  return cosVec;
}

// Gets the SUMPRODUCT of the sine/imaginary part
double autocorr_sin(vector<double> signalData, double digitalFreq)
{
  double sinVec = 0;
  for (int i = 0; i < signalData.size(); i++)
  {
    sinVec = sinVec + sin(digitalFreq * i) * sin(digitalFreq * i);
  }

  return sinVec;
}

// Gets the SUMPRODUCT of the signal input
double autocorr_signal(vector<double> signalData)
{
  double autocorr_sig = 0;
  for (int i = 0; i < signalData.size(); i++)
  {
    autocorr_sig = autocorr_sig + signalData[i] * signalData[i];
  }

  return autocorr_sig;
}

// Function for the real part
double cosineCoeff(vector<double> signalData, double digitalFreq)
{
  double cosine_autocorr = autocorr_cos(signalData, digitalFreq);
  double input_autocorr = autocorr_signal(signalData);
  double cosCoeff = 0;
  for (int i = 0; i < signalData.size(); i++)
  {
    cosCoeff += signalData[i] * cos(digitalFreq * i);
  }

  // Normalizing the real part
  return cosCoeff / sqrt(cosine_autocorr*input_autocorr);
}

// Function for the imaginary part
double sineCoeff(vector<double> signalData, double digitalFreq)
{
  double sine_autocorr = autocorr_sin(signalData, digitalFreq);
  double input_autocorr = autocorr_signal(signalData);
  double sinCoeff = 0;
  for (int i = 0; i < signalData.size(); i++)
  {
    sinCoeff += signalData[i] * sin(digitalFreq * i);
  }
  
  // Normalizing the imaginary part
  return sinCoeff / sqrt(sine_autocorr * input_autocorr);
}

// This function calculates how much is the size of every step
// based on the steps the user wants from the start to the
// end frequency
double getstepsize(int steps, double startFreq, double endFreq)
{
  double stepsize = 0;

  stepsize = (endFreq - startFreq) / steps;

  return stepsize;
}

// This is the main Discrete Fourier Transform function
// It takes in all the user's inputs, uses them to compute the DFT
// at specific steps and outputs the results to either one file
// or two files depending on the user's preferences
void DFT(char form, double startFreq, double endFreq, 
double samplingFreq, int steps, double stepsize, 
vector<double> signal)
{
  double real, imaginary, magnitude, phase;
  double analogFreq = startFreq;
  string fn1, fn2;
  ofstream f1; // File outputting the results to a file
  ofstream f2;
  
  // This determines whether one output file or two will be made
  if(form == 'R' || form == 'M')
  {
    cout << "\nPlease input your desired output filename for the";
    cout << "\noutput signals: ";
    cin >> fn1;
    f1.open(fn1);
  }
  else if(form == 'B')
  {
    cout << "\nPlease input your desired output filename for the";
    cout << "\nfirst output signal (Real-Imaginary): ";
    cin >> fn1;
    cout << "\nPlease input your desired output filename for the";
    cout << "\nsecond output signal (Magnitude-Phase): ";
    cin >> fn2;
    f1.open(fn1);
    f2.open(fn2);
  }

  // Main computation loop:
  for (int i = 0; i <= steps; i++)
  {
    double digitalFreq = digitalFrequency(analogFreq, samplingFreq);
    real = cosineCoeff(signal, digitalFreq); 
    imaginary = sineCoeff(signal, digitalFreq) * -1; 
    
    phase = atan2(imaginary, real) * (180/M_PI);
    magnitude = sqrt((imaginary * imaginary) + (real * real));

    // switch statement based on chosen option for file output
    switch(form)
    {
      case 'R':
        f1 << analogFreq << "\t " << real << "\t " << imaginary;
        f1 << endl;
        break;
      case 'M':
        f1 << analogFreq << "\t " << magnitude  << "\t " << phase;
        f1 << endl;
        break;
      case 'B':
        f1 << analogFreq << "\t" << real << "\t" << imaginary;
        f1 << endl;
        f2 << analogFreq << "\t" << magnitude  << "\t" << phase;
        f2 << endl;
        break;
      default:
        break;
    }
    analogFreq += stepsize;
  }
   
  if (form == 'R' || form == 'M')
  {
    cout << "\nOutput file " << fn1 << " created!";
  }
  else if (form == 'B')
  {
    cout << "\nOutput files " << fn1 << " and " << fn2;
    cout << " created!";
  }

  f1.close();
  f2.close();
}

int main()
{
  string signalFilename;
  string str_form;
  vector <double> signal;
  char form;
  int step;
  bool error = false;
  double samplingFreq, digitalFreq;
  double startFreq, endFreq;
  double real, imaginary, phase, magnitude;
  double stepsize;
  double double_check = 0;
  
  do
  {
    cout << "Please input filename:\t";
    cin >> signalFilename;
    double signalIndex = readFile(signalFilename, signal, error);
  }
  while (error == true);

  cout << "\nPlease type in only valid inputs. The program will not";
  cout << " begin unless a valid input is inputted and will keep "; 
  cout << "prompting you to re-type a valid input. Enjoy! \n";

  // do-while statements to ensure that valid inputs are typed;
  // cin.clear() and cin.ignore() to ensure that when an invalid
  // input is typed, cin is clear and the loop is repeated.
  do 
  {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cout << "Please input valid Sampling Frequency (in Hz):\t";
    cin >> samplingFreq;
  }
  while (samplingFreq <= 0 || !cin);

  do
  {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cout << "Please input valid Starting Frequency (in Hz):\t";
    cin >> startFreq;
  }
  while (!cin || startFreq < 0);

  do
  {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cout << "Please input valid End Frequency (in Hz):  \t";
    cin >> endFreq;
  } 
  while (!cin || startFreq >= endFreq);
  
  do
  {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cout << "Please input valid Steps from Start to End:\t";
    cin >> step;
  }
  while (step <= 0 || !cin);

  stepsize = getstepsize(step, startFreq, endFreq);
  cout << "Stepsize: \t" << stepsize << endl;

  cin.clear();
  cin.ignore(256,'\n');
  do
  {
    cout << "[R]eal-Imaginary or [M]agnitude-Phase or [B]oth?\n";
    cout << "Which form: R, M, B? \t";
    getline(cin, str_form);
    form = str_form[0];
    form = toupper(form);
  }
  while (form != 'R' && form != 'M' && form != 'B' 
  || str_form.length() > 1);

  DFT(form, startFreq, endFreq, samplingFreq, 
  step, stepsize, signal);
}