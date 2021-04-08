// This project was made by Kent Marc Kobe Bismark
// Gerson Cruz and Juan Glicerio Manlapaz
// in accordance with the requirements of ENGG 151.01

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;

// Class for the LTI system
class LTI
{
public:
  // Function for showing data in a vector
  void showData(vector<double> data)
  {
    for (int i = 0; i < data.size(); i++)
    {
      cout << data[i] << endl;
    }
  }

  // Function for adding zeroes to signal x where n < 0
  vector<double> addZeroesX(vector<double> LTICoeff)
  {
    vector<double> zeroes;
    double zero = 0.0;
    for (int i = 0; i < LTICoeff.size() - 1; i++)
    {
      zeroes.push_back(zero);
    }
    return zeroes;
  }

  // Function for adding zeroes to signal y where n < 0
  vector<double> addZeroesY(vector<double> LTICoeff)
  {
    vector<double> zeroes;
    double zero = 0.0;
    for (int i = 0; i < LTICoeff.size(); i++)
    {
      zeroes.push_back(zero);
    }
    return zeroes;
  }

  // Function for get the sum-product of two arrays; just like
  // the SUMPRODUCT() function in Excel
  double sumProduct(vector<double> array1, vector<double> array2)
  {
    double sumProductResult = 0;
    int n = 0;
    string errorMessage = "Error: Array lengths not equal";

    if (array1.size() == array2.size())
    {
      n = array1.size();
      for (int i = 0; i < n; i++)
      {
        sumProductResult += array1[i] * array2[i];
      }
      return sumProductResult;
    }
    return 0;
  }

  // Function to store the elements of a signal vector;
  // needed for the computation of the LTI output.
  vector<double> xTempVector(vector<double> signal,
                             vector<double> xNegatives)
  {
    vector<double> xTemp;

    for (int i = 0; i < xNegatives.size(); i++)
    {
      xTemp.push_back(xNegatives[i]);
    }
    xTemp.push_back(signal[0]);

    return xTemp;
  }

  // Function to store the elements of a LTI output vector;
  // needed for the computation of the LTI output.
  vector<double> yTempVector(vector<double> yNegatives)
  {
    vector<double> yTemp;

    for (int i = 0; i < yNegatives.size(); i++)
    {
      yTemp.push_back(yNegatives[i]);
    }

    return yTemp;
  }

  // Function to store LTI outputs in a vector.
  vector<double> LTIOutput(vector<double> NRC,
                           vector<double> RC,
                           vector<double> &signal)
  {
    vector<double> output;
    vector<double> xNegatives = addZeroesX(NRC);
    vector<double> yNegatives = addZeroesY(RC);
    vector<double> xTemp = xTempVector(signal, xNegatives);
    vector<double> yTemp = yTempVector(yNegatives);
    double outputValue;

    // Reverse NRC and RC arrays before LTI output computation
    reverse(RC.begin(), RC.end());
    reverse(NRC.begin(), NRC.end());

    for (int i = 1; i <= signal.size(); i++)
    {
      // outputvalue is computed, stored in vector variable output
      outputValue = -sumProduct(RC, yTemp) + sumProduct(NRC, xTemp);
      output.push_back(outputValue);

      // for loop for changing elements in vector xTemp.
      for (int j = 0; j < xTemp.size(); j++)
      {
        // All values in xTemp vector shifted to the left
        xTemp.at(j) = xTemp[j + 1];

        // Once j is now the last possible index for xTemp,
        // the element signal[i] is appended to xTemp
        if (j == xTemp.size() - 1)
        {
          xTemp.at(j) = signal[i];
        }
      }

      for (int j = 0; j < yTemp.size(); j++)
      {
        // All values in the yTemp vector are shifted to the left
        yTemp.at(j) = yTemp[j + 1];

        // Once j is now the last possible index for yTemp,
        // the element output[i-1] is appended to yTemp
        if (j == yTemp.size() - 1)
        {
          yTemp.at(j) = output[i - 1];
        }
      }
    }
    return output;
  }
};

// Function to return a boolean whether the input is a valid input
// without extra characters
bool isDouble(string input, double &number)
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

// readFile is responsible for reading the LTI system file
double readFile(string filename, vector<double> &inputData,
                bool &error, vector<double> &a_n,
                vector<double> &b_n, int &non_recur, int &recur)
{
  ifstream f;
  string line;
  string extraTest;
  double valid;
  double startingIndex = 0;
  error = false;

  int nrc = 0;
  int rc = 0;

  f.open(filename);
  // If file not found, prompt the user to re-type filename
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
  // Also checks if first entity is invalid
  ss >> extraTest;
  if (isDouble(extraTest, valid))
  {
    startingIndex = valid;
  }

  else // Prompts the user to check the file and re-input.
  {
    cout << "\nInvalid input. Please check your input file.\n";
    cout << "Here is the invalid line: " << line << endl;
    cout << "Check your file and input the filename again.\n";
    error = true;
    return 0;
  }

  // If there is still another element i.e. input 1 input 2
  // (e.g. -3 5), value is assigned to nrc
  if ((ss >> valid))
  {
    cout << valid << endl;
    nrc = valid;
  }
  // If none, then the first value is already the nrc
  // and startingIndex = 0
  else
  {
    nrc = startingIndex;
    startingIndex = 0;
  }

  // Get second line
  getline(f, line);
  stringstream s;
  s << line;

  // If valid, then the data in the second line is rc
  // the number of recursive coefficients
  if ((s >> valid))
  {
    rc = valid;
  }

  // For placing the non recursive coefficients into a vector
  for (int i = 0; i < nrc; i++)
  {
    getline(f, line);
    stringstream ss;
    ss << line;
    if (ss >> valid)
    {
      b_n.push_back(valid);
    }
  }

  // For placing recursive coefficients into a vector
  for (int i = 1; i <= rc; i++)
  {
    getline(f, line);
    stringstream ss;
    ss << line;
    if (ss >> valid)
    {
      a_n.push_back(valid);
    }
  }

  f.close();

  non_recur = nrc;
  recur = rc;

  return startingIndex;
}

// Function for reading a signal file
// separate for the LTI system file
double readSignalFile(string filename, vector<double> &signalData,
                      bool &error)
{
  ifstream f;
  string line;
  string extraTest;
  double valid;
  double startingIndex = 0;
  error = false;

  f.open(filename);
  // If file not found, prompt the user to re-type filename
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
  // and assign it to startIndex.
  // Also checks if first entity is invalid
  ss >> extraTest;
  if (isDouble(extraTest, valid))
  {
    startingIndex = valid;
  }

  else // Prompts the user to check the file and re-input.
  {
    cout << "\nInvalid input. Please check your input file.\n";
    cout << "Here is the invalid line: " << line << endl;
    cout << "Check your file and input the filename again.\n";
    error = true;
    return 0;
  }

  // If there is still another element
  // i.e. input 1 input 2 (e.g. -3 5)
  if ((ss >> valid))
  {
    signalData.push_back(valid);
  }
  // If none, then the first value is already
  // the first datapoint and startingIndex = 0
  else
  {
    signalData.push_back(startingIndex);
    startingIndex = 0;
  }

  // Read succeeding lines until EOF is reached
  while (getline(f, line))
  {
    stringstream ss;
    ss << line;
    if (ss >> valid)
    {
      signalData.push_back(valid);
    }
  }

  f.close();

  return startingIndex;
}

// -------------- MAIN FUNCTION -----------

int main()
{
  ofstream f;
  LTI ltisystem;
  string LTIFilename, SignalFilename, logFilename, opt_form, sigVal;
  bool error = false, exit = false, logfile = false, ltiread = false;
  char opt = ' ';
  vector<double> LTI, signal, a_n, b_n, LTIOutput, xtemp, ytemp;
  double signalValue;
  // rc: recursive, nrc: non recursive
  int rc = 0, nrc = 0, xcnt = 0, ycnt = 0;

  do
  {
    switch (opt)
    {
    case '1': // Load file containing specs for LTI system
    {
      // If 1 is pressed again, the LTI file will be overwritten.

      f.close();
      do
      {
        cout << "Please input LTI Filename:\t";
        // LTI system clear
        a_n.clear();
        b_n.clear();
        LTI.clear();
        signal.clear();
        xtemp.clear();
        ytemp.clear();
        LTIOutput.clear();

        ycnt = 0;
        xcnt = 0;
        cin >> LTIFilename;
        double startindexLTI = readFile(LTIFilename, LTI, error,
                                        a_n, b_n, nrc, rc);
      } while (error == true);

      cout << "LTI System Loaded from " << LTIFilename;

      if (logfile == true)
      {
        cout << "\nCurrent logfile has now been closed. ";
        cout << "Please enter a new log file if desired.\n";
        logfile = false;
      }

      ltiread = true;
      opt = ' ';
      break;
    }
    case '2': // readFile Log
    {
      do
      {
        cout << "Please input Log Filename:\t";
        cin >> logFilename;
      } while (error == true);
      f.open(logFilename);

      logfile = true;
      cout << "Log File Created: " << logFilename;

      opt = ' ';
      break;
    }
    case '3': // Display NRC, RC, inputs, and outputs
    {
      // If LTI file has not yet been read, cannot enter display
      if (!ltiread)
      {
        cout << "No LTI file read yet. Cannot display system ";
        cout << "info. Please enter an LTI file first. \n";
        opt = ' ';
        break;
      }

      cout << "\nLTI System Information:";
      cout << "\nNumber of Recursive Coef: " << rc << endl;

      // for displaying recursive coefs
      for (double n : a_n)
      {
        cout << n << "\t";
      }
      cout << "\nNumber of Non-recursive coef: " << nrc << endl;

      // for displaying non-recursive coefs
      for (double n : b_n)
      {
        cout << n << "\t";
      }

      // If empty, outputs and inputs are considered 0, 0.
      // Else, inputs and outputs depend on nrc and rc and
      // previous calculations found in xtemp and ytemp
      if (xtemp.size() != 0)
      {
        for (int n = 0; n < nrc; n++)
        {
          cout << "\nx(" << 0 - ((nrc - 1) - n) << ") = "
               << xtemp[xtemp.size() - (nrc - n)];
        }
        for (int n = 0; n < rc; n++)
        {
          cout << "\ny(" << 0 - ((rc - 1) - n) << ") = "
               << ytemp[ytemp.size() - (rc - n)];
        }
      }
      else
      {
        for (int n = 0; n < nrc; n++)
        {
          cout << "\nx(" << 0 - ((nrc - 1) - n) << ") = " << 0;
        }
        for (int n = 0; n < rc; n++)
        {
          cout << "\ny(" << 0 - ((rc - 1) - n) << ") = " << 0;
        }
        cout << "\nNo Detected Input";
      }

      opt = ' ';
      break;
    }
    case '4': // Interactive Mode
    {
      // If LTI file is not yet read, cannot enter interactive
      if (!ltiread)
      {
        cout << "No LTI file read yet. Cannot enter interactive ";
        cout << "mode. Please enter an LTI file first. \n";
        opt = ' ';
        break;
      }

      cout << "Entered Interactive Mode.";
      cout << " Press any letter to exit." << endl;
      bool truth = true;

      while (truth)
      {
        cout << "Input Value: ";
        cin >> sigVal;
        // If input is not double, enter any character, exits
        if (cin.fail() || !isDouble(sigVal, signalValue))
        {
          cout << "\nInvalid input. End of Interactive mode.";
          cin.clear();
          cin.ignore(10000, '\n');
          opt = ' ';
          break;
        }
        // xtemp and ytemp are running vectors
        // containing all previous inputs and outputs
        // xcnt and ycnt take note of where in the vector
        // the current input/output is
        signal.push_back(signalValue);
        xtemp.push_back(signalValue);
        xcnt++;

        LTIOutput = ltisystem.LTIOutput(b_n, a_n, xtemp);
        f << LTIOutput[ycnt] << endl;
        ytemp.push_back(LTIOutput[ycnt]);
        cout << "Output from LTI: " << ytemp[ycnt] << endl;
        ycnt++;
      }
      break;
    }
    case '5': // Specify a signal file from which the next inputs
    {         // to the system will be obtained

      signal.clear();
      // If LTI file has not yet been read, cannot solve any system
      if (!ltiread)
      {
        cout << "No LTI file read yet. Cannot load and solve ";
        cout << "signal file. Please enter an LTI file first. \n";
        opt = ' ';
        break;
      }

      do
      {
        cout << "Please input Signal Filename:\t";
        cin >> SignalFilename;
        double startindexSignal = readSignalFile(SignalFilename,
                                                 signal, error);
      } while (error == true);

      for (double n : signal)
      {
        xtemp.push_back(n);
        xcnt++;
      }
      cout << "Signal Input for LTI System: \n";
      ltisystem.showData(signal);

      // Calculate LTI Outputs for SignalFile
      LTIOutput = ltisystem.LTIOutput(b_n, a_n, signal);
      cout << "LTI Output: \n";
      ltisystem.showData(LTIOutput);
      for (double n : LTIOutput)
      {
        f << n << endl;
        ytemp.push_back(n);
        ycnt++;
      }

      opt = ' ';
      break;
    }
    case '6': // Clear Inputs and Outputs and Reset to 0
    {
      // Empties the vectors, the previous inputs and outputs
      // When displayed in '3' will be 0.

      f.close();
      signal.clear();
      xtemp.clear();
      ytemp.clear();
      LTIOutput.clear();

      ycnt = 0;
      xcnt = 0;

      cout << "Inputs and Outputs Cleared and Reset to 0.\n";

      if (logfile)
      {
        cout << "Current logfile has now been closed. ";
        cout << "Please enter a new log file if desired.\n";
        logfile = false;
      }

      opt = ' ';
      break;
    }
    case '7': // Exit Application
    {
      f.close();
      exit = true;

      if (exit)
      {
        cout << "Goodbye! Thank you for your patronage.\n";
        return 0;
      }
      break;
    }
    default: // Default View for LTI Menu Options
      do
      {
        cout << "\n\n    Options: \n";
        cout << "[1] Load File Containing LTI System Specs\n";
        cout << "[2] Specify a Log File\n";
        cout << "[3] Display LTI System Info\n";
        cout << "[4] Interactive Mode\n";
        cout << "[5] Load Signal File and Solve System\n";
        cout << "[6] Clear and Reset Inputs and Outputs\n";
        cout << "[7] Exit Application\n";
        cout << "Input Option:\t";
        cin >> opt_form;
        cin.clear();
        opt = opt_form[0];
        if (opt != '1' && opt != '2' && opt != '3' && opt != '4' 
            && opt != '5' && opt != '6' && opt != '7' 
            || opt_form.length() > 1)
        {
          cout << "\nInvalid Input. Please try again!";
        }
      } while (opt != '1' && opt != '2' && opt != '3' && opt != '4' 
               && opt != '5' && opt != '6' && opt != '7' 
               || opt_form.length() > 1);
      break;
    }
  } while (!exit);
}