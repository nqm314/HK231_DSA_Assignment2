#include "main.h"
#include "restaurant.cpp"

void simulate(string filename) {
  ifstream ss(filename);
  string str, name;
  int maxsize, num;
  ss >> str;
  if (str == "MAXSIZE") {
    ss >> maxsize;
    MAXSIZE = maxsize;
  }
  Solution s;
  int i = 2;
  while (ss >> str) {
    if (str == "LAPSE") // LAPSE <NAME>
    {
      ss >> name;
      // cout << "Line " << i << endl;
      i++;
      // cout << "LAPSE " << endl;
      s.LAPSE(name);
    } else if (str == "KOKUSEN") // KOKUSEN
    {
      // cout << "Line " << i << endl;
      i++;
      // cout << "KOKUSEN " << endl;
      s.KOKUSEN();
    } else if (str == "KEITEIKEN") // KEITEIKEN <NUM>
    {
      // cout << "Line " << i << endl;
      i++;
      ss >> num;
      // cout << "KEITEIKEN " << num << endl;
      s.KEITEIKEN(num);
    } else if (str == "HAND") // HAND
    {
      // cout << "Line " << i << endl;
      i++;
      // cout << "HAND " << endl;
      s.HAND();
    } else if (str == "LIMITLESS") // LIMITLESS <NUM>
    {
      // cout << "Line " << i << endl;
      i++;
      ss >> num;
      // cout << "LIMITLESS " << num << endl;
      s.LIMITLESS(num);
    } else if (str == "CLEAVE") // DOMAIN_EXPANSION
    {
      // cout << "Line " << i << endl;
      i++;
      ss >> num;
      // cout << "CLEAVE " << num << endl;
      s.CLEAVE(num);
    }
  }
  return;
}

int main(int argc, char *argv[]) {
  clock_t begin1 = clock();

  // Your code to measure runtime goes here
  string fileName = "testdsa.txt";
  simulate(fileName);

  clock_t end1 = clock();

  cout << "Execution time: " << (float)(end1 - begin1) / CLOCKS_PER_SEC << " s" << endl;

  return 0;
}
