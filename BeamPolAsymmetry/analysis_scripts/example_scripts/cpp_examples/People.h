#ifndef PEOPLE_H
#define PEOPLE_H

#include <string>
#include "Birthday.h"

using namespace std;

class Friend;

class People
{

  friend class Friend;
  
 public:
  People(string x, Birthday bo);
  void printInfo();
 private:
  string name;
  Birthday dateOfBirth;
  
};


#endif 
