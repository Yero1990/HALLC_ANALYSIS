#include <iostream>
#include "People.h"
#include "Birthday.h"
#include "Friend.h"
using namespace std;

int main() {

  //Birthday birthObj(11,19,1990);
  //birthObj.printDate();
  //People cyero("Carlos Yero", birthObj);
  
  //cyero.printInfo();
  Friend fr;
  Birthday b;
  //fr.printfInfo();
  //fr.printQuest();
  fr.printQuest(b);

  return 0;
}
