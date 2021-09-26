#include "Friend.h"
#include "Birthday.h"
Friend::Friend(){
  cout << "Initializing Constructor Friend " << endl;
  
}

void Friend::printfInfo(){
  cout << "Im a friendly class " << endl;

}

void Friend::printQuest(Birthday b){
  cout << "Thanks for being friend " << endl;
  //Birthday bd(11,19,1990);
  b.acceptFriend();
}
