#ifndef BIRTHDAY_H
#define BIRTHDAY_H

class Birthday
{

 public:
  Birthday(int m, int d, int y);  //constructor
  void printDate();
  static void acceptFriend();


 private:
  int month;
  int day;
  int year;

  //friend class Friend;


  
};


#endif 
