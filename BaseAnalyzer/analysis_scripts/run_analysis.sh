#!/bin/bash
 
# shell script to run main Hall C analysis


# user command line arguments 
# example: >> run_analysis.sh $1 $2 $3
#username=$1   # string
#age=$2        # int
#full_name=$3  # string

# create a variable to hold the input
read -p "Please enter username: " username
read -p "Please enter age: " age
read -p "Please enter full name: " full_name

echo "Username: $username";
echo "Age: $age";
echo "Full Name: $full_name";

# Check if string is empty using -z. For more 'help test'    
if [[ -z "$username" || -z "$age" ]]; then
   printf '%s\n' "missing input"
   exit 1

else
   # If userInput is not empty show what the user typed in and run ls -l
   printf "You entered %s " "$userInput"
fi

# define command
CMD="root -l -q -b \"test.C(\\\"${username}\\\", ${age}, \\\"${full_name}\\\")\""          

echo "executing command: $CMD" 
eval ${CMD}

