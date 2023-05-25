### Some Extra Helpful Commands 

To mass remove files. BE VERY CAREFULLY: 

    ls -1 */*/*/*.out | awk -F /Z '{print "rm " $1}' > BadFilesToDelete.txt