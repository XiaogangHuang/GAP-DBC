# GAP-DBC

The C++ code, binary code for "Fast Density-Based Clustering: Geometric Approach".

We uploaded the Household dataset for testing the efficiency of the algorithm. Other datasets could not be uploaded due to the restriction that the uploaded file cannot be larger than 25MB. 

**Remark**: Before downloading and using the binary code, you need to turn off the "Real-time Protection" option in the "Virus and Threat Protection Settings" on your windows system.

## DATA FORMAT
The input dataset needs to be preprocessed into a text file with the following format:

* In each line: its coordinates and the point's id, where the id is an integer in the range [0,  n-1].

* In total, there are n lines where the numbers at each line are separated by a space. 

For example, a 2-dimensional data set containing 4 points:
 
9.3498     56.7408     17.0527     0

9.3501     56.7406     17.6148     1

9.3505     56.7405     18.0835     2

9.3508     56.7404     18.2794     3

