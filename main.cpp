#include "Clustering.h"

/*
 *  Get next word from a string.
 */
void getNextWord(char* str, char* word);

int main(int argc, char** argv)
{
	double radius = 0.0, percent = 1.0;
	int minPts = 0;
	string rootfile = "";
	string suffixes = ".txt";
	string myfile = "";
	char str[20] = { '\0' };
	cout << "Usage:\n"
		<< "  .\\GAP-DBC-KDT [-r radius] [-m minPts] [-p percent] [-df data]\n"
		<< "  where:\n"
		<< "    radius	parameter of GAP-DBC-KDT\n"
		<< "    minPts	parameter of GAP-DBC-KDT\n"
		<< "    percent	parameter of GAP-DBC-KDT (default 1)\n"
		<< "    data	name of file containing data points\n"
		<< "Results are stored in \"data\\result.txt\".\n"
		<< "To run this demo use:\n"
		<< "  .\\GAP-DBC-KDT -r 5000 -m 50 -p 1 -df house\n";
	if (argc > 0)
	{
		int i = 1;
		while (i < argc) {							// read arguments
			if (!strcmp(argv[i], "-m")) {		// -nn option
				minPts = atoi(argv[++i]);				// get number of near neighbors
			}
			else if (!strcmp(argv[i], "-r")) {		// -r option
				sscanf_s(argv[++i], "%lf", &radius);		// get radius
			}
			else if (!strcmp(argv[i], "-p")) {		// -p option
				sscanf_s(argv[++i], "%lf", &percent);		// get percent
			}
			else if (!strcmp(argv[i], "-df")) {		// -df option
				getNextWord(argv[++i], str);
			}
			else {									// illegal syntax
				cerr << "Unrecognized option.\n";
				exit(1);
			}
			i++;
		}
		if (str[0] == '\0')
		{
			cerr << "name of file error.\n";
			exit(1);
		}
		else if (percent > 1 || percent < 0.5)
		{
			cerr << "percent must in [1/2, 1].\n";
			exit(1);
		}
		myfile = str;
	}
	string filename = rootfile + myfile + suffixes;
	printf("======================================================================\n");
	printf("The input parameters are:\n");
	cout << "file:       " << filename << endl;
	printf("radius:     %-5.2f\n", radius);
	printf("minPts:     %-5d\n", minPts);
	printf("percent:    %-5.2f\n", percent);
	printf("======================================================================\n");
	Clustering myClustering;       //Clustering algorithm object declaration.
	myClustering.Init((char*)filename.c_str(), radius, minPts, percent);
	printf("----------------------------------\n");
	printf("Start to run GAP-DBC-KDT.\n");
	printf("----------------------------------\n");
	myClustering.DoDBSCANRecursive();                    //Perform GAP-DBC.
	myClustering.PrintMessage();
	myClustering.WriteToFile("data\\result.txt");      //Save the result.

	return 0;
}

void getNextWord(char* str, char* word) {
	// Jump over all blanks
	while (*str == ' ') {
		str++;
	}

	while (*str != ' ' && *str != '\0') {
		*word = *str;
		str++;
		word++;
	}

	*word = '\0';
}