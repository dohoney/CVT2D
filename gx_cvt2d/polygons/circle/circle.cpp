#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;
#define NUM_EDGE 100



int main(int argc, char* argv[]) {
	if (argc != 2)
	{
		cout << "circle [radius]" << endl;
		return 0;
	}

	double radius = atof(argv[1]);

	for(int i=0; i<NUM_EDGE; ++i) {
		double angle = i * 2 * M_PI / NUM_EDGE;
		cout << "v " << radius * cos(angle) << " ";
		cout << radius * sin(angle) << endl;
	}

	for(int i=1; i<NUM_EDGE; ++i) {
		cout << "e " << i << " " << i+1 << endl;
	}
	cout << "e " << NUM_EDGE << " " << 1 << endl;

	return 0;
}