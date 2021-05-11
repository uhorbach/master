#include <iostream>
#include <cstdlib>
#include <string>
#include <time.h>
#include <fstream>
using namespace std;

int main()
{


	//Zufallszahl erzeugen
	srand(time(NULL));
	 int Z1= rand () %1000000;
	 int Z2= Z1 + rand () %1000000;
    fstream f;
    f.open("run1.mac");
    f.seekp(0, ios_base::beg);
    f <<"/random/setSeeds "<<Z1<<" "<<Z2<<endl;
    f.close();
		int Z3= rand () %1000000;
		int Z4= Z3 + rand () %1000000;
		 fstream g;
		 g.open("run2.mac");
		 g.seekp(0, ios_base::beg);
		 g <<"/random/setSeeds "<<Z3<<" "<<Z4<<endl;
		 g.close();
		int Z5= rand () %1000000;
		int Z6= Z5 + rand () %1000000;
			fstream h;
			h.open("run3.mac");
			h.seekp(0, ios_base::beg);
			h <<"/random/setSeeds "<<Z5<<" "<<Z6<<endl;
			h.close();



return 0;


}
