#include "TRint.h"
#include "digitalFilterOptimization.hh"
#include <iostream>

using namespace std;
/*void signal_handler(int signal_number){

	if(signal_number == SIGINT)
		std::cout<<"SIGINT signal = "<<signal_number<<" received!!!"<<std::endl;
	else if(signal_number == SIGTERM)
		std::cout<<"SIGTERM signal = "<<signal_number<<" received!!!"<<std::endl;

	GUser* a = GUser::g_instance;
	if(a){
		a->EndUser();
		delete (a);
		a = NULL;
	}
	gROOT->Reset();	
	exit(signal_number);

}
*/
int main(int argc, char** argv){
TRint *app = new TRint("",&argc, argv);
UShort_t par[6]={200,500, 50,50,200,50};

//std::signal(SIGINT,signal_handler);
//	std::signal(SIGTERM,signal_handler);
digitalFilterOptimization * a = new digitalFilterOptimization( par);

a->optimize_trapezoidal_filter_parameters(argv[1], stod(argv[2]), stod(argv[3]));

delete a;
app->Run();
return 0;
}
