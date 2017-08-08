#include "NGP.h"
#include <iostream>
#include <fstream>

int main(int argc, char*argv[])
{
	std::string source_mesh_name = "source_hand.ply";
	std::string target_mesh_name = "target_hand.ply";


	NRR::NGP test;
	//NRR::NGP test();
	test.initial(source_mesh_name,
				 target_mesh_name);
	test.optimization();

	//system("pause");
	return 0;
	
}