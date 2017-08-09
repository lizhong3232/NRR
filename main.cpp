#include "NGP.h"
#include <iostream>
#include <fstream>

int main(int argc, char*argv[])
{
	/*std::string source_mesh_name = "source_hand.ply";
	std::string target_mesh_name = "target_hand.ply";*/

	std::string source_mesh_name = "0070.ply";
	std::string target_mesh_name = "0071.ply";


	NRR::NGP test;
	//NRR::NGP test();
	test.initial(source_mesh_name,
				 target_mesh_name);
	test.optimization();

	//system("pause");
	return 0;
	
}