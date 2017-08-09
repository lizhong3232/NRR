#include "NGP.h"

namespace NRR
{
	void NGP::optimization() 
	{
		// vector initialization
		//VectorXd affineVector = initialAffineVector(int control_num)

		// optimization iteration
		int iter = 1;
		double coeff_rigid, coeff_smooth;

		while (1) 
		{
			coeff_rigid = coeff_rigid_tmp * pow(0.5 ,(iter - 1));
			coeff_smooth = coeff_smooth_tmp *pow(0.5, (iter - 1));
			// set break condition
			if (coeff_rigid <= 49)
				break;

			// optimization iteration main 
			Eigen::MatrixXd affineVector_out = gauss_newton(affineVector,
												control_vertex,
												source_vertex,
												target_vertex,
												target_normal,
												weight_smooth,
												weight_geodesic,
												coeff_rigid,
												coeff_smooth);
			affineVector = affineVector_out;
			++iter;

		}
	}

}