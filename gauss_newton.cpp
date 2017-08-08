#include "NGP.h"

namespace NRR
{
	// gauss_newton function
	VectorXd NGP::gauss_newton(
		VectorXd affine_vector,
		MatrixXd control_vertex,
		MatrixXd source_vertex,
		MatrixXd target_vertex,
		MatrixXd target_normal,
		MatrixXd weight_smooth,
		MatrixXd weight_Trans,
		float coeff_rigid,
		float coeff_smooth) 
	{
		VectorXd initial_sqrt = rsff_gv_square(
			affine_vector,
			control_vertex,
			 source_vertex,
			target_vertex,
			target_normal,
			weight_smooth,
			weight_Trans,
			coeff_rigid,
			coeff_smooth
		);

		double energy_previous = initial_sqrt.transpose()*initial_sqrt;
		std::cout << "Initial cost is " << energy_previous << std::endl;
		

		return affine_vector;
	}


	VectorXd NGP::initialAffineVector(int control_num)
	{
		Eigen::VectorXd affineVector(12 * control_num);

		#pragma omp parallel for
		for (int nodeIdx = 0; nodeIdx < control_num; ++nodeIdx)
		{
			affineVector(12 * nodeIdx + 0) = 1.0f;
			affineVector(12 * nodeIdx + 1) = 0.0f;
			affineVector(12 * nodeIdx + 2) = 0.5f;

			affineVector(12 * nodeIdx + 3) = 0.0f;
			affineVector(12 * nodeIdx + 4) = 1.0f;
			affineVector(12 * nodeIdx + 5) = 0.0f;

			affineVector(12 * nodeIdx + 6) = 0.0f;
			affineVector(12 * nodeIdx + 7) = 0.0f;
			affineVector(12 * nodeIdx + 8) = 1.0f;

			affineVector(12 * nodeIdx + 9) = 0.0f;
			affineVector(12 * nodeIdx + 10) = 0.0f;
			affineVector(12 * nodeIdx + 11) = 0.0f;
		}

		return affineVector;
	}

}