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
		// recontruct
		Map<MatrixXd> x_shape(affine_vector.data(), 12, control_vertex.rows());
		MatrixXd a1 = x_shape.middleRows<3>(0);
		MatrixXd a2 = x_shape.middleRows<3>(3);
		MatrixXd a3 = x_shape.middleRows<3>(6);
		MatrixXd b = x_shape.middleRows<3>(9);

		//MatrixXd v_RT,index_closest;

		// optimization iteration
		for (int i = 0; i < iter_optimization; ++i) 
		{
			clock_t optimization_time = clock();
			std::cout << "Estimzation jacoiba matrix...." << std::endl;

			MatrixXd v_RT = ApplyTrans(source_vertex, 
				a1, a2, a3, b, 
				control_vertex, weight_Trans);

			MatrixXd index_closest;
			KNearestSearch(v_RT, target_vertex, 1, index_closest);

			MatrixXd target_closest      = slice(target_vertex, index_closest);
			MatrixXd target_normal_index = slice(target_normal, index_closest);

			MatrixXd diff_v = v_RT - target_closest;

			MatrixXd symbol = (target_normal_index.cwiseProduct(diff_v)).rowwise().sum();

			//estimate jacobia matrix
			MatrixXd jac = Jacobia(
				a1, a2, a3, 
				control_vertex,
				source_vertex,
				weight_smooth,
				weight_Trans,
				coeff_rigid,
				coeff_smooth,
				symbol,
				target_normal_index);

			printf("estimate jacobia matrix done! Time taken: %.4fs\n", 
				(double)(clock() - optimization_time)
				/ CLOCKS_PER_SEC);
		}

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