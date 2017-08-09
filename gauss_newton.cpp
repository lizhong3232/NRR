#include "NGP.h"


namespace NRR
{
	template<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
	inline bool loadMatrix(std::string filename, Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>& m)
	{

		// General structure
		// 1. Read file contents into vector<double> and count number of lines
		// 2. Initialize matrix
		// 3. Put data in vector<double> into matrix

		std::ifstream input(filename.c_str());
		if (input.fail())
		{
			std::cerr << "ERROR. Cannot find file '" << filename << "'." << std::endl;
			m = Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>(0, 0);
			return false;
		}
		std::string line;
		Scalar d;

		std::vector<Scalar> v;
		int n_rows = 0;
		while (getline(input, line))
		{
			++n_rows;
			std::stringstream input_line(line);
			while (!input_line.eof())
			{
				input_line >> d;
				v.push_back(d);
			}
		}
		input.close();

		int n_cols = v.size() / n_rows;
		m = Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>(n_rows, n_cols);

		for (int i = 0; i<n_rows; i++)
			for (int j = 0; j<n_cols; j++)
				m(i, j) = v[i*n_cols + j];

		return true;
	}

	template<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
	inline bool saveMatrix(std::string filename, Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> matrix, bool overwrite)
	{
		//if (boost::filesystem::exists(filename))
		// {
		//   if (!overwrite)
		//    {
		//      // File exists, but overwriting is not allowed. Abort.
		//       std::cerr << "File '" << filename << "' already exists. Not saving data." << std::endl;
		//      return false;
		//    }
		// }


		std::ofstream file;
		file.open(filename.c_str());
		if (!file.is_open())
		{
			std::cerr << "Couldn't open file '" << filename << "' for writing." << std::endl;
			return false;
		}

		file << std::fixed;
		file << matrix;
		file.close();
		return true;

	}

	// gauss_newton function
	Eigen::VectorXd NGP::gauss_newton(
		Eigen::VectorXd affine_vector,
		Eigen::MatrixXd control_vertex,
		Eigen::MatrixXd source_vertex,
		Eigen::MatrixXd target_vertex,
		Eigen::MatrixXd target_normal,
		Eigen::MatrixXd weight_smooth,
		Eigen::MatrixXd weight_Trans,
		float coeff_rigid,
		float coeff_smooth) 
	{
		//loadMatrix("target_normal_matlab.txt", target_normal);

		Eigen::VectorXd initial_sqrt = rsff_gv_square(
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
		
		Eigen::VectorXd affine_vector_current = affine_vector;
		// recontruct
		Eigen::Map<Eigen::MatrixXd> x_shape(affine_vector_current.data(), 12, control_vertex.rows());
		Eigen::MatrixXd a1 = x_shape.middleRows<3>(0);
		Eigen::MatrixXd a2 = x_shape.middleRows<3>(3);
		Eigen::MatrixXd a3 = x_shape.middleRows<3>(6);
		Eigen::MatrixXd b = x_shape.middleRows<3>(9);

		//MatrixXd v_RT,index_closest;
		clock_t optimization_time;
		// optimization iteration

		for (int i = 0; i < iter_optimization; ++i) 
		{
			optimization_time = clock();
			std::cout << "Estimzation jacoiba matrix...." << std::endl;

			Eigen::MatrixXd v_RT = ApplyTrans(source_vertex,
				a1, a2, a3, b, 
				control_vertex, weight_Trans);

			/*std::cout << a1 << std::endl;
			std::cout << a2 << std::endl;
			std::cout << a3 << std::endl;
			std::cout << b << std::endl;*/


			Eigen::MatrixXd index_closest;
			KNearestSearch(v_RT, target_vertex, 1, index_closest);

			/*std::cout << index_closest.rows() << std::endl;
			std::cout << index_closest.row(0) << std::endl;
			std::cout << index_closest.row(1) << std::endl;
			std::cout << index_closest.row(2) << std::endl;
			std::cout << index_closest.row(3) << std::endl;
			std::cout << index_closest.row(4) << std::endl;*/

			Eigen::MatrixXd target_closest      = slice(target_vertex, index_closest);
			Eigen::MatrixXd target_normal_index = slice(target_normal, index_closest);

			Eigen::MatrixXd diff_v = v_RT - target_closest;

			Eigen::MatrixXd symbol = (target_normal_index.cwiseProduct(diff_v)).rowwise().sum();

			//std::cout << symbol.rows() << " " << symbol.cols() << std::endl;
			/*for (int i = 0; i < symbol.rows(); ++i) {
				std::cout << symbol.row(i) << std::endl;
			}
*/
			//estimate jacobia matrix
			Eigen::MatrixXd jac = Jacobia(
				a1, a2, a3, 
				control_vertex,
				source_vertex,
				weight_smooth,
				weight_Trans,
				coeff_rigid,
				coeff_smooth,
				symbol,
				target_normal_index);

			/*std::cout << a1 << std::endl;
			std::cout << a2 << std::endl;
			std::cout << a3 << std::endl;
			std::cout << b << std::endl;*/


			Eigen::SparseMatrix<double> jac_sparse = jac.sparseView();


			printf("estimate jacobia matrix done! Time taken: %.4fs\n", 
				(double)(clock() - optimization_time)
				/ CLOCKS_PER_SEC);

			optimization_time = clock();
			// solve least square
			std::cout << "solve least square...." << std::endl;
			VectorXd F = rsff_gv_square(
				affine_vector_current,
				control_vertex,
				source_vertex,
				target_vertex,
				target_normal,
				weight_smooth,
				weight_Trans,
				coeff_rigid,
				coeff_smooth
			);

			//std::cout << F.transpose() * F << std::endl;

			/*std::cout << F.rows() << std::endl;
			std::cout << jac.rows() << std::endl;*/

			SparseMatrix<double> A = jac_sparse.transpose() * jac_sparse;
			

			VectorXd b_vec = -1 * jac.transpose() * F;

			/*std::ofstream F1("F.txt");
			F1 << F;
			F1.close();
*/
			

			//VectorXd affine_delta;

			Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> CholeskySolver;
			CholeskySolver.compute(A);
			Eigen::VectorXd AffineVector_delta = CholeskySolver.solve(b_vec);

			printf("solve least square done! Time taken: %.4fs\n",
				(double)(clock() - optimization_time)
				/ CLOCKS_PER_SEC);

			// cacluate new x_output
			Eigen::VectorXd nextAffineVector = AffineVector_delta + affine_vector_current;
			VectorXd F_out = rsff_gv_square(
				nextAffineVector,
				control_vertex,
				source_vertex,
				target_vertex,
				target_normal,
				weight_smooth,
				weight_Trans,
				coeff_rigid,
				coeff_smooth
			);

			double energy_out = F_out.transpose() * F_out;
			std::cout << "coef_rigid is " << coeff_rigid << ", after num "<< i <<
				" iteration cost is " << energy_out << std::endl;

			// adjust ming
			double compare_value = abs(energy_out - energy_previous) / energy_previous;
			energy_previous = energy_out;

			std::cout << "coef_rigid is " << coeff_rigid << ", after num "<< i <<
				" iteration diff is " << compare_value << std::endl;

			if (compare_value < 0.005 && compare_value > 0.00005) 
			{
				std::cout << "stop because threshold set" << std::endl;
				break;
			}

			affine_vector_current = nextAffineVector;

			// test result after each iteration
			Eigen::Map<Eigen::MatrixXd> x_shape_inside(nextAffineVector.data(), 12, control_vertex.rows());
		    a1 = x_shape_inside.middleRows<3>(0);
			a2 = x_shape_inside.middleRows<3>(3);
			a3 = x_shape_inside.middleRows<3>(6);
			b  = x_shape_inside.middleRows<3>(9);

			/*std::cout << a1 << std::endl;
			std::cout << a2 << std::endl;
			std::cout << a3 << std::endl;
			std::cout << b << std::endl;*/

			v_RT = ApplyTrans(source_vertex, a1, a2, a3, b, control_vertex, weight_Trans);

			// write to file

		}

		return affine_vector;
	}


	Eigen::VectorXd NGP::initialAffineVector(int control_num)
	{
		Eigen::VectorXd affineVector(12 * control_num);

		#pragma omp parallel for
		for (int nodeIdx = 0; nodeIdx < control_num; ++nodeIdx)
		{
			affineVector(12 * nodeIdx + 0) = 1.0f;
			affineVector(12 * nodeIdx + 1) = 0.0f;
			affineVector(12 * nodeIdx + 2) = 0.0f;

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