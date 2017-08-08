#include "NGP.h"

namespace NRR 
{
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

	// Run jacobia matrix
	MatrixXd NGP::Jacobia(
		MatrixXd a1, 
		MatrixXd a2, 
		MatrixXd a3,
		MatrixXd control_point, 
		MatrixXd source_obj, 
		MatrixXd weight, 
		MatrixXd weightTrans,
		float coeff_rigid, 
		float coeff_smooth, 
		MatrixXd symbol, 
		MatrixXd target_normal_index) 
	{

		clock_t optimization_time = clock();
		unsigned long int m_size = control_point.rows() * 6 + pow(control_point.rows(), 2) * 3 + source_obj.rows() * 3 + source_obj.rows();

		std::cout << source_obj.rows() << std::endl;

		int variable_num = 12;

		// initial output matrix
		MatrixXd jac = MatrixXd::Zero(m_size, variable_num * control_point.rows());

		MatrixXd x1 = a1.row(0); MatrixXd x4 = a2.row(0); MatrixXd x7 = a3.row(0);
		MatrixXd x2 = a1.row(1); MatrixXd x5 = a2.row(1); MatrixXd x8 = a3.row(1);
		MatrixXd x3 = a1.row(2); MatrixXd x6 = a2.row(2); MatrixXd x9 = a3.row(2);

		///
		int index1, index2, index3, index4, index5, index6, index7, index8, index9;

		/**************************************************************************/
		/************************* E_rigid term ***********************************/
		/**************************************************************************/
		int seq = 1;
		int row_index;
		//---- - rigid part  1----------
#pragma omp parallel for
		for (int i = 0; i < control_point.rows(); ++i) {

			index1 = i * variable_num + 0;
			index2 = i * variable_num + 1;
			index3 = i * variable_num + 2;
			index4 = i * variable_num + 3;
			index5 = i * variable_num + 4;
			index6 = i * variable_num + 5;

			/*std::cout << index1 << std::endl;
			std::cout << index2 << std::endl;*/

			jac(i, index1) = x4(0, i);
			jac(i, index2) = x5(0, i);
			jac(i, index3) = x6(0, i);

			jac(i, index4) = x1(0, i);
			jac(i, index5) = x2(0, i);
			jac(i, index6) = x3(0, i);
		}
		/************* ---- rigid part 2 *********** ***********/
		//std::cout << jac.topRows(8) << std::endl;
		seq++;
#pragma omp parallel for
		for (int i = 0; i < control_point.rows(); ++i) {
			row_index = control_point.rows() * (seq - 1) + i;

			index1 = i * variable_num + 0;
			index2 = i * variable_num + 1;
			index3 = i * variable_num + 2;
			index7 = i * variable_num + 6;
			index8 = i * variable_num + 7;
			index9 = i * variable_num + 8;
			//%------ -
			jac(row_index, index1) = x7(i);
			jac(row_index, index2) = x8(i);
			jac(row_index, index3) = x9(i);

			jac(row_index, index7) = x1(i);
			jac(row_index, index8) = x2(i);
			jac(row_index, index9) = x3(i);
		}
		/************* ---- rigid part 3  ********** ***********/
		seq++;
#pragma omp parallel for
		for (int i = 0; i < control_point.rows(); ++i) {
			row_index = control_point.rows() * (seq - 1) + i;

			index4 = i * variable_num + 3;
			index5 = i * variable_num + 4;
			index6 = i * variable_num + 5;
			index7 = i * variable_num + 6;
			index8 = i * variable_num + 7;
			index9 = i * variable_num + 8;
			//%------ -
			jac(row_index, index4) = x7(i);
			jac(row_index, index5) = x8(i);
			jac(row_index, index6) = x9(i);

			jac(row_index, index7) = x4(i);
			jac(row_index, index8) = x5(i);
			jac(row_index, index9) = x6(i);
		}

		/**************** % ----- rigid part  2 ---------- ******************/
		seq++;
#pragma omp parallel for
		for (int i = 0; i < control_point.rows(); ++i) {
			row_index = control_point.rows() * (seq - 1) + i;

			index1 = i * variable_num + 0;
			index2 = i * variable_num + 1;
			index3 = i * variable_num + 2;
			//%------ -
			jac(row_index, index1) = -2 * x1(i);
			jac(row_index, index2) = -2 * x2(i);
			jac(row_index, index3) = -2 * x3(i);

		}
		/******************* rigid 2.2 ***********/
		seq++;
#pragma omp parallel for
		for (int i = 0; i < control_point.rows(); ++i) {
			row_index = control_point.rows() * (seq - 1) + i;

			index4 = i * variable_num + 3;
			index5 = i * variable_num + 4;
			index6 = i * variable_num + 5;
			//%------ -
			jac(row_index, index4) = -2 * x4(i);
			jac(row_index, index5) = -2 * x5(i);
			jac(row_index, index6) = -2 * x6(i);

		}
		/******************* rigid 2.3 ***********/
		seq++;
#pragma omp parallel for
		for (int i = 0; i < control_point.rows(); ++i) {
			row_index = control_point.rows() * (seq - 1) + i;

			index7 = i * variable_num + 6;
			index8 = i * variable_num + 7;
			index9 = i * variable_num + 8;
			//%------ -
			jac(row_index, index7) = -2 * x7(i);
			jac(row_index, index8) = -2 * x8(i);
			jac(row_index, index9) = -2 * x9(i);

		}

		jac.middleRows(0, 6 * control_point.rows()) *= sqrt(coeff_rigid);
		/**************************************************************************/
		/************************* E_smooth term ***********************************/
		/**************************************************************************/
		seq = seq + 1;

		MatrixXd control_trans = control_point.transpose();
		int control_point_num = control_point.rows();

		int row_smooth = control_point.rows() * (seq - 1);
		int smooth_start = control_point.rows() * (seq - 1);
		int smooth_end = control_point.rows() * (seq - 1) + control_point_num ^ 2 * 3;
		//
		int index1i, index2i, index3i, index4i, index5i, index6i, index7i, index8i, index9i, index10i, index11i, index12i;
		int	index1j, index2j, index3j, index4j, index5j, index6j, index7j, index8j, index9j, index10j, index11j, index12j;
		MatrixXd x_i, x_j, D;
		double Dx, Dy, Dz;
		float eff;
#pragma omp parallel for
		for (int i = 0; i < control_point_num; ++i) {
			index1i = i * variable_num + 0;
			index2i = i * variable_num + 1;
			index3i = i * variable_num + 2;
			index4i = i * variable_num + 3;
			index5i = i * variable_num + 4;
			index6i = i * variable_num + 5;
			index7i = i * variable_num + 6;
			index8i = i * variable_num + 7;
			index9i = i * variable_num + 8;
			index10i = i * variable_num + 9;
			index11i = i * variable_num + 10;
			index12i = i * variable_num + 11;
			for (int j = 0; j < control_point_num; ++j) {
				// preprocess initial variable
				x_i = control_trans.col(i);
				x_j = control_trans.col(j);
				D = x_j - x_i;
				Dx = D(0, 0);
				Dy = D(1, 0);
				Dz = D(2, 0);
				eff = sqrt(weight(i, j));
				// index
				index10j = j * variable_num + 9;
				index11j = j * variable_num + 10;
				index12j = j * variable_num + 11;
				for (int k = 0; k < 3; ++k) {
					if (k == 0) {
						jac(row_smooth, index1i) = eff * Dx;
						jac(row_smooth, index4i) = Dy * eff;
						jac(row_smooth, index7i) = Dz * eff;
						jac(row_smooth, index10i) = 1 * eff;
						jac(row_smooth, index10j) = -1 * eff;
					}

					else if (k == 1) {
						jac(row_smooth, index2i) = Dx * eff;
						jac(row_smooth, index5i) = Dy * eff;
						jac(row_smooth, index8i) = Dz * eff;
						jac(row_smooth, index11i) = 1 * eff;
						jac(row_smooth, index11j) = -1 * eff;
					}

					else if (k == 2) {
						jac(row_smooth, index3i) = Dx * eff;
						jac(row_smooth, index6i) = Dy * eff;
						jac(row_smooth, index9i) = Dz * eff;
						jac(row_smooth, index12i) = 1 * eff;
						jac(row_smooth, index12j) = -1 * eff;
					}
					row_smooth++;
				}
			}
		}

		jac.middleRows(smooth_start, pow(control_point.rows(), 2) * 3) *= sqrt(coeff_smooth);

		printf("estimate jacobia matrix done! Time taken: %.4fs\n",
			(double)(clock() - optimization_time)
			/ CLOCKS_PER_SEC);
		optimization_time = clock();
		/**************************************************************************/
		/************************* E_fit term ***********************************/
		/**************************************************************************/
		//%% E fit_1 term
		float alpha_point = 0.1;
		float alpha_plane = 1.0;
		float alpha;
		for (int i = 0; i < source_obj.rows(); ++i) {

			alpha = sqrt(alpha_point);
			for (int k = 0; k < 3; ++k) {
				for (int j = 0; j < control_point_num; ++j) {

					D = source_obj.row(i).transpose() - control_point.row(j).transpose();
					Dx = D(0, 0); Dy = D(1, 0); Dz = D(2, 0);
					index1j = j * variable_num + 0;
					index2j = j * variable_num + 1;
					index3j = j * variable_num + 2;
					index4j = j * variable_num + 3;
					index5j = j * variable_num + 4;
					index6j = j * variable_num + 5;
					index7j = j * variable_num + 6;
					index8j = j * variable_num + 7;
					index9j = j * variable_num + 8;
					index10j = j * variable_num + 9;
					index11j = j *variable_num + 10;
					index12j = j *variable_num + 11;
					///
					if (k == 0) {
						jac(row_smooth, index1j) = alpha*weightTrans(i, j) * Dx;
						jac(row_smooth, index4j) = alpha*weightTrans(i, j) * Dy;
						jac(row_smooth, index7j) = alpha*weightTrans(i, j) * Dz;
						jac(row_smooth, index10j) = alpha*weightTrans(i, j) * 1;
					}
					else if (k == 1) {
						jac(row_smooth, index2j) = alpha*weightTrans(i, j) * Dx;
						jac(row_smooth, index5j) = alpha*weightTrans(i, j) * Dy;
						jac(row_smooth, index8j) = alpha*weightTrans(i, j) * Dz;
						jac(row_smooth, index11j) = alpha*weightTrans(i, j) * 1;
					}
					else if (k == 2) {
						jac(row_smooth, index3j) = alpha*weightTrans(i, j) * Dx;
						jac(row_smooth, index6j) = alpha*weightTrans(i, j) * Dy;
						jac(row_smooth, index9j) = alpha*weightTrans(i, j) * Dz;
						jac(row_smooth, index12j) = alpha*weightTrans(i, j) * 1;
					}
				}
				row_smooth++;
			}
		}
		//%% E fit_2 term
		MatrixXd normal_current;
		double nx, ny, nz;
		int flag;
		for (int i = 0; i < source_obj.rows(); ++i) {

			alpha = sqrt(alpha_plane);
			normal_current = target_normal_index.row(i);
			nx = normal_current(0, 0); ny = normal_current(0, 1); nz = normal_current(0, 2);
			// flag to abs
			if (symbol(i, 0) < 0)
				flag = -1;
			else
				flag = 1;
			//
			for (int j = 0; j < control_point_num; ++j) {
				D = source_obj.row(i).transpose() - control_point.row(j).transpose();
				Dx = D(0, 0); Dy = D(1, 0); Dz = D(2, 0);
				/***********************************/
				index1j = j * variable_num + 0;
				index2j = j * variable_num + 1;
				index3j = j * variable_num + 2;
				index4j = j * variable_num + 3;
				index5j = j * variable_num + 4;
				index6j = j * variable_num + 5;
				index7j = j * variable_num + 6;
				index8j = j * variable_num + 7;
				index9j = j * variable_num + 8;
				index10j = j * variable_num + 9;
				index11j = j *variable_num + 10;
				index12j = j *variable_num + 11;
				/*********************/
				jac(row_smooth, index1j) = flag*nx*alpha*weightTrans(i, j) * Dx;
				jac(row_smooth, index2j) = flag*ny*alpha*weightTrans(i, j) * Dx;
				jac(row_smooth, index3j) = flag*nz*alpha*weightTrans(i, j) * Dx;

				jac(row_smooth, index4j) = flag*nx*alpha*weightTrans(i, j) * Dy;
				jac(row_smooth, index5j) = flag*ny*alpha*weightTrans(i, j) * Dy;
				jac(row_smooth, index6j) = flag*nz*alpha*weightTrans(i, j) * Dy;

				jac(row_smooth, index7j) = flag* nx*alpha*weightTrans(i, j) * Dz;
				jac(row_smooth, index8j) = flag*ny*alpha*weightTrans(i, j) * Dz;
				jac(row_smooth, index9j) = flag*nz*alpha*weightTrans(i, j) * Dz;

				jac(row_smooth, index10j) = flag*nx*alpha*weightTrans(i, j) * 1;
				jac(row_smooth, index11j) = flag*ny*alpha*weightTrans(i, j) * 1;
				jac(row_smooth, index12j) = flag*nz*alpha*weightTrans(i, j) * 1;
			}
			row_smooth++;
		}

		printf("estimate jacobia matrix done! Time taken: %.4fs\n",
			(double)(clock() - optimization_time)
			/ CLOCKS_PER_SEC);

		std::string output = "jac.txt";
		saveMatrix(output, jac, 1);
	}


}