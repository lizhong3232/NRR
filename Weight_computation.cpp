#include "NGP.h"

namespace NRR
{
	// matlab function equal B = A(index,:);
	MatrixXd NGP::slice(MatrixXd A, MatrixXd index) 
	{
		int index_num = index.rows();
		MatrixXd B(index_num, A.cols());
		B.setZero();
		//
		for (int i = 0; i < B.rows(); ++i) 
		{
			B.row(i) = A.row(index(i, 0));
		}

		return B;

	}

	MatrixXd NGP::Weight_smooth_ajacent_geodesic(
		MatrixXd control_vertex,
		MatrixXd r_distance,
		MatrixXd index_control,
		MatrixXd geodesic_table) 
	{
		int control_num = control_vertex.rows();
		MatrixXd weight_tmp(control_num, control_num);
		weight_tmp.setZero();
		//int d_2;
		double d_2,divend_tmp, divide_tmp;
		//
		//r = r(index_control,:);
		//smooth_table = geodesic_table(index_control, :);
		MatrixXd r            = slice(r_distance, index_control);
		MatrixXd smooth_table = slice(geodesic_table, index_control);
		//
		for (int i = 0; i < control_num; ++i) 
		{
			for (int j = 0; j < control_num; ++j)
			{
				if (i != j)
				{
				    divend_tmp = pow(smooth_table(i, j), 2);
					divide_tmp = pow((r(i, 0) + r(j, 0)), 2);
					d_2 = divend_tmp / divide_tmp;

					weight_tmp(i, j) = std::max(0.0, pow((1 - d_2), 3));
				}
			}
		}
		//
		MatrixXd div_tmp = weight_tmp.rowwise().sum();
		MatrixXd n_index_smooth;

		for (int i = 0; i < control_num; ++i) 
		{
			for (int j = 0; j < control_num; ++j) 
			{
				if (div_tmp(i, 0) == 0) 
				{
				  // if unreference vertex, maybe, cause expcetion
					KNearestSearch(control_vertex.row(i),
						control_vertex,
						3,
						n_index_smooth);

					weight_tmp(i, n_index_smooth(0, 0)) = 0.5;
					weight_tmp(i, n_index_smooth(0, 1)) = 0.3;
					weight_tmp(i, n_index_smooth(0, 2)) = 0.2;

					div_tmp(i, 0) = 1.0;
					std::cout << "overflow weight here" << std::endl;
				}
				weight_tmp(i, j) = weight_tmp(i, j) / div_tmp(i, 0);
			}
		}

		/*std::ofstream outfile("weight_smooth.txt");
		if (outfile.is_open())
		{
		outfile << weight_tmp;
		}
		outfile.close();*/

		return weight_tmp;

	}

	// check if is inside
	bool NGP::IsIn(MatrixXd self, double pointer, int& index)
	{
		for (int i = 0; i < self.rows(); ++i) 
		{
			if (self(i, 0) == pointer)
			{	
				index = i;
				return true;
			}
		}
		return false;
	}


	// compute geodesic weight
	MatrixXd NGP::WeightFunc_geodesic(MatrixXd source_vertex,
		MatrixXd control_vertex,
		MatrixXd r_distance,
		MatrixXd index_control,
		MatrixXd geodesic_table)
	{
		int vertex_num  = source_vertex.rows();
		int control_num = control_vertex.rows();

		bool flag_in;
		int index_geo;
		MatrixXd weight_tmp = MatrixXd::Zero(vertex_num, control_num);
		double d_2, divend_tmp, divide_tmp;
		//
		for (int i = 0; i < vertex_num; ++i) 
		{
			for (int j = 0; j < control_num; ++j) 
			{
				flag_in = IsIn(index_control, double(i), index_geo);

				if (flag_in) 
				{
					weight_tmp(i, index_geo) = 1.0;
					continue;
				}
				else 
				{
					divend_tmp = pow(geodesic_table(i,j), 2);
					divide_tmp = pow(r_distance(i, 0), 2);
					d_2 = divend_tmp / divide_tmp;

					weight_tmp(i, j) = std::max(0.0, pow((1 - d_2), 3));
				}
			}
		}

		// normlization
		//
		MatrixXd div_tmp = weight_tmp.rowwise().sum();
		MatrixXd n_index_weight;

		for (int i = 0; i < vertex_num; ++i)
		{
			for (int j = 0; j < control_num; ++j)
			{
				if (div_tmp(i, 0) == 0)
				{
					// if unreference vertex, maybe, cause expcetion
					KNearestSearch(source_vertex.row(i),
						control_vertex,
						3,
						n_index_weight);

					weight_tmp(i, n_index_weight(0, 0)) = 0.5;
					weight_tmp(i, n_index_weight(0, 1)) = 0.3;
					weight_tmp(i, n_index_weight(0, 2)) = 0.2;

					div_tmp(i, 0) = 1.0;
					std::cout << "overflow weight here" << std::endl;
				}
				weight_tmp(i, j) = weight_tmp(i, j) / div_tmp(i, 0);
			}
		}
		/*std::ofstream outfile("weight_geo.txt");
		if (outfile.is_open())
		{
			outfile << weight_tmp;
		}
		outfile.close();*/

		return weight_tmp;

	}

}