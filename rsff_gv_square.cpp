#include "NGP.h"

namespace NRR
{
	// non-rigid icp energy function
	VectorXd NGP::rsff_gv_square(
		VectorXd affine_vector,
		MatrixXd control_vertex,
		MatrixXd source_vertex,
		MatrixXd target_vertex,
		MatrixXd target_normal,
		MatrixXd weight,
		MatrixXd weightTrans,
		float coeff_rigid,
		float coeff_smooth
	)
	{
		int control_point_num = control_vertex.rows();
		//------------------------------------------------------------------------
		//	E_rigid
		//------------------------------------------------------------------------
		Map<MatrixXd> x_shape(affine_vector.data(), 12, control_vertex.rows());
		//std::cout << x_shape << std::endl;
		MatrixXd a1 = x_shape.middleRows<3>(0);
		MatrixXd a2 = x_shape.middleRows<3>(3);
		MatrixXd a3 = x_shape.middleRows<3>(6);
		MatrixXd b = x_shape.middleRows<3>(9);

		MatrixXd a1_square = (a1.cwiseProduct(a1)).colwise().sum();
		MatrixXd a2_square = (a2.cwiseProduct(a2)).colwise().sum();
		MatrixXd a3_square = (a3.cwiseProduct(a3)).colwise().sum();

		MatrixXd one_rigid = MatrixXd::Ones(1, a1.cols());

		MatrixXd rigid4 = (one_rigid - a1_square).transpose();
		MatrixXd rigid5 = (one_rigid - a2_square).transpose();
		MatrixXd rigid6 = (one_rigid - a3_square).transpose();

		MatrixXd rigid1 = ((a1.cwiseProduct(a2)).colwise().sum()).transpose();
		MatrixXd rigid2 = ((a1.cwiseProduct(a3)).colwise().sum()).transpose();
		MatrixXd rigid3 = ((a2.cwiseProduct(a3)).colwise().sum()).transpose();

		//cancanate matrix
		MatrixXd E_rigid_square(6 * control_vertex.rows(), rigid1.cols());
		E_rigid_square << rigid1,
			rigid2,
			rigid3,
			rigid4,
			rigid5,
			rigid6;

		E_rigid_square = sqrt(coeff_rigid)*E_rigid_square;

		//std::cout << E_rigid_square << std::endl;
		//------------------------------------------------------------------------
		//	E_smooth
		//------------------------------------------------------------------------
		MatrixXd E_smooth_square = MatrixXd::Zero(control_point_num*control_point_num, 3);
		int smooth_count = 0;
		MatrixXd A_i(3, 3), tmp_smooth;
		MatrixXd x_i(3, 1), x_j(3, 1), b_i(3, 1), b_j(3, 1);
		MatrixXd control_trans = control_vertex.transpose();

		for (int i = 0; i < control_point_num; ++i)
		{
			for (int j = 0; j < control_point_num; ++j)
			{
				A_i << a1.col(i), a2.col(i), a3.col(i);

				x_i = control_trans.col(i);
				x_j = control_trans.col(j);

				b_i = b.col(i);

				b_j = b.col(j);

				tmp_smooth = A_i * (x_j - x_i) + x_i + b_i - (x_j + b_j);

				E_smooth_square.row(smooth_count) = sqrt(coeff_smooth) *sqrt(weight(i,j))*tmp_smooth.transpose();

				++smooth_count;
			}
		}

		Map<MatrixXd> E_smooth_square_reshape(E_smooth_square.data(), control_point_num*control_point_num * 3, 1);
		//std::cout << E_smooth_square_reshape << std::endl;

		//------------------------------------------------------------------------
		//	E_fit
		//------------------------------------------------------------------------
		MatrixXd v_RT = ApplyTrans(
			source_vertex,
			a1,
			a2,
			a3,
			b,
			control_vertex,
			weightTrans);

		MatrixXd index_closest,D;
		bool flag = KNearestSearch_distance(v_RT, target_vertex, 1, index_closest,D);

		MatrixXd thres_ind = findValue(D, threshold, 1);
		/*for (int i = 0; i < thres_ind.rows(); ++i) 
		{
			std::cout << thres_ind.row(i) << std::endl;
		}*/

		//discard some outlier correspondense
		MatrixXd target_closest = slice(target_vertex, index_closest);
		MatrixXd target_normal_index = slice(target_normal,index_closest);

		CMeshO v_RT_mesh;
		vcg::tri::Append<CMeshO, CMeshO>::MeshCopy(v_RT_mesh, source_mesh_vcg);
		for (int i = 0; i < deform_graph_vcg.VN(); ++i) {
			v_RT_mesh.vert[i].P()[0]  = v_RT(i,0);
			v_RT_mesh.vert[i].P()[1]  = v_RT(i,1);
			v_RT_mesh.vert[i].P()[2]  = v_RT(i,2);
		}
		vcg::tri::Allocator<CMeshO>::CompactEveryVector(v_RT_mesh);
		tri::UpdateBounding<CMeshO>::Box(v_RT_mesh);
		MatrixXd v_RT_normal = computeNormal(v_RT_mesh);

		MatrixXd dot_normal = target_normal_index.cwiseProduct(v_RT_normal);
		MatrixXd dot_sum = dot_normal.rowwise().sum();
		MatrixXd thres_normal = findValue(dot_sum, 0.0, 0);

		// thres_normal protection
		if(thres_normal.rows() > (0.7*source_vertex.rows()) )
			thres_normal = findValue(dot_sum, 0.0, 1);

		MatrixXd diff_v = v_RT - target_closest;

		// earase threshold normal inconsistentcy
		//std::cout << diff_v.rows() << std::endl;
		erase_outlier(diff_v, thres_ind);
		erase_outlier(diff_v, thres_normal);

		// caculate Efit1 and Efit2
		//MatrixXd diff_v_transpose = diff_v.transpose();
		Map<MatrixXd> E_fit1_reshape(diff_v.transpose().data(), diff_v.rows() * 3, 1);
		
		MatrixXd E_fit_point = sqrt(alpha_point) * E_fit1_reshape;
		//std::cout << E_fit_point.rows() << " " << E_fit_point.cols() << std::endl;

		MatrixXd E_fit_plane = sqrt(alpha_plane) * ((target_normal_index.cwiseProduct(diff_v)).rowwise().sum()).cwiseAbs();
		/*std::cout << E_fit_plane.rows() << " " << E_fit_plane.cols() << std::endl;
		std::cout << E_rigid_square.rows() << " " << E_rigid_square.cols() << std::endl;
		std::cout << E_smooth_square.rows() << " " << E_smooth_square.cols() << std::endl;*/

		VectorXd energy((control_point_num * 6 + control_point_num*control_point_num * 3 + source_vertex.rows() * 4));

		energy << E_rigid_square,
			E_smooth_square_reshape,
			E_fit_point,
			E_fit_plane;

		/*std::cout << E_rigid_square.transpose()*E_rigid_square << std::endl;
		std::cout << E_smooth_square_reshape.transpose()*E_smooth_square_reshape << std::endl;

		std::cout << energy.transpose()*energy << std::endl;*/

		return energy;

	}

	// index is an vector, main is N*3 matrix,
	// output the index in vector of that main matrix is 
	// equal to zero
	void NGP::erase_outlier(MatrixXd& main, MatrixXd index) 
	{
		int current_index;
		for (int i = 0; i<index.rows(); ++i) 
		{
			current_index = index(i, 0);
			/*std::cout << main.row(i) << std::endl;
			std::cout << current_index << std::endl;*/
			
			main(current_index, 0) = 0.0;
			main(current_index, 1) = 0.0;
			main(current_index, 2) = 0.0;

			
		}
	}

	// find index where c_matrix element bigger or smaller to value, 
	// flag==0, means "<", flag==1 means ">"
	MatrixXd NGP::findValue(MatrixXd c_matrix, float value, int flag)
	{
		MatrixXd result = MatrixXd::Zero(c_matrix.rows(), 1);
		int f_count = 0;
		for (int i = 0; i < c_matrix.rows(); ++i) 
		{
			// if find element in cmatrix smaller than value
			if (flag == 0) 
			{	
				if (c_matrix(i, 0) < value) 
				{
					/*result.row(f_count) = c_matrix.row(i);*/
					result(f_count,0) = i;
					++f_count;
				}
				
			}

			// if find element in cmatrix bigger than value
			else if (flag == 1) 
			{
				if (c_matrix(i, 0) > value)
				{
					/*result.row(f_count) = c_matrix.row(i);*/
					result(f_count, 0) = i;
					++f_count;
				}
			}
		}

		MatrixXd result_final = result.topRows(f_count);

		return result_final;
	}


	// apply affine vector to source_vertex to new position
	MatrixXd NGP::ApplyTrans(
		MatrixXd source_vertex,
		MatrixXd a1,
		MatrixXd a2,
		MatrixXd a3,
		MatrixXd b,
		MatrixXd control_vertex,
		MatrixXd weightTrans) 
	{
		int source_length  = source_vertex.rows();
		int control_length = control_vertex.rows();

		MatrixXd source_trans  = source_vertex.transpose();
		MatrixXd control_trans = control_vertex.transpose();

		MatrixXd result = MatrixXd::Zero(3, source_length);
		MatrixXd A_i(3, 3);
		MatrixXd x_i,b_i;
		MatrixXd value;

		for (int i = 0; i < control_length; ++i) 
		{
			A_i << a1.col(i), a2.col(i), a3.col(i);
			x_i = control_trans.col(i).replicate(1, source_length);
			b_i = b.col(i).replicate(1, source_length);
			
			value = (weightTrans.col(i).transpose().replicate(3, 1)).
					cwiseProduct(A_i *(source_trans - x_i) + x_i + b_i);

			result += value;
			
		}

		MatrixXd result_final = result.transpose();
		//std::cout << result_final << std::endl;

		return result_final;
		
	}

	


}