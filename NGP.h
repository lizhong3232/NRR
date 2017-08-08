#ifndef _NGP_H_
#define _NGP_H_

#include <iostream>
#include <fstream>
#include "wrap\io_trimesh\import_ply.h"
#include "wrap\io_trimesh\export_ply.h"
#include <Eigen/Dense>
#include "quadric_tex_simp.h"
#include <algorithm>
#include <opencv2/opencv.hpp>

using namespace vcg;
using namespace Eigen;


namespace NRR 
{
	// non-rigid registration 
	class NGP 
	{
		public:
		struct NGP_mesh {// own define mesh structure

			MatrixXd vertice;
			MatrixXd face;
			MatrixXd uv;

		};

		public:
			NGP();
			void initial(std::string source_mesh_name, 
						 std::string target_mesh_name);
			//source or target
			void read_ply(std::string mesh_name, std::string type);
			void write_obj(MatrixXd vertex);

			void optimization();


		private:
			bool KNearestSearch(Eigen::MatrixXd source,
				Eigen::MatrixXd target,
				int k,
				Eigen::MatrixXd &n_index);

			// k nearest search with distance enable
			bool KNearestSearch_distance(Eigen::MatrixXd source,
				Eigen::MatrixXd target,
				int k,
				Eigen::MatrixXd &n_index, 
				Eigen::MatrixXd &distance_index);


			void extract_deform_nodes(CMeshO &test_mesh_vcg, int control_face_num);

			bool Simplify(CMeshO &cm, int targetfacenum, float QualityThr,
				float Extratcoordw, bool OptimalPlacement,
				bool PreserveBoundary, float BoundaryWeight,
				bool PlanarQuadric, bool PreserveNormal);

			void updateDataMask(CMeshO &cm, int neededDataMask);
			void QuadricTexSimplification(CMeshO &m, int  TargetFaceNum, 
				                           bool Selected, 
				                           tri::TriEdgeCollapseQuadricTexParameter &pp);
			// compute geodesic distance
			void ComputeGeodesic(std::string source, std::vector<int> index, int num, MatrixXd &geodesic_table,
				MatrixXd &r_distance, int nearest_geo);
			void NGP::Geodesic_pipeline(std::string source, MatrixXd index_control, int num_neighbor,
										MatrixXd& geodesic_table,MatrixXd& r_distance);

			// compute weight
			MatrixXd Weight_smooth_ajacent_geodesic(MatrixXd control_vertex,
												MatrixXd r_distance,
												MatrixXd index_control,
												MatrixXd geodesic_table);

			MatrixXd weight_smooth, weight_geodesic;

			MatrixXd WeightFunc_geodesic(MatrixXd source_vertex,
									 MatrixXd control_vertex,
								     MatrixXd r_distance,
									 MatrixXd index_control,
								     MatrixXd geodesic_table);
			// compute normal
			MatrixXd computeNormal(CMeshO& target_source);

			// initial affineVecotr
			VectorXd initialAffineVector(int control_num);
			VectorXd affineVector;

			// non-rigid icp function
			VectorXd rsff_gv_square(
				VectorXd affine_vector,
				MatrixXd control_vertex,
				MatrixXd source_vertex,
				MatrixXd target_vertex,
				MatrixXd target_normal,
				MatrixXd weight,
				MatrixXd weightTrans,
				float coeff_rigid,
				float coeff_smooth
			);

			// find index where c_matrix element bigger or smaller to value, 
			// flag==0, means "<", flag==1 means ">"
			MatrixXd findValue(MatrixXd c_matrix, float value, int flag);

			// index is an vector, main is N*3 matrix,
			// output the index in vector of that main matrix is 
			// equal to zero
			void erase_outlier(MatrixXd &main, MatrixXd index);

			// E_fit coefficient
			float alpha_point, alpha_plane;
			float threshold;// threshold for correspondense

			// apply affine vector to source_vertex to new position
			MatrixXd ApplyTrans(
				MatrixXd source_vertex,
				MatrixXd a1,
				MatrixXd a2,
				MatrixXd a3,
				MatrixXd b,
				MatrixXd control_vertex,
				MatrixXd weightTrans);

			// Gauss_newton function
			VectorXd gauss_newton(VectorXd affine_vector,
				MatrixXd control_vertex,
				MatrixXd source_vertex,
				MatrixXd target_vertex,
				MatrixXd target_normal,
				MatrixXd weight_smooth,
				MatrixXd weight_Trans,
				float coeff_rigid,
				float coeff_smooth);



			// slice function, like matlab B = A(index,:)
			MatrixXd slice(MatrixXd A, MatrixXd index);
			bool IsIn(MatrixXd self, double pointer, int& index);

			NGP_mesh source_mesh_format, target_mesh_format;
			MatrixXd source_vertex, target_vertex,control_vertex;
			// index vector
			MatrixXd index_control;

			// geodesic distance variable
			MatrixXd geodesic_table, r_distance;
			MatrixXd target_normal;
			float coeff_rigid_tmp, coeff_smooth_tmp;
			//int num_unknown;
			// unkonwn
			/*VectorXd affine_vector;*/

			CMeshO source_mesh_vcg, target_mesh_vcg,deform_graph_vcg;
			int control_face_num,num_neighbor;


	};
}
#endif