#ifndef _NGP_H_
#define _NGP_H_

#include <iostream>
#include <fstream>
#include "wrap\io_trimesh\import_ply.h"
#include "wrap\io_trimesh\export_ply.h"
#include <Eigen/Dense>
#include <Eigen\Eigen>
#include <Eigen\Sparse>
//#include <Eigen/Sparse>
//#include <Eigen/SparseCore>
//#include <Eigen/core>

#include "quadric_tex_simp.h"
#include <algorithm>
#include <opencv2/opencv.hpp>
#include <time.h>
//#include <Eigen/Sparse>
//#include <normals.h>     //class UpdateNormals

using namespace vcg;
using namespace Eigen;


namespace NRR 
{
	// non-rigid registration 
	class NGP 
	{
		public:
		struct NGP_mesh {// own define mesh structure

			Eigen::MatrixXd vertice;
			Eigen::MatrixXd face;
			Eigen::MatrixXd uv;

		};

		public:
			NGP();
			void initial(std::string source_mesh_name, 
						 std::string target_mesh_name);
			//source or target
			void read_ply(std::string mesh_name, std::string type);
			void write_obj(Eigen::MatrixXd vertex);

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
			void ComputeGeodesic(std::string source, 
				std::vector<int> index, 
				int num, 
				Eigen::MatrixXd &geodesic_table,
				Eigen::MatrixXd &r_distance, 
				int nearest_geo);

			void NGP::Geodesic_pipeline(
				std::string source, 
				Eigen::MatrixXd index_control, 
				int num_neighbor,
				Eigen::MatrixXd& geodesic_table,
				Eigen::MatrixXd& r_distance);

			// compute weight
			Eigen::MatrixXd Weight_smooth_ajacent_geodesic(
				Eigen::MatrixXd control_vertex,
				Eigen::MatrixXd r_distance,
				Eigen::MatrixXd index_control,
				Eigen::MatrixXd geodesic_table);

			Eigen::MatrixXd weight_smooth, weight_geodesic;

			Eigen::MatrixXd WeightFunc_geodesic(Eigen::MatrixXd source_vertex,
				Eigen::MatrixXd control_vertex,
				Eigen::MatrixXd r_distance,
				Eigen::MatrixXd index_control,
				Eigen::MatrixXd geodesic_table);
			// compute normal
			Eigen::MatrixXd computeNormal(CMeshO& target_source);

			// initial affineVecotr
			Eigen::VectorXd initialAffineVector(int control_num);
			Eigen::VectorXd affineVector;

			// non-rigid icp function
			Eigen::VectorXd rsff_gv_square(
				Eigen::VectorXd affine_vector,
				Eigen::MatrixXd control_vertex,
				Eigen::MatrixXd source_vertex,
				Eigen::MatrixXd target_vertex,
				Eigen::MatrixXd target_normal,
				Eigen::MatrixXd weight,
				Eigen::MatrixXd weightTrans,
				float coeff_rigid,
				float coeff_smooth
			);

			// find index where c_matrix element bigger or smaller to value, 
			// flag==0, means "<", flag==1 means ">"
			Eigen::MatrixXd findValue(Eigen::MatrixXd c_matrix, float value, int flag);

			// index is an vector, main is N*3 matrix,
			// output the index in vector of that main matrix is 
			// equal to zero
			void erase_outlier(Eigen::MatrixXd &main, Eigen::MatrixXd index);

			// E_fit coefficient
			float alpha_point, alpha_plane;
			float threshold;// threshold for correspondense

			// apply affine vector to source_vertex to new position
			Eigen::MatrixXd ApplyTrans(
				Eigen::MatrixXd source_vertex,
				Eigen::MatrixXd a1,
				Eigen::MatrixXd a2,
				Eigen::MatrixXd a3,
				Eigen::MatrixXd b,
				Eigen::MatrixXd control_vertex,
				Eigen::MatrixXd weightTrans);

			// Gauss_newton function
			Eigen::VectorXd gauss_newton(Eigen::VectorXd affine_vector,
				Eigen::MatrixXd control_vertex,
				Eigen::MatrixXd source_vertex,
				Eigen::MatrixXd target_vertex,
				Eigen::MatrixXd target_normal,
				Eigen::MatrixXd weight_smooth,
				Eigen::MatrixXd weight_Trans,
				float coeff_rigid,
				float coeff_smooth);

			// optimization iteration time
			int iter_optimization;

			// calculate jacbian matrix
			Eigen::MatrixXd Jacobia(
				Eigen::MatrixXd a1,
				Eigen::MatrixXd a2,
				Eigen::MatrixXd a3,
				Eigen::MatrixXd control_point,
				Eigen::MatrixXd source_obj,
				Eigen::MatrixXd weight,
				Eigen::MatrixXd weightTrans,
				float coeff_rigid,
				float coeff_smooth,
				Eigen::MatrixXd symbol,
				Eigen::MatrixXd target_normal_index);

			// slice function, like matlab B = A(index,:)
			Eigen::MatrixXd slice(Eigen::MatrixXd A, Eigen::MatrixXd index);
			bool IsIn(Eigen::MatrixXd self, double pointer, int& index);

			NGP_mesh source_mesh_format, target_mesh_format;
			Eigen::MatrixXd source_vertex, target_vertex,control_vertex;
			// index vector
			Eigen::MatrixXd index_control;

			// geodesic distance variable
			Eigen::MatrixXd geodesic_table, r_distance;
			Eigen::MatrixXd target_normal;
			float coeff_rigid_tmp, coeff_smooth_tmp;
			//int num_unknown;
			// unkonwn
			/*VectorXd affine_vector;*/

			CMeshO source_mesh_vcg, target_mesh_vcg,deform_graph_vcg;
			int control_face_num,num_neighbor;


	};
}
#endif