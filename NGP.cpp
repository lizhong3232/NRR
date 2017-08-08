#include "NGP.h"

namespace NRR
{
	NGP::NGP()
	{
		control_face_num = 50;
		num_neighbor = 4;

		coeff_rigid_tmp  = 100.0;
		coeff_smooth_tmp = 50.0;


		alpha_point = 0.1;
		alpha_plane = 1.0;
		threshold = 0.07;
		//num_unknown = 12;
	}

	// read mesh , extract deform graph
	void NGP::initial(std::string source_mesh_name, std::string target_mesh_name)
	{
		read_ply(source_mesh_name, "source");
		read_ply(target_mesh_name, "target");

		std::cout <<"source mesh has " << source_vertex.rows() <<" vertices" << std::endl;
		std::cout <<"target mesh has " << target_vertex.rows() <<" vertices" << std::endl;
		
		// extract deform graph
		extract_deform_nodes(source_mesh_vcg, control_face_num);

		// construct weight using geodesic distance
		Geodesic_pipeline(source_mesh_name, index_control, num_neighbor,
						  geodesic_table,r_distance);

		// initial weights
		 weight_smooth = Weight_smooth_ajacent_geodesic(
												control_vertex,
												r_distance,
												index_control,
												geodesic_table);
		// weight geodesic
		weight_geodesic = WeightFunc_geodesic(source_vertex,
														control_vertex,
														r_distance,
														index_control,
														geodesic_table);
		// compute normal
		target_normal = computeNormal(target_mesh_vcg);
		
		// initial unkonwn 12 * node size
		affineVector = initialAffineVector(control_vertex.rows());

	}

	/* extract deform nodes in matrixXd format */
	void NGP::extract_deform_nodes(CMeshO& test_mesh_vcg, int control_face_num)
	{
		// extract deform graph
		vcg::tri::Append<CMeshO, CMeshO>::MeshCopy(deform_graph_vcg, test_mesh_vcg);
		deform_graph_vcg.face.EnableWedgeTexCoord();
		Simplify(deform_graph_vcg, control_face_num, 0.3, 1.0, true, false, 1.0, false, false);

		// extract unreference part and save deform graph
		vcg::tri::Clean<CMeshO>::RemoveUnreferencedVertex(deform_graph_vcg);
		// update geometry
		vcg::tri::Allocator<CMeshO>::CompactEveryVector(deform_graph_vcg);
		tri::UpdateBounding<CMeshO>::Box(deform_graph_vcg);
		// debug
		std::cout << "deform graph has " << deform_graph_vcg.VN() << " nodes." << std::endl;
		int mask = tri::io::ExporterPLY<CMeshO>::GetExportMaskCapability();
		tri::io::ExporterPLY<CMeshO>::Save(deform_graph_vcg, "deform_graph.ply", mask, true);
		// convert to matrixXd version
		// load to structure
		MatrixXd vertex_tmp(deform_graph_vcg.VN(), 3);
		//MatrixXd face_tmp(mesh.FN(), 3);
		for (int i = 0; i < deform_graph_vcg.VN(); ++i) {
			vertex_tmp(i, 0) = deform_graph_vcg.vert[i].P()[0];
			vertex_tmp(i, 1) = deform_graph_vcg.vert[i].P()[1];
			vertex_tmp(i, 2) = deform_graph_vcg.vert[i].P()[2];
			/*std::cout << deform_graph_vcg.vert[i].P()[0]
				<< deform_graph_vcg.vert[i].P()[1]
				<< deform_graph_vcg.vert[i].P()[2] << std::endl;*/
		}
		//write_obj(vertex_tmp);
		//control_vertex = vertex_tmp;
		//std::cout << "deform graph has " << vertex_tmp.rows() << " nodes." << std::endl;
		
		// knnsearch the nearest deform-graph from deform_graph to source mesh
		bool flag = KNearestSearch(vertex_tmp, source_vertex,1, index_control);
		std::cout << index_control.rows() << std::endl;
		//deform_graph_vcg.vert.EnableMark();
		for (int i = 0; i < index_control.rows();++i)
		{
			int current_index = index_control(i, 0);
			/*std::cout << index_control(i, 0) << std::endl;
			std::cout << source_vertex.row(current_index)  << std::endl;
			std::cout << vertex_tmp.row(i) << std::endl<<std::endl;*/
			// rewrite deform_graph_vcg, bugg?
			deform_graph_vcg.vert[i].P()[0] = source_vertex(current_index, 0);
			deform_graph_vcg.vert[i].P()[1] = source_vertex(current_index, 1);
			deform_graph_vcg.vert[i].P()[2] = source_vertex(current_index, 2);
			// test eigen matrix
			vertex_tmp(i, 0) = source_vertex(current_index, 0);
			vertex_tmp(i, 1) = source_vertex(current_index, 1);
			vertex_tmp(i, 2) = source_vertex(current_index, 2);
		}

	   // give value to control_vertex
		control_vertex = vertex_tmp;
		// debug
		//tri::io::ExporterPLY<CMeshO>::Save(deform_graph_vcg, "deform_graph_revise.ply", mask, false);


	}

	void NGP::write_obj(MatrixXd vertex) {
		std::ofstream outfile("debug_graph.obj");
		for (int i = 0; i < vertex.rows(); i++) {
			outfile <<"v "<< vertex(i, 0) << " " << vertex(i, 1) 
				<< " " << vertex(i, 2) << std::endl;
		}
		outfile.close();
	}

	//read ply
	void NGP::read_ply(std::string mesh_name, std::string type)
	{
		int mask;
		CMeshO mesh;
		// import source ply
		tri::io::ImporterPLY<CMeshO>::LoadMask(mesh_name.c_str(), mask);
		int result = tri::io::ImporterPLY <CMeshO>::Open(mesh, mesh_name.c_str(), mask);
		if (result != 0) // all the importers return 0 on success
		{
			if (tri::io::ImporterPLY<CMeshO>::ErrorCritical(result))
			{
				if (tri::io::ImporterPLY<CMeshO>::ErrorCritical(result)) {
					std::cerr << "Load file error:  , " << 
					tri::io::ImporterPLY<CMeshO>::ErrorMsg(result);//<< flush;
					exit(0);
				}
			}
		}
		// load to structure
		MatrixXd vertex_tmp(mesh.VN(), 3);
		//MatrixXd face_tmp(mesh.FN(), 3);
		for (int i = 0; i < mesh.VN(); ++i) {
			vertex_tmp(i, 0) = mesh.vert[i].P()[0];
			vertex_tmp(i, 1) = mesh.vert[i].P()[1];
			vertex_tmp(i, 2) = mesh.vert[i].P()[2];

			//std::cout << vertex_tmp.row(i) << std::endl;
		}

		if (type == std::string("source")) {
			source_vertex = vertex_tmp;
			vcg::tri::Append<CMeshO, CMeshO>::MeshCopy(source_mesh_vcg, mesh);
		}
		else if (type == std::string("target")) {
			target_vertex = vertex_tmp;
			vcg::tri::Append<CMeshO, CMeshO>::MeshCopy(target_mesh_vcg, mesh);
		}

	}

	// k nearest search in opencv, is it fast?
	bool NGP::KNearestSearch(Eigen::MatrixXd source,
		Eigen::MatrixXd target,
		int k,
		Eigen::MatrixXd &n_index)
	{
		// k is k neighbor
		//std::cout << source.rows() << std::endl;

		/**** initiliaze points ***********/
		cv::Mat features(target.rows(), 3, CV_32F); // source_vertex
		//cv::Mat reference(source.rows(), 3, CV_32F);//vertex_tmp

		//#pragma omp parallel for
		for (int row = 0; row < features.rows; row++) {
			for (int col = 0; col < features.cols; col++) {
				features.at<float>(row, col) = target(row, col);
			}
		}
		// initial reference
	/*	for (int row = 0; row < reference.rows; row++) {
			for (int col = 0; col < reference.cols; col++) {
				reference.at<float>(row, col) = source(row, col);
			}
		}*/

		/***************** find K nearest neighbor ***************/
		// KdTree with 5 random trees
		// KdTree with 5 random trees
		Eigen::MatrixXd neighbors(source.rows(), k);// matrix load 10 neighbor for each point
		Eigen::MatrixXd distance(source.rows(), k);// matrix load 10 neighbor for each point


		cv::flann::LinearIndexParams indexParams;
		//create the Index
		cv::flann::Index kdtree(features, indexParams);//initilizate index

													 
		for (int i = 0; i < source.rows(); i++) {
			std::vector<float> singleQuery;
			std::vector<int> index(k);
			std::vector<float> dist(k);

			// Searching for the Mean
			singleQuery.push_back(source(i, 0));
			singleQuery.push_back(source(i, 1));
			singleQuery.push_back(source(i, 2));
			kdtree.knnSearch(singleQuery, index, dist, k, cv::flann::SearchParams(64));//get index for each query
																					   //
			for (int j = 0; j < k; j++) {
				neighbors(i, j) = index[j];
				distance(i, j) = sqrt(dist[j]);
			}			
		}
		
		n_index = neighbors;
		//r_distance = distance;

		return true;

	}


	// k nearest search with distance enable
	bool NGP::KNearestSearch_distance(Eigen::MatrixXd source,
		Eigen::MatrixXd target,
		int k,
		Eigen::MatrixXd &n_index, 
		Eigen::MatrixXd &distance_index)
	{
		// k is k neighbor
		//std::cout << source.rows() << std::endl;

		/**** initiliaze points ***********/
		cv::Mat features(target.rows(), 3, CV_32F); // source_vertex
													//cv::Mat reference(source.rows(), 3, CV_32F);//vertex_tmp

													//#pragma omp parallel for
		for (int row = 0; row < features.rows; row++) {
			for (int col = 0; col < features.cols; col++) {
				features.at<float>(row, col) = target(row, col);
			}
		}
	
		/***************** find K nearest neighbor ***************/
		// KdTree with 5 random trees
		// KdTree with 5 random trees
		Eigen::MatrixXd neighbors(source.rows(), k);// matrix load 10 neighbor for each point
		Eigen::MatrixXd distance(source.rows(), k);// matrix load 10 neighbor for each point


		cv::flann::LinearIndexParams indexParams;
		//create the Index
		cv::flann::Index kdtree(features, indexParams);//initilizate index


		for (int i = 0; i < source.rows(); i++) {
			std::vector<float> singleQuery;
			std::vector<int> index(k);
			std::vector<float> dist(k);

			// Searching for the Mean
			singleQuery.push_back(source(i, 0));
			singleQuery.push_back(source(i, 1));
			singleQuery.push_back(source(i, 2));
			kdtree.knnSearch(singleQuery, index, dist, k, cv::flann::SearchParams(64));//get index for each query
																					   //
			for (int j = 0; j < k; j++) {
				neighbors(i, j) = index[j];
				distance(i, j) = sqrt(dist[j]);
			}
		}

		n_index = neighbors;
		distance_index = distance;

		return true;

	}

	// computer normal
	MatrixXd NGP::computeNormal(CMeshO &target_source)
	{
		if (target_source.fn > 0) {
			/*tri::UpdateNormal<CMeshO>::PerFaceNormalized(target_source);
			tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(target_source);*/
			//tri::UpdateNormals<CMeshO>::PerVertexPerFace(target_source);
			tri::UpdateNormal<CMeshO>::PerVertexNormalized(target_source);
		}
		MatrixXd normal_tmp = MatrixXd::Zero(target_source.VN(), 3);

		for (int i = 0; i < target_source.VN(); ++i) 
		{
			normal_tmp(i, 0) = target_source.vert[i].N()[0];
			normal_tmp(i, 1) = target_source.vert[i].N()[1];
			normal_tmp(i, 2) = target_source.vert[i].N()[2];	
		}

		// debug
		/*std::ofstream outfile("normal.txt");
		if (outfile.is_open())
		{
		outfile << normal_tmp;
		}
		outfile.close();
*/
		return normal_tmp;
	}

	
}



