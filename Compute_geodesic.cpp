#include "NGP.h"
#include "geodesic.h"
#include <iomanip>
#include <numeric>    

namespace NRR
{	
	class MyEdge;
	class MyFace;
	class MyVertex;
	struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
		Use<MyEdge>     ::AsEdgeType,
		Use<MyFace>     ::AsFaceType> {};

	class MyVertex : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::Mark, vertex::VFAdj, vertex::Color4b, vertex::Qualityf, vertex::BitFlags  > {};
	class MyFace : public Face< MyUsedTypes, face::VFAdj, face::VertexRef, face::Normal3f, face::BitFlags > {};
	class MyEdge : public Edge<MyUsedTypes> {};
	class MyMesh : public tri::TriMesh< vector<MyVertex>, vector<MyFace>, vector<MyEdge>  > {};


	template <typename T>
	vector<size_t> sort_indexes(const vector<T> &v) {

		// initialize original index locations
		vector<size_t> idx(v.size());
		std::iota(idx.begin(), idx.end(), 0);

		// sort indexes based on comparing values in v
		sort(idx.begin(), idx.end(),
			[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

		return idx;
	}

	void NGP::Geodesic_pipeline(std::string source, MatrixXd index_control, int num_neighbor, 
								MatrixXd& geodesic_table, MatrixXd& r_distance)
	{
		std::vector<int> index;
		for (int i = 0; i < index_control.rows(); ++i) {
			index.push_back(index_control(i,0));
		}

		MatrixXd geodesic_table_tmp(source_vertex.rows(), index.size());
		MatrixXd r_distance_tmp(source_vertex.rows(), 1);// output

		/********* rest of geodesic pipline*******/
		auto step_size = 100ul;
		auto total_steps = source_vertex.rows() / step_size + 1;

		size_t steps_completed = 0;
		long double sum = 0;
		std::cout << "Compute geodesic distance for each node..." << std::endl;
#pragma omp parallel 
		{
			size_t local_count = 0;


			//#pragma omp for reduction(+:sum)
#pragma omp for
			for (int i = 0; i < index.size(); ++i)
			{
				ComputeGeodesic(source, index, i, geodesic_table_tmp, r_distance_tmp, num_neighbor);

				if (local_count++ % step_size == step_size - 1)
				{
#pragma omp atomic
					++steps_completed;

					if (steps_completed % 5 == 1)
					{
						//#pragma omp critical
						std::cout << "Progress: " << steps_completed << " of " << total_steps << " (" << std::fixed << std::setprecision(1) << (100.0*steps_completed / total_steps) << "%)\n";
					}
				}
			}
		}
		std::vector<size_t> index_sort;
		std::vector<double> distance;
		for (int i = 0; i < geodesic_table_tmp.rows(); ++i) {
			distance.clear();
			index_sort.clear();
			for (int j = 0; j < geodesic_table_tmp.cols(); ++j) {
				distance.push_back(geodesic_table_tmp(i, j));
			}
			index_sort = sort_indexes(distance);
			r_distance_tmp(i, 0) = distance[index_sort[num_neighbor]];
		}

		geodesic_table = geodesic_table_tmp;
		r_distance     = r_distance_tmp;

		std::cout << "..done" << std::endl;

		/*std::cout << geodesic_table.rows() << std::endl;*/

		/* debug   */
		//char filename[50];
		//sprintf_s(filename, "r_distance%.3d.txt", 1);
		//std::ofstream outfile(filename);
		///*for (int i = 0; i < r_distance.rows(); i++) {
		//outfile << r_distance(i, 0) << std::endl;
		//}*/
		//if (outfile.is_open())
		//{
		//	outfile << r_distance;
		//	//file << "m" << '\n' << colm(m) << '\n';
		//}
		//outfile.close();
		//std::cout << "R_distance complete!" << std::endl;

		//char geodesic_filename[100];
		//sprintf_s(geodesic_filename, "geodesic_table%.3d.txt", 1);
		//std::ofstream file(geodesic_filename);
		//if (file.is_open())
		//{
		//	file << geodesic_table;
		//	//file << "m" << '\n' << colm(m) << '\n';
		//}
		//file.close();
	}


	void NGP::ComputeGeodesic(std::string source, std::vector<int> index, int num, MatrixXd &geodesic_table,
		MatrixXd &r_distance, int nearest_geo) {

		MyMesh m;

		// initial matrix
		//MatrixXd geodesic_table_tmp(m.vert,)

		//  if(tri::io::ImporterPLY<MyMesh>::Open(m,"../../meshes/disk_irregular_1k.ply")!=0)
		if (tri::io::ImporterPLY<MyMesh>::Open(m, source.c_str()) != 0)
		{
			printf("Error reading file\n");
			exit(0);
		}
		//Point3f c = m.bbox.Center();
		Point3f c = m.vert[index[num]].P();

		MyVertex*closest = &*m.vert.begin();
		float minDist = Distance(closest->P(), c);
		for (MyMesh::VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
		{
			if (Distance(vi->P(), c)<minDist)
			{
				minDist = Distance(vi->P(), c);
				closest = &*vi;
			}
		}
		vector<MyVertex*> seedVec;
		seedVec.push_back(closest);
		tri::EuclideanDistance<MyMesh> ed;
		//tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
		tri::Allocator<MyMesh>::CompactEveryVector(m);
		tri::UpdateTopology<MyMesh>::VertexFace(m);

		//tri::Geodesic::DistanceFunctor distFunc;
		//int t0 = clock();
		tri::Geodesic<MyMesh>::Compute(m, seedVec, ed);
		//pair<float, float> minmax = tri::Stat<MyMesh>::ComputePerVertexQualityMinMax(m);
		//int t1 = clock();

		//printf("Geodesic distance %6.3f\n", float(t1 - t0) / CLOCKS_PER_SEC);
		//tri::UpdateColor<MyMesh>::PerVertexQualityRamp(m);

		//std::cout << m.vert[1].Q() << std::endl;
		//std::cout << m.vert[1].Q() << std::endl;
		//std::cout << m.vert[2].Q() << std::endl;
		// load into vector
		std::vector<double> distance;
		std::vector<size_t> index_sort;

		for (int i = 0; i < m.vert.size(); ++i) {
			/*distance.push_back(m.vert[i].Q());*/
			geodesic_table(i, num) = m.vert[i].Q();
		}
		
	}

}