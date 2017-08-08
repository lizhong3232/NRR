#include "NGP.h"

namespace NRR
{
	enum MeshElement {
		MM_NONE = 0x00000000,
		MM_VERTCOORD = 0x00000001,
		MM_VERTNORMAL = 0x00000002,
		MM_VERTFLAG = 0x00000004,
		MM_VERTCOLOR = 0x00000008,
		MM_VERTQUALITY = 0x00000010,
		MM_VERTMARK = 0x00000020,
		MM_VERTFACETOPO = 0x00000040,
		MM_VERTCURV = 0x00000080,
		MM_VERTCURVDIR = 0x00000100,
		MM_VERTRADIUS = 0x00000200,
		MM_VERTTEXCOORD = 0x00000400,
		MM_VERTNUMBER = 0x00000800,

		MM_FACEVERT = 0x00001000,
		MM_FACENORMAL = 0x00002000,
		MM_FACEFLAG = 0x00004000,
		MM_FACECOLOR = 0x00008000,
		MM_FACEQUALITY = 0x00010000,
		MM_FACEMARK = 0x00020000,
		MM_FACEFACETOPO = 0x00040000,
		MM_FACENUMBER = 0x00080000,

		MM_WEDGTEXCOORD = 0x00100000,
		MM_WEDGNORMAL = 0x00200000,
		MM_WEDGCOLOR = 0x00400000,

		// 	Selection
		MM_VERTFLAGSELECT = 0x00800000,
		MM_FACEFLAGSELECT = 0x01000000,

		// Per Mesh Stuff....
		MM_CAMERA = 0x08000000,
		MM_TRANSFMATRIX = 0x10000000,
		MM_COLOR = 0x20000000,
		MM_POLYGONAL = 0x40000000,
		MM_UNKNOWN = 0x80000000,

		MM_ALL = 0xffffffff
	};


	void NGP::updateDataMask(CMeshO &cm, int neededDataMask)
	{
		if ((neededDataMask & MM_FACEFACETOPO) != 0)
		{
			cm.face.EnableFFAdjacency();
			tri::UpdateTopology<CMeshO>::FaceFace(cm);
		}
		if ((neededDataMask & MM_VERTFACETOPO) != 0)
		{

			cm.vert.EnableVFAdjacency();
			cm.face.EnableVFAdjacency();

			tri::UpdateTopology<CMeshO>::VertexFace(cm);
		}

		if (((neededDataMask & MM_WEDGTEXCOORD) != 0))	cm.face.EnableWedgeTexCoord();
		if (((neededDataMask & MM_FACECOLOR) != 0))	cm.face.EnableColor();
		if (((neededDataMask & MM_FACEQUALITY) != 0))	cm.face.EnableQuality();
		if (((neededDataMask & MM_FACEMARK) != 0))		cm.face.EnableMark();
		if (((neededDataMask & MM_VERTMARK) != 0))		cm.vert.EnableMark();
		if (((neededDataMask & MM_VERTCURV) != 0))		cm.vert.EnableCurvature();
		if (((neededDataMask & MM_VERTCURVDIR) != 0))	cm.vert.EnableCurvatureDir();
		if (((neededDataMask & MM_VERTRADIUS) != 0))	cm.vert.EnableRadius();
		if (((neededDataMask & MM_VERTTEXCOORD) != 0))	cm.vert.EnableTexCoord();
	}


	bool NGP::Simplify(CMeshO &cm, int targetfacenum, float QualityThr, float Extratcoordw, bool OptimalPlacement,
		bool PreserveBoundary, float BoundaryWeight, bool PlanarQuadric, bool PreserveNormal) {


		updateDataMask(cm, MeshElement::MM_VERTFACETOPO | MeshElement::MM_VERTMARK);
		tri::UpdateFlags<CMeshO>::FaceBorderFromVF(cm);

		if (!tri::HasPerWedgeTexCoord(cm))
		{
			cerr << "Warning: nothing have been done. Mesh has no Texture." << flush;
			//return false;
		}
		//if (!tri::Clean<CMeshO>::HasConsistentPerWedgeTexCoord(cm)) {
		//	cerr << "Mesh has some inconsistent tex coords (some faces without texture)" << flush;;
		//	//return false; 
		//}

		int TargetFaceNum = targetfacenum;

		vcg::tri::TriEdgeCollapseQuadricTexParameter pp;

		pp.QualityThr = QualityThr;
		pp.ExtraTCoordWeight = Extratcoordw;
		pp.OptimalPlacement = OptimalPlacement;
		pp.PreserveBoundary = PreserveBoundary;
		pp.BoundaryWeight = pp.BoundaryWeight * BoundaryWeight;
		pp.QualityQuadric = PlanarQuadric;
		pp.NormalCheck = PreserveNormal;

		QuadricTexSimplification(cm, TargetFaceNum, false, pp);

		tri::UpdateBounding<CMeshO>::Box(cm);
		if (cm.fn > 0) {
			tri::UpdateNormal<CMeshO>::PerFaceNormalized(cm);
			tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(cm);
		}

		return true;
	}




	void NGP::QuadricTexSimplification(CMeshO &m, int  TargetFaceNum, bool Selected, tri::TriEdgeCollapseQuadricTexParameter &pp)
	{
		tri::UpdateNormal<CMeshO>::PerFace(m);
		math::Quadric<double> QZero;
		QZero.SetZero();
		QuadricTexHelper<CMeshO>::QuadricTemp TD3(m.vert, QZero);
		QuadricTexHelper<CMeshO>::TDp3() = &TD3;

		std::vector<std::pair<vcg::TexCoord2<float>, Quadric5<double> > > qv;

		QuadricTexHelper<CMeshO>::Quadric5Temp TD(m.vert, qv);
		QuadricTexHelper<CMeshO>::TDp() = &TD;

		if (Selected) // simplify only inside selected faces
		{
			// select only the vertices having ALL incident faces selected
			tri::UpdateSelection<CMeshO>::VertexFromFaceStrict(m);

			// Mark not writable un-selected vertices
			CMeshO::VertexIterator  vi;
			for (vi = m.vert.begin(); vi != m.vert.end(); ++vi) if (!(*vi).IsD())
			{
				if (!(*vi).IsS()) (*vi).ClearW();
				else (*vi).SetW();
			}
		}

		vcg::LocalOptimization<CMeshO> DeciSession(m, &pp);

		DeciSession.Init<MyTriEdgeCollapseQTex>();

		if (Selected)
			TargetFaceNum = m.fn - (m.sfn - TargetFaceNum);
		DeciSession.SetTargetSimplices(TargetFaceNum);
		DeciSession.SetTimeBudget(0.1f);
		//	int startFn=m.fn;

		int faceToDel = m.fn - TargetFaceNum;

		while (DeciSession.DoOptimization() && m.fn > TargetFaceNum)
		{
			char buf[256];
			sprintf(buf, "Simplifing heap size %i ops %i\n", int(DeciSession.h.size()), DeciSession.nPerfmormedOps);
		};

		DeciSession.Finalize<MyTriEdgeCollapseQTex>();

		if (Selected) // Clear Writable flags 
		{
			CMeshO::VertexIterator  vi;
			for (vi = m.vert.begin(); vi != m.vert.end(); ++vi)
				if (!(*vi).IsD()) (*vi).SetW();
		}
	}
}
