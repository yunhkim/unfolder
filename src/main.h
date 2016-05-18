//------------------------------------------------------------------------------
//  Copyright 2007-2008 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _BF_MINKOWSKI_SUM_H_
#define _BF_MINKOWSKI_SUM_H_

//
//
//
// This is the main header file is shared by m+3d.cpp (single thread) and pm+3d.cpp
// (multi-thread). Functions defined in this file are mainly for text parsing and
// OpenGL rendering 
//
//
//
#include <list>
#include <cfloat>
using namespace std;

#include "objReader.h"
#include "model.h"
#include "unfolder.h"
#include "config.h"
#include "util/DataHelper.h"

//using namespace masc::unfolding;


//-----------------------------------------------------------------------------
// INPUTS
vector<string> filenames;
Config config;

//-----------------------------------------------------------------------------
// Intermediate data
vector<Unfolder*> unfolders;

double R = 0;       //radius
Point3d COM;     //center of mass

//-----------------------------------------------------------------------------
//for random rotation
double rot_theta;
Vector3d rot_vec;
Quaternion current_rot;

//-----------------------------------------------------------------------------
//read M+ from file
bool readfromfile();
void computeCOM_R();

//-----------------------------------------------------------------------------
// (yh) the obj data
model* initobjmodel;
model* finalobjmodel;

//-----------------------------------------------------------------------------
//
//
//
//  Open GL stuff below
//
//
//-----------------------------------------------------------------------------

#include "draw.h"

//-----------------------------------------------------------------------------
bool parseArg(int argc, char ** argv) {
  if (argc == 1)
	return false;

  for (int i = 1; i < argc; i++) {
	auto const arg = string(argv[i]);
	if (arg == "-r") {
	  config.max_retries = std::stoi(argv[++i]);
	} else if (arg == "-h") {
	  auto const heuristic = string(argv[++i]);
	  if (heuristic[0] == 's')
		config.heuristic = CutHeuristic::STEEPEST_EDGE;
	  else if (heuristic[0] == 'f')
		config.heuristic = CutHeuristic::FLAT_TREE;
	  else if (heuristic[0] == 'p')
		config.heuristic = CutHeuristic::MINIMUM_PERIMETER;
	  else if (heuristic[0] == 'r')
		config.heuristic = CutHeuristic::RANDOM;
	  else if (heuristic[0] == 'b')
		config.heuristic = CutHeuristic::BRUTE_FORCE;
	  else if (heuristic[0] == '1')
		  config.heuristic = CutHeuristic::MY_UNFOLDER_01;
	  else if (heuristic[0] == '2')
		  config.heuristic = CutHeuristic::MY_UNFOLDER_02;
	  else {
		cerr << "!Error! Unknown heuristic type = " << heuristic << endl;
		return false;
	  }
	} else if (arg == "-g") {
	  config.disable_gui = true;
	} else if (arg == "-rb") {
	  config.random_baseface = true;
	} else if (arg == "-bf" || arg == "--base-face") {
	  config.baseface = std::stoi(argv[++i]);
	} else if (arg == "-lc") {
	  config.less_cuts = true;
	} else if (arg == "-s") {
	  config.seed = std::stoi(argv[++i]);
	} else if (arg == "-q") {
	  config.quite = true;
	} else if (arg == "-ns") {
	  config.shrink = false;
	} else if (arg == "-nfb") {
	  config.find_boundary = false;
	} else if (arg == "-k") {
	  config.k = std::stoi(argv[++i]);
	} else if (arg == "-i") {
	  config.max_iterations = std::stoi(argv[++i]);
	} else if (arg == "-run") {
	  config.run = std::stoi(argv[++i]);
	} else if (arg == "-weights") {
	  config.weights_filename = argv[++i];
	} else if (arg == "-sf") {
	  config.shrink_factor = stof(argv[++i]);
	} else if (arg == "-scale") {
	  config.scale = stof(argv[++i]);
	} else if (arg == "-nl") {
	  config.dump_labels = false;
	}
	else if (arg == "-lfs") {
		config.label_font_scale = stof(argv[++i]);
	}
	else if (arg == "-tab") {
		config.add_tabs = true;
	}
	else if (arg == "-gs"){
		//yh show the obj file
		config.onlydraw = true;
		config.save_obj = true;

		//config.inner_tri_scale = 0.0;
	}
	else if (arg[0] == '-') {
	  cerr << "!Error! Unknown arg = " << arg << endl;
	  return false;
	} else {
	  filenames.push_back(string(argv[i]));
	}
  }

  return true;
}

void printUsage(char * name) {

  Config default_config;

  //int offset=20;
  cerr << "Usage: " << name << " [options] *.obj/*.off \n\n";

  cerr << "Heuristic Methods\n";
  cerr << "  -h heuristic | use heuristic method\n";
  cerr << "      s        | STEEPEST_EDGE\n";
  cerr << "      f        | FLAT_TREE (default)\n";
  cerr << "      p        | MINIMUM_PERIMETER\n";
  cerr << "      r        | RANDOM\n";
  cerr << "      1        | Your Unfolder 01\n";
  cerr << "      2        | Your Unfolder 02\n";
  cerr << "\n";

  cerr << "Unfolding\n";
  cerr << "  -s seed      | specify random seed\n";

  cerr << "  -r times     | retry times, default is "
	  << default_config.max_retries << "\n";

  cerr << "  -weights fn  | using the specify weights to unfold the mesh.\n";

  cerr << "  -q           | quite mode.\n";

  cerr << "  -bf          | specify the base face.\n";

  cerr << "  -rb          | random base face.\n";

  cerr << "  -lc          | less cuts / don't cut flat edges.\n";

  cerr << "  -g           | disable GUI. dump outputs.\n";

  cerr << "  -tab         | add tabs in the net.\n";

  cerr << "  -gs          | (yh) show the obj file.\n";

  cerr << "Dumping SVG\n";
  cerr
	  << "  -scale       | scale factor applied. both *.svg and *.ori will be affected.\n";
  cerr << "  -nl          | do not dump labels in SVG file\n";
  cerr << "  -lfs         | label font size scale [default=1.0]\n";

  cerr << endl;
  cerr << "-- Complied on " << __DATE__ << endl;
  cerr << "-- Report bugs to: Jyh-Ming Lien jmlien@cs.gmu.edu" << endl;
}

//-----------------------------------------------------------------------------

bool readfromfiles() {
  long vsize = 0;
  long fsize = 0;

  int flattened = 0;

  uint id = 0;

  auto start = clock();

  for (const auto& filename : filenames) {

	// use the same seed for all the component
	mathtool::srand48(config.seed);
	srand(config.seed);

	cout << "- [" << ++id << "/" << filenames.size() << "] Start reading "
		<< filename << endl;

	config.filename = filename;

	model* m = new model();
	if (!m->build(config.filename)) {
	  delete m;
	  return false;
	}
	cout << "- Done reading " << m->v_size << " vertices and " << m->t_size << " facets" << endl;

	Unfolder* unfolder = new Unfolder(m, config);
	unfolder->measureModel();

	// using basic heuristic
	unfolder->buildUnfolding();

	unfolder->unfoldTo(0.0);

	if (unfolder->isFlattened())
	  flattened++;

	cout << string(40, '-') << endl;

	if (unfolder->getCheckOverlappingCalls() > 0) {
	  cout << "- Total CO calls  = " << unfolder->getCheckOverlappingCalls()
		  << endl;
	  cout << "- Total CO time   = "
		  << unfolder->getTotalCheckOverlappingTime() * 1.0 / CLOCKS_PER_SEC
		  << " s" << endl;
	  cout << "- Average CO time = "
		  << unfolder->getTotalCheckOverlappingTime() * 1.0
			  / unfolder->getCheckOverlappingCalls() / CLOCKS_PER_SEC << " s"
		  << endl;
	}

	unfolders.push_back(unfolder);

	cout << string(40, '-') << endl;
  }

  auto time_cost = (clock() - start) * 1.0 / CLOCKS_PER_SEC;

  cout << "Total flattened = " << flattened << "/" << filenames.size() << endl;
  cout << "Total time = " << time_cost << " s" << endl;

  computeCOM_R();

  return true;
}

void computeCOM_R() {
//compute a bbox
  double box[6] = { FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX };
//-------------------------------------------------------------------------
  for (auto& u : unfolders) {
	auto& m = *u->getModel();
	for (int j = 0; j < m.v_size; j++) {
	  Point3d& p = m.vertices[j].p;
	  if (p[0] < box[0])
		box[0] = p[0];
	  if (p[0] > box[1])
		box[1] = p[0];
	  if (p[1] < box[2])
		box[2] = p[1];
	  if (p[1] > box[3])
		box[3] = p[1];
	  if (p[2] < box[4])
		box[4] = p[2];
	  if (p[2] > box[5])
		box[5] = p[2];
	}				//j
  }				//i

//-------------------------------------------------------------------------
// compute center of mass and R...
  COM.set((box[1] + box[0]) / 2, (box[3] + box[2]) / 2, (box[5] + box[4]) / 2);

//-------------------------------------------------------------------------
  R = 0;
  for (auto& u : unfolders) {
	auto& m = *u->getModel();
	for (int j = 0; j < m.v_size; j++) {
	  Point3d& p = m.vertices[j].p;
	  double d = (p - COM).normsqr();
	  if (d > R)
		R = d;
	}    //j
  }    //i

  R = sqrt(R);
}

void computeCOM_R2(model& m) {
	//printf("123");
	//compute a bbox
	double box[6] = { FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX };
	auto M = m;

	//printf("1");
	//-------------------------------------------------------------------------
	//for (auto& u : unfolders) {
	//auto& m = *u->getModel();
	{
		//int i = 0;
		//printf("asdf\n");

		//printf("aaaaaaaaaaaaa: %d\n", M.v_size);

		for (int j = 0; j < M.v_size; j++) {
			//printf("ASDFASDF");
			//vertex& vp = M.vertices[j];
			Point3d& p = M.vertices[j].p;
			if (p[0] < box[0])
				box[0] = p[0];
			if (p[0] > box[1])
				box[1] = p[0];
			if (p[1] < box[2])
				box[2] = p[1];
			if (p[1] > box[3])
				box[3] = p[1];
			if (p[2] < box[4])
				box[4] = p[2];
			if (p[2] > box[5])
				box[5] = p[2];
			//i++;
			//cout<<i;
		}				//j

		//printf("end\n\n");
	}
	//}				//i

	//printf("A");
	//-------------------------------------------------------------------------
	// compute center of mass and R...
	COM.set((box[1] + box[0]) / 2, (box[3] + box[2]) / 2, (box[5] + box[4]) / 2);

	//printf("B");

	//-------------------------------------------------------------------------
	R = 0;

	{
		//for (auto& u : unfolders) {
		//auto& m = *u->getModel();
		//auto& m = *finalobjmodel;//->getModel();
		for (int j = 0; j < M.v_size; j++) {
			Point3d& p = M.vertices[j].p;
			double d = (p - COM).normsqr();
			if (d > R)
				R = d;
		}    //j
		//}    //i
	}
	//printf("C");

	R = sqrt(R);

	//printf("321");
}

// yh

void makeVolume(model& u){
	auto M = &u;
	model* final_M = new model();

	// Euler's formular(F+V-E=2)
	final_M->vertices = new vertex[(M->v_size) * 6];   //
	final_M->tris = new triangle[(M->t_size) * 26];    //
	final_M->edges = new edge[(M->v_size) * 6 + (M->t_size) * 26 - 2];
	
	//add volume
	//set values; padding thickness, object thickness, and inner triangle size(I will use ratio)
	config.inner_tri_scale = 0.0;
	config.pad_thickness = 0.0;
	config.obj_thickness = 0.4;


	if (config.inner_tri_scale >= 0 && config.inner_tri_scale <= 1){
		cout << "set the inner triangle ratio" << endl;
	}
	else{
		cout << "try it agian, the inner triangle ratio" << endl;
	}


	if (config.pad_thickness >= 0){
		cout << "set the padding thickness" << endl;
	}
	else{
		cout << "try it agian, the padding thickness" << endl;
	}


	if (config.obj_thickness >= 0){
		cout << "set the object thickness" << endl;
	}
	else{
		cout << "try it agian, the object thickness" << endl;
	}

	auto obj_height = config.obj_thickness / 2.0;
	auto in_tri_height = config.obj_thickness / 2.0 + config.pad_thickness;
	auto in_tri_ratio = config.inner_tri_scale;

	for (auto i = 0; i < M->t_size; i++){
		printf("2");

		triangle& t = M->tris[i];
		Vector3d t_n_vec = t.n.get();

		//Point<double, 3>* tri_v = new Point<double, 3>[3];
		Point<double, 3>* tri_v = new Point<double, 3>[3];

		tri_v[0] = M->vertices[t.v[0]].p.get();
		tri_v[1] = M->vertices[t.v[1]].p.get();
		tri_v[2] = M->vertices[t.v[2]].p.get();

		printf("3");

		//Point<double, 3 >* b = new Point<double, 3>();
		//Point<double, 3 > b = Point<double, 3>();


		//printf("\n%fl\n", tri_v[0][1]);
		//b = ;
			//; //+ tri_v[1] + tri_v[2]) / 3.0;

		//printf("3");
		if (config.obj_thickness>0){
			printf("4");

			int temp_vid = 0;
			//ot = obj thickness
			Point<double, 3>* ot_vertex = new Point<double, 3>[6];

			for (auto j = 0; j < 6; j++) {
				printf("5");

				if (j < 3){
					ot_vertex[j] = tri_v[j] + t_n_vec*obj_height;
						//final_M->v_size++;
				}
				else{
					ot_vertex[j] = tri_v[j] - t_n_vec*obj_height;
					
				}

				printf("6");

				//printf("\n  ot_vertex:  %fl %fl %fl \n", ot_vertex[j]);
				//printf("\n  ot_vertex2:  %fl %fl %fl \n", ot_vertex[j][0], ot_vertex[j][1], ot_vertex[j][2]);
				//printf("v_size: %d\n", final_M->v_size);

				//Point3d* temp_v = new Point3d();
				//temp_v[0] = &ot_vertex[j][0];
				//temp_v[1] = &ot_vertex[j][1];
				//temp_v[2] = &ot_vertex[j][2];
				//printf("\n  ot_vertex3:  %fl %fl %fl \n", temp_v->get()[0], temp_v->get()[1], temp_v->get()[2]);
				//printf("\n  ot_vertex4:  %fl %fl %fl \n", temp_v[0], temp_v[1], temp_v[2]);

				//vertex* temp_v = new vertex();
				//temp_v->p[0] = ot_vertex[j][0];
				//temp_v->p[1] = ot_vertex[j][1];
				//temp_v->p[2] = ot_vertex[j][2];
				//printf("\n  ot_vertex3:  %fl %fl %fl \n", temp_v->p[0], temp_v->p[1], temp_v->p[2]);
				//printf("\n  ot_vertex4:  %fl %fl %fl \n", temp_v->p.get()[0], temp_v->p.get()[1], temp_v->p.get()[2]);//[0], temp_v[1], temp_v[2]);

				//Point3d temp_v = Point3d();
				//temp_v.set(ot_vertex[j][0],ot_vertex[j][1],ot_vertex[j][2]);
				//printf("\n  ot_vertex3:  %fl %fl %fl \n", temp_v.get()[0], temp_v.get()[1], temp_v.get()[2]);
				
				//if (!final_M->vertices)
				final_M->vertices[final_M->v_size].p.set(ot_vertex[j][0], ot_vertex[j][1], ot_vertex[j][2]);
				final_M->vertices[final_M->v_size].bk_p = final_M->vertices[final_M->v_size].p;

				printf("vertex: ");
				cout <<final_M->vertices[final_M->v_size].p << endl;
				
				printf("7");
				//final_M->tris[final_M->t_size].v[temp_vid] = final_M->v_size;
				final_M->v_size++;

				printf("8");
				if (j == 2 || j == 5){
					temp_vid = 0;
					final_M->t_size++;
				}
			}

		}

		//	model* M = new model();
		//	
		//	M->v_size = u.v_size;
		//	M->t_size = u.t_size;
		//
		//	
		//	//vertices = new vertex[v_size];   //
		//	tris = new triangle[t_size];     //
		//	assert(vertices && tris);        //make sure enough memory
		//
		//	//copy vertices
		//	for (uint i = 0; i < v_size; i++) {
		//		vertices[i].p.set(&data.pts[i].x);
		//		vertices[i].bk_p = vertices[i].p;  //backup for modification
		//	}
		//
		//	cout << "vertices copied!" << endl;
		//
		//	//copy triangles
		//	int tid = 0;
		//	for (auto i = data.polys.begin(); i != data.polys.end(); i++) {
		//		const auto& ids = i->pts;
		//		//check if triangle
		//		if (ids.size() != 3) {
		//			cerr << "! Error: polygon " << tid << " is not a triangle, edge size is "
		//				<< ids.size() << "." << endl;
		//			return false;
		//		}
		//		int vid = 0;
		//		for (auto j = ids.begin(); j != ids.end(); j++) {
		//			assert(*j >= 0 && *j < v_size);
		//			tris[tid].v[vid++] = *j;
		//			vertices[*j].m_f.push_back(tid);
		//		}
		//		tid++;
		//	}
		//
		//	cout << "triangle copied!" << endl;

		//*finalobjmodel = M->build2(*initobjmodel, );
	//	M->tris;
	//	M->edges;
	//	M->vertices;
		//M->
		//printf("%fl", v1[0]);
		//printf("%fl", M->vertices[t.v[0]].p.get()[0]);
		

	}


	
	//finalobjmodel = M;
	finalobjmodel = final_M;


}
bool readfromfiles2() {
	long vsize = 0;
	long fsize = 0;

	int flattened = 0;

	uint id = 0;

	//auto start = clock();

	for (const auto& filename : filenames) {

		// use the same seed for all the component
		//mathtool::srand48(config.seed);
		//srand(config.seed);

		cout << "- [" << ++id << "/" << filenames.size() << "] Start reading "
			<< filename << endl;

		config.filename = filename;

		model* m = new model();
		if (!m->build(config.filename)) {
			delete m;
			return false;
		}

		initobjmodel = m;
		cout << "- Done reading " << m->v_size << " vertices and " << m->t_size << " facets" << endl;

		/*
		Unfolder* unfolder = new Unfolder(m, config);
		unfolder->measureModel();

		// using basic heuristic
		unfolder->buildUnfolding();

		unfolder->unfoldTo(0.0);

		if (unfolder->isFlattened())
			flattened++;

		cout << string(40, '-') << endl;

		if (unfolder->getCheckOverlappingCalls() > 0) {
			cout << "- Total CO calls  = " << unfolder->getCheckOverlappingCalls()
				<< endl;
			cout << "- Total CO time   = "
				<< unfolder->getTotalCheckOverlappingTime() * 1.0 / CLOCKS_PER_SEC
				<< " s" << endl;
			cout << "- Average CO time = "
				<< unfolder->getTotalCheckOverlappingTime() * 1.0
				/ unfolder->getCheckOverlappingCalls() / CLOCKS_PER_SEC << " s"
				<< endl;
		}

		unfolders.push_back(unfolder);
		*/
		cout << string(40, '-') << endl;
	}

	//*finalobjmodel = makeVolume(*initobjmodel);
	//computeCOM_R2(*finalobjmodel);

	makeVolume(*initobjmodel);
	computeCOM_R2(*finalobjmodel);
	
	return true;
}



#endif //_BF_MINKOWSKI_SUM_H_

