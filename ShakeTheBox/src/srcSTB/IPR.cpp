#include <NumDataIO.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <ctime>
#include <omp.h>
#include <ratio>
#include <chrono>

#include "IPR.h"
#include "Common.h"

#include <conio.h>

#include "Shlwapi.h"
#include <numeric>

using namespace std;
#define ALL_CAMS 100
#define REDUCED_CAMS 0
#define STBflag true
#define IPRflag false

//#define SAVE_INTERMEDIATE_PROJS

IPR::IPR(string& fname, int ncams, string& prjPath) : ncams(ncams)
{
	/*
	File contains
	--------------------------------
	# Path to cam calib file
	# Path to tiff files
	# Path to OTF .mat file
	# avg. particle size in pixels_orig (in 1D)
	# No. of outerloop iterations
	# No. of innerloop iterations
	# threshold for 2D finder
	# no. of bits/pixel
	# use reduce cams yes / no
	# no. of loops for each reduced camera combination
	-------------------------------
	*/
	// remove comments from the file 
	ifstream infile(fname.c_str(), ios::in);
	string line;
	stringstream parsed;
	while (getline(infile, line)) {
		size_t commentpos = line.find('#');
		if (commentpos > 0) {
			if (commentpos < string::npos) {
				line.erase(commentpos);
			}
			parsed << line << '\t';
		}
	}
	infile.close();

	int tri; triangulationOnly = false;									// Only use triangulation?
	parsed >> tri;
	if (tri != 0)
		triangulationOnly = true;

	int ipr; IPROnly = false;
	parsed >> ipr;
	if (ipr != 0)
		IPROnly = true;

	char buf[MAX_PATH];
	parsed >> calibfile;												// taking the path to calib file
	calibfile = PathCombineA(buf, prjPath.c_str(), calibfile.c_str());
	parsed >> tiffaddress;												// taking the path to tiff files	
	tiffaddress = PathCombineA(buf, prjPath.c_str(), tiffaddress.c_str());

	parsed >> otfFile;													// taking the path to .mat / .txt file with OTF parameters
	otfFile = PathCombineA(buf, prjPath.c_str(), otfFile.c_str());
	parsed >> psize;													// particle size in x or y direction

	parsed >> it_outerloop;												// # of outer and innerloop iterations
	parsed >> it_innerloop;
	cout << "# of outerloops: " << it_outerloop << endl;
	cout << "# of innerloops: " << it_innerloop << endl;

	parsed >> threshold;												// threshold for 2D finder
	int bitMax;															// max pixel value
	parsed >> bitMax;
	pixelMax = (int)pow(2, bitMax) - 1;

	parsed >> intensityLower;											// ghost particle intensity threshold
	parsed >> mindist_2D;												// reading in tolerences
	parsed >> mindist_3D;
	mindist_2D = mindist_2D * config.factor; // unit conversion
	mindist_3D = mindist_3D * config.factor;

	int reduced;														// reduced cams?
	parsed >> reduced;
	if (reduced == 1)
		reducedCams = true;

	parsed >> it_reducedCam;											// reading in tolerences for reduced cams
	parsed >> mindist_2Dr;
	parsed >> mindist_3Dr;
	mindist_2Dr = mindist_2Dr * config.factor;
	mindist_3Dr = mindist_3Dr * config.factor;

	ifstream infile2(calibfile.c_str(), ios::in);						// getting camera parameters
	string line2;
	parsed.str("");
	while (getline(infile2, line2)) {
		size_t commentpos = line2.find('#');
		if (commentpos > 0) {
			if (commentpos < string::npos) {
				line2.erase(commentpos);
			}
			parsed << line2 << '\t';
		}
	}
	infile.close();
	int dummy;
	parsed >> dummy;

	for (int i = 0; i < ncams; ++i)										// now read in the parameters for each camera
		camsAll.push_back(Camera(parsed));

	if (ncams >= 4)														// if there are more than 4 cams, reading the first 4 cams into cams4  TODO: Why?
		for (int i = 0; i < 4; ++i)
			cams4.push_back(camsAll[i]);
	else
		for (int i = 0; i < ncams; ++i)
			cams4.push_back(camsAll[i]);

	Npixh = camsAll[0].Get_Npixh();										// # of pixels in each dimension
	Npixw = camsAll[0].Get_Npixw();

	for (int n = 0; n < ncams; n++) {									// initializing original, residual and reprojected matrices
		pixels_orig.push_back(new int* [Npixh]);
		pixels_reproj.push_back(new int* [Npixh]);
		pixels_res.push_back(new int* [Npixh]);
		for (int i = 0; i < Npixh; i++) {
			pixels_orig[n][i] = new int[Npixw] {};
			pixels_reproj[n][i] = new int[Npixw] {};
			pixels_res[n][i] = new int[Npixw] {};
		}
	}
	m_particle_position_addr = tiffaddress + "ParticlePositions/";
	m_doing_STB = false;
}

void IPR::TestParticleGeneration(deque<int**>& _pixels_orig)
{
	size_t N = TEST_PARTICLES_GENERATION_N;

	if (m_bFirstTime)
	{
		m_bFirstTime = false;

		//auto get_random(std::bind(std::uniform_real_distribution<double>(0, 10), std::default_random_engine()));

		//rnd_points.resize(N);

		//for (size_t i = 0; i < N; ++i) {
		//	rnd_points[i] = Position(get_random(), get_random(), get_random());

		//}

		ifstream f("D:\\Users\\Mike\\calibpoints.txt");
		double NewX, NewY, NewZ, Newx1, Newy1, Newx2, Newy2, Newx3, Newy3, Newx4, Newy4;
		const int NN = 2688;
		Position p[NN];
		double ex1 = 0, ey1 = 0, ex2 = 0, ey2 = 0, ex3 = 0, ey3 = 0, ex4 = 0, ey4 = 0;
		for (int i = 0; i < NN; ++i)
		{
			f >> NewX;
			f >> NewY;
			f >> NewZ;
			f >> Newx1;
			f >> Newy1;
			f >> Newx2;
			f >> Newy2;
			f >> Newx3;
			f >> Newy3;
			f >> Newx4;
			f >> Newy4;
			p[i] = Position(NewX, NewY, NewZ, Newx1, Newy1, Newx2, Newy2, Newx3, Newy3, Newx4, Newy4, 0);
			//cout << p[i] << "\n";

			Position pos2D[4];
			for (int n = 0; n < ncams; n++) {
				pos2D[n] = camsAll[n].WorldToImage(p[i]); pos2D[n].Set_Z(0);
				pos2D[n] = camsAll[n].Distort(pos2D[n]);
			}

			ex1 += pos2D[0].X() - Newx1;
			ey1 += pos2D[0].Y() - Newy1;
			ex2 += pos2D[1].X() - Newx2;
			ey2 += pos2D[1].Y() - Newy2;
			ex3 += pos2D[2].X() - Newx3;
			ey3 += pos2D[2].Y() - Newy3;
			ex4 += pos2D[3].X() - Newx4;
			ey4 += pos2D[3].Y() - Newy4;
			//cout << pos2D[0].X() - Newx1 << " " << pos2D[0].Y() - Newy1 << " " << pos2D[1].X() - Newx2 << " " << pos2D[1].Y() - Newy2 << " " <<
			//	pos2D[2].X() - Newx3 << " " << pos2D[2].Y() - Newy3 << " " << pos2D[3].X() - Newx4 << " " << pos2D[3].Y() - Newy4 << "\n";
		}

		cout << ex1 / NN << " " << ey1 / NN << " " << ex2 / NN << " " << ey2 / NN << " " << ex3 / NN << " " << ey3 / NN << " " << ex4 / NN << " " << ey4 / NN << "\n";

		cin.get();
	}
	else
	{
		for (size_t i = 0; i < N; ++i) rnd_points[i] = Position(rnd_points[i].X(), rnd_points[i].Y() + 0.42, rnd_points[i].Z());
	}

	//int Counter = 0;
	//rnd_points.resize(11 * 21);
	//for (int j = -25; j <= 25; j += 5)
	//{
	//	for (int i = -50; i <= 50; i += 5)
	//	{
	//		rnd_points[Counter++] = Position(i, j, 5 + 10 * (frame - 1));
	//	}
	//}

	Frame my_estimate(rnd_points);
	ReprojImage(my_estimate, OTF(ncams, otfFile), _pixels_orig, STBflag);

#ifdef TEST_PARTICLES_GENERATION_SAVE
	NumDataIO<int> data_io;
	int* orig_pixel = new int[Npixh * Npixw];
	char buf[512];
	for (int n = 0; n < ncams; n++) {
		for (int i = 0; i < Npixh; i++) {
			for (int j = 0; j < Npixw; j++) {
				orig_pixel[i * Npixw + j] = _pixels_orig[n][i][j];
			}
		}
		data_io.SetTotalNumber(Npixh * Npixw);
		sprintf_s(buf, "%s\\orig_F%.03d_C%d.txt", tiffaddress.c_str(), frame, n);
		//data_io.SetFilePath(tiffaddress + "Tracks\\orig" + to_string(frame) + "_" + to_string(n) + ".txt");
		data_io.SetFilePath(buf);
		data_io.WriteData(orig_pixel);
	}
	delete[] orig_pixel;
	// change to write directly to an image according to https://research.cs.wisc.edu/graphics/Courses/638-f1999/libtiff_tutorial.htm
#endif //TEST_PARTICLES_GENERATION_SAVE
}

Frame IPR::FindPos3D(deque< deque<string> > imgNames, int frameNumber) {

	// clearing all the ipr variables
	pos3D.clear();
	ghost3D.clear();
	intensity3D.clear();
	/*
	 * Modified by Shiyong Tan, 3/30/2018
	 * iframes should be cleared before writing new data into it.
	 * Start:
	 */
	iframes.clear();
	// End

	frame = frameNumber;
	// creating a filename with all tiff image names at a particular frame
	/* make sure that the images are in the same sequence as the camera numbers*/
	for (int i = 0; i < ncams; i++)
		filename.push_back(imgNames[i][frame - 1]);

	// converting the images to 2D dynamic array and finding 2D particle centers
	Tiff2DFinder t(ncams, threshold, filename);
	t.FillPixels(pixels_orig);	// pixels_orig will be filled with pixel intensities for all cameras

#ifdef TEST_PARTICLES_GENERATION
	//Replacing pixels_orig
	TestParticleGeneration(pixels_orig);
#endif // TEST_PARTICLES_GENERATION

	for (int camID = 0; camID < ncams; camID++) {		// filling iframes with the 2D positions on each camera
		try {
			ParticleFinder p(pixels_orig[camID], Npixh, Npixw);//, t.Get_colors(), threshold);
			if (debug_mode == SKIP_IPR_2D_POSITION && frame - 1 < debug_frame_number) { // read 2D position directly
				iframes.push_back(p.ReadParticle2DCenter(imgNames[camID][frame - 1]));
				if (error == NO_FILE) {
					cout << "The file for reading particle 2D center can't be opened!";
					error = NONE;
					goto Find2DCenter;
				}
			}
			else {
			Find2DCenter:
				p.GetParticle2DCenter(t.Get_colors(), threshold);
				iframes.push_back(p.CreateFrame());
				p.SaveParticle2DCenter(imgNames[camID][frame - 1]);
			}
		}
		catch (out_of_range & e) {
			cerr << e.what() << endl;
			throw runtime_error("Caught out_of_range in ImageSequence::Particle2DList()");
		}
	}

	//Load_2Dpoints("S:/Projects/Bubble/09.28.17/Bubbles and Particles - 250fps/frame", frame, ALL_CAMS);

	int ncams4 = (ncams > 4) ? 4 : ncams;  //No idea why it was set like this.
//	int ncams4 = ncams;
	Calibration calib(cams4, mindist_2D, mindist_3D, ncams4);
	OTF OTFcalib(ncams, otfFile);

	double del = 0.01; //mm in 3D

	clock_t start0;
	start0 = clock();

	// iteratively correcting the position of each particle (SHAKING)

	// ALL (OR FIRST 4) CAMERA OUTERLOOP STARTS HERE 
	deque<int> camNums;
	for (int i = 0; i < ncams4; i++)
		camNums.push_back(i);

	for (int loopOuter = 0; loopOuter < it_outerloop; loopOuter++) {

		// time
		clock_t start;
		start = clock();

		// increasing the 2D threshold by 10% in each iteration
		calib.Set_min2D(pow(1.1, loopOuter) * mindist_2D);

		// update pos3Dnew with 3D particles from IPR of current outerloop
		Frame pos3Dnew = IPRLoop(calib, OTFcalib, camNums, ALL_CAMS, t.Get_colors(), pixels_orig, pixels_reproj, pixels_res, loopOuter);

		// add the current outerloop particles to pos3D
		pos3D.insert(pos3D.end(), pos3Dnew.begin(), pos3Dnew.end());

		double duration = clock() - start;
		cout << "\t# of particles detected in outerloop" << loopOuter << ": " << pos3Dnew.NumParticles() << endl;
		cout << "\t Time taken for outerloop " << loopOuter << ": " << duration / (double)CLOCKS_PER_SEC << "s" << endl;
		cout << "\tTotal particles (" << pos3D.size() << ")" << endl;

	}
	//	exit(0);

	//	_getch();

	// REDUCED CAMERA LOOP (ignoring 1 camera) 
	if (reducedCams) {
		m_reduce_cam_begin = 1;
		cout << "\t\t\t\t Applying IPR for reduced cams" << endl;

		Calibration calibReduced(calibfile, REDUCED_CAMS, mindist_2Dr, mindist_3Dr, ncams4); // new calibration file for reduced cameras	

		for (int loopOuter = 0; loopOuter < it_reducedCam; loopOuter++) {
			// increasing the 2D threshold by 10% in each iteration
			calibReduced.Set_min2D(pow(1.1, loopOuter) * mindist_2Dr);

			// running IPR by ignoring cameras one by one
			for (int ignoreCam = 0; ignoreCam < ncams4; ignoreCam++) {
				Frame pos3Dnew = IPRLoop(calibReduced, OTFcalib, camNums, ignoreCam, t.Get_colors(), pixels_orig, pixels_reproj, pixels_res, loopOuter);
				cout << "\t # of particles detected in outerloop" << loopOuter << ", ignoring cam" << ignoreCam << ": " << pos3Dnew.NumParticles() << endl;

				pos3D.insert(pos3D.end(), pos3Dnew.begin(), pos3Dnew.end());
			}

			cout << "\t Total particles (" << pos3D.size() << ")" << endl;
		}
		m_reduce_cam_begin = 0;
	}

	double duration0 = clock() - start0;
	cout << "\t Total # of particles detected: " << pos3D.size() << endl;
	cout << "\t Total time taken by IPR: " << duration0 / (double)CLOCKS_PER_SEC << "s" << endl;

	// FOR TESTING 
		// saving pos3D as a .mat file
		//if (triangulationOnly || IPROnly) 
	{
		//string reproj = "ReprojectedImage", res = "ResidualImage" + to_string(frame);
		stringstream mat_name1; mat_name1 << tiffaddress << "pos3Dframe" << frameNumber << ".txt";
		//stringstream mat_name2; mat_name2 << tiffaddress << "cleanlistpos2Dframe" << frameNumber;
		SaveParticlePositions(pos3D, mat_name1.str());
		//SaveParticlePositions(calib.good2Dpos, mat_name2.str());
		//MatfilePositions(ghost3D, "ghost3D");
	}

	// saving the final residual and reprojected images as .mat files
	/*ReprojImage(pos3D, OTFcalib);
	// shifting the pixels of reprojected image: 1,1 --> 0,0
	for (int n = 0; n < ncams; n++) {
		for (int i = 1; i < Npixh; i++) {
			for (int j = 1; j < Npixw; j++) {
				pixels_reproj[n][i - 1][j - 1] = pixels_reproj[n][i][j];
			}
		}
	}*/
	// shifting the pixels of residual image: 1,1 --> 0,0
	/*for (int n = 0; n < ncams; n++) {
		for (int i = 1; i < Npixh; i++) {
			for (int j = 1; j < Npixw; j++) {
				pixels_res[n][i - 1][j - 1] = pixels_res[n][i][j];
			}
		}
	}*/

	//MatfileImage(pixels_reproj, reproj);
	//MatfileImage(pixels_res, res);

#ifdef TEST_PARTICLES_GENERATION
	//Check acuracy by synthetic images
	auto dist = [](const Position& a, const Position& b)
	{
		return std::hypot(std::hypot(a.X() - b.X(), a.Y() - b.Y()), a.Z() - b.Z());
	};
#undef max
	double error = 0;
	double error_x = 0;
	double error_y = 0;
	double error_z = 0;

	for (const auto& pos : pos3D)
	{
		double _min = std::numeric_limits<double>::max();
		double _dist = 0;
		Position p_nearest;
		for (const auto& rnd_p : rnd_points)
		{
			//_min = min(_min, dist(pos, rnd_p));
			_dist = dist(pos, rnd_p);
			if (_min > _dist)
			{
				_min = _dist;
				p_nearest = rnd_p;
			}
		}
		error += _min;
		error_x += pos.X() - p_nearest.X();
		error_y += pos.Y() - p_nearest.Y();
		error_z += pos.Z() - p_nearest.Z();
	}
	printf("TEST_PARTICLES_GENERATION mean error = %lf\n", error / pos3D.size());
	printf("TEST_PARTICLES_GENERATION mean error_x = %lf\n", error_x / pos3D.size());
	printf("TEST_PARTICLES_GENERATION mean error_y = %lf\n", error_y / pos3D.size());
	printf("TEST_PARTICLES_GENERATION mean error_z = %lf\n", error_z / pos3D.size());
	cin.get();
#endif // TEST_PARTICLES_GENERATION
	return pos3D;
}


//############################################################################# SUB-FUNCTIONS ##########################################################################

// a function that gives the 3D position by applying a single outerloop of IPR
Frame IPR::IPRLoop(Calibration& calib, OTF& OTFcalib, deque<int> camNums, int ignoreCam, double colors,
	deque<int**>& orig, deque<int**>& reproj, deque<int**>& res, int loop_time) {
	Frame pos3Dnew;
	double del;
	int id = 0;
	// getting 2D particles from the updated original image (residual after removing identified particles) on reduced cams
	// for the first loop we should not clear iframes
	// for the following loop of IPR in initial phase as well as every loop in the convergence phase, 
	// it should be refreshed by using residual images.
	if (loop_time >= 1 || m_reduce_cam_begin || m_doing_STB) {
		iframes.clear();
		for (int camID = 0; camID < camNums.size(); camID++)
			if (camNums[camID] != ignoreCam) {
				try {
					ParticleFinder p(orig[camNums[camID]], Npixh, Npixw);//, colors, threshold);
					p.GetParticle2DCenter(colors, threshold);
					//cout << "\nQQQQQQQQQQQQQQ = " << threshold << endl;
					iframes.push_back(p.CreateFrame());
					//p.SaveParticle2DCenter("/home/tanshiyong/Documents/Data/Single-Phase/11.03.17/Run1/frame100_" + to_string(camID) + ".txt");
				}
				catch (out_of_range & e) {
					cerr << e.what() << endl;
					throw runtime_error("Caught out_of_range in ImageSequence::Particle2DList()");
				}
			}
	}


	//	Load_2Dpoints("S:/Projects/Bubble/Cam_Config_of_10.22.17/10.29.17/BubblesNParticlesHigh_4000fps/BubblesNParicleswithBreakup/Bubble_Reconstruction_Corrected/Bubble_2D_centers", frame, ignoreCam);
	cout << "\tThe 2D particle number in a image is: ";
	int cNum = camNums.size();
	if (ignoreCam != ALL_CAMS) cNum--;
	for (int i = 0; i < cNum; i++)
		cout << iframes[i].NumParticles() << " ";
	cout << endl;
	//	cout << iframes[1].NumParticles() << endl;
	//	cout << iframes[2].NumParticles() << endl;
	//	cout << iframes[3].NumParticles() << endl;
		// stereomatching to give 3D positions from triangulation
	if ((debug_mode == SKIP_IPR_TRIANGULATION || debug_mode == SKIP_IPR_SHAKING)
		&& frame - 1 < debug_frame_number) {
		if (!m_reduce_cam_begin) {
			pos3Dnew = ReadParticlePositions(m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + ".txt");
		}
		else {
			pos3Dnew = ReadParticlePositions(m_particle_position_addr + "frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam)
				+ "Loop" + to_string(loop_time) + ".txt");
		}
		if (error == NO_FILE) {
			string message;
			if (!m_reduce_cam_begin) {
				message = "The triangulation file for frame" + to_string(frame) + "Loop" + to_string(loop_time) + "can't be opened!\n";
			}
			else {
				message = "The triangulation file for frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam)
					+ "Loop" + to_string(loop_time) + "can't be opened!\n";
			}
			cout << message;
			error = NONE;
			goto StereoMatch;
		}
		else {
			printf("Read in %i triangulated particles\n", pos3Dnew.NumParticles());
		}
	}
	else {
	StereoMatch:
		//double start = clock();
		pos3Dnew = calib.Stereomatch(iframes, frame, ignoreCam);
		//cout << "StereoMatchTime = " << (clock() - start) / (double)CLOCKS_PER_SEC << "\n";
		//cin.get();
		if (frame < 4) {  // we save the data only for the first 4 frames.
			if (!m_reduce_cam_begin) {
				SaveParticlePositions(pos3Dnew.Get_PosDeque(), m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + ".txt");
			}
			else {
				SaveParticlePositions(pos3Dnew.Get_PosDeque(),
					m_particle_position_addr + "frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam) + "Loop" + to_string(loop_time) + ".txt");
			}
		}

	}

	if (((!triangulationOnly) && IPROnly) || !(triangulationOnly || IPROnly)) {
		// initializing the 3D intensity
		deque<double> intensity3Dnew;
		if (debug_mode == SKIP_IPR_SHAKING && frame - 1 < debug_frame_number) {
			string file_path, file_path1;
			if (!m_reduce_cam_begin) {
				file_path = m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + "Shaking.txt";
				file_path1 = m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + "ShakingIntensity.txt";
			}
			else {
				file_path = m_particle_position_addr + "frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam) + "Loop" + to_string(loop_time) + "Shaking.txt";
				file_path1 = m_particle_position_addr + "frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam) + "Loop" + to_string(loop_time) + "ShakingIntensity.txt";
			}
			pos3Dnew = ReadParticlePositions(file_path);
			intensity3Dnew = ReadParticleIntensity(file_path1);
			if (error == NO_FILE) {
				string message;
				if (!m_reduce_cam_begin) {
					message = "The shaking file for frame" + to_string(frame) + "Loop" + to_string(loop_time) + "can't be opened!\n";
				}
				else {
					message = "The shaking file for frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam)
						+ "Loop" + to_string(loop_time) + "can't be opened!\n";
				}
				cout << message;
				error = NONE;
				goto Shaking;
			}
			else {
				printf("Read in %i Shaked particles\n", pos3Dnew.NumParticles());
			}
		}
		else {
		Shaking:
			for (unsigned int i = 0; i < pos3Dnew.NumParticles(); i++)
				intensity3Dnew.push_back(1.0);

			Frame estPosOld = pos3Dnew;

			// ######## Innerloop starts here ########
			for (int loopInner = 0; loopInner < it_innerloop; loopInner++) {

				//pair<int, int> bad = Rem(pos3Dnew, intensity3Dnew, mindist_3D);
//						double start = clock();

				if (pos3Dnew.NumParticles() == 0)
					break;

				// Creating the reprojected images (pixels_reproj) by reprojecting the 3D particles onto all cameras using Gaussian ellipse.
				ReprojImage(pos3Dnew, OTFcalib, reproj, IPRflag);
				//						ReprojImage(pos3Dnew, OTFcalib, reproj, 0.5);

										// residual image
#ifdef SAVE_INTERMEDIATE_PROJS
				NumDataIO<int> data_io;
				int* res_pixel = new int[Npixh * Npixw];
				int* orig_pixel = new int[Npixh * Npixw];
				int* proj_pixel = new int[Npixh * Npixw];
#endif
				for (int n = 0; n < camNums.size(); n++) {
					for (int i = 0; i < Npixh; i++) {
						for (int j = 0; j < Npixw; j++) {
							int residual = (orig[camNums[n]][i][j] - reproj[camNums[n]][i][j]);
							res[camNums[n]][i][j] = residual;// (residual < 0) ? 0 : residual;
#ifdef SAVE_INTERMEDIATE_PROJS
							res_pixel[i * Npixw + j] = residual;
							orig_pixel[i * Npixw + j] = orig[camNums[n]][i][j];
							proj_pixel[i * Npixw + j] = reproj[camNums[n]][i][j];
#endif
						}
					}
#ifdef SAVE_INTERMEDIATE_PROJS
					data_io.SetTotalNumber(Npixh* Npixw);

					data_io.SetFilePath(tiffaddress + "Tracks\\orig" + to_string(loopInner) + "_" + to_string(frame) + "_" + to_string(n) + ".txt");
					data_io.WriteData(orig_pixel);

					//data_io.SetFilePath(tiffaddress + "Tracks\\proj" + to_string(loopInner) + "_" + to_string(frame) + "_" + to_string(n) + ".txt");
					//data_io.WriteData(proj_pixel);
					
					//data_io.SetFilePath(tiffaddress + "Tracks\\res" + to_string(loopInner) + "_" + to_string(frame) + "_" + to_string(n) + ".txt");
					//data_io.WriteData(res_pixel);
#endif
				}
#ifdef SAVE_INTERMEDIATE_PROJS
				delete[] res_pixel;
				delete[] orig_pixel;
				delete[] proj_pixel;
#endif

				// updating the 3D position and intensity field by shaking
				Frame::const_iterator pIDend = pos3Dnew.end();
				//						int index = 0;

				vector<Position> shakeLims(pos3Dnew.NumParticles());
				double mindist_2D_ = 0.04;

				if (loopInner < 2)  del = mindist_2D_;
				else if (loopInner < 5)  del = mindist_2D_ / pow(2, loopInner - 1);
				else  del = mindist_2D_ / 100;
				//							if (loopInner < 1)  del = config.shaking_shift;			// initial shake TODO
				//							else if (loopInner < 5)  del = config.shaking_shift / pow(2,loopInner - 1);//_ipr.mindist_2D/10;	// normal shakes TODO
				//							else  del = config.shaking_shift/100;

				//cout << "del = " << del << "\n";

				// shaking
//						Frame::const_iterator pID = pos3Dnew.begin();
//						auto start = std::chrono::system_clock::now();
#ifndef SAVE_INTERMEDIATE_SHAKES
#pragma omp parallel //num_threads(8)
#endif
				{//2 4 6 8 10 20 15
//							int TID = omp_get_thread_num();
//							printf("Thread %d is runing\n", TID);
#ifndef SAVE_INTERMEDIATE_SHAKES
#pragma omp for
#endif
					for (int i = 0; i < pos3Dnew.NumParticles(); ++i) {
						//						for (Frame::const_iterator pID = pos3Dnew.begin(); pID != pIDend; ++pID) {
						Frame::const_iterator pID = pos3Dnew.begin() + i;
						OTF otf_calib(OTFcalib);
						Shaking s(ncams, ignoreCam, otf_calib, Npixw, Npixh, psize, del, *pID, camsAll, res, intensity3Dnew[i], loop_time, loopInner);

						pos3Dnew[i] = s.Get_posnew();
						intensity3Dnew[i] = s.Get_int();
						shakeLims[i] = s.Get_lims();
						//							index++;
						//							pID++;
					}
				}
				//						double duration = (clock() - start) / (double) CLOCKS_PER_SEC;
				//						cout<<"Time taken for shaking in each loop:"<<duration<<endl;
				//						auto end = std::chrono::system_clock::now();
				//						auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
				//						cout << elapsed.count() << '\n';

				Position totLims = accumulate(shakeLims.begin(), shakeLims.end(), Position(0, 0, 0));

				if (pos3Dnew.NumParticles() != 0) totLims /= pos3Dnew.NumParticles();
				else totLims = Position(0, 0, 0);

				cout << "Shake limits percentage: " << totLims.X() << " " << totLims.Y() << " " << totLims.Z() << "\n";

			} // ############# Innerloop ends here #############


			//ofstream f(tiffaddress + "Tracks\\testShake3.dat", ios::out);
			//for (int index = 0; index < pos3Dnew.NumParticles(); index++) {						// adding 2D image centers and intensity data to the estimates 
			//	Position pos3Dold = estPosOld[index];
			//	Position pos3Dnew_ = pos3Dnew[index];
			//	f << pos3Dold.X() - pos3Dnew_.X() << "\t" << pos3Dold.Y() - pos3Dnew_.Y() << "\t" << pos3Dold.Z() - pos3Dnew_.Z() << endl;
			//}

			// Save the data
			if (to_save_data) {
				string file_path, file_path1;
				if (!m_reduce_cam_begin) {
					file_path = m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + "Shaking.txt";
					file_path1 = m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + "ShakingIntensity.txt";
				}
				else {
					file_path = m_particle_position_addr + "frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam) + "Loop" + to_string(loop_time) + "Shaking.txt";
					file_path1 = m_particle_position_addr + "frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam) + "Loop" + to_string(loop_time) + "ShakingIntensity.txt";
				}
				SaveParticlePositions(pos3Dnew.Get_PosDeque(), file_path);
				SaveParticleIntensity(intensity3Dnew, file_path1);
			}
		}

		if (pos3Dnew.NumParticles() != 0) {
			//			cout <<  "Number of particles after shaking:" <<pos3Dnew.NumParticles()<<endl;
						// save the 3D position, 3D intensity and their correspoding 2D positions
			for (int index = 0; index < pos3Dnew.NumParticles(); index++)
				FullData(pos3Dnew[index], intensity3Dnew[index], camNums.size(), ignoreCam);

			// removing ambiguous and ghost particles
			Rem(pos3Dnew, intensity3Dnew, mindist_3D);

			// updating the reprojected image
			ReprojImage(pos3Dnew, OTFcalib, reproj, IPRflag);
			//			ReprojImage(pos3Dnew, OTFcalib, reproj, 1.5);

						// updating the original image by removing correctly identified 3D particles
			for (int n = 0; n < camNums.size(); n++)
				for (int i = 0; i < Npixh; i++) {
					for (int j = 0; j < Npixw; j++) {
						int residual = (orig[camNums[n]][i][j] - config.fpt * reproj[camNums[n]][i][j]);
						orig[camNums[n]][i][j] = (residual < 0) ? 0 : residual;
					}
				}
		}
	}
	filename.clear();
	return pos3Dnew;
}



void IPR::ReprojImage(Frame matched3D, OTF& OTFcalib, deque<int**>& pixels_reproj, bool STB) {
	int size = psize;
	// doubling the area of reprojection for STB
	// also double for IPR for consistency with Shaking.
//	if (STB)
	size = 2 * psize;

	// intializing pixel_reproj to 0
	for (int camID = 0; camID < ncams; camID++) {
		for (int i = 0; i < Npixh; i++) {
			for (int j = 0; j < Npixw; j++) {
				pixels_reproj[camID][i][j] = 0;
			}
		}
	}

	Frame::const_iterator pIDend = matched3D.end();
	int p = 0;
	for (Frame::const_iterator pID = matched3D.begin(); pID != pIDend; ++pID) {

		// stores the particle's center on each camera image
		deque<Position> particle2Dcenters;

		for (int n = 0; n < ncams; n++) {
			// Use OTFparam matrix and interpolate (trilinear) to get the parameters at the current 3D particle location
			vector <double> otfParam = OTFcalib.OTFgrid(n, matched3D[p]); // otfParam contains a,b,c and alpha for camera 'n' at 3D position 'pos'
			// finding the 2D center
			Position pos2Dmm = camsAll[n].WorldToImage(*pID); pos2Dmm.Set_Z(0);
			particle2Dcenters.push_back(camsAll[n].Distort(pos2Dmm));

			// *Reporjecting* //
			// pixel range for each particle
			int xmin = max(1, (int)floor(particle2Dcenters[n].X() - size / 2 + 1));
			int ymin = max(1, (int)floor(particle2Dcenters[n].Y() - size / 2 + 1));
			int xmax = min(Npixw, (int)floor(particle2Dcenters[n].X() + size / 2 + 1));
			int ymax = min(Npixh, (int)floor(particle2Dcenters[n].Y() + size / 2 + 1));

			for (int x = xmin; x < xmax; x++) {
				for (int y = ymin; y < ymax; y++) {
					// reprojecting the particle using Gaussian ellipse
					int proj = round(PixelReproj(particle2Dcenters[n], otfParam, x, y));
#ifndef TEST_PARTICLES_GENERATION
					pixels_reproj[n][y][x] = max(pixels_reproj[n][y][x], proj);// important comment: Not sure max is the right thing to use here for overlapping particles
#else
					pixels_reproj[n][y][x] = max(0, min(65535, pixels_reproj[n][y][x] + proj));
#endif //TEST_PARTICLES_GENERATION

				}
			}
		}
		p++;
	}

}

void IPR::ReprojImage(Frame matched3D, OTF& OTFcalib, deque<int**>& pixels_reproj, double projsize) {
	//	int size = psize;
	//	// doubling the area of reprojection for STB
	//	// also double for IPR for consistency with Shaking.
	////	if (STB)
	//		size = 2 * psize;
	projsize = projsize * psize;

	// intializing pixel_reproj to 0
	for (int camID = 0; camID < ncams; camID++) {
		for (int i = 0; i < Npixh; i++) {
			for (int j = 0; j < Npixw; j++) {
				pixels_reproj[camID][i][j] = 0;
			}
		}
	}

	Frame::const_iterator pIDend = matched3D.end();
	int p = 0;
	for (Frame::const_iterator pID = matched3D.begin(); pID != pIDend; ++pID) {

		// stores the particle's center on each camera image
		deque<Position> particle2Dcenters;

		for (int n = 0; n < ncams; n++) {
			// Use OTFparam matrix and interpolate (trilinear) to get the parameters at the current 3D particle location
			vector <double> otfParam = OTFcalib.OTFgrid(n, matched3D[p]); // otfParam contains a,b,c and alpha for camera 'n' at 3D position 'pos'
			// finding the 2D center
			Position pos2Dmm = camsAll[n].WorldToImage(*pID); pos2Dmm.Set_Z(0);
			particle2Dcenters.push_back(camsAll[n].Distort(pos2Dmm));

			// *Reporjecting* //
			// pixel range for each particle
			int xmin = max(1, (int)floor(particle2Dcenters[n].X() - projsize / 2));
			int ymin = max(1, (int)floor(particle2Dcenters[n].Y() - projsize / 2));
			int xmax = min(Npixw, (int)floor(particle2Dcenters[n].X() + projsize / 2));
			int ymax = min(Npixh, (int)floor(particle2Dcenters[n].Y() + projsize / 2));

			for (int x = xmin; x < xmax; x++) {
				for (int y = ymin; y < ymax; y++) {
					// reprojecting the particle using Gaussian ellipse
					int proj = round(PixelReproj(particle2Dcenters[n], otfParam, x, y));
					pixels_reproj[n][y][x] = max(pixels_reproj[n][y][x], proj);
					// important comment: Not sure max is the right thing to use here for overlapping particles
				}
			}
		}
		p++;
	}
}

// a function for Gaussian ellipse reprojection at position (x,y)
double IPR::PixelReproj(Position& particle2Dcenter, vector <double>& otfParam, int x, int y) {
	double xx = ((double)x - particle2Dcenter.X()) * cos(otfParam[3]) + ((double)y - particle2Dcenter.Y()) * sin(otfParam[3]);
	double yy = -((double)x - particle2Dcenter.X()) * sin(otfParam[3]) + ((double)y - particle2Dcenter.Y()) * cos(otfParam[3]);
	double value = otfParam[0] * exp(-(otfParam[1] * pow(xx, 2) + otfParam[2] * pow(yy, 2)));
	return(value);
}

// removing ghost and ambiguous particles
pair<int, int> IPR::Rem(Frame& pos3D, deque<double>& int3D, double mindist_3D) {
	int ambiguous = 0;

	// deleting particles that are very close to each other
	for (int i = 0; i < pos3D.NumParticles(); i++) {
		for (int j = i + 1; j < pos3D.NumParticles();) {
			if (Distance(pos3D[i], pos3D[j]) < 9 * mindist_3D * mindist_3D) {
				pos3D.Delete(j); int3D.erase(int3D.begin() + j);
				ambiguous++;
			}
			else
				j++;
		}
	}

	int ghost = 0;
	// deleting based on intensity
	double avgInt = 0;
	for (int i = 0; i < int3D.size(); i++) {
		if (int3D[i] > 0) {
			avgInt = avgInt + int3D[i];
		}
	}
	avgInt = avgInt / int3D.size();
	//	cout<<"average intensity:"<<avgInt<<endl;

	for (int index = int3D.size() - 1; index >= 0; index--) {
		// removing a particle if its intensity falls below a % of the avg intensity
		if (int3D[index] < intensityLower * avgInt) {
			ghost3D.push_back(pos3D[index]);
			pos3D.Delete(index); int3D.erase(int3D.begin() + index);
			ghost++;
		}
	}

	return make_pair(ambiguous, ghost);
}

void IPR::FullData(Position& pos, double intensity, int Cams, int ignoreCam) {
	deque< deque<double> > pos2D(4);
	for (int camID = 0; camID < 4; camID++) {
		if (camID < Cams) {
			deque<double> tmp2D(2);
			Position tmp = camsAll[camID].WorldToImage(pos); tmp.Set_Z(0);
			tmp = camsAll[camID].Distort(tmp);
			tmp2D[0] = tmp.X(); tmp2D[1] = tmp.Y();
			pos2D[camID] = tmp2D;
		}
		else {
			/*
			 * Modified by Shiyong Tan, 1/31/18
			 * Illegal initialization with {...}
			 * Start:
			 */
			 //			deque<double> tmp = { 0.0, 0.0 };
			deque<double> tmp(2);
			tmp[0] = 0.0; tmp[1] = 0.0;
			// End
			pos2D[camID] = tmp;
		}
	}
	Position temp(pos.X(), pos.Y(), pos.Z(), pos2D[0][0], pos2D[0][1], pos2D[1][0], pos2D[1][1], pos2D[2][0], pos2D[2][1], pos2D[3][0], pos2D[3][1], intensity);
	pos = temp;
}

// creates a matfile of positions
void IPR::SaveParticlePositions(deque<Position> pos, string file_path) {
	/*
	 * Modified by Shiyong Tan, 2/6/18
	 * Discard using matio, use DataIO instead
	 * Start:
	 */
	 //	// Create a .mat file with pos3D
	 //	size_t sizeofpos3D = pos.size();
	 //	mat_t    *matfp;
	 //	matvar_t *cell_array, *cell_element;
	 //	size_t dims[2] = { sizeofpos3D, 1 };
	 //	stringstream s; s << "pos3Dframe" << frame;
	 //	string mat_name = name + ".mat";
	 //	matfp = Mat_CreateVer(mat_name.c_str(), NULL, MAT_FT_DEFAULT);
	 //	switch (NULL == matfp) {
	 //		fprintf(stderr, "Error creating MAT file \"pos3D.mat\"!\n");
	 //		break;
	 //	}
	 //
	 //	cell_array = Mat_VarCreate(s.str().c_str(), MAT_C_CELL, MAT_T_CELL, 2, dims, NULL, 0);
	 //	if (NULL == cell_array) {
	 //		fprintf(stderr, "Error creating variable for 'pos3D'\n");
	 //	}
	 //	else {
	 //		for (int i = 0; i < sizeofpos3D; i++) {
	 //			dims[0] = 1;
	 //			dims[1] = 12;
	 //			double temp[12] = { pos[i].X(),pos[i].Y(),pos[i].Z(),pos[i].X1(),pos[i].Y1(),pos[i].X2(),pos[i].Y2(),pos[i].X3(),pos[i].Y3(),pos[i].X4(),pos[i].Y4(),pos[i].Info() };
	 //			cell_element = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, temp, 0);
	 //			switch (NULL == cell_element) {
	 //				fprintf(stderr, "Error creating cell element variable\n");
	 //				Mat_VarFree(cell_array);
	 //				Mat_Close(matfp);
	 //				break;
	 //			}
	 //			Mat_VarSetCell(cell_array, i, cell_element);
	 //		}
	 //	}
	 //
	 //	Mat_VarWrite(matfp, cell_array, MAT_COMPRESSION_NONE);
	 //	Mat_VarFree(cell_array);
	 //	Mat_Close(matfp);

		 // TODO: Check whether it works. Shiyong Tan
	size_t sizeofpos3D = pos.size();
	vector<double> array(sizeofpos3D * 12);
	//Convert Position into array
	Position2Array(pos, &array[0]);

	NumDataIO<double> data_io;
	data_io.SetFilePath(file_path);  // setting the path to save the data.
	data_io.SetTotalNumber(sizeofpos3D * 12);
	data_io.WriteData((double*)&array[0]);
	// END
}

Frame IPR::ReadParticlePositions(string file_path) {
	NumDataIO<double> data_io;
	data_io.SetFilePath(file_path);
	int num = data_io.GetTotalNumber(); //Get total data number
	vector<double> array(num);
	data_io.ReadData((double*)&array[0]);
	deque<Position> pos = Array2Position(num / 12, &array[0]);
	Frame frame(pos);
	return frame;
}

// creates a matfile for images
void IPR::MatfileImage(deque<int**>& pix, string name) {

	/*
	 * Modified by Shiyong Tan, 2/7/18
	 * Matio is discared. Use DataIO instead.
	 * Start:
	 */
	 //	size_t sizeofpos3D = pix.size();
	 //	mat_t    *matfp;
	 //	size_t dims[2] = { Npixw, 1 };
	 //	string mat_name = name + ".mat";
	 //	matfp = Mat_CreateVer(mat_name.c_str(), NULL, MAT_FT_DEFAULT);
	 //	switch (NULL == matfp) {
	 //		fprintf(stderr, "Error creating MAT file \"pos3D.mat\"!\n");
	 //		break;
	 //	}
	 //
	 //	for (int i = 0; i < sizeofpos3D; i++) {
	 //		dims[0] = Npixw;
	 //		dims[1] = 1;
	 //		stringstream cam;
	 //		cam << "cam" << i;
	 //		matvar_t **cell_element = new matvar_t*[Npixw];
	 //		matvar_t *cell_array = Mat_VarCreate(cam.str().c_str(), MAT_C_CELL, MAT_T_CELL, 2, dims, NULL, 0);
	 //		if (NULL == cell_array) {
	 //			fprintf(stderr, "Error creating variable for 'pos3D'\n");
	 //		}
	 //		else {
	 //			for (int j = 0; j < Npixw; j++) {
	 //				dims[0] = 1;
	 //				dims[1] = Npixh;
	 //				int* temp = new int[Npixh];
	 //				temp = pix[i][j];
	 //
	 //				cell_element[j] = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims, temp, 0);
	 //
	 //				switch (NULL == cell_element) {
	 //					fprintf(stderr, "Error creating cell element variable\n");
	 //					Mat_VarFree(cell_array);
	 //					Mat_Close(matfp);
	 //					break;
	 //				}
	 //				Mat_VarSetCell(cell_array, j, cell_element[j]);
	 //				//Mat_VarFree(cell_element);
	 //				//delete[] temp;
	 //			}
	 //		}
	 //		Mat_VarWrite(matfp, cell_array, MAT_COMPRESSION_NONE);
	 //		for (int j = 0; j < Npixw; j++) {
	 //			Mat_VarFree(cell_element[j]);
	 //		}
	 //		delete[] cell_element;
	 //		//Mat_VarFree(cell_array);
	 //	}
	 //	Mat_Close(matfp);
		 // TODO: check whether it works.
	size_t sizeofpos3D = pix.size();
	// get all the values of pixels into 3D matrix
	vector<int> tmp(sizeofpos3D * Npixw * Npixh);
	for (int i = 0; i < sizeofpos3D; i++) {
		for (int j = 0; j < Npixw; j++) {
			for (int k = 0; k < Npixh; k++) tmp[Npixw * Npixh * i + Npixh * j + k] = pix[i][j][k];
		}
	}
	NumDataIO<int> data_io;
	data_io.SetFilePath(name + ".txt");
	data_io.SetTotalNumber(sizeofpos3D * Npixw * Npixh);
	data_io.WriteData((int*)&tmp[0]);
	// End
}

void IPR::Load_2Dpoints(string path, int frame, int ignoreCam) {
	/*
	 * Modified by Shiyong Tan, 2/7/18
	 * Discard using matio, use DataIo istead.
	 * Start:
	 */
	 //	stringstream s; s << path << ".mat";
	 //	string file = s.str();
	 //
	 //	const char *fileName = file.c_str();
	 //	mat_t *mat = Mat_Open(fileName, MAT_ACC_RDONLY);
	 //
	 //	if (mat == NULL) {
	 //		cout << " error in reading the 2D pos mat file" << endl;
	 //	}
	 //
	 //	for (int id = 0; id < ncams; id++) {
	 //		if (id != ignoreCam) {
	 //			int k = id;
	 //
	 //			//if (id == 2)
	 //			//	k = 3;
	 //			//if (id == 3)
	 //			//	k = 2;
	 //			stringstream s; s << "frame_" << frame;
	 //
	 //			string varName = s.str();
	 //			Frame pos2D;
	 //
	 //			if (mat) {
	 //				//std::cout << "Open file to read\n\tmat == " << mat << "\n";
	 //
	 //				matvar_t *matVar = 0;
	 //				matVar = Mat_VarRead(mat, (char*)varName.c_str());
	 //
	 //				if (matVar) {
	 //					int rows;	int cols;
	 //					rows = matVar->dims[0]; cols = matVar->dims[1];
	 //					unsigned namesize = matVar->nbytes / matVar->data_size;
	 //					double *namedata = static_cast<double*>(matVar->data);
	 //					for (int i = 0; i < rows; i++) {
	 //						Position pos(namedata[(2*id)*rows + i], namedata[(2*id + 1)*rows + i], 0);
	 //						pos2D.Add(pos);
	 //					}
	 //				}
	 //				else
	 //					cout << "Cannot open mat variable\n";
	 //			}
	 //			else
	 //				cout << "Cannot open mat file\n";
	 //
	 //			iframes.push_back(pos2D);
	 //		}
	 //	}
	 //
	 //	Mat_Close(mat);
		 //TODO: to check whether it works.
	NumDataIO<double> data_io;
	data_io.SetFilePath(path + ".txt");
	int total_number = data_io.GetTotalNumber(); //get the total number of elements in txt file
	// suppose the data format is: X1, Y1, X2, Y2, X3, Y3, X4, Y4
	int rows = total_number / 8;
	//double points_array[rows][8];
	vector<double> points_array(rows * 8);
	data_io.ReadData((double*)&points_array[0]);
	for (int id = 0; id < ncams; id++) {
		Frame pos2D;
		if (id != ignoreCam) {
			for (int i = 0; i < rows; i++) {
				Position pos(points_array[8 * i + 2 * id], points_array[8 * i + 2 * id + 1], 0);
				pos2D.Add(pos);
			}
		}
		iframes.push_back(pos2D);
	}
	// End
}

void IPR::Position2Array(deque<Position> pos, double* array) {
	size_t sizeofpos3D = pos.size();
	for (int i = 0; i < sizeofpos3D; i++) {
		array[12 * i] = pos[i].X(); array[12 * i + 1] = pos[i].Y(); array[12 * i + 2] = pos[i].Z();
		array[12 * i + 3] = pos[i].X1(); array[12 * i + 4] = pos[i].Y1(); array[12 * i + 5] = pos[i].X2(); array[12 * i + 6] = pos[i].Y2();
		array[12 * i + 7] = pos[i].X3(); array[12 * i + 8] = pos[i].Y3(); array[12 * i + 9] = pos[i].X4(); array[12 * i + 10] = pos[i].Y4();
		array[12 * i + 11] = pos[i].Info();
	}

}

deque<Position> IPR::Array2Position(int num_particle, double* array) {
	deque<Position> pos;
	for (int i = 0; i < num_particle; i++) {
		Position point(array[12 * i], array[12 * i + 1], array[12 * i + 2], array[12 * i + 3], array[12 * i + 4], array[12 * i + 5], array[12 * i + 6],
			array[12 * i + 7], array[12 * i + 8], array[12 * i + 9], array[12 * i + 10], array[12 * i + 11]);
		pos.push_back(point);
	}
	return pos;
}

//Save Particle intensity
void IPR::SaveParticleIntensity(deque<double> intensity, string file_path) {
	NumDataIO<double> data_io;
	int num = intensity.size();
	data_io.SetFilePath(file_path);
	vector<double> intensity_array(num); // Convert deque into double
	for (int i = 0; i < num; i++) {
		intensity_array[i] = intensity[i];
	}
	data_io.SetTotalNumber(num);
	data_io.WriteData((double*)&intensity_array[0]);
}

// Read the particle intensity
deque<double> IPR::ReadParticleIntensity(string file_path) {
	NumDataIO<double> data_io;
	data_io.SetFilePath(file_path);
	int num = data_io.GetTotalNumber();
	vector<double> intensity_array(num);
	data_io.ReadData((double*)&intensity_array[0]);
	deque<double> intensity;
	for (int i = 0; i < num; i++) {
		intensity.push_back(intensity_array[i]);
	}
	return intensity;
}
