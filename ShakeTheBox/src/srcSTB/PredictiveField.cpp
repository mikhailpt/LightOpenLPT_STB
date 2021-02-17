
/*
 * Modified by Shiyong Tan, 2/1/18
 * cstring is used for C++ instead of string
 * Start:
 */
//#include <string>
#include <cstring>
// End
#include <deque>
#include <stdexcept>
#include <utility>

#include <linterp.h>
#include <PredictiveField.h>
#include "NumDataIO.h"
#include "Common.h"
#include "../../inc/Kriging/OrdinaryKriging.hpp"	//TODO set propper include in settings

#define PREV_FRAME 0
#define CURR_FRAME 1

using namespace std;

void PredictiveField::GetPredictiveField(Frame prevFramePos, Frame currFramePos, std::string& fname, int frame) {
	matchedPrev = prevFramePos;
	matchedCurr = currFramePos;
//TODO: report err when cannot open file.
	// remove comments from the file
	ifstream infile(fname.c_str(), ios::in);
	string line;
	stringstream parsed;
	while (getline(infile, line)) {
		size_t commentpos = line.find('#');
		if (commentpos > 0) {
			if (commentpos < string::npos) 
				line.erase(commentpos);
			
			parsed << line << '\t';
		}
	}
	infile.close();

	
	viewAreaLimits[0] = config.x_lower_limit; viewAreaLimits[1] = config.x_upper_limt;
	viewAreaLimits[2] = config.y_lower_limit; viewAreaLimits[3] = config.y_upper_limt;
	viewAreaLimits[4] = config.z_lower_limit; viewAreaLimits[5] = config.z_upper_limt;

	for (int i = 0; i < 3; i++) {
		parsed >> gridSize[i];
		gridSize[i] = gridSize[i] * config.factor;
	}
	
	parsed >> radius;
	radius = radius * config.factor;

	int mat; bool matlabFlag = false;
	parsed >> mat;
	if (mat != 0)
		matlabFlag = true;

	getPredictiveField = true;

	// saving the grid points
	GridPoints();

	if (matlabFlag) {
		string path;
		parsed >> path;
		cout << "\t\tLoading predictive field from matfile" << endl;
		// geting the field from matfile
		// TODO: to check whether it works.
		Load_field(path);
	}
	else {
		cout << "\t\tCalculating predictive field" << endl;
		// calculating the field
		Field();
//		cout<<field[0][1];
	}



	//########### FOR TESTING ###########
	// getting the displacement of particles
	/*for (Frame::const_iterator pID = matchedPrev.begin(); pID != matchedPrev.end(); ++pID) {
		Position pVel(ParticleInterpolation(*pID));
		vector<double> pDisplacement;
		pDisplacement.push_back(pVel.X()); pDisplacement.push_back(pVel.Y()); pDisplacement.push_back(pVel.Z());
		pDisplacement.push_back(pID->X()); pDisplacement.push_back(pID->Y()); pDisplacement.push_back(pID->Z());
		particleDisplacements.push_back(pDisplacement);
	}*/

	//saving the field
	//string grid_name = "grids", field_name = "field", pDisp_name = "particleDisplacement";
	//cout << "\t\tsaving field.mat" << endl;
	////MatfileSave(gridPoints, grid_name);
	//MatfileSave(field, field_name);
	//MatfileSave(particleDisplacements, pDisp_name);

}

void PredictiveField::myLoad_field(string path, string& fname) {

	//TODO: report err when cannot open file.
		// remove comments from the file
	ifstream infile(fname.c_str(), ios::in);
	string line;
	stringstream parsed;
	while (getline(infile, line)) {
		size_t commentpos = line.find('#');
		if (commentpos > 0) {
			if (commentpos < string::npos)
				line.erase(commentpos);

			parsed << line << '\t';
		}
	}
	infile.close();

	viewAreaLimits[0] = config.x_lower_limit; viewAreaLimits[1] = config.x_upper_limt;
	viewAreaLimits[2] = config.y_lower_limit; viewAreaLimits[3] = config.y_upper_limt;
	viewAreaLimits[4] = config.z_lower_limit; viewAreaLimits[5] = config.z_upper_limt;

	for (int i = 0; i < 3; i++) {
		parsed >> gridSize[i];
		gridSize[i] = gridSize[i] * config.factor;
	}

	parsed >> radius;
	radius = radius * config.factor;

	// saving the grid points
	GridPoints();

	// load data
	string fnameTot = (path + "_af1.dat");
	ifstream f(fnameTot);

	if (!f.is_open()) {
		error = NO_FILE;
		return;
	}

	int i, status;

	for (int z = 0; z < numGrids[2]; z++) {
		for (int y = 0; y < numGrids[1]; y++) {
			for (int x = 0; x < numGrids[0]; x++) {
				i = x * numGrids[2] * numGrids[1] + y * numGrids[2] + z;
				f >> gridPoints[0][i];
				f >> gridPoints[1][i];
				f >> gridPoints[2][i];
				f >> field[0][i];
				f >> field[1][i];
				f >> field[2][i];
				f >> field[3][i];
				//cout << gridPoints[0][i] << "\t" << gridPoints[1][i] << "\t" << gridPoints[2][i] << "\t" << field[0][i] << "\t" << field[1][i] << "\t" << field[2][i] << "\t" << field[3][i] << "\n";
				//cin >> status;
			}
		}
	}

	cout << "Predictive field from file " << fnameTot << " was loaded!" << endl;
}

void PredictiveField::GridPoints() {

	for (int i = 0; i < 3; i++) {
		numGrids[i] = ((viewAreaLimits[2 * i + 1] - viewAreaLimits[2 * i]) / gridSize[i]) + 1;
	}

	totalGridPoints = numGrids[0] * numGrids[1] * numGrids[2];
	gridPoints = new vector<double>[3];

	for (int n = 0; n < 4; n++) {
		field.push_back(new double[totalGridPoints]);
	}

	for (double x = viewAreaLimits[0]; x <= viewAreaLimits[1]; x = x + gridSize[0]) {
		for (double y = viewAreaLimits[2]; y <= viewAreaLimits[3]; y = y + gridSize[1]) {
			for (double z = viewAreaLimits[4]; z <= viewAreaLimits[5]; z = z + gridSize[2]) {
				gridPoints[0].push_back(x);
				gridPoints[1].push_back(y);
				gridPoints[2].push_back(z);
			}
		}
	}

}

void PredictiveField::Field() {

	double rsqr = pow(radius, 2);
	int dispMapRes = 10;
	double radius_ = 0.5;// maximum displacement in mm
	Size = 4 * dispMapRes * radius_ + 1;
	// converting the Map index (i) to displacement (dx) as i = m*dx + c;
	m = dispMapRes;
	c = (Size - 1) / 2;

#pragma omp parallel //num_threads(8)
						{
#pragma omp for
	for (int i = 0; i < totalGridPoints; i++) {

		double*** dispMap = new double**[Size];
		for (int j = 0; j < Size; j++) {
			dispMap[j] = new double*[Size];
				for (int k = 0; k < Size; k++) {
					dispMap[j][k] = new double[Size];
				}
		}

		deque<Frame::const_iterator> prevFrame;
		deque<Frame::const_iterator> currFrame;
		deque<deque<double>> displacements;
		
		// getting the points within the interrogation sphere around the grid point
		// previous frame
		SearchVolumeParticles(prevFrame, rsqr, gridPoints, i, PREV_FRAME);
		// current frame
		SearchVolumeParticles(currFrame, rsqr, gridPoints, i, CURR_FRAME);

		// getting the displacement vectors (if there are particles in the search volume for both the frames)
		if (currFrame.size() != 0 && prevFrame.size() != 0) {
			for (int curr = 0; curr < currFrame.size(); curr++) {
				for (int prev = 0; prev < prevFrame.size(); prev++) {
					deque<double> disp(4); // dx,dy,dz,I_curr*I_prev
					disp[0] = currFrame[curr]->X() - prevFrame[prev]->X();
					disp[1] = currFrame[curr]->Y() - prevFrame[prev]->Y();
					disp[2] = currFrame[curr]->Z() - prevFrame[prev]->Z();
					disp[3] = 1;// currFrame[curr]->Info() * prevFrame[prev]->Info();
					if ( fabs(disp[0]) < radius_ && fabs(disp[1]) < radius_ && fabs(disp[2]) < radius_ )
						displacements.push_back(disp);
				}
			}

			// intializing the displacement map to 0
			for (int j = 0; j < Size; j++) {
				for (int k = 0; k < Size; k++) {
					memset(dispMap[j][k], 0, sizeof(dispMap[0][0][0])*Size);
				}
			}

			// getting the 3D displacement correlation map
			DisplacementMap(dispMap, displacements);

			// finding the peak of displacement map
			deque<double> peak(3);
			peak = DisplacementMapPeak(dispMap);
			// saving the peak dx,dy,dz to field
			field[0][i] = peak[0];
			field[1][i] = peak[1];
			field[2][i] = peak[2];
			field[3][i] = displacements.size();
		}
		else {
			deque<double> peak(3, 0);
			field[0][i] = peak[0];
			field[1][i] = peak[1];
			field[2][i] = peak[2];
			field[3][i] = 0;
		}
		for (int j = 0; j < Size; j++) {
			for (int k = 0; k < Size; k++) {
				delete[] dispMap[j][k];
			}
			delete[] dispMap[j];
		}
		delete[] dispMap;
	}
						}


	
}

void PredictiveField::SearchVolumeParticles(deque<Frame::const_iterator>& Frame, double rsqr, vector<double>* gridPoint, int grid, int frame) {
	double x = gridPoint[0][grid], y = gridPoint[1][grid], z = gridPoint[2][grid];
	Position Grid(gridPoint[0][grid], gridPoint[1][grid], gridPoint[2][grid]);
	if (frame == PREV_FRAME) {
		for (Frame::const_iterator pID = matchedPrev.begin(); pID < matchedPrev.end(); ++pID) {
			if (pID->X() < x + radius && pID->X() > x - radius && pID->Y() < y + radius && pID->Y() > y - radius && pID->Z() < z + radius && pID->Z() > z - radius) {
				double distsqr = Distance(*pID, Grid);
				if (distsqr < rsqr) {
					Frame.push_back(pID);
				}
			}
		}
	}
	else if (frame == CURR_FRAME) {
		for (Frame::const_iterator pID = matchedCurr.begin(); pID < matchedCurr.end(); ++pID) {
			if (pID->X() < x + radius && pID->X() > x - radius && pID->Y() < y + radius && pID->Y() > y - radius && pID->Z() < z + radius && pID->Z() > z - radius) {
				double distsqr = pow((pID->X() - x), 2) + pow((pID->Y() - y), 2) + pow((pID->Z() - z), 2);
				if (distsqr < rsqr) {
					Frame.push_back(pID);
				}
			}
		}
	}
}

void PredictiveField::DisplacementMap(double*** &dispMap, deque<deque<double>> displacements) {
	
	for (int d = 0; d < displacements.size(); d++) {
		double dx = m*displacements[d][0] + c, dy = m * displacements[d][1] + c, dz = m * displacements[d][2] + c, I = displacements[d][3];
		int x = round(dx), y = round(dy), z = round(dz);
		for (int i = max(0, x - 1); i <= min(Size - 1, x + 1); i++) {
			double xx = i - dx;
			for (int j = max(0, y - 1); j <= min(Size - 1, y + 1); j++) {
				double yy = j - dy;
				for (int k = max(0, z - 1); k <= min(Size - 1, z + 1); k++) {
					double zz = k - dz;
					dispMap[i][j][k] = dispMap[i][j][k] + I*exp(-(pow(xx, 2) + pow(yy, 2) + pow(zz, 2)));
				}
			}
		}	
	}
}

deque<double> PredictiveField::DisplacementMapPeak(double*** &dispMap) {
	deque<int> index(3);
	// finding the index of largest element
	index = IndexofLargestElement(dispMap);
	int x = index[0], y = index[1], z = index[2];

	deque<double> peak(3);
	// finding the peak using 1D Gaussian in x, y & z direction
	if (x == 0)
		peak[0] = 0;
	else if (x == Size - 1)
		peak[0] == Size - 1;
	else
		peak[0] = Gaussian1DPeak(x - 1, dispMap[x - 1][y][z], x, dispMap[x][y][z], x + 1, dispMap[x + 1][y][z]);

	if (y == 0)
		peak[1] = 0;
	else if (y == Size - 1)
		peak[1] == Size - 1;
	else
		peak[1] = Gaussian1DPeak(y - 1, dispMap[x][y - 1][z], y, dispMap[x][y][z], y + 1, dispMap[x][y + 1][z]);
	
	if (z == 0)
		peak[2] = 0;
	else if (z == Size - 1)
		peak[2] == Size - 1;
	else
		peak[2] = Gaussian1DPeak(z - 1, dispMap[x][y][z - 1], z, dispMap[x][y][z], z + 1, dispMap[x][y][z + 1]);
	
	return peak;
	
}

// returns index of largest element
deque<int> PredictiveField::IndexofLargestElement(double*** &array) {

	deque<int> maxInd(3);
	double maximum = array[0][0][0];
	for (int x = 0; x < Size; x++) {
		for (int y = 0; y < Size; y++) {
			for (int z = 0; z < Size; z++) {
				if (array[x][y][z] > maximum) {
					maximum = array[x][y][z];
					maxInd[0] = x; maxInd[1] = y; maxInd[2] = z;
				}
			}
		}
	}

	return maxInd;
}

// returns 1D Gaussian peak
double PredictiveField::Gaussian1DPeak(double y1, double v1, double y2, double v2, double y3, double v3) {
	double lnz1, lnz2, lnz3;
	if (v1 == 0) {
			lnz1 = log(0.0001);
		}
		else {
			lnz1 = log(v1);
		}
		if (v2 == 0) {
			lnz2 = log(0.0001);
		}
		else {
			lnz2 = log(v2);
		}
		if (v3 == 0) {
			lnz3 = log(0.0001);
		}
		else {
			lnz3 = log(v3);
		}
	
	

	double yc = -0.5 * ((lnz1 * ((y2 * y2) - (y3 * y3))) - (lnz2 * ((y1 * y1) - (y3 * y3))) +
		(lnz3 * ((y1 * y1) - (y2 * y2))))
		/ ((lnz1 * (y3 - y2)) - (lnz3 * (y1 - y2)) + (lnz2 * (y1 - y3)));

	//cout << "xc: " << xc << " yc: " << yc << endl;
	return (yc - c)/m;
}

/*void PredictiveField::GridInterpolation() {

	vector< vector<double>::iterator > grid_iter_list;
	grid_iter_list.push_back(gridPoints[0].begin());
	grid_iter_list.push_back(gridPoints[1].begin());
	grid_iter_list.push_back(gridPoints[2].begin());

	// the size of grid in each dimension
	array<int, 3> grid_sizes;
	grid_sizes[0] = numGrids[0];
	grid_sizes[1] = numGrids[1];
	grid_sizes[2] = numGrids[2];

	// total number of elements
	int num_elements = grid_sizes[0] * grid_sizes[1] * grid_sizes[2];

	// construct the interpolator, the last two arguments are pointers to the underlying data
	InterpMultilinear<3, double> field_x(grid_iter_list.begin(), grid_sizes.begin(), field[0], field[0] + num_elements);
	InterpMultilinear<3, double> field_y(grid_iter_list.begin(), grid_sizes.begin(), field[1], field[1] + num_elements);
	InterpMultilinear<3, double> field_z(grid_iter_list.begin(), grid_sizes.begin(), field[2], field[2] + num_elements);


	vector<double> displacement(3);
	array<double, 3> pos1 = { pos3D.X(), pos3D.Y(), pos3D.Z() };
	displacement[0] = field_x.interp(pos1.begin());
	displacement[1] = field_y.interp(pos1.begin());
	displacement[2] = field_z.interp(pos1.begin());


	return displacement;
}*/

Position PredictiveField::TrilinearInterpolation(Position pos3D) {

	std::vector<double> gridx = linspace(viewAreaLimits[0], viewAreaLimits[1], numGrids[0]);
	std::vector<double> gridy = linspace(viewAreaLimits[2], viewAreaLimits[3], numGrids[1]);
	std::vector<double> gridz = linspace(viewAreaLimits[4], viewAreaLimits[5], numGrids[2]);
	vector< vector<double>::iterator > grid_iter_list;
	grid_iter_list.push_back(gridx.begin());
	grid_iter_list.push_back(gridy.begin());
	grid_iter_list.push_back(gridz.begin());


	// the size of grid in each dimension
	array<int, 3> grid_sizes;
	grid_sizes[0] = numGrids[0];
	grid_sizes[1] = numGrids[1];
	grid_sizes[2] = numGrids[2];

	// total number of elements
	int num_elements = grid_sizes[0] * grid_sizes[1] * grid_sizes[2];

	// construct the interpolator, the last two arguments are pointers to the underlying data
	InterpMultilinear<3, double> field_x(grid_iter_list.begin(), grid_sizes.begin(), field[0], field[0] + num_elements);
	InterpMultilinear<3, double> field_y(grid_iter_list.begin(), grid_sizes.begin(), field[1], field[1] + num_elements);
	InterpMultilinear<3, double> field_z(grid_iter_list.begin(), grid_sizes.begin(), field[2], field[2] + num_elements);


	
	array<double, 3> pos1 = { pos3D.X(), pos3D.Y(), pos3D.Z() };
	Position displacement(field_x.interp(pos1.begin()), field_y.interp(pos1.begin()), field_z.interp(pos1.begin()));

	return displacement;
}

std::vector<Position> PredictiveField::KrigingInterpolation(const std::vector<Position>& estimateCord_)
{

	std::vector<double> gridx = linspace(viewAreaLimits[0], viewAreaLimits[1], numGrids[0]);
	std::vector<double> gridy = linspace(viewAreaLimits[2], viewAreaLimits[3], numGrids[1]);
	std::vector<double> gridz = linspace(viewAreaLimits[4], viewAreaLimits[5], numGrids[2]);

	// the size of grid in each dimension
	array<size_t, 3> grid_sizes;
	grid_sizes[0] = numGrids[0];
	grid_sizes[1] = numGrids[1];
	grid_sizes[2] = numGrids[2];

	// total number of elements
	size_t num_elements = grid_sizes[0] * grid_sizes[1] * grid_sizes[2];

	auto unpackPositionCord = [gridx, gridy, gridz, grid_sizes](size_t index)
	{
		size_t x[3];
		for (size_t i = 2; i < 3; --i)
		{	//� �� ������ � ����� ������� ���������� ������� (� ���� �������� � field)
			x[i] = index % grid_sizes[i];
			index /= grid_sizes[i];
		}
		assert(index == 0);
		return Position(gridx[x[0]], gridy[x[1]], gridz[x[2]]);
	};

	auto unpackPositionVel = [this](size_t index)
	{
		return Position(field[0][index], field[1][index], field[2][index]);
	};

	using Element = std::pair<Position, Position>;	//cordinates and velocity
	auto distance = [](const Element& a, const Element& b)
	{
		auto delta = a.first - b.first;
		return std::hypot(std::hypot(delta.X(), delta.Y()), delta.Z());	//Linear Kriging
	};

	auto access = [](Element& a)->Element::second_type&
	{
		return a.second;
	};

	std::vector<Element> estimate(estimateCord_.size());
	for (size_t i = 0; i < estimate.size(); ++i)
	{
		estimate[i] = std::make_pair(estimateCord_[i], Position());
	}

	std::vector<Element> data;
	data.reserve(num_elements);
	for (size_t i = 0; i < num_elements; ++i)
	{
		Element e = std::make_pair(unpackPositionCord(i), unpackPositionVel(i));
		if (e.second == Position(0,0,0))
		{
			continue;
		}
		bool close = false;
		for (size_t j = 0; j < estimate.size(); ++j)
		{
			if (distance(e, estimate[j]) < 3.0)
			{
				close = true;
				break;
			}
		}
		if (close)
		{
			data.push_back(e);
		}
	}

	Kriging::OrdinaryKriging(data.cbegin(), data.cend(), estimate.begin(), estimate.end(), distance, access);

	std::vector<Position> result(estimate.size());
	for (size_t i = 0; i < result.size(); ++i)
	{
		result[i] = estimate[i].second;
	}

	return result;
}

Position PredictiveField::KrigingInterpolation(const Position& p)
{
	std::vector<Position> v(1);
	v[0] = p;
	return KrigingInterpolation(v)[0];
}

// return an evenly spaced 1-corrSize grid of doubles.
std::vector<double> PredictiveField::linspace(double first, double last, int len) {
	std::vector<double> result(len);
	double step = (last - first) / (len - 1);
	for (int i = 0; i<len; i++) { result[i] = first + i*step; }
	return result;
}


// ########################################### MAT FILES #########################################
void PredictiveField::MatfileSave(vector<double*> pos, string name) {
	/*
	 * Modified by Shiyong Tan, 2/7/18
	 * Discard using matio, use DataIO instead.
	 */
//	// Create a .mat file with pos3D
//	size_t sizeofpos3D = totalGridPoints;
//	mat_t    *matfp;
//	matvar_t *cell_array, *cell_element;
//	size_t dims[2] = { sizeofpos3D, 1 };
//	string mat_name = name + ".mat";
//	matfp = Mat_CreateVer(mat_name.c_str(), NULL, MAT_FT_DEFAULT);
//	switch (NULL == matfp) {
//		fprintf(stderr, "Error creating MAT file \"pos3D.mat\"!\n");
//		break;
//	}
//
//	cell_array = Mat_VarCreate(name.c_str(), MAT_C_CELL, MAT_T_CELL, 2, dims, NULL, 0);
//	if (NULL == cell_array) {
//		fprintf(stderr, "Error creating variable for 'pos3D'\n");
//	}
//	else {
//		for (int i = 0; i < sizeofpos3D; i++) {
//			dims[0] = 1;
//			dims[1] = 3;
//			double temp[3] = { pos[0][i],pos[1][i],pos[2][i] };
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
// TODO: to check whether it works
	// Convert pos into 2D matrix
	size_t sizeofpos3D = totalGridPoints;
	vector<double> point_array(sizeofpos3D*3);
	for (int i = 0; i < sizeofpos3D; i++) {
		point_array[3 * i] = pos[0][i];
		point_array[3 * i + 1] = pos[1][i];
		point_array[3 * i + 2] = pos[2][i];
	}

	NumDataIO<double> data_io;
	data_io.SetFilePath(name + ".txt");
	data_io.SetTotalNumber(sizeofpos3D * 3);
	data_io.WriteData((double*) &point_array[0]);
	// END
}

void PredictiveField::MatfileSave(vector<double> pos[3], string name) {
	/*
	 * Modified by Shiyong Tan, 2/7/18
	 * Discard using matio, use DataIO instead.
	 */
//	// Create a .mat file with pos3D
//	size_t sizeofpos3D = totalGridPoints;
//	mat_t    *matfp;
//	matvar_t *cell_array, *cell_element;
//	size_t dims[2] = { sizeofpos3D, 1 };
//	string mat_name = name + ".mat";
//	matfp = Mat_CreateVer(mat_name.c_str(), NULL, MAT_FT_DEFAULT);
//	switch (NULL == matfp) {
//		fprintf(stderr, "Error creating MAT file \"pos3D.mat\"!\n");
//		break;
//	}
//
//	cell_array = Mat_VarCreate(name.c_str(), MAT_C_CELL, MAT_T_CELL, 2, dims, NULL, 0);
//	if (NULL == cell_array) {
//		fprintf(stderr, "Error creating variable for 'pos3D'\n");
//	}
//	else {
//		for (int i = 0; i < sizeofpos3D; i++) {
//			dims[0] = 1;
//			dims[1] = 3;
//			double temp[3] = { pos[0][i],pos[1][i],pos[2][i] };
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
	// Convert pos into 2D matrix
	size_t sizeofpos3D = totalGridPoints;
	vector<double> point_array(sizeofpos3D*3);
	for (int i = 0; i < sizeofpos3D; i++) {
		point_array[3 * i] = pos[0][i];
		point_array[3 * i + 1] = pos[1][i];
		point_array[3 * i + 2] = pos[2][i];
	}

	NumDataIO<double> data_io;
	data_io.SetFilePath(name + ".txt");
	data_io.SetTotalNumber(sizeofpos3D * 3);
	data_io.WriteData((double*) &point_array[0]);
	// END
}
//
void PredictiveField::MatfileSave(vector<vector<double>> pos, string name) {
	/*
	 * Modified by Shiyong Tan, 2/7/18
	 * Discard using matio, use DataIO instead.
	 */
//	// Create a .mat file with pos3D
//	size_t sizeofpos3D = pos.size();
//	mat_t    *matfp;
//	matvar_t *cell_array, *cell_element;
//	size_t dims[2] = { sizeofpos3D, 1 };
//	string mat_name = name + ".mat";
//	matfp = Mat_CreateVer(mat_name.c_str(), NULL, MAT_FT_DEFAULT);
//	switch (NULL == matfp) {
//		fprintf(stderr, "Error creating MAT file \"pos3D.mat\"!\n");
//		break;
//	}
//
//	cell_array = Mat_VarCreate(name.c_str(), MAT_C_CELL, MAT_T_CELL, 2, dims, NULL, 0);
//	if (NULL == cell_array) {
//		fprintf(stderr, "Error creating variable for 'pos3D'\n");
//	}
//	else {
//		for (int i = 0; i < sizeofpos3D; i++) {
//			dims[0] = 1;
//			dims[1] = 6;
//			double temp[6] = { pos[i][0],pos[i][1],pos[i][2], pos[i][3],pos[i][4],pos[i][5] };
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
	size_t sizeofpos3D = totalGridPoints;
	vector<double> point_array(sizeofpos3D*6);
	for (int i = 0; i < sizeofpos3D; i++) {
		point_array[6 * i] = pos[i][0];
		point_array[6 * i + 1] = pos[i][1];
		point_array[6 * i + 2] = pos[i][2];
		point_array[6 * i + 3] = pos[i][3];
		point_array[6 * i + 4] = pos[i][4];
		point_array[6 * i + 5] = pos[i][5];
	}

	NumDataIO<double> data_io;
	data_io.SetFilePath(name + ".txt");
	data_io.SetTotalNumber(sizeofpos3D * 6);
	data_io.WriteData((double*) &point_array[0]);
	// END
}

//// Load field from field.mat
void PredictiveField::Load_field(string path) {
	/*
	 * Modified by Shiyong Tan, 2/7/18
	 * Discard using matio, use DataIO instead.
	 */
//	string file = path + to_string(frame) + ".mat";
//	const char *fileName = file.c_str();
//	mat_t *mat = Mat_Open(fileName, MAT_ACC_RDONLY);
//
//	if (mat == NULL)
//		cout << " error in reading the mat file at" << path << frame << ".mat" << endl;
//
//	string cam;
//	string name;
//
//	if (mat) {
//		//std::cout << "Open file to read\n\tmat == " << mat << "\n";
//
//		matvar_t *matVar = 0;
//		name = path;
//		matVar = Mat_VarRead(mat, (char*)name.c_str());
//
//		if (matVar) {
//			int rows;	int cols;
//			rows = matVar->dims[0]; cols = matVar->dims[1];
//			//totalGridPoints = rows;
//
//			unsigned namesize = matVar->nbytes / matVar->data_size;
//			double *namedata = static_cast<double*>(matVar->data);
//			for (int i = 0; i < rows; i++) {
//				field[0][i] = namedata[i];
//				field[1][i] = namedata[rows + i];
//				field[2][i] = namedata[2 * rows + i];
//			}
//		}
//	}
//	else {
//		cout << "Cannot open pfield mat file\n";
//	}
//
//	Mat_Close(mat);
//	string file = path + to_string(frame) + ".txt";

	NumDataIO<double> data_io;
	data_io.SetFilePath(path);
	int total_num = data_io.GetTotalNumber();
	vector<double> field_data(total_num);  // data format: rows * 3
	data_io.SetTotalNumber(total_num);
	data_io.ReadData((double*) &field_data[0]);
	for (int i = 0; i < total_num /3; i++) {
		field[0][i] = field_data[3 * i];
		field[1][i] = field_data[3 * i + 1];
		field[2][i] = field_data[3 * i + 2];
	}
//END
}

//void PredictiveField::SaveField(string file_path) {
//	NumDataIO<double> data_io;
//	data_io.SetFilePath(file_path);
//	data_io.SetTotalNumber(totalGridPoints * 6);
//	vector<double> field_data(totalGridPoints*6);
//	for (int i = 0; i < totalGridPoints; i++) {
//		field_data[3 * i]     = gridPoints[0][i];
//		field_data[3 * i + 1] = gridPoints[1][i];
//		field_data[3 * i + 2] = gridPoints[2][i];
//		field_data[3 * i + 3] = field[0][i];
//		field_data[3 * i + 4] = field[1][i];
//		field_data[3 * i + 5] = field[2][i];
//	}
//	data_io.WriteData((double*)&field_data[0]);
//}

void PredictiveField::SaveField(string file_path)
{
	int zeroCount = 0;
	for (int i = 0; i < totalGridPoints; i++)
		if (field[0][i] == 0 || field[1][i] == 0 || field[2][i] == 0)
			zeroCount++;

	int nonZeroCount = totalGridPoints - zeroCount;

	float ranges[3];
	for (int i = 0; i < 3; i++) {
		ranges[i] = viewAreaLimits[2 * i + 1] - viewAreaLimits[2 * i];
	}

	float meanDist = pow(ranges[0] * ranges[1] * ranges[2] / nonZeroCount, 1.0/3);

	string str = string("Percent of nonzeros is " + to_string(100*float(nonZeroCount)/ totalGridPoints) + "\% or " + to_string(nonZeroCount) + " of " + to_string(totalGridPoints) + ", with mean vector distance " + to_string(meanDist) + "mm");

	std::ofstream f(file_path);
	f  << "TITLE     = \"" + str + "\"\nVARIABLES = \"x [mm]\"\n\"y [mm]\"\n\"z [mm]\"\n\"dx [mm]\"\n\"dy [mm]\"\n\"dz [mm]\"\n\"particlesNum [units]\"\nZONE T=\"Rectangular zone\"\nI=" << numGrids[0] << ", J=" << numGrids[1] << ", K=" << numGrids[2] << ", ZONETYPE=Ordered\nDATAPACKING=POINT\n";

	std::ofstream f2(file_path + "_nonzerro");
	f2 << "TITLE     = \"" + str + "\"\nVARIABLES = \"x [mm]\"\n\"y [mm]\"\n\"z [mm]\"\n\"dx [mm]\"\n\"dy [mm]\"\n\"dz [mm]\"\n\"particlesNum [units]\"\nZONE T=\"Rectangular zone\"\nI=" << nonZeroCount << ", J=1, K=1, ZONETYPE=Ordered\nDATAPACKING=POINT\n";

	std::ofstream f3(file_path + "_af");
	f3 << "x[pix]\ty[pix]\tz[pix]\tvx[pix/dt]\tvy[pix/dt]\tvz[pix/dt]\ts[UNKNOWN]\n";

	int i, status;

	for (int z = 0; z < numGrids[2]; z++) {
		for (int y = 0; y < numGrids[1]; y++) {
			for (int x = 0; x < numGrids[0]; x++) {
				i = x * numGrids[2] * numGrids[1] + y * numGrids[2] + z;
				f << gridPoints[0][i] << "\t" << gridPoints[1][i] << "\t" << gridPoints[2][i] << "\t" << field[0][i] << "\t" << field[1][i] << "\t" << field[2][i] << "\t" << field[3][i] << "\n";
//Save to file
				//Position p = KrigingInterpolation(Position(gridPoints[0][i], gridPoints[1][i] + 0.625, gridPoints[2][i] + 0.625));
				//f << gridPoints[0][i] << "\t" << gridPoints[1][i] + 0.625 << "\t" << gridPoints[2][i] + 0.625 << "\t" << p.X() << "\t" << p.Y() << "\t" << p.Z() << "\n";

				status = 3;

				if ( field[0][i] == 0 && field[1][i] == 0 && field[2][i] == 0 ) status = 2;

				// this works if the field has been already filtered
				if ( true == m_bFiltered && 6 == field[3][i] ) status = 6;

				if ( 3 == status || 6 == status ) f2 << gridPoints[0][i] << "\t" << gridPoints[1][i] << "\t" << gridPoints[2][i] << "\t" << field[0][i] << "\t" << field[1][i] << "\t" << field[2][i] << "\t" << field[3][i] << "\n";

				f3 << gridPoints[0][i] << "\t" << gridPoints[1][i] << "\t" << gridPoints[2][i] << "\t" << field[0][i] << "\t" << field[1][i] << "\t" << field[2][i] << "\t" << status << "\n";
			}
		}
	}
}

bool is_inside(double pos, double center, double boders) {
	if ((pos - center) > boders || (center - pos) > boders)
		return false;
	return true;
}

int PredictiveField::ParticleFiltration() {

	int n = 5, m = 5, p = 5;
	int nd2 = (n - 1) / 2;
	int md2 = (m - 1) / 2;
	int pd2 = (p - 1) / 2;

	int N = 3;
	double acpt_fact = 0.5;

	int counter = 0, ind;

	for (int q = 0; q < N; q++) {

		double k = 0;
		double cur_k = 0;

		for (int l = 0; l < numGrids[2]; l++) {
			for (int j = 0; j < numGrids[1]; j++) {
				for (int i = 0; i < numGrids[0]; i++) {

					// calc proper mean
					double x = 0, y = 0, z = 0;
					counter = 0;

					for (int c = -pd2; c <= pd2; c++) {
						for (int b = -md2; b <= md2; b++) {
							for (int a = -nd2; a <= nd2; a++) {
								//���� ������ ������� �� ������� ������� - �� �������������
								if (l + c < 0 || l + c >= numGrids[2] || j + b < 0 || j + b >= numGrids[1] || i + a < 0 || i + a >= numGrids[0])
									continue;

								ind = (i + a) * numGrids[2] * numGrids[1] + (j + b) * numGrids[2] + (l + c);

								if (field[0][ind] == 0 && field[1][ind] == 0 && field[2][ind] == 0) continue;

								x += field[0][ind];
								y += field[1][ind];
								z += field[2][ind];

								counter++;
							}
						}
					}

					ind = i * numGrids[2] * numGrids[1] + j * numGrids[2] + l;
					if (field[0][ind] == 0 && field[1][ind] == 0 && field[2][ind] == 0) {
						field[3][ind] = 2;
					}

					if (counter == 0) continue;
					Position avrU(x, y, z);
					avrU /= counter;
					
					Position dU(field[0][ind], field[1][ind], field[2][ind]);
					dU -= avrU;

					//std::cout << dU.X() << " " << dU.Y() << " " << dU.Z() << endl;
					if ((cur_k = dU.Magnitude()) > k) k = cur_k;

				}
			}

		}

		k *= acpt_fact;

		//std::cout << k << endl;
		for (int l = 0; l < numGrids[2]; l++) {
			for (int j = 0; j < numGrids[1]; j++) {
				for (int i = 0; i < numGrids[0]; i++) {

					// calc proper mean
					double x = 0, y = 0, z = 0;
					counter = 0;

					for (int c = -pd2; c <= pd2; c++) {
						for (int b = -md2; b <= md2; b++) {
							for (int a = -nd2; a <= nd2; a++) {
								//���� ������ ������� �� ������� ������� - �� �������������
								if (l + c < 0 || l + c >= numGrids[2] || j + b < 0 || j + b >= numGrids[1] || i + a < 0 || i + a >= numGrids[0])
									continue;

								ind = (i + a) * numGrids[2] * numGrids[1] + (j + b) * numGrids[2] + (l + c);

								if (field[0][ind] == 0 && field[1][ind] == 0 && field[2][ind] == 0)
									continue;

								x += field[0][ind];
								y += field[1][ind];
								z += field[2][ind];

								counter++;
							}
						}
					}

					if (counter == 0) continue;
					Position avrU(x, y, z);
					avrU /= counter;

					ind = i * numGrids[2] * numGrids[1] + j * numGrids[2] + l;
					Position dU(field[0][ind], field[1][ind], field[2][ind]);

					dU -= avrU;

					if (dU.Magnitude() > k || field[3][ind] == 2) {
						field[0][ind] = avrU.X();
						field[1][ind] = avrU.Y();
						field[2][ind] = avrU.Z();
						field[3][ind] = 6;
					}

				}
			}
		}

	}

	m_bFiltered = true;

	return 0;
}
