#include <iostream>
#include <vector>
#include <map>
#include <chrono>
#include "csv.h"
#include "csv.cpp"
#include "catalog.h"

using namespace std;
using namespace chrono;

int32_t row = 0;
int32_t column = 0;
vector<float> centroid_x;
vector<float> centroid_y;

void findStar(const char* filename, int32_t sigma, int32_t clusterMinSize, int32_t clusterMaxSize, bool saveCSV, bool debugMode);
void getAttitude(int32_t pairMaxNumber, double fmm, double err, double slope, double yintercept, bool saveCSV, bool debugMode);

int main(void)
{
	//findStar("11 vf 1.8 100.bmp", 4, 14, 250, false, true);
	findStar("11 c 1.8 100.bmp", 5, 14, 250, false, true);
	//findStar("test.bmp", 1, 5, 40000, true, true);

	getAttitude(20, 16.27, 1.5e-5, 4.96309e-7, 0.942007877, true, false);

	return 0;
}

void getAttitude(int32_t pairMaxNumber, double fmm, double err, double slope, double yintercept, bool saveCSV, bool debugMode)
{
	// 5. selection
	cout << "5. selection" << endl;
	steady_clock::time_point begin = steady_clock::now();
	vector<int32_t> pairEuclidean;
	vector<int32_t> pairIndex1;
	vector<int32_t> pairIndex2;
	int32_t numberofCluster = centroid_x.size();
	int32_t index = 0;
	for (int32_t i = 0; i < numberofCluster - 1; i++)
	{
		for (int32_t j = i + 1; j < numberofCluster; j++)
		{
			float deltaX = centroid_x[i] - centroid_x[j];
			float deltaY = centroid_y[i] - centroid_y[j];
			pairEuclidean.push_back((int32_t)sqrtf(deltaX*deltaX + deltaY*deltaY));
			pairIndex1.push_back(i);
			pairIndex2.push_back(j);
		}
	}
	// sort in ascending order (bubble sort)
	int32_t numberofPair = pairEuclidean.size();
	for (int32_t i = 0; i < pairMaxNumber && i < numberofPair; i++)
	{
		for (int32_t j = numberofPair - 1; j > i; j--)
		{
			if (pairEuclidean[j] < pairEuclidean[j - 1])
			{
				int32_t swapEuclidean = pairEuclidean[j - 1];
				int32_t swapIndex1 = pairIndex1[j - 1];
				int32_t swapIndex2 = pairIndex2[j - 1];
				pairEuclidean[j - 1] = pairEuclidean[j];
				pairIndex1[j - 1] = pairIndex1[j];
				pairIndex2[j - 1] = pairIndex2[j];
				pairEuclidean[j] = swapEuclidean;
				pairIndex1[j] = swapIndex1;
				pairIndex2[j] = swapIndex2;
			}
		}
	}
	// find all possible pyramids
	vector<int32_t> frequency(numberofCluster, 0);
	vector<vector<int32_t>> pyramid;
	for (int32_t i = 0; i < pairMaxNumber && i < numberofPair; i++)
	{
		frequency[pairIndex1[i]]++;
		frequency[pairIndex2[i]]++;
	}

	for (int32_t i = 0; i < numberofCluster; i++)
	{
		if (frequency[i] > 2) // 3 pairs (or more) with a commmon element i forms a pyramid (or more)
		{
			vector<int32_t> pyramidElement;
			for (int32_t j = 0; j < pairMaxNumber && j < numberofPair; j++)
			{
				if (i == pairIndex1[j]) pyramidElement.push_back(pairIndex2[j]);
				else if (i == pairIndex2[j]) pyramidElement.push_back(pairIndex1[j]);
			}
			for (int32_t j = 0; j < frequency[i] - 2; j++)
			{
				for (int32_t k = j + 1; k < frequency[i] - 1; k++)
				{
					for (int32_t l = k + 1; l < frequency[i]; l++)
					{
						vector<int32_t> entry;
						entry.push_back(i); entry.push_back(pyramidElement[j]);
						entry.push_back(pyramidElement[k]); entry.push_back(pyramidElement[l]);
						for (int32_t l1 = 0; l1 < 4; l1++) // entry element in ascending order
						{
							for (int32_t l2 = 3; l2 > l1; l2--)
							{
								if (entry[l2] < entry[l2 - 1])
								{
									int32_t swap = entry[l2 - 1];
									entry[l2 - 1] = entry[l2];
									entry[l2] = swap;
								}
							}
						}
						pyramid.push_back(entry);
					}
				}
			}
		}
	}
	// post processing
	vector<vector<int32_t>> pyramid_sorted;
	for (int32_t i = 0; i < numberofCluster - 3; i++)
	{
		for (int32_t j = i + 1; j < numberofCluster - 2; j++)
		{
			for (int32_t k = j + 1; k < numberofCluster - 1; k++)
			{
				for (int32_t l = k + 1; l < numberofCluster; l++)
				{
					int32_t m = 0;
					bool foundMatch = false;
					while (m < pyramid.size() && !foundMatch)
					{
						if (pyramid[m][0] == i && pyramid[m][1] == j && pyramid[m][2] == k && pyramid[m][3] == l)
						{
							pyramid_sorted.push_back(vector<int32_t>());
							pyramid_sorted.back().push_back(i);
							pyramid_sorted.back().push_back(j);
							pyramid_sorted.back().push_back(k);
							pyramid_sorted.back().push_back(l);
							foundMatch = true;
						}
						m++;
					}
				}
			}
		}
	}
	// find the best pyramid (smallest pyramid)
	int32_t pyramidEuclidMin = row*column;
	int32_t pyramidEuclidMinIndex = -1;
	for (int32_t i = 0; i < pyramid_sorted.size(); i++)
	{
		int32_t pyramidEuclidSum = 0;
		for (int32_t j = 0; j < 3; j++)
		{
			for (int32_t k = j + 1; k < 4; k++)
			{
				int32_t myID1 = pyramid_sorted[i][j];
				int32_t myID2 = pyramid_sorted[i][k];
				float deltaX = centroid_x[myID1] - centroid_x[myID2];
				float deltaY = centroid_y[myID1] - centroid_y[myID2];
				pyramidEuclidSum += (int32_t)sqrtf(deltaX*deltaX + deltaY*deltaY);
			}
		}
		if (pyramidEuclidSum < pyramidEuclidMin)
		{
			pyramidEuclidMin = pyramidEuclidSum;
			pyramidEuclidMinIndex = i;
		}
		pyramid_sorted[i].push_back(pyramidEuclidSum);
	}
	steady_clock::time_point end = steady_clock::now();
	cout << "elapsed time: " << duration_cast<milliseconds>(end - begin).count() << "ms" << endl;
	if (saveCSV)
	{
		cout << "saving results in csv format..." << endl;
		csv<int32_t> csvout("5pyramid_candidates_myID.csv");
		vector<string> header;
		for (int32_t i = 0; i < 4; i++) header.push_back("myID " + to_string(i + 1));
		header.push_back("pyramid Euclid sum");
		csvout.write(header, pyramid_sorted);
		cout << "saved 5pyramid_candidates_myID.csv" << endl;
	}

	// 6. catalog matching (k-vector search)
	cout << "6. catalog matching" << endl;
	begin = steady_clock::now();
	double f = column / 5.7*fmm; // focal length in pixels
	double pairSpherical[6];
	int32_t lowerIndex[6];
	int32_t upperIndex[6];
	int32_t count = 0;
	vector<int32_t> pyramid_myID = pyramid_sorted[pyramidEuclidMinIndex];
	// find index range
	for (int32_t i = 0; i < 3; i++)
	{
		for (int32_t j = i + 1; j < 4; j++)
		{
			double x1 = centroid_x[pyramid_myID[i]];
			double y1 = centroid_y[pyramid_myID[i]];
			double x2 = centroid_x[pyramid_myID[j]];
			double y2 = centroid_y[pyramid_myID[j]];
			pairSpherical[count] = (x1*x2 + y1*y2 + f*f) / (sqrt(x1*x1 + y1*y1 + f*f)*sqrt(x2*x2 + y2*y2 + f*f));
			upperIndex[count] = kvector[(int32_t)((pairSpherical[count] + err - yintercept) / slope) - 1][2]; // upper bound (excluded)
			lowerIndex[count] = kvector[(int32_t)((pairSpherical[count] - err - yintercept) / slope) + 1][2]; // lower bound (included)
			count++;
		}
	}
	// find pyramid candidates
	bool pairID[2647][6] = { false, };
	int32_t pairIDcount[2647] = { 0, };
	for (int32_t i = 0; i < 6; i++)
	{
		for (int32_t j = lowerIndex[i]; j < upperIndex[i]; j++)
		{
			if (!pairID[kvector[j][0]][i])
			{
				pairID[kvector[j][0]][i] = true;
				pairIDcount[kvector[j][0]]++;
			}
			if (!pairID[kvector[j][1]][i])
			{
				pairID[kvector[j][1]][i] = true;
				pairIDcount[kvector[j][1]]++;
			}
		}
	}
	int32_t pairPattern[4][3] =
	{
		0, 1, 2,
		0, 3, 4,
		1, 3, 5,
		2, 4, 5
	};
	vector<vector<int32_t>> pyramid_starID;
	for (int32_t i = 0; i < 2647; i++)
	{
		int32_t cond = -1;
		if (pairID[i][0] && pairID[i][1] && pairID[i][2] && !pairID[i][3] && !pairID[i][4] && !pairID[i][5]) cond = 0;
		else if (pairID[i][0] && !pairID[i][1] && !pairID[i][2] && pairID[i][3] && pairID[i][4] && !pairID[i][5]) cond = 1;
		else if (!pairID[i][0] && pairID[i][1] && !pairID[i][2] && pairID[i][3] && !pairID[i][4] && pairID[i][5]) cond = 2;
		else if (!pairID[i][0] && !pairID[i][1] && pairID[i][2] && !pairID[i][3] && pairID[i][4] && pairID[i][5]) cond = 3;
		if (!(cond < 0))
		{
			vector<int32_t> entry(4);
			entry[cond] = i; // common ID
			bool existDuplicate = false;
			for (int32_t j = 0; j < 3; j++)
			{
				int32_t count = 0;
				for (int32_t k = lowerIndex[pairPattern[cond][j]]; k < upperIndex[pairPattern[cond][j]]; k++)
				{
					if (kvector[k][0] == i)
					{
						if (cond > j) entry[j] = kvector[k][1];
						else entry[j + 1] = kvector[k][1];
						count++;
					}
					else if (kvector[k][1] == i)
					{
						if (cond > j) entry[j] = kvector[k][0];
						else entry[j + 1] = kvector[k][0];
						count++;
					}
				}
				if (count > 1) existDuplicate = true;
			}
			if (!existDuplicate) pyramid_starID.push_back(entry);
		}
	}
	end = steady_clock::now();
	cout << "elapsed time: " << duration_cast<milliseconds>(end - begin).count() << "ms" << endl;
	if (saveCSV)
	{
		cout << "saving results in csv format..." << endl;
		csv<int32_t> csvout("6pyramid_starID.csv");
		vector<string> header;
		for (int32_t i = 0; i < 4; i++) header.push_back("starID " + to_string(i + 1));
		csvout.write(header, pyramid_starID);
		cout << "saved 6pyramid_starID.csv" << endl;
	}
	if (saveCSV)
	{
		cout << "saving results in csv format..." << endl;
		csv<string> csvout("6pair_candidates_starID.csv");
		vector<string> header;
		vector<vector<string>> fields;
		header.push_back("pyramid pair"); header.push_back("spherical distance");
		header.push_back("starID 1"); header.push_back("starID 2");
		for (int32_t i = 0; i < 6; i++)
		{
			for (int32_t j = lowerIndex[i]; j < upperIndex[i]; j++)
			{
				fields.push_back(vector<string>());
				fields.back().push_back(to_string(i + 1));
				fields.back().push_back(to_string(pairSpherical[i]));
				fields.back().push_back(to_string(kvector[j][0]));
				fields.back().push_back(to_string(kvector[j][1]));
			}
		}
		csvout.write(header, fields);
		cout << "saved 6pair_candidates_starID.csv" << endl;
	}

	// 7. TRIAD
	cout << "7. TRIAD" << endl;
	begin = steady_clock::now();
	// find skyCoord & pixelCoord
	double skyCoord1[3];
	double skyCoord2[3];
	double pixelCoord1[3];
	double pixelCoord2[3];
	for (int32_t i = 0; i < 3; i++)
	{
		skyCoord1[i] = mcf[pyramid_starID[0][0] - 1][i];
		skyCoord2[i] = mcf[pyramid_starID[0][1] - 1][i];
	}
	double sfx1 = atan(centroid_x[pyramid_myID[0]]/f);
	double sfy1 = atan(centroid_y[pyramid_myID[0]]/f*cos(sfx1));
	double sfx2 = atan(centroid_x[pyramid_myID[1]]/f);
	double sfy2 = atan(centroid_y[pyramid_myID[1]]/f*cos(sfx2));
	pixelCoord1[0] = -sin(sfx1)*cos(sfy1); pixelCoord1[1] = cos(sfx1)*cos(sfy1); pixelCoord1[2] = -sin(sfy1);
	pixelCoord2[0] = -sin(sfx2)*cos(sfy2); pixelCoord2[1] = cos(sfx2)*cos(sfy2); pixelCoord2[2] = -sin(sfy2);
	// find ifn & sfn
	double ifn[3][3];
	double sfn[3][3];
	for (int32_t i = 0; i < 3; i++) ifn[0][i] = skyCoord1[i];
	ifn[1][0] = skyCoord1[1]*skyCoord2[2] - skyCoord2[1]*skyCoord1[2];
	ifn[1][1] = skyCoord2[0]*skyCoord1[2] - skyCoord1[0]*skyCoord2[2];
	ifn[1][2] = skyCoord1[0]*skyCoord2[1] - skyCoord2[0]*skyCoord1[1];
	ifn[2][0] = -(skyCoord1[1]*ifn[1][2] - ifn[1][1]*skyCoord1[2]);
	ifn[2][1] = -(ifn[1][0]*skyCoord1[2] - skyCoord1[0]*ifn[1][2]);
	ifn[2][2] = -(skyCoord1[0]*ifn[1][1] - ifn[1][0]*skyCoord1[1]);
	for (int32_t i = 1; i < 3; i++) // normalize
	{
		double norm = 0;
		for (int32_t j = 0; j < 3; j++) norm += ifn[i][j] * ifn[i][j];
		norm = sqrt(norm);
		for (int32_t j = 0; j < 3; j++) ifn[i][j] /= norm;
	}
	for (int32_t i = 0; i < 3; i++) sfn[0][i] = pixelCoord1[i];
	sfn[1][0] = pixelCoord1[1]*pixelCoord2[2] - pixelCoord2[1]*pixelCoord1[2];
	sfn[1][1] = pixelCoord2[0]*pixelCoord1[2] - pixelCoord1[0]*pixelCoord2[2];
	sfn[1][2] = pixelCoord1[0]*pixelCoord2[1] - pixelCoord2[0]*pixelCoord1[1];
	sfn[2][0] = -(pixelCoord1[1]*sfn[1][2] - sfn[1][1]*pixelCoord1[2]);
	sfn[2][1] = -(sfn[1][0]*pixelCoord1[2] - pixelCoord1[0]*sfn[1][2]);
	sfn[2][2] = -(pixelCoord1[0]*sfn[1][1] - sfn[1][0]*pixelCoord1[1]);
	for (int32_t i = 1; i < 3; i++) // normalize
	{
		double norm = 0;
		for (int32_t j = 0; j < 3; j++) norm += sfn[i][j] * sfn[i][j];
		norm = sqrt(norm);
		for (int32_t j = 0; j < 3; j++) sfn[i][j] /= norm;
	}
	double rotationMatrix[3][3] = {0, };
	for (int32_t i = 0; i < 3; i++)
	{
		for (int32_t j = 0; j < 3; j++)
		{
			for (int32_t k = 0; k < 3; k++) rotationMatrix[i][j] += sfn[i][k]*ifn[j][k];
		}
	}
	end = steady_clock::now();
	cout << "elapsed time: " << duration_cast<milliseconds>(end - begin).count() << "ms" << endl;
	if (saveCSV)
	{
		cout << "saving results in csv format..." << endl;
		csv<double> csvout("7rotation_Matrix.csv");
		vector<string> header;
		vector<vector<double>> fields;
		for (int32_t i = 0; i < 3; i++) header.push_back("C" + to_string(i + 1));
		for (int32_t i = 0; i < 3; i++) fields.push_back(vector<double>(3));
		for (int32_t i = 0; i < 3; i++)
		{
			for (int32_t j = 0; j < 3; j++) fields[i][j] = rotationMatrix[i][j];
		}
		csvout.write(header, fields);
		cout << "saved 7rotation_Matrix.csv" << endl;
	}
}

void findStar(const char* filename, int32_t sigma, int32_t clusterMinSize, int32_t clusterMaxSize, bool saveCSV, bool debugMode)
{
	// 1. read image
	cout << "1. read image" << endl;
	steady_clock::time_point begin = steady_clock::now();
	FILE *in = fopen(filename, "rb");
	// read header & extract row/column size
	uint8_t header[54];
	fread(header, sizeof(uint8_t), 54, in);
	row = *(int32_t*)(header + 22);
	column = *(int32_t*)(header + 18);
	// read raw image 
	int32_t size = row*column;
	uint8_t* bgr = new uint8_t[3*size];
	uint8_t* data = new uint8_t[size];
	fread(bgr, sizeof(uint8_t), 3*size, in);
	for (int32_t i = 0; i < row; i++) // (image stored left to right from bottom row to top row)
	{
		for (int32_t j = 0; j < column; j++)
		{
			data[(row - 1 - i)*column + j] = (bgr[3*(i*column + j)] + bgr[3*(i*column + j) + 1] + bgr[3*(i*column + j) + 2])/3;
		}
	}
	steady_clock::time_point end = steady_clock::now();
	cout << "elapsed time: " << duration_cast<milliseconds>(end - begin).count() << "ms" << endl;
	if (saveCSV)
	{
		cout << "saving results in csv format..." << endl;
		csv<int32_t> csvout("1original_image.csv");
		vector<string> header;
		vector<vector<int32_t>> fields;
		for (int32_t i = 0; i < column; i++) header.push_back("C" + to_string(i + 1));
		for (int32_t i = 0; i < row; i++)
		{
			fields.push_back(vector<int32_t>());
			for (int32_t j = 0; j < column; j++) fields.back().push_back((int32_t)data[i*column + j]);
		}
		csvout.write(header, fields);
		cout << "saved 1original_image.csv" << endl;
	}
	
	// 2. thresholding (global)
	cout << "2. thresholding" << endl;
	begin = steady_clock::now();
	int32_t sum = 0;
	int64_t ssum = 0;
	float mean = 0;
	float stddev = 0;
	for (int32_t i = 0; i < size; i++)
	{
		sum += (int32_t)data[i];
		ssum += ((int64_t)data[i])*((int64_t)data[i]);
	}
	mean = (float)sum/size;
	stddev = sqrtf((float)ssum/size - mean*mean);
	for (int32_t i = 0; i < size; i++) if (data[i] < sigma*stddev) data[i] = 0;
	end = steady_clock::now();
	cout << "elapsed time: " << duration_cast<milliseconds>(end - begin).count() << "ms" << endl;
	if (saveCSV)
	{
		cout << "saving results in csv format..." << endl;
		csv<int32_t> csvout("2threshold_image.csv");
		vector<string> header;
		vector<vector<int32_t>> fields;
		for (int32_t i = 0; i < column; i++) header.push_back("C" + to_string(i + 1));
		for (int32_t i = 0; i < row; i++)
		{
			fields.push_back(vector<int32_t>());
			for (int32_t j = 0; j < column; j++) fields.back().push_back((int32_t)data[i*column + j]);
		}
		csvout.write(header, fields);
		cout << "saved 2threshold_image.csv" << endl;
	}

	// 3. find and label connected components
	int32_t* nonZeroCluster = new int32_t[size];
	int32_t* nonZeroClusterSize = new int32_t[size];
	int32_t numberofCluster = 0;
	map<int32_t, int32_t> nonZeroCluster_dict;
	map<int32_t, int32_t> nonZeroClusterSize_dict;
	int32_t numberofCluster_dict = 0;
	int32_t connectivity8[8] = { -column - 1, -column, -column + 1, -1, 1, column - 1, column, column + 1 };
	
	// crude implementation (with arrays)
	cout << "3. clustering & labelling (with arrays)" << endl;
	begin = steady_clock::now();
	// initialize array
	for (int32_t i = 0; i < size; i++)
	{
		nonZeroCluster[i] = -1;
		nonZeroClusterSize[i] = 0;
	}
	// find non-zero connected component (excluding single pixel clusters) and assign a parent label
	for (int32_t i = 0; i < size; i++)
	{
		if (data[i] > 0) // for non-zero pixels
		{
			int32_t minimumNeighborID = size;
			int32_t myID = i;
			bool isolated = true;
			// find minimum neighbor ID (excluding self)
			for (int32_t j = 0; j < 8; j++)
			{
				// for admissible indices
				if (-1 < i + connectivity8[j] && i + connectivity8[j] < size)
				{
					// for ID assigned non-zero pixels
					if (data[i + connectivity8[j]] > 0)
					{
						if (nonZeroCluster[i + connectivity8[j]] > -1)
						{
							int32_t rootID = nonZeroCluster[i + connectivity8[j]];
							while (rootID != nonZeroCluster[rootID]) rootID = nonZeroCluster[rootID]; // find root
							if (minimumNeighborID > rootID) minimumNeighborID = rootID;
						}
						isolated = false;
					}
				}
			}
			if (!isolated) // exclude single pixel clusters
			{
				// update self
				if (minimumNeighborID >= myID) minimumNeighborID = i;
				nonZeroCluster[i] = minimumNeighborID;
				nonZeroClusterSize[minimumNeighborID]++;
				// update neighbors
				for (int32_t j = 0; j < 8; j++)
				{
					if (-1 < i + connectivity8[j] && i + connectivity8[j] < size) // for admissible indices
					{
						if (data[i + connectivity8[j]] > 0) // for non-zero pixels
						{
							if (nonZeroCluster[i + connectivity8[j]] < 0) // unassinged non-zero pixels
							{
								nonZeroCluster[i + connectivity8[j]] = minimumNeighborID;
							}
							else // ID assigned non-zero pixels
							{
								int32_t rootID = nonZeroCluster[i + connectivity8[j]];
								while (rootID != nonZeroCluster[rootID]) rootID = nonZeroCluster[rootID]; // find root
								if (rootID > minimumNeighborID) nonZeroCluster[rootID] = minimumNeighborID; // update root
								nonZeroCluster[i + connectivity8[j]] = minimumNeighborID; // update neighbor
							}
						}
					}
				}
			}
		}
	}
	// post processing
	for (int32_t i = 0; i < size; i++)
	{
		if (nonZeroCluster[i] != -1 && nonZeroCluster[i] != i) // for non-zero, non-root, non-single pixel cluster pixels
		{
			int32_t rootID = nonZeroCluster[i];
			while (rootID != nonZeroCluster[rootID]) rootID = nonZeroCluster[rootID]; // find root
			nonZeroCluster[i] = rootID; // update cluster ID
			if (nonZeroClusterSize[i] != 0)
			{
				nonZeroClusterSize[nonZeroCluster[i]] += nonZeroClusterSize[i]; // update cluster size
				nonZeroClusterSize[i] = 0;
			}
		}
		if (nonZeroCluster[i] == i) numberofCluster++;
	}
	end = steady_clock::now();
	cout << "elapsed time: " << duration_cast<milliseconds>(end - begin).count() << "ms" << endl;
	cout << "memory usage: " << 2*4*size/1000 << "Kbyte" << endl;
	if (saveCSV)
	{
		cout << "saving results in csv format..." << endl;
		csv<int32_t> csvout("3non-zero_pixel_map.csv");
		vector<string> header;
		vector<vector<int32_t>> fields;
		for (int32_t i = 0; i < column; i++) header.push_back("C" + to_string(i + 1));
		for (int32_t i = 0; i < row; i++)
		{
			fields.push_back(vector<int32_t>());
			for (int32_t j = 0; j < column; j++) fields.back().push_back(nonZeroCluster[i*column + j]);
		}
		csvout.write(header, fields);
		cout << "saved 3non-zero_pixel_map.csv" << endl;
	}
	if (saveCSV)
	{
		cout << "saving results in csv format..." << endl;
		csv<int32_t> csvout("3non-zero_cluster_size.csv");
		vector<string> header;
		vector<vector<int32_t>> fields;
		for (int32_t i = 0; i < column; i++) header.push_back("C" + to_string(i + 1));
		for (int32_t i = 0; i < row; i++)
		{
			fields.push_back(vector<int32_t>());
			for (int32_t j = 0; j < column; j++) fields.back().push_back(nonZeroClusterSize[i*column + j]);
		}
		csvout.write(header, fields);
		cout << "saved 3non-zero_cluster_size.csv" << endl;
	}

	// memory optimized implementation (with map STL)
	cout << "3. clustering & labelling (with map STL)" << endl;
	begin = steady_clock::now();
	// find non-zero connected component (excluding single pixel clusters) and assign a parent label
	for (int32_t i = 0; i < size; i++)
	{
		if (data[i] > 0) // for non-zero pixels
		{
			int32_t minimumNeighborID = size;
			int32_t myID = i;
			bool isolated = true;
			// find minimum neighbor ID (excluding self)
			for (int32_t j = 0; j < 8; j++)
			{
				// for admissible indices
				if (-1 < i + connectivity8[j] && i + connectivity8[j] < size)
				{
					/* (little bit slower)
					// for ID assigned non-zero pixels
					if (data[i + connectivity8[j]] > 0)
					{
						if (nonZeroCluster_dict.count(i + connectivity8[j]))
						{
							int32_t rootID = nonZeroCluster_dict[i + connectivity8[j]];
							while (rootID != nonZeroCluster_dict[rootID]) rootID = nonZeroCluster_dict[rootID]; // find root
							if (minimumNeighborID > rootID) minimumNeighborID = rootID;
						}
						isolated = false;
					}
					*/

					if (data[i + connectivity8[j]] > 0)
					{
						if (j < 4) // already registered in nonZeroCluster_dict
						{
							int32_t rootID = nonZeroCluster_dict[i + connectivity8[j]];
							int32_t rootID_ = rootID + 1;
							while (rootID != rootID_) // find root
							{
								rootID_ = rootID;
								rootID = nonZeroCluster_dict[rootID];
							}
							if (minimumNeighborID > rootID) minimumNeighborID = rootID;
						}
						else if (j < 7) // may or may not be registered in nonZeroCluster_dict
						{
							if (nonZeroCluster_dict.find(i + connectivity8[j]) != nonZeroCluster_dict.end())
							{
								int32_t rootID = nonZeroCluster_dict[i + connectivity8[j]];
								int32_t rootID_ = rootID + 1;
								while (rootID != rootID_) // find root
								{
									rootID_ = rootID;
									rootID = nonZeroCluster_dict[rootID];
								}
								if (minimumNeighborID > rootID) minimumNeighborID = rootID;
							}
						}
						else // not registered in nonZeroCluster_dict
						{
							// do nothing
						}
						isolated = false;
					}

				}
			}
			if (!isolated) // exclude single pixel clusters
			{
				// update self
				if (minimumNeighborID >= myID) minimumNeighborID = i;
				nonZeroCluster_dict[i] = minimumNeighborID;
				nonZeroClusterSize_dict[minimumNeighborID]++;
				// update neighbors
				for (int32_t j = 0; j < 8; j++)
				{
					if (-1 < i + connectivity8[j] && i + connectivity8[j] < size) // for admissible indices
					{
						if (data[i + connectivity8[j]] > 0) // for non-zero pixels
						{
							/* (little bit slower)
							if (nonZeroCluster_dict.count(i + connectivity8[j]) == 0) // unassinged non-zero pixels
							{
								nonZeroCluster_dict[i + connectivity8[j]] = minimumNeighborID;
							}
							else // ID assigned non-zero pixels
							{
								int32_t rootID = nonZeroCluster_dict[i + connectivity8[j]];
								int32_t rootID_ = rootID + 1;
								while (rootID != rootID_) // find root
								{
									rootID_ = rootID;
									rootID = nonZeroCluster_dict[rootID];
								}
								if (rootID > minimumNeighborID) nonZeroCluster_dict[rootID] = minimumNeighborID; // update root
								nonZeroCluster_dict[i + connectivity8[j]] = minimumNeighborID; // update neighbor
							}
							*/

							if (j < 4) // already registered in nonZeroCluster_dict
							{
								int32_t rootID = nonZeroCluster_dict[i + connectivity8[j]];
								int32_t rootID_ = rootID + 1;
								while (rootID != rootID_) // find root
								{
									rootID_ = rootID;
									rootID = nonZeroCluster_dict[rootID];
								}
								if (rootID > minimumNeighborID) nonZeroCluster_dict[rootID] = minimumNeighborID; // update root
								nonZeroCluster_dict[i + connectivity8[j]] = minimumNeighborID; // update neighbor
							}
							else if (j < 7) // may or may not be registered in nonZeroCluster_dict
							{
								if (nonZeroCluster_dict.find(i + connectivity8[j]) == nonZeroCluster_dict.end()) // unassinged non-zero pixels
								{
									nonZeroCluster_dict[i + connectivity8[j]] = minimumNeighborID;
								}
								else // ID assigned non-zero pixels
								{
									int32_t rootID = nonZeroCluster_dict[i + connectivity8[j]];
									int32_t rootID_ = rootID + 1;
									while (rootID != rootID_) // find root
									{
										rootID_ = rootID;
										rootID = nonZeroCluster_dict[rootID];
									}
									if (rootID > minimumNeighborID) nonZeroCluster_dict[rootID] = minimumNeighborID; // update root
									nonZeroCluster_dict[i + connectivity8[j]] = minimumNeighborID; // update neighbor
								}
							}
							else // not registered in nonZeroCluster_dict
							{
								nonZeroCluster_dict[i + connectivity8[j]] = minimumNeighborID;
							}

						}
					}
				}
			}
		}
	}
	// post processing
	for (auto it = nonZeroCluster_dict.begin(); it != nonZeroCluster_dict.end(); it++) // for non-zero, non-single pixel cluster pixels
	{
		if (it->first != it->second) // for non-root cluster pixels
		{
			int32_t rootID = it->second;
			while (rootID != nonZeroCluster_dict[rootID]) rootID = nonZeroCluster_dict[rootID]; // find root
			it->second = rootID; // update cluster ID
			if (nonZeroClusterSize_dict.find(it->first) != nonZeroClusterSize_dict.end())
			{
				nonZeroClusterSize_dict[it->second] += nonZeroClusterSize_dict[it->first]; // update cluster size
				nonZeroClusterSize_dict.erase(it->first);
			}
		}
		else numberofCluster_dict++; // for root cluster pixels
	}
	end = steady_clock::now();
	cout << "elapsed time: " << duration_cast<milliseconds>(end - begin).count() << "ms" << endl;
	cout << "memory usage: " << (nonZeroCluster_dict.size() + nonZeroClusterSize_dict.size())*8/1000 << "Kbyte" << endl;
	if (saveCSV)
	{
		cout << "saving results in csv format..." << endl;
		csv<int32_t> csvout("3non-zero_pixel_map_STL.csv");
		vector<string> header;
		vector<vector<int32_t>> fields;
		for (int32_t i = 0; i < column; i++) header.push_back("C" + to_string(i + 1));
		for (int32_t i = 0; i < row; i++) fields.push_back(vector<int32_t>(column, -1)); // initialize with -1
		for (auto it = nonZeroCluster_dict.begin(); it != nonZeroCluster_dict.end(); it++)
		{
			fields[(it->first)/column][(it->first)%column] = it->second;
		}
		csvout.write(header, fields);
		cout << "saved 3non-zero_pixel_map_STL.csv" << endl;
	}
	if (saveCSV)
	{
		cout << "saving results in csv format..." << endl;
		csv<int32_t> csvout("3non-zero_cluster_size_STL.csv");
		vector<string> header;
		vector<vector<int32_t>> fields;
		for (int32_t i = 0; i < column; i++) header.push_back("C" + to_string(i + 1));
		for (int32_t i = 0; i < row; i++) fields.push_back(vector<int32_t>(column, 0)); // initialize with 0
		for (auto it = nonZeroClusterSize_dict.begin(); it != nonZeroClusterSize_dict.end(); it++)
		{
			fields[(it->first)/column][(it->first)%column] = it->second;
		}
		csvout.write(header, fields);
		cout << "saved 3non-zero_cluster_size_STL.csv" << endl;
	}
		
	// 4. find centroid
	cout << "4. get cendtroid" << endl;
	begin = steady_clock::now();
	vector<float> centroid_pixelx;
	vector<float> centroid_pixely;
	// array for clustering
	/*for (int32_t i = 0; i < size; i++)
	{
		if (nonZeroClusterSize[i] > clusterMinSize && nonZeroClusterSize[i] < clusterMaxSize)
		{
			int32_t window = i + connectivity8[7] + 1;
			int32_t sumx = 0;
			int32_t sumy = 0;
			int32_t sum = 0;
			for (int32_t j = i; j < window && j < size; j++)
			{
				if (nonZeroCluster[j] == i)
				{
					sumx += (j%column)*data[j];
					sumy += (j/column)*data[j];
					sum += data[j];
					window = j + connectivity8[7] + 1;
				}
			}
			centroid_x.push_back(((float)sumx)/sum + 0.5f - 0.5f*column);
			centroid_y.push_back(0.5f*row - ((float)sumy)/sum - 0.5f);
			centroid_pixelx.push_back(((float)sumx)/sum);
			centroid_pixely.push_back(((float)sumy)/sum);
		}
	}*/
	// map for clustering
	for (auto it = nonZeroClusterSize_dict.begin(); it != nonZeroClusterSize_dict.end(); it++)
	{
		if (it->second > clusterMinSize && it->second < clusterMaxSize)
		{
			int32_t window = it->first + connectivity8[7] + 1;
			int32_t sumx = 0;
			int32_t sumy = 0;
			int32_t sum = 0;
			int32_t count = 0;
			for (auto it2 = nonZeroCluster_dict.find(it->first); it2 != nonZeroCluster_dict.end() && count < it->second; it2++)
			{
				if (it2->second == it->first)
				{
					int32_t index = it2->first;
					sumx += (index%column)*data[index];
					sumy += (index/column)*data[index];
					sum += data[index];
					count++;
				}
			}
			centroid_x.push_back(((float)sumx)/sum + 0.5f - 0.5f*column);
			centroid_y.push_back(0.5f*row - ((float)sumy)/sum - 0.5f);
			centroid_pixelx.push_back(((float)sumx)/sum);
			centroid_pixely.push_back(((float)sumy)/sum);
		}
	}
	end = steady_clock::now();
	cout << "elapsed time: " << duration_cast<milliseconds>(end - begin).count() << "ms" << endl;
	if (saveCSV)
	{
		cout << "saving results in csv format..." << endl;
		csv<float> csvout("4centroid.csv");
		vector<string> header;
		vector<vector<float>> fields;
		header.push_back("myID"); header.push_back("X"); header.push_back("Y");
		header.push_back("X'"); header.push_back("Y'");
		for (uint32_t i = 0; i < centroid_x.size(); i++)
		{
			fields.push_back(vector<float>());
			fields.back().push_back(i);
			fields.back().push_back(centroid_x[i]); fields.back().push_back(centroid_y[i]);
			fields.back().push_back(centroid_pixelx[i]); fields.back().push_back(centroid_pixely[i]);
		}
		csvout.write(header, fields);
		cout << "saved 4centroid.csv" << endl;
	}

	// for debugging purposes
	if (debugMode)
	{
		cout << endl;
		cout << "width(columns): " << column << ", height(rows): " << row << endl;
		cout << "sum: " << sum << ", ssum: " << ssum << endl;
		cout << "mean: " << mean << ", stddev: " << stddev << endl;
		cout << "number of non-zero clusters: " << numberofCluster << endl;
		cout << "number of stars: " << centroid_x.size() << endl;
		cout << endl;
	}

	// return
	fclose(in);
	delete[] bgr;
	delete[] data;
	delete[] nonZeroCluster;
	delete[] nonZeroClusterSize;
}