#include <iostream>
#include <vector>
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
	//findStar("test2.bmp", 1, 5, 40000, true, true);

	getAttitude(15, 16.27, 1.5e-5, 4.96309e-7, 0.942007877, false, false);

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
			pairEuclidean.push_back(sqrtf(deltaX*deltaX + deltaY*deltaY));
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
			//
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
	steady_clock::time_point end = steady_clock::now();
	cout << "elapsed time: " << duration_cast<milliseconds>(end - begin).count() << "ms" << endl;

	// 6. catalog matching (k-vector search)
	cout << "6. catalog matching" << endl;
	begin = steady_clock::now();
	double f = column/5.7*fmm; // focal length in pixels
	double pairSpherical[6];
	int32_t lowerIndex[6];
	int32_t upperIndex[6];
	int32_t count = 0;
	vector<int32_t> pyramid_myID = pyramid_sorted[0];
	// find index range
	for (int32_t i = 0; i < 3; i++)
	{
		for (int32_t j = i + 1; j < 4; j++)
		{
			double x1 = centroid_x[pyramid_myID[i]];
			double y1 = centroid_y[pyramid_myID[i]];
			double x2 = centroid_x[pyramid_myID[j]];
			double y2 = centroid_y[pyramid_myID[j]];
			pairSpherical[count] = (x1*x2 + y1*y2 + f*f)/(sqrt(x1*x1 + y1*y1 + f*f)*sqrt(x2*x2 + y2*y2 + f*f));
			upperIndex[count] = kvector[(int32_t)((pairSpherical[count] + err - yintercept)/slope) - 1][2]; // upper bound (excluded)
			lowerIndex[count] = kvector[(int32_t)((pairSpherical[count] - err - yintercept)/slope) + 1][2]; // lower bound (included)
			count++;
		}
	}
	// find pyramid candidates
	bool pairID[2647][6] = {false, };
	int32_t pairIDcount[2647] = {0, };
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
				if (count > 1) cout << i << endl;
			}
			pyramid_starID.push_back(entry);
		}
	}
	end = steady_clock::now();

	for (int32_t i = 0; i < pyramid_starID.size(); i++)
	{
		for (int32_t j = 0; j < 4; j++) cout << pyramid_starID[i][j] << " ";
		cout << endl;
	}

	cout << "elapsed time: " << duration_cast<milliseconds>(end - begin).count() << "ms" << endl;

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
	double sfx2 = atan(centroid_x[pyramid_myID[1]] / f);
	double sfy2 = atan(centroid_y[pyramid_myID[1]] / f*cos(sfx2));
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
	sfn[1][0] = pixelCoord1[1] * pixelCoord2[2] - pixelCoord2[1] * pixelCoord1[2];
	sfn[1][1] = pixelCoord2[0] * pixelCoord1[2] - pixelCoord1[0] * pixelCoord2[2];
	sfn[1][2] = pixelCoord1[0] * pixelCoord2[1] - pixelCoord2[0] * pixelCoord1[1];
	sfn[2][0] = -(pixelCoord1[1] * sfn[1][2] - sfn[1][1] * pixelCoord1[2]);
	sfn[2][1] = -(sfn[1][0] * pixelCoord1[2] - pixelCoord1[0] * sfn[1][2]);
	sfn[2][2] = -(pixelCoord1[0] * sfn[1][1] - sfn[1][0] * pixelCoord1[1]);
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

	for (int32_t i = 0; i < 3; i++)
	{
		for (int32_t j = 0; j < 3; j++) cout << rotationMatrix[i][j] << " ";
		cout << endl;
	}

	double eye[3][3] = {0, };
	for (int32_t i = 0; i < 3; i++)
	{
		for (int32_t j = 0; j < 3; j++)
		{
			for (int32_t k = 0; k < 3; k++) eye[i][j] += rotationMatrix[i][k]*rotationMatrix[j][k];
		}
	}

	for (int32_t i = 0; i < 3; i++)
	{
		for (int32_t j = 0; j < 3; j++) cout << eye[i][j] << " ";
		cout << endl;
	}

	end = steady_clock::now();
	cout << "elapsed time: " << duration_cast<milliseconds>(end - begin).count() << "ms" << endl;
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
		csv<int32_t> csvout("original_image.csv");
		vector<string> header;
		vector<vector<int32_t>> fields;
		for (int32_t i = 0; i < column; i++) header.push_back("C" + to_string(i + 1));
		for (int32_t i = 0; i < row; i++)
		{
			fields.push_back(vector<int32_t>());
			for (int32_t j = 0; j < column; j++) fields.back().push_back((int32_t)data[i*column + j]);
		}
		csvout.write(header, fields);
		cout << "saved original_image.csv" << endl;
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
		csv<int32_t> csvout("threshold_image.csv");
		vector<string> header;
		vector<vector<int32_t>> fields;
		for (int32_t i = 0; i < column; i++) header.push_back("C" + to_string(i + 1));
		for (int32_t i = 0; i < row; i++)
		{
			fields.push_back(vector<int32_t>());
			for (int32_t j = 0; j < column; j++) fields.back().push_back((int32_t)data[i*column + j]);
		}
		csvout.write(header, fields);
		cout << "saved threshold_image.csv" << endl;
	}

	// 3. find and label connected components (crude implementation)
	cout << "3. clustering & labelling" << endl;
	begin = steady_clock::now();
	int32_t* nonZeroCluster = new int32_t[size];
	int32_t* nonZeroClusterSize = new int32_t[size];
	int32_t numberofCluster = 0;
	int32_t connectivity8[8] = { -column - 1, -column, -column + 1, -1, 1, column - 1, column, column + 1 };
	// initialize array
	for (int32_t i = 0; i < size; i++)
	{
		nonZeroCluster[i] = -1;
		nonZeroClusterSize[i] = 0;
	}
	// find non-zero connected component (of size > 1) and assign a parent label
	for (int32_t i = 0; i < size; i++)
	{
		if (data[i] > 0) // for non-zero pixels
		{
			int32_t minimumNeighborID = size;
			int32_t myID = i;
			// find minimum neighbor ID (excluding self)
			for (int32_t j = 0; j < 8; j++)
			{
				// for admissible indices
				if (-1 < i + connectivity8[j] && i + connectivity8[j] < size)
				{
					// for ID assigned non-zero pixels
					if (data[i + connectivity8[j]] > 0 && nonZeroCluster[i + connectivity8[j]] > -1)
					{
						int32_t rootID = nonZeroCluster[i + connectivity8[j]];
						while (rootID != nonZeroCluster[rootID]) rootID = nonZeroCluster[rootID]; // find root
						if (minimumNeighborID > rootID) minimumNeighborID = rootID;
					}
				}
			}
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
							nonZeroCluster[i + connectivity8[j]] = rootID; // update neighbor
						}
					}
				}
			}
		}
	}
	// post processing
	for (int32_t i = 0; i < size; i++)
	{
		if (data[i] > 0 && nonZeroCluster[i] != i) // for non-zero, non-root pixels
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
	if (saveCSV)
	{
		csv<int32_t> csvout("non-zero_pixel_map.csv");
		vector<string> header;
		vector<vector<int32_t>> fields;
		for (int32_t i = 0; i < column; i++) header.push_back("C" + to_string(i + 1));
		for (int32_t i = 0; i < row; i++)
		{
			fields.push_back(vector<int32_t>());
			for (int32_t j = 0; j < column; j++) fields.back().push_back(nonZeroCluster[i*column + j]);
		}
		csvout.write(header, fields);
		cout << "saved non-zero_pixel_map.csv" << endl;
	}
	if (saveCSV)
	{
		csv<int32_t> csvout("non-zero_cluster_size.csv");
		vector<string> header;
		vector<vector<int32_t>> fields;
		for (int32_t i = 0; i < column; i++) header.push_back("C" + to_string(i + 1));
		for (int32_t i = 0; i < row; i++)
		{
			fields.push_back(vector<int32_t>());
			for (int32_t j = 0; j < column; j++) fields.back().push_back(nonZeroClusterSize[i*column + j]);
		}
		csvout.write(header, fields);
		cout << "saved non-zero_cluster_size.csv" << endl;
	}
		
	// 4. find centroid
	cout << "4. get cendtroid" << endl;
	begin = steady_clock::now();
	for (int32_t i = 0; i < size; i++)
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
			centroid_x.push_back(((float)sumx) / sum + 0.5f - 0.5f*column);
			centroid_y.push_back(0.5f*row - ((float)sumy) / sum - 0.5f);
		}
	}
	end = steady_clock::now();
	cout << "elapsed time: " << duration_cast<milliseconds>(end - begin).count() << "ms" << endl;
	if (saveCSV)
	{
		csv<float> csvout("centroid.csv");
		vector<string> header;
		vector<vector<float>> fields;
		header.push_back("X"); header.push_back("Y");
		for (uint32_t i = 0; i < centroid_x.size(); i++)
		{
			fields.push_back(vector<float>());
			fields.back().push_back(centroid_x[i]); fields.back().push_back(centroid_y[i]);
		}
		csvout.write(header, fields);
		cout << "saved centroid.csv" << endl;
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