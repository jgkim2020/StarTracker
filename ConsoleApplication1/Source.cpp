#include <iostream>
#include <vector>
#include <chrono>
#include "csv.h"
#include "csv.cpp"

using namespace std;
using namespace chrono;

void findStar(const char* filename, int32_t sigma, int32_t clusterMinSize, int32_t clusterMaxSize, bool saveCSV, bool debugMode);
void getAttitude();

int main(void)
{
	//findStar("11 c 1.8 100.bmp", 5, 14, 250, false, true);
	findStar("test2.bmp", 1, 5, 40000, true, true);
	getAttitude();
	return 0;
}

void getAttitude()
{
	// 5. 
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
	int32_t row = *(int32_t*)(header + 22);
	int32_t column = *(int32_t*)(header + 18);
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
	vector<float> centroid_x;
	vector<float> centroid_y;
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
	}

	// return
	fclose(in);
	delete[] bgr;
	delete[] data;
	delete[] nonZeroCluster;
	delete[] nonZeroClusterSize;
}