#ifndef csv_cpp
#define csv_cpp

#include "csv.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

template <class T>
void csv<T>::read()
{
	ifstream in(filename);
	if (in)
	{
		string line;
		bool isheader = true;
		while (getline(in, line))
		{
			stringstream sep(line);
			string field;
			if (!isheader)
			{
				fields.push_back(vector<T>());
				while (getline(sep, field, ','))
				{
					fields.back().push_back((T)stod(field));
				}
			}
			else
			{
				while (getline(sep, field, ','))
				{
					header.push_back(field);
				}
				isheader = false;
			}
		}
	}
}

template <class T>
void csv<T>::write(vector<string>& wheader, vector<vector<T>>& wfields)
{
	ofstream out(filename);
	for (auto it = wheader.begin(); it != wheader.end(); it++)
	{
		out << *it << ",";
	}
	out << endl;
	for (auto it = wfields.begin(); it != wfields.end(); it++)
	{
		for (auto it2 = it->begin(); it2 != it->end(); it2++)
		{
			out << *it2 << ",";
		}
		out << endl;
	}
}

template <class T>
vector<string> csv<T>::getheader()
{
	if (!quietMode)
	{
		for (auto it = header.begin(); it != header.end(); ++it)
		{
			cout << *it;
		}
	}
	return header;
}

template <class T>
vector<vector<T>> csv<T>::getfields()
{
	if (!quietMode)
	{
		for (auto row : fields)
		{
			for (auto field : row)
			{
				cout << field << ' ';
			}
			cout << '\n';
		}
	}
	return fields;
}

template <class T>
void csv<T>::setquietMode(bool arg)
{
	quietMode = arg;
}

#endif