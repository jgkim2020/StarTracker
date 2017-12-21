#ifndef csv_h
#define csv_h

#include <vector>
using namespace std;

template <class T>
class csv
{
public:
	csv<T>(const char* filename)
		: filename(filename)
	{}
	void read();
	void write(vector<string>& wheader, vector<vector<T>>& wfields);
	vector<string> getheader();
	vector<vector<T>> getfields();
	void setquietMode(bool arg);
	
private:
	const char* filename;
	vector<string> header;
	vector<vector<T>> fields;
	bool quietMode = true;
};

#endif