#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>

using namespace std;


#define LEFT_UP 4
#define LEFT 1
#define UP 5

template<typename T>
T ** new2D(size_t max_r, size_t max_c, T initVal)
{
	T** mat = new T*[max_r];
	for (size_t i = 0; i < max_r; i++)
	{
		mat[i] = new T[max_c];
		memset(mat[i], initVal, max_c * sizeof(T));
	}
	return mat;
}

template<typename T>
void delete2D(T ** mat, size_t max_r)
{
	for (size_t i = 0; i < max_r; i++)
	{
		delete[] mat[i];
	}
	delete[] mat;
}

int No_Of_Ways(vector<int>& value, int V) {
	int** ways = new2D<int>(value.size(), V + 1, 0);

	for (int i = 0; i < value.size(); i++)
	{
		ways[i][0] = 1;
	}
	for (int v = 1; v <= V; v++)
	{
		ways[0][v] = 1;
	}

	for (int i = 1; i < value.size(); i++)
	{
		for (int v = 1; v <= V; v++)
		{
			ways[i][v] = ways[i - 1][v] + (v >= value[i] ? ways[i][v - value[i]] : 0);
			//cout << value[i] << " " << v << " --- " << ways[i][v] << endl;
		}
	}

	int n_ways = ways[value.size() - 1][V];

	delete2D<int>(ways, value.size());

	return n_ways;
}

int main() {
	/*freopen("C:/Users/suban/Documents/in.txt", "r", stdin);
	freopen("C:/Users/suban/Documents/out.txt", "w", stdout);
	*/
	vector<int> vals = { 1,5,10,25,50 };


	int V = 0;
	while (cin >> V)
	{
		cout << No_Of_Ways(vals, V) << endl;
	}


	return 0;
}

