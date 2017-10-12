#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>

using namespace std;


#define LEFT_UP 4
#define LEFT 1
#define UP 5


int LCS_LENGTH(const char X[], size_t len_X,const char Y[], size_t len_Y)
{
	size_t m = len_X;
	size_t n = len_Y;

	if (m == 0 || n == 0) return 0;

	int** b = new int*[m + 1];
	int**c = new int*[m + 1];

	for (size_t i = 0; i <= m; ++i)
	{
		c[i] = new int[n + 1];
		b[i] = new int[n + 1];
		c[i][0] = 0;
		//memset(c[i], 0, n + 1);
		//memset(b[i], 0, n + 1);
	}
	for (size_t j = 0; j <= n; ++j)
		c[0][j] = 0;

	for (size_t i = 1; i <= m; ++i) {
		for (size_t j = 1; j <= n; ++j) {
			if (X[i - 1] == Y[j - 1]) {
				c[i][j] = c[i - 1][j - 1] + 1;
				b[i][j] = LEFT_UP;
			}
			else if (c[i - 1][j] >= c[i][j - 1]) {
				c[i][j] = c[i - 1][j];
				b[i][j] = UP;
			}
			else {
				c[i][j] = c[i][j - 1];
				b[i][j] = LEFT;
			}
		}
	}

	int lcs_len = c[len_X][len_Y];


	for (size_t i = 0; i <= m; ++i)
	{
		delete[] c[i];
		delete[] b[i];
	}

	delete[] b, c;

	return lcs_len;
}



int main() {
	//freopen("C:/Users/suban/Documents/in.txt", "r", stdin);
	//freopen("C:/Users/suban/Documents/out.txt", "w", stdout);

	string X, Y;
	while (getline(cin,X) && getline(cin,Y))
	{
		size_t len_X = X.length();
		size_t len_Y = Y.length();

		int lcs_len = LCS_LENGTH(X.data(), len_X, Y.data(), len_Y);

		cout << lcs_len << endl;
	}
	//cout << X << endl;
	//cout << Y << endl;


	return 0;
}

