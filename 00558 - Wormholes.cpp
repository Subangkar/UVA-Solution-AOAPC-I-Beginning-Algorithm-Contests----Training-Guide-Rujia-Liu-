#define ACTIVE
#define _CRT_SECURE_NO_WARNINGS
#ifdef ACTIVE


#include<iostream>
#include<algorithm>
#include<vector>
#include<list>
#include<queue>
#include<cstring>
#include<cstdio>


using namespace std;


#define NULL_VALUE -999999
#define INFINITY 999999

#define NIL -1









//****************Heap start************************
#pragma region MINHEAP



#define MAX_HEAP_SIZE 100000+1

template <typename T> class HeapItem
{
public:
	T data;
	double key;
};


template <typename T> class MinHeap
{

public:
	MinHeap();
	~MinHeap();

	T MIN();
	T Extract_MIN();




	void Insert(T data, double key);

	bool IncreaseKey(int i, double key);
	bool DecreaseKey(int i, double key);

	size_t size();

	void BuildHeap(T Array[], double key[], int size);

	bool hasItem(T data);
	int getPos(T data);

	bool isEmpty();

	double getKey(T data);



	void PrintHeap()
	{
		for (int i = 1; i <= heap_size; i++)
		{
			cout << "---" << i << " === " << map[A[i].data] << "---" << A[i].data << "---" << A[i].key << endl;
		}
	}
protected:
	HeapItem<T>* A = NULL;
	int *map = NULL;
	size_t heap_size;

	void heapify_UpBottom(int i);
	void heapify_BottomUp(int i);

	int left(int i);
	int right(int i);
	int parent(int i);

	void swap(HeapItem<T>& x, HeapItem<T>& y);
	//void swap(int& x, int& y);

	bool smaller(HeapItem<T>& x, HeapItem<T>& y);
};


template <typename T>
MinHeap<T>::MinHeap()
{
	A = new HeapItem<T>[MAX_HEAP_SIZE];
	map = new int[MAX_HEAP_SIZE];
	memset(map, NULL, MAX_HEAP_SIZE * sizeof(int));
	heap_size = 0;
}


template <typename T>
MinHeap<T>::~MinHeap()
{
	if (A) delete[] A;
	if (map) delete[] map;
	heap_size = 0;
	A = NULL;
	map = NULL;
}



template <typename T>
void MinHeap<T>::Insert(T data, double key)
{
	heap_size = heap_size + 1;
	A[heap_size].data = data;
	A[heap_size].key = key;
	map[data] = heap_size;
	heapify_BottomUp(heap_size);
}

template<typename T>
T MinHeap<T>::Extract_MIN()
{
	if (!heap_size)
		return T(NULL);
	else {
		//cout << "Before Xtract ::: " << endl;
		//PrintHeap();
		T x = A[1].data;
		//map[A[1].data] = NULL;
		swap(A[1], A[heap_size]);
		heap_size = heap_size - 1;
		heapify_UpBottom(1);
		map[x] = NULL;
		//cout << "After Xtract ::: " << endl;
		//PrintHeap();
		return x;
	}
}

template<typename T>
T MinHeap<T>::MIN()
{
	if (isEmpty()) return T(NULL);
	else return A[1].data;
}


template<typename T>
void MinHeap<T>::heapify_BottomUp(int i)
{
	if (i > heap_size)	return;

	while (i > 1 && smaller(A[i], A[parent(i)]))
	{
		swap(A[i], A[parent(i)]);
		i = parent(i);
	}
}

template<typename T>
void MinHeap<T>::heapify_UpBottom(int i)
{
	if (i<1 || i > heap_size) return;

	int smallest = 0;

	int l = left(i), r = right(i);

	if (l <= heap_size && smaller(A[l], A[i]))
		smallest = l;
	else
		smallest = i;

	if (r <= heap_size && smaller(A[r], A[smallest]))
		smallest = r;

	if (smallest != i) {
		swap(A[i], A[smallest]);
		heapify_UpBottom(smallest);
	}
}


template<typename T>
bool MinHeap<T>::IncreaseKey(int i, double key)
{
	HeapItem<T> updated;
	updated.data = A[i].data;
	updated.key = key;
	if (smaller(updated, A[i]))
		return false;
	else
	{
		A[i].key = key;
		heapify_UpBottom(i);
		return true;
	}
}

template<typename T>
bool MinHeap<T>::DecreaseKey(int i, double key)
{
	HeapItem<T> updated;
	updated.data = A[i].data;
	updated.key = key;
	if (!smaller(updated, A[i]))
		return false;
	else
	{
		A[i].key = key;
		heapify_BottomUp(i);
		return true;
	}
}















template<typename T>
size_t MinHeap<T>::size()
{
	return heap_size;
}

template<typename T>
void MinHeap<T>::BuildHeap(T Array[], double key[], int size)
{
	heap_size = size;
	for (size_t i = 1; i <= heap_size; i++)
	{
		A[i].data = Array[i];
		A[i].key = key[i];

		map[A[i].data] = i;
	}

	for (size_t i = heap_size / 2; i > 0; i--)
	{
		heapify_UpBottom(i);
	}
}

template<typename T>
bool MinHeap<T>::hasItem(T data)
{
	return map[data];
}

template<typename T>
int MinHeap<T>::getPos(T data)
{
	return map[data];
}

template<typename T>
bool MinHeap<T>::isEmpty()
{
	return !heap_size;
}

template<typename T>
double MinHeap<T>::getKey(T data)
{
	//cout << "----"<<data<<"----" << A[map[data]].key << endl;
	if (!hasItem(data)) return NULL;
	return A[map[data]].key;
}

template<typename T>
int MinHeap<T>::left(int i)
{
	return 2 * i;
}

template<typename T>
int MinHeap<T>::right(int i)
{
	return 2 * i + 1;
}

template<typename T>
int MinHeap<T>::parent(int i)
{
	return i / 2;
}

template<typename T>
void MinHeap<T>::swap(HeapItem<T> & x, HeapItem<T> & y)
{
	T tempData = x.data;
	x.data = y.data;
	y.data = tempData;


	//swap(map[x.data], map[y.data]);
	int t = map[x.data];
	map[x.data] = map[y.data];
	map[y.data] = t;



	double tmpKey = x.key;
	x.key = y.key;
	y.key = tmpKey;
}

//template<typename T>
//void MinHeap<T>::swap(int & x, int & y)
//{
//	int t = x;
//	x = y;
//	y = t;
//}






template<typename T>
bool MinHeap<T>::smaller(HeapItem<T>& x, HeapItem<T>& y)
{
	if (x.key < y.key) return true;
	else if (x.key == y.key && x.data < y.data) return true;
	return false;
}


#pragma endregion
//****************Heap end************************






class Edge
{
public:
	int v = -1;
	int w = 0;
	bool operator==(Edge edg2)
	{
		return v == edg2.v;
	}

	Edge(int v = -1) { this->v = v; }
	Edge(int v, int w) { this->v = v; this->w = w; }
};

















#pragma region GRAPH


#define WHITE 1 // NOT VISITED
#define GREY 2 // DISCOVERED
#define BLACK 3 // VISITED



#define EDGE_Y 1
#define EDGE_N NIL


#define ADJ_LIST_TYPE vector<Edge>
#define ADJLISTPOS vector<Edge>::iterator 

//******************Graph class starts here**************************
class Graph
{
	int nVertices, nEdges;
	bool directed;
	bool isMatrix;
	bool startAt1;
	int W_total;

	int ** matrix; //adjacency matrix to store the graph
				   //define other variables required for bfs such as color, parent, and dist
				   //you must use pointers and dynamic allocation
	ADJ_LIST_TYPE* adjList;


	int *color;
	int *dist;
	int *parent;

	int *time_discovery;
	int *time_finishing;
	int time;


	double *key;
	int minCost;

	ADJ_LIST_TYPE::iterator findVInListOfU(int u, int v) {
		return find(adjList[u].begin(), adjList[u].end(), Edge(v));
	}

	bool isNULL(ADJ_LIST_TYPE::iterator it, int u) {
		return it == adjList[u].end();
	}


	template<typename T>
	T** new2D(size_t max_r, size_t max_c, T initVal = NULL);

	template<typename T>
	void delete2D(T** mat, size_t max_r);

	void relax(int u,int v,int w) {
		if (dist[v] > dist[u] + w) {
			dist[v] = dist[u] + w;
			parent[v] = u;
		}	
	}
public:
	Graph(bool dir = false, bool isMatrix = false, bool startAt1 = true);
	Graph(const Graph& g);
	~Graph();
	void setnVertices(int n);
	void addEdge(int u, int v, int w);
	void removeEdge(int u, int v);
	bool isEdge(int u, int v);
	bool hasVertex(int u) { return !(!nVertices || u < 0 || u > nVertices); }
	int getWeight(int u, int v);
	void setWeight(int u, int v, int w);


	int getDegree(int u);
	bool hasCommonAdjacent(int u, int v);
	bool isConnected(int source);


	bool ShortestPath_BellmanFord(int src);
	void printShortestPathInfo(int source);

	void printEdges();
	void printGraph();

};

template<typename T>
T ** Graph::new2D(size_t max_r, size_t max_c, T initVal)
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
void Graph::delete2D(T ** mat, size_t max_r)
{
	for (size_t i = 0; i < max_r; i++)
	{
		delete[] mat[i];
	}
	delete[] mat;
}

int Graph::getWeight(int u, int v)
{
	if (hasVertex(u) && hasVertex(v) && isEdge(u, v)) {
		if (isMatrix)
			return matrix[u][v];
		else {
			ADJLISTPOS i = findVInListOfU(u, v);
			if (!isNULL(i, u)) return i->w;
		}
	}
	return 0;
}

Graph::Graph(bool dir, bool isMatrix, bool start1)
{
	nVertices = 0;
	nEdges = 0;
	W_total = 0;

	matrix = NULL;
	adjList = NULL;

	directed = dir; //set direction of the graph
					//define other variables to be initialized
	startAt1 = start1;
	this->isMatrix = isMatrix;

	color = NULL;
	dist = NULL;
	parent = NULL;

	time_discovery = NULL;
	time_finishing = NULL;

	time = 0;

	minCost = 0;
	key = NULL;
}

Graph::Graph(const Graph & g)
{
	isMatrix = g.isMatrix;
	directed = g.directed;
	startAt1 = g.startAt1;

	setnVertices(g.nVertices);

	memcpy(color, g.color, nVertices + 1);
	memcpy(dist, g.dist, nVertices + 1);
	memcpy(parent, g.parent, nVertices + 1);

	memcpy(time_discovery, g.time_discovery, nVertices + 1);
	memcpy(time_finishing, g.time_finishing, nVertices + 1);

	memcpy(key, g.key, nVertices + 1);


}

Graph::~Graph()
{
	if (nVertices) {
		if (isMatrix) {
			delete2D<int>(matrix, nVertices + 1);
		}
		else
		{
			if (adjList) delete[] adjList; //delete previous list
		}

		delete[] color, dist, parent;

		delete[] time_discovery, time_finishing;

		delete[] key;
	}

	matrix = NULL;
	adjList = NULL;

	color = NULL;
	dist = NULL;
	parent = NULL;

	time_discovery = NULL;
	time_finishing = NULL;
	time = 0;


	nVertices = 0;
	nEdges = 0;
	W_total = 0;

	key = NULL;
}



// Set Vertices Numbered from 0 to n
void Graph::setnVertices(int n)
{
	if (nVertices) this->~Graph();
	nVertices = n + 1;


	if (isMatrix) {
		//allocate space for the matrix
		matrix = new2D<int>(nVertices, nVertices, EDGE_N);
	}
	else
	{
		adjList = new ADJ_LIST_TYPE[nVertices];
	}


	color = new int[nVertices];
	dist = new int[nVertices];
	parent = new int[nVertices];

	time_discovery = new int[nVertices];
	time_finishing = new int[nVertices];


	key = new double[nVertices];

	nVertices = n;
	W_total = 0;
}

// O(1)
void Graph::addEdge(int u, int v, int w)
{
	if (!hasVertex(u) || !hasVertex(v)) return;
	nEdges++;


	if (isMatrix) {
		matrix[u][v] = w;
		if (!directed) matrix[v][u] = w;
	}
	else
	{
		adjList[u].push_back(Edge(v, w));
		if (!directed) adjList[v].push_back(Edge(u, w));
	}

	W_total += w;
}

// MAT-O(1) & LIST-O(deg(V))
void Graph::removeEdge(int u, int v)
{
	if (!hasVertex(u) || !hasVertex(v) || !isEdge(u, v)) return; //vertex out of range
	--nEdges;
	W_total -= getWeight(u, v);

	if (isMatrix) {
		matrix[u][v] = EDGE_N;
		if (!directed) matrix[v][u] = EDGE_N;
	}
	else
	{
		ADJLISTPOS i = findVInListOfU(u, v);
		if (isNULL(i, u)) return;
		adjList[u].erase(i);
		if (!directed) adjList[v].erase(findVInListOfU(v, u));
	}
}

//returns true if (u,v) is an edge, otherwise should return false
// MAT-O(1) & LIST-O(deg(V))
bool Graph::isEdge(int u, int v)
{
	if (isMatrix)
		return hasVertex(u) && hasVertex(v) && matrix[u][v] != EDGE_N;
	else
		return hasVertex(u) && hasVertex(v) && !isNULL(findVInListOfU(u, v), u);
}




// MAT-O(deg(V)) & LIST-O(1)
int Graph::getDegree(int u)
{
	//returns the degree of vertex u
	if (!hasVertex(u)) return NIL;

	size_t noOfDeg = 0;
	if (isMatrix)for (int v = 0; v <= nVertices; ++v) if (isEdge(u, v)) ++noOfDeg;
	else noOfDeg = adjList[u].size();
	return noOfDeg;
}

bool Graph::hasCommonAdjacent(int u, int v)
{
	if (isMatrix) {
		//returns true if vertices u and v have common adjacent vertices
		for (int i = 0; i <= nVertices; ++i) {
			if (isEdge(u, i) && isEdge(v, i)) return true;
		}
	}
	else
	{
		for (size_t iu = 0; iu < adjList[u].size(); ++iu)
			for (size_t iv = 0; iv < adjList[v].size(); ++iv)
				if (*findVInListOfU(v, iv) == *findVInListOfU(u, iu)) return true;
	}
	return false;
}




void Graph::printGraph()
{
	printf("\nNumber of vertices: %d, Number of edges: %d\n", nVertices, nEdges);

	if (isMatrix) {
		for (int i = 0; i <= nVertices; i++)
		{
			cout << i << ":";
			for (int j = 0; j <= nVertices; j++)
			{
				if (isEdge(i, j))
					cout << " " << j << " ( " << getWeight(i, j) << " )";
			}
			cout << endl;
		}
	}

	else {
		for (int i = 0; i <= nVertices; i++)
		{
			cout << i << ":";
			for (int j = 0; j < adjList[i].size(); j++)
			{
				cout << " " << adjList[i][j].v << " ( " << getWeight(i, adjList[i][j].v) << " )";
			}
			cout << endl;
		}
	}



}

bool Graph::ShortestPath_BellmanFord(int src)
{
	for (int u = 0; u <= nVertices; u++)
	{
		dist[u] = INFINITY;
		parent[u] = NIL;
	}

	dist[src] = 0;
	parent[src] = NIL;

	for (int i = 1; i < nVertices; i++)
	{
		for (int u = 1; u <= nVertices; u++)
		{
			if (isMatrix) {
				for (int v = 1; v <= nVertices; v++)
				{
					if (isEdge(u, v)) relax(u, v, getWeight(u,v));
				}
			}
			else
			{
				int v = 0;
				for (int c = 0; c < adjList[u].size() ; c++)
				{
					v = adjList[u][c].v;
					relax(u, v, getWeight(u, v));
				}
			}
		}
	}


	for (int u = 1; u <= nVertices; u++)
	{
		if (isMatrix) {
			for (int v = 1; v <= nVertices; v++)
			{
				if (isEdge(u, v) && dist[v] > dist[u]+getWeight(u,v)) return false;
			}
		}
		else
		{
			int v = 0;
			for (int c = 0; c < adjList[u].size(); c++)
			{
				v = adjList[u][c].v;
				if (dist[v] > dist[u] + getWeight(u, v)) return false;
			}
		}
	}


	return true;
}

void Graph::printShortestPathInfo(int source)
{
	bool hasNegEdgeCycle = !ShortestPath_BellmanFord(source);

	if (!hasNegEdgeCycle) {
		cout << "Source : " << source << endl;
		//cout << "v <> v.d <> v.p" << endl;
		printf("%3s <> %3s <> %3s\n", "v", "v.d", "v.p");
		for (int u = 1; u <= nVertices; u++)
		{
			//cout << u << " <> " << dist[u] << " <> " << parent[u] << endl;
			printf("%3d <> %3d <> %3d\n", u, dist[u], parent[u]);
		}
	}

	cout << "Negative Cycle? " << (hasNegEdgeCycle ? "YES" : "NO") << endl;
}



void Graph::printEdges()
{
	bool **printed = new2D<bool>(nVertices + 1, nVertices + 1, false);

	for (size_t u = 0; u <= nVertices; u++)
	{
		for (int v = 0; v <= nVertices; ++v)
		{
			if (!printed[u][v] && isEdge(u, v))
			{
				cout << u << " " << v << endl;
				printed[u][v] = true;
				if (!directed) printed[v][u] = true;
			}
		}
	}

	delete2D<bool>(printed, nVertices + 1);
}

void Graph::setWeight(int u, int v, int w)
{
	W_total -= getWeight(u, v);

	if (isMatrix) {
		matrix[u][v] = w;
		if (!directed) matrix[v][u] = w;
	}
	else
	{
		ADJLISTPOS i = findVInListOfU(u, v);
		if (isNULL(i, u)) {
			addEdge(u, v, w);
		}
		else
		{
			i->w = w;
			if (!directed) {
				findVInListOfU(v, u)->w = w;
			}
		}
	}

	W_total += getWeight(u, v);
}


#pragma endregion




int main()
{
	int t = 0, cs = 0;

	cin >> t;


	while (cs<t)
	{
		cs++;
		bool startFrom1 = true;
		int source = 1;
		bool matrix = false;
		bool directed = true;


		Graph g(directed, matrix, startFrom1);

		int n, m;
		n = m = 0;

		int u, v, w;
		u = v = w = 0;

		cin >> n >> m;

		g.setnVertices(n);

		for (size_t i = 0; i < m; i++)
		{
			cin >> u >> v >> w;
			g.addEdge(u + 1, v + 1, w);
		}
		int src = 1;

		if (!g.ShortestPath_BellmanFord(src)) cout << "possible" << endl;
		else cout << "not possible" << endl;

	}
	return 0;
}



#endif // ACTIVE

