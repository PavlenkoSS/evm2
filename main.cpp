#include "headMat.h"
using namespace std;


int main(int nargs, char** args)
{
	if (nargs < 5)
	{
		cout << "Not enought args: launch as following..." << endl;
		return -1;
	}
	int n = 0;
	int m = 0;
	int k = 1;
	double eps = 1;
	char* filename = new char[128];
	//string filn = "matfil.txt";
	if (sscanf(args[1], "%d", &n) != 1 || sscanf(args[2], "%d", &m) != 1 || sscanf(args[3], "%lf", &eps) != 1 ||sscanf(args[4], "%d", &k) != 1)
	{
		cout << "bad params" << endl;
		return -1;
	}
	if (n == 0)
	{
		cout << "bad param n" << endl;
		return -1;
	}
	if (k == 0 && (args[5] == nullptr|| sscanf(args[5], "%s", filename) != 1))
	{
		cout << "No filename" << endl;
		return -1;
	}
	double* M;
	double* a;
	double* c;
	double* d;
	double* lambdas;
	double* l;
	M = new double[n * n];
	a = new double[n];
	c = new double[n - 1];
	d = new double[n - 1];
	lambdas = new double[n];
	l = new double[n];
	if (fulMat(M, n, k, filename) == -1)
	{
		delete[]M;
		delete[]a;
		delete[]c;
		delete[]d;
		delete[]lambdas;
		delete[]l;
		return -1;
	}
	cout << "M = " << endl;
	outMat1(M, n, m);

	double nrm = normByMaxMat(M, n); // перенести в один метод 
	//
	//threeMat(M, n);
	//outMat1(M, n, m);

	// threeMat(M,n);
	// matToVecs(M, a, c, d, n); //
	// eigenBisection(lambdas, l, a, c, d, n, nrm, eps); //
	double start_time = clock();
	superFun(M,n, a,c,d, lambdas,l,nrm,eps);
	double end_time = clock();

	fulMat(M, n, k, filename);
	cout << setprecision(3);
	cout << "First invariant " << firstInvariant(M, lambdas, n) << endl;
	cout << "Second invariant " << secondInvariant(M, lambdas, n) << endl;
	outMat1(lambdas, 1, n, m);
	
	cout << "Time of solving = " << (-start_time + end_time) / CLOCKS_PER_SEC << endl;
	delete[]M;
	delete[]a;
	delete[]c;
	delete[]d;
	delete[]lambdas;
	delete[]l;
	return 0;
}
