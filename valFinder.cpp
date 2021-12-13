#include "headMat.h"

using namespace std;



int uptriangleMat(double(*mat), double(*tam), double(*m), double(*t), int n, int I, int J, double eps)
{
	double x = mat[I * n + I];
	double y = mat[J * n + I];
	if ((abs(x) > 1e-12 * eps) || (abs(y) > 1e-12 * eps))
	{
		double Cos = x / sqrt(x * x + y * y);
		double Sin = -y / sqrt(x * x + y * y);
		for (int j = I; j < n; j++)
		{

			m[j] = Cos * mat[I * n + j] - Sin * mat[J * n + j];
			m[n + j] = Sin * mat[I * n + j] + Cos * mat[J * n + j];
		}
		for (int j = I; j < n; j++)
		{
			mat[I * n + j] = m[j];
			mat[J * n + j] = m[n + j];
		}
		for (int j = 0; j < n; j++)
		{
			t[j] = Cos * tam[I * n + j] - Sin * tam[J * n + j];
			t[n + j] = Sin * tam[I * n + j] + Cos * tam[J * n + j];
		}
		for (int j = 0; j < n; j++)
		{
			tam[I * n + j] = t[j];
			tam[J * n + j] = t[n + j];
		}

		return 0;
	}

	return 1;
}

double normMat(double(*mat), double(*tam), int n)
{
	double a = 0, A = 0;
	for (int j = 0; j < n; j++)
	{
		a += abs(mat[j] - tam[j]);
	}
	A = a;
	for (int i = 1; i < n; i++)
	{
		a = 0;
		for (int j = 0; j < n; j++)
		{
			a += abs(mat[i * n + j] - tam[i * n + j]);
		}
		if (a > A)
		{
			A = a;
		}
	}
	return A;
}
/*
int rotMat(double(*mat), double(*T), int n, int I, int J) //����� ������ ������ �����
{
	double x = mat[J * n + J];
	double y = mat[I * n + J];

	idMat(T, n);
	if ((abs(x)  < 1e-16) || (abs(y)  < 1e-16))
	{
		double Cos = x / sqrt(x * x + y * y);
		double Sin = -y / sqrt(x * x + y * y);
		T[J * n + J] = Cos;
		T[J * n + I] = -Sin;
		T[I * n + J] = Sin;
		T[I * n + I] = Cos;
		return 0;
	}
	return 1;
}
*/
int threeMat(double(*mat), int n)
{
	double cos, sin;
	for (int i = 1; i < n - 1; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			double x = mat[i * n + i - 1];
			double y = mat[j * n + i - 1];
			double sq = sqrt(x * x + y * y);
			if ((abs(y) < 1e-15)|| (sq < 1e-15))
			{
				continue;
			}
			cos = x / sq;
			sin = -y / sq;
			mat[i * n + i - 1] = sq; 
			mat[(i - 1) * n + i] = sq;
			mat[j * n + i - 1] = 0; 
			mat[(i - 1) * n + j] = 0;
			for (int k = i + 1; k < n; k++)
			{
				if (k == j)
				{
					continue;
				}
				x = mat[i * n + k];
				y = mat[j * n + k];
				mat[k * n + i] = x * cos - y * sin;
				mat[i * n + k] = x * cos - y * sin;
				mat[k * n + j] = x * sin + y * cos;
				mat[j * n + k] = x * sin + y * cos;
			}
			double a_ii = mat[i * n + i] * cos - mat[j * n + i] * sin;
			double a_ji = mat[i * n + i] * sin + mat[j * n + i] * cos;
			double a_ij = mat[i * n + j] * cos - mat[j * n + j] * sin;
			double a_jj = mat[i * n + j] * sin + mat[j * n + j] * cos;
			mat[i * n + i] = a_ii * cos - a_ij * sin;
			mat[j * n + i] = a_ii * sin + a_ij * cos;
			mat[i * n + j] = mat[j * n + i];
			mat[j * n + j] = a_ji * sin + a_jj * cos;
		}
	}
	return 0;
}

void matToVecs(double(*mat), double(*a), double(*c), double(*d), int n)
{
	for (int i = 0; i < n-1; i++)
	{
		a[i] = mat[i * n + i];
		c[i] = mat[i * n + i + 1];
		d[i] = mat[(i + 1) * n + i];
	}
	a[n - 1] = mat[(n - 1) * n + n - 1];
}
int signChanges(double(*l), double(*a), double(*c), double(*d), int n, double lambd, double eps)
{
	int N=0;
	double L;
	//double L1;
	int i;
	for (i = 0; i < n; i++)
	{
		l[i] = a[i] - lambd;
	}
	i = 0;
	if (abs(l[i]) < eps*1e-15)
	{
		return -1000000;
	}
	L = l[i];
	//L1 = l[i];
	if (l[i] < 0)
	{
		N++;
	}
	for (i = 1; i < n; i++)
	{
		if (abs(l[i] - d[i-1] * c[i-1] / L) < eps*1e-15)
		{
			return -1000000;	
		}
		L = (l[i] - d[i-1] * c[i-1] / L);
		if (L < 0)
		{
			N++;
		}
	}
	return N;
}

void eigenBisection(double(*lambdas), double(*l), double(*a), double(*c), double(*d), int n, double rbord, double eps)
{
	rbord = rbord + 1;
	double RBORD = rbord;
	double lbord = -rbord;
	double centre = 0;
	int N = 0;
	int sgn = 0;
	int sgnr = 0;
	int sgnl = 0;
	int K = 0;
	while (N+K < n)
	{
		while (rbord - lbord > eps)
		{
			centre = (lbord + rbord) / 2;
			sgn = signChanges(l, a, c, d, n, centre, eps);
			if (sgn == -1000000)
			{
				//cout << "warp" << endl;
				sgnr = signChanges(l, a, c, d, n, centre + 3 * eps, eps);
				sgnl = signChanges(l, a, c, d, n, centre - 3 * eps, eps);
				//int rsg = signChanges(l, a, c, d, n, RBORD, eps);
				//subEigenBisection(lambdas, l, a, c, d, n, lbord, centre-3*eps, sgnl, 0,eps);
				//cout << "warp 2" << endl;
				//subEigenBisection(lambdas, l, a, c, d, n, centre + 3 * eps, RBORD , rsg-sgnr, sgnr, eps);
				//cnt = true;
				
				
				//cout << centre <<  ' ' << K << ' ' << N << endl;
				int current_N = sgnr - sgnl;
				for (int j = 0; j < current_N; j++)
				{
				//	cout << centre << ' ' << K << ' ' << N << ' ' << n - K - current_N << endl;
					lambdas[n-K-1-j] = centre;
				}
				K += current_N;
				if (sgnl < N + 1)
				{
					lbord = centre- 3*eps;
				}
				else
				{
					rbord = centre + 3*eps;
				}
				continue;
			}
			if (sgn < N + 1)
			{
				//cout << sgn << ' ' << N << " l" << endl;
				lbord = centre;
			}
			else
			{
				//cout << sgn << ' ' << N << " r" << endl;
				rbord = centre;
			}
		}
		centre = (lbord + rbord) / 2;

		int current_N = signChanges(l, a, c, d, n, rbord, eps) - signChanges(l, a, c, d, n, lbord, eps);
		//cout << lbord << ' ' << rbord << ' ' << current_N << endl;
		//if ((N + current_N) / n > itr)
		//	{
		//		cout << N << ' ' << current_N << ' ' << (N + current_N) / n << ' ';
		//		itr += 0.1;
		//	}
		//if (current_N != 0)
		//{
			for (int j = 0; j < current_N; j++)
			{
				//cout << "norm";
				lambdas[N + j] = centre;
			}
			N += current_N;
			//cout << lbord << ' ' << rbord << endl;
			lbord = centre;
			rbord = RBORD;
			///cout << lbord << ' ' << rbord << endl;
	//	}
		//else
		//{
		//	lbord = sub_centre;
			//rbord = RBORD;
		//}
	}
}
//void subEigenBisection(double(*lambdas), double(*l), double(*a), double(*c), double(*d), int n, double lbord, double rbord, int numEig, int L, double eps)
//{
//	rbord = rbord + 1.5*eps;
//	double RBORD = rbord;
//	int N = 0;
//	int sgn = 0;
//	bool cnt = false;
//	int sgnr = 0;
//	int sgnl = 0;
//	int K = 0;
//	double centre;
//	while (N < numEig)
//	{
//		double sub_centre = (lbord + rbord) / 2;
//		while (rbord - lbord > eps)
//		{
//			centre = (lbord + rbord) / 2;
//			sgn = signChanges(l, a, c, d, n, centre, eps)-L;
//			if (sgn == -1000000)
//			{
//				sgnr = signChanges(l, a, c, d, n, centre + 3 * eps, eps)-L;
//				sgnl = signChanges(l, a, c, d, n, centre - 3 * eps, eps)-L;
//				int current_N = sgnr - sgnl;
//				for (int j = 0; j < current_N; j++)
//				{
//					//	cout << centre << ' ' << K << ' ' << N << ' ' << n - K - current_N << endl;
//					//lambdas[n - K - 1 - j] = centre;
//					lambdas[N + L + j] = centre;
//				}
//				int rsg = signChanges(l, a, c, d, n, RBORD, eps);
//				int lsg = signChanges(l, a, c, d, n, lbord, eps);
//				subEigenBisection(lambdas, l, a, c, d, n, lbord, centre - 3 * eps, sgnl+L-lsg, lsg, eps);
//				subEigenBisection(lambdas, l, a, c, d, n, centre + 3 * eps, RBORD, rsg-sgnr+L, sgnr+L, eps);
//				cnt = true;
//				//cout << centre <<  ' ' << K << ' ' << N << endl;
//
//	
//				//K += current_N;
//				//if (sgnl < N + 1)
//				//{
//				//	lbord = centre- 3*eps;
//				//}
//				//else
//				//{
//				//	rbord = centre + 3*eps;
//				//}
//				break;
//			}
//			if (sgn < N + 1)
//			{
//				lbord = centre;
//			}
//			else
//			{
//				rbord = centre;
//			}
//		}
//		if (cnt)
//		{
//			break;
//		}
//		centre = (lbord + rbord) / 2;
//		//cout << lbord << ' ' << rbord << endl;
//		int current_N = signChanges(l, a, c, d, n, centre + 5 * eps, eps) - signChanges(l, a, c, d, n, centre - 5 * eps, eps);
//		cout << lbord << ' ' << rbord << ' ' << current_N << endl;
//		//if ((N + current_N) / n > itr)
//		//	{
//		//		cout << N << ' ' << current_N << ' ' << (N + current_N) / n << ' ';
//		//		itr += 0.1;
//		//	}
//		//if (current_N != 0)
//		//{
//		for (int j = 0; j < current_N; j++)
//		{
//			//cout << "norm";
//			lambdas[N +L+ j] = centre;
//		}
//		N += current_N;
//		cout << lbord << ' ' << rbord << endl;
//		lbord = centre;
//		rbord = RBORD;
//		cout << lbord << ' ' << rbord << endl;
//		//	}
//			//else
//			//{
//			//	lbord = sub_centre;
//				//rbord = RBORD;
//			//}
//	}
//}

double normByMaxMat(double(*mat), int n)
{
	double max = 0;
	double a;
	for (int i = 0; i < n; i++)
	{
		a = 0;
		for (int j = 0; j < n; j++)
		{
			a = a + abs(mat[i * n + j]);

		}
		if (a > max)
		{
			max = a;
		}
	}
	return max;
}

double firstInvariant(double(*mat), double(*lambdas), int n)
{
	double s=0;
	for (int i = 0; i < n; i++)
	{
		s = s+mat[i * n + i] - lambdas[i];
	}
	return abs(s);
}
double secondInvariant(double(*mat), double(*lambdas), int n)
{
	double s1 = 0;
	double s2 = 0;
	for (int i = 0; i < n; i++)
	{
		s1 = s1 + lambdas[i]* lambdas[i];
	}
	s1 = sqrt(s1);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			s2 = s2 + mat[i * n + j]* mat[i * n + j];
		}
	}
	s2 = sqrt(s2);
	return abs(s2 - s1);
}