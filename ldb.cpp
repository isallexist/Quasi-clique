#include<iostream>
#include<fstream>
#include<sstream>
#include<string.h>

#include<vector>
#include<set>

#include<algorithm>
#include<time.h>
#include<math.h>

#include<omp.h>

#include "mkl.h"

using namespace std;

typedef vector<string> strv;

typedef set<int> ints;			
typedef vector<ints> adj;

typedef vector<double> douv;
typedef vector<douv> douvv;

typedef vector<int> intv;
typedef vector<intv> intvv;

const double eps = 1e-8;
const double INF = 1e8;
const int MAT_FORMAT = LAPACK_ROW_MAJOR;
const double GAMMA[] = {0.95,0.85};

//-----------------------------------------------
// input & output, not changing often
const string INPUT = "sparse"; //"dimacs";	
const string OUTPUT = "output/results.csv";
//------------------------------------------------

// Marco ------------------------------------------------------------------------
//#define INWIN
#define HASOMP
//-------------------------------------------------------------------------------

void splitstr(strv & splits, const string& str, const string& deli = " ")
{	// Note: splits are never cleared inside!!!
	size_t start, end;
	
	start = str.find_first_not_of(deli);
	while(start!=string::npos)
	{
		end = str.find_first_of(deli,start+1);

		if(end!=string::npos)
		{
			splits.push_back(str.substr(start,end-start));
			start = str.find_first_not_of(deli,end+1);
		}
		else 
		{
			splits.push_back(str.substr(start));	
			break;
		}
	}
}

class graph
{	
private:
	int thres;
	
	void ParseDimacs(const string& struct_str)
	{
		strv splits;
		splitstr(splits,struct_str);
		vcnt = atoi(splits[2].c_str());
		ecnt = 0;	//atoi(splits[3].c_str());
	}
	void AddDimacsEdge(const string& edge_str)
	{
		// we dont care about the weights
		if(edge_str.size()<1||edge_str.at(0)!='e')
		return;

		int e1,e2;
		strv splits;
		splitstr(splits,edge_str);

		e1=atoi(splits[1].c_str())-1; //atoi is changing char to int, e1,e2 are int.
		e2=atoi(splits[2].c_str())-1; //Be careful. If -1 will to be max int value

		if(e1!=e2)	// take care of loops
		{
			AdjList[e1].insert(e2);
			AdjList[e2].insert(e1);
		}
	}
	void Dimacs(const string& path_str)
	{
		string str;
		ifstream fs(path_str.c_str());

		if(!fs)
		{
			cout<<"Cannot open!"<<endl;
			return;
		}

		do{getline(fs,str);} 
		while(str.at(0)!='p');

		ParseDimacs(str);

		AdjList=adj(vcnt);

		while(!fs.eof())
		{
			getline(fs,str);
			AddDimacsEdge(str);
		}

		fs.close();

		for(int i=0;i<vcnt;i++)
			ecnt += (int)AdjList[i].size();

		ecnt = ecnt/2;
	}

public:
	graph(void){}
	~graph(void){}

	graph(const string& path_str, int graph_type = 1, int thres_value = 0)
	{
		thres = thres_value;

		switch(graph_type)
		{
			case 1:
				Dimacs(path_str);
				break;
	
			default:
				cout<<"Constructing...";
		}
	}

	adj AdjList;
	int vcnt, ecnt;

	bool IsAdjTo(int a, int b)
	{
		return (AdjList[a].find(b)!=AdjList[a].end());
	}
};

//-------------------------------------------------------------------------------------//
// Matrix functions; complicated computation is based on MKL LAPACKE
void print_mat(int n, double* mat)
{
	for(int i =0; i<n; i++)
	{
		for(int j=0; j<n; j++)
			cout<<mat[i*n+j]<<"\t";
		cout<<endl;			
	}		
	cout<<endl;
}	
void MxV(double* M, douv& V, douv& prod)
{ // matrix times vector
	int n = (int)V.size(), i, j; // n is the dimesion 
								 // of the subproblem
	prod.assign(n, 0);

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (j >= i)
				prod[i] += M[i * n + j] * V[j];
			else
				prod[i] += M[j * n + i] * V[j];
		}
	}
}
double VxV(douv &V1, douv &V2)
{ // vector times vector
	int n = (int)V1.size(), i;
	double prod = 0;

	for(i=0;i<n;i++)
	{
		prod += V1[i]*V2[i];
	}

	return prod;
}
double mkl_min_eig(int n, double* mat)
{
	int m, info;
	double w[1];
	double *mat_cpy = new double[n*n];
	memcpy(mat_cpy, mat, sizeof(double)*n*n);	// keep mat as lapacke change mat	
	
	info = LAPACKE_dsyevr( MAT_FORMAT, 'N', 'I', 'U', n, mat_cpy, n,
		NULL, NULL, 1, 1, NULL, &m, w, NULL, n, NULL );
	
	/* Check for convergence */
	if (info > 0)
	{
		cout << "The algorithm failed to compute eigenvalues.\n";
		exit(1);
	}
	
	if(w[0] >= -eps)
	{
		cout<<"Caution: The min eigenvalues is non-negative!!!"<<endl;
		exit(1);
	}
	
	double r = w[0];	
	delete[] mat_cpy;
	
	return r;
}
void mkl_inv(int n, double* mat, double* mat_out)
{		
	int info = LAPACKE_dpotrf(MAT_FORMAT, 'U', n, mat_out, n);
	int info2 = LAPACKE_dpotri (MAT_FORMAT, 'U', n , mat_out, n);	
	
	/* Check for convergence */
	if (info > 0 || info2 > 0)
	{
		cout << "The algorithm failed to compute eigenvalues.\n";
		exit(1);
	}
}
void lift_mat(int n, double* mat, double tau, double* mat_out)
{
	memcpy(mat_out, mat, sizeof(double)*n*n);
	for(int i=0;i<n;i++)
		mat_out[i*n+i] += (2*tau);
}
/////////////////////////////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------------------//
// the qcp formulation, and its math information
class QCP
{
private:
	double* qbar;

public:
	int n; // total number of variables
	int ucnt; // unfixed variables in BB

	// tau is d in the new notations
	// hcoef is the nominator of the head  
	// coefficient of the nasty function
	double tau, lmd, hcoef, eta;

	double* q, * M;

	douv a, c;

	QCP(double gamma, graph& g)
	{
		int i;
		ints::iterator j;

		n = g.vcnt;
		q = new double[n * n];
		qbar = new double[n * n];
		M = new double[n * n];

		//---------------------------//
		// Parts of the formula other
		// than Q
		ucnt = n;
		eta = 0;
		a.assign(ucnt, 0);
		c.assign(ucnt, 1);
		///////////////////////////////

		for (i = 0; i < n * n; i++)
		{
			q[i] = -gamma;
			qbar[i] = gamma;
		}

		for (i = 0; i < n; i++)
		{
			q[i * n + i] = 0;
			qbar[i * n + i] = 0;
		}

		for (i = 0; i < n; i++)
		{
			for (j = g.AdjList[i].begin(); j != g.AdjList[i].end(); j++)
			{
				q[i * n + (*j)]++;
				qbar[i * n + (*j)]--;
			}
		}
	}
	~QCP(void)
	{
		delete[] q;
		delete[] qbar;
		delete[] M;
	}
	void EigM(double mueps)
	{
		lmd = mkl_min_eig(n, qbar);
		tau = -(1 + mueps) * lmd / 2;

		lift_mat(n, qbar, tau, M);
		mkl_inv(n, qbar, M);

		hcoef = lmd + 2 * tau;
	}
	void CalM(double mueps)
	{
		tau = -(1 + mueps) * lmd / 2;

		lift_mat(n, qbar, tau, M);
		mkl_inv(n, qbar, M);

		hcoef = lmd + 2 * tau;
	}
};
/////////////////////////////////////////////////////////////////////////////////////////


//-------------------------------------------------------------------------------------//
// LDB
//------------------------------------------------------------//
// double comparison
bool eq(double x, double y)
{
	return (x <= y+eps && x >= y-eps);
}
bool lgr(double x,double y)	
{
	return (x > y+eps);
}
bool les(double x,double y)	
{
	return (x < y-eps);
}
bool ge(double x,double y)	
{
	return (x >= y-eps);
}
bool le(double x,double y)	
{
	return (x <= y+eps);
}
//------------------------------------------------------------//
// Algorithm starts ...
struct triple
{   // a * u^-1 + b * u + c
	double a, b, c;
};
struct asc_by_phi
{	// sort phi
	const douv & phi;
	asc_by_phi(const douv &phi_input):phi(phi_input){}

	bool operator()(int a, int b)
	{
		return phi[a] < phi[b];
	}
};
class LDB //it is a calculator, create just onece
{
private:
	int n, m; // m is the number of brpts
	double ub_u_eq_0, hcoef, eta;

	douv alpha, beta; // They are Mc, Md respectively

	vector<triple> dsqr; // square of distance 
	triple g; // triple of function g_\mu

	intv brptci;	// current indecies of break points among m, used for sorting
	intv brpti;		// original indecies of break points among n	
	douv phi;		// break points, phi		

	void set_piece()
	{	// strictly speaking, should be distinct 
		// points, but we dont worry here

		// phi_ele is psi[i] = alpha[i]/(0.5-beta[i])
		// beta-1 since beta>=0.5
		double phi_ele, beta_minus_1;

		for (int i = 0; i < n; i++)
		{
			if (eq(alpha[i], 0))
			{ // case that alpha is 0				
				if (ge(beta[i], 0.5))
				{
					beta_minus_1 = beta[i] - 1;
					dsqr[i].a = 0;
					dsqr[i].b = hcoef * beta_minus_1 * beta_minus_1;
					dsqr[i].c = 0;
				}
				else
				{
					dsqr[i].a = 0;
					dsqr[i].b = hcoef * beta[i] * beta[i];
					dsqr[i].c = 0;
				}
			}
			else if (eq(beta[i], 0.5) || ((phi_ele = alpha[i] / (0.5 - beta[i])) < 0))
			{ // careful!!! we dont worry that sigma in [-eps,eps]
				if (lgr(alpha[i], 0))
				{
					beta_minus_1 = beta[i] - 1;
					dsqr[i].a = hcoef * alpha[i] * alpha[i];
					dsqr[i].b = hcoef * beta_minus_1 * beta_minus_1;
					dsqr[i].c = hcoef * 2 * alpha[i] * beta_minus_1;
				}
				else
				{
					dsqr[i].a = hcoef * alpha[i] * alpha[i];
					dsqr[i].b = hcoef * beta[i] * beta[i];
					dsqr[i].c = hcoef * 2 * alpha[i] * beta[i];
				}
			}
			else
			{
				brptci.push_back(m);
				m++;

				brpti.push_back(i);
				phi.push_back(phi_ele);
			}
		}

		sort(brptci.begin(), brptci.end(), asc_by_phi(phi));
	}

	//-----------------------------------------------------//
	// solve max h, the function h is now xi
	// return -infty if odd happens !!!
	double maxh_wo_piece()
	{	// LDB when there is no scenario
		// i.e. for m = 0
		// if return INF, infeasible

		double r;
		double a = 0, b = 0, c = 0; // the final expression

		vector<triple>::iterator i;
		for (i = dsqr.begin(); i != dsqr.end(); i++)
		{
			a += i->a;
			b += i->b;
			c += i->c;
		}

		a += g.a;
		b += (g.b + 2 * eta); // divide 2 later
		c += g.c;

		// totally 9 cases !!!
		if (les(a, 0) && les(b, 0)) // a < 0 and b < 0, concave
		{
			double spt = sqrt(a / b);	// stationary point		
			r = (a / spt + b * spt + c) / 2;
			// devided by 2 since first two parts have 1/2
			// moved it here
		}
		else if (eq(a, 0) && le(b, 0)) // a = 0 and b <= 0 
		{
			r = c / 2;	//sup when b<0, use u->0
		}
		else if (les(a, 0) && eq(b, 0)) // a < 0 and b = 0
		{
			r = c / 2;	// sup as u->inf
			cout << "A little odd 1"; // a->inf => 0
		}
		else if ((lgr(a, 0) && lgr(b, 0)) || (lgr(a, 0) && eq(b, 0)))
		{	// a>0, b>0 convex, rho->inf or rho->0 to max to inf
			// a>0, b<=0, u->0 to max to inf
			cout << "Odd 1"; // u->0 to max shall not happen
			r = -INF;	// set result to be negative infty
		}
		else {
			r = INF;
			//a<=0, b>0, u->inf to max to inf
		}

		return r;
	}
	double maxh_each_piece(double e1, double e2)
	{	// -1 is infinity end point
		// e2>0, then not the last piece
		// If return INF, infeasible
		double r, t1, t2, spt;
		double a = 0, b = 0, c = 0; // the final expression

		bool is1stpiece = false;
		if (e1 < -1)
		{
			is1stpiece = true;
			e1 = 0;
		}

		vector<triple>::iterator i;
		for (i = dsqr.begin(); i != dsqr.end(); i++)
		{
			a += i->a;
			b += i->b;
			c += i->c;
		}

		a += g.a;
		b += (g.b + 2 * eta); // divide 2 later
		c += g.c;

		if (lgr(a, 0) && le(b, 0))	// a>0, b<=0, decreasing curve
		{
			if (!is1stpiece)
				r = (a / e1 + b * e1 + c) / 2;
			else
			{ // rhi->0, max to infinity
				cout << "Odd 2";
				r = -INF;
			}
		}
		else if (eq(a, 0) && le(b, 0))	// a=0, b<=0, decreaing line
			r = (b * e1 + c) / 2;
		else if (le(a, 0) && lgr(b, 0)) // a<=0, b>0, increasing curve or line
		{
			if (e2 > 0)
				r = (a / e2 + b * e2 + c) / 2;
			else
				r = INF;//rho->infty to max to infty, infeasible
		}
		else if (les(a, 0) && eq(b, 0))	// a<0, b=0, increasing curve
		{
			if (e2 > 0)
				r = (a / e2 + c) / 2;
			else
			{
				r = c / 2;
				cout << "A little odd 2"; // a->inf => 0
			}
		}
		else if (lgr(a, 0) && lgr(b, 0))	// a>0, b>0 convex
		{
			if (!is1stpiece && e2 > 0) // either the first or the last piece
			{
				t1 = (a / e1 + b * e1 + c) / 2;
				t2 = (a / e2 + b * e2 + c) / 2;

				r = max(t1, t2);
			}
			else if (is1stpiece)
			{
				cout << "Odd 3"; // e1 is 0, positive infinity
				r = -INF;
			}
			else
				r = INF; //u->infty to max to infty, infeasible
		}
		else	// a<0, b<0 concave
		{
			spt = sqrt(a / b);	// stationary point
			if (spt < e1)
				r = (a / e1 + b * e1 + c) / 2;
			else if (spt > e2 && e2 > 0)
				r = (a / e2 + b * e2 + c) / 2;
			else
				r = (a / spt + b * spt + c) / 2;
		}

		return r;
	}
	double maxh(douv & subhs)
	{	// get the max h for all pieces
		// no need for w/o piece case
		douv::iterator i;
		double max_val = subhs[0];
		for (i = subhs.begin(); i != subhs.end(); i++)
		{
			if (*i > max_val)
				max_val = *i;
		}
		return max_val;
	}
	void update_obj_func(QCP & Q)
	{ // call before computing the bound once in line search
		m = 0;
		hcoef = Q.hcoef;
		douv d(n, Q.tau);

		MxV(Q.M, Q.c, alpha);
		MxV(Q.M, d, beta);

		dsqr.resize(n);
		brpti.clear();
		brptci.clear();
		phi.clear();

		g.a = -VxV(Q.c, alpha);
		g.b = -VxV(d, beta);
		g.c = -VxV(d, alpha) - VxV(Q.c, beta);
	}
	bool check_solve_one(double& bound_result)
	{	// if the return value is false infeasible for sure, if ture means 
		// not sure but upper bound is calculable
		// The upper bound has decimals !!!!

		// maxh functions return -infty if odd happens i.e. rho->0 => ->infty
		// so if ->infty must rho->0 and infeasible
		// we also take care of the case that "u = 0"

		int i, j, k;

		//-----------------------------------------------------------------//
		// note that o'= Md = beta; g = -d^TMd/2 = g.b/2
		// hcoef, eta are the same
		// we will devide by 2 later as we have been doing always
		double xi = 0, o_minus_1;

		for (i = 0; i < n; i++)
		{
			if (ge(beta[i], 0.5))
			{
				o_minus_1 = beta[i] - 1;
				xi += o_minus_1 * o_minus_1;
			}
			else {
				xi += beta[i] * beta[i];
			}
		}

		xi = (hcoef * xi + g.b) / 2;

		// see the second sufficient 
		// condition for infeasibility
		if (lgr(xi + eta, 0))
			return false;
		///////////////////////////////////////////////////////////////////////

		double beta_minus_1, t;
		douv subhs;

		set_piece();

		if (m == 0)
		{
			t = maxh_wo_piece();

			if (ge(t, INF))
				return false; // if ->infty, infeasible

			bound_result = min(-t, ub_u_eq_0);
			return true;
		}

		for (i = m; i >= 0; i--)
		{ // from iterval_m to interval_0
		  // interval_0 = ( 0, phi_0 ]
		  // phi_0 is the 1st one in sorted phi set
			for (j = 0; j < m; j++)
			{
				k = brpti[brptci[j]];
				if (j >= i)
				{
					if (les(alpha[k], 0))
					{
						dsqr[k].a = hcoef * alpha[k] * alpha[k];
						dsqr[k].b = hcoef * beta[k] * beta[k];
						dsqr[k].c = hcoef * 2 * alpha[k] * beta[k];
					}
					else
					{
						beta_minus_1 = beta[k] - 1;
						dsqr[k].a = hcoef * alpha[k] * alpha[k];
						dsqr[k].b = hcoef * beta_minus_1 * beta_minus_1;
						dsqr[k].c = hcoef * 2 * alpha[k] * beta_minus_1;
					}
				}
				else
				{
					if (les(alpha[k], 0))
					{
						beta_minus_1 = beta[k] - 1;
						dsqr[k].a = hcoef * alpha[k] * alpha[k];
						dsqr[k].b = hcoef * beta_minus_1 * beta_minus_1;
						dsqr[k].c = hcoef * 2 * alpha[k] * beta_minus_1;
					}
					else
					{
						dsqr[k].a = hcoef * alpha[k] * alpha[k];
						dsqr[k].b = hcoef * beta[k] * beta[k];
						dsqr[k].c = hcoef * 2 * alpha[k] * beta[k];
					}
				}
			}

			if (i == 0)
			{
				t = maxh_each_piece(-2, phi[brptci[i]]);
				// although use -2, but it is actually 0
				// always worry about the double 0 becomes eps
			}
			else if (i == m)
			{
				t = maxh_each_piece(phi[brptci[m - 1]], -1);	//-1 is the infinity end point
				if (ge(t, INF))
					return false;
			}
			else {
				t = maxh_each_piece(phi[brptci[i - 1]], phi[brptci[i]]);
			}

			subhs.push_back(t);
		}

		bound_result = min(-maxh(subhs), ub_u_eq_0);
		return true;
	}
public:
	LDB() {}
	~LDB() {}

	// Do not use bisection search, give a mueps 
	// also we dont care about the status of matrix
	// as we always assume that they are not positive definite
	// used to conduct experiment on large sparse graphs !!!
	double simple_solve(QCP & Q, int fix, double mueps, bool caleig)
	{
		//-------------------------------------------------//
		// Handle matrix case usiingline search
		// set unchanged variables in line seach first
		n = Q.ucnt; ub_u_eq_0 = n; eta = Q.eta;
		/////////////////////////////////////////////////////

		double result;
		bool IsFeasible;

		if (caleig)
			Q.EigM(mueps);
		else
			Q.CalM(mueps);

		update_obj_func(Q);
		IsFeasible = check_solve_one(result);
		if (!IsFeasible)
			return -1;
		else
		{
			result += fix;
			result = ceil((result * pow(10.0, 6)) - 0.49) / pow(10.0, 6); // round to 6th decimal for safety
			return result;
		}
	}

	// here we use a bisection search
	// and take the floor to make it integer !!!!
	// return -1 if infeasible, OK for this problem
	int bisec_solve(QCP & Q, int fix, int matstat)
	{
		//--------------------------------------//
		// deal with 1x1 and empty matrix cases
		if (matstat == -2)
		{
			if (lgr(Q.eta, 0))
				return -1;
			else
				return fix;
		}
		else if (matstat == -3)
		{
			if (lgr(Q.eta, Q.a[0]))
				return -1;
			else
				return Q.ucnt + fix;
		}
		/////////////////////////////////////////

		//-------------------------------------------------//
		// Handle matrix case usiingline search
		// set unchanged variables in line seach first
		n = Q.ucnt; ub_u_eq_0 = n; eta = Q.eta;
		/////////////////////////////////////////////////////

		double a = 0, b = 2, mid, c, d;
		double fc, fd;
		double result;

		int dir0_count = 0;
		bool _1st = true;
		bool IsFeasible;

		while (b - a >= 0.005)
		{
			mid = (a + b) / 2;
			c = mid - 0.002;
			d = mid + 0.002;

			if (_1st)
			{
				Q.EigM(c);
				//---------------------------------------//
				// handle 0 matrix!!!
				if (Q.tau < -0.5)
				{
					if (lgr(Q.eta, 0))
						return -1;
					else
						return Q.ucnt + fix;
				}
				///////////////////////////////////////////
				_1st = false;
			}
			else { Q.CalM(c); }

			update_obj_func(Q);
			IsFeasible = check_solve_one(fc);
			if (!IsFeasible) {
				return -1;
			}

			Q.CalM(d);
			update_obj_func(Q);
			IsFeasible = check_solve_one(fd);
			if (!IsFeasible) {
				return -1;
			}

			if (fabs(fc - fd) < 0.01)
			{
				dir0_count++;
				if (dir0_count > 3) { break; }
			}

			if (le(fc, fd))
				b = mid;
			else
				a = mid;
		}

		Q.CalM(mid);
		update_obj_func(Q);
		IsFeasible = check_solve_one(result);
		if (!IsFeasible)
			return -1;
		else
		{
			result += fix;
			result = ceil((result * pow(10.0, 6)) - 0.49) / pow(10.0, 6); // round to 6th decimal for safety
			return (int)result; // truncate to return integer floor
		}
	}
};
///////////////////////////////////////////////////////////////////////////////////////////////


void GetInput(strv& g_names)
{
	string str, ful_file;
	ful_file = "input/" + INPUT + ".txt";

	ifstream fs(ful_file.c_str());

	if (!fs)
	{
		cout << "Cannot open!" << endl;
		return;
	}

	while (!fs.eof())
	{
		getline(fs, str);
		if (str.size() > 0)
			g_names.push_back(str);
	}

	fs.close();
}
void batch_solve(strv& g_names)
{
	int gamma_cnt = sizeof(GAMMA) / sizeof(double);
	ofstream fout(OUTPUT.c_str());

	for (int k = 0; k < 5; k++)	//the 1st run is warm up then run 4 times to take average
	{
		cout << "Round " << k << endl;
		fout << "Round " << k << ":" << endl;
		for (int i = 0; i < gamma_cnt; i++)
		{
			fout << "Gamma = " << GAMMA[i] << ",ub,runtime" << endl;
			for (strv::iterator j = g_names.begin(); j != g_names.end(); j++)
			{
				graph g("input/" + INPUT + "/" + *j, 1);
				QCP qcp(GAMMA[i], g);	LDB ldb;

				double t = omp_get_wtime();
				double ub = ldb.simple_solve(qcp, 0, 0.01, true);
				//int ub = ldb.bisec_solve(qcp, 0, 1);
				t = omp_get_wtime() - t;

				fout << *j << "," << ub << "," << t << endl;
			}
			fout << endl;
		}
		fout << endl;
	}

	fout.close();
}

int main()
{
	strv gnames;
	GetInput(gnames);
	batch_solve(gnames);

#ifdef INWIN
	system("pause");
#endif

	return 0;
}
