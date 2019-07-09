#include<iostream>
#include<fstream>
#include<string.h>

#include<vector>
#include<set>

#include<algorithm>
#include<time.h>
#include<math.h>

#include "omp.h"

#include "gurobi_c++.h"

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
const double GAMMA[] = {0.85};

//-----------------------------------------------
// input & output, not changing often
const string INPUT = "dimacs";	//"sparse";
const string OUTPUT = "output/results.csv";
const string LDBINPUT = "input/bounds.csv";
const double TIMEL = 7200;
//------------------------------------------------


//---------------------------------------------------
// IPF parameters
//const char TYPE = GRB_CONTINUOUS;	
const char TYPE = GRB_BINARY;
const int WHICHF = 2;			// can be 2,3,4
//---------------------------------------------------


// Marco ------------------------------------------------------------------------
//#define INWIN
#define HASOMP
//#define USELDBCPM

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
// grb auxiliary functions
GRBVar* GVars(int n, GRBModel &model, double lb = 0, double ub = GRB_INFINITY, 
				 char type = GRB_CONTINUOUS, double obj_coeff = 0) 
{
	// delete after use
	GRBVar* vars = new GRBVar[n];

	for(int i=0;i<n;i++){
		vars[i] = model.addVar(lb, ub, obj_coeff, type);
	}

	return vars;
}
GRBLinExpr GSum(int n, GRBVar* vars)
{
	GRBLinExpr expr = 0;

	for(int i=0;i<n;i++){
		expr += vars[i];		
	}

	return expr;
}
GRBLinExpr GSum(int n, GRBVar* vars, douv &coeff)
{
	GRBLinExpr expr = 0;

	for(int i=0;i<n;i++){
		expr += (vars[i]*coeff[i]);	
	}

	return expr;
}
//-----------------------------------------------------------------------------------------------//

struct EdgePt
{
	int i,j;
};

void F2(graph &g, double gamma, double tl, double &obj, double &bound,  double &mtime, double &stime, double &ncnt)
{
	double t = clock();

	int n = g.vcnt, i, k;
	ints::iterator j;

	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag,0);
	env.set(GRB_DoubleParam_TimeLimit,tl);
	GRBModel model = GRBModel(env);	

	GRBVar* x = GVars(n, model, 0.0, 1.0, TYPE);
	GRBVar* y = GVars(n, model, -GRB_INFINITY);
	model.update();

	model.setObjective(GSum(n,x),GRB_MAXIMIZE);

	douv a(n), u(n), l(n);
	for(i=0;i<n;i++)
	{
		a[i] = -gamma;
		u[i] = (1-gamma) * ( g.AdjList[i].size() );
		l[i] = ( n - 1 - (g.AdjList[i].size()) )*(-gamma);
	}

	model.addConstr(GSum(n,y)>=0);
	for(i=0;i<n;i++)
	{
		model.addConstr( y[i] <= u[i]*x[i] );
		model.addConstr( y[i] >= l[i]*x[i] );
	}
	for(i=0;i<n;i++)
	{
		for(j=g.AdjList[i].begin();j!=g.AdjList[i].end();j++)
			a[*j]++;	

		model.addConstr( y[i] >= gamma*x[i] + GSum(n,x,a) - u[i]*(1-x[i]) );
		model.addConstr( y[i] <= gamma*x[i] + GSum(n,x,a) - l[i]*(1-x[i]) );

		for(k=0;k<n;k++){
			a[k] = -gamma;
		}
	}	
		
	//------------------------------------------//
	// Release the memory for variables
	delete[] y;
	delete[] x;
	//////////////////////////////////////////////

	mtime = (clock()-t)/CLOCKS_PER_SEC;

	model.optimize();
	
	obj = model.get(GRB_DoubleAttr_ObjVal);
	stime = model.get(GRB_DoubleAttr_Runtime);	

	bound = -1;
	if(TYPE == GRB_BINARY){
		bound = model.get(GRB_DoubleAttr_ObjBound);
	}
	
	ncnt = model.get(GRB_DoubleAttr_NodeCount);
}
void F3(graph &g, double gamma, int lb, int ub, double tl, double &obj, double &bound,  double &mtime, double &stime, double &ncnt)
{
#ifdef HASOMP
	double st = omp_get_wtime();
#else
	double st = clock();
#endif

	int i, n = g.vcnt, m = g.ecnt;
	int k = ub - lb + 1;

	vector<EdgePt> ends;
	ints::iterator j;

	for(i=0;i<n;i++)
	{
		for(j=g.AdjList[i].begin();
			j!=g.AdjList[i].end();j++)
		{
			if(*j>i)
			{
				EdgePt e;
				e.i = i;
				e.j = *j;
				ends.push_back(e);
			}
		}
	}

	GRBEnv *env = 0;
	env = new GRBEnv();
	env->set(GRB_IntParam_OutputFlag,0);
	env->set(GRB_DoubleParam_TimeLimit,tl);

    GRBModel model = GRBModel(*env);
	GRBVar* x = GVars(n, model, 0.0, 1.0, TYPE, -1);
	GRBVar* y = GVars(m, model);
	GRBVar* z = GVars(k, model);
	model.update();	
	
	//---------------------------------------------------//
	douv coefs;
	for(i=lb;i<=ub;i++)
	{
		double t = i*(i-1)/2;
		coefs.push_back(t);
	}
	GRBLinExpr lhs = GSum(m,y);
	GRBLinExpr rhs = GSum(k,z,coefs);
	model.addConstr(lhs >= gamma*rhs);

	for(i=0;i<m;i++)
	{
		model.addConstr(y[i]<=x[ends[i].i]);
		model.addConstr(y[i]<=x[ends[i].j]);
	}

	lhs =  GSum(n,x);
	for(i=0;i<k;i++){
		coefs[i] = lb+i;
	}
	rhs = GSum(k,z,coefs);
	model.addConstr(lhs == rhs);

	lhs =  GSum(k,z);
	model.addConstr(lhs == 1);

#ifdef HASOMP
	mtime = omp_get_wtime() - st;
#else
	mtime = (clock()-st)/CLOCKS_PER_SEC;
#endif

	model.optimize();

	stime = model.get(GRB_DoubleAttr_Runtime);
	
	obj = -model.get(GRB_DoubleAttr_ObjVal);
	bound = -1;
	if(TYPE == GRB_BINARY){
		bound = -model.get(GRB_DoubleAttr_ObjBound);
	}

	ncnt = model.get(GRB_DoubleAttr_NodeCount);

	delete[] x,y,z;
	delete env;
}
void F4(graph &g, double gamma, int lb, int ub, double tl, double &obj, double &bound,  double &mtime, double &stime, double &ncnt)
{
#ifdef HASOMP
	double st = omp_get_wtime();
#else
	double st = clock();
#endif

	int i, n = g.vcnt;	
	int k = ub - lb + 1;

	ints::iterator j;

	GRBEnv *env = 0;
	env = new GRBEnv();
	env->set(GRB_IntParam_OutputFlag,0);
	env->set(GRB_DoubleParam_TimeLimit,tl);

    GRBModel model = GRBModel(*env);
	GRBVar* x = GVars(n, model, 0.0, 1.0, TYPE, -1);
	GRBVar* v = GVars(n, model);
	GRBVar* z = GVars(k, model);
	model.update();

	//---------------------------------------------------//
	douv coefs;
	for(i=lb;i<=ub;i++)
	{
		double t = i*(i-1);
		coefs.push_back(t);
	}
	GRBLinExpr lhs = GSum(n,v);
	GRBLinExpr rhs = GSum(k,z,coefs);
	model.addConstr(lhs >= gamma*rhs);

	for(i=0;i<n;i++)
	{
		ints &cur = g.AdjList[i];
		double deg = (double)cur.size();
		model.addConstr( v[i] <= deg*x[i] );

		rhs = 0;
		for(j=cur.begin();j!=cur.end();j++){
			rhs += x[*j];
		}
		model.addConstr( v[i] <= rhs );
	}

	lhs =  GSum(n,x);
	for(i=0;i<k;i++){
		coefs[i] = lb+i;
	}
	rhs = GSum(k,z,coefs);
	model.addConstr(lhs == rhs);

	lhs =  GSum(k,z);
	model.addConstr(lhs == 1);

#ifdef HASOMP
	mtime = omp_get_wtime() - st;
#else
	mtime = (clock()-st)/CLOCKS_PER_SEC;
#endif

	model.optimize();

	stime = model.get(GRB_DoubleAttr_Runtime);
	
	obj = -model.get(GRB_DoubleAttr_ObjVal);
	bound = -1;
	if(TYPE == GRB_BINARY){
		bound = -model.get(GRB_DoubleAttr_ObjBound);
	}	

	ncnt = model.get(GRB_DoubleAttr_NodeCount);

	delete[] x,v,z;
	delete env;
}
///////////////////////////////////////////////////////////////////////////////////////////////////


void GetInput(strv& g_names)
{
	string str,ful_file;
	ful_file = "input/" + INPUT + ".txt";

	ifstream fs(ful_file.c_str());

	if(!fs)
	{
		cout<<"Cannot open!"<<endl;
		return;
	}

	while(!fs.eof())
	{
		getline(fs,str);
		if(str.size()>0)
			g_names.push_back(str);
	}

	fs.close();
}

//-----------------------------------------------------------------------------------------------------
void SolveOne(graph &g, double gamma, double tl, ofstream &fout)
{	// with analytical UB, and LB=2
	double obj = -1, bound = -1, mtime = -1, stime = -1, ncnt = 0;
	
	double t = clock();

	double dou_ub = 0.5 + 0.5 * sqrt(1 + 8.0 * g.ecnt / gamma); 
	int int_ub = (int)floor(dou_ub);

	t = (clock()-t)/CLOCKS_PER_SEC;
	
	if(WHICHF == 2)
	{
		F2(g,gamma,tl,obj,bound,mtime,stime,ncnt);
	}
	else if(WHICHF == 3)
	{		
		F3(g,gamma,2,int_ub,tl,obj,bound,mtime,stime,ncnt);
	}
	else if(WHICHF == 4)
	{
		F4(g,gamma,2,int_ub,tl,obj,bound,mtime,stime,ncnt);
	}
	else
	{
		cout<<"No other formulations!!!"<<endl;
	}

	fout<<2<<","<<0<<","<<int_ub<<","<<t<<","
		<<obj<<","<<bound<<","<<mtime<<","<<stime<<","<<ncnt<<endl;
}
void SolveOne(graph &g, double gamma, double tl, ofstream &fout, int int_lb, int int_ub)
{	// Use ldb as UB, and cpm as LB
	double obj = -1, bound = -1, mtime = -1, stime = -1, ncnt = 0;
	
	if(WHICHF == 2)
	{
		F2(g,gamma,tl,obj,bound,mtime,stime,ncnt);
	}
	else if(WHICHF == 3)
	{		
		F3(g,gamma,int_lb,int_ub,tl,obj,bound,mtime,stime,ncnt);
	}
	else if(WHICHF == 4)
	{
		F4(g,gamma,int_lb,int_ub,tl,obj,bound,mtime,stime,ncnt);
	}
	else
	{
		cout<<"No other formulations!!!"<<endl;
	}

	fout<<int_lb<<",0,"<<int_ub<<",0,"<<obj<<","<<bound<<","<<mtime<<","<<stime<<","<<ncnt<<endl;
}
//------------------------------------------------------------------------------------------------------

void loadbounds(intvv &bvals)
{	//lb ub appear alternatively!!!
	ifstream fs(LDBINPUT.c_str());

	if(!fs)
	{
		cout<<"Cannot open!"<<endl;
		return;
	}

	bvals.clear();

	int	i = -1; //indicate which 0.95 or 0.85 line
	string str;
	while(!fs.eof())
	{
		getline(fs,str);
		
		if(str.size()<1)
			continue;

		if(str[0]=='*')
		{
			bvals.push_back(intv());
			i++;
		}
		else
		{
			strv splits;
			splitstr(splits,str,",");

			int t = atoi(splits[0].c_str());
			bvals[i].push_back(t);

			t = atoi(splits[1].c_str());
			bvals[i].push_back(t);
		}
	}

	fs.close();
}
void SolveBatch(double tl)
{
#ifdef USELDBCPM
	intvv bvals;
	loadbounds(bvals);
#endif

	strv gnames;
	GetInput(gnames);

	ofstream fout(OUTPUT.c_str());

	if(WHICHF == 2)
	{
		fout<<"This is F2"<<endl;
	}
	else if(WHICHF == 3)
	{		
		fout<<"This is F3"<<endl;
	}
	else if(WHICHF == 4)
	{
		fout<<"This is F4"<<endl;
	}
	fout<<"gamma,graph,n,m,ini_lb,inilb_t,ini_ub,iniub_t,obj,bound,mtime,stime,ncnt"<<endl;

	for(int i=0;i<sizeof(GAMMA)/sizeof(double);i++)
	{	
		int k = 0;
		fout<<GAMMA[i]<<endl;
		cout<<GAMMA[i]<<endl;
		for(strv::iterator j=gnames.begin();j!=gnames.end();j++)
		{
			fout<<","<<*j;
			cout<<" "<<*j<<" solving ..."<<endl;
			graph g("input/" + INPUT + "/" + *j,1);				
			fout<<","<<g.vcnt<<","<<g.ecnt<<",";


	#ifdef USELDBCPM
			SolveOne(g,GAMMA[i],tl,fout,bvals[i][k],bvals[i][k+1]);
	#else
			SolveOne(g,GAMMA[i],tl,fout);
	#endif

			k+=2;
		}

		fout<<endl;
	}

	fout.close();
}

int main()
{
	SolveBatch(TIMEL);	

	//system("pause");

	return 0;
}
