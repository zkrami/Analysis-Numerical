#pragma once
#include "Poly.h"
typedef long long ll ;






using namespace std;

namespace Project2 {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	using namespace MetroFramework::Drawing; 
	using namespace MetroFramework;
	using namespace MetroFramework::Forms ; 
	using namespace MetroFramework::Controls;
	using namespace MetroFramework::Interfaces;
	using namespace System::Collections::Generic;




	/// <summary>
	/// Summary for Interpolation
	/// </summary>
	public ref class Interpolation : public MetroFramework::Forms::MetroForm
	{
	public:

		// Functions 
		ll factorial(int n){
			ll fact[1000];	
			fact[0] = fact[1] = 1;
			for (int i = 2; i < 50; i++)
				fact[i] = i * fact[i - 1];
			return fact[n];
		}
		spline_solution Spline(vector<pair<ld, ld>>& points) {
			Poly X;
			X.add(1, 1);
			spline_solution ans(points.size() - 1);
			for (int i = 0; i < (int)points.size() - 1; i++) {
				ans[i].p = (X - points[i].first) * ((points[i + 1].second - points[i].second) / (points[i + 1].first - points[i].first)) + points[i].second;
				ans[i].interval = make_pair(points[i].first, points[i + 1].first);
			}
			return ans;
		}

		bool is_const(vector<pair<ld, ld> >& points)
		{
			if (points.size() < 2) return 0;
			ld h = points[1].first - points[0].first;
			for (int i = 1; i < points.size(); i++)
			{
				if (abs(points[i].first - points[i - 1].first - h) > error)
					return 0;
			}				
			return 1;
		}
		vector<ld> gaussian_elimination(vector<vector<ld>>&aug) {
			int i, j, k, l;
			ld t;

			for (i = 0; i < (int)aug.size() - 1; i++) {
				l = i;
				for (j = i + 1; j < (int)aug.size(); j++)
					if (fabs(aug[j][i]) > fabs(aug[l][i]))
						l = j;
				for (k = i; k <= (int)aug.size(); k++)
					t = aug[i][k], aug[i][k] = aug[l][k], aug[l][k] = t;
				for (j = i + 1; j < (int)aug.size(); j++)
					for (k = aug.size(); k >= i; k--)
						aug[j][k] -= aug[i][k] * aug[j][i] / aug[i][i];
			}

			vector<ld> Ans(aug.size());
			for (j = aug.size() - 1; j >= 0; j--) {
				for (t = 0.0, k = j + 1; k < (int)aug.size(); k++) t += aug[j][k] * Ans[k];
				Ans[j] = (aug[j][aug.size()] - t) / aug[j][j];
			}
			return Ans;
		}


		Poly general_interpolation(vector<pair<ld, ld>> interpol)
		{
			vector<vector<ld> > ans((int)(interpol.size()), vector<ld>(1 + (int)(interpol.size())));
			for (int i = 0; i < interpol.size(); i++)
			{
				for (int j = 0; j < interpol.size(); j++)
				{
					ans[j][i] = pow(interpol[j].first, i);
				}
			}
			for (int i = 0; i < interpol.size(); i++)
			{
				ans[i][interpol.size()] = interpol[i].second;
			}
			return Poly(gaussian_elimination(ans));
		}
		Poly Polys_forward(vector<pair<ld, ld> >& points, int i)
		{

			Poly res = 1;
			for (int j = 0; j < i;j++)
			{
				Poly temp;
				temp.add(1, 1);
				temp -= points[j].first;
				res *= temp;
			}
			return res;
		}

		Poly Polys_backward(vector<pair<ld, ld> >& points, int i)
		{
			int n = points.size();

			Poly res = 1;
			for (int j = 0; j < i; j++)
			{
				Poly temp;
				temp.add(1, 1);
				temp -= points[n-j-1].first;
				res *= temp;
			}
			return res;
		}
		pair<Poly , vector<vector<ld>> >  Newton_Divided_Difference_forward(vector<pair<ld, ld>>& points)
		{
			int n = points.size();
			vector<vector<ld>> ans(points.size(), vector<ld>(points.size()));
			for (int i = 0; i < points.size(); i++)
			{
				ans[i][0] = points[i].second;
			}
			for (int j = 1; j < n; j++)
			{
				for (int i = j; i < n; i++)
				{

					ans[i][j] = (ans[i][j - 1] - ans[i - 1][j - 1]) / (points[i].first - points[i-j].first);
				}
			}

			Poly res;

			for (int i = 0; i < n; i++)
			{
				res += Polys_forward(points, i)*(ans[i][i]);
			}
			return make_pair(res,ans);
		}
		pair< Poly , pair < vector<vector<ld> >  , vector<vector<ld> >  > > squares( vector<pair<ld,ld> > & v , int n ) { 


			vector< vector<ld> > xk(v.size()+1 , vector<ld> (2*n -1 , 1) );
			vector<vector<ld> > yk(v.size() + 1 , vector<ld> (n));
			for(int j =0 ; j < 2*n-1; j++){
				if(j)
				{
					ld s = 0 ;
					for(int i = 0; i < v.size() ; i++ ){
						xk[i][j] = xk[i][j-1]*v[i].first ;
						s += xk[i][j];
					}
					xk[v.size()][j] = s ;

				}else {
					xk[v.size()][j] = v.size();
				}

			}
			for(int j = 0; j < n ; j++){

				ld s =0 ;
				for(int i =0 ; i < v.size() ; i++){
					yk[i][j] = v[i].second*xk[i][j] ;
					s += yk[i][j];
				}
				yk[v.size()][j] = s;
			}
			vector<vector<ld> > solve(n , vector<ld>(n+1));
			for(int i =0 ; i < n ; i++){

				for(int j = 0 ; j < n ; j++){
					solve[i][j] = xk[v.size()][i+j] ;
				}
				solve[i][n] = yk[v.size()][i] ;
			}

			Poly res = Poly(gaussian_elimination(solve)) ; 

			return make_pair(res , make_pair(xk,yk)) ; 

		}
		pair< Poly , vector<vector<ld>> >  Newton_Divided_Difference_backward(vector<pair<ld, ld>>& points)
		{
			int n = points.size();
			vector<vector<ld>> ans(points.size(), vector<ld>(points.size()));
			for (int i = 0; i < points.size(); i++)
			{
				ans[i][0] = points[i].second;
			}
			for (int j = 1; j < n; j++)
			{
				for (int i = j; i < n; i++)
				{

					ans[i][j] = (ans[i][j - 1] - ans[i - 1][j - 1]) / (points[i].first - points[i - j].first);
				}
			}

			Poly res;

			for (int i = 0; i < n; i++)
			{
				res += Polys_backward(points, i)*(ans[n-1][i]);
			}
			return make_pair(res,ans);
		}





		Poly Permutation(Poly p, int n) {
			Poly res;
			res.add(1, 0);
			for (int i = 0; i < n; i++) {
				res *= (p - i);
			}
			return res;
		}
		Poly Multi(Poly p, int n)
		{
			Poly res;
			res.add(1, 0);
			for (int i = 0; i < n; i++) {
				res *= (p + i);
			}
			return res;
		}
		pair<Poly,vector<vector<ld>> > Newton_backward(vector<pair<ld, ld>> points)
		{
			int n = points.size();
			vector<vector<ld>> ans(points.size(), vector<ld>(points.size()));
			for (int i = 0; i < points.size(); i++)
			{
				ans[i][0] = points[i].second;
			}
			for (int j = 1; j < n; j++)
			{
				for (int i = j; i < n; i++)
				{

					ans[i][j] = ans[i][j - 1] - ans[i - 1][j - 1];
				}
			}

			double h = points[1].first - points[0].first;
			Poly S;
			S.add(1, 1);
			S -= points[n-1].first;
			S /= h;

			Poly res;

			for (int i = 0; i < n; i++)
			{
				res += Multi(S, i)*(ans[n-1][i] / factorial(i));
			}
			return make_pair(res,ans);
		}


		pair<Poly,vector<vector<ld>> > Newton_forward(vector<pair<ld, ld>> points) {
			int n = points.size();
			vector<vector<ld>> ans(points.size(), vector<ld>(points.size()));
			for (int i = 0; i < points.size(); i++)
			{
				ans[i][0] = points[i].second;
			}
			for (int j = 1; j < n; j++)
			{
				for (int i = j; i < n; i++)
				{

					ans[i][j] = ans[i][j - 1] - ans[i - 1][j - 1];
				}
			}

			double h = points[1].first - points[0].first;
			Poly P;
			P.add(1, 1);
			P -= points[0].first;
			P /= h;

			Poly res;

			for (int i = 0; i < n; i++)
			{
				res += Permutation(P, i)*(ans[i][i] / factorial(i));
			}
			return make_pair(res,ans);
		}
		Poly Lagrange(vector<pair<ld, ld>> points) {
			Poly res, X;
			X.add(1,1);

			for (int i = 0; i < (int)points.size(); i++) {
				Poly nom(1);
				ld dom = 1;
				for (int j = 0; j < (int)points.size(); j++) {
					if (i != j) {
						nom *= (X - points[j].first);
						dom *= (points[i].first - points[j].first);
					}
				}
				nom /= dom;
				nom *= points[i].second;
				res += nom;
			}
			return res;
		}


		Interpolation(void)
		{
			InitializeComponent();

			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~Interpolation()
		{
			if (components)
			{
				delete components;
			}
		}
	private: MetroFramework::Controls::MetroGrid^  pointsGrid;
	protected: 
	private: System::Windows::Forms::DataGridViewTextBoxColumn^  x;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^  y;
	private: System::Windows::Forms::NumericUpDown^  pointsCountUpDown;
	private: MetroFramework::Controls::MetroTabControl^  tabControl;
	private: MetroFramework::Controls::MetroTabPage^  nuetonProgressiveTab;
	private: MetroFramework::Controls::MetroTabPage^  newtonBackword;

	private: MetroFramework::Controls::MetroTabPage^  NewtonDividedDifferencebackward;




	private: MetroFramework::Controls::MetroGrid^  nuetonProgressiveGrid;
	private: MetroFramework::Controls::MetroTabPage^  helpTab;


	private: MetroFramework::Controls::MetroLabel^  neutoProgressiveFail;
	private: MetroFramework::Controls::MetroLabel^  neutonProgressiveResult;

	private: MetroFramework::Controls::MetroLabel^  metroLabel1;
	private: MetroFramework::Controls::MetroTextBox^  neutonProgressiveValue;
	private: MetroFramework::Controls::MetroLabel^  nuetonProgreesivePolyLabel;

	private: MetroFramework::Controls::MetroLabel^  newtonBackwordResult;

	private: MetroFramework::Controls::MetroLabel^  metroLabel5;
	private: MetroFramework::Controls::MetroTextBox^  newtonBackwordValue;
	private: MetroFramework::Controls::MetroLabel^  newtonBackwordPolyLabel2;


	private: MetroFramework::Controls::MetroLabel^  newtonBackwordFail;
	private: MetroFramework::Controls::MetroGrid^  newtonBackwordGrid;
	private: MetroFramework::Controls::MetroGrid^  newtonDividedBackwordGrid;

	private: MetroFramework::Controls::MetroLabel^  metroLabel3;
	private: MetroFramework::Controls::MetroLabel^  metroLabel2;
	private: MetroFramework::Controls::MetroLabel^  metroLabel4;
	private: MetroFramework::Controls::MetroLabel^  newtonDividedBackwordResult;
	private: MetroFramework::Controls::MetroLabel^  metroLabel7;
	private: MetroFramework::Controls::MetroTextBox^  newtonDividedBackwordValue;
	private: MetroFramework::Controls::MetroLabel^  newtonDividedBackwordPolyLabel;
	private: MetroFramework::Controls::MetroTabPage^  NewtonDividedDifferenceforward;
	private: MetroFramework::Controls::MetroLabel^  metroLabel6;
	private: MetroFramework::Controls::MetroLabel^  newtonDividedForwardResult;
	private: MetroFramework::Controls::MetroLabel^  metroLabel9;
	private: MetroFramework::Controls::MetroTextBox^  newtonDividedForwardValue;
	private: MetroFramework::Controls::MetroLabel^  newtonDividedForwardPolyLabel;
	private: MetroFramework::Controls::MetroGrid^  newtonDividedForwardGrid;
	private: MetroFramework::Controls::MetroTabPage^  LagrangeTab;
	private: MetroFramework::Controls::MetroLabel^  metroLabel8;
	private: MetroFramework::Controls::MetroLabel^  lagrangeResult;
	private: MetroFramework::Controls::MetroLabel^  metroLabel11;
	private: MetroFramework::Controls::MetroTextBox^  lagrangeValue;
	private: MetroFramework::Controls::MetroLabel^  lagrangePolyLabel;
	private: MetroFramework::Controls::MetroTabPage^  generalInterpolationTab;
	private: MetroFramework::Controls::MetroLabel^  metroLabel10;
	private: MetroFramework::Controls::MetroLabel^  generalInterpolationResult;
	private: MetroFramework::Controls::MetroLabel^  metroLabel13;
	private: MetroFramework::Controls::MetroTextBox^  generalInterpolationValue;
	private: MetroFramework::Controls::MetroLabel^  generalInterpolationPolyLabel;
	private: System::Windows::Forms::TabPage^  leastSquareTab;
	private: MetroFramework::Controls::MetroLabel^  metroLabel14;
	private: System::Windows::Forms::NumericUpDown^  leastSquareDegree;
	private: MetroFramework::Controls::MetroLabel^  metroLabel12;
	private: MetroFramework::Controls::MetroLabel^  leastSquareResult;
	private: MetroFramework::Controls::MetroLabel^  metroLabel15;
	private: MetroFramework::Controls::MetroTextBox^  leastSquareValue;
	private: MetroFramework::Controls::MetroLabel^  leastSquarePolyLabel;
	private: MetroFramework::Controls::MetroGrid^  leastSquareYk;
	private: MetroFramework::Controls::MetroGrid^  leastSquareXk;
private: System::Windows::Forms::TabPage^  subLineTab;
private: MetroFramework::Controls::MetroLabel^  metroLabel16;
private: MetroFramework::Controls::MetroLabel^  subLineResult;

private: MetroFramework::Controls::MetroLabel^  metroLabel18;
private: MetroFramework::Controls::MetroTextBox^  subLineValue;

private: MetroFramework::Controls::MetroGrid^  subLineGrid;

private: System::Windows::Forms::DataGridViewTextBoxColumn^  column1;
private: MetroFramework::Controls::MetroLabel^  metroLabel17;
private: MetroFramework::Controls::MetroLabel^  metroLabel19;
private: MetroFramework::Controls::MetroPanel^  panel;
























	protected: 



	protected: 

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle1 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle2 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle3 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle4 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle5 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle6 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle7 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle8 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle9 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle10 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle11 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle12 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle13 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle14 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle15 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle16 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle17 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle18 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle19 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle20 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle21 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle22 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle23 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle24 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			this->pointsGrid = (gcnew MetroFramework::Controls::MetroGrid());
			this->x = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->y = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->pointsCountUpDown = (gcnew System::Windows::Forms::NumericUpDown());
			this->tabControl = (gcnew MetroFramework::Controls::MetroTabControl());
			this->helpTab = (gcnew MetroFramework::Controls::MetroTabPage());
			this->nuetonProgressiveTab = (gcnew MetroFramework::Controls::MetroTabPage());
			this->metroLabel3 = (gcnew MetroFramework::Controls::MetroLabel());
			this->neutonProgressiveResult = (gcnew MetroFramework::Controls::MetroLabel());
			this->metroLabel1 = (gcnew MetroFramework::Controls::MetroLabel());
			this->neutonProgressiveValue = (gcnew MetroFramework::Controls::MetroTextBox());
			this->nuetonProgreesivePolyLabel = (gcnew MetroFramework::Controls::MetroLabel());
			this->neutoProgressiveFail = (gcnew MetroFramework::Controls::MetroLabel());
			this->nuetonProgressiveGrid = (gcnew MetroFramework::Controls::MetroGrid());
			this->newtonBackword = (gcnew MetroFramework::Controls::MetroTabPage());
			this->metroLabel2 = (gcnew MetroFramework::Controls::MetroLabel());
			this->newtonBackwordResult = (gcnew MetroFramework::Controls::MetroLabel());
			this->metroLabel5 = (gcnew MetroFramework::Controls::MetroLabel());
			this->newtonBackwordValue = (gcnew MetroFramework::Controls::MetroTextBox());
			this->newtonBackwordPolyLabel2 = (gcnew MetroFramework::Controls::MetroLabel());
			this->newtonBackwordFail = (gcnew MetroFramework::Controls::MetroLabel());
			this->newtonBackwordGrid = (gcnew MetroFramework::Controls::MetroGrid());
			this->NewtonDividedDifferencebackward = (gcnew MetroFramework::Controls::MetroTabPage());
			this->metroLabel4 = (gcnew MetroFramework::Controls::MetroLabel());
			this->newtonDividedBackwordResult = (gcnew MetroFramework::Controls::MetroLabel());
			this->metroLabel7 = (gcnew MetroFramework::Controls::MetroLabel());
			this->newtonDividedBackwordValue = (gcnew MetroFramework::Controls::MetroTextBox());
			this->newtonDividedBackwordPolyLabel = (gcnew MetroFramework::Controls::MetroLabel());
			this->newtonDividedBackwordGrid = (gcnew MetroFramework::Controls::MetroGrid());
			this->NewtonDividedDifferenceforward = (gcnew MetroFramework::Controls::MetroTabPage());
			this->metroLabel6 = (gcnew MetroFramework::Controls::MetroLabel());
			this->newtonDividedForwardResult = (gcnew MetroFramework::Controls::MetroLabel());
			this->metroLabel9 = (gcnew MetroFramework::Controls::MetroLabel());
			this->newtonDividedForwardValue = (gcnew MetroFramework::Controls::MetroTextBox());
			this->newtonDividedForwardPolyLabel = (gcnew MetroFramework::Controls::MetroLabel());
			this->newtonDividedForwardGrid = (gcnew MetroFramework::Controls::MetroGrid());
			this->LagrangeTab = (gcnew MetroFramework::Controls::MetroTabPage());
			this->metroLabel8 = (gcnew MetroFramework::Controls::MetroLabel());
			this->lagrangeResult = (gcnew MetroFramework::Controls::MetroLabel());
			this->metroLabel11 = (gcnew MetroFramework::Controls::MetroLabel());
			this->lagrangeValue = (gcnew MetroFramework::Controls::MetroTextBox());
			this->lagrangePolyLabel = (gcnew MetroFramework::Controls::MetroLabel());
			this->generalInterpolationTab = (gcnew MetroFramework::Controls::MetroTabPage());
			this->metroLabel10 = (gcnew MetroFramework::Controls::MetroLabel());
			this->generalInterpolationResult = (gcnew MetroFramework::Controls::MetroLabel());
			this->metroLabel13 = (gcnew MetroFramework::Controls::MetroLabel());
			this->generalInterpolationValue = (gcnew MetroFramework::Controls::MetroTextBox());
			this->generalInterpolationPolyLabel = (gcnew MetroFramework::Controls::MetroLabel());
			this->leastSquareTab = (gcnew System::Windows::Forms::TabPage());
			this->metroLabel14 = (gcnew MetroFramework::Controls::MetroLabel());
			this->leastSquareDegree = (gcnew System::Windows::Forms::NumericUpDown());
			this->metroLabel12 = (gcnew MetroFramework::Controls::MetroLabel());
			this->leastSquareResult = (gcnew MetroFramework::Controls::MetroLabel());
			this->metroLabel15 = (gcnew MetroFramework::Controls::MetroLabel());
			this->leastSquareValue = (gcnew MetroFramework::Controls::MetroTextBox());
			this->leastSquarePolyLabel = (gcnew MetroFramework::Controls::MetroLabel());
			this->leastSquareYk = (gcnew MetroFramework::Controls::MetroGrid());
			this->leastSquareXk = (gcnew MetroFramework::Controls::MetroGrid());
			this->subLineTab = (gcnew System::Windows::Forms::TabPage());
			this->metroLabel16 = (gcnew MetroFramework::Controls::MetroLabel());
			this->subLineResult = (gcnew MetroFramework::Controls::MetroLabel());
			this->metroLabel18 = (gcnew MetroFramework::Controls::MetroLabel());
			this->subLineValue = (gcnew MetroFramework::Controls::MetroTextBox());
			this->subLineGrid = (gcnew MetroFramework::Controls::MetroGrid());
			this->column1 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->metroLabel17 = (gcnew MetroFramework::Controls::MetroLabel());
			this->metroLabel19 = (gcnew MetroFramework::Controls::MetroLabel());
			this->panel = (gcnew MetroFramework::Controls::MetroPanel());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->pointsGrid))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->pointsCountUpDown))->BeginInit();
			this->tabControl->SuspendLayout();
			this->helpTab->SuspendLayout();
			this->nuetonProgressiveTab->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->nuetonProgressiveGrid))->BeginInit();
			this->newtonBackword->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->newtonBackwordGrid))->BeginInit();
			this->NewtonDividedDifferencebackward->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->newtonDividedBackwordGrid))->BeginInit();
			this->NewtonDividedDifferenceforward->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->newtonDividedForwardGrid))->BeginInit();
			this->LagrangeTab->SuspendLayout();
			this->generalInterpolationTab->SuspendLayout();
			this->leastSquareTab->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->leastSquareDegree))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->leastSquareYk))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->leastSquareXk))->BeginInit();
			this->subLineTab->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->subLineGrid))->BeginInit();
			this->SuspendLayout();
			// 
			// pointsGrid
			// 
			this->pointsGrid->AllowUserToAddRows = false;
			this->pointsGrid->AllowUserToDeleteRows = false;
			this->pointsGrid->AllowUserToResizeRows = false;
			this->pointsGrid->BackgroundColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(250)), 
				static_cast<System::Int32>(static_cast<System::Byte>(250)), static_cast<System::Int32>(static_cast<System::Byte>(250)));
			this->pointsGrid->BorderStyle = System::Windows::Forms::BorderStyle::None;
			this->pointsGrid->CellBorderStyle = System::Windows::Forms::DataGridViewCellBorderStyle::None;
			this->pointsGrid->ColumnHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle1->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle1->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle1->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle1->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle1->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle1->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle1->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->pointsGrid->ColumnHeadersDefaultCellStyle = dataGridViewCellStyle1;
			this->pointsGrid->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->pointsGrid->Columns->AddRange(gcnew cli::array< System::Windows::Forms::DataGridViewColumn^  >(2) {this->x, this->y});
			dataGridViewCellStyle2->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle2->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle2->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle2->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(136)), 
				static_cast<System::Int32>(static_cast<System::Byte>(136)), static_cast<System::Int32>(static_cast<System::Byte>(136)));
			dataGridViewCellStyle2->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle2->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle2->WrapMode = System::Windows::Forms::DataGridViewTriState::False;
			this->pointsGrid->DefaultCellStyle = dataGridViewCellStyle2;
			this->pointsGrid->EnableHeadersVisualStyles = false;
			this->pointsGrid->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			this->pointsGrid->GridColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)));
			this->pointsGrid->Location = System::Drawing::Point(23, 91);
			this->pointsGrid->Name = L"pointsGrid";
			this->pointsGrid->RowHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle3->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle3->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle3->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle3->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle3->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle3->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle3->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->pointsGrid->RowHeadersDefaultCellStyle = dataGridViewCellStyle3;
			this->pointsGrid->RowHeadersWidthSizeMode = System::Windows::Forms::DataGridViewRowHeadersWidthSizeMode::DisableResizing;
			this->pointsGrid->SelectionMode = System::Windows::Forms::DataGridViewSelectionMode::FullRowSelect;
			this->pointsGrid->Size = System::Drawing::Size(192, 278);
			this->pointsGrid->TabIndex = 0;
			this->pointsGrid->CellValueChanged += gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &Interpolation::pointsGrid_CellValueChanged);
			// 
			// x
			// 
			this->x->HeaderText = L"X";
			this->x->Name = L"x";
			this->x->Width = 75;
			// 
			// y
			// 
			this->y->HeaderText = L"Y";
			this->y->Name = L"y";
			this->y->Width = 75;
			// 
			// pointsCountUpDown
			// 
			this->pointsCountUpDown->BorderStyle = System::Windows::Forms::BorderStyle::None;
			this->pointsCountUpDown->Location = System::Drawing::Point(24, 69);
			this->pointsCountUpDown->Minimum = System::Decimal(gcnew cli::array< System::Int32 >(4) {2, 0, 0, 0});
			this->pointsCountUpDown->Name = L"pointsCountUpDown";
			this->pointsCountUpDown->Size = System::Drawing::Size(191, 16);
			this->pointsCountUpDown->TabIndex = 1;
			this->pointsCountUpDown->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) {2, 0, 0, 0});
			this->pointsCountUpDown->ValueChanged += gcnew System::EventHandler(this, &Interpolation::pointsCountUpDown_ValueChanged);
			// 
			// tabControl
			// 
			this->tabControl->Controls->Add(this->helpTab);
			this->tabControl->Controls->Add(this->nuetonProgressiveTab);
			this->tabControl->Controls->Add(this->newtonBackword);
			this->tabControl->Controls->Add(this->NewtonDividedDifferencebackward);
			this->tabControl->Controls->Add(this->NewtonDividedDifferenceforward);
			this->tabControl->Controls->Add(this->LagrangeTab);
			this->tabControl->Controls->Add(this->generalInterpolationTab);
			this->tabControl->Controls->Add(this->leastSquareTab);
			this->tabControl->Controls->Add(this->subLineTab);
			this->tabControl->Location = System::Drawing::Point(230, 69);
			this->tabControl->Name = L"tabControl";
			this->tabControl->SelectedIndex = 0;
			this->tabControl->Size = System::Drawing::Size(709, 304);
			this->tabControl->TabIndex = 0;
			this->tabControl->UseSelectable = true;
			this->tabControl->SelectedIndexChanged += gcnew System::EventHandler(this, &Interpolation::metroTabControl1_SelectedIndexChanged);
			// 
			// helpTab
			// 
			this->helpTab->Controls->Add(this->metroLabel19);
			this->helpTab->Controls->Add(this->metroLabel17);
			this->helpTab->HorizontalScrollbarBarColor = true;
			this->helpTab->HorizontalScrollbarHighlightOnWheel = false;
			this->helpTab->HorizontalScrollbarSize = 10;
			this->helpTab->Location = System::Drawing::Point(4, 38);
			this->helpTab->Name = L"helpTab";
			this->helpTab->Size = System::Drawing::Size(701, 262);
			this->helpTab->TabIndex = 3;
			this->helpTab->Text = L"About";
			this->helpTab->VerticalScrollbarBarColor = true;
			this->helpTab->VerticalScrollbarHighlightOnWheel = false;
			this->helpTab->VerticalScrollbarSize = 10;
			// 
			// nuetonProgressiveTab
			// 
			this->nuetonProgressiveTab->Controls->Add(this->metroLabel3);
			this->nuetonProgressiveTab->Controls->Add(this->neutonProgressiveResult);
			this->nuetonProgressiveTab->Controls->Add(this->metroLabel1);
			this->nuetonProgressiveTab->Controls->Add(this->neutonProgressiveValue);
			this->nuetonProgressiveTab->Controls->Add(this->nuetonProgreesivePolyLabel);
			this->nuetonProgressiveTab->Controls->Add(this->neutoProgressiveFail);
			this->nuetonProgressiveTab->Controls->Add(this->nuetonProgressiveGrid);
			this->nuetonProgressiveTab->HorizontalScrollbarBarColor = true;
			this->nuetonProgressiveTab->HorizontalScrollbarHighlightOnWheel = false;
			this->nuetonProgressiveTab->HorizontalScrollbarSize = 10;
			this->nuetonProgressiveTab->Location = System::Drawing::Point(4, 38);
			this->nuetonProgressiveTab->Name = L"nuetonProgressiveTab";
			this->nuetonProgressiveTab->Size = System::Drawing::Size(701, 262);
			this->nuetonProgressiveTab->TabIndex = 0;
			this->nuetonProgressiveTab->Text = L"Newton Forward Method";
			this->nuetonProgressiveTab->VerticalScrollbarBarColor = true;
			this->nuetonProgressiveTab->VerticalScrollbarHighlightOnWheel = false;
			this->nuetonProgressiveTab->VerticalScrollbarSize = 10;
			this->nuetonProgressiveTab->Click += gcnew System::EventHandler(this, &Interpolation::nuetonProgressiveTab_Click);
			// 
			// metroLabel3
			// 
			this->metroLabel3->AutoSize = true;
			this->metroLabel3->Location = System::Drawing::Point(127, 240);
			this->metroLabel3->Name = L"metroLabel3";
			this->metroLabel3->Size = System::Drawing::Size(59, 19);
			this->metroLabel3->TabIndex = 9;
			this->metroLabel3->Text = L"P(X) => ";
			// 
			// neutonProgressiveResult
			// 
			this->neutonProgressiveResult->AutoSize = true;
			this->neutonProgressiveResult->Location = System::Drawing::Point(192, 240);
			this->neutonProgressiveResult->Name = L"neutonProgressiveResult";
			this->neutonProgressiveResult->Size = System::Drawing::Size(13, 19);
			this->neutonProgressiveResult->TabIndex = 8;
			this->neutonProgressiveResult->Text = L" ";
			// 
			// metroLabel1
			// 
			this->metroLabel1->AutoSize = true;
			this->metroLabel1->Location = System::Drawing::Point(6, 236);
			this->metroLabel1->Name = L"metroLabel1";
			this->metroLabel1->Size = System::Drawing::Size(34, 19);
			this->metroLabel1->TabIndex = 6;
			this->metroLabel1->Text = L"X = ";
			// 
			// neutonProgressiveValue
			// 
			this->neutonProgressiveValue->Lines = gcnew cli::array< System::String^  >(1) {L"0"};
			this->neutonProgressiveValue->Location = System::Drawing::Point(46, 236);
			this->neutonProgressiveValue->MaxLength = 32767;
			this->neutonProgressiveValue->Name = L"neutonProgressiveValue";
			this->neutonProgressiveValue->PasswordChar = '\0';
			this->neutonProgressiveValue->ScrollBars = System::Windows::Forms::ScrollBars::None;
			this->neutonProgressiveValue->SelectedText = L"";
			this->neutonProgressiveValue->Size = System::Drawing::Size(75, 23);
			this->neutonProgressiveValue->TabIndex = 5;
			this->neutonProgressiveValue->Text = L"0";
			this->neutonProgressiveValue->UseSelectable = true;
			this->neutonProgressiveValue->TextChanged += gcnew System::EventHandler(this, &Interpolation::neutonProgressiveValue_TextChanged);
			this->neutonProgressiveValue->Click += gcnew System::EventHandler(this, &Interpolation::neutonProgressiveValue_Click);
			// 
			// nuetonProgreesivePolyLabel
			// 
			this->nuetonProgreesivePolyLabel->AutoSize = true;
			this->nuetonProgreesivePolyLabel->Location = System::Drawing::Point(6, 210);
			this->nuetonProgreesivePolyLabel->Name = L"nuetonProgreesivePolyLabel";
			this->nuetonProgreesivePolyLabel->Size = System::Drawing::Size(48, 19);
			this->nuetonProgreesivePolyLabel->TabIndex = 4;
			this->nuetonProgreesivePolyLabel->Text = L"P(x) = ";
			// 
			// neutoProgressiveFail
			// 
			this->neutoProgressiveFail->AutoSize = true;
			this->neutoProgressiveFail->Location = System::Drawing::Point(189, 99);
			this->neutoProgressiveFail->Name = L"neutoProgressiveFail";
			this->neutoProgressiveFail->Size = System::Drawing::Size(250, 19);
			this->neutoProgressiveFail->TabIndex = 3;
			this->neutoProgressiveFail->Text = L"the difference is not const between points";
			// 
			// nuetonProgressiveGrid
			// 
			this->nuetonProgressiveGrid->AllowUserToAddRows = false;
			this->nuetonProgressiveGrid->AllowUserToDeleteRows = false;
			this->nuetonProgressiveGrid->AllowUserToResizeRows = false;
			this->nuetonProgressiveGrid->BackgroundColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(250)), 
				static_cast<System::Int32>(static_cast<System::Byte>(250)), static_cast<System::Int32>(static_cast<System::Byte>(250)));
			this->nuetonProgressiveGrid->BorderStyle = System::Windows::Forms::BorderStyle::None;
			this->nuetonProgressiveGrid->CellBorderStyle = System::Windows::Forms::DataGridViewCellBorderStyle::None;
			this->nuetonProgressiveGrid->ColumnHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle4->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle4->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle4->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle4->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle4->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle4->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle4->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->nuetonProgressiveGrid->ColumnHeadersDefaultCellStyle = dataGridViewCellStyle4;
			this->nuetonProgressiveGrid->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			dataGridViewCellStyle5->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle5->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle5->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle5->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(136)), 
				static_cast<System::Int32>(static_cast<System::Byte>(136)), static_cast<System::Int32>(static_cast<System::Byte>(136)));
			dataGridViewCellStyle5->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle5->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle5->WrapMode = System::Windows::Forms::DataGridViewTriState::False;
			this->nuetonProgressiveGrid->DefaultCellStyle = dataGridViewCellStyle5;
			this->nuetonProgressiveGrid->EnableHeadersVisualStyles = false;
			this->nuetonProgressiveGrid->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			this->nuetonProgressiveGrid->GridColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			this->nuetonProgressiveGrid->Location = System::Drawing::Point(6, 15);
			this->nuetonProgressiveGrid->Name = L"nuetonProgressiveGrid";
			this->nuetonProgressiveGrid->ReadOnly = true;
			this->nuetonProgressiveGrid->RowHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle6->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle6->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle6->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle6->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle6->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle6->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle6->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->nuetonProgressiveGrid->RowHeadersDefaultCellStyle = dataGridViewCellStyle6;
			this->nuetonProgressiveGrid->RowHeadersWidthSizeMode = System::Windows::Forms::DataGridViewRowHeadersWidthSizeMode::DisableResizing;
			this->nuetonProgressiveGrid->SelectionMode = System::Windows::Forms::DataGridViewSelectionMode::FullRowSelect;
			this->nuetonProgressiveGrid->Size = System::Drawing::Size(695, 192);
			this->nuetonProgressiveGrid->TabIndex = 2;
			// 
			// newtonBackword
			// 
			this->newtonBackword->Controls->Add(this->metroLabel2);
			this->newtonBackword->Controls->Add(this->newtonBackwordResult);
			this->newtonBackword->Controls->Add(this->metroLabel5);
			this->newtonBackword->Controls->Add(this->newtonBackwordValue);
			this->newtonBackword->Controls->Add(this->newtonBackwordPolyLabel2);
			this->newtonBackword->Controls->Add(this->newtonBackwordFail);
			this->newtonBackword->Controls->Add(this->newtonBackwordGrid);
			this->newtonBackword->HorizontalScrollbarBarColor = true;
			this->newtonBackword->HorizontalScrollbarHighlightOnWheel = false;
			this->newtonBackword->HorizontalScrollbarSize = 10;
			this->newtonBackword->Location = System::Drawing::Point(4, 38);
			this->newtonBackword->Name = L"newtonBackword";
			this->newtonBackword->Size = System::Drawing::Size(701, 229);
			this->newtonBackword->TabIndex = 1;
			this->newtonBackword->Text = L"Newton Backward Method";
			this->newtonBackword->VerticalScrollbarBarColor = true;
			this->newtonBackword->VerticalScrollbarHighlightOnWheel = false;
			this->newtonBackword->VerticalScrollbarSize = 10;
			// 
			// metroLabel2
			// 
			this->metroLabel2->AutoSize = true;
			this->metroLabel2->Location = System::Drawing::Point(124, 235);
			this->metroLabel2->Name = L"metroLabel2";
			this->metroLabel2->Size = System::Drawing::Size(59, 19);
			this->metroLabel2->TabIndex = 15;
			this->metroLabel2->Text = L"P(X) => ";
			// 
			// newtonBackwordResult
			// 
			this->newtonBackwordResult->AutoSize = true;
			this->newtonBackwordResult->Location = System::Drawing::Point(189, 235);
			this->newtonBackwordResult->Name = L"newtonBackwordResult";
			this->newtonBackwordResult->Size = System::Drawing::Size(0, 0);
			this->newtonBackwordResult->TabIndex = 14;
			// 
			// metroLabel5
			// 
			this->metroLabel5->AutoSize = true;
			this->metroLabel5->Location = System::Drawing::Point(3, 235);
			this->metroLabel5->Name = L"metroLabel5";
			this->metroLabel5->Size = System::Drawing::Size(34, 19);
			this->metroLabel5->TabIndex = 13;
			this->metroLabel5->Text = L"X = ";
			// 
			// newtonBackwordValue
			// 
			this->newtonBackwordValue->Lines = gcnew cli::array< System::String^  >(1) {L"0"};
			this->newtonBackwordValue->Location = System::Drawing::Point(43, 235);
			this->newtonBackwordValue->MaxLength = 32767;
			this->newtonBackwordValue->Name = L"newtonBackwordValue";
			this->newtonBackwordValue->PasswordChar = '\0';
			this->newtonBackwordValue->ScrollBars = System::Windows::Forms::ScrollBars::None;
			this->newtonBackwordValue->SelectedText = L"";
			this->newtonBackwordValue->Size = System::Drawing::Size(75, 23);
			this->newtonBackwordValue->TabIndex = 12;
			this->newtonBackwordValue->Text = L"0";
			this->newtonBackwordValue->UseSelectable = true;
			this->newtonBackwordValue->TextChanged += gcnew System::EventHandler(this, &Interpolation::newtonBackwordValue_TextChanged);
			// 
			// newtonBackwordPolyLabel2
			// 
			this->newtonBackwordPolyLabel2->AutoSize = true;
			this->newtonBackwordPolyLabel2->Location = System::Drawing::Point(3, 209);
			this->newtonBackwordPolyLabel2->Name = L"newtonBackwordPolyLabel2";
			this->newtonBackwordPolyLabel2->Size = System::Drawing::Size(48, 19);
			this->newtonBackwordPolyLabel2->TabIndex = 11;
			this->newtonBackwordPolyLabel2->Text = L"P(x) = ";
			// 
			// newtonBackwordFail
			// 
			this->newtonBackwordFail->AutoSize = true;
			this->newtonBackwordFail->Location = System::Drawing::Point(186, 98);
			this->newtonBackwordFail->Name = L"newtonBackwordFail";
			this->newtonBackwordFail->Size = System::Drawing::Size(250, 19);
			this->newtonBackwordFail->TabIndex = 10;
			this->newtonBackwordFail->Text = L"the difference is not const between points";
			// 
			// newtonBackwordGrid
			// 
			this->newtonBackwordGrid->AllowUserToAddRows = false;
			this->newtonBackwordGrid->AllowUserToDeleteRows = false;
			this->newtonBackwordGrid->AllowUserToResizeRows = false;
			this->newtonBackwordGrid->BackgroundColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(250)), 
				static_cast<System::Int32>(static_cast<System::Byte>(250)), static_cast<System::Int32>(static_cast<System::Byte>(250)));
			this->newtonBackwordGrid->BorderStyle = System::Windows::Forms::BorderStyle::None;
			this->newtonBackwordGrid->CellBorderStyle = System::Windows::Forms::DataGridViewCellBorderStyle::None;
			this->newtonBackwordGrid->ColumnHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle7->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle7->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle7->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle7->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle7->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle7->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle7->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->newtonBackwordGrid->ColumnHeadersDefaultCellStyle = dataGridViewCellStyle7;
			this->newtonBackwordGrid->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			dataGridViewCellStyle8->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle8->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle8->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle8->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(136)), 
				static_cast<System::Int32>(static_cast<System::Byte>(136)), static_cast<System::Int32>(static_cast<System::Byte>(136)));
			dataGridViewCellStyle8->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle8->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle8->WrapMode = System::Windows::Forms::DataGridViewTriState::False;
			this->newtonBackwordGrid->DefaultCellStyle = dataGridViewCellStyle8;
			this->newtonBackwordGrid->EnableHeadersVisualStyles = false;
			this->newtonBackwordGrid->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			this->newtonBackwordGrid->GridColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			this->newtonBackwordGrid->Location = System::Drawing::Point(3, 14);
			this->newtonBackwordGrid->Name = L"newtonBackwordGrid";
			this->newtonBackwordGrid->ReadOnly = true;
			this->newtonBackwordGrid->RowHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle9->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle9->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle9->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle9->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle9->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle9->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle9->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->newtonBackwordGrid->RowHeadersDefaultCellStyle = dataGridViewCellStyle9;
			this->newtonBackwordGrid->RowHeadersWidthSizeMode = System::Windows::Forms::DataGridViewRowHeadersWidthSizeMode::DisableResizing;
			this->newtonBackwordGrid->SelectionMode = System::Windows::Forms::DataGridViewSelectionMode::FullRowSelect;
			this->newtonBackwordGrid->Size = System::Drawing::Size(695, 192);
			this->newtonBackwordGrid->TabIndex = 9;
			// 
			// NewtonDividedDifferencebackward
			// 
			this->NewtonDividedDifferencebackward->Controls->Add(this->metroLabel4);
			this->NewtonDividedDifferencebackward->Controls->Add(this->newtonDividedBackwordResult);
			this->NewtonDividedDifferencebackward->Controls->Add(this->metroLabel7);
			this->NewtonDividedDifferencebackward->Controls->Add(this->newtonDividedBackwordValue);
			this->NewtonDividedDifferencebackward->Controls->Add(this->newtonDividedBackwordPolyLabel);
			this->NewtonDividedDifferencebackward->Controls->Add(this->newtonDividedBackwordGrid);
			this->NewtonDividedDifferencebackward->HorizontalScrollbarBarColor = true;
			this->NewtonDividedDifferencebackward->HorizontalScrollbarHighlightOnWheel = false;
			this->NewtonDividedDifferencebackward->HorizontalScrollbarSize = 10;
			this->NewtonDividedDifferencebackward->Location = System::Drawing::Point(4, 38);
			this->NewtonDividedDifferencebackward->Name = L"NewtonDividedDifferencebackward";
			this->NewtonDividedDifferencebackward->Size = System::Drawing::Size(701, 229);
			this->NewtonDividedDifferencebackward->TabIndex = 2;
			this->NewtonDividedDifferencebackward->Text = L"Newton Backward-Divided-Difference Method";
			this->NewtonDividedDifferencebackward->VerticalScrollbarBarColor = true;
			this->NewtonDividedDifferencebackward->VerticalScrollbarHighlightOnWheel = false;
			this->NewtonDividedDifferencebackward->VerticalScrollbarSize = 10;
			// 
			// metroLabel4
			// 
			this->metroLabel4->AutoSize = true;
			this->metroLabel4->Location = System::Drawing::Point(125, 243);
			this->metroLabel4->Name = L"metroLabel4";
			this->metroLabel4->Size = System::Drawing::Size(59, 19);
			this->metroLabel4->TabIndex = 20;
			this->metroLabel4->Text = L"P(X) => ";
			// 
			// newtonDividedBackwordResult
			// 
			this->newtonDividedBackwordResult->AutoSize = true;
			this->newtonDividedBackwordResult->Location = System::Drawing::Point(190, 243);
			this->newtonDividedBackwordResult->Name = L"newtonDividedBackwordResult";
			this->newtonDividedBackwordResult->Size = System::Drawing::Size(0, 0);
			this->newtonDividedBackwordResult->TabIndex = 19;
			// 
			// metroLabel7
			// 
			this->metroLabel7->AutoSize = true;
			this->metroLabel7->Location = System::Drawing::Point(4, 243);
			this->metroLabel7->Name = L"metroLabel7";
			this->metroLabel7->Size = System::Drawing::Size(34, 19);
			this->metroLabel7->TabIndex = 18;
			this->metroLabel7->Text = L"X = ";
			// 
			// newtonDividedBackwordValue
			// 
			this->newtonDividedBackwordValue->Lines = gcnew cli::array< System::String^  >(1) {L"0"};
			this->newtonDividedBackwordValue->Location = System::Drawing::Point(44, 243);
			this->newtonDividedBackwordValue->MaxLength = 32767;
			this->newtonDividedBackwordValue->Name = L"newtonDividedBackwordValue";
			this->newtonDividedBackwordValue->PasswordChar = '\0';
			this->newtonDividedBackwordValue->ScrollBars = System::Windows::Forms::ScrollBars::None;
			this->newtonDividedBackwordValue->SelectedText = L"";
			this->newtonDividedBackwordValue->Size = System::Drawing::Size(75, 23);
			this->newtonDividedBackwordValue->TabIndex = 17;
			this->newtonDividedBackwordValue->Text = L"0";
			this->newtonDividedBackwordValue->UseSelectable = true;
			this->newtonDividedBackwordValue->TextChanged += gcnew System::EventHandler(this, &Interpolation::newtonDividedBackwordValue_TextChanged);
			// 
			// newtonDividedBackwordPolyLabel
			// 
			this->newtonDividedBackwordPolyLabel->AutoSize = true;
			this->newtonDividedBackwordPolyLabel->Location = System::Drawing::Point(4, 217);
			this->newtonDividedBackwordPolyLabel->Name = L"newtonDividedBackwordPolyLabel";
			this->newtonDividedBackwordPolyLabel->Size = System::Drawing::Size(48, 19);
			this->newtonDividedBackwordPolyLabel->TabIndex = 16;
			this->newtonDividedBackwordPolyLabel->Text = L"P(x) = ";
			// 
			// newtonDividedBackwordGrid
			// 
			this->newtonDividedBackwordGrid->AllowUserToAddRows = false;
			this->newtonDividedBackwordGrid->AllowUserToDeleteRows = false;
			this->newtonDividedBackwordGrid->AllowUserToResizeRows = false;
			this->newtonDividedBackwordGrid->BackgroundColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(250)), 
				static_cast<System::Int32>(static_cast<System::Byte>(250)), static_cast<System::Int32>(static_cast<System::Byte>(250)));
			this->newtonDividedBackwordGrid->BorderStyle = System::Windows::Forms::BorderStyle::None;
			this->newtonDividedBackwordGrid->CellBorderStyle = System::Windows::Forms::DataGridViewCellBorderStyle::None;
			this->newtonDividedBackwordGrid->ColumnHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle10->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle10->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle10->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle10->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle10->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle10->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle10->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->newtonDividedBackwordGrid->ColumnHeadersDefaultCellStyle = dataGridViewCellStyle10;
			this->newtonDividedBackwordGrid->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			dataGridViewCellStyle11->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle11->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle11->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle11->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(136)), 
				static_cast<System::Int32>(static_cast<System::Byte>(136)), static_cast<System::Int32>(static_cast<System::Byte>(136)));
			dataGridViewCellStyle11->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle11->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle11->WrapMode = System::Windows::Forms::DataGridViewTriState::False;
			this->newtonDividedBackwordGrid->DefaultCellStyle = dataGridViewCellStyle11;
			this->newtonDividedBackwordGrid->EnableHeadersVisualStyles = false;
			this->newtonDividedBackwordGrid->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			this->newtonDividedBackwordGrid->GridColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			this->newtonDividedBackwordGrid->Location = System::Drawing::Point(3, 12);
			this->newtonDividedBackwordGrid->Name = L"newtonDividedBackwordGrid";
			this->newtonDividedBackwordGrid->ReadOnly = true;
			this->newtonDividedBackwordGrid->RowHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle12->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle12->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle12->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle12->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle12->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle12->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle12->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->newtonDividedBackwordGrid->RowHeadersDefaultCellStyle = dataGridViewCellStyle12;
			this->newtonDividedBackwordGrid->RowHeadersWidthSizeMode = System::Windows::Forms::DataGridViewRowHeadersWidthSizeMode::DisableResizing;
			this->newtonDividedBackwordGrid->SelectionMode = System::Windows::Forms::DataGridViewSelectionMode::FullRowSelect;
			this->newtonDividedBackwordGrid->Size = System::Drawing::Size(695, 192);
			this->newtonDividedBackwordGrid->TabIndex = 10;
			// 
			// NewtonDividedDifferenceforward
			// 
			this->NewtonDividedDifferenceforward->Controls->Add(this->metroLabel6);
			this->NewtonDividedDifferenceforward->Controls->Add(this->newtonDividedForwardResult);
			this->NewtonDividedDifferenceforward->Controls->Add(this->metroLabel9);
			this->NewtonDividedDifferenceforward->Controls->Add(this->newtonDividedForwardValue);
			this->NewtonDividedDifferenceforward->Controls->Add(this->newtonDividedForwardPolyLabel);
			this->NewtonDividedDifferenceforward->Controls->Add(this->newtonDividedForwardGrid);
			this->NewtonDividedDifferenceforward->HorizontalScrollbarBarColor = true;
			this->NewtonDividedDifferenceforward->HorizontalScrollbarHighlightOnWheel = false;
			this->NewtonDividedDifferenceforward->HorizontalScrollbarSize = 10;
			this->NewtonDividedDifferenceforward->Location = System::Drawing::Point(4, 38);
			this->NewtonDividedDifferenceforward->Name = L"NewtonDividedDifferenceforward";
			this->NewtonDividedDifferenceforward->Size = System::Drawing::Size(701, 346);
			this->NewtonDividedDifferenceforward->TabIndex = 4;
			this->NewtonDividedDifferenceforward->Text = L"Newton Forward-Divided-Difference Method";
			this->NewtonDividedDifferenceforward->VerticalScrollbarBarColor = true;
			this->NewtonDividedDifferenceforward->VerticalScrollbarHighlightOnWheel = false;
			this->NewtonDividedDifferenceforward->VerticalScrollbarSize = 10;
			// 
			// metroLabel6
			// 
			this->metroLabel6->AutoSize = true;
			this->metroLabel6->Location = System::Drawing::Point(125, 240);
			this->metroLabel6->Name = L"metroLabel6";
			this->metroLabel6->Size = System::Drawing::Size(59, 19);
			this->metroLabel6->TabIndex = 26;
			this->metroLabel6->Text = L"P(X) => ";
			// 
			// newtonDividedForwardResult
			// 
			this->newtonDividedForwardResult->AutoSize = true;
			this->newtonDividedForwardResult->Location = System::Drawing::Point(190, 240);
			this->newtonDividedForwardResult->Name = L"newtonDividedForwardResult";
			this->newtonDividedForwardResult->Size = System::Drawing::Size(0, 0);
			this->newtonDividedForwardResult->TabIndex = 25;
			// 
			// metroLabel9
			// 
			this->metroLabel9->AutoSize = true;
			this->metroLabel9->Location = System::Drawing::Point(4, 240);
			this->metroLabel9->Name = L"metroLabel9";
			this->metroLabel9->Size = System::Drawing::Size(34, 19);
			this->metroLabel9->TabIndex = 24;
			this->metroLabel9->Text = L"X = ";
			// 
			// newtonDividedForwardValue
			// 
			this->newtonDividedForwardValue->Lines = gcnew cli::array< System::String^  >(1) {L"0"};
			this->newtonDividedForwardValue->Location = System::Drawing::Point(44, 240);
			this->newtonDividedForwardValue->MaxLength = 32767;
			this->newtonDividedForwardValue->Name = L"newtonDividedForwardValue";
			this->newtonDividedForwardValue->PasswordChar = '\0';
			this->newtonDividedForwardValue->ScrollBars = System::Windows::Forms::ScrollBars::None;
			this->newtonDividedForwardValue->SelectedText = L"";
			this->newtonDividedForwardValue->Size = System::Drawing::Size(75, 23);
			this->newtonDividedForwardValue->TabIndex = 23;
			this->newtonDividedForwardValue->Text = L"0";
			this->newtonDividedForwardValue->UseSelectable = true;
			this->newtonDividedForwardValue->TextChanged += gcnew System::EventHandler(this, &Interpolation::newtonDividedForwardValue_TextChanged);
			// 
			// newtonDividedForwardPolyLabel
			// 
			this->newtonDividedForwardPolyLabel->AutoSize = true;
			this->newtonDividedForwardPolyLabel->Location = System::Drawing::Point(4, 214);
			this->newtonDividedForwardPolyLabel->Name = L"newtonDividedForwardPolyLabel";
			this->newtonDividedForwardPolyLabel->Size = System::Drawing::Size(48, 19);
			this->newtonDividedForwardPolyLabel->TabIndex = 22;
			this->newtonDividedForwardPolyLabel->Text = L"P(x) = ";
			// 
			// newtonDividedForwardGrid
			// 
			this->newtonDividedForwardGrid->AllowUserToAddRows = false;
			this->newtonDividedForwardGrid->AllowUserToDeleteRows = false;
			this->newtonDividedForwardGrid->AllowUserToResizeRows = false;
			this->newtonDividedForwardGrid->BackgroundColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(250)), 
				static_cast<System::Int32>(static_cast<System::Byte>(250)), static_cast<System::Int32>(static_cast<System::Byte>(250)));
			this->newtonDividedForwardGrid->BorderStyle = System::Windows::Forms::BorderStyle::None;
			this->newtonDividedForwardGrid->CellBorderStyle = System::Windows::Forms::DataGridViewCellBorderStyle::None;
			this->newtonDividedForwardGrid->ColumnHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle13->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle13->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle13->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle13->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle13->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle13->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle13->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->newtonDividedForwardGrid->ColumnHeadersDefaultCellStyle = dataGridViewCellStyle13;
			this->newtonDividedForwardGrid->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			dataGridViewCellStyle14->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle14->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle14->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle14->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(136)), 
				static_cast<System::Int32>(static_cast<System::Byte>(136)), static_cast<System::Int32>(static_cast<System::Byte>(136)));
			dataGridViewCellStyle14->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle14->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle14->WrapMode = System::Windows::Forms::DataGridViewTriState::False;
			this->newtonDividedForwardGrid->DefaultCellStyle = dataGridViewCellStyle14;
			this->newtonDividedForwardGrid->EnableHeadersVisualStyles = false;
			this->newtonDividedForwardGrid->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			this->newtonDividedForwardGrid->GridColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			this->newtonDividedForwardGrid->Location = System::Drawing::Point(3, 9);
			this->newtonDividedForwardGrid->Name = L"newtonDividedForwardGrid";
			this->newtonDividedForwardGrid->ReadOnly = true;
			this->newtonDividedForwardGrid->RowHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle15->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle15->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle15->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle15->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle15->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle15->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle15->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->newtonDividedForwardGrid->RowHeadersDefaultCellStyle = dataGridViewCellStyle15;
			this->newtonDividedForwardGrid->RowHeadersWidthSizeMode = System::Windows::Forms::DataGridViewRowHeadersWidthSizeMode::DisableResizing;
			this->newtonDividedForwardGrid->SelectionMode = System::Windows::Forms::DataGridViewSelectionMode::FullRowSelect;
			this->newtonDividedForwardGrid->Size = System::Drawing::Size(695, 192);
			this->newtonDividedForwardGrid->TabIndex = 21;
			// 
			// LagrangeTab
			// 
			this->LagrangeTab->Controls->Add(this->metroLabel8);
			this->LagrangeTab->Controls->Add(this->lagrangeResult);
			this->LagrangeTab->Controls->Add(this->metroLabel11);
			this->LagrangeTab->Controls->Add(this->lagrangeValue);
			this->LagrangeTab->Controls->Add(this->lagrangePolyLabel);
			this->LagrangeTab->HorizontalScrollbarBarColor = true;
			this->LagrangeTab->HorizontalScrollbarHighlightOnWheel = false;
			this->LagrangeTab->HorizontalScrollbarSize = 10;
			this->LagrangeTab->Location = System::Drawing::Point(4, 38);
			this->LagrangeTab->Name = L"LagrangeTab";
			this->LagrangeTab->Size = System::Drawing::Size(701, 229);
			this->LagrangeTab->TabIndex = 5;
			this->LagrangeTab->Text = L"Lagrange Method";
			this->LagrangeTab->VerticalScrollbarBarColor = true;
			this->LagrangeTab->VerticalScrollbarHighlightOnWheel = false;
			this->LagrangeTab->VerticalScrollbarSize = 10;
			// 
			// metroLabel8
			// 
			this->metroLabel8->AutoSize = true;
			this->metroLabel8->Location = System::Drawing::Point(3, 74);
			this->metroLabel8->Name = L"metroLabel8";
			this->metroLabel8->Size = System::Drawing::Size(59, 19);
			this->metroLabel8->TabIndex = 25;
			this->metroLabel8->Text = L"P(X) => ";
			// 
			// lagrangeResult
			// 
			this->lagrangeResult->AutoSize = true;
			this->lagrangeResult->Location = System::Drawing::Point(68, 74);
			this->lagrangeResult->Name = L"lagrangeResult";
			this->lagrangeResult->Size = System::Drawing::Size(0, 0);
			this->lagrangeResult->TabIndex = 24;
			// 
			// metroLabel11
			// 
			this->metroLabel11->AutoSize = true;
			this->metroLabel11->Location = System::Drawing::Point(3, 39);
			this->metroLabel11->Name = L"metroLabel11";
			this->metroLabel11->Size = System::Drawing::Size(34, 19);
			this->metroLabel11->TabIndex = 23;
			this->metroLabel11->Text = L"X = ";
			// 
			// lagrangeValue
			// 
			this->lagrangeValue->Lines = gcnew cli::array< System::String^  >(1) {L"0"};
			this->lagrangeValue->Location = System::Drawing::Point(43, 39);
			this->lagrangeValue->MaxLength = 32767;
			this->lagrangeValue->Name = L"lagrangeValue";
			this->lagrangeValue->PasswordChar = '\0';
			this->lagrangeValue->ScrollBars = System::Windows::Forms::ScrollBars::None;
			this->lagrangeValue->SelectedText = L"";
			this->lagrangeValue->Size = System::Drawing::Size(75, 23);
			this->lagrangeValue->TabIndex = 22;
			this->lagrangeValue->Text = L"0";
			this->lagrangeValue->UseSelectable = true;
			this->lagrangeValue->TextChanged += gcnew System::EventHandler(this, &Interpolation::lagrangeValue_TextChanged);
			// 
			// lagrangePolyLabel
			// 
			this->lagrangePolyLabel->AutoSize = true;
			this->lagrangePolyLabel->Location = System::Drawing::Point(3, 13);
			this->lagrangePolyLabel->Name = L"lagrangePolyLabel";
			this->lagrangePolyLabel->Size = System::Drawing::Size(48, 19);
			this->lagrangePolyLabel->TabIndex = 21;
			this->lagrangePolyLabel->Text = L"P(x) = ";
			// 
			// generalInterpolationTab
			// 
			this->generalInterpolationTab->Controls->Add(this->metroLabel10);
			this->generalInterpolationTab->Controls->Add(this->generalInterpolationResult);
			this->generalInterpolationTab->Controls->Add(this->metroLabel13);
			this->generalInterpolationTab->Controls->Add(this->generalInterpolationValue);
			this->generalInterpolationTab->Controls->Add(this->generalInterpolationPolyLabel);
			this->generalInterpolationTab->HorizontalScrollbarBarColor = true;
			this->generalInterpolationTab->HorizontalScrollbarHighlightOnWheel = false;
			this->generalInterpolationTab->HorizontalScrollbarSize = 10;
			this->generalInterpolationTab->Location = System::Drawing::Point(4, 38);
			this->generalInterpolationTab->Name = L"generalInterpolationTab";
			this->generalInterpolationTab->Size = System::Drawing::Size(701, 229);
			this->generalInterpolationTab->TabIndex = 6;
			this->generalInterpolationTab->Text = L"General Interpolation Method";
			this->generalInterpolationTab->VerticalScrollbarBarColor = true;
			this->generalInterpolationTab->VerticalScrollbarHighlightOnWheel = false;
			this->generalInterpolationTab->VerticalScrollbarSize = 10;
			// 
			// metroLabel10
			// 
			this->metroLabel10->AutoSize = true;
			this->metroLabel10->Location = System::Drawing::Point(3, 70);
			this->metroLabel10->Name = L"metroLabel10";
			this->metroLabel10->Size = System::Drawing::Size(59, 19);
			this->metroLabel10->TabIndex = 30;
			this->metroLabel10->Text = L"P(X) => ";
			// 
			// generalInterpolationResult
			// 
			this->generalInterpolationResult->AutoSize = true;
			this->generalInterpolationResult->Location = System::Drawing::Point(68, 70);
			this->generalInterpolationResult->Name = L"generalInterpolationResult";
			this->generalInterpolationResult->Size = System::Drawing::Size(0, 0);
			this->generalInterpolationResult->TabIndex = 29;
			// 
			// metroLabel13
			// 
			this->metroLabel13->AutoSize = true;
			this->metroLabel13->Location = System::Drawing::Point(3, 35);
			this->metroLabel13->Name = L"metroLabel13";
			this->metroLabel13->Size = System::Drawing::Size(34, 19);
			this->metroLabel13->TabIndex = 28;
			this->metroLabel13->Text = L"X = ";
			// 
			// generalInterpolationValue
			// 
			this->generalInterpolationValue->Lines = gcnew cli::array< System::String^  >(1) {L"0"};
			this->generalInterpolationValue->Location = System::Drawing::Point(43, 35);
			this->generalInterpolationValue->MaxLength = 32767;
			this->generalInterpolationValue->Name = L"generalInterpolationValue";
			this->generalInterpolationValue->PasswordChar = '\0';
			this->generalInterpolationValue->ScrollBars = System::Windows::Forms::ScrollBars::None;
			this->generalInterpolationValue->SelectedText = L"";
			this->generalInterpolationValue->Size = System::Drawing::Size(75, 23);
			this->generalInterpolationValue->TabIndex = 27;
			this->generalInterpolationValue->Text = L"0";
			this->generalInterpolationValue->UseSelectable = true;
			this->generalInterpolationValue->TextChanged += gcnew System::EventHandler(this, &Interpolation::generalInterpolationValue_TextChanged);
			// 
			// generalInterpolationPolyLabel
			// 
			this->generalInterpolationPolyLabel->AutoSize = true;
			this->generalInterpolationPolyLabel->Location = System::Drawing::Point(3, 9);
			this->generalInterpolationPolyLabel->Name = L"generalInterpolationPolyLabel";
			this->generalInterpolationPolyLabel->Size = System::Drawing::Size(48, 19);
			this->generalInterpolationPolyLabel->TabIndex = 26;
			this->generalInterpolationPolyLabel->Text = L"P(x) = ";
			// 
			// leastSquareTab
			// 
			this->leastSquareTab->BackColor = System::Drawing::Color::White;
			this->leastSquareTab->Controls->Add(this->metroLabel14);
			this->leastSquareTab->Controls->Add(this->leastSquareDegree);
			this->leastSquareTab->Controls->Add(this->metroLabel12);
			this->leastSquareTab->Controls->Add(this->leastSquareResult);
			this->leastSquareTab->Controls->Add(this->metroLabel15);
			this->leastSquareTab->Controls->Add(this->leastSquareValue);
			this->leastSquareTab->Controls->Add(this->leastSquarePolyLabel);
			this->leastSquareTab->Controls->Add(this->leastSquareYk);
			this->leastSquareTab->Controls->Add(this->leastSquareXk);
			this->leastSquareTab->Location = System::Drawing::Point(4, 38);
			this->leastSquareTab->Name = L"leastSquareTab";
			this->leastSquareTab->Size = System::Drawing::Size(701, 229);
			this->leastSquareTab->TabIndex = 7;
			this->leastSquareTab->Text = L"Least Square Method";
			// 
			// metroLabel14
			// 
			this->metroLabel14->AutoSize = true;
			this->metroLabel14->Location = System::Drawing::Point(590, 16);
			this->metroLabel14->Name = L"metroLabel14";
			this->metroLabel14->Size = System::Drawing::Size(52, 19);
			this->metroLabel14->TabIndex = 33;
			this->metroLabel14->Text = L"Degree";
			// 
			// leastSquareDegree
			// 
			this->leastSquareDegree->BorderStyle = System::Windows::Forms::BorderStyle::None;
			this->leastSquareDegree->Location = System::Drawing::Point(592, 38);
			this->leastSquareDegree->Minimum = System::Decimal(gcnew cli::array< System::Int32 >(4) {2, 0, 0, 0});
			this->leastSquareDegree->Name = L"leastSquareDegree";
			this->leastSquareDegree->Size = System::Drawing::Size(85, 16);
			this->leastSquareDegree->TabIndex = 32;
			this->leastSquareDegree->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) {2, 0, 0, 0});
			this->leastSquareDegree->ValueChanged += gcnew System::EventHandler(this, &Interpolation::leastSquareDegree_ValueChanged);
			// 
			// metroLabel12
			// 
			this->metroLabel12->AutoSize = true;
			this->metroLabel12->Location = System::Drawing::Point(141, 243);
			this->metroLabel12->Name = L"metroLabel12";
			this->metroLabel12->Size = System::Drawing::Size(59, 19);
			this->metroLabel12->TabIndex = 31;
			this->metroLabel12->Text = L"P(X) => ";
			// 
			// leastSquareResult
			// 
			this->leastSquareResult->AutoSize = true;
			this->leastSquareResult->Location = System::Drawing::Point(206, 248);
			this->leastSquareResult->Name = L"leastSquareResult";
			this->leastSquareResult->Size = System::Drawing::Size(0, 0);
			this->leastSquareResult->TabIndex = 30;
			// 
			// metroLabel15
			// 
			this->metroLabel15->AutoSize = true;
			this->metroLabel15->Location = System::Drawing::Point(20, 243);
			this->metroLabel15->Name = L"metroLabel15";
			this->metroLabel15->Size = System::Drawing::Size(34, 19);
			this->metroLabel15->TabIndex = 29;
			this->metroLabel15->Text = L"X = ";
			// 
			// leastSquareValue
			// 
			this->leastSquareValue->Lines = gcnew cli::array< System::String^  >(1) {L"0"};
			this->leastSquareValue->Location = System::Drawing::Point(60, 243);
			this->leastSquareValue->MaxLength = 32767;
			this->leastSquareValue->Name = L"leastSquareValue";
			this->leastSquareValue->PasswordChar = '\0';
			this->leastSquareValue->ScrollBars = System::Windows::Forms::ScrollBars::None;
			this->leastSquareValue->SelectedText = L"";
			this->leastSquareValue->Size = System::Drawing::Size(75, 23);
			this->leastSquareValue->TabIndex = 28;
			this->leastSquareValue->Text = L"0";
			this->leastSquareValue->UseSelectable = true;
			this->leastSquareValue->TextChanged += gcnew System::EventHandler(this, &Interpolation::leastSquareValue_TextChanged);
			// 
			// leastSquarePolyLabel
			// 
			this->leastSquarePolyLabel->AutoSize = true;
			this->leastSquarePolyLabel->Location = System::Drawing::Point(20, 217);
			this->leastSquarePolyLabel->Name = L"leastSquarePolyLabel";
			this->leastSquarePolyLabel->Size = System::Drawing::Size(48, 19);
			this->leastSquarePolyLabel->TabIndex = 27;
			this->leastSquarePolyLabel->Text = L"P(x) = ";
			// 
			// leastSquareYk
			// 
			this->leastSquareYk->AllowUserToAddRows = false;
			this->leastSquareYk->AllowUserToDeleteRows = false;
			this->leastSquareYk->AllowUserToResizeRows = false;
			this->leastSquareYk->BackgroundColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(250)), 
				static_cast<System::Int32>(static_cast<System::Byte>(250)), static_cast<System::Int32>(static_cast<System::Byte>(250)));
			this->leastSquareYk->BorderStyle = System::Windows::Forms::BorderStyle::None;
			this->leastSquareYk->CellBorderStyle = System::Windows::Forms::DataGridViewCellBorderStyle::None;
			this->leastSquareYk->ColumnHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle16->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle16->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle16->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle16->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle16->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle16->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle16->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->leastSquareYk->ColumnHeadersDefaultCellStyle = dataGridViewCellStyle16;
			this->leastSquareYk->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			dataGridViewCellStyle17->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle17->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle17->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle17->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(136)), 
				static_cast<System::Int32>(static_cast<System::Byte>(136)), static_cast<System::Int32>(static_cast<System::Byte>(136)));
			dataGridViewCellStyle17->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle17->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle17->WrapMode = System::Windows::Forms::DataGridViewTriState::False;
			this->leastSquareYk->DefaultCellStyle = dataGridViewCellStyle17;
			this->leastSquareYk->EnableHeadersVisualStyles = false;
			this->leastSquareYk->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			this->leastSquareYk->GridColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)));
			this->leastSquareYk->Location = System::Drawing::Point(304, 16);
			this->leastSquareYk->Name = L"leastSquareYk";
			this->leastSquareYk->ReadOnly = true;
			this->leastSquareYk->RowHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle18->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle18->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle18->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle18->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle18->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle18->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle18->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->leastSquareYk->RowHeadersDefaultCellStyle = dataGridViewCellStyle18;
			this->leastSquareYk->RowHeadersWidthSizeMode = System::Windows::Forms::DataGridViewRowHeadersWidthSizeMode::DisableResizing;
			this->leastSquareYk->SelectionMode = System::Windows::Forms::DataGridViewSelectionMode::FullRowSelect;
			this->leastSquareYk->Size = System::Drawing::Size(276, 192);
			this->leastSquareYk->TabIndex = 23;
			// 
			// leastSquareXk
			// 
			this->leastSquareXk->AllowUserToAddRows = false;
			this->leastSquareXk->AllowUserToDeleteRows = false;
			this->leastSquareXk->AllowUserToResizeRows = false;
			this->leastSquareXk->BackgroundColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(250)), 
				static_cast<System::Int32>(static_cast<System::Byte>(250)), static_cast<System::Int32>(static_cast<System::Byte>(250)));
			this->leastSquareXk->BorderStyle = System::Windows::Forms::BorderStyle::None;
			this->leastSquareXk->CellBorderStyle = System::Windows::Forms::DataGridViewCellBorderStyle::None;
			this->leastSquareXk->ColumnHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle19->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle19->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle19->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle19->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle19->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle19->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle19->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->leastSquareXk->ColumnHeadersDefaultCellStyle = dataGridViewCellStyle19;
			this->leastSquareXk->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			dataGridViewCellStyle20->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle20->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle20->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle20->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(136)), 
				static_cast<System::Int32>(static_cast<System::Byte>(136)), static_cast<System::Int32>(static_cast<System::Byte>(136)));
			dataGridViewCellStyle20->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle20->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle20->WrapMode = System::Windows::Forms::DataGridViewTriState::False;
			this->leastSquareXk->DefaultCellStyle = dataGridViewCellStyle20;
			this->leastSquareXk->EnableHeadersVisualStyles = false;
			this->leastSquareXk->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			this->leastSquareXk->GridColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)));
			this->leastSquareXk->Location = System::Drawing::Point(14, 16);
			this->leastSquareXk->Name = L"leastSquareXk";
			this->leastSquareXk->ReadOnly = true;
			this->leastSquareXk->RowHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle21->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle21->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle21->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle21->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle21->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle21->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle21->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->leastSquareXk->RowHeadersDefaultCellStyle = dataGridViewCellStyle21;
			this->leastSquareXk->RowHeadersWidthSizeMode = System::Windows::Forms::DataGridViewRowHeadersWidthSizeMode::DisableResizing;
			this->leastSquareXk->SelectionMode = System::Windows::Forms::DataGridViewSelectionMode::FullRowSelect;
			this->leastSquareXk->Size = System::Drawing::Size(276, 192);
			this->leastSquareXk->TabIndex = 22;
			// 
			// subLineTab
			// 
			this->subLineTab->BackColor = System::Drawing::Color::White;
			this->subLineTab->Controls->Add(this->metroLabel16);
			this->subLineTab->Controls->Add(this->subLineResult);
			this->subLineTab->Controls->Add(this->metroLabel18);
			this->subLineTab->Controls->Add(this->subLineValue);
			this->subLineTab->Controls->Add(this->subLineGrid);
			this->subLineTab->Location = System::Drawing::Point(4, 38);
			this->subLineTab->Name = L"subLineTab";
			this->subLineTab->Size = System::Drawing::Size(701, 229);
			this->subLineTab->TabIndex = 8;
			this->subLineTab->Text = L"Spline Method(Linear)";
			// 
			// metroLabel16
			// 
			this->metroLabel16->AutoSize = true;
			this->metroLabel16->Location = System::Drawing::Point(21, 239);
			this->metroLabel16->Name = L"metroLabel16";
			this->metroLabel16->Size = System::Drawing::Size(59, 19);
			this->metroLabel16->TabIndex = 35;
			this->metroLabel16->Text = L"P(X) => ";
			// 
			// subLineResult
			// 
			this->subLineResult->AutoSize = true;
			this->subLineResult->Location = System::Drawing::Point(86, 244);
			this->subLineResult->Name = L"subLineResult";
			this->subLineResult->Size = System::Drawing::Size(0, 0);
			this->subLineResult->TabIndex = 34;
			// 
			// metroLabel18
			// 
			this->metroLabel18->AutoSize = true;
			this->metroLabel18->Location = System::Drawing::Point(18, 209);
			this->metroLabel18->Name = L"metroLabel18";
			this->metroLabel18->Size = System::Drawing::Size(34, 19);
			this->metroLabel18->TabIndex = 33;
			this->metroLabel18->Text = L"X = ";
			// 
			// subLineValue
			// 
			this->subLineValue->Lines = gcnew cli::array< System::String^  >(1) {L"0"};
			this->subLineValue->Location = System::Drawing::Point(58, 209);
			this->subLineValue->MaxLength = 32767;
			this->subLineValue->Name = L"subLineValue";
			this->subLineValue->PasswordChar = '\0';
			this->subLineValue->ScrollBars = System::Windows::Forms::ScrollBars::None;
			this->subLineValue->SelectedText = L"";
			this->subLineValue->Size = System::Drawing::Size(75, 23);
			this->subLineValue->TabIndex = 32;
			this->subLineValue->Text = L"0";
			this->subLineValue->UseSelectable = true;
			this->subLineValue->TextChanged += gcnew System::EventHandler(this, &Interpolation::subLineValue_TextChanged);
			// 
			// subLineGrid
			// 
			this->subLineGrid->AllowUserToAddRows = false;
			this->subLineGrid->AllowUserToDeleteRows = false;
			this->subLineGrid->AllowUserToResizeRows = false;
			this->subLineGrid->BackgroundColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(250)), 
				static_cast<System::Int32>(static_cast<System::Byte>(250)), static_cast<System::Int32>(static_cast<System::Byte>(250)));
			this->subLineGrid->BorderStyle = System::Windows::Forms::BorderStyle::None;
			this->subLineGrid->CellBorderStyle = System::Windows::Forms::DataGridViewCellBorderStyle::None;
			this->subLineGrid->ColumnHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle22->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle22->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle22->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle22->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle22->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle22->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle22->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->subLineGrid->ColumnHeadersDefaultCellStyle = dataGridViewCellStyle22;
			this->subLineGrid->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->subLineGrid->Columns->AddRange(gcnew cli::array< System::Windows::Forms::DataGridViewColumn^  >(1) {this->column1});
			dataGridViewCellStyle23->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle23->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle23->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle23->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(136)), 
				static_cast<System::Int32>(static_cast<System::Byte>(136)), static_cast<System::Int32>(static_cast<System::Byte>(136)));
			dataGridViewCellStyle23->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle23->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle23->WrapMode = System::Windows::Forms::DataGridViewTriState::False;
			this->subLineGrid->DefaultCellStyle = dataGridViewCellStyle23;
			this->subLineGrid->EnableHeadersVisualStyles = false;
			this->subLineGrid->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			this->subLineGrid->GridColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)));
			this->subLineGrid->Location = System::Drawing::Point(14, 12);
			this->subLineGrid->Name = L"subLineGrid";
			this->subLineGrid->ReadOnly = true;
			this->subLineGrid->RowHeadersBorderStyle = System::Windows::Forms::DataGridViewHeaderBorderStyle::None;
			dataGridViewCellStyle24->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle24->BackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(174)), 
				static_cast<System::Int32>(static_cast<System::Byte>(219)));
			dataGridViewCellStyle24->Font = (gcnew System::Drawing::Font(L"Segoe UI", 11, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Pixel));
			dataGridViewCellStyle24->ForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(255)), 
				static_cast<System::Int32>(static_cast<System::Byte>(255)), static_cast<System::Int32>(static_cast<System::Byte>(255)));
			dataGridViewCellStyle24->SelectionBackColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), 
				static_cast<System::Int32>(static_cast<System::Byte>(198)), static_cast<System::Int32>(static_cast<System::Byte>(247)));
			dataGridViewCellStyle24->SelectionForeColor = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(17)), 
				static_cast<System::Int32>(static_cast<System::Byte>(17)), static_cast<System::Int32>(static_cast<System::Byte>(17)));
			dataGridViewCellStyle24->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->subLineGrid->RowHeadersDefaultCellStyle = dataGridViewCellStyle24;
			this->subLineGrid->RowHeadersWidthSizeMode = System::Windows::Forms::DataGridViewRowHeadersWidthSizeMode::DisableResizing;
			this->subLineGrid->SelectionMode = System::Windows::Forms::DataGridViewSelectionMode::FullRowSelect;
			this->subLineGrid->Size = System::Drawing::Size(680, 192);
			this->subLineGrid->TabIndex = 24;
			// 
			// column1
			// 
			this->column1->HeaderText = L"Si(x)";
			this->column1->MinimumWidth = 620;
			this->column1->Name = L"column1";
			this->column1->ReadOnly = true;
			this->column1->Width = 620;
			// 
			// metroLabel17
			// 
			this->metroLabel17->AutoSize = true;
			this->metroLabel17->Location = System::Drawing::Point(279, 36);
			this->metroLabel17->Name = L"metroLabel17";
			this->metroLabel17->Size = System::Drawing::Size(118, 19);
			this->metroLabel17->TabIndex = 2;
			this->metroLabel17->Text = L"Numerical Analysis";
			this->metroLabel17->Click += gcnew System::EventHandler(this, &Interpolation::metroLabel17_Click);
			// 
			// metroLabel19
			// 
			this->metroLabel19->AutoSize = true;
			this->metroLabel19->Location = System::Drawing::Point(279, 94);
			this->metroLabel19->Name = L"metroLabel19";
			this->metroLabel19->Size = System::Drawing::Size(118, 19);
			this->metroLabel19->TabIndex = 3;
			this->metroLabel19->Text = L"Numerical Analysis";
			// 
			// panel
			// 
			this->panel->HorizontalScrollbarBarColor = true;
			this->panel->HorizontalScrollbarHighlightOnWheel = false;
			this->panel->HorizontalScrollbarSize = 10;
			this->panel->Location = System::Drawing::Point(24, 380);
			this->panel->Name = L"panel";
			this->panel->Size = System::Drawing::Size(911, 160);
			this->panel->TabIndex = 2;
			this->panel->VerticalScrollbarBarColor = true;
			this->panel->VerticalScrollbarHighlightOnWheel = false;
			this->panel->VerticalScrollbarSize = 10;
			this->panel->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &Interpolation::panel_Paint);
			// 
			// Interpolation
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(953, 557);
			this->Controls->Add(this->panel);
			this->Controls->Add(this->tabControl);
			this->Controls->Add(this->pointsCountUpDown);
			this->Controls->Add(this->pointsGrid);
			this->MaximizeBox = false;
			this->Name = L"Interpolation";
			this->Text = L"Interpolation";
			this->Load += gcnew System::EventHandler(this, &Interpolation::Interpolation_Load);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->pointsGrid))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->pointsCountUpDown))->EndInit();
			this->tabControl->ResumeLayout(false);
			this->helpTab->ResumeLayout(false);
			this->helpTab->PerformLayout();
			this->nuetonProgressiveTab->ResumeLayout(false);
			this->nuetonProgressiveTab->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->nuetonProgressiveGrid))->EndInit();
			this->newtonBackword->ResumeLayout(false);
			this->newtonBackword->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->newtonBackwordGrid))->EndInit();
			this->NewtonDividedDifferencebackward->ResumeLayout(false);
			this->NewtonDividedDifferencebackward->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->newtonDividedBackwordGrid))->EndInit();
			this->NewtonDividedDifferenceforward->ResumeLayout(false);
			this->NewtonDividedDifferenceforward->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->newtonDividedForwardGrid))->EndInit();
			this->LagrangeTab->ResumeLayout(false);
			this->LagrangeTab->PerformLayout();
			this->generalInterpolationTab->ResumeLayout(false);
			this->generalInterpolationTab->PerformLayout();
			this->leastSquareTab->ResumeLayout(false);
			this->leastSquareTab->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->leastSquareDegree))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->leastSquareYk))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->leastSquareXk))->EndInit();
			this->subLineTab->ResumeLayout(false);
			this->subLineTab->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->subLineGrid))->EndInit();
			this->ResumeLayout(false);

		}

#pragma endregion
	private: System::Void changeRows(MetroGrid^ grid , int x){


				 this->pointsGrid->CellValueChanged -= gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &Interpolation::pointsGrid_CellValueChanged);

				 int d = grid->Rows->Count ; 

				 List<Tuple<String^ , String^ >^ > ^v = gcnew List<Tuple<String^ , String^ >^ > (d); 

				 for(int i =0 ; i < d ; i++) { 
					 v->Add( gcnew Tuple<String^ , String ^> (  grid->Rows[i]->Cells[0]->Value->ToString() , grid->Rows[i]->Cells[1]->Value->ToString()));
				 }
				 grid->Rows->Clear(); 

				 for(int i =0 ; i < x ; i ++) { 
					 grid->Rows->Add(); 

					 grid->Rows[i]->Cells[0]->Value = ( i < d ?    v[i]->Item1 : "0"); 
					 grid->Rows[i]->Cells[1]->Value = ( i < d ?    v[i]->Item2 : "0"); 

				 }
				 this->pointsGrid->CellValueChanged += gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &Interpolation::pointsGrid_CellValueChanged);
				 makeChange();
			 }
	private: System::Void metroButton1_Click(System::Object^  sender, System::EventArgs^  e) {
			 }
	private: System::Void Interpolation_Load(System::Object^  sender, System::EventArgs^  e) {
				 this->changeRows(this->pointsGrid , Convert::ToInt32(pointsCountUpDown->Value));
			 }

	private: System::Void pointsCountUpDown_ValueChanged(System::Object^  sender, System::EventArgs^  e) {
				 this->changeRows(this->pointsGrid , Convert::ToInt32( pointsCountUpDown->Value)) ; 

			 }

			 double tryConvertDouble(String^ s){

				 double  f= 0 ; 
				 try{
					 f= Convert::ToDouble(s) ;
				 }
				 catch(...){
					 f= 0 ; 
				 }
				 return f; 
			 }
	private: System::Void pointsGrid_CellValueChanged(System::Object^  sender, System::Windows::Forms::DataGridViewCellEventArgs^  e) {

				 int x = e->RowIndex ; 
				 int y = e->ColumnIndex ; 
				 if(x <= -1 || y <= -1 ) return ;
				 if(x >= pointsGrid->Rows->Count || y >= pointsGrid->Columns->Count ) return ; 
				 double f = 0 ; 
				 try{
					 f = tryConvertDouble(pointsGrid[y,x]->Value->ToString());
				 }catch(...){
					 f =0 ; 
				 }
				 pointsGrid[y,x]->Value = f; 
				 makeChange();




			 }
			 void fillGridFromVector(MetroGrid^ grid , vector<vector<ld>> v , String^ header ,vector< pair<ld ,ld > > & points ) {

				 grid->Rows->Clear();
				 grid->Columns->Clear();


				 grid->Columns->Add( "X" , "X") ; 

				 for(int i =0 ; i < v.size() ; i++){
					 grid->Columns->Add( i.ToString() , i ?  header + Convert::ToString(i) : "Y" ) ; 
				 }
				 for(int i =0 ; i < v.size() ; i++){
					 grid->Rows->Add();
					 for(int j =0 ; j < v[i].size() ; j++){
						 grid[j+1,i]->Value = v[i][j]; 
					 }
				 }			 
				 for(int i =0 ; i < v.size() ; i++){
					 grid[0,i]->Value = points[i].first ; 
				 }


			 }
			 void fillGridFromVector(MetroGrid^ grid , vector<vector<ld>> v , String^ header ) {

				 grid->Rows->Clear();
				 grid->Columns->Clear();


				 for(int i =0 ; i < v[0].size() ; i++){
					 grid->Columns->Add( i.ToString() ,   header + Convert::ToString(i)  ) ; 
				 }
				 for(int i =0 ; i < v.size() ; i++){
					 grid->Rows->Add();
					 for(int j =0 ; j < v[i].size() ; j++){
						 grid[j,i]->Value = v[i][j]; 
					 }
				 }			 


			 }

			 void makeChange(){

				 int x =tabControl->SelectedIndex ; 
				 switch (x)
				 {
				 case 0 : 
					 return ; 
					 break;
				 case 1:
					 fillNeutonProgressive();
					 break;
				 case 2: 
					 fillNewtonBackword();
					 break;
				 case 3:
					 fillNewtonDividedBackword();
					 break;
				 case 4 :
					 fillNewtonDividedForword();
					 break;
				 case 5:
					 fillLagrange();
					 break;
				 case 6 :
					 fillGeneral();
					 break;
				 case 7:
					 fillLeastSquare();
					 break;
				 case 8:
					 fillSubLine();
				break;
				 default:
					 break;
				 }

			 }
			 String^ stringToString(string& s){
				 String^ res = "" ; 
				 for(int i =0 ; i < s.length() ; i++){
					 res += Convert::ToChar( s[i]);
				 }
				 return res; 
			 }

			 void changeLabelPoly(MetroLabel^ l , Poly& p , double x){
				 l->Text = p.value(x).ToString();
			 }


			 void fillNewtonBackword(){

				 vector<pair<ld , ld > > points = gridToPonts(this->pointsGrid);
				 pair<Poly,vector<vector<ld>> >  res = Newton_backward(points);
				 if(!is_const(points)){
					 newtonBackwordGrid->Visible = false ; 
					 newtonBackwordPolyLabel2->Visible = false; 				 
					 metroLabel5->Visible = false ;
					 metroLabel2->Visible = false ; 
					 newtonBackwordValue->Visible = false; 

					 newtonBackwordResult->Visible = false ; 

					 newtonBackwordFail->Visible= true; 

					 return ; 

				 }else {

					 newtonBackwordGrid->Visible =! false ; 
					 newtonBackwordPolyLabel2->Visible =! false; 				 
					 metroLabel5->Visible =! false ;
					 newtonBackwordValue->Visible =! false; 
					 metroLabel2->Visible = true;
					 newtonBackwordResult->Visible = !false ; 

					 newtonBackwordFail->Visible= !true; 


				 }
				 newtonBackwordPolyLabel2->Text  = stringToString(res.first.to_string());
				 fillGridFromVector(this->newtonBackwordGrid , res.second , gcnew String("Delta") , points) ; 
				 changeLabelPoly(newtonBackwordResult , res.first ,  tryConvertDouble(newtonBackwordValue->Text));
				 draw(res.first);
			 }
			 void fillNewtonDividedBackword(){

				 vector<pair<ld , ld > > points = gridToPonts(this->pointsGrid);
				 pair<Poly,vector<vector<ld>> >  res = Newton_Divided_Difference_backward(points);

				 newtonDividedBackwordPolyLabel->Text  = stringToString(res.first.to_string());
				 fillGridFromVector(this->newtonDividedBackwordGrid , res.second , gcnew String("Delta") , points) ; 
				 changeLabelPoly(newtonDividedBackwordResult , res.first ,  tryConvertDouble(newtonDividedBackwordValue->Text));
				 draw(res.first);

			 }
			 void fillLeastSquare(){
				 vector<pair<ld , ld > > points = gridToPonts(this->pointsGrid);

				 auto res = squares(points , (int)leastSquareDegree->Value);
				 leastSquarePolyLabel->Text  = stringToString(res.first.to_string());

				 fillGridFromVector(this->leastSquareXk , res.second.first , gcnew String("Xk^") ) ;
				 fillGridFromVector(this->leastSquareYk , res.second.second , gcnew String("Yk*Xk^") ) ;
				 changeLabelPoly(leastSquareResult , res.first ,  tryConvertDouble(leastSquareValue->Text));
				 draw(res.first);
			 }
			 void fillNewtonDividedForword(){

				 vector<pair<ld , ld > > points = gridToPonts(this->pointsGrid);
				 pair<Poly,vector<vector<ld>> >  res = Newton_Divided_Difference_forward(points);

				 newtonDividedForwardPolyLabel->Text  = stringToString(res.first.to_string());
				 fillGridFromVector(this->newtonDividedForwardGrid , res.second , gcnew String("Delta") , points) ; 
				 changeLabelPoly(newtonDividedForwardResult , res.first ,  tryConvertDouble(newtonDividedForwardValue->Text));
				 draw(res.first);

			 }
			 void fillGeneral(){
				 vector<pair<ld , ld > > points = gridToPonts(this->pointsGrid);
				 Poly res = general_interpolation(points);
				 generalInterpolationPolyLabel->Text = stringToString(res.to_string());
				 changeLabelPoly(generalInterpolationResult , res ,  tryConvertDouble(generalInterpolationValue->Text));
				 draw(res);
			 }
			 void fillLagrange(){
				 vector<pair<ld , ld > > points = gridToPonts(this->pointsGrid);
				 Poly res = Lagrange(points);
				 lagrangePolyLabel->Text = stringToString(res.to_string());
				 changeLabelPoly(lagrangeResult , res ,  tryConvertDouble(lagrangeValue->Text));
				 draw(res);
			 }
			 void fillSubLine(){
				  vector<pair<ld , ld > > points = gridToPonts(this->pointsGrid);
				  spline_solution res = Spline(points);
				  subLineGrid->Rows->Clear();
				  for(int i = 0 ; i < (int) res.v.size() ; i++){
					  subLineGrid->Rows->Add();
					  subLineGrid->Rows[i]->Cells[0]->Value = stringToString(res[i].to_string());
				  }

				  try{
					  subLineResult->Text =   res.evaluate(tryConvertDouble(subLineValue->Text)).ToString();
				  }catch(exception e){
					  subLineResult->Text = stringToString(string(e.what()));
				  }
				  
				 

			 }
			 void draw(Poly& P ){

				 
				 ld w = this->panel->Width; 
				 ld h = this->panel->Height; 
				 Graphics^e =  panel->CreateGraphics() ; 
				 e->Clear(Color::White);
				 vector<pair<long double, long double > > res = P.Draws(w,h ,20 , 1);

				 try{
				 for (int i = 1; i < res.size(); ++i)
				 {
					 cout << res[i].first << " " << res[i].second << endl; 
					 e->DrawLine(System::Drawing::Pens::Black, (float)(res[i - 1].first + w / 2), -2 * (float)res[i - 1].second + (float)h / 2, (float)(res[i].first + (float)w / 2), -2 * (float)res[i].second + (float)h / 2);
				 }
				 }
				 catch(...){

				 }

			 }
			 void fillNeutonProgressive(){
				 vector< pair<ld ,ld> > v ; 
				 vector<pair<ld , ld > > points = gridToPonts(this->pointsGrid);
				 pair<Poly,vector<vector<ld>> >  res = Newton_forward(points);
				 if(!is_const(points)){

					 nuetonProgressiveGrid->Visible = false; 
					 nuetonProgreesivePolyLabel->Visible = false; 
					 metroLabel1->Visible = false ;
					 neutonProgressiveValue->Visible = false; 
					 metroLabel3->Visible = false; 
					 neutonProgressiveResult->Visible = false ; 

					 neutoProgressiveFail->Visible = true;  
					 return ; 

				 }else {

					 nuetonProgressiveGrid->Visible = !false; 
					 nuetonProgreesivePolyLabel->Visible = !false; 
					 metroLabel1->Visible = !false ;
					 neutonProgressiveValue->Visible = !false; 
					 metroLabel3->Visible = !false; 
					 neutonProgressiveResult->Visible = !false ; 


					 neutoProgressiveFail->Visible = false ;

				 }
				 draw(res.first);
				 nuetonProgreesivePolyLabel->Text  = stringToString(res.first.to_string());
				 fillGridFromVector(this->nuetonProgressiveGrid , res.second , gcnew String("Delta") , points) ; 
				 changeLabelPoly(neutonProgressiveResult , res.first ,  tryConvertDouble(neutonProgressiveValue->Text));
			 }

			 vector< pair<ld ,ld> > gridToPonts(MetroGrid^ grid){

				 vector< pair<ld ,ld> > res ;
				 for(int i =0 ; i < grid->Rows->Count ; i++){
					 res.push_back(make_pair(tryConvertDouble(grid->Rows[i]->Cells[0]->Value->ToString() ), tryConvertDouble(grid->Rows[i]->Cells[1]->Value->ToString() )));
				 }
				 sort(res.begin() , res.end()) ; 
				 return res;

			 }
	private: System::Void metroTabControl1_SelectedIndexChanged(System::Object^  sender, System::EventArgs^  e) {
				 makeChange();
			 }



	private: System::Void nuetonProgressiveTab_Click(System::Object^  sender, System::EventArgs^  e) {

			 }
	private: System::Void neutonProgressiveValue_Click(System::Object^  sender, System::EventArgs^  e) {
			 }

	private: System::Void neutonProgressiveValue_TextChanged(System::Object^  sender, System::EventArgs^  e) {

				 makeChange();

			 }
	private: System::Void T(System::Object^  sender, System::EventArgs^  e) {
			 }
	private: System::Void newtonBackwordValue_TextChanged(System::Object^  sender, System::EventArgs^  e) {
				 makeChange();
			 }
	private: System::Void newtonDividedBackwordValue_TextChanged(System::Object^  sender, System::EventArgs^  e) {
				 makeChange();
			 }
	private: System::Void newtonDividedForwardValue_TextChanged(System::Object^  sender, System::EventArgs^  e) {
				 makeChange();
			 }
	private: System::Void generalInterpolationValue_TextChanged(System::Object^  sender, System::EventArgs^  e) {
				 makeChange();
			 }
	private: System::Void lagrangeValue_TextChanged(System::Object^  sender, System::EventArgs^  e) {
				 makeChange();
			 }
	private: System::Void leastSquareDegree_ValueChanged(System::Object^  sender, System::EventArgs^  e) {
				 makeChange();
			 }
	private: System::Void leastSquareValue_TextChanged(System::Object^  sender, System::EventArgs^  e) {
				 makeChange();
			 }
	private: System::Void subLineValue_TextChanged(System::Object^  sender, System::EventArgs^  e) {
				 makeChange();
			 }
private: System::Void metroLabel17_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void panel_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e) {

		 }
};
}
