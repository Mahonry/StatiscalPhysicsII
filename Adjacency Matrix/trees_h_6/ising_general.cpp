#include<cstdlib>
#include<fstream> 
#include<iostream>
#include<vector>
#include<random>
#include<cmath>
using namespace std;
#define Nv 1093


int RandomChoice(vector<double> Freqs,mt19937 &gen)
{

	discrete_distribution<> d(Freqs.begin(),Freqs.end());
	return d(gen);	
}

double Magnetization_spin(int State[Nv])
{
	double m_spin = 0;
	for(int x=0; x<Nv;x++)
	{
		m_spin+=State[x];
	}
	return (1.0*m_spin)/(1.0*Nv);
}

double Energy_spin(int State[Nv], vector<vector<int>> Neighbords)
{
	double energy = 0;
	int count = 0;

	for(int x=0; x<Nv;x++)
	{
		for(int j = 0; j<Neighbords[count].size(); j++)
		{
			energy += -State[x]*State[Neighbords[count][j]];			
		}
		count += 1;
	}

	return 1.0*energy/(2.0*Nv);
}

void Flips(double T, int delta_t, int State[Nv],  vector<vector<int>> Neighbords, uniform_int_distribution<> &S_choice,mt19937 &gen)
{
	for(int t = 0; t<delta_t;t++)
	{
		int k = S_choice(gen);

		int DeltaE = 0;

		for(int j = 0; j < Neighbords[k].size(); j ++)
		{
			DeltaE += 2*State[k]*State[Neighbords[k][j]];
		}


		if(DeltaE>0)
		{
			double expval = exp(-DeltaE/T);
			State[k] = (2*RandomChoice({expval, 1 - expval},gen)-1)*State[k];
			
		}
		else
		{
			State[k] = -1*State[k];
		}
	}
}










int main()
{
	random_device rd;
	mt19937 gen(rd());
	vector<vector<int>> Neighbords(Nv);
	vector<double> Mspin;
	int State[Nv];
	int count = 0;	
	int k;
	uniform_int_distribution<> S_choice(0, Nv - 1);
	double T = 5.0;
	int a, b;
	string str;

	
	//Cargamos la red 
	ifstream file("Lines.dat");

	while(getline(file,str))	
	{
		a = stoi(str.substr(0,str.find("  ")));
		b = stoi(str.substr(str.find("  ")+2));		
		Neighbords[a].push_back(b);
		Neighbords[b].push_back(a);
	}

	//Inicializamos los spines
	for(int i = 0; i <Nv; i++)
	{
		State[i] = 2*RandomChoice({0.5,0.5}, gen)-1;
	}
	

	//Thermalize
	cout<<"Termalizando"<<endl;
	Flips(T,1000000*Nv, State, Neighbords, S_choice, gen);	




//Metropolis
	
	int realizations=10000;
	ofstream Results ("Results_ensemble.dat");
	
	while(T>0.1)
	{
		string name = "./Results/T_" + to_string(T).substr(0,3) + ".dat";
		ofstream State_T (name);
		for(int n = 0; n<realizations; n++)
		{
			
			Flips(T, 100*Nv, State,  Neighbords, S_choice, gen);
			Results<<T<<" "<<Magnetization_spin(State)<<" "<<Energy_spin(State, Neighbords)<<endl;
			for (int i = 0; i < Nv; i++)
			{
				State_T<<State[i]<<" ";
			}
			State_T<<endl;
		}
		State_T.close();
		T-=0.1;
		cout<<" For "<<T<<endl;
	}	

	Results.close();
	return 0;
}
