//###########################################//
//### Hoping Mobility Kinetic Monte-Carlo ###//
//###########################################//

#include "DMS.h"
#include <iomanip>
#include <omp.h>
#include <string>
#include <math.h>

using namespace std;
using namespace std::chrono;
using namespace DMS;

int    next_hop(std::vector<std::pair<int, double>> Hopping_Probs);				// determines next hop < index, rate >
void   relaxation(DisorderedSystem & DS, const int Hopps);						// relxation procedure
double mobility_calculation(DisorderedSystem & DS, const int Hopps);			// hopping procedure to determine mobility

int main()
{

	// Inputs
	const int    relax_hop_num  = 1e6;
	const int    mob_hop_num    = 1e7;
	const double N				= 1e4;
	const double _E0			= 4;
	const double _Ef			= -4;
	const double _F_i			= 0.01;
	const double _F_f			= 3.5;
	const double _dF			= 0.5;
	const double _traps_conc	= 2.7e-2;
	const double cutoff			= 5 / cbrt(_traps_conc);						 // [ Avg. Dist. units ]

	const int	 NMC			= 32;
	// ------

	cout << " ############################################################## " << endl;
	cout << " ### --- Hoping Mobility Kinetic Monte-Carlo Simulation --- ### " << endl;
	cout << " ###  (c) Vitalii Valkovskii, University of Marburg, 2017   ### " << endl;
	cout << " ############################################################## " << endl;
	cout << endl;
	
	ofstream output_file("Mobility.dat");
	
	cout << " Simulation in progress ... " << endl;
	auto t1 = high_resolution_clock::now();
	for (double _F = _F_i; _F <= _F_f; _F += _dF)
	{
		double mobility = 0;
		double exp_factor;
#pragma omp parallel for default(shared) reduction(+:mobility) schedule(dynamic)  
		for (int n = 0; n < NMC; n++)
		{
			std::random_device rd;
			std::mt19937 generator(rd());
			std::uniform_int_distribution <int> distribution(0, N);
			
			DisorderedSystem DS(N, _E0, _traps_conc, cutoff, _F); 
			DS.Calc_Effective_T(0.67);
			DS.Set_dls_E_Fermi(_Ef);
			DS.Calc_Transport_Energy_ET(0.01);
			DS.Set_Fermi_Occupation_Probs_ET();
			DS.Precalculate_Hopping_Rates_NA_BC();
			DS.Precalculate_Hopping_Probs();
			DS.Add_Carrier(0, distribution(generator), 0, 1);
			relaxation(DS, relax_hop_num);
			mobility += mobility_calculation(DS, mob_hop_num);
			exp_factor = (DS.Get_dls_E_Transport() - _Ef) / DS.Get_T_Eff_Factor();
		}
		
		output_file << _F << '	' << exp_factor <<'	' << mobility / NMC << endl;	   
	}
	output_file.close();
	auto t2 = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(t2 - t1).count();
	cout << " Simulation is finished. " << endl << " Time: " << duration << " seconds" << endl;
	cin.get();
	return 0;
}

int		next_hop(std::vector<std::pair<int, double>> Hopping_Probs)
{
	random_device rd1;
	mt19937 generator1(rd1());
	uniform_real_distribution <double> distribution(0, 1);
	double p = distribution(generator1);
	int i = 0;
	while (Hopping_Probs[i].second < p) i++;
	return Hopping_Probs[i].first;
}
void	relaxation(DisorderedSystem & DS, const int Hopps)
{
 	int current_pos = DS.Carriers[0].Get_Loc_State_Index();
	while (DS.Carriers[0].Get_Hops_Counter() < Hopps)
	{
		current_pos = next_hop(DS.Hopping_Probs_Table[current_pos]); 
		DS.Move_Carrier(0, current_pos);
		DS.Carriers[0].Increase_Hops_Counter(1);
	}
	DS.Carriers[0].Set_Hops_Counter(0);
}
double	mobility_calculation(DisorderedSystem & DS, const int Hopps)
{
	random_device rd1;
	mt19937 generator1(rd1());
	uniform_real_distribution <double> distribution(0, 1);
		
	double mobility = 0;
	double carr_path = 0;
	int current_pos = DS.Carriers[0].Get_Loc_State_Index();
	int new_pos = 0;
	while (DS.Carriers[0].Get_Hops_Counter() < Hopps)
	{
		new_pos = next_hop(DS.Hopping_Probs_Table[current_pos]);
		double CummRate = DS.Hopping_Rates_Table[current_pos][DS.Hopping_Rates_Table[current_pos].size() - 1].second;
		double t = distribution(generator1);
		DS.Carriers[0].Increase_Timer(-log(t) / CummRate);
		DS.Carriers[0].Increase_Hops_Counter(1);
		double p = distribution(generator1);
		if (p > DS.LS_Fermi_Occupation_Probs[new_pos])
		{
			DS.Move_Carrier(0, new_pos);
			int dx = DS.LS[new_pos].Get_X() - DS.LS[current_pos].Get_X();
			if (abs(dx) < DS.Get_Sys_Size() / 2) carr_path += dx;
			else
			{
				dx = ((dx < 0) - (dx > 0)) * (DS.Get_Sys_Size() - abs(dx));
				carr_path += dx;
			}
			current_pos = new_pos;
		}
		else continue;
	}
	return carr_path / (DS.Carriers[0].Get_Timer_Value() * DS.Get_dls_Field());
}
