#include "DMS.h"

namespace DMS
{
	//LOCALIZED STATE CLASS IMPLEMENTATION
	LocState::LocState()
	{
		// Default constructor
		// Setting all of LS parameters to "0" / "False"
		index = 0;
		type = 0;
		X = 0;
		Y = 0;
		Z = 0;
		E = 0;
		Loc_Rad = 0;
		is_occupied = false;
		occupation_prob = 0;
		Loc_Carrier_index = -1;
		is_radiative = false; 
	}
	LocState::LocState(int LS_index, int LS_type, double LS_x_coord, double LS_y_coord,
					   double LS_z_coord, double LS_energy, double LS_Loc_Rad, 
					   bool occupation_flag, bool radiative_flag)
	{
		// Parametrized constuctor
		index = LS_index;
		type = LS_type;
		X = LS_x_coord;
		Y = LS_y_coord;
		Z = LS_z_coord;
		E = LS_energy;
		Loc_Rad = LS_Loc_Rad;
		is_occupied = occupation_flag;
		occupation_prob = 0;
		Loc_Carrier_index = -1;
		is_radiative = radiative_flag;
	}

	void LocState::Set_Index(int LS_index)
	{
		// LS index seter
		index = LS_index;
	}
	void LocState::Set_Type(int LS_type)
	{
		// LS type seter
		type = LS_type;
	}
	void LocState::Set_X(double LS_x_coord)
	{
		// LS x-coordinate seter
		X = LS_x_coord;
	}
	void LocState::Set_Y(double LS_y_coord)
	{
		// LS y-coordinate seter
		Y = LS_y_coord;
	}
	void LocState::Set_Z(double LS_z_coord)
	{
		// LS z-coordinate seter
		Z = LS_z_coord;
	}
	void LocState::Set_E(double LS_energy)
	{
		// LS energy seter
		E = LS_energy;
	}
	void LocState::Set_Loc_Rad(double LS_Loc_Rad)
	{
		Loc_Rad = LS_Loc_Rad;
	}
	void LocState::Set_Occupied_Flag(bool occupation_flag)
	{
		// LS occupied-flag seter
		is_occupied = occupation_flag;
	}
	void LocState::Set_Occupation_Probability(double occupation_probability)
	{
		occupation_prob = occupation_probability;
	}
	void LocState::Set_Loc_Carrier_index(int LC_index)
	{
		Loc_Carrier_index = LC_index;
	}
	void LocState::Set_Radiative_Flag(bool radiative_flag)
	{
		// LS radiative-flag seter
		is_radiative = radiative_flag;
	}

	int    LocState::Get_Index() const noexcept
	{
		return index;
	}
	int    LocState::Get_Type() const noexcept
	{
		return type;
	}
	double LocState::Get_X() const noexcept
	{
		return X;
	}
	double LocState::Get_Y() const noexcept
	{
		return Y;
	}
	double LocState::Get_Z() const noexcept
	{
		return Z;
	}
	double LocState::Get_E() const noexcept
	{
		return E;
	}
	double LocState::Get_Loc_Rad() const noexcept
	{
		return Loc_Rad;
	}
	bool   LocState::Get_Occupied_Flag() const noexcept
	{
		return is_occupied;
	}
	double LocState::Get_Occupation_Probability() const noexcept
	{
		return occupation_prob;
	}
	int    LocState::Get_Loc_Carrier_index() const noexcept
	{
		return Loc_Carrier_index;
	}
	bool   LocState::Get_Radiative_Flag() const noexcept
	{
		return is_radiative;
	}


	// CARRIER CLASS IMPLEMENTATION
	Carrier::Carrier()
	{
		// Default constructor
		// Creates delocalized electron
		index = 0;
		type = 0;
		loc_state_index = -1;
		life_time = 0;
		is_localized = 0;
		timer = 0;
		hops_counter = 0;
	}
	Carrier::Carrier(int Carr_index, int Carr_type, int C_LS_index, 
					 double Carr_life_time, bool localized_flag)
	{
		//Parametrized constructor
		index = Carr_index;
		type = Carr_type;
		loc_state_index = C_LS_index;
		life_time = Carr_life_time;
		is_localized = localized_flag;
		timer = 0;
		hops_counter = 0;
	}

	void Carrier::Set_Index(int Carr_index)
	{
		index = Carr_index;
	}
	void Carrier::Set_Type(int Carr_type)
	{
		type = Carr_type;
	}
	void Carrier::Set_Loc_State_Index(int C_LS_index)
	{
		loc_state_index = C_LS_index;
	}
	void Carrier::Set_Life_Time(double Carr_life_time)
	{
		life_time = Carr_life_time;
	}
	void Carrier::Set_Localized_Flag(bool localized_flag)
	{
		is_localized = localized_flag;
	}
	void Carrier::Set_Timer_Value(double Carr_time)
	{
		timer = Carr_time;
	}
	void Carrier::Set_Hops_Counter(int Carr_hops_counter)
	{
		hops_counter = Carr_hops_counter;
	}

	int    Carrier::Get_Index() const noexcept
	{
		return index;
	}
	int    Carrier::Get_Type() const noexcept
	{
		return type;
	}
	int    Carrier::Get_Loc_State_Index() const noexcept
	{
		return loc_state_index;
	}
	double Carrier::Get_Life_Time() const noexcept
	{
		return life_time;
	}
	bool   Carrier::Get_Localized_Flag() const noexcept
	{
		return is_localized;
	}
	double Carrier::Get_Timer_Value() const noexcept
	{
		return timer;
	}
	int    Carrier::Get_Hops_Counter() const noexcept
	{
		return hops_counter;
	}

	void Carrier::Increase_Timer(double Carr_timer_incr) 
	{
		timer += Carr_timer_incr;
	}
	void Carrier::Decrease_Timer(double Carr_timer_decr)
	{
		timer -= Carr_timer_decr;
	}

	void Carrier::Increase_Hops_Counter(int Carr_hops_counter_incr)
	{
		hops_counter += Carr_hops_counter_incr;
	}
	void Carrier::Decrease_Hops_Counter(int Carr_hops_counter_decr)
	{
		hops_counter -= Carr_hops_counter_decr;
	}

	// DISORDERED SYSTEM CLASS IMPLEMENTATION
	DisorderedSystem::DisorderedSystem()
	{
		// Default constructor
		// Creates empty system with no LS, zero-temperature [ALL ]
		Dim							= 0;
		N							= 0;
		Loc_Carr_Num				= 0;
		Free_Carr_Num				= 0;
		LS_index_counter			= 0;
		carrier_index_counter		= 0;
		E0						    = 0;
		E1						    = 0;
		Traps_Conc					= 0;
		Sys_Size					= 0;
		Sys_Loc_Rad					= 0;
		Cutoff_Rad					= 0;
		Niu0						= 0;
		T							= 0;
		kT							= 0;
		F							= 0;
		timer					    = 0;
		_E0							= 0;
		_E1							= 0;
		_Traps_Conc					= 0;
		_F							= 0;
		T_Eff_Factor			    = 0;
		T_eff                       = 0;
		kT_eff                      = 0;
	}
	DisorderedSystem::DisorderedSystem(int LS_Num, double DS_energy_scale_1, 
									   double DS_Loc_Rad, double DS_Cutoff, double DS_Traps_Concentration, 
		                               double DS_Temp, double DS_Field)
	{
		// parametrized constructor for 3d-cubic system, empty, room-temperature, Gauss-distrib
		//seting up systems parameters
		Dim							= 3;
		N							= LS_Num;
		Loc_Carr_Num				= 0;
		Free_Carr_Num			    = 0;
		LS_index_counter			= 0;
		carrier_index_counter		= 0;
		E0							= DS_energy_scale_1;
		E1						    = 0;
		Traps_Conc					= DS_Traps_Concentration;
		Sys_Size					= std::cbrt(N / Traps_Conc);
		Sys_Loc_Rad					= DS_Loc_Rad * a_B;
		Cutoff_Rad					= DS_Cutoff;
		Niu0					    = 1e12;
		T							= DS_Temp;
		kT							= k*T;
		F							= DS_Field;
		
		_E0							= E0 / kT;
		_E1							= E1 / kT;
		_Traps_Conc					= Traps_Conc;  
		_F							= (F * Sys_Loc_Rad / cbrt(Traps_Conc)) / E0;

		Calc_Effective_T(0.67);

		timer = 0;

		std::random_device rd;
		std::mt19937 generator(rd());
		std::uniform_real_distribution <double> uni_distribution(-0.5, 0.5);
		std::normal_distribution <double> norm_distribution(0, _E0);

		// generating LS's
		LS.resize(N);
		Carriers.reserve(N);
		LS_Occupation_Status.resize(N);
		double E, x, y, z;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			E = norm_distribution(generator);	
			x = uni_distribution(generator) * Sys_Size;
			y = uni_distribution(generator) * Sys_Size;
			z = uni_distribution(generator) * Sys_Size;
			LS[_LS_count].LocState::LocState(LS_index_counter, 0, x, y, z, E, 
										     Sys_Loc_Rad, false, false);
			LS_Occupation_Status[_LS_count] = 0;
			LS_index[LS_index_counter] = &LS[_LS_count]; LS_index_counter++;
		}
	}
	DisorderedSystem::DisorderedSystem(int LS_Num, double DS_energy_scale_1,
									   double DS_Traps_Concentration, double DS_Cutoff, double DS_Field) 
	{
		// parametrized constructor for 3d-cubic system, empty, room-temperature, Gauss-distrib
		//seting up systems parameters [DIMENSIONLESS PARAMETERS]
		Dim = 3;
		N = LS_Num;
		Loc_Carr_Num = 0;
		Free_Carr_Num = 0;
		LS_index_counter = 0;
		carrier_index_counter = 0;
		Niu0 = 1e12;
		T = 300;
		kT = k*T;

		_E0 = DS_energy_scale_1;
		_E1 = 0;
		_Traps_Conc = DS_Traps_Concentration;
		_F = DS_Field;

		E0 = _E0*kT;
		E1 = 0;
		Traps_Conc = _Traps_Conc;
		Sys_Size = std::cbrt(N / Traps_Conc);
		Sys_Loc_Rad = a_B;
		Cutoff_Rad = DS_Cutoff;
		F = (_F * E0) * cbrt(_Traps_Conc) / Sys_Loc_Rad;

		Calc_Effective_T(0.67);

		timer = 0;

		std::random_device rd;
		std::mt19937 generator(rd());
		std::uniform_real_distribution <double> uni_distribution(0, 1);
		std::normal_distribution <double> norm_distribution(0, _E0);

		// generating LS's
		LS.resize(N);
		Carriers.reserve(N);
		LS_Occupation_Status.resize(N);
		double E, x, y, z;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			E = norm_distribution(generator);
			x = uni_distribution(generator) * Sys_Size;
			y = uni_distribution(generator) * Sys_Size;
			z = uni_distribution(generator) * Sys_Size;
			LS[_LS_count].LocState::LocState(LS_index_counter, 0, x, y, z, E,
				Sys_Loc_Rad, false, false);
			LS_Occupation_Status[_LS_count] = 0;
			LS_index[LS_index_counter] = &LS[_LS_count]; LS_index_counter++;
		}
	}

	void DisorderedSystem::Set_Dim(int DS_dimensions_num)
	{
		Dim = DS_dimensions_num;
	}
	void DisorderedSystem::Set_N(int LS_Num)
	{
		// LS number seter
		N = LS_Num;
	}
	void DisorderedSystem::Set_Sys_Size(double DS_size)
	{
		// System size seter
		Sys_Size = DS_size;
	}
	void DisorderedSystem::Set_E0(double DS_energy_scale_exp)
	{
		// Disordere energy scale seter
		E0 = DS_energy_scale_exp;
	}
	void DisorderedSystem::Set_Sys_Loc_Rad(double DS_Loc_Rad)
	{
		Sys_Loc_Rad = DS_Loc_Rad;
	}
	void DisorderedSystem::Set_Cutoff_Rad(double DS_Cutoff)
	{
		Cutoff_Rad = DS_Cutoff;
	}
	void DisorderedSystem::Set_Traps_Conc(double DS_Traps_concentration)
	{
		Traps_Conc = DS_Traps_concentration;
	}
	void DisorderedSystem::Set_Escape_Rate(double LS_Escape_rate)
	{
		Niu0 = LS_Escape_rate;
	}
	void DisorderedSystem::Set_Sys_Temperature(double DS_Temperature)
	{
		T = DS_Temperature;
		kT = k*T;
	}
	void DisorderedSystem::Set_kT(double DS_kT)
	{
		kT = DS_kT;
		T = kT / k;
	}
	void DisorderedSystem::Set_Field(double DS_Field)
	{
		F = DS_Field;
	}
	void DisorderedSystem::Set_Timer_Value(double DS_time)
	{
		timer = DS_time;
	}
	void DisorderedSystem::Set_dls_E0(double _DS_energy_scale_exp)
	{
		_E0 = _DS_energy_scale_exp;
		E0  = _DS_energy_scale_exp*kT;
	}
	void DisorderedSystem::Set_dls_Traps_Conc(double _DS_Traps_concentration)
	{
		_Traps_Conc = _DS_Traps_concentration;
		Traps_Conc  = _DS_Traps_concentration;
	}
	void DisorderedSystem::Set_dls_Field(double _DS_Field)
	{
		_F = _DS_Field;
		F  = (_F * E0) * pow(Traps_Conc, 1/Dim) / Sys_Loc_Rad;
	}
	void DisorderedSystem::Set_E_Fermi(double DS_Fermi_level)
	{
		E_Fermi  = DS_Fermi_level;
		_E_Fermi = E_Fermi / kT;
	}
	void DisorderedSystem::Set_dls_E_Fermi(double _DS_Fermi_level)
	{
		_E_Fermi = _DS_Fermi_level;
		E_Fermi  = _E_Fermi * kT;
	}

	int	   DisorderedSystem::Get_Dim() const noexcept
	{
		return Dim;
	}
	int    DisorderedSystem::Get_N() const noexcept
	{
		// LS number geter 
		return N;
	}
	double DisorderedSystem::Get_Sys_Size() const noexcept
	{
		// System syze geter
		return Sys_Size;
	}
	double DisorderedSystem::Get_E0() const noexcept
	{
		// Disordere energy scale geter
		return E0;
	}
	double DisorderedSystem::Get_Sys_Loc_Rad() const noexcept
	{
		return Sys_Loc_Rad;
	}
	double DisorderedSystem::Get_Cutoff_Rad() const noexcept
	{
		return Cutoff_Rad;
	}
	double DisorderedSystem::Get_Traps_Conc() const noexcept
	{
		return Traps_Conc;
	}
	int    DisorderedSystem::Get_Loc_Carr_Num() const noexcept
	{
		return Loc_Carr_Num;
	}
	int    DisorderedSystem::Get_Free_Carr_Num() const noexcept
	{
		return Free_Carr_Num;
	}
	double DisorderedSystem::Get_Escape_Rate() const noexcept
	{
		return Niu0;
	}
	double DisorderedSystem::Get_Sys_Temperature() const noexcept
	{
		return T;
	}
	double DisorderedSystem::Get_kT() const noexcept
	{
		return kT;
	}
	double DisorderedSystem::Get_Field() const noexcept
	{
		return F;
	}
	double DisorderedSystem::Get_Timer_Value() const noexcept
	{
		return timer;
	}
	double DisorderedSystem::Get_dls_E0() const noexcept
	{
		return _E0;
	}
	double DisorderedSystem::Get_dls_Traps_Conc() const noexcept
	{
		return _Traps_Conc;
	}
	double DisorderedSystem::Get_dls_Field() const noexcept
	{
		return _F;
	}
	double DisorderedSystem::Get_E_Fermi() const noexcept
	{
		return E_Fermi;
	}
	double DisorderedSystem::Get_dls_E_Fermi() const noexcept
	{
		return _E_Fermi;
	}
	double DisorderedSystem::Get_T_Eff_Factor() const noexcept
	{
		return T_Eff_Factor;
	}
	double DisorderedSystem::Get_T_Eff() const noexcept
	{
		return T_eff;
	}
	double DisorderedSystem::Get_kT_Eff_Factor() const noexcept
	{
		kT_eff;
	}
	double DisorderedSystem::Get_E_Transport() const noexcept
	{
		return E_Transport;
	}
	double DisorderedSystem::Get_dls_E_Transport() const noexcept
	{
		return _E_Transport;
	}

	void DisorderedSystem::Calc_Effective_T(double gamma)
	{
		T_Eff_Factor = sqrt(1 + pow((_F * E0) * cbrt(_Traps_Conc) * gamma / kT, 2));
		T_eff = T * T_Eff_Factor;
		kT_eff = k * T_eff;
	}

	void DisorderedSystem::Calc_Transport_Energy_ET(double precision)
	{
		const double fermi_level = _E_Fermi;
		const double effective_temp_factor = T_Eff_Factor;
		const double traps_concentration = _Traps_Conc;
		const double sigma = _E0;
		const double B_c = 2.7;		// percolation criterion
		const double prefactor = (2.0/3.0) * pow((4.0 * 3.1415) / (3.0 * B_c), - 1.0 / 3.0);
		
		// this lamda returns value of Fermi-function f(E) WITH EFFECTIVE TEMPERATURE
		auto fermi_function = [fermi_level, effective_temp_factor] (double energy) -> double {
			return 1 / (1 + exp((energy - fermi_level) / effective_temp_factor));
		};
		// this lamda returns value of DOS(E)
		auto gaussian_dos = [traps_concentration, sigma](double energy) {
			return (traps_concentration / (sigma * sqrt(2 * 3.1415))) * exp(-energy * energy / (2 * sigma * sigma));
		};
		// this lambda returns value of integral [1 - f(E)DOS(E)]dE over [ - inf(~5*sigma) , transport_energy ]
		auto free_traps_integral = [sigma, fermi_function, gaussian_dos](double transport_energy, int grid_steps_count) {
			double result = 0;
			double step = (transport_energy - (-5 * sigma)) / grid_steps_count;
			for (int i = 0; i < grid_steps_count; i++)
				result += (1 - fermi_function((-5 * sigma) + step * (i + 0.5))) * gaussian_dos((-5 * sigma) + step * (i + 0.5));
			result *= step;
			return result;
		};
		// calculating transport energy using method of dihotomy
		double left_TE = -5 * sigma;
		double right_TE = 5 * sigma;
		while ((right_TE - left_TE) > precision) {
			double median_TE = (left_TE + right_TE) / 2;
			double left_value = prefactor * pow(free_traps_integral(left_TE, 1e5), -4.0 / 3.0) 
				                          * (1 - fermi_function(left_TE)) * gaussian_dos(left_TE) - 1;
			double right_value = prefactor * pow(free_traps_integral(right_TE, 1e5), -4.0 / 3.0)
									       * (1 - fermi_function(right_TE)) * gaussian_dos(right_TE) - 1;
			double median_value = prefactor * pow(free_traps_integral(median_TE, 1e5), -4.0 / 3.0)
				                            * (1 - fermi_function(median_TE)) * gaussian_dos(median_TE) - 1;
			if (median_value * right_value < 0) left_TE = median_TE;
			else right_TE = median_TE;
		}
		_E_Transport = (right_TE + left_TE) / 2;
		E_Transport = _E_Transport * kT;
	}

	void DisorderedSystem::Increase_Timer(double DS_timer_incr)
	{
		timer += DS_timer_incr;
	}
	void DisorderedSystem::Decrease_Timer(double DS_timer_decr)
	{
		timer += DS_timer_decr;
	}

	void DisorderedSystem::Gen_Coords_Distrib_3D_Random()
	{
		std::random_device rd;
		std::mt19937 generator(rd());
		std::uniform_real_distribution <double> uni_distribution(0, 1);
		for (auto& _LS : LS)
		{
			_LS.Set_X(uni_distribution(generator) * Sys_Size);
			_LS.Set_Y(uni_distribution(generator) * Sys_Size);
			_LS.Set_Z(uni_distribution(generator) * Sys_Size);
		}
	}
	void DisorderedSystem::Gen_Coords_Distrib_3D_Cubic_Lattice()
	{
		const int count = std::cbrt(N);
		const double dl = Sys_Size / count;
		int i = 0;
		for (int x = 0; x < count; x++)
		{
			for (int y = 0; y < count; y++)
			{
				for (int z = 0; z < count; z++)
				{
					if (i < N)
					{
						LS[i].Set_X(x*dl);
						LS[i].Set_Y(y*dl);
						LS[i].Set_Z(z*dl);
						i++;
					}
				}
			}
		}
	}
	void DisorderedSystem::Gen_Energies_Distrib_Exp()
	{
		std::random_device rd;
		std::mt19937 generator(rd());
		std::exponential_distribution <double> exp_distribution(1 / E0);
		for (auto& _LS : LS)
			_LS.Set_E(-exp_distribution(generator));
	}
	void DisorderedSystem::Gen_Energies_Distrib_Gauss()
	{
		std::random_device rd;
		std::mt19937 generator(rd());
		std::normal_distribution <double> norm_distribution(0, E0);
		for (auto& _LS : LS)
			_LS.Set_E(-abs(norm_distribution(generator)));
	}

	void DisorderedSystem::Set_LS_Index_Ptrs()
	{
		// Setup map of LS-pointinters for index-access
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
			LS_index[LS[_LS_count].Get_Index()] = &LS[_LS_count];
	}
	int  DisorderedSystem::Get_LS_Position(int _LS_index)
	{
		// Return position of the LS in LS-vector by LS-index
		int _LS_pos = 0;
		while ((LS_index[_LS_index] != &LS[_LS_pos]) && (_LS_pos < LS.size()))
			_LS_pos++;

		return _LS_pos;
	}
	void DisorderedSystem::Set_Carriers_Index_Ptrs()
	{
		// Setup of Carriers-pointinters for index-acces
		for (int _LC_count = 0; _LC_count < Loc_Carr_Num + Free_Carr_Num; _LC_count++)
			Carriers_index[Carriers[_LC_count].Get_Index()] = &Carriers[_LC_count];
	}
	int  DisorderedSystem::Get_Carriers_Position(int _Carr_index)
	{
		// Return position of the Carrier in LocCarriers-vector by LocCarr-index
		int _LC_pos = 0;
		while ((Carriers_index[_Carr_index] != &Carriers[_LC_pos]) && (_LC_pos < Carriers.size()))
			_LC_pos++;

		return _LC_pos;
	}
	
	void DisorderedSystem::Add_Trap(int LS_type, double LS_x_coord, double LS_y_coord, 
									double LS_z_coord, double LS_energy, 
									double LS_Loc_Rad, bool radiative_flag)
	{
	// Add new trap to the end of LS - vector, index = max__current_index + 1
		N++;
		LS.resize(N);
		LS_Occupation_Status.resize(N);
		LS_Occupation_Status[N - 1] = 0;
		LS[N - 1].LocState::LocState(LS_index_counter, LS_type, LS_x_coord, LS_y_coord, 
			                         LS_z_coord, LS_energy, LS_Loc_Rad, false, radiative_flag);
		LS_index[LS_index_counter] = &LS[N - 1]; LS_index_counter++;
	}
	void DisorderedSystem::Remove_Trap(int _LS_index)
	{
		//removes by postition _LS_index 
		N--;
		LS_index.erase(LS[_LS_index].Get_Index());
		LS.erase(LS.begin() + _LS_index);
		LS_Occupation_Status.erase(LS_Occupation_Status.begin() + _LS_index);
		Set_LS_Index_Ptrs();
	}
	void DisorderedSystem::Remove_Trap_IA(int _LS_index)
	{
		//removes by index_<_LS_index>
		N--;
		Set_LS_Index_Ptrs();
		int _LS_pos = Get_LS_Position(_LS_index);
		LS_index.erase(_LS_index);
		LS.erase(LS.begin() + _LS_pos);
		LS_Occupation_Status.erase(LS_Occupation_Status.begin() + _LS_pos);
		Set_LS_Index_Ptrs();
	}

	void DisorderedSystem::Add_Carrier(int Carr_type, int C_LS_index, 
							double Carr_life_time, bool localized_flag)
	{
		// Add new carrier localized / 
		//free  <->  localized_flag = 0 / 1 
		if (localized_flag)
		{
			Loc_Carr_Num++;
			LS[C_LS_index].Set_Occupied_Flag(true);
			LS_Occupation_Status[C_LS_index] = true;
			LS[C_LS_index].Set_Loc_Carrier_index(carrier_index_counter); 
			Carriers.emplace_back(Carrier::Carrier(carrier_index_counter, Carr_type, 
								                   C_LS_index, Carr_life_time, localized_flag)); 
		}
		else
		{
			Free_Carr_Num++;
			Carriers.emplace_back(Carrier::Carrier(carrier_index_counter, Carr_type, -1, 
								                   Carr_life_time, localized_flag)); 
		}
		Carriers_index[carrier_index_counter] = (&Carriers[Loc_Carr_Num + Free_Carr_Num - 1]); 
		carrier_index_counter++;
	}
	void DisorderedSystem::Add_Carrier_IA(int Carr_type, int C_LS_index, 
							double Carr_life_time, bool localized_flag)
	{
		// Add new carrier by index_<LS_index> localized / 
		// free  <->  localized_flag = 0 / 1 
		if (localized_flag)
		{
			Loc_Carr_Num++;
			LS_index[C_LS_index]->Set_Occupied_Flag(true);
			LS_Occupation_Status[Get_LS_Position(C_LS_index)] = true;
			LS_index[C_LS_index]->Set_Loc_Carrier_index(carrier_index_counter);
			Carriers.emplace_back(Carrier::Carrier(carrier_index_counter, Carr_type, 
				                                   Get_LS_Position(C_LS_index), 
												   Carr_life_time, localized_flag));
		}
		else
		{
			Free_Carr_Num++;
			Carriers.emplace_back(Carrier::Carrier(carrier_index_counter, Carr_type, -1, 
				                                   Carr_life_time, localized_flag));
		}
		Carriers_index[carrier_index_counter] = (&Carriers[Loc_Carr_Num + Free_Carr_Num - 1]);  
		carrier_index_counter++;
	}
	void DisorderedSystem::Remove_Carrier(int Carr_index)
	{
		// Removre carreir by pos Carr_index
		if (Carriers[Carr_index].Get_Localized_Flag())
		{
			Loc_Carr_Num--;
			int indx = Carriers[Carr_index].Get_Loc_State_Index();
			LS[indx].Set_Occupied_Flag(false);
			LS_Occupation_Status[indx] = false;
			LS[indx].Set_Loc_Carrier_index(-1);
			Carriers_index.erase(Carriers[Carr_index].Get_Index());
			Carriers.erase(Carriers.begin() + Carr_index);
			Set_Carriers_Index_Ptrs();
		}
		else
		{
			Free_Carr_Num--;
			Carriers.erase(Carriers.begin() + Carr_index);
			Carriers_index.erase(Carriers[Carr_index].Get_Index());
			Set_Carriers_Index_Ptrs();
		}
	}
	void DisorderedSystem::Remove_Carrier_IA(int Carr_index)
	{
		// Removre carreir by index_<Carr_index> 
		if (Carriers_index[Carr_index]->Get_Localized_Flag())
		{
			Loc_Carr_Num--;
			int _Carr_pos = Get_Carriers_Position(Carr_index);
			int indx = Carriers[_Carr_pos].Get_Loc_State_Index();
			LS[indx].Set_Occupied_Flag(false);
			LS_Occupation_Status[indx] = false;
			LS[indx].Set_Loc_Carrier_index(-1);
			Carriers_index.erase(Carr_index);
			Carriers.erase(Carriers.begin() + _Carr_pos);
			Set_Carriers_Index_Ptrs();
		}
		else
		{
			Set_Carriers_Index_Ptrs();
			int _Carr_pos = Get_Carriers_Position(Carr_index);
			Carriers.erase(Carriers.begin() + _Carr_pos);
			Carriers_index.erase(Carr_index);
			Set_Carriers_Index_Ptrs();
		}
	}

	double DisorderedSystem::Calc_Distance(int i, int j)
	{
		// returns distance between <i> and <j> sites
		double d_x = LS[j].Get_X() - LS[i].Get_X();
		double d_y = LS[j].Get_Y() - LS[i].Get_Y();
		double d_z = LS[j].Get_Z() - LS[i].Get_Z();

		return sqrt(pow(d_x, 2) + pow(d_y, 2) + pow(d_z, 2));
	}
	double DisorderedSystem::Calc_Distance_BC(int i, int j)
	{
		// returns distance between <i> and <j> sites assuming boundary conditions
		double d_x = LS[j].Get_X() - LS[i].Get_X();
		double d_y = LS[j].Get_Y() - LS[i].Get_Y();
		double d_z = LS[j].Get_Z() - LS[i].Get_Z();
		// setting up boundary conditions
		if (abs(d_x) > Sys_Size / 2) d_x = ((d_x < 0) - (d_x > 0)) * (Sys_Size - abs(d_x));
		if (abs(d_y) > Sys_Size / 2) d_y = ((d_y < 0) - (d_y > 0)) * (Sys_Size - abs(d_y));
		if (abs(d_z) > Sys_Size / 2) d_z = ((d_z < 0) - (d_z > 0)) * (Sys_Size - abs(d_z));
		 
		return sqrt(pow(d_x, 2) + pow(d_y, 2) + pow(d_z, 2));
	}
	double DisorderedSystem::Calc_Energy_Diff(int i, int j)
	{
		// returns energy difference between sites <j> and <i>
		double d_x = LS[j].Get_X() - LS[i].Get_X();
		return LS[j].Get_E() - LS[i].Get_E() - _F * _E0 * pow(_Traps_Conc, 1 / double(Dim)) * d_x;;
	}

	double DisorderedSystem::Calc_Hopping_Rate(int i, int j)
	{
		// returns hopping rate from LS[i] to LS[j] or 
		// activation rate from LS[i] (if i == j)
		if (i == j) return Niu0 * exp(LS[i].Get_E());
		else
		{
			const double d_x = LS[j].Get_X() - LS[i].Get_X();
			const double d_y = LS[j].Get_Y() - LS[i].Get_Y();
			const double d_z = LS[j].Get_Z() - LS[i].Get_Z();
			const double R = sqrt(pow(d_x, 2) + pow(d_y, 2) + pow(d_z, 2));
			const double  d_E = LS[j].Get_E() - LS[i].Get_E() - _F*_E0*pow(Traps_Conc, 1 / Dim)*d_x;
			if ((Cutoff_Rad != 0) && (R > Cutoff_Rad)) return 0;
			else if (d_E > 0) return Niu0 * exp(-2 * R - d_E);
			else return Niu0*exp(-2 * R);
		}
	}
	double DisorderedSystem::Calc_Hopping_Rate_BC(int i, int j)
	{
		// returns hopping rate from LS[i] to LS[j] or 
		// activation rate from LS[i] (if i == j)
		if (i == j) return Niu0 * exp(LS[i].Get_E());
		else
		{
			double d_x = LS[j].Get_X() - LS[i].Get_X();
			double d_y = LS[j].Get_Y() - LS[i].Get_Y();
			double d_z = LS[j].Get_Z() - LS[i].Get_Z();
			// setting up boundary conditions
			if (abs(d_x) > Sys_Size / 2) d_x = ((d_x < 0) - (d_x > 0)) * (Sys_Size - abs(d_x));
			if (abs(d_y) > Sys_Size / 2) d_y = ((d_y < 0) - (d_y > 0)) * (Sys_Size - abs(d_y));
			if (abs(d_z) > Sys_Size / 2) d_z = ((d_z < 0) - (d_z > 0)) * (Sys_Size - abs(d_z));
			// calcualting hopping rates
			const double R = sqrt(pow(d_x, 2) + pow(d_y, 2) + pow(d_z, 2));
			const double d_E = LS[j].Get_E() - LS[i].Get_E() - _F * _E0 * pow(_Traps_Conc, 1/double(Dim)) * d_x;
			if ((Cutoff_Rad != 0) && (R > Cutoff_Rad)) return 0;
			else if (d_E > 0) return Niu0 * exp(-2 * R - d_E);
			else return Niu0 * exp(-2 * R);
		}
	}
	double DisorderedSystem::Calc_Hopping_Rate_IA(int i, int j)
	{
		// returns hopping rate from LS_index[i] to LS_index[j] or 
		// activation rate from LS_index[i] (if i == j)
		if (i == j) return Niu0 * exp(LS_index[i]->Get_E());
		else
		{
			const double d_x = LS_index[j]->Get_X() - LS_index[i]->Get_X();
			const double d_y = LS_index[j]->Get_Y() - LS_index[i]->Get_Y();
			const double d_z = LS_index[j]->Get_Z() - LS_index[i]->Get_Z();
			const double R   = sqrt(pow(d_x, 2) + pow(d_y, 2) + pow(d_z, 2));
			const double d_E = LS_index[j]->Get_E() - LS_index[i]->Get_E() - _F*_E0*pow(Traps_Conc, 1 / double(Dim))*d_x;
			if ((Cutoff_Rad != 0) && (R > Cutoff_Rad)) return 0;
			else if (d_E > 0) return Niu0 * exp(-2 * R - d_E);
			else return Niu0*exp(-2 * R);
		}
	}
	double DisorderedSystem::Calc_Hopping_Rate_BC_IA(int i, int j)
	{
		// returns hopping rate from LS_index[i] to LS_index[j] or 
		// activation rate from LS_index[i] (if i == j)
		if (i == j) return Niu0 * exp(LS_index[i]->Get_E());
		else
		{
			double d_x = LS_index[j]->Get_X() - LS_index[i]->Get_X();
			double d_y = LS_index[j]->Get_Y() - LS_index[i]->Get_Y();
			double d_z = LS_index[j]->Get_Z() - LS_index[i]->Get_Z();
			// setting up boundary conditions
			if (abs(d_x) > Sys_Size / 2) d_x = ((d_x < 0) - (d_x > 0))*(Sys_Size - abs(d_x));
			if (abs(d_y) > Sys_Size / 2) d_y = ((d_y < 0) - (d_y > 0))*(Sys_Size - abs(d_y));
			if (abs(d_z) > Sys_Size / 2) d_z = ((d_z < 0) - (d_z > 0))*(Sys_Size - abs(d_z));
			const double R = sqrt(pow(d_x, 2) + pow(d_y, 2) + pow(d_z, 2));
			const double d_E = LS_index[j]->Get_E() - LS_index[i]->Get_E() - _F*_E0*pow(Traps_Conc, 1 / double(Dim))*d_x;
			if ((Cutoff_Rad != 0) && (R > Cutoff_Rad)) return 0;
			else if (d_E > 0) return Niu0 * exp(-2 * R - d_E);
			else return Niu0*exp(-2 * R);
		}
	}
	double DisorderedSystem::Calc_Hopping_Rate_ZeroT(int i, int j)
	{
		// returns hopping rate from LS[i] to LS[j] for zero T case
		if (i == j) return Niu0 * exp(LS[i].Get_E());
		else
		{
			const double d_x = LS[j].Get_X() - LS[i].Get_X();
			const double d_y = LS[j].Get_Y() - LS[i].Get_Y();
			const double d_z = LS[j].Get_Z() - LS[i].Get_Z();
			const double R = sqrt(pow(d_x, 2) + pow(d_y, 2) + pow(d_z, 2));
			const double d_E = LS[j].Get_E() - LS[i].Get_E() - _F*_E0*pow(Traps_Conc, 1 / double(Dim))*d_x;
			if ((Cutoff_Rad != 0) && (R > Cutoff_Rad) || (d_E > 0)) return 0;
			else  return Niu0 * exp(-2 * R);
		}
	}
	double DisorderedSystem::Calc_Hopping_Rate_ZeroT_BC(int i, int j)
	{
		// returns hopping rate from LS[i] to LS[j] for zero T case
		if (i == j) return Niu0 * exp(LS[i].Get_E());
		else
		{
			double d_x = LS[j].Get_X() - LS[i].Get_X();
			double d_y = LS[j].Get_Y() - LS[i].Get_Y();
			double d_z = LS[j].Get_Z() - LS[i].Get_Z();
			// setting up boundary conditions
			if (abs(d_x) > Sys_Size / 2) d_x = ((d_x < 0) - (d_x > 0))*(Sys_Size - abs(d_x));
			if (abs(d_y) > Sys_Size / 2) d_y = ((d_y < 0) - (d_y > 0))*(Sys_Size - abs(d_y));
			if (abs(d_z) > Sys_Size / 2) d_z = ((d_z < 0) - (d_z > 0))*(Sys_Size - abs(d_z));
			const double R = sqrt(pow(d_x, 2) + pow(d_y, 2) + pow(d_z, 2));
			const double d_E = LS[j].Get_E() - LS[i].Get_E() - _F*_E0*pow(Traps_Conc, 1 / double(Dim))*d_x;
			if ((Cutoff_Rad != 0) && (R > Cutoff_Rad) || (d_E > 0)) return 0;
			else  return Niu0 * exp(-2 * R);
		}
	}
	double DisorderedSystem::Calc_Hopping_Rate_ZeroT_IA(int i, int j)
	{
		// returns hopping rate from LS_index[i] to LS_index[j] for zero T case
		if (i == j) return Niu0 * exp(LS[i].Get_E());
		else
		{
			const double d_x = LS_index[j]->Get_X() - LS_index[i]->Get_X();
			const double d_y = LS_index[j]->Get_Y() - LS_index[i]->Get_Y();
			const double d_z = LS_index[j]->Get_Z() - LS_index[i]->Get_Z();
			const double R = sqrt(pow(d_x, 2) + pow(d_y, 2) + pow(d_z, 2));
			const double d_E = LS_index[j]->Get_E() - LS_index[i]->Get_E() - _F*_E0*pow(Traps_Conc, 1 / double(Dim))*d_x;
			if ((Cutoff_Rad != 0) && (R > Cutoff_Rad) || (d_E > 0)) return 0;
			else  return Niu0 * exp(-2 * R);
		}
	}
	double DisorderedSystem::Calc_Hopping_Rate_ZeroT_BC_IA(int i, int j)
	{
		// returns hopping rate from LS_index[i] to LS_index[j] for zero T case
		if (i == j) return Niu0 * exp(LS[i].Get_E());
		else
		{
			double d_x = LS_index[j]->Get_X() - LS_index[i]->Get_X();
			double d_y = LS_index[j]->Get_Y() - LS_index[i]->Get_Y();
			double d_z = LS_index[j]->Get_Z() - LS_index[i]->Get_Z();
			// setting up boundary conditions
			if (abs(d_x) > Sys_Size / 2) d_x = ((d_x < 0) - (d_x > 0))*(Sys_Size - abs(d_x));
			if (abs(d_y) > Sys_Size / 2) d_y = ((d_y < 0) - (d_y > 0))*(Sys_Size - abs(d_y));
			if (abs(d_z) > Sys_Size / 2) d_z = ((d_z < 0) - (d_z > 0))*(Sys_Size - abs(d_z));
			const double R = sqrt(pow(d_x, 2) + pow(d_y, 2) + pow(d_z, 2));
			const double d_E = LS_index[j]->Get_E() - LS_index[i]->Get_E() - _F*_E0*pow(Traps_Conc, 1 / double(Dim))*d_x;
			if ((Cutoff_Rad != 0) && (R > Cutoff_Rad) || (d_E > 0)) return 0;
			else  return Niu0 * exp(-2 * R);
		}
	}

	std::vector<std::pair<int, double>> DisorderedSystem::Calc_Hopping_Rates(int i)
	{
		// Returns  vector of pairs <j-trap pos, hopping rate> to hop 
		// from trap <i> to trap <j> or activate to extended states
		std::pair<int, double> Pos_Rate;
		std::vector<std::pair<int, double> > Hopp_Probs;
		int cur_pos = i;
		double cum_rate = 0;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			if ((_LS_count != cur_pos) & (LS[_LS_count].Get_Occupied_Flag())) continue;
			else
			{
				Pos_Rate.first = _LS_count;
				Pos_Rate.second = Calc_Hopping_Rate(cur_pos, _LS_count);
				cum_rate += Pos_Rate.second;
				if (Pos_Rate.second != 0) Hopp_Probs.emplace_back(Pos_Rate);
			}
		}
		Pos_Rate.first = N;
		Pos_Rate.second = cum_rate;
		Hopp_Probs.emplace_back(Pos_Rate);
		return Hopp_Probs;
	}
	std::vector<std::pair<int, double>> DisorderedSystem::Calc_Hopping_Rates_BC(int i)
	{
		// Returns  vector of pairs <j-trap pos, hopping rate> to hop 
		// from trap <i> to trap <j> or activate to extended states
		std::pair<int, double> Pos_Rate;
		std::vector<std::pair<int, double> > Hopp_Probs;
		int cur_pos = i;
		double cum_rate = 0;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			if ((_LS_count != cur_pos) & (LS[_LS_count].Get_Occupied_Flag())) continue;
			else
			{
				Pos_Rate.first = _LS_count;
				Pos_Rate.second = Calc_Hopping_Rate_BC(cur_pos, _LS_count);
				cum_rate += Pos_Rate.second;
				if (Pos_Rate.second != 0) Hopp_Probs.emplace_back(Pos_Rate);
			}
		}
		Pos_Rate.first = N;
		Pos_Rate.second = cum_rate;
		Hopp_Probs.emplace_back(Pos_Rate);
		return Hopp_Probs;
	}
	std::vector<std::pair<int, double>> DisorderedSystem::Calc_Hopping_Rates_NA(int i)
	{
		// Returns  vector of pairs <j-trap pos, hopping rate> to hop 
		// from trap <i> to trap <j> , NO-ACTIVATION [NA] allowed
		std::pair<int, double> Pos_Rate;
		std::vector<std::pair<int, double> > Hopp_Probs;
		int cur_pos = i;
		double cum_rate = 0;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			if (LS[_LS_count].Get_Occupied_Flag()) continue;
			else
			{
				Pos_Rate.first = _LS_count;
				Pos_Rate.second = Calc_Hopping_Rate(cur_pos, _LS_count);
				cum_rate += Pos_Rate.second;
				if (Pos_Rate.second != 0) Hopp_Probs.emplace_back(Pos_Rate);
			}
		}
		Pos_Rate.first = N;
		Pos_Rate.second = cum_rate;
		Hopp_Probs.emplace_back(Pos_Rate);
		return Hopp_Probs;
	}
	std::vector<std::pair<int, double>> DisorderedSystem::Calc_Hopping_Rates_NA_BC(int i)
	{
		// Returns  vector of pairs <j-trap pos, hopping rate> to hop 
		// from trap <i> to trap <j> , NO-ACTIVATION [NA] allowed
		std::pair<int, double> Pos_Rate;
		std::vector<std::pair<int, double> > Hopp_Probs;
		int cur_pos = i;
		double cum_rate = 0;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{			
			if ((LS[_LS_count].Get_Occupied_Flag()) || 
				(Calc_Distance_BC(cur_pos, _LS_count) > Cutoff_Rad)) continue;
			else
			{
				Pos_Rate.first = _LS_count;
				Pos_Rate.second = Calc_Hopping_Rate_BC(cur_pos, _LS_count);
				cum_rate += Pos_Rate.second;
				Hopp_Probs.emplace_back(Pos_Rate);
			}
		}
		Pos_Rate.first = N;
		Pos_Rate.second = cum_rate;
		Hopp_Probs.emplace_back(Pos_Rate);
		return Hopp_Probs;
	}
	std::vector<std::pair<int, double>> DisorderedSystem::Calc_Hopping_Rates_IA(int i)
	{	
		// Returns  vector of pairs <j-trap pos, hopping rate> to hop 
		// from trap index_<i> to trap index_<j> [IA] or activate to extended states
		std::pair<int, double> Pos_Rate;
		std::vector<std::pair<int, double> > Hopp_Probs;
		int cur_pos = i;
		double cum_rate = 0;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			if ((_LS_count != cur_pos) & (LS[_LS_count].Get_Occupied_Flag()))  continue;
			else
			{
				Pos_Rate.first = LS[_LS_count].Get_Index();
				Pos_Rate.second = Calc_Hopping_Rate_IA(cur_pos, LS[_LS_count].Get_Index());
				cum_rate += Pos_Rate.second;
				if (Pos_Rate.second != 0) Hopp_Probs.emplace_back(Pos_Rate);
			}
		}
		Pos_Rate.first = N;
		Pos_Rate.second = cum_rate;
		Hopp_Probs.emplace_back(Pos_Rate);
		return Hopp_Probs;
	}
	std::vector<std::pair<int, double>> DisorderedSystem::Calc_Hopping_Rates_BC_IA(int i)
	{
		// Returns  vector of pairs <j-trap pos, hopping rate> to hop 
		// from trap index_<i> to trap index_<j> [IA] or activate to extended states
		std::pair<int, double> Pos_Rate;
		std::vector<std::pair<int, double> > Hopp_Probs;
		int cur_pos = i;
		double cum_rate = 0;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			if ((_LS_count != cur_pos) & (LS[_LS_count].Get_Occupied_Flag()))  continue;
			else
			{
				Pos_Rate.first = LS[_LS_count].Get_Index();
				Pos_Rate.second = Calc_Hopping_Rate_BC_IA(cur_pos, LS[_LS_count].Get_Index());
				cum_rate += Pos_Rate.second;
				if (Pos_Rate.second != 0) Hopp_Probs.emplace_back(Pos_Rate);
			}
		}
		Pos_Rate.first = N;
		Pos_Rate.second = cum_rate;
		Hopp_Probs.emplace_back(Pos_Rate);
		return Hopp_Probs;
	}
	std::vector<std::pair<int, double>> DisorderedSystem::Calc_Hopping_Rates_NA_IA(int i)
	{
		// Returns  vector of pairs <j-trap pos, hopping rate> to hop 
		// from trap index_<i> to trap index_<j> , NO-ACTIVATION [NA] allowed
		std::pair<int, double> Pos_Rate;
		std::vector<std::pair<int, double> > Hopp_Probs;
		int cur_pos = i;
		double cum_rate = 0;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			if (LS[_LS_count].Get_Occupied_Flag()) continue;
			else
			{
				Pos_Rate.first = LS[_LS_count].Get_Index();
				Pos_Rate.second = Calc_Hopping_Rate_IA(cur_pos, LS[_LS_count].Get_Index());
				cum_rate += Pos_Rate.second;
				if (Pos_Rate.second != 0) Hopp_Probs.emplace_back(Pos_Rate);
			}
		}
		Pos_Rate.first = N;
		Pos_Rate.second = cum_rate;
		Hopp_Probs.emplace_back(Pos_Rate);
		return Hopp_Probs;
	}
	std::vector<std::pair<int, double>> DisorderedSystem::Calc_Hopping_Rates_NA_BC_IA(int i)
	{
		// Returns  vector of pairs <j-trap pos, hopping rate> to hop 
		// from trap index_<i> to trap index_<j> , NO-ACTIVATION [NA] allowed
		std::pair<int, double> Pos_Rate;
		std::vector<std::pair<int, double> > Hopp_Probs;
		int cur_pos = i;
		double cum_rate = 0;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			if (LS[_LS_count].Get_Occupied_Flag()) continue;
			else
			{
				Pos_Rate.first = LS[_LS_count].Get_Index();
				Pos_Rate.second = Calc_Hopping_Rate_BC_IA(cur_pos, LS[_LS_count].Get_Index());
				cum_rate += Pos_Rate.second;
				if (Pos_Rate.second != 0) Hopp_Probs.emplace_back(Pos_Rate);
			}
		}
		Pos_Rate.first = N;
		Pos_Rate.second = cum_rate;
		Hopp_Probs.emplace_back(Pos_Rate);
		return Hopp_Probs;
	}
	std::vector<std::pair<int, double>> DisorderedSystem::Calc_Hopping_Rates_ZeroT(int i)
	{
		// Returns  vector of pairs <j-trap pos, hopping rate> to hop 
		// from trap <i> to trap <j>, Zero-T case
		std::pair<int, double> Pos_Rate;
		std::vector<std::pair<int, double> > Hopp_Probs;
		int cur_pos = i;
		double cum_rate = 0;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			if ((_LS_count != cur_pos) & (LS[_LS_count].Get_Occupied_Flag())) continue;
			else
			{
				Pos_Rate.first = _LS_count;
				Pos_Rate.second = Calc_Hopping_Rate_ZeroT(cur_pos, _LS_count);
				cum_rate += Pos_Rate.second;
				if (Pos_Rate.second != 0) Hopp_Probs.emplace_back(Pos_Rate);
			}
		}
		Pos_Rate.first = N;
		Pos_Rate.second = cum_rate;
		Hopp_Probs.emplace_back(Pos_Rate);
		return Hopp_Probs;
	}
	std::vector<std::pair<int, double>> DisorderedSystem::Calc_Hopping_Rates_ZeroT_BC(int i)
	{
		// Returns  vector of pairs <j-trap pos, hopping rate> to hop 
		// from trap <i> to trap <j>, Zero-T case
		std::pair<int, double> Pos_Rate;
		std::vector<std::pair<int, double> > Hopp_Probs;
		int cur_pos = i;
		double cum_rate = 0;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			if ((_LS_count != cur_pos) & (LS[_LS_count].Get_Occupied_Flag())) continue;
			else
			{
				Pos_Rate.first = _LS_count;
				Pos_Rate.second = Calc_Hopping_Rate_ZeroT_BC(cur_pos, _LS_count);
				cum_rate += Pos_Rate.second;
				if (Pos_Rate.second != 0) Hopp_Probs.emplace_back(Pos_Rate);
			}
		}
		Pos_Rate.first = N;
		Pos_Rate.second = cum_rate;
		Hopp_Probs.emplace_back(Pos_Rate);
		return Hopp_Probs;
	}
	std::vector<std::pair<int, double>> DisorderedSystem::Calc_Hopping_Rates_ZeroT_IA(int i)
	{
		// Returns  vector of pairs <j-trap pos, hopping rate> to hop 
		// from trap index_<i> to trap index_<j> [IA], Zero-T case
		std::pair<int, double> Pos_Rate;
		std::vector<std::pair<int, double> > Hopp_Probs;
		int cur_pos = i;
		double cum_rate = 0;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			if ((_LS_count != cur_pos) & (LS[_LS_count].Get_Occupied_Flag()))  continue;
			else
			{
				Pos_Rate.first = LS[_LS_count].Get_Index();
				Pos_Rate.second = Calc_Hopping_Rate_ZeroT_IA(cur_pos, LS[_LS_count].Get_Index());
				cum_rate += Pos_Rate.second;
				if (Pos_Rate.second != 0) Hopp_Probs.emplace_back(Pos_Rate);
			}
		}
		Pos_Rate.first = N;
		Pos_Rate.second = cum_rate;
		Hopp_Probs.emplace_back(Pos_Rate);
		return Hopp_Probs;
	}
	std::vector<std::pair<int, double>> DisorderedSystem::Calc_Hopping_Rates_ZeroT_BC_IA(int i)
	{
		// Returns  vector of pairs <j-trap pos, hopping rate> to hop 
		// from trap index_<i> to trap index_<j> [IA], Zero-T case
		std::pair<int, double> Pos_Rate;
		std::vector<std::pair<int, double> > Hopp_Probs;
		int cur_pos = i;
		double cum_rate = 0;
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			if ((_LS_count != cur_pos) & (LS[_LS_count].Get_Occupied_Flag()))  continue;
			else
			{
				Pos_Rate.first = LS[_LS_count].Get_Index();
				Pos_Rate.second = Calc_Hopping_Rate_ZeroT_BC_IA(cur_pos, LS[_LS_count].Get_Index());
				cum_rate += Pos_Rate.second;
				if (Pos_Rate.second != 0) Hopp_Probs.emplace_back(Pos_Rate);
			}
		}
		Pos_Rate.first = N;
		Pos_Rate.second = cum_rate;
		Hopp_Probs.emplace_back(Pos_Rate);
		return Hopp_Probs;
	}

	void DisorderedSystem::Precalculate_Hopping_Rates_NA_BC()
	{
	// Pre-calculating all the hopping rates for the system	
		for (int i = 0; i < N; i++)
		{
			Hopping_Rates_Table.emplace_back(Calc_Hopping_Rates_NA_BC(i));
			std::sort(Hopping_Rates_Table[i].begin(), Hopping_Rates_Table[i].end()-1, [](auto &left, auto &right) {
				return left.second > right.second;
			});
		}
	}
	void DisorderedSystem::Precalculate_Hopping_Probs()
	{
		// Pre-calculating all the hopping probabilities for the system	
		Hopping_Probs_Table.resize(Hopping_Rates_Table.size());
		for (int i = 0; i < N; i++)
		{
			double cum_rate = Hopping_Rates_Table[i][Hopping_Rates_Table[i].size() - 1].second;
			Hopping_Probs_Table[i].resize(Hopping_Rates_Table[i].size() - 1);
			Hopping_Probs_Table[i][0].first = Hopping_Rates_Table[i][0].first;
			Hopping_Probs_Table[i][0].second = Hopping_Rates_Table[i][0].second / cum_rate;
			for (int j = 1; j < Hopping_Probs_Table[i].size(); j++)
			{
				Hopping_Probs_Table[i][j].first = Hopping_Rates_Table[i][j].first;
				Hopping_Probs_Table[i][j].second = Hopping_Probs_Table[i][j - 1].second + Hopping_Rates_Table[i][j].second / cum_rate;
			}
		}
	}
	
	void DisorderedSystem::Acticvate_Carrier(int i)
	{
		// Carrier <i> activation 
		Loc_Carr_Num--;
		Free_Carr_Num++;
		int indx = Carriers[i].Get_Loc_State_Index();
		Carriers[i].Set_Localized_Flag(false);
		LS[indx].Set_Occupied_Flag(false);
		LS_Occupation_Status[indx] = false;
		LS[indx].Set_Loc_Carrier_index(-1);	
	}
	void DisorderedSystem::Acticvate_Carrier_IA(int i)
	{
		// Carrier index_<i> activation 
		Loc_Carr_Num--;
		Free_Carr_Num++;
		int indx = Carriers_index[i]->Get_Loc_State_Index();
		Carriers_index[i]->Set_Localized_Flag(false);
		LS[indx].Set_Occupied_Flag(false);
		LS_Occupation_Status[indx] = false;
		LS[indx].Set_Loc_Carrier_index(-1);	
	}
	void DisorderedSystem::Capture_Carrier(int i, int j)
	{
		// Capture carrier <i> to trap <j> 
		Free_Carr_Num--;
		Loc_Carr_Num++;
		Carriers[i].Set_Localized_Flag(true);
		Carriers[i].Set_Loc_State_Index(j);
		LS[j].Set_Occupied_Flag(true);
		LS_Occupation_Status[j] = true;
		LS[j].Set_Loc_Carrier_index(i);
	}
	void DisorderedSystem::Capture_Carrier_IA(int i, int j)
	{
		// Capture carrier index_<i> to trap inedx_<j> 
		Free_Carr_Num--;
		Loc_Carr_Num++;
		int j_indx = Get_LS_Position(j);
		Carriers_index[i]->Set_Localized_Flag(true);
		Carriers_index[i]->Set_Loc_State_Index(j_indx);
		LS_index[j]->Set_Occupied_Flag(true);
		LS_Occupation_Status[j_indx] = true;
		LS_index[j]->Set_Loc_Carrier_index(Get_Carriers_Position(i));
	}
	void DisorderedSystem::Move_Carrier(int Carr_i, int LS_j)
	{
		// Move carrier <i> to trap <j> 
		int LS_indx = Carriers[Carr_i].Get_Loc_State_Index();
		Carriers[Carr_i].Set_Loc_State_Index(LS_j);
		LS[LS_indx].Set_Occupied_Flag(false);
		LS_Occupation_Status[LS_indx] = false;
		LS[LS_indx].Set_Loc_Carrier_index(-1);
		LS[LS_j].Set_Occupied_Flag(true);
		LS_Occupation_Status[LS_j] = true;
		LS[LS_j].Set_Loc_Carrier_index(Carr_i);
	}
	void DisorderedSystem::Move_Carrier_IA(int Carr_i, int LS_j)
	{
		// Move carrier index_<i> to LS index_<j> 
		int LS_indx = Carriers_index[Carr_i]->Get_Loc_State_Index();
		Carriers_index[Carr_i]->Set_Loc_State_Index(LS_j);
		LS_index[LS_indx]->Set_Occupied_Flag(false);
		LS_Occupation_Status[Get_LS_Position(Carr_i)] = false;
		LS_index[LS_indx]->Set_Loc_Carrier_index(-1);
		LS_index[LS_j]->Set_Occupied_Flag(true);
		LS_Occupation_Status[Get_LS_Position(LS_j)] = true;
		LS_index[LS_j]->Set_Loc_Carrier_index(Get_Carriers_Position(Carr_i));
	}

	void DisorderedSystem::Set_Static_Fermi_Occupation()
	{
	// Marks LS occupied according to Fermi distribution
		std::random_device rd;
		std::mt19937 generator(rd());
		std::uniform_real_distribution <double> uni_distribution(0, 1);
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			double p = uni_distribution(generator);
			if (p <= 1 / (exp(LS[_LS_count].Get_E() - _E_Fermi) + 1))
			{
				LS[_LS_count].Set_Occupied_Flag(true);
				LS_Occupation_Status[_LS_count] = true;
			}
		}
	}
	void DisorderedSystem::Set_Static_Fermi_Occupation_ET()
	{
	// Marks LS occupied according to Fermi distribution 
	//with effective temperature
		std::random_device rd;
		std::mt19937 generator(rd());
		std::uniform_real_distribution <double> uni_distribution(0, 1);
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			double p = uni_distribution(generator);
			if (p <= 1 / (exp((LS[_LS_count].Get_E() - _E_Fermi) / T_Eff_Factor) + 1))
			{
				LS[_LS_count].Set_Occupied_Flag(true);
				LS_Occupation_Status[_LS_count] = true;
			}
		}
	}
	void DisorderedSystem::Set_Fermi_Occupation_Probs()
	{
	// Calculates Probailities of LS being occupied according to Fermidistribution
		LS_Fermi_Occupation_Probs.resize(N);
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			LS[_LS_count].Set_Occupation_Probability(1 / (exp(LS[_LS_count].Get_E() - _E_Fermi) + 1));
			LS_Fermi_Occupation_Probs[_LS_count] = 1 / (exp(LS[_LS_count].Get_E() - _E_Fermi) + 1);
		}

	}
	void DisorderedSystem::Set_Fermi_Occupation_Probs_ET()
	{
	// Calculates Probailities of LS being occupied according to Fermidistribution 
	// with effective temperature
		LS_Fermi_Occupation_Probs.resize(N);
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			LS[_LS_count].Set_Occupation_Probability(1 / ((exp(LS[_LS_count].Get_E() - _E_Fermi) / T_Eff_Factor) + 1));
			LS_Fermi_Occupation_Probs[_LS_count] = 1 / ((exp(LS[_LS_count].Get_E() - _E_Fermi) / T_Eff_Factor) + 1);
		}
	}
	void DisorderedSystem::Reset_Static_Occupation()
	{
		// Drops all LS <is_occupied> to false 
		// for trpas previously marked as occupied
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
		{
			if ((LS[_LS_count].Get_Occupied_Flag()) & (LS[_LS_count].Get_Loc_Carrier_index() == -1))
			{
				LS[_LS_count].Set_Occupied_Flag(true);
				LS_Occupation_Status[_LS_count] = false;
			}
		}
	}
	
	// OUTPUT METHODS
	void DisorderedSystem::Output_Energies(std::string filename)
	{
		// Output LS energies distribution (DOS) to file
		std::ofstream output_file(filename);
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
			output_file << LS[_LS_count].Get_E() << std::endl;
		output_file.close();
	}
	void DisorderedSystem::Output_Occupied_Energies(std::string filename)
	{
		// Output LS energies distribution (DOS) to file
		std::ofstream output_file(filename);
		for (int _LS_count = 0; _LS_count < N; _LS_count++)
			if (LS[_LS_count].Get_Occupied_Flag()) 
				output_file << LS[_LS_count].Get_E() << std::endl;
		output_file.close();
	}

}
