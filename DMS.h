#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <iterator>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>  
#include <functional> 
#include <iterator>


namespace DMS
{
	// PHYSICAL CONSTANTS
	const double k	 = 8.625e-5;		// Boltzman constant [ ~ Ev/K ]	
	const double a_B = 5.3e-9;		    // Bohr radius		 [ ~ cm   ]
	const double Ry  = 13.6;			// Rydberg energy	 [ ~ Ev   ]	

	// LOCALIZED STATE CLASS DECLARATION
	class LocState
	{
	// Localized State [LS] Class
	private:
		// LS parameters 
		bool is_occupied;		// flag to check if the LS is occupied by carrier [ LS is occupied if "True"									]
		bool is_radiative;		// flag determines if the LS is of radiative type [ LS is radiative if "True"									] 
								// (for models which assume radiative and nonradiative channels of exciton recombination)
		
		int Loc_Carrier_index;	// index(POSSITION) of the carrier at the LS      [ -1 if LS is not occuoied									]
		int index;				// index of the LS
		int type;				// type of the LS
		
		double X;				// x - coordinate
		double Y;				// y - coordinate
		double Z;				// z - coordinate
		double E;				// energy potential of the LS
		double Loc_Rad;			// localization radius of the LS				  [ possible to set various Loc_Rad different from Sys_Loc_Rad  ]
								//												  [ ~ Bohr - radius units                                       ]
		double occupation_prob; // probability of being occupied                  [ from 0 to 1                                                 ]
	
	public:
		// default constructor
		LocState();
		// parametrized constructor
		LocState(int LS_index, int LS_type, double LS_x_coord, double LS_y_coord, double LS_z_coord, 
					double LS_energy, double LS_Loc_Rad, bool occupation_flag, bool radiative_flag);

		// Set-methods
		void Set_Index(int LS_index);
		void Set_Type(int LS_type);
		void Set_X(double LS_x_coord);
		void Set_Y(double LS_y_coord);
		void Set_Z(double LS_z_coord);
		void Set_E(double LS_energy);
		void Set_Loc_Rad(double LS_Loc_Rad);
		void Set_Occupied_Flag(bool occupation_flag);
		void Set_Occupation_Probability(double occupation_probability);
		void Set_Loc_Carrier_index(int LC_index);
		void Set_Radiative_Flag(bool radiative_flag);
		
		// Get-methods
		int Get_Index() const noexcept;
		int Get_Type() const noexcept;
		double Get_X() const noexcept;
		double Get_Y() const noexcept;
		double Get_Z() const noexcept;
		double Get_E() const noexcept;
		double Get_Loc_Rad() const noexcept;
		bool Get_Occupied_Flag() const noexcept;
		double Get_Occupation_Probability() const noexcept;
		int Get_Loc_Carrier_index() const noexcept;
		bool Get_Radiative_Flag() const noexcept;
	};

	// CARRIER CLASS DECLARATION
	class Carrier
	{
	// Carrier Class
	private:
		// Carrier parameters
		bool is_localized;			// localization flag								[ "True" if carrier is localized             ]
		int index;				    // carriers index					
		int type;					// type												[ ~ 0 = electron, 1 = exciton                ]
		int loc_state_index;		// position of the LS								[ carrier localized at this LS				 ] 
									//													[ = -1 for delocalized carrier				 ]	
		int hops_counter;			// to count total number of hops made by carrier
		double life_time;			// carriers lifetime(for esxitons)					[ ~ Nu0 * Tau				                 ]
		double timer;				// timer to count  time since carrier was generated	[ ~ sec										 ]

	public:
		// default constructor
		Carrier();
		// parametrizedconstructor
		Carrier(int Carr_index, int Carr_type, int C_LS_index, double Carr_life_time, bool localized_flag);

		// Set-methods
		void Set_Index(int Carr_index);
		void Set_Type(int Carr_type);
		void Set_Loc_State_Index(int C_LS_index);
		void Set_Life_Time(double Carr_life_time);
		void Set_Localized_Flag(bool localized_flag);
		void Set_Timer_Value(double Carr_time);
		void Set_Hops_Counter(int Carr_hops_counter);

		// Get-methods
		int Get_Index() const noexcept;
		int Get_Type() const noexcept;
		int Get_Loc_State_Index() const noexcept;
		double Get_Life_Time() const noexcept;
		bool Get_Localized_Flag() const noexcept;
		double Get_Timer_Value() const noexcept;
		int Get_Hops_Counter() const noexcept;

		// Carrier-timer methods
		void Increase_Timer(double Carr_timer_incr);
		void Decrease_Timer(double Carr_timer_decr);

		// Hops-counter methods
		void Increase_Hops_Counter(int Carr_hops_counter_incr);
		void Decrease_Hops_Counter(int Carr_hops_counter_decr);
	};

	// DISORDERED SYSTEM CLASS DECLARATION 
	class DisorderedSystem
	{
	// Disordered system [DS] class
	private:
		// Disordered system parameters
		int Dim;							// sys's dimensionality
		int N;								// total number of traps(LS's) in the system
		int Loc_Carr_Num;					// number of localized carriers in the system
		int Free_Carr_Num;					// number of free carriers in the system
		
		int LS_index_counter;				// counter for LS indexes									[ counter - 1 = max index                    ]
		int carrier_index_counter;			// counter for carriers indexes								[ counter - 1 = max index					 ]

		double E0;							// first characteristic scale of disorder potential			[ for exp or Gauss density of states         ]
		double E1;							// second characteristic scale of disorder potential		[ for complex density of states				 ]
		double Traps_Conc;					// concentration of traps per loc-rad-volume				[ N0 x (loc_rad ^ dim) ~ dimensionless		 ]
		double Sys_Size;					// characteristic length of cubic cell with LS's			[ ~ localization radius units				 ]
		double Sys_Loc_Rad;					// localization radius at LS's of the system				[ ~ Bohr-radius units						 ]
		double Cutoff_Rad;					// cutoff radius for hopping transitions					[ ~ Bohr-radius units						 ]

		double Niu0;						// attempt to escape frequency (= capture rate)             [ ~ sec ^ -1								 ]
		
		double T;							// system temperature										[ ~ K                                        ]			
		double kT;							// system temperature * Boltzman's constant					[ ~ Ev                                       ]
		double F;							// electric field applied to system							[ ~ V/cm									 ]
		
		double T_Eff_Factor;				// Effective temperature factor                             [ T_eff = T_Eff_Factor * T (real)			 ]
		double T_eff;						// Effective temperature value                              [ ~ K                                        ]
		double kT_eff;                      // Effective temperature * Boltzman's constant              [ ~ Ev                                       ]
		
		double E_Fermi;						// Fermi level											    [ ~ Ev										 ]
		double _E_Fermi;					// Dimensionless Fermi level                                [ E_Fermi / kT								 ]

		double E_Transport;					// Transport energy                                         [ ~ Ev                                       ] 
		double _E_Transport;				// Dimensionless transport energy                           [ E_Transport / kT                           ]

		double timer;						// timer to count global system time						[ ~ sec										 ]
		
		// Dimensionless parameters
		double _E0;							// dimesnionless E0											[  E0 / kT									 ]
		double _E1;							// dimensionless E1											[  E1 / kT                                   ]	
		double _Traps_Conc;					// dimensionless concentrarion of traps						[  N0 x (loc_rad ^ dim)						 ]
		double _F;							// dimensionless electric field								[  e x F x (N0 ^ 1/dim)/E0   			     ]

	public:
		std::vector<DMS::LocState> LS;						// vector to store LS's
		std::map<int, DMS::LocState*> LS_index;				// vector of pointers to LS's				[ to access by index						 ]
		std::vector<DMS::Carrier> Carriers;					// vector to store Carriers
		std::map<int, DMS::Carrier*> Carriers_index;		// vector to store pointers to Carriers		[ to access by index						 ]	
		std::vector<bool> LS_Occupation_Status;				// vector to store LS occupation status
		std::vector<double> LS_Fermi_Occupation_Probs;		// vector to store LS occupation probabilities
		std::vector<std::vector<std::pair<int, double>>> Hopping_Rates_Table; // vector to store precalculated hopping rates
		std::vector<std::vector<std::pair<int, double>>> Hopping_Probs_Table; // vector to store precalculated hopping probabilities
		
        // default constructor
		DisorderedSystem();
		// parametrized constructor	[default 3d-cubic case / Gaussian DOS]
		DisorderedSystem(int LS_Num, double DS_energy_scale_1 , double DS_Loc_Rad, double DS_Cutoff,
						double DS_Traps_Concentration, double DS_Temp, double DS_Field);
		// parametrized constructor	[default 3d-cubic case / Gaussian DOS/ dimesionless parameters]
		DisorderedSystem(int LS_Num, double DS_energy_scale_1, double DS_Traps_Concentration, 
						 double DS_Cutoff, double DS_Field);

		// Set-methods
		void Set_Dim(int DS_dimensions_num);
		void Set_N(int LS_Num);
		void Set_E0(double DS_energy_scale_exp);
		void Set_Sys_Size(double DS_size);
		void Set_Sys_Loc_Rad(double DS_Loc_Rad);
		void Set_Cutoff_Rad(double DS_Cutoff);
		void Set_Traps_Conc(double DS_Traps_concentration);
		void Set_Escape_Rate(double LS_Escape_rate);
		void Set_Sys_Temperature(double DS_Temperature);
		void Set_kT(double DS_kT);
		void Set_Field(double DS_Field);
		void Set_Timer_Value(double DS_time);
		void Set_dls_E0(double DS_energy_scale_exp);
		void Set_dls_Traps_Conc(double DS_Traps_concentration);
		void Set_dls_Field(double DS_Field);
		void Set_E_Fermi(double DS_Fermi_level);
		void Set_dls_E_Fermi(double _DS_Fermi_level);

		// Get-methods
		int Get_Dim() const noexcept;
		int Get_N() const noexcept;
		double Get_E0() const noexcept;
		double Get_Sys_Size() const noexcept;
		double Get_Sys_Loc_Rad() const noexcept;
		double Get_Cutoff_Rad() const noexcept;
		double Get_Traps_Conc() const noexcept;
		int Get_Loc_Carr_Num() const noexcept;
		int Get_Free_Carr_Num() const noexcept;
		double Get_Escape_Rate() const noexcept;
		double Get_Sys_Temperature() const noexcept;
		double Get_kT() const noexcept;
		double Get_Field() const noexcept;
		double Get_Timer_Value() const noexcept;
		double Get_dls_E0() const noexcept;
		double Get_dls_Traps_Conc() const noexcept;
		double Get_dls_Field() const noexcept;
		double Get_E_Fermi() const noexcept;
		double Get_dls_E_Fermi() const noexcept;
		double Get_T_Eff_Factor() const noexcept;
		double Get_T_Eff() const noexcept;
		double Get_kT_Eff_Factor() const noexcept;
		double Get_E_Transport() const noexcept;
		double Get_dls_E_Transport() const noexcept;

		// Effective temperature calculation
		void Calc_Effective_T(double gamma);

		// Transport Energy Calculation
		void Calc_Transport_Energy_ET(double precision);
		
		// System-timer methods
		void Increase_Timer(double DS_timer_incr);
		void Decrease_Timer(double DS_timer_decr);

		// Energies end coordinates distribution methods
		void Gen_Coords_Distrib_3D_Random();
		void Gen_Coords_Distrib_3D_Cubic_Lattice();
		void Gen_Energies_Distrib_Exp();
		void Gen_Energies_Distrib_Gauss();

		// Index access methods
		void Set_LS_Index_Ptrs();
		int  Get_LS_Position(int _LS_index);
		void Set_Carriers_Index_Ptrs();
		int Get_Carriers_Position(int _Carr_index);

		// Add/Remove traps and carriers
		// here and further if IA (index_access) -> traps are vec[i]/vec[j] 
		// ,else -> traps are vec_index[i]/vec_index[j] 
		// where vec - vector of systems object (LS, Carriers , etc... )
		// ..._IA methods for Idex Access!
		void Add_Trap(int LS_type, double LS_x_coord, double LS_y_coord, double LS_z_coord, 
						double LS_energy, double LS_Loc_Rad, bool radiative_flag);
		void Remove_Trap(int LS_index);
		void Remove_Trap_IA(int LS_index);
		void Add_Carrier(int Carr_type, int C_LS_index, double Carr_life_time, bool localized_flag);
		void Add_Carrier_IA(int Carr_type, int C_LS_index, double Carr_life_time, bool localized_flag);
		void Remove_Carrier(int Carr_index);
		void Remove_Carrier_IA(int Carr_index);

		// LS comparison methods
		double Calc_Distance(int i, int j);
		double Calc_Distance_BC(int i, int j);
		double Calc_Energy_Diff(int i, int j);

		// Calculate hopping rates methods											
		// assume i , j are the initial and final sites
		// ..._NA-methods to suppress activation
		// ..._BC-methods for enabled boundatirs conditions ( infinit system ) 
		// ..._IA-methods for Idex Access
		double Calc_Hopping_Rate(int i, int j);
		double Calc_Hopping_Rate_BC(int i, int j);
		double Calc_Hopping_Rate_IA(int i, int j);
		double Calc_Hopping_Rate_BC_IA(int i, int j);
		double Calc_Hopping_Rate_ZeroT(int i, int j);
		double Calc_Hopping_Rate_ZeroT_BC(int i, int j);
		double Calc_Hopping_Rate_ZeroT_IA(int i, int j);
		double Calc_Hopping_Rate_ZeroT_BC_IA(int i, int j);

		std::vector<std::pair<int, double>> Calc_Hopping_Rates(int i);
		std::vector<std::pair<int, double>> Calc_Hopping_Rates_BC(int i);
		std::vector<std::pair<int, double>> Calc_Hopping_Rates_NA(int i);
		std::vector<std::pair<int, double>> Calc_Hopping_Rates_NA_BC(int i);
		std::vector<std::pair<int, double>> Calc_Hopping_Rates_IA(int i);
		std::vector<std::pair<int, double>> Calc_Hopping_Rates_BC_IA(int i);
		std::vector<std::pair<int, double>> Calc_Hopping_Rates_NA_IA(int i);
		std::vector<std::pair<int, double>> Calc_Hopping_Rates_NA_BC_IA(int i);
		std::vector<std::pair<int, double>> Calc_Hopping_Rates_ZeroT(int i);
		std::vector<std::pair<int, double>> Calc_Hopping_Rates_ZeroT_BC(int i);
		std::vector<std::pair<int, double>> Calc_Hopping_Rates_ZeroT_IA(int i);
		std::vector<std::pair<int, double>> Calc_Hopping_Rates_ZeroT_BC_IA(int i);

		void Precalculate_Hopping_Rates_NA_BC();
		void Precalculate_Hopping_Probs();

		// Carrier transition methods - activation/capture/hopping
		void Acticvate_Carrier(int i);			// here i is carrier position/index
		void Acticvate_Carrier_IA(int i);
		void Capture_Carrier(int i, int j);		// here i is carriers position/index and j is is LS position/index 	
		void Capture_Carrier_IA(int i, int j);
		void Move_Carrier(int i, int j);		// i , j are the initial and final sites position/index				
		void Move_Carrier_IA(int i, int j);

		// Occupation methods
		// ET - effective temperature
		void Set_Static_Fermi_Occupation();		// staticaly sets random LS to be occupied according to Fermi distribution
		void Set_Static_Fermi_Occupation_ET();  // staticaly sets random LS to be occupied according to Fermi distribution with effective temperature [ET]
		void Set_Fermi_Occupation_Probs();		// sets probablities of being occupied for LS according Fermi distribution
		void Set_Fermi_Occupation_Probs_ET();   // sets probablities of being occupied for LS according Fermi distribution with effective temperature [ET]
		void Reset_Static_Occupation();			// resets static occupation


		// Output Energies distribution (DOS) to file
		void Output_Energies(std::string filename);
		void Output_Occupied_Energies(std::string filename);
	};

}
