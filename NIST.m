classdef NIST
    %NIST Fundamental physical constants.
    % For more information about NIST physical constants, visit: https://physics.nist.gov/cuu/Constants/
    %
    % See Also QE, VPHMAT, MATDYN

    properties (Constant)
        elementary_charge = 1.602176634e-19; %C
        electron_mass = 9.1093837015e-31; %kg, atomic unit of mass
        proton_mass = 1.67262192595e-27; %kg, m_p
        deuteron_mass = 3.3435837768e-27; %kg, m_d
        reduced_Planck_constant = 1.054571817e-34; %J s
        Planck_constant = 6.62607015e-34; %J Hz^{-1}
        Boltzmann_constant = 1.380649e-23; %J K^{-1}
        Boltzmann_constant_in_eV_K = 8.617333262e-5; % eV/K
        Bohr_radius = 5.29177210903e-11; %m
        atomic_mass_constant = 1.66053906660e-27; %kg
        Rydberg_constant = 10973731.568157; %m^{-1}
        Rydberg_constant_in_eV = 13.605693122990; %eV, Rydberg constant times hc in eV
        Rydberg_constant_in_J = 2.1798723611030e-18; %J =electron_mass*elementary_charge^4/(8*vacuum_electric_permittivity^2*Plank_constant^2)
        Hartree_energy = 4.3597447222060e-18; %J
        Hartree_energy_in_eV = 27.211386245981; %eV
        atomic_unit_of_time = 2.4188843265864e-17; %s, \hbar/E_h
        speed_of_light = 299792458; %m/s
        vacuum_electric_permittivity = 8.8541878188e-12; %F/m
        vacuum_magnetic_permeability = 1.25663706127e-6; %N A^{-2} \mu_0
        Newtonian_constant_of_gravitation = 6.67430e-11; %m^3 kg^{-1} s^{-2}
        Avogadro_constant=6.02214076e23; %mol^{-1} N_A
    end

    properties (Access = private, Constant)
        convmat=NIST.initconvmat; % consistent judgement matrix
    end

    methods (Static)
        function coef = convertunit(dataunit, unit)
            % unit can be 'mev','ev','ry','thz','cmm1'
            unit2ind = containers.Map({'mev','ev','ry','thz','cmm1'}, 1:5);
            convmat = NIST.convmat;
            coef = convmat(unit2ind(dataunit), unit2ind(unit));
        end
    end

    methods (Access=private, Static)
        function convmat = initconvmat()
            n = 5;
            convmat = ones(n, n);
            convmat(1,2) = 1e-3; % mev_to_ev
            convmat(2,3) = 1/NIST.Rydberg_constant_in_eV; % ev_to_ry
            convmat(3,4) = 1e-12*NIST.Rydberg_constant_in_J/NIST.Planck_constant; % ry_to_thz
            convmat(4,5) = 1e10/NIST.speed_of_light; % thz_to_cmm1

            for i = 1:n
                prod_val = 1;
                for j = i+1:n
                    prod_val = prod_val * convmat(j-1, j);
                    convmat(i, j) = prod_val;
                end
            end
            for i = 2:n
                for j = 1:i-1
                    convmat(i, j) = 1 / convmat(j, i);
                end
            end
        end
    end
end