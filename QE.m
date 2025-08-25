classdef QE < handle
    %QE   Parse and store Quantum ESPRESSO structure data.
    %   QE provides a container for structural and atomic data compatible with
    %   Quantum ESPRESSO, with utilities for converting from VASP POSCAR files.
    %
    %   qe = QE.fromPOSCAR('POSCAR') reads a VASP POSCAR file and returns a QE object
    %   populated with structure and atomic information.
    %
    %   qe = QE.fromXML('prefix.xml') reads a Quantum ESPRESSO XML output file and returns
    %   a QE object with the parsed geometry and atomic information.
    %
    %   This class stores cell vectors, atomic coordinates, species, and other
    %   fundamental data required for Quantum ESPRESSO calculations, and provides
    %   utility functions such as element symbol to atomic number conversion.
    %
    %   Properties (accessible after construction):
    %     prefix      - (char) System prefix or comment (from POSCAR line 1 or QE input)
    %     alat        - (double) Lattice constant in Bohr
    %     at          - (3x3 double) Lattice vectors (columns: a1, a2, a3, unit: alat)
    %     tpiba       - (double) 2*pi/alat
    %     bg          - (3x3 double) Reciprocal lattice vectors (columns: b1, b2, b3, unit: tpiba)
    %     atm         - (cellstr) Names of atomic species
    %     amass       - (vector) Mass of atomic species
    %     atnum       - (vector) Atomic number of each atom
    %     na          - (row vector) Number of atoms of each species
    %     nsp         - (int) Number of species
    %     nat         - (int) Total number of atoms
    %     ityp        - (row vector) Species index for each atom (1...nsp)
    %     tau         - (3xnat double) Atomic positions (Bohr)
    %     if_pos      - (3xnat logical) Whether each atomic coordinate is fixed
    %     fileType    - (string) File type: 'qe' (default) or 'vasp' (if from POSCAR)
    %     poscarData  - (struct) Original POSCAR data as unpacked struct (if constructed from POSCAR)
    %
    %   Static Methods:
    %     qe = QE.fromPOSCAR(filePath)
    %         Read a VASP POSCAR file and construct a QE object.
    %
    %     qe = QE.fromXML(filePath)
    %         Read a Quantum ESPRESSO XML output file and construct a QE object.
    %
    %     atnum = QE.get_atomic_number(symbols)
    %         Return atomic number(s) for given element symbol(s).
    %
    %   Example:
    %     qe = QE.fromPOSCAR('POSCAR');
    %     disp(qe.nat);          % total number of atoms
    %     disp(qe.at);           % lattice vectors
    %
    %   For more information about Quantum ESPRESSO, see:
    %     https://www.quantum-espresso.org/
    %
    %   See also VPHMAT, MATDYN, NIST
    %
    %   Copyright (c) 2025 Lizhou Yang
    %   MIT License

    properties
        prefix             % prefix or comment of the system
        alat               % lattice constant (Bohr)
        at                 % lattice vectors (3×3, each column is a_i, unit: alat)
        tpiba              % 2π/alat
        bg                 % reciprocal lattice vectors (3×3, each column is b_i, unit: tpiba)
        atm                % names of atomic species (cell)
        amass              % mass of atomic species (vector)
        atnum              % atomic number of each atom
        na                 % number of atoms of each species (vector)
        nsp                % number of species
        nat                % total number of atoms
        ityp               % type index of each atom 1...nsp
        tau                % atomic positions (3×nat, unit: Bohr)
        if_pos             % whether position is fixed (3×nat, logical)
        fileType = "qe"    % file type
        poscarData         % original POSCAR unpacked struct
    end

    methods(Static)
        function qe = fromPOSCAR(filePath)
            % Initialize a QE object from a VASP POSCAR file
            qe = QE;
            qe.fileType = "vasp";
            poscar.rawPOSCARContent = fileread(filePath);

            fid = fopen(filePath, 'r');
            if fid == -1
                error('QE:OpenFileError', ...
                    'Cannot open this POSCAR file: %s', filePath);
            end
            cleaner = onCleanup(@() fclose(fid));

            line1 = strtrim(fgetl(fid));
            if isvarname(line1)
                qe.prefix = line1;
            else
                qe.prefix = "pwscf";
            end
            poscar.comment = line1;
            scale = fscanf(fid, ' %f ', 1);
            poscar.scalingFactor = scale;
            lat = scale * fscanf(fid, ' %f %f %f \n', [3 3]).';
            poscar.lattice = lat;
            qe.alat = 1; % alat=1 Bohr
            qe.at = lat.' / (NIST.Bohr_radius * 1e10); % convert to alat units
            qe.tpiba = 2 * pi / qe.alat;
            qe.bg = inv(qe.at.'); % reciprocal lattice

            poscar.speciesNames = split(strtrim(fgetl(fid)));
            qe.atm = poscar.speciesNames;
            poscar.ionsPerSpecies = str2num(fgetl(fid));
            qe.na = poscar.ionsPerSpecies;
            qe.nsp = numel(qe.na);
            qe.nat = sum(qe.na);

            ityp = arrayfun(@(n, idx) repmat(idx, 1, n), qe.na, 1:qe.nsp, 'UniformOutput', false);
            qe.ityp = [ityp{:}];
            atm_num=QE.get_atomic_number(qe.atm);
            qe.atnum=atm_num(qe.ityp);
            qe.atnum=qe.atnum(:);

            s = strtrim(fgetl(fid));
            if startsWith(s, 's', 'IgnoreCase', true) % Selective dynamics
                poscar.selectiveDynamics = true;
                s = strtrim(fgetl(fid));
                poscar.direct = ~startsWith(s, {'k', 'c'}, 'IgnoreCase', true);
            else
                poscar.selectiveDynamics = false;
                poscar.direct = ~startsWith(s, {'k', 'c'}, 'IgnoreCase', true);
            end

            if poscar.selectiveDynamics
                tmp = fscanf(fid, ' %f %f %f %s %s %s \n', [6 qe.nat]).';
                poscar.ionPositions = tmp(:, 1:3);
                poscar.ionPositionsDynamics = tmp(:, 4:end) == 84; % 'T'->84, 'F'->70
            else
                coords = textscan(fid, '%f %f %f %s', qe.nat);
                poscar.ionPositions = [coords{1}, coords{2}, coords{3}];
                poscar.ionPositionsDynamics = true(qe.nat, 3);
            end

            qe.if_pos = ~poscar.ionPositionsDynamics.';
            if poscar.direct
                qe.tau = (poscar.ionPositions * poscar.lattice).' / (NIST.Bohr_radius * 1e10);
            else
                qe.tau = scale * poscar.ionPositions.' / (NIST.Bohr_radius * 1e10);
            end

            qe.poscarData = poscar;
        end

        function qe = fromXML(filePath)
            % Initialize a QE object from prefix.xml
            qe = QE;
            qe.fileType = "qe";
            s=readstruct(filePath);
            a1=str2num(s.output.atomic_structure.cell.a1);
            a2=str2num(s.output.atomic_structure.cell.a2);
            a3=str2num(s.output.atomic_structure.cell.a3);
            qe.at=[a1;a2;a3].'; % Bohr
            qe.alat=s.output.atomic_structure.alatAttribute;
            qe.tpiba=2*pi/qe.alat;
            qe.nat=s.output.atomic_structure.natAttribute;
            qe.at=qe.at/qe.alat; % alat
            b1=str2num(s.output.basis_set.reciprocal_lattice.b1);
            b2=str2num(s.output.basis_set.reciprocal_lattice.b2);
            b3=str2num(s.output.basis_set.reciprocal_lattice.b3);
            qe.bg=[b1;b2;b3].'; % reciprocal lattice, tpiba
            qe.atm=cellstr([s.input.atomic_species.species.nameAttribute])';
            qe.amass=[s.input.atomic_species.species.mass].';
            qe.prefix=s.input.control_variables.prefix;
            sp_num=QE.get_atomic_number(qe.atm);
            tau=zeros(qe.nat,3);
            atnum=zeros(qe.nat,1);
            for ia=1:qe.nat
                tau(ia,:)=str2num(s.output.atomic_structure.atomic_positions.atom(ia).Text);
                atnum(ia)=QE.get_atomic_number(s.output.atomic_structure.atomic_positions.atom(ia).nameAttribute);
            end
            qe.tau=tau.';
            qe.atnum=atnum;
            tmphash=zeros(max(sp_num),1);
            tmphash(sp_num)=1:size(sp_num,1);
            qe.ityp=tmphash(qe.atnum);
            qe.nsp=max(qe.ityp);
            tmphash=zeros(1,qe.nsp);
            for ia=1:qe.nat
                tmphash(qe.ityp(ia))=tmphash(qe.ityp(ia))+1;
            end
            qe.na=tmphash;
        end
    end

    methods
        function this = QE()
            % Empty constructor
        end
    end

    methods(Static)
        function atnum = get_atomic_number(symbols)
            % get_atomic_number Returns atomic number(s) given element symbol(s)
            % Supports string, char, cellstr, or string array as input
            persistent symbol2num
            if isempty(symbol2num)
                elems = { ...
                    'H', 1;   'He', 2;  ...
                    'Li', 3;  'Be', 4;  'B', 5;  'C', 6;  'N', 7;  'O', 8;  'F', 9;  'Ne', 10; ...
                    'Na', 11; 'Mg', 12; 'Al', 13; 'Si', 14; 'P', 15; 'S', 16; 'Cl', 17; 'Ar', 18; ...
                    'K', 19;  'Ca', 20; 'Sc', 21; 'Ti', 22; 'V', 23; 'Cr', 24; 'Mn', 25; 'Fe', 26; ...
                    'Co', 27; 'Ni', 28; 'Cu', 29; 'Zn', 30; 'Ga', 31; 'Ge', 32; 'As', 33; 'Se', 34; 'Br', 35; 'Kr', 36; ...
                    'Rb', 37; 'Sr', 38; 'Y', 39;  'Zr', 40; 'Nb', 41; 'Mo', 42; 'Tc', 43; 'Ru', 44; 'Rh', 45; 'Pd', 46; ...
                    'Ag', 47; 'Cd', 48; 'In', 49; 'Sn', 50; 'Sb', 51; 'Te', 52; 'I', 53;  'Xe', 54; ...
                    'Cs', 55; 'Ba', 56; 'La', 57; 'Ce', 58; 'Pr', 59; 'Nd', 60; 'Pm', 61; 'Sm', 62; 'Eu', 63; 'Gd', 64; ...
                    'Tb', 65; 'Dy', 66; 'Ho', 67; 'Er', 68; 'Tm', 69; 'Yb', 70; 'Lu', 71; ...
                    'Hf', 72; 'Ta', 73; 'W', 74;  'Re', 75; 'Os', 76; 'Ir', 77; 'Pt', 78; 'Au', 79; 'Hg', 80; ...
                    'Tl', 81; 'Pb', 82; 'Bi', 83; 'Po', 84; 'At', 85; 'Rn', 86; ...
                    'Fr', 87; 'Ra', 88; 'Ac', 89; 'Th', 90; 'Pa', 91; 'U', 92;  'Np', 93; 'Pu', 94; 'Am', 95; ...
                    'Cm', 96; 'Bk', 97; 'Cf', 98; 'Es', 99; 'Fm', 100; 'Md', 101; 'No', 102; 'Lr', 103; ...
                    'Rf', 104; 'Db', 105; 'Sg', 106; 'Bh', 107; 'Hs', 108; 'Mt', 109; 'Ds', 110; 'Rg', 111; ...
                    'Cn', 112; 'Fl', 114; 'Lv', 116; 'Ts', 117; 'Og', 118; ...
                    };
                symbol2num = containers.Map(elems(:,1), cell2mat(elems(:,2)));
            end

            if ischar(symbols) || isstring(symbols)
                symbols = cellstr(symbols);
            end

            atnum = zeros(size(symbols));
            for i = 1:numel(symbols)
                symb = char(symbols{i});
                if isKey(symbol2num, symb)
                    atnum(i) = symbol2num(symb);
                else
                    error('Unknown element symbol: %s', symb);
                end
            end
        end
    end

end
