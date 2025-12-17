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
        zv
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

        function toPOSCAR(alat,at,tau,atnum,filePath)
            fid=fopen(filePath,'w');
            if fid == -1
                error('QE:OpenFileError', ...
                    'Cannot open this POSCAR file: %s', filePath);
            end
            cleaner = onCleanup(@() fclose(fid));
            % fprintf(fid,['Generated by VPHMAT\n %.16f\n' ...
            %     '%.16f %.16f %.16f\n' ...
            %     '%.16f %.16f %.16f\n' ...
            %     '%.16f %.16f %.16f\n'], ...
            %     alat*NIST.Bohr_radius*1e10,at);
            fprintf(fid,['Generated by VPHMAT\n %.16f\n' ...
                '%.16f %.16f %.16f\n' ...
                '%.16f %.16f %.16f\n' ...
                '%.16f %.16f %.16f\n'], ...
                1.0,alat*NIST.Bohr_radius*1e10*at);
            [atm_num,~,ityp]=unique(atnum,'stable');
            [~,I]=sort(ityp);
            atm=QE.get_atomic_sym(atm_num);
            nat=numel(atnum);
            nsp=numel(atm_num);
            tmphash=zeros(nsp,1);
            for ia=1:nat
                tmphash(ityp(ia))=tmphash(ityp(ia))+1;
            end
            na=tmphash;
            fprintf(fid,'%s ',atm{:});
            fprintf(fid,'\n');
            fprintf(fid,'%d ',na);
            fprintf(fid,'\nDirect\n');
            fprintf(fid,'%f %f %f\n',(alat*at)\tau(:,I));
        end

        function tocif(prefix,alat,at,tau,atnum,filePath)
            % refer to https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/
            fid=fopen(filePath,'w');
            if fid == -1
                error('QE:OpenFileError', ...
                    'Cannot open this *.cif file: %s', filePath);
            end
            cleaner = onCleanup(@() fclose(fid));
            nat=numel(atnum);
            if isvarname(prefix) && ~strcmp(prefix,'pwscf')
                fprintf(fid,'data_%s\n',prefix);
            else
                [atnum_u,~,ic]=unique(atnum);
                atm=QE.get_atomic_sym(atnum_u);
                nsp=numel(atnum_u);
                na=zeros(nsp,1);
                for i=1:nat
                    na(ic(i))=na(ic(i))+1;
                end
                fprintf(fid,'data_');
                for i=1:nsp
                    fprintf(fid,'%s%d',atm{i},na(i));
                end
                fprintf(fid,'\n');
            end
            dt = datetime('today');
            dateStr = string(dt, 'yyyy-MM-dd');
            audit_lab = {
                '_audit_creation_date'
                '_audit_creation_method'};
            audit_val={dateStr,'''VPHMAT'''};
            at=alat*at;
            tau=at\tau;
            at=at*NIST.Bohr_radius*1e10;
            a=at(:,1);
            b=at(:,2);
            c=at(:,3);
            an=vecnorm(a);
            bn=vecnorm(b);
            cn=vecnorm(c);
            alpha=acosd(dot(b,c)/(bn*cn));
            beta=acosd(dot(a,c)/(an*cn));
            gamma=acosd(dot(a,b)/(an*bn));
            cell_volume=abs(det(at));
            cell_lab={'_cell_length_a','_cell_length_b','_cell_length_c', ...
                '_cell_angle_alpha','_cell_angle_beta','_cell_angle_gamma',...
                '_cell_volume'};
            cell_val={an,bn,cn,alpha,beta,gamma,cell_volume};
            space_group_lab={'_space_group_name_H-M_alt', ...
                '_space_group_IT_number'};
            space_group_val={'''P 1''','1'};
            atom_site_lab={'_atom_site_label'
                '_atom_site_type_symbol'
                '_atom_site_fract_x'
                '_atom_site_fract_y'
                '_atom_site_fract_z'
                '_atom_site_U_iso_or_equiv'
                '_atom_site_adp_type'
                '_atom_site_occupancy'};
            colwidth=35;
            for i = 1:length(audit_lab)
                fprintf(fid,'%-*s %s\n', colwidth, audit_lab{i}, audit_val{i});
            end
            for i = 1:length(cell_lab)
                fprintf(fid,'%-*s %.8f\n', colwidth, cell_lab{i}, cell_val{i});
            end
            for i = 1:length(space_group_lab)
                fprintf(fid,'%-*s %s\n', colwidth, space_group_lab{i}, space_group_val{i});
            end
            fprintf(fid,'loop_\n_space_group_symop_operation_xyz\n   ''x, y, z''\n');
            fprintf(fid,'loop_\n');
            for i=1:length(atom_site_lab)
                fprintf(fid,'%s\n',atom_site_lab{i});
            end
            atm=QE.get_atomic_sym(atnum);
            for i=1:nat
                fprintf(fid,'%-8s %-5s %.8f   %.8f   %.8f   0.00000  Uiso   1.00\n', ...
                    sprintf('%s%d',atm{i},i),atm{i},tau(1,i),tau(2,i),tau(3,i));
            end
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

        function qe = fromCIF(filePath)
            % Initialize a QE object from *.cif
            % please read https://www.iucr.org/resources/cif
            qe=QE;
            fid=fopen(filePath);
            if fid == -1
                error('QE:OpenFileError', ...
                    'Cannot open this *.cif file: %s', filePath);
            end
            cleaner = onCleanup(@() fclose(fid));
            linet = strtrim(fgetl(fid));
            atom_site_type=containers.Map({'_atom_site_type_symbol','_atom_site_label', ...
                '_atom_site_symmetry_multiplicity','_atom_site_fract_x', ...
                '_atom_site_fract_y','_atom_site_fract_z',...
                '_atom_site_occupancy','_atom_site_U_iso_or_equiv', ...
                '_atom_site_adp_type'}, ...
                {'%s','%s','%d','%f','%f','%f','%d','%f','%s'});
            while ~feof(fid)
                lines=strsplit(linet);
                switch lines{1}
                    case '_cell_length_a'
                        a=str2double(lines{2});
                    case '_cell_length_b'
                        b=str2double(lines{2});
                    case '_cell_length_c'
                        c=str2double(lines{2});
                    case '_cell_angle_alpha'
                        alpha=str2double(lines{2});
                    case '_cell_angle_beta'
                        beta=str2double(lines{2});
                    case '_cell_angle_gamma'
                        gamma=str2double(lines{2});
                    case 'loop_'
                        loopn=0;
                        keys{10,1}='tmp';
                        for i=1:1e2
                            tmp=fgetl(fid);
                            if feof(fid)
                                break
                            end
                            tmp=strtrim(tmp);
                            if startsWith(tmp,'_')
                                loopn=loopn+1;
                                keys{loopn}=tmp;
                                ftmp=ftell(fid);
                            else
                                fseek(fid,ftmp,-1);
                                keys1=keys{1};
                                if startsWith(keys1,'_atom_site')
                                    atom_site_idx=containers.Map(keys(1:loopn),1:loopn);
                                    vals=values(atom_site_type,keys(1:loopn));
                                    ftmp=ftell(fid);
                                    for j=1:1e4
                                            tmp=fgetl(fid);
                                            if feof(fid)
                                                j=j+1;
                                                break
                                            elseif startsWith(strtrim(tmp),'_') || startsWith(strtrim(tmp),'loop_')
                                                break
                                            end
                                    end
                                    fseek(fid,ftmp,-1);
                                    atom_site=textscan(fid,strjoin(vals,' '),j-1);
                                    type_symbol=atom_site{:,atom_site_idx('_atom_site_type_symbol')};
                                    fract_x=atom_site{:,atom_site_idx('_atom_site_fract_x')};
                                    fract_y=atom_site{:,atom_site_idx('_atom_site_fract_y')};
                                    fract_z=atom_site{:,atom_site_idx('_atom_site_fract_z')};                                    
                                    tau_cryst_ori=[fract_x fract_y fract_z].';
                                elseif startsWith(keys1,'_symmetry_equiv')
                                    if strcmp(keys1,'_symmetry_equiv_pos_site_id') ...
                                            && strcmp(keys{2},'_symmetry_equiv_pos_as_xyz')
                                        pos_site=textscan(fid,'%d %s %s %s', ...
                                            'Delimiter',' ,''', ...
                                            'MultipleDelimsAsOne',true);
                                        xop=pos_site{2};
                                        yop=pos_site{3};
                                        zop=pos_site{4};
                                    elseif strcmp(keys1,'_symmetry_equiv_pos_as_xyz') ...
                                            && isempty(keys{2})
                                        ftmp=ftell(fid);                                        
                                        for j=1:1e4
                                            tmp=fgetl(fid);
                                            if feof(fid) || startsWith(strtrim(tmp),'_') || startsWith(strtrim(tmp),'loop_')
                                                break
                                            end
                                        end
                                        fseek(fid,ftmp,-1);
                                        pos_site=textscan(fid,'%s %s %s',j-1, ...
                                            'Delimiter',' ,''', ...
                                            'MultipleDelimsAsOne',true);
                                        xop=pos_site{1};
                                        yop=pos_site{2};
                                        zop=pos_site{3};
                                    end
                                end
                                break
                            end
                        end
                end
                if feof(fid)
                    break
                end
                linet=strtrim(fgetl(fid));
            end
            nop=size(xop,1);
            tau_cryst_tmp=zeros(3,size(tau_cryst_ori,2),nop);
            for iop=1:nop
                fx=str2func(['@(x,y,z)' xop{iop}]);
                fy=str2func(['@(x,y,z)' yop{iop}]);
                fz=str2func(['@(x,y,z)' zop{iop}]);
                tau_cryst_tmp(:,:,iop)=[fx(tau_cryst_ori(1,:),tau_cryst_ori(2,:),tau_cryst_ori(3,:));
                    fy(tau_cryst_ori(1,:),tau_cryst_ori(2,:),tau_cryst_ori(3,:));
                    fz(tau_cryst_ori(1,:),tau_cryst_ori(2,:),tau_cryst_ori(3,:))];
            end
            tau_cryst_tmp=mod(permute(tau_cryst_tmp,[1 3 2]),1);
            nat=0;
            atm=regexp(type_symbol, '^[A-Za-z]+', 'match');
            atm=[atm{:}].';
            atnumtmp=QE.get_atomic_number(atm);
            tau_cryst=zeros(3,1e4);
            atnum=zeros(1e4,1);
            for i=1:size(tau_cryst_tmp,3)
                tau_cryst_t=tau_cryst_tmp(:,:,i);
                tau_cryst_tu=QE.uniquetolcir(tau_cryst_t',1e-6,1);
                natp=nat+size(tau_cryst_tu,1);
                tau_cryst(:,nat+1:natp)=tau_cryst_tu.';
                atnum(nat+1:natp)=atnumtmp(i);
                nat=natp;
            end
            qe.nat=nat;            
            tau_cryst=tau_cryst(:,1:nat);
            atnum=atnum(1:nat);
            [atnumu,~,qe.ityp]=unique(atnum,'stable');
            qe.atm=QE.get_atomic_sym(atnumu);
            qe.atnum=atnum;
            qe.nsp=numel(qe.atm);
            na=zeros(qe.nsp,1);
            for ia=1:nat
                na(qe.ityp(ia))=na(qe.ityp(ia))+1;
            end
            qe.na=na;
            at=[a b*cosd(gamma) c*cosd(beta);
                0 b*sind(gamma) c*(cosd(alpha)-cosd(beta)*cosd(gamma))/sind(gamma)
                0 0 c*sqrt(1-cosd(alpha)^2-cosd(beta)^2-cosd(gamma)^2+2*cosd(alpha)...
                *cosd(beta)*cosd(gamma))/sind(gamma)];
            qe.at=at/(NIST.Bohr_radius*1e10);
            qe.alat=1;
            qe.tpiba = 2 * pi / qe.alat;
            qe.bg = inv(qe.at.'); % reciprocal lattice
            qe.tau=qe.alat*qe.at*tau_cryst;
        end
    end

    methods
        function this = QE()
            % Empty constructor
        end

        function exportPOSCAR(this,filePath)
            QE.toPOSCAR(this.alat,this.at,this.tau,this.atnum,filePath)
        end

        function exportCIF(this,filePath)
            QE.tocif(this.prefix,this.alat,this.at,this.tau,this.atnum,filePath)
        end

        function pwinput(this,kpr)
        end
        function cp2kinput(this)
            fid=fopen('cp2k.inp','w+');
            if fid == -1
                error('QE:OpenFileError', ...
                    'Cannot open this cp2k.inp file: %s', 'cp2k.inp');
            end
            cleaner = onCleanup(@() fclose(fid));
            sglobal=['&GLOBAL\n' ...
                '    PROJECT %s\n' ...
                '    RUN_TYPE GEO_OPT\n' ...
                '    PRINT_LEVEL LOW\n' ...
                '&END GLOBAL\n\n'];
            fprintf(fid,sglobal,this.prefix);
            smotion=['&MOTION\n' ...
                ' &GEO_OPT\n' ...
                '    MAX_ITER 400\n' ...
                '    MAX_FORCE 6.0E-4\n' ...
                '    OPTIMIZER CG\n' ...
                ' &END GEO_OPT\n' ...
                '&END MOTION\n\n'];
            fprintf(fid,smotion);
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

        function atsym = get_atomic_sym(atnum)
            % get_atomic_sym Returns element symbol(s) given atomic number(s)
            % Input: atnum - scalar or array of atomic numbers

            symcell={'H','He','Li','Be','B','C','N','O','F','Ne', ...
                'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca', ...
                'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn', ...
                'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr', ...
                'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', ...
                'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd', ...
                'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', ...
                'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg', ...
                'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', ...
                'Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', ...
                'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', ...
                'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'}.';
            atsym=symcell(atnum);
        end

        function [C,IA,IC]=uniquetolcir(A,tol,period)
            % A is col vector or matrix 'ByRows'
            % mod(a-b,period)<tol | mod(a-b,period)>period-tol: a,b are the same
            % period is a scalar or row vector (1,size(A,2))
            % C = A(IA,:)
            % A~C(IC)
            n=size(A,1);
            if n==1
                C=A;
                IA=1;
                IC=1;
                return
            end
            idx=nchoosek(1:n,2);
            dA=A(idx(:,1),:)-A(idx(:,2),:);
            lA=all(mod(dA,period)<tol | mod(dA,period)>period-tol,2);
            idx=idx(lA,:);
            % hsame=eye(n,n,'logical');
            % hsame(sub2ind([n n],idx(:,1),idx(:,2)))=true;
            % hsame=hsame | hsame.';
            tmp=1:n;
            for i=1:n
                tmp(idx(i,2))=tmp(idx(i,1));
            end
            [~,IA,IC]=unique(tmp);
            C=A(IA,:);
        end
    end
end