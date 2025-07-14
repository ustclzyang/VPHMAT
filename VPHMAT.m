classdef VPHMAT < handle
    %VPHMAT   Visualize Quantum ESPRESSO crystal structures and phonon modes.
    %   VPHMAT(QE) creates a VPHMAT object for visualizing crystal structures parsed
    %   from a Quantum ESPRESSO (QE) object. The VPHMAT class provides various
    %   visualization styles and customization options for atoms, bonds, and cells.
    %
    %   V = VPHMAT(QE) constructs a visualizer for the QE structure object QE.
    %
    %   V.plot() displays the atomic structure, including atoms, bonds, and the
    %   unit cell, in a 3D figure. The visualization can be customized via the
    %   'style', 'view', and 'supercell' properties.
    %
    %   VPHMAT also supports interactive visualization of phonon modes, allowing
    %   users to animate and explore vibrational eigenmodes obtained from phonon
    %   calculations. Users can select a phonon branch and q-point interactively,
    %   and visualize atomic displacements in real time.
    %
    %   Supported visualization styles include:
    %     'Ball-and-stick'    Atoms as spheres, bonds as cylinders (default)
    %     'Space-filling'     Atoms shown at their van der Waals radii
    %     'Wireframe'         Bonds as lines, atoms represented minimally
    %     'Stick'             Bonds as cylinders, small atoms
    %     'Polyhedral'        (If implemented) Polyhedral representation
    %
    %   The supercell property allows expansion of the primitive cell in the
    %   a, b, and c directions, e.g. V.supercell = [2 2 2] for a 2x2x2 supercell.
    %
    %   VPHMAT automatically assigns CPK or Jmol colors to atoms based on their
    %   atomic numbers, and determines bond connectivity using empirical or
    %   tabulated bond lengths. Bonds can be customized and edited interactively.
    %
    %   For interactive phonon visualization:supercell
    %      V.animate() enables branch and q-point selection and animates atomic
    %      displacements according to the chosen phonon eigenmode.
    %
    %   Example:
    %      qe = QE.fromPOSCAR('POSCAR'); % Parse POSCAR to QE object
    %      v = VPHMAT(qe);                % Create VPHMAT visualizer
    %      v.style = 'Space-filling';    % Set visualization style
    %      v.supercell = [2 2 1];        % Expand cell
    %      v.plot();                     % Visualize structure
    %
    %      % Interactive phonon mode visualization (requires phonon data)
    %      v.matdyn = MATDYN('matdyn.modes'); % Load phonon dispersions
    %      v.animate();                      % Start interactive phonon animation
    %
    %   For more information about Quantum ESPRESSO, visit:
    %     https://www.quantum-espresso.org/
    %   For the original VESTA visualization software, see:
    %     https://jp-minerals.org/vesta/en/
    %   For the original tools-phonon-dispersion, see:
    %     https://github.com/materialscloud-org/tools-phonon-dispersion
    %
    %   See also: QE, NIST, MATDYN
    %
    %   Copyright (c) 2025 Lizhou Yang
    %   Licensed under the MIT License.

    properties
        qe
        matdyn
        general
        data
        atoms
        bonds
        polyhedra
        isosurfaces
        ax
        phonon
        light
        hatom     % ploted atoms' handle
        hbond     % ploted bonds' handle
    end

    properties (SetObservable)
        style
        view
        supercell
        boundary
    end

    methods
        function this = VPHMAT(qe)
            this.qe = qe;
            this.data.prefix=qe.prefix;
            this.data.atnum=qe.atnum;
            this.data.at=qe.at*qe.alat*NIST.Bohr_radius*1e10; % in \AA unit
            this.data.bg=qe.bg;
            this.data.nat=qe.nat;
            this.data.tau=qe.tau*NIST.Bohr_radius*1e10; % in \AA unit
            % this.boundary='none'; % 'none': Do not search atoms beyond the boundary
            % this.boundary='A1'; % 'A1': Search additional atoms within 1 \AA near each face, see initA1
            this.boundary='B1'; % 'B1': Search additional atoms that can form bonds with atoms in unitcell and its atomic number is smaller
            
            this.bonds = this.initBonds();
            this.atoms = this.initAtoms();
            switch this.boundary
                case 'A1'
                    this.initA1();
            end
            this.initSupercell();
            this.general=this.initGeneral();
            this.style = 'Ball-and-stick'; % 'Space-filling' 'Polyhedral' 'Wireframe' 'Stick'
            this.initStyle();
            addlistener(this, 'style', 'PostSet', @(src, evt)this.onStyleChanged());
            addlistener(this, 'view', 'PostSet', @(src, evt)this.onViewChanged());
            addlistener(this, 'supercell', 'PostSet', @(src, evt)this.onSupercellChanged());
            addlistener(this, 'boundary', 'PostSet', @(src, evt)this.onBoundaryChanged());
            this.phonon.amplitude=0.5;
            this.phonon.speed=1;
        end

        function onStyleChanged(this)
            this.initStyle();
        end

        function onViewChanged(this)
            if isempty(this.ax) || ~isvalid(this.ax)
                return;
            end
            if ischar(this.view) || isstring(this.view)
                switch char(this.view)
                    case 'a'
                        view(this.ax, this.data.at(:,1));
                    case 'a*'
                        view(this.ax, -this.data.at(:,1));
                    case 'b'
                        view(this.ax, this.data.at(:,2));
                    case 'b*'
                        view(this.ax, -this.data.at(:,2));
                    case 'c'
                        view(this.ax, this.data.at(:,3));
                    case 'c*'
                        view(this.ax, -this.data.at(:,3));
                    case '3'
                        view(this.ax, 3)
                end
            elseif isnumeric(this.view) && numel(this.view) == 2
                view(this.ax, this.view(1), this.view(2));
            end
        end

        function onSupercellChanged(this)
            n1=this.supercell(1);
            n2=this.supercell(2);
            n3=this.supercell(3);
            [cx,cy,cz]=meshgrid(2:n1+1,2:n2+1,2:n3+1);
            x=cx(:);
            y=cy(:);
            z=cz(:);

            % change atoms ---------------------------------------
            atomshash=false(this.data.nat,n1+2,n2+2,n3+2);
            celltmp=this.atoms.cell;
            for i=1:size(x,1)
                xp=x(i)+celltmp(2,:);
                yp=y(i)+celltmp(3,:);
                zp=z(i)+celltmp(4,:);
                ind = sub2ind(size(atomshash), celltmp(1,:), xp, yp, zp);
                atomshash(ind) = true;
            end
            [idx1, idx2, idx3, idx4] = ind2sub(size(atomshash), find(atomshash));
            % this.atoms.hash=atomshash;
            idx2=idx2-2; idx3=idx3-2; idx4=idx4-2;
            this.atoms.supercell=[idx1 idx2 idx3 idx4].';
            this.atoms.nat=size(idx1,1);
            this.atoms.tau=this.data.tau(:,idx1)+this.data.at*this.atoms.supercell(2:4,:);
            this.atoms.tau_bal=this.atoms.tau;
            this.atoms.atnum=this.data.atnum(idx1);
            this.atoms.color=VPHMAT.CPK_coloring(this.atoms.atnum,'J');
            this.atoms.radii=this.data.radii(idx1);

            % change bonds ----------------------------------------
            this.changeBondsFromAtoms();
        end

        function onBoundaryChanged(this)
            switch this.boundary
                case 'A1'
                    this.initA1();
                % case 'B1'
                %     this.initB1();
                case 'none'
                    this.initnone();
            end
            this.changeBondsFromAtoms();
        end

        function initStyle(this)
            switch this.style
                case 'Ball-and-stick'
                    this.bonds.style='bicolorCylinder';
                    this.atoms.radiiScale=0.4;
                case 'Space-filling'
                    this.bonds.style='none';
                    this.atoms.radiiScale=1;
                case 'Polyhedral'
                    % this.bonds.style='bicolorCylinder';
                case 'Wireframe'
                    this.bonds.style='gradientLine';
                    this.atoms.radiiScale=NaN;
                case 'Stick'
                    this.bonds.style='bicolorCylinder';
                    this.atoms.radiiScale=-Inf;
            end
        end

        function bonds=initBonds(this)
            % Reference: https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/chem.200800987

            % in the unit cell --------------------------------
            bondind=nchoosek(1:this.data.nat,2);
            bondindn=size(bondind,1);
            tmp=zeros(bondindn,7);
            tmp(:,1:2)=bondind;                 % all posible atom indices that form bonds
            tmp(:,3:4)=this.data.atnum(bondind); % two chemical componds (atomic number) that form bonds
            tmp(:,5)=vecnorm(this.data.tau(:,tmp(:,1))-this.data.tau(:,tmp(:,2))).'; % bond length, bohr
            % Reference: https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
            % Covalent (single bond) in pm units
            bonds.radiiType='wikipedia';
            coval=[32 46 133 102 85 75 71 63 64 67 155 139 126 116 ...
                111 103 99 96 196 171 148 136 134 122 119 116 111 110 112 118 124 ...
                121 121 116 114 117 210 185 163 154 147 138 128 125 125 120 128 136 ...
                142 140 140 136 133 131 232 196 180 163 176 174 173 172 168 169 168 ...
                167 166 165 164 170 162 152 146 137 131 129 122 123 124 133 144 144 ...
                151 145 147 142 NaN 201 186 175 169 170 171 172 166 166 NaN NaN NaN ...
                NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
            coval(isnan(coval))=120;
            bonds.radii=coval/100; % in \AA units

            % Determine critical bond length
            filename = fullfile('config', 'bondlengthdic.mat');
            if exist(filename, 'file')
                S = load(filename, 'bondlengthdic');
                bondlengthdic = S.bondlengthdic;
                if anynan(bondlengthdic)
                    bondlengthdictmp=bonds.radii+bonds.radii.';
                    bondlengthdic(isnan(bondlengthdic))=bondlengthdictmp(isnan(bondlengthdic));
                    save(filename, 'bondlengthdic');
                end
            else
                if ~exist('config', 'dir')
                    mkdir('config');
                end
                bondlengthdic=bonds.radii+bonds.radii.';
                save(filename, 'bondlengthdic');
            end
            for i = 1:bondindn
                tmp(i,6) = bondlengthdic(tmp(i,3),tmp(i,4)); % corresponding critical bond length
            end
            % tmp(:,6)=bonds.radii(tmp(:,3))+bonds.radii(tmp(:,4)); % bond length=radii(Z1) + radii(Z2) + delta

            tmp(:,7)=tmp(:,6)>tmp(:,5); % determin whether to plot this bond
            bonds.unitcell=tmp;

            % in the nearby cells --------------------------------
            [cx,cy,cz]=meshgrid(-1:1);
            bondind=combinations(1:this.data.nat,1:this.data.nat);
            bondind=[bondind.Var1 bondind.Var2];
            bondindn=size(bondind,1);
            nearcellnum=size(cx(:),1); % including unit cell itself
            tmp=zeros(bondindn,5+2*nearcellnum);
            tmp(:,1:2)=bondind;                  % col1 represents atoms in unicell, col2 represents atoms in nearby cells
            tmp(:,3:4)=this.data.atnum(bondind); % two chemical componds (atomic number) that form bonds
            for i = 1:bondindn
                tmp(i,5) = bondlengthdic(tmp(i,3),tmp(i,4)); % corresponding critical bond length
            end
            centerind=14;
            for ni=1:27
                if ni~=centerind % ni==centerind is unit cell
                    vnear=this.data.at*[cx(ni); cy(ni); cz(ni)];
                    tmp(:,2*ni+4)=vecnorm(this.data.tau(:,tmp(:,1))-this.data.tau(:,tmp(:,2))-vnear).';
                    tmp(:,2*ni+5)=tmp(:,5)>tmp(:,2*ni+4); % determin whether to plot this bond
                end
            end
            bonds.nearcell=tmp; % nothing to do with plot style, but depends on bond length

            formedindex=bonds.unitcell(bonds.unitcell(:,7)==1,1:2);
            formed=zeros(size(formedindex,1),10); % ia, ib, ia's cord, ib's cord, ia's meshgrid, ib's meshgrid
            formed(:,1:2)=formedindex;
            formed(:,3:5)=this.data.tau(:,formed(:,1)).';
            formed(:,6:8)=this.data.tau(:,formed(:,2)).';
            formed(:,9:10)=centerind;

            switch this.boundary
                case 'none'
                    bonds.formed=formed;
                case 'A1' % A1 style is special, has bonds that bonds.nearcell not includes, so I handle it after initAtoms
                    bonds.formed=formed;
                case 'B1'
                    nearcell=bonds.nearcell;
                    for ni=1:nearcellnum
                        nearcell(:,2*ni+5)= logical(nearcell(:,2*ni+5)) & this.data.atnum(nearcell(:,1))>=this.data.atnum(nearcell(:,2)); % B1 style, formed bands between centercell and outercell
                    end
                    formed_near=zeros(sum(nearcell(:,7:2:end),"all"),10);
                    i=1;
                    for ni=1:nearcellnum
                        if ni~=centerind % ni==centerind is unit cell
                            tmp=nearcell(logical(nearcell(:,2*ni+5)),1:2); % col1 represents atoms in unicell, col2 represents atoms in nearby cells
                            ii=i+size(tmp,1)-1;
                            vnear=this.data.at*[cx(ni); cy(ni); cz(ni)];
                            formed_near(i:ii,1:2)=tmp;
                            formed_near(i:ii,3:5)=this.data.tau(:,tmp(:,1)).';
                            formed_near(i:ii,6:8)=(this.data.tau(:,tmp(:,2))+vnear).';
                            formed_near(i:ii,9)=centerind;
                            formed_near(i:ii,10)=ni;
                            i=ii+1;
                        end
                    end
                    bonds.formed=[formed; formed_near];
            end

            if size(bonds.formed,1) > 2000
                bonds.style='gradientLine';
            elseif isfield(this.bonds,'style')
                bonds.style=this.bonds.style;
            else
                bonds.style='bicolorCylinder';
            end
            bonds.radius=0.10; % only avail for cylinder
            bonds.width=2.0; % only avail for line
        end

        function initA1(this)
            [cx,cy,cz]=meshgrid(-1:1);
            cx=cx(:);
            cy=cy(:);
            cz=cz(:);
            tautmp=repmat(this.qe.at\this.qe.tau,[1 27]); % direct
            iatmp=repmat(1:this.qe.nat,[1 27]);
            Rtmp=repelem([cx cy cz].',1,this.qe.nat);
            tautmp=tautmp+Rtmp;
            d=1e-10/NIST.Bohr_radius; % Bohr
            originmovd=VPHMAT.expandat(this.qe.alat,this.qe.at,this.qe.bg,d); % direct
            ind=all(tautmp>=originmovd & tautmp<=(1-originmovd));  % ind(qe.nat*14-qe.nat+1:qe.nat*14) is true
            this.atoms.nat=sum(ind);
            this.atoms.tau=1e10*NIST.Bohr_radius*this.qe.alat*this.qe.at*tautmp(:,ind); % \AA
            this.atoms.tau_bal=this.atoms.tau;
            this.atoms.cell=[iatmp(ind);Rtmp(:,ind)];
            this.atoms.atnum=this.qe.atnum(this.atoms.cell(1,:));
            this.atoms.color=VPHMAT.CPK_coloring(this.atoms.atnum,'J');
            this.atoms.radii=this.data.radii(this.atoms.cell(1,:));
            this.bonds.boudary='A1';
        end

        function initnone(this)
            this.atoms.nat=this.data.nat;
            this.atoms.tau=this.data.tau;
            this.atoms.tau_bal=this.atoms.tau;
            this.atoms.cell=[1:this.atoms.nat;zeros(3,this.atoms.nat)];
            this.atoms.atnum=this.data.atnum;
            this.atoms.color=VPHMAT.CPK_coloring(this.atoms.atnum,'J');
            this.atoms.radii=this.data.radii(this.atoms.cell(1,:));
        end

        function changeBondsFromAtoms(this)
            bondind=nchoosek(1:this.atoms.nat,2);
            bondindn=size(bondind,1);
            tmp=zeros(bondindn,7);
            tmp(:,1:2)=bondind;                 % all posible atom indices that form bonds
            tmp(:,3:4)=this.atoms.atnum(bondind); % two chemical componds (atomic number) that form bonds
            tmp(:,5)=vecnorm(this.atoms.tau(:,tmp(:,1))-this.atoms.tau(:,tmp(:,2))).'; % bond length, bohr

            % Determine critical bond length
            filename = fullfile('config', 'bondlengthdic.mat');
            if exist(filename, 'file')
                S = load(filename, 'bondlengthdic');
                bondlengthdic = S.bondlengthdic;
            else
                error('VPHMAT:OpenFileError', ...
                    '%s does not exist', filename);
            end
            for i = 1:bondindn
                tmp(i,6) = bondlengthdic(tmp(i,3),tmp(i,4)); % corresponding critical bond length
            end
            tmp(:,7)=tmp(:,6)>tmp(:,5);

            formedindex=tmp(tmp(:,7)==1,1:2);
            formed=zeros(size(formedindex,1),10); % ia, ib, ia's cord, ib's cord, ia's meshgrid, ib's meshgrid
            formed(:,1:2)=formedindex;
            formed(:,3:5)=this.atoms.tau(:,formed(:,1)).';
            formed(:,6:8)=this.atoms.tau(:,formed(:,2)).';
            formed(:,9:10)=0;
            this.bonds.formed=formed;
        end

        function atoms=initAtoms(this)
            % Reference: https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
            % empirical in pm units
            empirical=[25 120 145 105 85 70 65 60 50 160 180 150 125 110 ...
                100 100 100 71 220 180 160 140 135 140 140 140 135 135 135 135 ...
                130 125 115 115 115 NaN 235 200 180 155 145 145 135 130 135 140 ...
                160 155 155 145 145 140 140 NaN 260 215 195 185 185 185 185 185 ...
                185 180 175 175 175 175 175 175 175 155 145 135 135 130 135 135 ...
                135 150 190 180 160 190 NaN NaN NaN 215 195 180 180 175 175 175 ...
                175 176 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
            empirical(isnan(empirical))=150;

            nearcellnum=(size(this.bonds.nearcell,2)-5)/2; % including center cell
            atomsbondsnum=zeros(this.data.nat,nearcellnum); % center cell col index is this.bonds.formed(1,end-1)
            atomsneartau=nan(this.data.nat,nearcellnum,3);

            if size(this.bonds.formed,1)
                centerind=this.bonds.formed(1,end-1);
                for i=1:size(this.bonds.formed,1)
                    atomsbondsnum(this.bonds.formed(i,1),this.bonds.formed(i,end-1))=atomsbondsnum(this.bonds.formed(i,1),this.bonds.formed(i,end-1))+1;
                    atomsbondsnum(this.bonds.formed(i,2),this.bonds.formed(i,end))=atomsbondsnum(this.bonds.formed(i,2),this.bonds.formed(i,end))+1;
                    atomsneartau(this.bonds.formed(i,2),this.bonds.formed(i,end),:)=this.bonds.formed(i,6:8);
                end
                nat_cent=this.data.nat;
                atomsbondsnump=atomsbondsnum;
                atomsbondsnump(:,centerind)=0;
                nat_near_ind=atomsbondsnump~=0;
                [row,col] = find(nat_near_ind); % row: atom index in unit cell, col: cell index
                [cx,cy,cz]=meshgrid(-1:1);
                cx=cx(:);
                cy=cy(:);
                cz=cz(:);
                nat_near=sum(nat_near_ind,"all");
                atoms.nat=nat_cent+nat_near;
                atoms.tau=zeros(3,atoms.nat);
                atoms.cell=zeros(4,atoms.nat); % each col is atom index in unit cell; cell vector in data.at, dependends on boundary style
                atoms.cell(1,1:this.data.nat)=1:this.data.nat;
                atoms.tau(:,1:nat_cent)=this.data.tau;
                for i=1:nat_near
                    atoms.tau(:,i+nat_cent)=atomsneartau(row(i),col(i),:);
                    atoms.cell(:,i+nat_cent)=[row(i);cx(col(i));cy(col(i));cz(col(i))];
                end
            else
                atoms.nat=this.data.nat;
                atoms.cell=zeros(4,atoms.nat); % each col is atom index in unit cell; cell vector in data.at, dependends on boundary style
                atoms.cell(1,1:this.data.nat)=1:this.data.nat;
                atoms.tau=this.data.tau;
            end

            atoms.tau_bal=atoms.tau;
            atoms.atnum=this.data.atnum(atoms.cell(1,:));
            atoms.radiiType='wikipedia';
            atoms.radii=empirical(atoms.atnum)/100; % in \AA units
            this.data.radii=empirical(atoms.atnum)/100;
            atoms.radiiScale=0.4;
            if atoms.nat > 1000
                atoms.radiiScale=NaN;
            end
            atoms.color=VPHMAT.CPK_coloring(atoms.atnum,'J');
        end

        function initSupercell(this)
            this.supercell=[1 1 1];
            this.changeBondsFromAtoms();
            this.atoms.supercell=this.atoms.cell;
        end

        function general=initGeneral(this)
            general.unitCell.line='allUnitCells';
            general.unitCell.lineStyle='-'; % '-' ':' '--'
            general.unitCell.lineWidth=0.5;
            general.unitCell.color='k';
            general.axes.showCompass=true;
            general.axes.showAxisLabels=true;
        end

        function plot(this)
            % plot: Visualize structure including atoms, bonds, and cell
            % figure
            hold on
            plotBonds(this);
            plotAtoms(this);
            plotCell(this);
            axis equal
            axis off
            if isempty(this.light) || ~isvalid(this.light)
                this.light = light;
                this.light.Color = [1 1 1];
                this.light.Position = [1 0 1];
                lighting GOURAUD
            end
            this.ax=gca;
            this.view=[70 30];
        end

        function plotBonds(this)
            switch this.bonds.style
                case 'dottedLine'
                    this.hbond=gobjects(size(this.bonds.formed,1),1);
                    for i=1:size((this.bonds.formed),1)
                        this.hbond(i)=plot3(this.bonds.formed(i,[3 6]),...
                            this.bonds.formed(i,[4 7]),this.bonds.formed(i,[5 8]),...
                            'LineWidth',this.bonds.width,'LineStyle',':','Color',[0.5 0.5 0.5]);
                    end
                case 'bicolorCylinder'
                    this.hbond=gobjects(size(this.bonds.formed,1),2);
                    for i=1:size((this.bonds.formed),1)
                        v1=this.bonds.formed(i,3:5).';
                        v2=this.bonds.formed(i,6:8).';
                        ia=this.bonds.formed(i,1);
                        ib=this.bonds.formed(i,2);
                        c1=this.atoms.color(ia,:);
                        c2=this.atoms.color(ib,:);
                        r1=this.atoms.radii(ia)*this.atoms.radiiScale;
                        r2=this.atoms.radii(ib)*this.atoms.radiiScale;
                        this.hbond(i,:)=VPHMAT.plotBicolorCylinder(v1,v2,c1,c2,this.bonds.radius,r1,r2); % v1 is col vector, c1 is the row color
                    end
                case 'gradientLine'
                    this.hbond=gobjects(size(this.bonds.formed,1),1);
                    for i=1:size((this.bonds.formed),1)
                        v1=this.bonds.formed(i,3:5).';
                        v2=this.bonds.formed(i,6:8).';
                        ia=this.bonds.formed(i,1);
                        ib=this.bonds.formed(i,2);
                        c1=this.atoms.color(ia,:);
                        c2=this.atoms.color(ib,:);
                        this.hbond(i)=VPHMAT.plotGradientLine(v1,v2,c1,c2,2);
                    end
            end
        end

        function plotBonds2(this)
            % change bond positions ploted by poltBonds
            switch this.bonds.style
                case 'dottedLine'
                    for i=1:size((this.bonds.formed),1)
                        set(this.hbond(i), 'XData', this.bonds.formed(i,[3 6]), ...
                            'YData', this.bonds.formed(i,[4 7]), ...
                            'ZData', this.bonds.formed(i,[5 8]));
                    end
                case 'bicolorCylinder'
                    for i=1:size((this.bonds.formed),1)
                        v1=this.bonds.formed(i,3:5).';
                        v2=this.bonds.formed(i,6:8).';
                        ia=this.bonds.formed(i,1);
                        ib=this.bonds.formed(i,2);
                        c1=this.atoms.color(ia,:);
                        c2=this.atoms.color(ib,:);
                        r1=this.atoms.radii(ia)*this.atoms.radiiScale;
                        r2=this.atoms.radii(ib)*this.atoms.radiiScale;
                        VPHMAT.plotBicolorCylinder2(v1,v2,c1,c2,this.bonds.radius,r1,r2,this.hbond(i,:)); % v1 is col vector, c1 is the row color
                    end
                case 'gradientLine'
                    for i=1:size((this.bonds.formed),1)
                        v1=this.bonds.formed(i,3:5).';
                        v2=this.bonds.formed(i,6:8).';
                        ia=this.bonds.formed(i,1);
                        ib=this.bonds.formed(i,2);
                        c1=this.atoms.color(ia,:);
                        c2=this.atoms.color(ib,:);
                        VPHMAT.plotGradientLine2(v1,v2,c1,c2,2,this.hbond(i));
                    end
            end
        end

        function plotAtoms(this)
            [X,Y,Z] = sphere;
            this.hatom = gobjects(this.atoms.nat,1);
            if ~isnan(this.atoms.radiiScale)
                if this.atoms.radiiScale==-Inf
                    for ia=1:this.atoms.nat
                        r=this.bonds.radius;
                        x=this.atoms.tau(1,ia);
                        y=this.atoms.tau(2,ia);
                        z=this.atoms.tau(3,ia);
                        this.hatom(ia)=surf(X*r+x,Y*r+y,Z*r+z, ...
                            'EdgeColor', 'none', 'FaceColor', this.atoms.color(ia,:));
                    end
                else
                    for ia=1:this.atoms.nat
                        r=this.atoms.radii(ia)*this.atoms.radiiScale;
                        x=this.atoms.tau(1,ia);
                        y=this.atoms.tau(2,ia);
                        z=this.atoms.tau(3,ia);
                        this.hatom(ia)=surf(X*r+x,Y*r+y,Z*r+z, ...
                            'EdgeColor', 'none', 'FaceColor', this.atoms.color(ia,:));
                    end
                end
            end
        end

        function plotAtoms2(this)
            % change atom positions ploted by poltAtoms
            [X,Y,Z] = sphere;
            if ~isnan(this.atoms.radiiScale)
                if this.atoms.radiiScale==-Inf
                    for ia=1:this.atoms.nat
                        r=this.bonds.radius;
                        x=this.atoms.tau(1,ia);
                        y=this.atoms.tau(2,ia);
                        z=this.atoms.tau(3,ia);
                        set(this.hatom(ia), 'XData', X*r+x, 'YData', Y*r+y, 'ZData', Z*r+z);
                    end
                else
                    for ia=1:this.atoms.nat
                        r=this.atoms.radii(ia)*this.atoms.radiiScale;
                        x=this.atoms.tau(1,ia);
                        y=this.atoms.tau(2,ia);
                        z=this.atoms.tau(3,ia);
                        set(this.hatom(ia), 'XData', X*r+x, 'YData', Y*r+y, 'ZData', Z*r+z);
                    end
                end
            end
        end

        function plotCell(this)
            switch this.general.unitCell.line
                case 'allUnitCells'
                    plot3([0 this.data.at(1,1) this.data.at(1,2)+this.data.at(1,1) this.data.at(1,2) 0], ...
                        [0 this.data.at(2,1) this.data.at(2,2)+this.data.at(2,1) this.data.at(2,2) 0], ...
                        [0 this.data.at(3,1) this.data.at(3,2)+this.data.at(3,1) this.data.at(3,2) 0], ...
                        'Color',this.general.unitCell.color,'LineWidth',this.general.unitCell.lineWidth,'LineStyle',this.general.unitCell.lineStyle);
                    plot3([0 this.data.at(1,1) this.data.at(1,2)+this.data.at(1,1) this.data.at(1,2) 0]+this.data.at(1,3), ...
                        [0 this.data.at(2,1) this.data.at(2,2)+this.data.at(2,1) this.data.at(2,2) 0]+this.data.at(2,3), ...
                        [0 this.data.at(3,1) this.data.at(3,2)+this.data.at(3,1) this.data.at(3,2) 0]+this.data.at(3,3), ...
                        'Color',this.general.unitCell.color,'LineWidth',this.general.unitCell.lineWidth,'LineStyle',this.general.unitCell.lineStyle);
                    plot3([0 this.data.at(1,3)], [0 this.data.at(2,3)], [0 this.data.at(3,3)], ...
                        'Color',this.general.unitCell.color,'LineWidth',this.general.unitCell.lineWidth,'LineStyle',this.general.unitCell.lineStyle);
                    plot3(this.data.at(1,1)+[0 this.data.at(1,3)], this.data.at(2,1)+[0 this.data.at(2,3)], this.data.at(3,1)+[0 this.data.at(3,3)], ...
                        'Color',this.general.unitCell.color,'LineWidth',this.general.unitCell.lineWidth,'LineStyle',this.general.unitCell.lineStyle);
                    plot3(this.data.at(1,2)+this.data.at(1,1)+[0 this.data.at(1,3)], this.data.at(2,2)+this.data.at(2,1)+[0 this.data.at(2,3)], this.data.at(3,2)+this.data.at(3,1)+[0 this.data.at(3,3)], ...
                        'Color',this.general.unitCell.color,'LineWidth',this.general.unitCell.lineWidth,'LineStyle',this.general.unitCell.lineStyle);
                    plot3(this.data.at(1,2)+[0 this.data.at(1,3)], this.data.at(2,2)+[0 this.data.at(2,3)], this.data.at(3,2)+[0 this.data.at(3,3)], ...
                        'Color',this.general.unitCell.color,'LineWidth',this.general.unitCell.lineWidth,'LineStyle',this.general.unitCell.lineStyle);
            end
            if this.general.axes.showCompass
                plot3([0 this.data.at(1,1)], [0 this.data.at(2,1)], [0 this.data.at(3,1)], ...
                    'Color','r','LineWidth',this.general.unitCell.lineWidth*2,'LineStyle','-');
                plot3([0 this.data.at(1,2)], [0 this.data.at(2,2)], [0 this.data.at(3,2)], ...
                    'Color','g','LineWidth',this.general.unitCell.lineWidth*2,'LineStyle','-');
                plot3([0 this.data.at(1,3)], [0 this.data.at(2,3)], [0 this.data.at(3,3)], ...
                    'Color','b','LineWidth',this.general.unitCell.lineWidth*2,'LineStyle','-');
            end
            if this.general.axes.showAxisLabels
                text(this.data.at(1,1)/2, this.data.at(2,1)/2, this.data.at(3,1)/2,'a');
                text(this.data.at(1,2)/2, this.data.at(2,2)/2, this.data.at(3,2)/2,'b');
                text(this.data.at(1,3)/2, this.data.at(2,3)/2, this.data.at(3,3)/2,'c');
            end
        end

        function animate(this)
            this.matdyn.plot();
            hlines = this.matdyn.hLines;
            for ib = 1:numel(hlines)
                set(hlines(ib), ...
                    'ButtonDownFcn', @(src, evt) this.onLineClicked(src, evt,ib), ...
                    'PickableParts', 'all', ...
                    'HitTest', 'on');
            end
            fig=figure;
            fig.Position=fig.Position+[400 0 0 0];
            this.plot();
        end

        function onLineClicked(this, src, evt, ib)
            clickPoint = evt.IntersectionPoint; % 3x1 [x y z]
            X = get(src, 'XData');
            Y = get(src, 'YData');
            iq=VPHMAT.findNearestPoint(this.matdyn.ax,X,Y,clickPoint(1),clickPoint(2));
            if ~isempty(this.matdyn.hpoint) && isvalid(this.matdyn.hpoint)
                delete(this.matdyn.hpoint); % remove last set hight light point
            end
            this.matdyn.hpoint=plot(X(iq), Y(iq),"Color","#0072BD",'Marker','.','MarkerSize',10);


            fprintf('\nbranch = %d\nq = %f %f %f\n',ib,this.matdyn.q(:,iq));
            % z=this.matdyn.z(:,ib,iq);
            % z=reshape(z,[3,this.matdyn.nat]);
            % disp(z)
            axes(this.ax)
            fig=gcf;
            xl = xlim(this.ax);
            yl = ylim(this.ax);
            zl = zlim(this.ax);
            for i=0:1000
                if ~isvalid(fig) || ~isequal(gca, this.ax)
                    break;
                end
                t=i;
                if isfield(this.phonon,'save') && ~isempty(this.phonon.save)
                    this.phonon.speed=1;
                end
                this.tAtomsTau(ib, iq, t);
                this.plotAtoms2();
                this.plotBonds2();
                xlim(this.ax, xl);
                ylim(this.ax, yl);
                zlim(this.ax, zl);
                drawnow
                if ~isfield(this.phonon,'save') || isempty(this.phonon.save)
                    pause(0.05)
                else
                    frame = getframe(this.ax);
                    switch this.phonon.save
                        case 'gif'                            
                            im = frame2im(frame);
                            [A,map] = rgb2ind(im,256);
                            if i==0
                                filename=sprintf('%s-%d-%d.gif',this.qe.prefix,ib,iq);
                                imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",0);
                            else
                                imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",0);
                            end
                            if i==8*2-1
                                break;
                            end
                        case 'mp4'
                            if i==0
                                filename=sprintf('%s-%d-%d.mp4',this.qe.prefix,ib,iq);
                                v = VideoWriter(filename, 'MPEG-4');
                                v.FrameRate = 10;
                                open(v);
                            end
                            writeVideo(v, frame);
                            if i==8*6-1
                                close(v);
                                break;
                            end
                    end
                end
            end
        end

        function tAtomsTau(this, ib, iq, t)
            % time dependent atoms position

            % omega=this.phonon.speed*this.matdyn.freq(ib,iq)/max(this.matdyn.freq,[],"all")*pi/3; % speed set to be related with freq
            omega=this.phonon.speed*pi/8; % speed set not to be related with freq
            z=this.matdyn.z(:,ib,iq);
            z=reshape(z,[3,this.matdyn.nat]);
            q=this.matdyn.q(:,iq)*this.qe.tpiba; % 1/Bohr
            for ia=1:this.atoms.nat
                R=this.qe.alat*this.qe.at*this.atoms.supercell(2:4,ia);
                this.atoms.tau(:,ia)=this.atoms.tau_bal(:,ia)+...
                    this.phonon.amplitude*real(z(:,this.atoms.supercell(1,ia))...
                    *exp(1i*(omega*t-q'*R)));
            end

            % update the bonds positions
            this.bonds.formed(:,3:5)=this.atoms.tau(:,this.bonds.formed(:,1)).';
            this.bonds.formed(:,6:8)=this.atoms.tau(:,this.bonds.formed(:,2)).';
        end

    end

    methods (Static, Hidden)
        function idx=findNearestPoint(ax, X, Y, clickX, clickY)
            ax_pos = getpixelposition(ax, true);
            xlim_ = xlim(ax);
            ylim_ = ylim(ax);
            px = (X - xlim_(1)) / (xlim_(2)-xlim_(1));
            py = (Y - ylim_(1)) / (ylim_(2)-ylim_(1));
            Xpix = ax_pos(1) + px * ax_pos(3);
            Ypix = ax_pos(2) + py * ax_pos(4);
            click_px = (clickX - xlim_(1)) / (xlim_(2)-xlim_(1));
            click_py = (clickY - ylim_(1)) / (ylim_(2)-ylim_(1));
            click_Xpix = ax_pos(1) + click_px * ax_pos(3);
            click_Ypix = ax_pos(2) + click_py * ax_pos(4);
            dists = hypot(Xpix - click_Xpix, Ypix - click_Ypix);
            [~, idx] = min(dists);
        end
    end

    methods (Static)
        function setBondLength(a, b, len)
            if ischar(a) || isstring(a)
                a = QE.get_atomic_number(a);
            end
            if ischar(b) || isstring(b)
                b = QE.get_atomic_number(b);
            end
            if ~exist('config', 'dir')
                mkdir('config');
            end
            filename = fullfile('config', 'bondlengthdic.mat');
            if exist(filename, 'file')
                S = load(filename, 'bondlengthdic');
                bondlengthdic = S.bondlengthdic;
            else
                bondlengthdic = nan(118, 118);
            end
            bondlengthdic(a, b) = len;
            bondlengthdic(b, a) = len;
            save(filename, 'bondlengthdic');
        end

        function originmovd=expandat(alat,at,bg,d)
            % d is the expanding thickness, scalar or 3*1 vector
            % at'*bg=I
            d=d(:).';
            if size(d,1)==1
                d=[d d d];
            end
            bgnormed=bg./vecnorm(bg);
            movtmp=-bgnormed.*d;
            interpoint=(movtmp')\(vecnorm(movtmp).^2)'; % movtmp(:,i)'*(interpoint-movtmp(:,i))==0
            originmovd=(at\interpoint)/alat; % direct
        end

        function h=plotBicolorCylinder(v1,v2,c1,c2,r,r1,r2)
            [X,Y,Z] = cylinder(r);
            a = v2 - v1;
            z = a / norm(a);
            if isnan(r1) || isinf(r1)
                v3=(v1+v2)/2;
            else
                v1p=v1+z*r1;
                v2p=v2-z*r2;
                v3=(v1p+v2p)/2;
            end
            normv1=norm(v1-v3);
            normv2=norm(v2-v3);
            % R*[0 0 1]'+T=v2
            % R*[0 0 0]'+T=v3
            % where R must be a rotation matrix, can be represented by R = rotz(alpha)*roty(beta)
            T=v3;

            % tmp can not be parallel with z
            if abs(z(1)) < 0.99
                tmp = [1;0;0];
            else
                tmp = [0;1;0];
            end
            y = cross(z, tmp); y = y/norm(y);
            x = cross(y, z); x = x/norm(x);
            R = [x, y, z];
            n=size(X,2);
            XYZ=[X(1,:) X(2,:); Y(1,:) Y(2,:); Z(1,:) normv2*Z(2,:)];
            xyz=R*XYZ+T;
            x=[xyz(1,1:n);xyz(1,n+1:end)];
            y=[xyz(2,1:n);xyz(2,n+1:end)];
            z=[xyz(3,1:n);xyz(3,n+1:end)];
            h(1)=surf(x,y,z,'FaceColor',c2,'EdgeColor','none');

            XYZ=[X(1,:) X(2,:); Y(1,:) Y(2,:); Z(1,:) -normv1*Z(2,:)];
            xyz=R*XYZ+T;
            x=[xyz(1,1:n);xyz(1,n+1:end)];
            y=[xyz(2,1:n);xyz(2,n+1:end)];
            z=[xyz(3,1:n);xyz(3,n+1:end)];
            h(2)=surf(x,y,z,'FaceColor',c1,'EdgeColor','none');
        end

        function plotBicolorCylinder2(v1,v2,c1,c2,r,r1,r2,h)
            [X,Y,Z] = cylinder(r);
            a = v2 - v1;
            z = a / norm(a);
            if isnan(r1) || isinf(r1)
                v3=(v1+v2)/2;
            else
                v1p=v1+z*r1;
                v2p=v2-z*r2;
                v3=(v1p+v2p)/2;
            end
            normv1=norm(v1-v3);
            normv2=norm(v2-v3);
            % R*[0 0 1]'+T=v2
            % R*[0 0 0]'+T=v3
            % where R must be a rotation matrix, can be represented by R = rotz(alpha)*roty(beta)
            T=v3;

            % tmp can not be parallel with z
            if abs(z(1)) < 0.99
                tmp = [1;0;0];
            else
                tmp = [0;1;0];
            end
            y = cross(z, tmp); y = y/norm(y);
            x = cross(y, z); x = x/norm(x);
            R = [x, y, z];
            n=size(X,2);
            XYZ=[X(1,:) X(2,:); Y(1,:) Y(2,:); Z(1,:) normv2*Z(2,:)];
            xyz=R*XYZ+T;
            x=[xyz(1,1:n);xyz(1,n+1:end)];
            y=[xyz(2,1:n);xyz(2,n+1:end)];
            z=[xyz(3,1:n);xyz(3,n+1:end)];
            set(h(1),'XData',x,'YData',y,'ZData',z,'FaceColor',c2);

            XYZ=[X(1,:) X(2,:); Y(1,:) Y(2,:); Z(1,:) -normv1*Z(2,:)];
            xyz=R*XYZ+T;
            x=[xyz(1,1:n);xyz(1,n+1:end)];
            y=[xyz(2,1:n);xyz(2,n+1:end)];
            z=[xyz(3,1:n);xyz(3,n+1:end)];
            set(h(2),'XData',x,'YData',y,'ZData',z,'FaceColor',c1);
        end

        function h=plotGradientLine(v1,v2,c1,c2,lineWidth)
            c=[c1;c2;NaN NaN NaN];
            h=patch([v1(1) v2(1) NaN],[v1(2) v2(2) NaN],[v1(3) v2(3) NaN],'r','FaceVertexCData',c,'LineWidth',lineWidth,'EdgeColor','interp');
        end

        function h=plotGradientLine2(v1,v2,c1,c2,lineWidth,h)
            c=[c1;c2;NaN NaN NaN];
            set(h, 'XData', [v1(1) v2(1) NaN], 'YData', [v1(2) v2(2) NaN], 'ZData', [v1(3) v2(3) NaN],'FaceVertexCData',c,'LineWidth',lineWidth);
        end

        function RBGtriplet=CPK_coloring(atnum,theme)
            % Reference: https://en.wikipedia.org/wiki/CPK_coloring
            % https://jmol.sourceforge.net/jscolors/
            % Convert atomic number to RGB triplet
            % theme can be 'J'
            switch theme
                case 'J'
                    jmol = [
                        255, 255, 255;
                        217, 255, 255;
                        204, 128, 255;
                        194, 255, 0;
                        255, 181, 181;
                        144, 144, 144;
                        48, 80, 248;
                        255, 13, 13;
                        144, 224, 80;
                        179, 227, 245;
                        171, 92, 242;
                        138, 255, 0;
                        191, 166, 166;
                        240, 200, 160;
                        255, 128, 0;
                        255, 255, 48;
                        31, 240, 31;
                        128, 209, 227;
                        143, 64, 212;
                        61, 255, 0;
                        230, 230, 230;
                        191, 194, 199;
                        166, 166, 171;
                        138, 153, 199;
                        156, 122, 199;
                        224, 102, 51;
                        240, 144, 160;
                        80, 208, 80;
                        200, 128, 51;
                        125, 128, 176;
                        194, 143, 143;
                        102, 143, 143;
                        189, 128, 227;
                        255, 161, 0;
                        166, 41, 41;
                        92, 184, 209;
                        112, 46, 176;
                        0, 255, 0;
                        148, 255, 255;
                        148, 224, 224;
                        115, 194, 201;
                        84, 181, 181;
                        59, 158, 158;
                        36, 143, 143;
                        10, 125, 140;
                        0, 105, 133;
                        192, 192, 192;
                        255, 217, 143;
                        166, 117, 115;
                        102, 128, 128;
                        158, 99, 181;
                        212, 122, 0;
                        148, 0, 148;
                        66, 158, 176;
                        87, 23, 143;
                        0, 201, 0;
                        112, 212, 255;
                        255, 255, 199;
                        217, 255, 199;
                        199, 255, 199;
                        163, 255, 199;
                        143, 255, 199;
                        97, 255, 199;
                        69, 255, 199;
                        48, 255, 199;
                        31, 255, 199;
                        0, 255, 156;
                        0, 230, 117;
                        0, 212, 82;
                        0, 191, 56;
                        0, 171, 36;
                        77, 194, 255;
                        77, 166, 255;
                        33, 148, 214;
                        38, 125, 171;
                        38, 102, 150;
                        23, 84, 135;
                        208, 208, 224;
                        255, 209, 35;
                        184, 184, 208;
                        166, 84, 77;
                        87, 89, 97;
                        158, 79, 181;
                        171, 92, 0;
                        117, 79, 69;
                        66, 130, 150;
                        66, 0, 102;
                        0, 125, 0;
                        112, 171, 250;
                        0, 186, 255;
                        0, 161, 255;
                        0, 143, 255;
                        0, 128, 255;
                        0, 107, 255;
                        84, 92, 242;
                        120, 92, 227;
                        138, 79, 227;
                        161, 54, 212;
                        179, 31, 212;
                        179, 31, 186;
                        179, 13, 166;
                        189, 13, 135;
                        199, 0, 102;
                        204, 0, 89;
                        209, 0, 79;
                        217, 0, 69;
                        224, 0, 56;
                        230, 0, 46;
                        235, 0, 38
                        ]/255;
                    jmol(110:118,1)=jmol(109,1);
                    jmol(110:118,2)=jmol(109,2);
                    jmol(110:118,3)=jmol(109,3);
                    RBGtriplet=jmol(atnum,:);
            end
        end
    end
end