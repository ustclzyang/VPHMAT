classdef MATDYN < handle
    %MATDYN Parse and store matdyn.f90's data
    %
    %   For more information about Quantum ESPRESSO, visit:
    %     https://www.quantum-espresso.org/
    %   For more information about matdyn.f90, visit:
    %     https://gitlab.com/QEF/q-e/-/blob/develop/PHonon/PH/matdyn.f90
    %   See also: QE, NIST, VPHMAT
    %
    %   Copyright (c) 2025 Lizhou Yang
    %   Licensed under the MIT License.

    properties
        nat
        nq
        q         % 3*nq tpiba
        freq      % nat3*nq Ry
        freq_thz  % nat3*nq THz
        freq_cmm1 % nat3*nq cm^{-1}
        unit      % 'ry' 'thz' 'cmm1'
        klabels   % {"\Gamma","X","S","Y","\Gamma"}
        eigvec
        z
        ax        % phonon dispersion gca handle
    end

    properties (Hidden)
        hLines    % nbnd*1 ploted lines' handle
        hpoint    % hight light point's handle
    end

    methods
        function this = MATDYN(flvec)
            % read matdyn.modes
            % see matdyn.f90: "The normalized phonon displacements
            %      are the eigenvectors divided by the square root of the mass,
            %      then normalized. As such they are not orthogonal."
            % see write_eigenvectors.f90: subroutine writemodes (nat,q,w2,z,iout)

            fid = fopen(flvec, 'r');
            if fid == -1
                error('MATDYN:OpenFileError', ...
                    'Cannot open this matdyn.modes file: %s', flvec);
            end
            cleaner = onCleanup(@() fclose(fid));

            nq=0;
            nat3=0;
            while ~feof(fid)
                tmp=fgetl(fid);
                if contains(tmp,'q =')
                    nq=nq+1;
                end
                if nq==1 && contains(tmp,'freq')
                    nat3=nat3+1;
                end
            end
            frewind(fid);
            nat=nat3/3;
            this.nat=nat;
            this.nq=nq;
            q=zeros(3,nq);
            ztmp=zeros(6,nat,nat3,nq);
            z=zeros(nat3,nat3,nq); % complex, normalized, the first index is atoms and its direction, the second index is phonon branch
            freq_ind=zeros(nat3,nq);
            freq_thz=zeros(nat3,nq);
            freq_cmm1=zeros(nat3,nq);
            for iq=1:nq
                fgetl(fid);
                fgetl(fid);
                q(:,iq)=fscanf(fid,' q = %f %f %f\n',3);
                fgetl(fid);
                for i=1:nat3
                    tmp=fscanf(fid,' freq (  %f) = %f [THz] = %f [cm-1] \n',3);
                    freq_ind(i,iq)=tmp(1);
                    freq_thz(i,iq)=tmp(2);
                    freq_cmm1(i,iq)=tmp(3);
                    ztmp(:,:,i,iq)=fscanf(fid,' ( %f %f %f %f %f %f ) \n',[6,nat]);
                    tmp=complex(ztmp([1 3 5],:,i,iq),ztmp([2 4 6],:,i,iq));
                    z(:,i,iq)=tmp(:);
                end
                fgetl(fid);
            end
            this.z=z;
            this.q=q;
            ry_to_thz=1e-12*NIST.Rydberg_constant_in_J/NIST.Planck_constant;
            % ry_to_cmm1=1e-2*NIST.Rydberg_constant_in_J/NIST.Planck_constant/NIST.speed_of_light;
            this.freq=freq_thz/ry_to_thz; % Ry
            % this.freq=freq_cmm1/ry_to_cmm1; % Ry
            this.freq_thz=freq_thz; % THz
            this.freq_cmm1=freq_cmm1; % cm^{-1}
            this.unit='cmm1';
            this.klabels={};


        end

        function plot(this)
            switch this.unit
                case 'ry'
                    [this.ax, this.hLines] = MATDYN.plotband(this.q,this.freq,this.klabels);
                    ylabel('Ry');
                case 'thz'
                    [this.ax, this.hLines] = MATDYN.plotband(this.q,this.freq_thz,this.klabels);
                    ylabel('THz');
                case 'cmm1'
                    [this.ax, this.hLines] = MATDYN.plotband(this.q,this.freq_cmm1,this.klabels);
                    ylabel('cm^{-1}');
            end
        end
    end

    methods(Static)
        function [ax, hLines]=plotband(q,freq,klabels)
            %PLOTBAND Plot band or phonon dispersion.
            % PLOTBAND(q,freq,klabels) plot freq(nat*3,nq) versus q(3,nq),
            % klabels (optional) can be {"\Gamma","X","S","Y","\Gamma"}
            if nargin < 3
                klabels={};
            end
            hold on

            nks=size(q,2);
            nbnd=size(freq,1);
            ks=q.';

            ksd=[0 0 0;diff(ks)]; % vector formed by nearby k points
            ksdd=diff(ksd);
            xksdd=~all(abs(ksdd)<0.001,2);
            xksdd=[xksdd;true]; % high symmetry points index
            xksdd_tmp=[false; xksdd(1:end-1)];
            xksdd_discon=xksdd & xksdd_tmp; % this one is discontinue with last one
            n_discon=sum(xksdd_discon); % number of discontinue points
            high_symmetry_point=ks(xksdd,:);
            % disp(high_symmetry_point)
            xksd=vecnorm(ksd,2,2);
            xksd(xksdd_discon)=0;
            xks=cumsum(xksd,1);

            high_symmetry_point_x=xks(xksdd);

            label_index=cumsum(xksdd,1);
            i_discon=label_index(xksdd_discon); % among them, i_discon(i) and i_discon(i)-1 is discon，should merge their klabels
            high_symmetry_point_xd=high_symmetry_point_x(1:end-n_discon);
            j=1;
            for i=1:size(high_symmetry_point_x,1)
                if any(i==i_discon-1)
                    high_symmetry_point_xd(j)=high_symmetry_point_x(i);
                    j=j+1;
                elseif ~any(i==i_discon)
                    high_symmetry_point_xd(j)=high_symmetry_point_x(i);
                    j=j+1;
                end
            end
            xticks(high_symmetry_point_xd)

            if ~isempty(klabels)
                if size(high_symmetry_point_x,1)~=size(klabels,2)
                    error('there are %d high symmetry points, but %d klabels provided', size(high_symmetry_point_x,1),size(klabels,2));
                end
                klabelsd=klabels(1:end-n_discon);
                j=1;
                for i=1:size(high_symmetry_point_x,1)
                    if any(i==i_discon-1)
                        klabelsd{j}=join([klabels{i},'/',klabels{i+1}],"");
                        j=j+1;
                    elseif ~any(i==i_discon)
                        klabelsd{j}=klabels{i};
                        j=j+1;
                    end
                end
                xticklabels(klabelsd)
            else
                xticklabels({})
            end

            hLines = gobjects(nbnd,1);
            for ib=1:nbnd
                hLines(ib) = plot(xks,freq(ib,:),'Color','k');
            end
            axis tight
            xlim([0 xks(end)+1e-5])
            ax=gca;
        end
    end
end