function [po] = runPOsolver(Const, Solver_setup, ~, ~, refIsol)

    %runPOsolver
    %   Usage:
    %       [po] = runPOsolver(Const)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details    
    %       zMatrices
    %           The Z-matrices data. This can be from FEKO (extracted from the *.mat file, or internally
    %           calculated).
    %       yVectors
    %           The Yrhs-vector data
    %       refIsol
    %           The reference solution-vector data (e.g. MoM solution of FEKO or SUN-EM)
    %
    %   Output Arguments:
    %       po
    %           Structs containing PO solution and timing data
    %
    %   Description:
    %       Runs the PO solution based on the Z and Y data that was read / parsed
    %       from the FEKO *.out, *.mat, *.str and *.rhs files
    %
    %   =======================
    %   Written by Cullen Stewart-Burger on 27 July 2020
    %   Adapted from runMoMsolver.m by Danie Ludick
    %   Stellenbosch University
    %   Email: 20751028@sun.ac.za
    
    %   indev notes (28/07/2020):
    %       Current implementation assumes fully illuminated target
    %       Excitation source is set to a unit phi-polarised plane wave
    %       travelling in -r direction.

    narginchk(5,5);

    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Running PO solver'));

    % Initialise the return values
    po  = [];
    po.name = 'po';
    Npo = Solver_setup.num_mom_basis_functions;   % Total number of basis functions for whole problem
    %numSols = xVectors.numSols;                   % The number of solutions configurations
    po.numSols = 1; %numSols;                     % For now, set to 1. (TO-DO: Update)
    numFreq = Solver_setup.frequencies.freq_num;   % The number of frequency points to process
    numRHSperFreq = po.numSols / numFreq;         % The number of solutions per frequency point.
                                                   % For now, should be 1 (TO-DO: Update)
                                                   
    % Some info about the solution configurations
    message_fc(Const,sprintf('  numSols : %d', po.numSols));
    message_fc(Const,sprintf('  numFreq : %d', numFreq));
    message_fc(Const,sprintf('  numRHSperFreq : %d', numRHSperFreq));

    % Calculate the solution vector (all frequency points, all RHSes)
    po.Isol = complex(zeros(Npo,po.numSols));

    % The timing calculations also need to take into account that there is a
    % frequency loop
    po.setupTime = zeros(1,numFreq);
    % Zero also the total times (for all frequency iterations)
    po.totsetupTime = 0.0;
    po.totfactorisationTime = 0.0;
    po.totsolTime = 0.0;

    % Start the frequency loop now
    for freq=1:numFreq

        % Start timing (PO setup)
        tic

        % Reset each solution per frequency point here
        solStart = 1;
        solEnd   = numRHSperFreq;

%         % Allocate here space for the MoM matrix so that it can be filled
%         ObservRWGs = [1:Npo];
%         SourceRWGs = [1:Npo];
%         % Note: Since 2017-06-25, we are also passing a freq. variable here to
%         % indicate for which frequency point we are extracting the matrix
%         Zmom = (calcZmn(Const, zMatrices, freq, 1, 1, ObservRWGs, SourceRWGs));

        % End timing (calculating the impedance matrix)
        po.setupTime(freq) = toc;

        % Start timing (MoM factorisation)
        tic
        % LU-decomposition of the Z-matrix
        %[L,U] = lu(Zmom);

        % Loop over the RHs vectors (numSols) and calculate each of the currents.
        for index=solStart:Npo

            % Take care where we are extracting the values from (out of Yrhs)
            % and also where we will be storing these values again (in Xsol)
            %index = solNum + (freq-1)*numRHSperFreq;
            % DJdbg --> remove
            %message_fc(Const,sprintf('  index: %d', index));

            % Back-wards substitution
            %b = L\yVectors.values(:,index);
            rn = Solver_setup.rwg_basis_functions_shared_edge_centre(index, :);
            shared_nodes = Solver_setup.rwg_basis_functions_shared_edge_nodes(index, :);
            ln = Solver_setup.nodes_xyz(shared_nodes(1), :) - Solver_setup.nodes_xyz(shared_nodes(2), :);
            ln = ln/Solver_setup.rwg_basis_functions_length_m(index);
            %We do not know if this direction for ln is correct according
            %to our reference. This is checked and possibly corrected below.
            
            side = Solver_setup.nodes_xyz(shared_nodes(2), :)- Solver_setup.nodes_xyz(Solver_setup.rwg_basis_functions_trianglePlusFreeVertex(index), :);
            normTest = cross(side, ln);
            reverse = dot(normTest, Solver_setup.triangle_normal_vector(Solver_setup.rwg_basis_functions_trianglePlus(index), :));
            if(reverse < 0)
                ln = -ln;
            end
            
            
            delta = 1;  %currently assume model is fully illuminated
            k = 2*pi*500E6/3E8; %set 1000MHz freq
            [kx, ky, kz] = sph2cart(Solver_setup.phi*Const.DEG2RAD, (90-Solver_setup.theta)*Const.DEG2RAD, -k);
            k_vec = [kx, ky, kz];
            H = (1/Const.ETA_0)*exp(j*dot(k_vec, rn));    %Find impressed H field (r directed plane wave with E field theta-polarised)
            a_phi = [-sin(Solver_setup.phi*Const.DEG2RAD), cos(Solver_setup.phi*Const.DEG2RAD), 0];
            H_vec = -H*a_phi;
            po.Isol(index) = dot(2*delta*H_vec, ln);
        end%for

        % End timing (MoM factorisation)
        po.factorisationTime(freq) = toc;

        % Total time (MoM matrix setup + factorisation)
        po.solTime(freq) = po.setupTime(freq) + po.factorisationTime(freq);

        % Memory usage remains constant between frequency iterations
        po.memUsage = 0;%byteSize(Zmom);

        % Calculate the total MoM times
        po.totsetupTime = po.totsetupTime + po.setupTime(freq);
        po.totfactorisationTime = po.totfactorisationTime + po.factorisationTime(freq);
        po.totsolTime = po.totsolTime + po.solTime(freq);

    end%for freq=1:numFreq

    message_fc(Const,sprintf('Finished MoM solver in %f sec.',po.totsolTime));
    message_fc(Const,sprintf('(Times for Z-setup : %f sec. and LU-fact. : %f)',po.totsetupTime, ...
        po.totfactorisationTime));
    message_fc(Const,sprintf('Memory usage of PO %s',po.memUsage));

    % Compare the MoM solution obtained with MATLAB, with that obtained by FEKO
    % that was stored in xVectors.values (for each frequency iteration (and each solution within the frequency iteration)
    % Calculate also space for the relative error here
    po.relError = zeros(1,po.numSols);
    for freq=1:numFreq
        for solNum=1:numRHSperFreq
            index = solNum + (freq-1)*numRHSperFreq;
            po.relError(index) = calculateErrorNormPercentage(refIsol.Isol(1:Solver_setup.num_metallic_edges,index), po.Isol(:,index));
            message_fc(Const,sprintf('Rel. error norm. for Sol. %d of %d of freq. %d of %d compared to reference sol. %f percent',solNum, ...
                numRHSperFreq, freq, numFreq, po.relError(index)));
        end
    end
    
    % Write the MoM solution to a ASCII str file, so that it can be read
    % again by FEKO (for plotting in POSTFEKO) - only if requested (i.e. if the filename is defined)
    if (~isempty(Const.SUNEMmomstrfilename))
        writeSolToFile(Const, po);
    end%if

    