function HfieldAtPointCartesian =  calculateHfieldAtPointRWGCart(Const, Px, Py, Pz, ...
    Solver_setup, Isol)
    %calculateEfieldAtPointRWG
    %   Usage:
    %       [EfieldAtPointCartesian] =  calculateHfieldAtPointRWGCart(Const, Px, Py, Pz, ...
    %           Solver_setup, Isol)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Px, Py, Pz
    %           Cartesian co-ordinate of the observation point.
    %       Solver_setup (struct)
    %           Struct containing the frequency range, nodal co-ordinates of the 
    %           triangles, as well as information on the shared edges, RWG basis functions, 
    %           etc. 
    %       Isol
    %           The solution vector, i.e. the MoM expansion coefficients
    %   Output Arguments:
    %       EfieldAtPointCartesian
    %           The H-field strength calculated at the point (Px, Py, Pz)
    %           EfieldAtPoint = A x^ + B y^ + C z^
    %
    %   Description:
    %       Calculates the H-field value at a certain point
    %       (r,theta,phi). The H-field value is based on the dipole
    %       model, as discussed in [1] and [2].
    %
    %   TO-DO: Z-displaced atennas is not working correctly
    %
    %   References:
    %   [1] S. Makarov, "MoM Antenna Simulations with Matlab: RWG Basis Functions"
    %   [2] Balanis, "Antenna Theory: Analysis and Design (3rd Edition)"
    %   =======================
    %   Adapted from 'calculateEfieldAtPointRWG.m' written by Danie Ludick on 2018.05.04
    %   Adapted by Cullen Stewart-Burger (20751028@sun.ac.za)
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    % Initialisations
    % Note: for each frequency, we will have a different field
    HfieldAtPointCartesian = zeros(Solver_setup.frequencies.freq_num, 3);

    for freq_index = 1:Solver_setup.frequencies.freq_num

        % Wavelength
        lambda = Const.C0/Solver_setup.frequencies.samples(freq_index);
        k = 2*pi/lambda;
        K = 1j*k/(4*pi);
        %K = 1.0;

        % The E-field at the point is calculated as the superposition of the individual 
        % RWG (dipole) contributions. See definitions for m, M and C (Eqs. 9a and 9b)
        for rwg_bf_index = 1:Solver_setup.num_metallic_edges

            % --------------------------------
            % calculate m (the dipole moment vector) associate with the mth edge: 
            %     m = lm*Im*[ (r_mc^-) - (r_mc^+)]
            % where "c" denotes the triangle midpoint.
            % --------------------------------
            lengthM = Solver_setup.rwg_basis_functions_length_m(rwg_bf_index);
            Im = Isol(rwg_bf_index);
            if (Im == 0)
                continue
            end

            triangle_plus = Solver_setup.rwg_basis_functions_trianglePlus(rwg_bf_index);
            triangle_minus = Solver_setup.rwg_basis_functions_triangleMinus(rwg_bf_index);

            rMplusX  = Solver_setup.triangle_centre_point(triangle_plus,1);
            rMminusX = Solver_setup.triangle_centre_point(triangle_minus,1);
            
            rMplusY  = Solver_setup.triangle_centre_point(triangle_plus,2);
            rMminusY = Solver_setup.triangle_centre_point(triangle_minus,2);
            
            rMplusZ  = Solver_setup.triangle_centre_point(triangle_plus,3);
            rMminusZ = Solver_setup.triangle_centre_point(triangle_minus,3);

            rwgDipoleCentreXYZ  = lineCentre([rMplusX, rMplusY, rMplusZ],[rMminusX, rMminusY, rMminusZ]);
            rwgDipoleCentre = rwgDipoleCentreXYZ;
            
            rL = sqrt( (rwgDipoleCentre(1) - Px)^2 + (rwgDipoleCentre(2) - Py)^2 + ...
                (rwgDipoleCentre(3) - Pz)^2 );
            
            rLvecX = Px - rwgDipoleCentreXYZ(1);
            rLvecY = Py - rwgDipoleCentreXYZ(2);
            rLvecZ = Pz - rwgDipoleCentreXYZ(3);
            rLvec = [rLvecX, rLvecY, rLvecZ];        
                                                            
            C = (1/(rL*rL)) * (1 + 1/(1j*k*rL));

            mX = lengthM*Im*(rMminusX - rMplusX);
            mY = lengthM*Im*(rMminusY - rMplusY);
            mZ = lengthM*Im*(rMminusZ - rMplusZ);
            m = [mX, mY, mZ];


            % calculate the H-field (x,y,z) due to this RWG Basis Function. Note, the phase is referenced 
            % relative to the centre of the dipole. Account for that below.
            phaseTerm = exp(-1j*k*rL);  
            
            HfieldCurrentRWG = K * cross(m, rLvec) * C * phaseTerm;

            HfieldAtPointCartesian(freq_index, :) = HfieldAtPointCartesian(freq_index, :) + HfieldCurrentRWG;
        end

    end%for freq_index = 1:Solver_setup.freq_num


