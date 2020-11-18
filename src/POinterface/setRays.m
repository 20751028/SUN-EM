function [rays] = setRays(Solver_setup, theta_degrees, phi_degrees)
%setRays
%   Usage:
%       [rays] = setRays(Solver_setup, theta_degrees, phi_degrees)
%   Input Arguments:
%       Solver_setup
%           Solver specific struct, e.g. frequency range, basis function details, geometry details
%       theta_degres, phi_degrees
%           spherical coordinates for the observation angle
%
%   Output Arguments:
%       rays
%           Structure containing the origin and direction of each ray
%
%   Description:
%       Sets up a uniform grid of rays at the viewport directed from
%       (theta_degreed, phi_degrees).
%
%   =======================
%   Written by Cullen Stewart-Burger on 1 Sept 2020
%   Stellenbosch University
%   Email: 20751028@sun.ac.za

%For now we project all the centre points onto the plane
pt = Solver_setup.triangle_centre_point;
[norm(1), norm(2), norm(3)] = sph2cart(phi_degrees*pi/180, pi/2-theta_degrees*pi/180, 1);
norm = -norm;
N = size(pt,1);
norm_mat = repmat(norm,[N 1]);
dist = dot(pt, norm_mat, 2)+Solver_setup.r_bounding+0.1;
rays.origin = pt - dist.*norm;
rays.dir = norm_mat;

end