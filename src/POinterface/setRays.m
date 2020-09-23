function [rays] = setRays(Solver_setup, r)
%setRays
%   Usage:
%       [rays] = calcRCS(Const, Solver_setup, theta_grid, phi_grid, xVectors)
%   Input Arguments:
%       Const
%           A global struct, containing general data
%       Solver_setup
%           Solver specific struct, e.g. frequency range, basis function details, geometry details
%       r
%           Distance to viewport
%
%   Output Arguments:
%       rays
%           Structure containing the origin and direction of each ray
%
%   Description:
%       Sets up a uniform grid of rays at the viewport directed in the
%       incident direction
%
%   =======================
%   Written by Cullen Stewart-Burger on 1 Sept 2020
%   Stellenbosch University
%   Email: 20751028@sun.ac.za

%get U and V vectors for incident plane
%a_phi = [-sin(Solver_setup.phi*Const.DEG2RAD), cos(Solver_setup.phi*Const.DEG2RAD), 0];
%a_theta = [cosd(Solver_setup.phi)*cosd(Solver_setup.theta), cosd(Solver_setup.phi)*sind(Solver_setup.theta), -sind(Solver_setup.theta)];

%For now we project all the centre points onto the plane
pt = Solver_setup.triangle_centre_point;
[norm(1), norm(2), norm(3)] = sph2cart(Solver_setup.phi, pi/2-Solver_setup.theta*pi/810, r);
norm = -norm;
N = size(pt,1);
norm_mat = repmat(norm,[N 1]);
dist = dot(pt, norm_mat, 2)+r;
rays.origin = pt - dist.*norm;
rays.dir = norm_mat;

end