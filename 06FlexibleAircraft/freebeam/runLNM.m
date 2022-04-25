%
% runLNM.m: Compute mode shapes in quasi-coordinates.
%
% Copyright, Rafael Palacios, June 2018
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all, close all
 global invM0 Kss m IB rcgcross Mrx Mpx Mxx F 
 global Nflexmodes Tmax omega Fmodal
 global invMrigid
 
% General problem parameters.
L= 10;          % Length of the beam.
EI=2;           % Bending stiffness
rho_s=1;        % Mass density per unit length (assume A=1).
N= 25;          % Number of elements.
rBP=0:L/N:L;    % Coordinates of the beam nodes.
adofs=[3:2*(N+1)]; % Active degrees of freedom (after BCs are enforced).
Nflexmodes=4;

% Define the problem.
mod_setup

% Linear normal modes
mod_lnm

% eof