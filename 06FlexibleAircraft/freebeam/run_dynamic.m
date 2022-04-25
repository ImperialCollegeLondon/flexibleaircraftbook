%
% run_dynamic.m: Run dynamic simulation in quasi-coordinates.
%
% Copyright, Rafael Palacios, June 2018
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all
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
Fmax=0.25;      % Value of the applied force.
Tmax=5;         % Time for linear ramp.
tfinal=20;      % Integration time.
dt=0.01;

% Define the problem.
mod_setup

% Linear normal modes
mod_lnm

% Dynamic simulations -- small elastic deformations.
mod_dynamic

% eof
