%%
% Author:   MJ Schnieders
% Date:     April 19th, 2013
% Function: vanDerWaalsForces.m
% Purpose:  Compute the van der Waals energy and forces for Argon atoms.
% Input:    
%           x - a vector of x-coordinates (Angstroms)
%           y - a vector of y-coordinates (Angstroms)
% Output:  
%           energy - the potential energy (Kcal/mole)
%           fx     - the x-component of the force (Kcal/mole/Angstrom)
%           fy     - the y-component of the force (Kcal/mole/Angstrom)
%
function [energy, fx, fy] = vanderWaals(x, y)
    eps4 = 0.2824 * 4.0;
    rmin = 3.361;
    rmin6 = rmin^6;
    rmin12 = rmin^12;
    nAtoms = length(x);
    % Zero out energy.
    energy = 0;
    % Zero out forces.
    fx = zeros(1,nAtoms);
    fy = zeros(1,nAtoms);
    % Loop over atoms to accumulate energy and forces.
    for i=1:nAtoms
        xi = x(i);
        yi = y(i);
        % Loop over atoms whose index is greater than i.
        for j=i+1:nAtoms
            % Compute separation distance (Angstroms).
            dx = xi - x(j);
            dy = yi - y(j);
            r = sqrt(dx*dx + dy*dy);
            % Compute interaction energy.
            ir = 1.0/r;
            ir6 = ir^6;
            ir12 = ir6*ir6;
            e = eps4*(rmin12*ir12 - rmin6*ir6);
            energy = energy + e;
            % Compute equal and opposite forces, includeing chain rule term.
            ir7 = ir6*ir;
            ir13 = ir12*ir;
            de = eps4*(6.0*rmin6*ir7 - 12.0*rmin12*ir13);
            dxr = de*dx/r;
            dyr = de*dy/r;
            fx(i) = fx(i) - dxr;
            fy(i) = fy(i) - dyr;
            fx(j) = fx(j) + dxr;
            fy(j) = fy(j) + dyr;
        end
    end
end