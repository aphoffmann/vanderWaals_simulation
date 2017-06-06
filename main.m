 %************************************************************  
 % Programmer: Alex Hoffmann                                *
 % Date: 5/2/17                                             *
 % Name: main.m                                             *
 % Description: 1. Initialize 2d lattice of atoms and plot, *
 % compute potential energy and kinetic energy, simulate    *
 % atomic movement, summarizes results                      *
 %***********************************************************
 function main()
    %Contains all code outlined above. No input or output.
    
    %% Part 1 
    % Allocate memory for x and y
    x = zeros(1,25);
    y = zeros(1,25);
    
    %Initialize x and y values with 4loop
    for i = 1:5
       x((5*i-4):(5*i)) = ones(1, 5)*3.8*i-11.4;
       y((5*i-4):(5*i)) = -7.6:3.8:7.6;
    end
    
    %Plot atoms' location
    figure(1);
    scatter(x,y)
    title('2D Argon Atom System')
    axis([-10 10 -10 10])
    xlabel('X-coordinate (Angstroms)')
    ylabel('Y-coordinate (Angstroms)')
    
    %% Part 2 %%
    %Compute potential energy and vector forces using vanderWaals
    [energy, Fx, Fy] = vanderWaals(x,y);
    fprintf('Initial Potential Energy: %f\n', energy)
    
    %Initialize random velocites
    x_vel = 0.05*randn(1,25);
    y_vel = 0.05*randn(1,25);
    
    
    %% Part 3
    %Allocate memory for potential and kinetic energy
    potential = zeros(1,5000);
    kinetic = zeros(1,5000);
    
    for t = 1:5000
        %3.1 compute new x and y values
        x = x + x_vel*.01 +.5*Fx*.01^2;
        y = y + y_vel*.01 +.5*Fy*.01^2;
        
        %3.2 Compute new potential, velocity vectors, and vector forces
        tempFx = Fx;
        tempFy = Fy;
        [energy, Fx, Fy] = vanderWaals(x,y);
        x_vel = x_vel + (tempFx + Fx).*.005;
        y_vel = y_vel + (tempFy + Fy).*.005;
        
        % 3.3 Record Potential and Kinetic Energies
        kinetic(t) = .5*sum(x_vel.^2+y_vel.^2);
        potential(t) = energy; 
        
        
    end
    total_energy = kinetic+potential;
    fprintf('Energy: %f\n', energy);
    
    %% Part 4
    % Plot new argon atom formation
    figure(2);
    scatter(x,y)
    title('2D Argon Atom System (Final)')
    axis([-10 10 -10 10])
    xlabel('X-coordinate (Angstroms)')
    ylabel('Y-coordinate (Angstroms)')
    
    % Calculate means and standard deviations
    potential_mean = mean(potential(end/2:end));
    potential_std = std(potential(end/2:end));
    
    kinetic_mean = mean(kinetic(end/2:end));
    kinetic_std = std(kinetic(end/2:end));
    
    total_mean = mean(total_energy(end/2:end));
    total_std = std(total_energy(end/2:end));
    
    %Plot Kinetic and potential energy variation over 5000 steps
    figure(3)
    plot(1:5000, potential);
    title('Argon Energy') % Labels
    ylabel('Energy (Kcal/mole)')
    xlabel('MD Step)')
    axis([0 5000 -20 5]) %Set Axes
    hold on;
    plot(1:5000, kinetic) % Add other data to plot
    plot(1:5000, total_energy)
    %Legend! using sprintf
    legend(sprintf('Potential Energy: %f +/- %f', potential_mean, potential_std), ...
        sprintf('Kinetic Energy: %f +/- %f', kinetic_mean,kinetic_std),...
        sprintf('Total Energy: %f +/- %f', total_mean,total_std))
    
 end
 
 
