% Aaron Wold
% 3-D Elastic Collision Simulator [Fission]

clc; clear;

% All distances are measured in pm
% All velocities are measured in pm/s
% Mass is in amu

particles = [];

% Boundary Properties
boundary = 1000;
xboundary.max = boundary;
xboundary.min = -boundary;

yboundary.max = boundary;
yboundary.min = -boundary;

zboundary.max = boundary;
zboundary.min = -boundary;

% Neutron Properties

%{
neutron.radius = 100; % Neutron radius in picometers
[nx, ny, nz] = sphere(num_faces); % Generate coords for sphere

neutron_initial_pos = [500 0 0]; % Neutron initial position
neutron.pos = neutron_initial_pos; % Neutron current position

neutron_initial_velocity = [-500 0 -150]; % Neutron initial velocity
neutron.velocity = neutron_initial_velocity; % Neutron current velocity

neutron.mass = 1; % Neutron mass (unitless for now)
%}

neutron = generate_particle("neutron", 5, 1.008, [randi([xboundary.min xboundary.max]) randi([yboundary.min yboundary.max]) randi([zboundary.min zboundary.max])], [randi([2500 10000]) randi([2500 10000]) randi([2500 10000])]);
particles = [particles neutron]; % Added neutron to list of particles

% U-235 Properties

%{
uranium.radius = 100; % Uranium radius in picometers
[ux, uy, uz] = sphere(num_faces); 

uranium_initial_pos = [0 0 0]; % Uranium initial position
uranium.pos = uranium_initial_pos; % Uranium current position

uranium_initial_velocity = [0 0 0]; % Uranium initial velocity
uranium.velocity = uranium_initial_velocity; % Uranium current velocity

uranium.mass = 2;
%}

uranium = generate_particle("uranium", 156, 235, [0, 0, 0], [0, 0, 0]);
particles = [particles uranium];

% Time
dt = 0.00001; % Time step
t = 0;

% Figure Initialization
figure(Name="Neutron Collision")
frame_counter = 0;

% Simulation loop
while 1
    t = t + dt;
    % Calculate distance between centers of each sphere
    d = uranium.pos - neutron.pos;
    dmag = sqrt(sum(d.^2));
    l.mag = dmag - neutron.radius - uranium.radius; % Distance between sphere surfaces

    n = d/dmag; % Normal vector
    

    % Figure Animation
    if mod(frame_counter, 3000) == 1
         neutron_final = surf(neutron.px*neutron.radius+neutron.pos(1), neutron.py*neutron.radius+neutron.pos(2), neutron.pz*neutron.radius+neutron.pos(3), FaceColor= [0 0 1]);
         hold on
         uranium_final = surf(uranium.px*uranium.radius+uranium.pos(1), uranium.py*uranium.radius+uranium.pos(2), uranium.pz*uranium.radius+uranium.pos(3), FaceColor= [0 1 0]);
         hold off
        
        title(sprintf("t = %.1f", t))
        %view(-45, 45);
        
        view(frame_counter/10000,45)

        xlim([xboundary.min xboundary.max])
        ylim([yboundary.min yboundary.max])
        zlim([zboundary.min zboundary.max])

        %xlim([neutron.pos(1)-200 neutron.pos(1)+200])
        %ylim([neutron.pos(2)-200 neutron.pos(2)+200])
        %zlim([neutron.pos(3)-200 neutron.pos(3)+200])

        drawnow
    end

    % Check for Collisions

    % w/ Walls
    %{
    for i = 1:size(particles, 2)
        if particles(i).pos(1) > xboundary.max-particles(i).radius
            particles(i).velocity(1) = -particles(i).velocity(1);
        elseif particles(i).pos(1) < xboundary.min+particles(i).radius
            particles(i).velocity(1) = -particles(i).velocity(1);
        end

        if particles(i).pos(2) > yboundary.max-particles(i).radius
            particles(i).velocity(2) = -particles(i).velocity(2);
        elseif particles(i).pos(2) < yboundary.min+particles(i).radius
            particles(i).velocity(2) = -particles(i).velocity(2);
        end

        if particles(i).pos(3) > zboundary.max-particles(i).radius
            particles(i).velocity(3) = -particles(i).velocity(3);
        elseif particles(i).pos(3) < zboundary.min+particles(i).radius
            particles(i).velocity(3) = -particles(i).velocity(3);
        end
    end
    %}
    
    if uranium.pos(1) > xboundary.max-uranium.radius || uranium.pos(1) < xboundary.min+uranium.radius
        uranium.velocity(1) = -uranium.velocity(1);
    end
    if uranium.pos(2) > yboundary.max-uranium.radius || uranium.pos(2) < yboundary.min+uranium.radius 
        uranium.velocity(2) = -uranium.velocity(2);
    end
    if uranium.pos(3) > zboundary.max-uranium.radius || uranium.pos(3) < zboundary.min+uranium.radius
        uranium.velocity(3) = -uranium.velocity(3);
    end

    if neutron.pos(1) > xboundary.max-neutron.radius || neutron.pos(1) < xboundary.min+neutron.radius
        neutron.velocity(1) = -neutron.velocity(1);
    end
    if neutron.pos(2) > yboundary.max-neutron.radius || neutron.pos(2) < yboundary.min+neutron.radius 
        neutron.velocity(2) = -neutron.velocity(2);
    end
    if neutron.pos(3) > zboundary.max-neutron.radius || neutron.pos(3) < zboundary.min+neutron.radius
        neutron.velocity(3) = -neutron.velocity(3);
    end
    

    % w/ Neutrons
    if  l.mag <= 1E-5 % With neutrons
        fprintf("COLLISION @ t = " + t + "s \n")
        
        % Account for sphere overlap
        for i = 1:size(d,2)
            if d(i) ~= 0
                neutron.pos(i) = neutron.pos(i) + l.mag;
            end
        end

        % Conservation of momentum and kinetic energy system
        
        temp = uranium.velocity;

        v_relative = uranium.velocity - neutron.velocity;
        v_normal = (dot(n, v_relative))*n;

    %    uranium.velocity = (uranium.mass-neutron.mass)/(uranium.mass+neutron.mass)*uranium.velocity + (2*neutron.mass)/(uranium.mass+neutron.mass)*neutron.velocity; %v_normal;
    %    neutron.velocity = (neutron.mass-uranium.mass)/(neutron.mass+uranium.mass)*neutron.velocity + (2*uranium.mass)/(uranium.mass+neutron.mass)*temp; %v_normal;

        uranium.velocity = uranium.velocity - (2*neutron.mass)/(uranium.mass + neutron.mass)*v_normal;
        neutron.velocity = neutron.velocity + (2*uranium.mass)/(uranium.mass + neutron.mass)*v_normal;

        fprintf("Neutron Velocity: [" + num2str(neutron.velocity) + "]\n");
        fprintf("Uranium Velocity: [" + num2str(uranium.velocity) + "]\n");
        
    end

    % Calculate new sphere positions
    uranium.pos = uranium.pos + dt*uranium.velocity;
    neutron.pos = neutron.pos + dt*neutron.velocity;

   
    frame_counter = frame_counter + 1;
end

% Output data
fprintf("Neutron Center Position: [" + num2str(neutron.pos) + "]\n");
fprintf("Uranium Center Position: [" + num2str(uranium.pos) + "]\n");


