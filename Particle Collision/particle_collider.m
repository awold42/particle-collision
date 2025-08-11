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

for i = 1:20
    neutron = generate_particle("neutron", 5, 1.008, [randi([-750 750]) randi([-750 750]) randi([-750 750])], [randi([0 2000]) randi([0 2000]) randi([0 2000])]);
    particles = [particles neutron]; % Added neutron to list of particles
end

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

for i = 1:10
    uranium = generate_particle("uranium", 20, 235, [randi([-750 750]) randi([-750 750]) randi([-750 750])], [0 0 0]); % r = 156
    particles = [particles uranium];
end
% Time
dt = 0.01; % Time step
t = 0;

% Figure Initialization
figure(Name="Neutron Collision")
frame_counter = 0;

collision = [];

% Simulation loop
while 1
    t = t + dt;

    % Calculate distance between centers of each sphere
    for i = 1:size(particles, 2)
        for j = 1:size(particles, 2)
            d = particles(i).pos - particles(j).pos;
            if d == 0
                break;
            end
            dmag = sqrt(sum(d.^2));
            lmag = dmag - particles(i).radius - particles(j).radius; % Distance between sphere surfaces

            % Check for collision with another particle
            if  lmag <= 1E-4
                
                collision = particles(i).name + particles(j).name;
                fprintf("COLLISION @ t = " + t + "s \n")
                
                n = d/dmag; % Normal vector
        
                
        
                % Conservation of momentum and kinetic energy system
        
                v_relative = particles(i).velocity - particles(j).velocity;
                v_normal = (dot(n, v_relative))*n;
        
                particles(i).velocity = particles(i).velocity - (2*particles(j).mass)/(particles(i).mass + particles(j).mass)*v_normal;
                particles(j).velocity = particles(j).velocity + (2*particles(i).mass)/(particles(i).mass + particles(j).mass)*v_normal;
                
                fprintf("P1 Name: [" + particles(i).name + "]\n");
                fprintf("P2 Name: [" + particles(j).name + "]\n");
                
                if collision == "uraniumneutron" || collision == "neutronuranium"
                    fprintf("FISSION\n");
                    
                    for k = 1:3
                        neutron = generate_particle("neutron", 5, 1.008, particles(j).pos+10*k, [randi([0 2000]) randi([0 2000]) randi([0 2000])]);
                        particles = [particles neutron]; % Added neutron to list of particles
                    end
                    particles(i) = [];
                    particles(j) = [];
                end

                fprintf("P1 Velocity: [" + num2str(particles(i).velocity) + "]\n");
                fprintf("P2 Velocity: [" + num2str(particles(j).velocity) + "]\n");
        
            end

        end
    end

    % Check for Collisions w/ Walls

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
      

    % Calculate new sphere positions

    for i = 1:size(particles, 2)
        particles(i).pos = particles(i).pos + dt*particles(i).velocity;
    end

    % Figure Animation
    if mod(frame_counter, 5) == 1

       for i = 1:size(particles, 2)
         neutron_final = surf(particles(i).px*particles(i).radius+particles(i).pos(1), particles(i).py*particles(i).radius+particles(i).pos(2), particles(i).pz*particles(i).radius+particles(i).pos(3), FaceColor= "interp", EdgeColor="none");
         hold on
       end

        hold off

        title(sprintf("t = %.1f", t))
        %view(-45, 45);
        
        view(frame_counter/10,45)
        %view(90,90)

        xlim([xboundary.min xboundary.max])
        ylim([yboundary.min yboundary.max])
        zlim([zboundary.min zboundary.max])

        %xlim([neutron.pos(1)-200 neutron.pos(1)+200])
        %ylim([neutron.pos(2)-200 neutron.pos(2)+200])
        %zlim([neutron.pos(3)-200 neutron.pos(3)+200])

        drawnow
    end
   
    frame_counter = frame_counter + 1;
end


