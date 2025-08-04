function particle = generate_particle(name, radius, mass, pos, velocity)
    
    particle.name = name;
    particle.radius = radius;
    particle.mass = mass;
    particle.pos = pos;
    particle.velocity = velocity;

    num_faces = 10;

    [x, y, z] = sphere(num_faces);
    particle.px = x;
    particle.py = y;
    particle.pz = z;


end