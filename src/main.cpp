/*******************************************************************************
 *
 * Description: Main class that uses the IISPH solver to simulate a fluid.
 * 				The result is output in .geo format after every frame.
 *
 ******************************************************************************/

#include "IISPH.h"
#include "SPHParticleHoudiniIO.cpp"

#include <string>
#include <math.h>
#include <iostream>

//--------------------------------------------------------------------------------------------------
// initParticlesBox
void initParticlesBox(const Vec3d& 				boxMin,
                      const Vec3d& 				boxMax,
                      double 					spacing,
                      std::vector<SPHParticle>&	particles)
{
    int resX = static_cast<int>(floor((boxMax.x-boxMin.x) / spacing));
    int resY = static_cast<int>(floor((boxMax.y-boxMin.y) / spacing));
    int resZ = static_cast<int>(floor((boxMax.z-boxMin.z) / spacing));

    // Clear previous particles
    particles.clear();

    // Fill box with particles
    particles.reserve(resX*resY*resZ);
    for (int x=0; x<resX; ++x)
    {
        for (int y=0; y<resY; ++y)
        {
            for (int z=0; z<resZ; ++z)
            {
                SPHParticle particle;

                // Init position
                particle.pos.x = static_cast<double>(x) * spacing + boxMin.x;
                particle.pos.y = static_cast<double>(y) * spacing + boxMin.y;
                particle.pos.z = static_cast<double>(z) * spacing + boxMin.z;

                // Init other properties
                particle.vel = Vec3d(0,0,0);
				particle.pressure = 0.0;
				particle.sigma = 0.0;
				particle.Psi = 0.0;

                // Add particle
                particles.push_back(particle);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------
//init a static Cubic
void initCubic(const Vec3d& 					boxMin,
                      const Vec3d& 				boxMax,
                      double 					spacing,
                      std::vector<SPHParticle>&	particles)
{
    int resX = static_cast<int>(floor((boxMax.x-boxMin.x) / spacing));
    int resY = static_cast<int>(floor((boxMax.y-boxMin.y) / spacing));
    int resZ = static_cast<int>(floor((boxMax.z-boxMin.z) / spacing));

    // Clear previous particles
    particles.clear();
	particles.reserve(resX*resY*resZ);
    for (int x=0; x<=resX; ++x)	
    {
		for (int y=0; y<resY; ++y)
        {
            for (int z=0; z<resZ; ++z)
            {
                SPHParticle particle;

                // Init position
                particle.pos.x = static_cast<double>(x) * spacing + boxMin.x;
                particle.pos.y = static_cast<double>(y) * spacing + boxMin.y;
                particle.pos.z = static_cast<double>(z) * spacing + boxMin.z;

                // Init other properties
                particle.vel = Vec3d(0,0,0);
				particle.accel = Vec3d(0,0,0);
				particle.veladv = Vec3d(0,0,0);
				particle.density = 0.0;
				particle.oneOverDensity = 0.0;
				particle.pressure = 0.0;
				particle.sigma = 0.0;
				particle.Psi = 0.0;
                // Add particle
                particles.push_back(particle);
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------
// initParticlesBox
void initBoundary(const Vec3d&						boxMin,
				  const Vec3d&						boxMax,
				  double							spacing,
				  std::vector<SPHParticle>&			particles)
{
    int resX = static_cast<int>(floor((boxMax.x-boxMin.x) / spacing));
    int resY = static_cast<int>(floor((boxMax.y-boxMin.y) / spacing));
    int resZ = static_cast<int>(floor((boxMax.z-boxMin.z) / spacing));

    // Clear previous particles
    particles.clear();

    // Fill the surface of rigid body with particles
    particles.reserve(resY*resZ*2+resX*resY*2+resX*resZ*2);
    for (int x=0; x<resX; x += resX-1)
    {
        for (int y=0; y<resY; ++y)
        {
            for (int z=0; z<resZ; ++z)
            {
                SPHParticle particle;

                // Init position
                particle.pos.x = static_cast<double>(x) * spacing + boxMin.x;
                particle.pos.y = static_cast<double>(y) * spacing + boxMin.y;
                particle.pos.z = static_cast<double>(z) * spacing + boxMin.z;

                 //Init other properties
                particle.vel = Vec3d(0,0,0);
				particle.accel = Vec3d(0,0,0);
				particle.density = 0.0;
				particle.oneOverDensity = 0.0;
				particle.pressure = 0.0;
				particle.sigma = 0.0;
				particle.Psi = 0.0;
				particle.veladv = Vec3d(0,0,0);
                // Add particle
                particles.push_back(particle);
            }
        }
    }
	for (int z=0; z<resZ; z += resZ-1)
    {
        for (int y=0; y<resY; ++y)
        {
            for (int x=1; x<resX-1; ++x)
            {
                SPHParticle particle;

                // Init position
                particle.pos.x = static_cast<double>(x) * spacing + boxMin.x;
                particle.pos.y = static_cast<double>(y) * spacing + boxMin.y;
                particle.pos.z = static_cast<double>(z) * spacing + boxMin.z;

                // Init other properties
                particle.vel = Vec3d(0,0,0);
				particle.accel = Vec3d(0,0,0);
				particle.density = 0.0;
				particle.oneOverDensity = 0.0;
				particle.pressure = 0.0;
				particle.sigma = 0.0;
				particle.Psi = 0.0;
				particle.veladv = Vec3d(0,0,0);
                // Add particle
                particles.push_back(particle);
            }
        }
    }
    for (int y=0; y<resY; y += resY-1)	
    {
		for (int z=1; z<resZ-1; ++z)
        {
            for (int x=1; x<resX-1; ++x)
            {
                SPHParticle particle;

                // Init position
                particle.pos.x = static_cast<double>(x) * spacing + boxMin.x;
                particle.pos.y = static_cast<double>(y) * spacing + boxMin.y;
                particle.pos.z = static_cast<double>(z) * spacing + boxMin.z;

                // Init other properties
                particle.vel = Vec3d(0,0,0);
				particle.accel = Vec3d(0,0,0);
				particle.density = 0.0;
				particle.oneOverDensity = 0.0;
				particle.pressure = 0.0;
				particle.sigma = 0.0;
				particle.Psi = 0.0;
				particle.veladv = Vec3d(0,0,0);
                // Add particle
                particles.push_back(particle);
            }
        }
    }
}
//--------------------------------------------------------------------------------------------------
// MAIN
int main(int argc, char* argv[])
{
    // Validate arguments
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " OutputDirectory" << std::endl;
        return -1;
    }

    std::string outputPath(argv[1]);


    // Init simulator
    std::cout << "Initializing simulator..." << std::endl;
    //const Vec3d volumeMin(-5,  -1, -2);
    //const Vec3d volumeMax( 5,  5,  2);
    const Vec3d volumeMin(-9, -3, -6);
    const Vec3d volumeMax( 9, 6,  6);
    const double mass = 3.1;
    const double restDensity = 998.23;
    const double h = 0.3;
    const double k = 100.0;
    const double dt = 0.011;
    IISPH sph(volumeMin, volumeMax, mass, restDensity, h, k, dt);
    //sph.enableAdaptiveTimeStep();

    // Init particles
    std::cout << "Initializing particles" << std::endl;
    Vec3d boxMin(-9, 0, -2.25);
    Vec3d boxMax(-3, 4.5, 2.25);
    initParticlesBox(boxMin, boxMax, h/2.0, sph.particles());

	// Init boundary particles
	std::cout << "Initializing boundary particles" << std::endl;
	std::vector<SPHParticle> cubicparticles;
	cubicparticles.clear();
    Vec3d cubicMin( 2.4, 0, -0.75);
    Vec3d cubicMax( 3.9, 3.0,  0.75);
    initCubic(cubicMin, cubicMax, h/2.0, cubicparticles);
	std::vector<SPHParticle> faceparticles1;
    Vec3d boundaryMin1( -9.15, -0.15, -2.4);
    Vec3d boundaryMax1( 9.15, 6.15, 2.4);
    initBoundary(boundaryMin1, boundaryMax1, h/2.0, faceparticles1);
    sph.boundaryparticles().insert(sph.boundaryparticles().begin(),faceparticles1.begin(),faceparticles1.end());
    sph.boundaryparticles().insert(sph.boundaryparticles().begin(),cubicparticles.begin(),cubicparticles.end());
	sph.searchBoundaryNeighbors();
	sph.computePsi();

    // Output first frame (initial frames)
    const std::string filePrefix("particles_");
	// Merge two vectors
	std::vector<SPHParticle> allparticles;
	allparticles.clear();
    //allparticles.insert(allparticles.begin(), sph.boundaryparticles().begin(), sph.boundaryparticles().end());
    allparticles.insert(allparticles.begin(),cubicparticles.begin(),cubicparticles.end());
	allparticles.insert(allparticles.begin(),sph.particles().begin(),sph.particles().end());
    //
    SPHParticleHoudiniIO::outputParticles(allparticles, outputPath, filePrefix, 1);

    // Run simulation and output a frame every 1/24 second
    std::cout << "Running simulation!" << std::endl;
    const double frameTime = 1.0/24.0;
    const double totalSimulationTime = 10.0;
    double time = 0.0;
    int currentFrame = 2;
    while (time < totalSimulationTime)
    {
        std::cout << "Simulating frame " << currentFrame << " (" << (frameTime+time) << "s)";
        std::cout << std::endl;
        // Run simulation
        sph.run(frameTime);
        // Output particles
		// Merge two vectors
		std::vector<SPHParticle> allparticles;
		allparticles.clear();
        allparticles.insert(allparticles.begin(),cubicparticles.begin(),cubicparticles.end());
		allparticles.insert(allparticles.begin(),sph.particles().begin(),sph.particles().end());
		//
        SPHParticleHoudiniIO::outputParticles(allparticles, outputPath, filePrefix, currentFrame);
        // Update simulation time
        time += frameTime;
        ++currentFrame;
    }

   std::cout << "Done!" << std::endl;
   return 0;
}
