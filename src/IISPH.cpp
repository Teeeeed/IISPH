/*******************************************************************************
 *
 * Description:	IISPH fluid solver based on the paper "Implicit Incompressible
 * 				SPH" (Ihmsen and Cornelis 2013). Surface tension forces are 
 *				not handled in this version, but they will be in the 
 *				future. We use a M4 Spline kernel in our calculations.
 * 				The particles are integrated using the semi-implicit Euler
 * 				method and we handle boundaries by using boundary particles 
 *				with repulsive forces. 
 *
 ******************************************************************************/

#include "IISPH.h"
#include "SpatialGrid.h"

#include <iostream>

namespace
{
    // This is much faster than calling pow(val, exponent)
    inline double pow2(double val) { return val*val; }
    inline double pow3(double val) { return val*val*val; }
    inline double pow7(double val) { return val*val*val*val*val*val*val; }

    inline void checknan(const double& v, int id, const char* name, bool dothrow = true)
    {
        if (v!=v)
        {
            std::cout << name << "[" << id << "] = " << v << std::endl;
            if (dothrow)
            {
                throw(std::string("NAN"));
            }
        }
    }

    inline void checknan(const Vec3d& v, int id, const char* name, bool dothrow = true)
    {
        checknan(v.x, id, name, dothrow);
        checknan(v.y, id, name, dothrow);
        checknan(v.z, id, name, dothrow);
    }

    inline void checkpositive(const double& v, int id, const char* name)
    {
        if (v<0.0)
        {
            std::cout << name << "[" << id << "] = " << v << std::endl;
            throw(std::string("NEGATIVE"));
        }
    }
}

//--------------------------------------------------------------------------------------------------
// Constructor / Destructor
//--------------------------------------------------------------------------------------------------
IISPH::IISPH(Vec3d volumeMin, Vec3d volumeMax, double mass, double restDensity, double h,
                   double k, double dt)
    : _volumeMin(volumeMin),
      _volumeMax(volumeMax),
      _mass(mass),
      _restDensity(restDensity),
      _k(k),
      _dt(dt),
      _h(h),
      _bulkViscosity(0.5),
      _shearViscosity(0.0),
      _useAdaptiveTimeStep(false),
      _maxuij(0.0)
{
    precomputeKernelCoefficients();
}

IISPH::~IISPH()
{}

//--------------------------------------------------------------------------------------------------
// Public functions
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
// enableAdaptiveTimeStep
void IISPH::enableAdaptiveTimeStep(double tolerance, double min, double max)
{
    _useAdaptiveTimeStep = true;
    _tolerance = tolerance;
    _minDT = min;
    _maxDT = max;
}

//--------------------------------------------------------------------------------------------------
// run
void IISPH::run(double time)
{
    double oldTimeStep = _dt;
    // Run simulation!
    double timeLeft = time;
    while (timeLeft > 0.0)
    {
        // Run simulation steps

        //computeDensityAndPressure();
        //addExternalForces();
        //computeArtificialViscosityForces();
        //computePressureForces();
         //Compute time step
        //if (_useAdaptiveTimeStep)
        //    computeTimeStep();

        // Limit timestep to the time left in the simulation
        if (timeLeft < _dt)
        {
            _dt = timeLeft;
        }
		//searchBoundaryNeighbors();//build boundary particles' boundary neighbors
		//computePsi();
        buildFluidNeighbors();//build fluid particles' fluid neighbors
		buildBoundaryNeighbors();//build fluid particles' boundary neighbors
		predictadvection();
        // Update particles
        integrate();

        // Update time
        timeLeft -= _dt;

        std::cout << "Substep done! With timestep = " << _dt << std::endl;
    }
    // Restore old time step
    _dt = oldTimeStep;
}

//--------------------------------------------------------------------------------------------------
// Simulation stepsrhoij
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
// searchneighbors
void IISPH::searchBoundaryNeighbors()
{
	// Reserve space and initialize neighbors' data
    _rigidneighbors.clear();
    _rigidneighbors.resize(_boundaryparticles.size());

    // Init spatial grid
	/*int h = 0.02;*/
    Vec3d borders(_h, _h, _h);
	Vec3d gridMin = _volumeMin;
    gridMin -= borders;
    Vec3d gridMax = _volumeMax;
    gridMax += borders;
    SpatialGrid<long> grid(_h, gridMin, gridMax);

    // Insert particles into grid
    for (long p=0; p<_boundaryparticles.size(); ++p)
    {
        grid.insert(p, _boundaryparticles[p].pos);
    }

    // Use grid to retrieve neighbor particles
    double h2 = _h*_h;
    std::vector<long*> nearbyParticles;
    for (long p=0; p<_boundaryparticles.size(); ++p)
    {
        const SPHParticle &particle = _boundaryparticles[p];
		nearbyParticles.clear();
        // Get nearby particles
        grid.getElements(particle.pos, _h, nearbyParticles);

        // Find particles that are within smoothing radius
        _rigidneighbors[p].reserve(100);
        for (long i=0; i<nearbyParticles.size(); ++i)
        {
            long nID = *nearbyParticles[i];
            const SPHParticle &neighborParticle = _boundaryparticles[nID];

            // Skip current particle
            if (nID==p)
                continue;

            Vec3d xij = particle.pos - neighborParticle.pos;

            // Check if distance is lower than smoothing radius
            double dist2 = xij.dot(xij);
            if (dist2 < h2)
            {
                // Yup! Add the particle to the neighbors list along with
                // some precomputed informations
                _rigidneighbors[p].push_back(Neighbor(nID, xij, sqrt(dist2)));
            }
        }	
    }
}

//--------------------------------------------------------------------------------------------------
// computePsi
void IISPH::computePsi()
{
    const double Wii = getKernelValue(0.0);

    // Iterate particles
    for (long p=0; p<_boundaryparticles.size(); ++p)
    {
        SPHParticle &particle = _boundaryparticles[p];

        // Reinitialize particle properties
		particle.vel = Vec3d(0,0,0);
        particle.density = 0.0;
        particle.oneOverDensity = 0.0;
        particle.accel = Vec3d(0,0,0);
        particle.pressure = 0.0;
		particle.sigma = 0.0;
		particle.Psi = 0.0;

        // Compute Psi from neighbors contributions
        // sigma_i = SUM_j ( Wij)
        for (long n=0; n<_rigidneighbors[p].size(); ++n)
        {
            const Neighbor &neighbor = _rigidneighbors[p][n];

            // Add contribution
			// Compute Boundary particles' sigma
            particle.sigma += getKernelValue(neighbor.dist);
        }

        // Add current particle's contribution
        particle.sigma += Wii;

        // Compute Boundary particles' Psi
        if (particle.sigma != 0.0)
        {
            particle.Psi= 1.0/particle.sigma;
            particle.Psi *= _restDensity;
        }
    }
}

//--------------------------------------------------------------------------------------------------
// buildFluidNeighbors
void IISPH::buildFluidNeighbors()
{
	int max=0;
	int max1=0;
    // Reserve space and initialize neighbors' data
    _neighbors.clear();
    _neighbors.resize(_particles.size());

    // Init spatial grid
    Vec3d borders(_h, _h, _h);
    Vec3d gridMin = _volumeMin;
    gridMin -= borders;
    Vec3d gridMax = _volumeMax;
    gridMax += borders;
    SpatialGrid<long> grid(_h, gridMin, gridMax);

    // Insert particles into grid
    for (long p=0; p<_particles.size(); ++p)
    {
        grid.insert(p, _particles[p].pos);
    }

    // Use grid to retrieve neighbor particles
    double h2 = _h*_h;
    std::vector<long*> nearbyParticles;
    for (long p=0; p<_particles.size(); ++p)
    {
        const SPHParticle &particle = _particles[p];

        // Get nearby particles
		nearbyParticles.clear();
        grid.getElements(particle.pos, _h, nearbyParticles);
		if(nearbyParticles.size()>max1){max1=nearbyParticles.size();}
        // Find particles that are within smoothing radius
        _neighbors[p].clear();
        for (long i=0; i<nearbyParticles.size(); ++i)
        {
            long nID = *nearbyParticles[i];
            const SPHParticle &neighborParticle = _particles[nID];

            // Skip current particle
            if (nID==p)
                continue;

            Vec3d xij = particle.pos - neighborParticle.pos;

            // Check if distance is lower than smoothing radius
            double dist2 = xij.dot(xij);
            if (dist2 < h2)
            {
				// Yup! Add the particle to the neighbors list along with
				// some precomputed informations
				_neighbors[p].push_back(Neighbor(nID, xij, sqrt(dist2)));
            }
        }
		if(_neighbors[p].size()>max){max=_neighbors[p].size();}
    }
	std::cout<<max<<std::endl;
	std::cout<<max1<<std::endl;
}

//--------------------------------------------------------------------------------------------------
// buildBoundaryNeighbors
void IISPH::buildBoundaryNeighbors()
{
	int max=0;
	int max1=0;
    // Reserve space and initialize neighbors' data
    _boundaryneighbors.clear();
    _boundaryneighbors.resize(_particles.size());

    // Init spatial grid
    Vec3d borders(_h, _h, _h);
    Vec3d gridMin = _volumeMin;
    gridMin -= borders;
    Vec3d gridMax = _volumeMax;
    gridMax += borders;
    SpatialGrid<long> grid(_h, gridMin, gridMax);

    // Insert particles into grid
    for (long p=0; p<_boundaryparticles.size(); ++p)
    {
        grid.insert(p, _boundaryparticles[p].pos);
    }

    // Use grid to retrieve neighbor particles
    double h2 = _h*_h;
    std::vector<long*> nearbyParticles;
    for (long p=0; p<_particles.size(); ++p)
    {
        const SPHParticle &particle = _particles[p];
		_boundaryneighbors[p].clear();
		nearbyParticles.clear();
        // Get nearby particles
        grid.getElements(particle.pos, _h, nearbyParticles);
		if(nearbyParticles.size()>max1){max1=nearbyParticles.size();}
        // Find particles that are within smoothing radius
        //_boundaryneighbors[p].reserve(100);
        for (long i=0; i<nearbyParticles.size(); ++i)
        {
            long nID = *nearbyParticles[i];
            const SPHParticle &neighborParticle = _boundaryparticles[nID];

            // Skip current particle
            //if (nID==p)
            //    continue;

            Vec3d xij = particle.pos - neighborParticle.pos;

            // Check if distance is lower than smoothing radius
            double dist2 = xij.dot(xij);
            if (dist2 < h2)
            {
                // Yup! Add the particle to the neighbors list along with
                // some precomputed informations
                _boundaryneighbors[p].push_back(Neighbor(nID, xij, sqrt(dist2)));
            }
        }
		if(_boundaryneighbors[p].size()>max){max=_boundaryneighbors[p].size();}
    }
	std::cout<<max<<std::endl;
	std::cout<<max1<<std::endl;
}

//IISPH
void IISPH::viadv(SPHParticle &particle,int i)
{

    const double alpha = 1.0;	// Bulk viscosity
    const double beta = 0.0;	// Shear viscosity

    // Precompute coefficients
    const double speedOfSound = sqrt(_k);
    const double h2 = _h*_h;

	// Reset acceleration
    particle.accel = Vec3d(0,0,0);
    // Get neighbors contributions
#if 1
#if 1
    for (long n=0; n<_neighbors[i].size(); ++n)
    {
        const Neighbor &neighbor = _neighbors[i][n];
        const SPHParticle &neighborParticle = _particles[neighbor.id];
        // Compute contribution
        Vec3d vij = particle.vel - neighborParticle.vel;
        double vijxij = vij.dot(neighbor.xij);
        double dij = neighbor.dist;
        double uij = _h*vijxij / (dij*dij + 0.01*h2);
        if (uij < 0)
        {
            // Compute contribution
            double avgDensity = 0.5 * (particle.density + neighborParticle.density);

            double IIij = alpha*uij*speedOfSound / avgDensity;
            Vec3d contribution = getKernelGradient(neighbor.dist, neighbor.xij);
            contribution *= IIij;
            contribution *= _mass;
            particle.accel += contribution;
        }
    }

	//viscosity force between fluid particles and boundary particles
#if 1
    for (long n=0; n<_boundaryneighbors[i].size(); ++n)
    {
		const Neighbor &boundaryneighbor = _boundaryneighbors[i][n];
        const SPHParticle &boundaryneighborparticle = _boundaryparticles[boundaryneighbor.id];
        Vec3d vij = particle.vel - boundaryneighborparticle.vel;
        double vijxij = vij.dot(boundaryneighbor.xij);
        double dij = boundaryneighbor.dist;
        double uij = _h*vijxij / (dij*dij + 0.01*h2);
        if (uij < 0)
        {
            // Compute contribution
            double avgDensity = 2.0 * particle.density;
            double IIij = (2.0*uij*speedOfSound) / avgDensity;
            Vec3d contribution = getKernelGradient(boundaryneighbor.dist, boundaryneighbor.xij);
            contribution *= IIij;
            contribution *= boundaryneighborparticle.Psi;

            particle.accel += contribution;
        }
    }
#endif
#endif
#if 1
    particle.accel += Vec3d(0.0, -9.81, 0.0);
#endif

#endif
    particle.veladv = particle.vel + (particle.accel * _dt);
}

//predict advection
void IISPH::predictadvection()
{
    double trho=0,finalterm=0,trho1=0;
    double t=0,t1=0,t2=0;
    Vec3d gradwij=Vec3d(0,0,0);
    Vec3d temp=Vec3d(0,0,0);
	Vec3d temp1=Vec3d(0,0,0);
    Vec3d term=Vec3d(0,0,0);
    Vec3d term2=Vec3d(0,0,0);
    Vec3d term3=Vec3d(0,0,0);
    Vec3d dji=Vec3d(0,0,0);

    std::vector<double> rhoadv;
    std::vector<double> prevtpi;
    std::vector<double> tpi;
    std::vector<double> aii;
    std::vector<Vec3d> dii;
    std::vector<Vec3d> term1;

    rhoadv.resize(_particles.size());
    prevtpi.resize(_particles.size());
    tpi.resize(_particles.size());
    aii.resize(_particles.size());
    dii.resize(_particles.size());
    term1.resize(_particles.size());

    const double Wii = _mass*getKernelValue(0.0);

    for(int i=0;i<_particles.size();++i)
    {
        trho = 0;
		trho1 = 0;
        SPHParticle &p=_particles[i];
        for(int j=0;j<_neighbors[i].size();++j)
        {
            Neighbor &n=_neighbors[i][j];
            trho += _mass * getKernelValue(n.dist);
        }

        trho += Wii;

		//boundary part
#if 1
        for(int j=0;j<_boundaryneighbors[i].size();++j)
		{
			Neighbor &n=_boundaryneighbors[i][j];
			SPHParticle &b=_boundaryparticles[n.id];
			trho1 += b.Psi * getKernelValue(n.dist);
        }
#endif
        if(trho != 0 || trho1 != 0)
        {
            p.density=trho+trho1;
        }
        temp=Vec3d(0,0,0);
        for(int j=0;j<_neighbors[i].size();++j)
        {
            gradwij=Vec3d(0,0,0);
            Neighbor &n=_neighbors[i][j];
            gradwij = getKernelGradient(n.dist,n.xij);
            t = -_mass/(p.density*p.density);
            temp += gradwij * t;
        }
		//boundary part
#if 1
        temp1=Vec3d(0,0,0);
        for(int j=0;j<_boundaryneighbors[i].size();++j)
        {
            gradwij=Vec3d(0,0,0);
            Neighbor &n=_boundaryneighbors[i][j];
            SPHParticle &b=_boundaryparticles[n.id];
            gradwij = getKernelGradient(n.dist,n.xij);
            t = -_boundaryparticles[n.id].Psi/(p.density*p.density);
            temp1 += gradwij * t;
        }
#endif
        dii[i] = (temp + temp1) * _dt*_dt;
        viadv(p,i);
    }
    for(int i=0;i<_particles.size();++i)
    {
        t=0;
		t1=0;
        t2=0;
        SPHParticle &p=_particles[i];

        for(int j=0;j<_neighbors[i].size();++j)
        {
            gradwij=Vec3d(0,0,0);
            temp=Vec3d(0,0,0);
            dji=Vec3d(0,0,0);
            Neighbor &n=_neighbors[i][j];
			SPHParticle &particle=_particles[n.id];
			//compute rhoadv
            gradwij = getKernelGradient(n.dist,n.xij);
            temp=p.veladv - _particles[n.id].veladv;
            t += gradwij.dot(temp)*_mass*_dt;

			//compute aii
            gradwij =  getKernelGradient(n.dist,n.xij);
            dji = gradwij * (_mass*_dt*_dt/(p.density*p.density));
            temp=dii[i] - dji;
            t2 += (temp.dot(gradwij)*_mass);

        }
		//boundary part
#if 1
        for(int j=0;j<_boundaryneighbors[i].size();++j)
        {
            gradwij=Vec3d(0,0,0);
            temp=Vec3d(0,0,0);
            Neighbor &n=_boundaryneighbors[i][j];
			SPHParticle &b=_boundaryparticles[n.id];

            gradwij = getKernelGradient(n.dist,n.xij);
            temp=p.veladv - b.veladv;
            t += gradwij.dot(temp)*b.Psi*_dt;

            // contribution to aii
            t2 += dii[i].dot(gradwij) * b.Psi;
        }
#endif
        if(t2==0) t2=1;
        aii[i] = t2;
        rhoadv[i] = p.density + t;
        tpi[i] = (0.5*p.pressure);
        prevtpi[i] = (0.5*p.pressure);
    }
#if 1   //iterate pi
    int l=0;
    double averageDensityError = 1e10;
    double eta = 0.001*_restDensity;
    while(((averageDensityError > eta) || l<2) && (l<200))
    //while (l<10)
    {
        double sumDensities = 0.0;
        long nbParticlesInSummation = 0;
        gradwij=Vec3d(0,0,0);
        term=Vec3d(0,0,0);
        term2=Vec3d(0,0,0);
        term3=Vec3d(0,0,0);
        temp=Vec3d(0,0,0);
		finalterm=0;
        for(int i=0;i<_particles.size();i++)
        {
            temp=Vec3d(0,0,0);
            for(int j=0;j<_neighbors[i].size();++j)
            {
                gradwij=Vec3d(0,0,0);
                Neighbor &n=_neighbors[i][j];
                gradwij = getKernelGradient(n.dist,n.xij);
                t=_particles[n.id].density;
                t2 = -(_mass*prevtpi[n.id])/(t*t);
                temp += gradwij * t2;
            }
			//term1.resize(particles.size());
            term1[i] = (temp * (_dt*_dt));
        }
        for(int i=0;i<_particles.size();++i)
        {
            term=Vec3d(0,0,0);
            term3=Vec3d(0,0,0);
            term2=Vec3d(0,0,0);
            finalterm=0;
            for(int j=0;j<_neighbors[i].size();++j)
            {
                Neighbor &n=_neighbors[i][j];

                term2=dii[n.id] * prevtpi[n.id];
                term3=Vec3d(0,0,0);
                temp=Vec3d(0,0,0);

                // Compute dji
                Vec3d dji = getKernelGradient(n.dist, n.xij);
                dji *= _dt*_dt * _mass / (_particles[i].density * _particles[i].density);

                // Compute sum(djk*pk)-dji*pi
                term3 = term1[n.id] - dji*prevtpi[i];


                /*for(int k=0;k<_neighbors[n.id].size();++k)
                {

                    Neighbor &n2=_neighbors[n.id][k];

                    if(n2.id==i) continue;

                    gradwij=getKernelGradient(n2.dist,n2.xij);
                    t=_particles[n2.id].density;
                    t2=-(_mass*prevtpi[n2.id])/(t*t);
                    temp += gradwij * t2;

                }

                term3=temp * (_dt*_dt);*/

                gradwij=getKernelGradient(n.dist,n.xij);

                temp=Vec3d(0,0,0);

                temp = term1[i] - term2 ;

                term = temp - term3;//term1[id];

                finalterm=finalterm+(_mass * term.dot(gradwij));

            }
			//boundary part
#if 1
            for(int j=0;j<_boundaryneighbors[i].size();++j)
            {
				gradwij=Vec3d(0,0,0);
                Neighbor &n=_boundaryneighbors[i][j];
				SPHParticle &b=_boundaryparticles[n.id];
				gradwij=getKernelGradient(n.dist,n.xij);
                finalterm=finalterm + (b.Psi * term1[i].dot(gradwij));
            }
#endif
            t=0.5*prevtpi[i] + (0.5/aii[i]) * (_restDensity- rhoadv[i]- finalterm);
            if (t<0.0)
            {
                t = 0.0;
            }
            else
            {
                // Estimate rho_i
                // Note: rho_i is computed using the previous pressure values, which
                //		 will stop the loop one iteration later. However, this way is
                //		 more efficient, since we use previously computed terms, instead
                //		 of using lengthy computations every iteration.
                sumDensities += rhoadv[i] + aii[i]*prevtpi[i] + finalterm;
                ++nbParticlesInSummation;
            }
            tpi.at(i)=t;

        }			
        prevtpi=tpi;
        term1.clear();
        term1.resize(_particles.size());
        l=l+1;

        // Compute the average density error rho_avg-rho_o (used by the stopping
        // criterion)
        if (nbParticlesInSummation > 0)
        {
            averageDensityError = (sumDensities / nbParticlesInSummation) -_restDensity;
        }
        else
        {
            averageDensityError = 0.0;
        }
    }

    std::cout << "Pressure solver ended after " << l
              << " iterations... (rho_avg-rho_o = "
              << averageDensityError << ")" << std::endl;
#endif

    for(int i=0; i<_particles.size(); i++) _particles[i].density = rhoadv[i];
    for(int i=0; i<_particles.size(); i++) _particles[i].pressure = tpi.at(i);
    tpi.clear();
    prevtpi.clear();
    rhoadv.clear();
    aii.clear();
    dii.clear();
    term1.clear();
}

void IISPH::integrate()
{

    double t=0,t1=0,t2=0,t3=0;
    Vec3d fpi=Vec3d(0,0,0);
	Vec3d bpi=Vec3d(0,0,0);
    Vec3d gradwij=Vec3d(0,0,0);
    // Update particles velocity and position
    for (long p=0; p<_particles.size(); ++p)
    {
		Vec3d contribution=Vec3d(0,0,0);
        SPHParticle &particle = _particles[p];

        t1=particle.pressure/(particle.density*particle.density);
        fpi=Vec3d(0,0,0);
		bpi=Vec3d(0,0,0);
        for(int j=0;j<_neighbors[p].size();j++)
        {
			gradwij=Vec3d(0,0,0);
            Neighbor &n=_neighbors[p][j];
            t2=_particles[n.id].pressure/(_particles[n.id].density*_particles[n.id].density);
            gradwij=getKernelGradient(n.dist,n.xij);
            t=t1+t2;
            fpi += gradwij * t;
        }
		//boundary pressure force
#if 1
        for(int j=0;j<_boundaryneighbors[p].size();++j)
		{
			gradwij=Vec3d(0,0,0);
			Neighbor &n=_boundaryneighbors[p][j];
			gradwij=getKernelGradient(n.dist,n.xij);
			bpi += gradwij*t1*_boundaryparticles[n.id].Psi;
        }
#endif
        particle.vel = particle.veladv + (fpi * (-_mass*_dt)) - bpi*_dt;
        // Update position
        particle.pos += (particle.vel * _dt);
        //particle.vel = Vec3d(0,0,0);	// TODO: REMOVE ME!
        //particle.veladv = Vec3d(0,0,0);	// TODO: ME TOO!

        // Apply boundary conditions

        if (particle.pos.x < _volumeMin.x)
        {
            particle.pos.x = _volumeMin.x;
            particle.vel.x = 0.0;
        }
       else if (particle.pos.x > _volumeMax.x)
        {
            particle.pos.x = _volumeMax.x;
            particle.vel.x = 0.0;
        }
        if (particle.pos.y < _volumeMin.y)
        {
            particle.pos.y = _volumeMin.y;
            particle.vel.y = 0.0;
        }
        else if (particle.pos.y > _volumeMax.y)
        {
            particle.pos.y = _volumeMax.y;
            particle.vel.y = 0.0;
        }
        if (particle.pos.z < _volumeMin.z)
        {
            particle.pos.z = _volumeMin.z;
            particle.vel.z = 0.0;
        }
        else if (particle.pos.z > _volumeMax.z)
        {
            particle.pos.z = _volumeMax.z;
            particle.vel.z = 0.0;
        }

    }

}
//IISPH


//--------------------------------------------------------------------------------------------------
// computeTimeStep
void IISPH::computeTimeStep()
{
    // Find maximum acceleration
    double maxAccel2 = 0.0;
    for (long p=0; p<_particles.size(); ++p)
    {
        const Vec3d &accel = _particles[p].accel;

        // Test squared acceleration length
        double accelLength2 = accel.x*accel.x + accel.y*accel.y + accel.z*accel.z;
        if (accelLength2 > maxAccel2)
        {
            maxAccel2 = accelLength2;
        }
    }

    // Compute force
    double maxForce = _mass * sqrt(maxAccel2);	// f = mass * a

    // Compute timestep (Based on paper "Smoothed Particles Hydrodynamics" (Monaghan 1992))
    // dt = tolerance * min(dta, dtcv)
    // dta = min_i(h/|f_i|)
    // dtcv = min_i( h/(c + 0.6*(alpha * c + beta * max_j(uij))) )
    double alpha = _bulkViscosity;
    double beta = _shearViscosity;
    double speedOfSound = sqrt(_k);	// c
    double tf = sqrt(_h / maxForce);
    double tcv = _h / (speedOfSound + 0.6*(alpha*speedOfSound + beta*_maxuij));
    _dt = (tf < tcv) ? tf : tcv;
    _dt *= _tolerance;

    // Clamp time step
    if (_dt < _minDT) _dt = _minDT;
    if (_dt > _maxDT) _dt = _maxDT;
}

//--------------------------------------------------------------------------------------------------
// Kernel functions
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
// precomputeKernelCoefficients
void IISPH::precomputeKernelCoefficients()
{
    const double PI = 3.14159265359;

    _halfH = _h/2.0;	// In Monaghan2005, h=half of smoothing radius

    // Precompute value coefficient (Identical for part A and B)
    _kernelValueCoeff = 1.0 / (4.0*PI*pow(_halfH,3));

    // Precompute gradient coefficients
    //_kernelGradientCoeffA = 3.0 / (4.0*PI*pow(_halfH,5));
    //_kernelGradientCoeffB = -3.0 / (4.0*PI*pow(_halfH,4));
    _kernelGradientCoeffA = 3.0 / (4.0*PI*pow(_halfH,4));
    _kernelGradientCoeffB = -3.0 / (4.0*PI*pow(_halfH,4));

    // Precompute laplacian coefficients
    _kernelLaplacianCoeffA = -9.0 / (PI*pow(_halfH,5));
    _kernelLaplacianCoeffB = 3.0 / (PI*pow(_halfH,5));

}

//--------------------------------------------------------------------------------------------------
// getKernelValue
inline
double IISPH::getKernelValue(double dist) const
{
    double q = dist/_halfH;
    if (q<1.0)
    {
        return _kernelValueCoeff * ( pow3(2.0-q)-4*pow3(1.0-q) );
    }
    else
    {
        return _kernelValueCoeff * pow3(2.0-q);
    }
}

//--------------------------------------------------------------------------------------------------
// getKernelGradient
inline
Vec3d IISPH::getKernelGradient(double dist, const Vec3d& xij) const
{
    double q = dist/_halfH;
    Vec3d gradient = xij;
    if (q<= 0.0)
    {
        gradient = Vec3d(0,0,0);
    }
    else if (q<1.0)
    {
        gradient *= _kernelGradientCoeffA * (4.0 * pow2(1.0-q) - pow2(2.0-q)) / dist;
    }
    else
    {
        gradient *= (_kernelGradientCoeffB * pow2(2.0 - q)) / dist;
    }

    return gradient;

}

//--------------------------------------------------------------------------------------------------
// getKernelLaplacian
inline
double IISPH::getKernelLaplacian(double dist) const
{
    double q = dist/_halfH;
    if (q<=1.0)
    {
        return _kernelLaplacianCoeffA * (1.0-q);
    }
    else
    {
        return _kernelLaplacianCoeffB * (3.0-q-(2.0/q));
    }
}
