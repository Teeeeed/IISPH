/*******************************************************************************
 *
 * Description: Adapter class between BasicSPH and HoudiniFileDumpHelper
 *
 ******************************************************************************/

#ifndef SPHPARTICLEHOUDINIIO_H
#define SPHPARTICLEHOUDINIIO_H

#include <vector>
#include <string>

class SPHParticle;

namespace SPHParticleHoudiniIO
{
    void outputParticles(const std::vector<SPHParticle>&	particles,
                         const std::string&					outputPath,
                         const std::string&					prefix,
                         int								frameNo);
}

#endif // SPHPARTICLEHOUDINIIO_H
