//includes
#include "CalibTreeRecord.h"
namespace CalibTreeRecord {

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 
  void CalibTreeRecord::Clear() 
  {
    eves.clear();
    hits.clear();
    ophits.clear();
    run=0;
    subrun=0;
    event=0;
  }

  void CalibTreeRecord::stdout_dump() 
  {
    std::cout<<"N eves is "<<eves.size()<<"\n";
    for(auto eve : eves){
      std::cout<<"N Particles is : "<<eve.particles.size()<<"\n";
      std::cout<<"Eve position is X:"<<eve.x_pos<<" Y:"<<eve.y_pos<<" Z:"<<eve.z_pos<<" T:"<<eve.t_pos<<"\n";
      for(auto particle: eve.particles){
        std::cout<<"Particle position is X:"<<particle.x_pos<<" Y:"<<particle.y_pos<<" Z:"<<particle.z_pos<<" t:"<<particle.t_pos<<"\n";
        for(auto hit : particle.partial_hits){
          std::cout<<"Partial Hit. Charge:"<<hit.charge<<" ne:"<<hit.num_electrons<<" en:"<<hit.energy<<" t:"<<hit.width<<"\n";
        }
        for(auto hit :particle.partial_ophits){
          std::cout<<"ParticlOphit. pes: "<<hit.pes<<" np:"<<hit.num_photons<<" t:"<<hit.time<<"\n";
        }
      }
    }
  }
}

