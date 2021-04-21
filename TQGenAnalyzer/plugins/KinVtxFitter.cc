#include "KinVtxFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h" // MIGHT be useful for Phi->KK?

KinVtxFitter::KinVtxFitter(const std::vector<reco::TransientTrack> tracks, 
                           const std::vector<double> masses, 
                           std::vector<float> sigmas):
  n_particles_{masses.size()} {
  
  KinematicParticleFactoryFromTransientTrack factory;
  std::vector<RefCountedKinematicParticle> particles;
  std::cout<<"size: "<<tracks.size()<<std::endl;
  std::cout<<"======================="<<std::endl;
  for(size_t i = 0; i < tracks.size(); ++i) {
    std::cout<<"m: "<<masses.at(i)<<" chi2: "<<kin_chi2_<<" nof: "<<kin_ndof_<<" sigma: "<<sigmas[i]<<std::endl;
    particles.emplace_back(
      factory.particle(
        tracks.at(i), masses.at(i), kin_chi2_, 
        kin_ndof_, sigmas[i]
        )
      );
  }

  KinematicConstrainedVertexFitter kcv_fitter;    
  RefCountedKinematicTree vtx_tree = kcv_fitter.fit(particles);
  std::cout<<" empty? "<<vtx_tree->isEmpty()<<std::endl;
  std::cout<<" valid? "<<vtx_tree->isValid()<<std::endl;
  //  std::cout<<" consistent? "<<vtx_tree->isConsistent()<<std::endl;
  if (vtx_tree->isEmpty() || !vtx_tree->isValid() || !vtx_tree->isConsistent()) {
    success_ = false; 
    return;
  }

  vtx_tree->movePointerToTheTop(); 
  fitted_particle_ = vtx_tree->currentParticle();
  fitted_vtx_ = vtx_tree->currentDecayVertex();
  if (!fitted_particle_->currentState().isValid() || !fitted_vtx_->vertexIsValid()){ 
    success_ = false; 
    return;
  }
  fitted_state_ = fitted_particle_->currentState();
  fitted_children_ = vtx_tree->finalStateParticles();
  if(fitted_children_.size() != n_particles_) { 
    success_=false; 
    return;
  }
  fitted_track_ = fitted_particle_->refittedTransientTrack();
  success_ = true;
}
