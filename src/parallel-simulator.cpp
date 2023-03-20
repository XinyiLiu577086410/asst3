#include "quad-tree.h"
#include "world.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <omp.h>

// TASK 2

// NOTE: You may modify this class definition as you see fit, as long as the
// class name, and type of simulateStep and buildAccelerationStructure remain
// the same. You may modify any code outside this class unless otherwise
// specified.

const int QuadTreeLeafSize = 8;
class ParallelNBodySimulator : public INBodySimulator {
public:
  // TODO: implement a function that builds and returns a quadtree containing
  // particles. You do not have to preserve this function type.
  std::unique_ptr<QuadTreeNode> buildQuadTree(std::vector<Particle> &particles,
                                              Vec2 bmin, Vec2 bmax) {
    // TODO: implement a function that builds and returns a quadtree containing
    // particles.
    // BEGIN EDIT
    std::unique_ptr<QuadTreeNode> node(new QuadTreeNode);
    node->particles.clear();
    for(int i = 0; i < 4; i++) node->children[i] = nullptr;
    if(particles.size()<=QuadTreeLeafSize){
      node->isLeaf = true;
      node->particles.swap(particles);
    }
    else{
      node->isLeaf = false;
      std::vector<Particle> Subparticles[4];
      Vec2 bmid = (bmin + bmax) * 0.5f;
      Vec2 borders[4][2] = {
        {{bmin.x,bmin.y}, {bmid.x,bmid.y}},//0
        {{bmid.x,bmin.y}, {bmax.x,bmid.y}},//1
        {{bmin.x,bmid.y}, {bmid.x,bmax.y}},//2
        {{bmid.x,bmid.y}, {bmax.x,bmax.y}},//3
        };
      for(int j = 0;j < 4;j++){
        for(size_t i = 0; i < particles.size();){
          auto p = particles[i];
          if(p.position.x >= borders[j][0].x && p.position.x <= borders[j][1].x
            && p.position.y >= borders[j][0].y && p.position.y <= borders[j][1].y) {
            Subparticles[j].push_back(p);
            if(particles[i].id != particles.back().id) {
              std::swap(particles[i], particles.back());
              particles.pop_back();
              continue;
            }
            else{
              particles.pop_back();
              continue;
            }
          }
          i++;
        }
      }
      for(int i = 0; i < 4; i++){
        node->children[i] = buildQuadTree(Subparticles[i], borders[i][0], borders[i][1]);
      }
    }
    // END EDIT
    return node;
  }
  
  // Do not modify this function type.
  virtual std::unique_ptr<AccelerationStructure>
  buildAccelerationStructure(std::vector<Particle> &particles) {
    // build quad-tree
    auto quadTree = std::make_unique<QuadTree>();

    // find bounds
    Vec2 bmin(1e30f, 1e30f);
    Vec2 bmax(-1e30f, -1e30f);

    for (auto &p : particles) {
      bmin.x = fminf(bmin.x, p.position.x);
      bmin.y = fminf(bmin.y, p.position.y);
      bmax.x = fmaxf(bmax.x, p.position.x);
      bmax.y = fmaxf(bmax.y, p.position.y);
    }

    quadTree->bmin = bmin;
    quadTree->bmax = bmax;

    // build nodes
    std::vector<Particle> tmp (particles);
    quadTree->root = buildQuadTree(particles, bmin, bmax);
    particles.swap(tmp);
    
    if (!quadTree->checkTree()) {
      std::cout << "Your Tree has Error!" << std::endl;
    }

    return quadTree;
  }

  // Do not modify this function type.
  virtual void simulateStep(AccelerationStructure *accel,
                            std::vector<Particle> &particles,
                            std::vector<Particle> &newParticles,
                            StepParameters params) override {
    // TODO: implement parallel version of quad-tree accelerated n-body
    // BEGIN EDIT
    #pragma opm parallel for 
    for (int i = 0; i < (int)particles.size(); i++) {
      Vec2 force = Vec2(0.0f,0.0f);
      std::vector<Particle> attractors;
      auto pi = particles[i];
      accel->getParticles(attractors, pi.position, params.cullRadius);
      for (size_t j = 0; j < attractors.size(); j++) {
          if(pi.id == attractors[j].id) continue;
          force += computeForce(pi, attractors[j], params.cullRadius);
      }
      newParticles[i] = updateParticle(pi, force, params.deltaTime);
    }
    // END EDIT
    // simulation here, using quadTree as acceleration structure
  }
};

// Do not modify this function type.
std::unique_ptr<INBodySimulator> createParallelNBodySimulator() {
  return std::make_unique<ParallelNBodySimulator>();
}
