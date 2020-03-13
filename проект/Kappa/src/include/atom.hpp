//
//  atom.hpp
//
// The Atom class loads and stores spectroscopic data for atoms


#ifndef atom_hpp
#define atom_hpp

#include <armadillo>
#include <map>
#include <string>
#include <vector>
#include "constants.h"
#include "particle.hpp"

namespace kappa {

    class Atom : public Particle {
    public:
        // <constructor_docs>
        Atom(const std::string &name, const std::string &filename = "particles.yaml");
        /*
        <cdoc>
            <doc lang="rus">
                Создает объект типа Atom, загружая данные об атоме из базы данных
                
                **Параметры**:
                
                    * ``const std::string &name`` - название атома

                    * ``const std::string &filename`` - путь к файлу базы данных
            </doc>
        </cdoc>
        */
        // </constructor_docs>
    };
}

#endif /* atom_hpp */
