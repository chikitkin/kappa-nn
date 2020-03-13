//
//  particle.hpp
//
// The Particle class loads and stores basic data for particles; the only particle that can be loaded only as a Particle object is the electron

#ifndef particle_hpp
#define particle_hpp

#include <armadillo>
#include <string>

namespace kappa {

    class Particle {
        // base class for particles
    public:
        // <constructor_docs>
        Particle(const std::string &name, const std::string &filename = "particles.yaml");
        /*
        <cdoc>
            <doc lang="rus">
                Создает объект типа Particle, загружая данные о частице из базы данных
                
                **Параметры**:
                
                    * ``const std::string &name`` - название молекулы

                    * ``const std::string &filename`` - путь к файлу базы данных
            </doc>
        </cdoc>
        */

        Particle();
        /*
        <cdoc>
            <doc lang="rus">
                Создает пустой объект типа Particle с нулевыми значениями всех параметров
            </doc>
        </cdoc>
        */
        // </constructor_docs>

        // <field_docs>
        std::string name;
        /*
        <ddoc>
            <doc lang="rus">
                Название частицы
            </doc>
        </ddoc>
        */

        double mass = 0.0;
        /*
        <ddoc>
            <doc lang="rus">
                Масса частицы
            </doc>
        </ddoc>
        */

        double diameter = 0.0;
        /*
        <ddoc>
            <doc lang="rus">
                Диаметр частицы
            </doc>
        </ddoc>
        */

        int charge = 0;
        /*
        <ddoc>
            <doc lang="rus">
                Заряд частицы (выраженный в элементарных электрических зарядах)
            </doc>
        </ddoc>
        */

        double formation_energy = 0.0;
        /*
        <ddoc>
            <doc lang="rus">
                Энергия образования частицы
            </doc>
        </ddoc>
        */
        
        double LennardJones_epsilon = 0.0;
        /*
        <ddoc>
            <doc lang="rus">
                Глубина потенциальной ямы в потенциале Леннарда-Джонса (если частица является электроном, считается равным 0)
            </doc>
        </ddoc>
        */

        double ionization_potential = 0.0;
        /*
        <ddoc>
            <doc lang="rus">
                Потенциал ионизации частицы (если частица является электроном, считается равным 0)
            </doc>
        </ddoc>
        */

        int num_electron_levels = 0;
        /*
        <ddoc>
            <doc lang="rus">
                Число электронных уровней частицы (если частица является электроном, считается равным 0)
            </doc>
        </ddoc>
        */

        arma::vec electron_energy;
        /*
        <ddoc>
            <doc lang="rus">
                Вектор энергий электронных уровней частицы (если частица является электроном, не инициализируется)
            </doc>
        </ddoc>
        */

        arma::Col<unsigned long> statistical_weight;
        /*
        <ddoc>
            <doc lang="rus">
                Вектор стат. весов электронных уровней частицы (если частица является электроном, не инициализируется)
            </doc>
        </ddoc>
        */
        // </field_docs>
        

    private:
        void readData(const std::string &name, const std::string &filename);
    };

}

#endif /* particle_hpp */
