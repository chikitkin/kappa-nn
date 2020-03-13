//
//  molecule.hpp
//
// The Molecule class loads and stores spectroscopic data for diatomic molecules

#ifndef molecule_hpp
#define molecule_hpp

#include <armadillo>
#include <map>
#include <string>
#include <vector>
#include "constants.h"
#include "particle.hpp"

namespace kappa {

    class Molecule : public Particle {
    public:
        // <constructor_docs>
        Molecule(const std::string &name, bool anharmonic_spectrum=true, bool rigid_rotator = true, const std::string &filename = "particles.yaml");
        /*
        <cdoc>
            <doc lang="rus">
                Создает объект типа Molecule, загружая данные о двухатомной молекуле из базы данных и производя расчет вращательных и колебательных спектров для каждого электронного состояния
                
                **Параметры**:
                
                    * ``const std::string &name`` - название молекулы

                    * ``bool anharmonic_spectrum`` - имеет ли молекула ангармонический колебательный спектр (значение по умолчанию ``true``)

                    * ``bool rigid_rotator`` - является ли молекула жестким ротатором (значение по умолчанию ``true``)

                    * ``const std::string &filename`` - путь к файлу базы данных
            </doc>
        </cdoc>
        */
        // </constructor_docs>

        // <field_docs>
        bool anharmonic_spectrum;
        /*
        <ddoc>
            <doc lang="rus">
                Является ли спектр ангармоническим
            </doc>
        </ddoc>
        */

        bool rigid_rotator;
        /*
        <ddoc>
            <doc lang="rus">
                Является ли молекула жестким ротатором
            </doc>
        </ddoc>
        */

        double reduced_osc_mass;
        /*
        <ddoc>
            <doc lang="rus">
                Приведенная масса осциллятора, равная <imath>m_{A}m_{B} / (m_{A} + m_{B})</imath>, где <imath>m_{A}</imath>, <imath>m_{B}</imath> - массы атомов молекулы
            </doc>
        </ddoc>
        */

        double mA_mAB;
        /*
        <ddoc>
            <doc lang="rus">
                Отношение массы первого атома молекулы к массе молекулы, равное <imath>m_{A} / (m_{A} + m_{B})</imath>, где <imath>m_{A}</imath>, <imath>m_{B}</imath> - массы атомов молекулы
            </doc>
        </ddoc>
        */

        double mB_mAB;
        /*
        <ddoc>
            <doc lang="rus">
                Отношение массы второго атома молекулы к массе молекулы, равное <imath>m_{A} / (m_{A} + m_{B})</imath>, где <imath>m_{A}</imath>, <imath>m_{B}</imath> - массы атомов молекулы
            </doc>
        </ddoc>
        */

        double rot_inertia;
        /*
        <ddoc>
            <doc lang="rus">
                Вращательный момент инерции молекулы
            </doc>
        </ddoc>
        */

        double internuclear_distance;
        /*
        <ddoc>
            <doc lang="rus">
                Межъядерное расстояние
            </doc>
        </ddoc>
        */

        int rot_symmetry;
        /*
        <ddoc>
            <doc lang="rus">
                Фактор симметрии молекул (равен 2 для гомоядерных и 1 для гетероядерных молекул)
            </doc>
        </ddoc>
        */

        arma::vec vibr_frequency;
        /*
        <ddoc>
            <doc lang="rus">
                Вектор колебательных частот (для каждого электронного состояния)
            </doc>
        </ddoc>
        */

        arma::vec vibr_we_xe;
        /*
        <ddoc>
            <doc lang="rus">
                Вектор параметра ангармоничности <imath>\omega_{e}x_{e}</imath> (для каждого электронного состояния)
            </doc>
        </ddoc>
        */

        arma::vec vibr_we_ye;
        /*
        <ddoc>
            <doc lang="rus">
                Вектор параметра ангармоничности <imath>\omega_{e}y_{e}</imath> (для каждого электронного состояния)
            </doc>
        </ddoc>
        */

        arma::vec vibr_we_ze;
        /*
        <ddoc>
            <doc lang="rus">
                Вектор параметра ангармоничности <imath>\omega_{e}z_{e}</imath> (для каждого электронного состояния)
            </doc>
        </ddoc>
        */

        arma::vec rot_be;
        /*
        <ddoc>
            <doc lang="rus">
                Вектор параметра <imath>B_{e}</imath>, описывающего зависимость вращательной энергии от вращательного уровня (для каждого электронного состояния)
            </doc>
        </ddoc>
        */

        arma::vec rot_ae;
        /*
        <ddoc>
            <doc lang="rus">
                Вектор параметра <imath>\alpha_{e}</imath>, описывающего зависимость вращательной энергии от вращательного и колебательного уровня (для каждого электронного состояния)
            </doc>
        </ddoc>
        */
        
        arma::vec characteristic_vibr_temperatures;
        /*
        <ddoc>
            <doc lang="rus">
                Вектор характеристических колебательных температур (для каждого электронного состояния)
            </doc>
        </ddoc>
        */
        arma::vec diss_energy;
        /*
        <ddoc>
            <doc lang="rus">
                Вектор энергий диссоциации (для каждого электронного состояния)
            </doc>
        </ddoc>
        */

        std::vector<std::vector<int> > num_rot_levels; 
        /*
        <ddoc>
            <doc lang="rus">
                Массив числа вращательных уровней (1й индекс - номер электронного уровня, 2й индекс - номер колебательного уровня; в случае, если молекула является жестким ротатором, для любого фиксированного 1го индекса
                все значения одинаков (число вращательных уровней не зависит от колебательного состояния молекулы))
            </doc>
        </ddoc>
        */
        
        std::vector<int> num_vibr_levels;
        /*
        <ddoc>
            <doc lang="rus">
                Массив числа колебательных уровней (для каждого электронного состояния)
            </doc>
        </ddoc>
        */

        std::vector<std::vector<arma::vec> > rot_energy;
        /*
        <ddoc>
            <doc lang="rus">
                Массив векторов энергий вращательных уровней (1й индекс - номер электронного уровня, 2й индекс - номер колебательного уровня; в случае, если молекула является жестким ротатором, для любого фиксированного 1го индекса
                все значения одинаковы (вращательный спектр не зависит от колебательного состояния молекулы))
            </doc>
        </ddoc>
        */

        std::vector<arma::vec> vibr_energy;
        /*
        <ddoc>
            <doc lang="rus">
                Вектор векторов энергий колебательных уровней (индекс - номер электронного уровня)
            </doc>
        </ddoc>
        */

        double parker_const;
        /*
        <ddoc>
            <doc lang="rus">
                Значения постоянной в формуле Паркера, используемой для расчета времен вращательной релаксации
            </doc>
        </ddoc>
        */
        // </field_docs>
        
    private:
        void readData(const std::string &name, const std::string &filename);
    };

}

#endif /* molecule_hpp */
