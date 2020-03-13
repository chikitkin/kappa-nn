//
//  interaction.hpp
//
// The Interaction class loads and stores data for particle-particle interactions

#ifndef kappa_interaction_hpp
#define kappa_interaction_hpp

#include <string>
#include <unordered_map>
#include "particles.hpp"
#include "../exceptions.hpp"

namespace kappa {

    enum class interaction_types {interaction_neutral_neutral, interaction_neutral_ion, interaction_neutral_electron, interaction_charged_charged}; 

    class Interaction {
        // class for interactions
    private:
        std::unordered_map<std::string, double> data;
        
        void readData(const std::string &name, const std::string &filename);
        
        
    public:
        // <constructor_docs>
        Interaction(const Particle &Particle1, const Particle &Particle2, const std::string &filename = "interaction.yaml");
        /*
        <cdoc>
            <doc lang="rus">
                Создает объект типа Interaction, загружая данные о взаимодействии двух частиц из базы данных
                
                **Параметры**:
                
                    * ``const Particle &Particle1`` - первая частица, участвующая во взаимодействии

                    * ``const Particle &Particle2`` - вторая частица, участвующая во взаимодействии

                    * ``const std::string &filename`` - путь к файлу базы данных
            </doc>
        </cdoc>
        */
        // </constructor_docs>

        // <field_docs>
        std::string particle1_name;
        /*
        <ddoc>
            <doc lang="rus">
                Название первой частицы, участвующей во взаимодействии
            </doc>
        </ddoc>
        */

        std::string particle2_name;
        /*
        <ddoc>
            <doc lang="rus">
                Название второй частицы, участвующей во взаимодействии
            </doc>
        </ddoc>
        */

        int charge1;
        /*
        <ddoc>
            <doc lang="rus">
                Заряд первой частицы, участвующей во взаимодействии (выраженный в элементарных электрических зарядах)
            </doc>
        </ddoc>
        */

         
        int charge2;
        /*
        <ddoc>
            <doc lang="rus">
                Заряд второй частицы, участвующей во взаимодействии (выраженный в элементарных электрических зарядах)
            </doc>
        </ddoc>
        */

        double collision_mass;
        /*
        <ddoc>
            <doc lang="rus">
                Приведеннная масса сталкивающихся частиц
            </doc>
        </ddoc>
        */

        double collision_diameter;
        /*
        <ddoc>
            <doc lang="rus">
                Диаметр столкновения частиц
            </doc>
        </ddoc>
        */

        double epsilon;
        /*
        <ddoc>
            <doc lang="rus">
                Глубина потенциальной ямы в потенциале Леннарда-Джонса 
            </doc>
        </ddoc>
        */

        bool vss_data;
        /*
        <ddoc>
            <doc lang="rus">
                Присутствуют ли данные потенциала VSS для данного взаимодействия (``true`` если присутствуют, ``false`` если отсутствуют)
            </doc>
        </ddoc>
        */

        double vss_dref;
        /*
        <ddoc>
            <doc lang="rus">
                Эталонный диаметр <imath>d_{ref}</imath> в потенциале VSS
            </doc>
        </ddoc>
        */

        double vss_omega;
        /*
        <ddoc>
            <i1><math>d = d_{ref} (g_{ref} / g)^(\omega - 0.5)</math></i1>

            <doc lang="rus">
                Параметр <imath>\omega</imath> в потенциале VSS. Диаметр столкновения выражается через относительную скорость сталкивающихся частиц <imath>g</imath> как <insert>i1</insert>
                где <imath>g_{ref}</imath> - эталонная скорость
            </doc>
        </ddoc>
        */

        double vss_alpha;
        /*
        <ddoc>
            <doc lang="rus">
                Значение параметра <imath>\alpha</imath>, описывающего рассеяние частиц, в потенциале VSS
            </doc>
        </ddoc>
        */


        double vss_Tref;
        /*
        <ddoc>
            <doc lang="rus">
                Эталонная температура <imath>T_{ref}</imath> в потенциале VSS
            </doc>
        </ddoc>
        */

        double vss_c_d;
        /*
        <ddoc>
            <i1><math>d = vss_{c,d} g^(0.5 - \omega)</math></i1>

            <doc lang="rus">
                Вспомогательная величина <imath>vss_{c,d}</imath>, диаметр столкновения выражается через нее по формуле <insert>i1</insert>
            </doc>
        </ddoc>
        */

        double vss_c_cs;
        /*
        <ddoc>
            <i1><math>\pi d^2 = vss_{c,cs} g^(1 - 2\omega)</math></i1>

            <doc lang="rus">
                Вспомогательная величина <imath>vss_{c,cs}</imath>, сечение столкновения выражается через нее по формуле <insert>i1</insert>
            </doc>
        </ddoc>
        */

        interaction_types interaction_type;
        /*
        <ddoc>
            <doc lang="rus">
                Тип взаимодействия, возможные значения:

                    * ``kappa::interaction_types::interaction_neutral_neutral`` - взаимодействие незаряженных частиц

                    * ``kappa::interaction_types::interaction_neutral_ion`` - взаимодействие нейтральной частиц и иона

                    * ``kappa::interaction_types::interaction_neutral_electron`` - взаимодействие нейтральной частиц и электрона

                    * ``kappa::interaction_types::interaction_charged_charged`` - взаимодействие двух заряженных частиц
            </doc>
        </ddoc>
        */
        // </field_docs>

        // arma::mat neutral_electron_elastic_coeffs;
        
        // <method_docs>
        const double& operator[](const std::string &name) const;
        /*
        <fdoc>
            <doc lang="rus">
                Доступ к скалярному параметру взаимодействия по названию параметра
            </doc>
        </fdoc>
        */

        const double& operator[](const char* name) const;
        /*
        <fdoc>
            <doc lang="rus">
                Доступ к скалярному параметру взаимодействия по названию параметра
            </doc>
        </fdoc>
        */
        // </method_docs>
    };
}

#endif /* interaction_hpp */
