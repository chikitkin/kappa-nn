////
////  approximation.hpp
////
//
///*
//<topdoc>
//    <label>approximation-label</label>
//    <type>classdoc</type>
//    <doc lang="eng">
//        The Mixture class contains
//    </doc>
//    <doc lang="rus">
//        Класс Mixture содержит методы
//    </doc>
//</topdoc>
//*/
//
//
//#ifndef kappa_mixture_hpp
//#define kappa_mixture_hpp
//
//#include <armadillo>
//#include <vector>
//#include "approximation.hpp"
//#include "numeric.hpp"
//#include "models.h"
//
//
//namespace kappa {
//
//    class Mixture : public Approximation  {
//    using Approximation::c_tr;
//    using Approximation::c_rot;
//    using Approximation::debye_length;
//    protected:
//        
//        int inter_index(int i, int j);
//        std::vector<std::string> split_string(std::string input);
//
//
//        arma::vec molecule_charges;
//        arma::vec atom_charges;
//        arma::vec molecule_charges_sq;
//        arma::vec atom_charges_sq; // squared charges
//
//        arma::mat omega_11;
//        arma::mat omega_12;
//        arma::mat omega_13;
//        arma::mat omega_22;
//        arma::mat rot_rel_times;
//        arma::vec c_rot_arr;
//        arma::vec c_rot_rigid_rot_arr;
//        
//        arma::mat shear_viscosity_LHS;
//        arma::mat thermal_conductivity_LHS;
//        arma::mat thermal_conductivity_rigid_rot_LHS;
//        arma::mat diffusion_LHS;
//        arma::mat diffusion_rigid_rot_LHS;
//        arma::mat bulk_viscosity_LHS;
//        arma::mat bulk_viscosity_rigid_rot_LHS;
//
//        arma::vec thermal_conductivity_RHS;
//        arma::vec thermal_conductivity_rigid_rot_RHS;
//        arma::vec diffusion_RHS;
//        arma::vec diffusion_rigid_rot_RHS;
//        arma::vec shear_viscosity_RHS;
//        arma::vec bulk_viscosity_RHS;
//        arma::vec bulk_viscosity_rigid_rot_RHS;
//
//        arma::vec thermal_conductivity_rigid_rot_coeffs;
//        arma::vec thermal_conductivity_coeffs;
//        arma::vec diffusion_coeffs;
//        arma::vec diffusion_rigid_rot_coeffs;
//        arma::vec shear_viscosity_coeffs;
//        arma::vec bulk_viscosity_coeffs;
//        arma::vec bulk_viscosity_rigid_rot_coeffs;
//
//        arma::vec empty_n_atom; // in case our mixture has no atoms, we pass this as a dummy parameter to functions which expect numeric density of atomic species
//        std::vector<arma::vec> empty_n_vl_molecule;
//
//        arma::vec this_n_atom;
//        std::vector<arma::vec> this_n_vl_mol;
//        arma::vec this_n_molecules;
//        double this_total_n;
//        double this_total_dens;
//        double this_ctr;
//        double this_crot;
//        double this_n_electrons;
//
//        std::vector<kappa::Molecule> molecules;
//        std::vector<kappa::Atom> atoms;
//        kappa::Particle electron;
//
//        int num_molecules;
//        int num_atoms;
//        int n_vibr_levels_total;
//        int n_particles;
//        bool cache_on;
//        bool all_rigid_rotators = true;
//        bool is_ionized;
//
//        std::vector<int> vl_offset; // vl_offset[0] = 0; vl_offset[1] = molecules[0].num_vibr_levels[0]; vl_offset[1] = molecules[0].num_vibr_levels[0] + molecules[1].num_vibr_levels[0], etc.
//        // this allows us to easily calculate which elements of state-to-state matrix we need to fill, given the molecule indices and the vibrational levels
//
//        std::vector<kappa::Interaction> interactions;
//        std::map<std::string, int> molecule_name_map;
//        std::map<std::string, int> atom_name_map;
//
//        void init_matrices(const std::string &particles_filename);
//        void add_interactions(const std::string &filename);
//
//        void check_n_vl_molecule(const std::vector<arma::vec> &n_vl_molecule); // check sizes of arrays and test values for non-negativity
//        void check_n_molecule(const arma::vec &n_molecule); // check sizes of arrays and test values for non-negativity
//        void check_n_atom(const arma::vec &n_atom);
//        void check_n(const arma::vec &n);
//
//        void compute_omega11(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//        void compute_omega12(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//        void compute_omega13(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//        void compute_omega22(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//
//        // Particle + e- interactions
//        void compute_omega11(double T, double debye_length, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//        void compute_omega12(double T, double debye_length, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//        void compute_omega13(double T, double debye_length, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//        void compute_omega22(double T, double debye_length, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//
//        void compute_c_rot_rigid_rot(double T);
//        void compute_c_rot(double T); // compute c_rot arrays
//        void compute_full_crot_rigid_rot(double T);
//        void compute_full_crot(double T); // compute c_rot of mixture
//        void compute_rot_rel_times(double T, double n, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//        void inplace_compute_n_molecule(const std::vector<arma::vec> &n_vl_molecule);
//
//        double th_cond, sh_visc, b_visc;
//        arma::vec th_diff;
//
//
//		const arma::mat &compute_bulk_viscosity_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//        const arma::vec &compute_bulk_viscosity_RHS(double T);
//        const arma::vec &compute_bulk_viscosity_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//
//        const arma::mat &compute_bulk_viscosity_rigid_rot_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//        const arma::vec &compute_bulk_viscosity_rigid_rot_RHS(double T);
//        const arma::vec &compute_bulk_viscosity_rigid_rot_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//
//        const arma::mat &compute_thermal_conductivity_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//        const arma::vec &compute_thermal_conductivity_RHS(double T);
//        const arma::vec &compute_thermal_conductivity_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//
//        const arma::mat &compute_thermal_conductivity_rigid_rot_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//		const arma::vec &compute_thermal_conductivity_rigid_rot_RHS(double T);
//		const arma::vec &compute_thermal_conductivity_rigid_rot_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//
//        const arma::mat &compute_diffusion_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//
//        const arma::mat &compute_shear_viscosity_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//        const arma::vec &compute_shear_viscosity_RHS(double T);
//        const arma::vec &compute_shear_viscosity_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//
//        double bulk_viscosity(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//		double thermal_conductivity(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);    
//        void thermodiffusion(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//		double shear_viscosity(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
//        
//    public:
//        // <constructor_docs>
//        Mixture(const std::vector<kappa::Molecule> &i_molecules, const std::vector<kappa::Atom> &i_atoms, const std::string &interactions_filename="interaction.yaml", const std::string &particles_filename="particles.yaml");
//        /*
//        <cdoc>
//            <doc lang="rus">
//                Создает объект класса Mixture с указанными молекулами и атомами. В случае, если среди молекул и атомов есть заряженные частицы, считается, что смесь также содержит электроны.
//
//                **Параметры**:
//
//                    * ``const std::vector<kappa::Molecule> &i_molecules`` - Вектор молекул, которые будут в смеси
//
//                    * ``const std::vector<kappa::Atom> &i_atoms`` - Вектор атомов, которые будут в смеси
//
//                    * ``const std::string &interactions_filename`` - путь к файлу базы данных взаимодействий (значение по умолчанию ``"interaction.yaml"``)
//
//                    * ``const std::string &particles_filename`` - путь к файлу базы данных частиц, необходимо для загрузки данных об электроне (значение по умолчанию ``"particles.yaml"``)
//            </doc>
//        </cdoc>
//        */
//
//        Mixture(const std::vector<kappa::Molecule> &i_molecules, const std::string &interactions_filename="interaction.yaml", const std::string &particles_filename="particles.yaml");
//        /*
//        <cdoc>
//            <doc lang="rus">
//                Создает объект класса Mixture с указанными молекулами. В случае, если среди молекул есть заряженные частицы, считается, что смесь также содержит электроны.
//
//                **Параметры**:
//
//                    * ``const std::vector<kappa::Molecule> &i_molecules`` - Вектор молекул, которые будут в смеси
//
//                    * ``const std::string &interactions_filename`` - путь к файлу базы данных взаимодействий (значение по умолчанию ``"interaction.yaml"``)
//
//                    * ``const std::string &particles_filename`` - путь к файлу базы данных частиц, необходимо для загрузки данных об электроне (значение по умолчанию ``"particles.yaml"``)
//            </doc>
//        </cdoc>
//        */
//
//        Mixture(const std::vector<kappa::Atom> &i_atoms, const std::string &interactions_filename="interaction.yaml", const std::string &particles_filename="particles.yaml");
//        /*
//        <cdoc>
//            <doc lang="rus">
//                Создает объект класса Mixture с указанными атомами. В случае, если среди атомов есть заряженные частицы, считается, что смесь также содержит электроны.
//
//                **Параметры**:
//
//                    * ``const std::vector<kappa::Atom> &i_atoms`` - Вектор атомов, которые будут в смеси
//
//                    * ``const std::string &interactions_filename`` - путь к файлу базы данных взаимодействий (значение по умолчанию ``"interaction.yaml"``)
//
//                    * ``const std::string &particles_filename`` - путь к файлу базы данных частиц, необходимо для загрузки данных об электроне (значение по умолчанию ``"particles.yaml"``)
//            </doc>
//        </cdoc>
//        */
//
//        Mixture(const kappa::Molecule &molecule, const kappa::Atom &atom, const std::string &interactions_filename="interaction.yaml", const std::string &particles_filename="particles.yaml");
//        /*
//        <cdoc>
//            <doc lang="rus">
//                Создает объект класса Mixture с указанными молекулой и атомом. В случае, если среди них есть заряженные частицы, считается, что смесь также содержит электроны.
//
//                **Параметры**:
//
//                    * ``const kappa::Molecule &molecule`` - Молекула, котоая будет в смеси
//
//                    * ``const kappa::Atom &atom`` - Атом, который будет в смеси
//
//                    * ``const std::string &interactions_filename`` - путь к файлу базы данных взаимодействий (значение по умолчанию ``"interaction.yaml"``)
//
//                    * ``const std::string &particles_filename`` - путь к файлу базы данных частиц, необходимо для загрузки данных об электроне (значение по умолчанию ``"particles.yaml"``)
//            </doc>
//        </cdoc>
//        */
//
//        Mixture(const kappa::Molecule &molecule, const std::string &interactions_filename="interaction.yaml", const std::string &particles_filename="particles.yaml");
//        /*
//        <cdoc>
//            <doc lang="rus">
//                Создает объект класса Mixture с указанной молекулой. В случае, если она заряжена, считается, что смесь также содержит электроны.
//
//                **Параметры**:
//
//                    * ``const kappa::Molecule &molecule`` - Молекула, котоая будет в смеси
//
//                    * ``const std::string &interactions_filename`` - путь к файлу базы данных взаимодействий (значение по умолчанию ``"interaction.yaml"``)
//
//                    * ``const std::string &particles_filename`` - путь к файлу базы данных частиц, необходимо для загрузки данных об электроне (значение по умолчанию ``"particles.yaml"``)
//            </doc>
//        </cdoc>
//        */
//
//        Mixture(const kappa::Atom &atom, const std::string &interactions_filename="interaction.yaml", const std::string &particles_filename="particles.yaml");
//        /*
//        <cdoc>
//            <doc lang="rus">
//                Создает объект класса Mixture с указанным атомом. В случае, если он заряжен, считается, что смесь также содержит электроны.
//
//                **Параметры**:
//
//                    * ``const kappa::Atom &atom`` - Атом, который будет в смеси
//
//                    * ``const std::string &interactions_filename`` - путь к файлу базы данных взаимодействий (значение по умолчанию ``"interaction.yaml"``)
//
//                    * ``const std::string &particles_filename`` - путь к файлу базы данных частиц, необходимо для загрузки данных об электроне (значение по умолчанию ``"particles.yaml"``)
//            </doc>
//        </cdoc>
//        */
//
//        Mixture(const std::string particle_names, const std::string &interactions_filename="interaction.yaml", const std::string &particles_filename="particles.yaml", bool anharmonic=true, bool rigid_rotators=true);
//        /*
//        <cdoc>
//            <doc lang="rus">
//                Создает объект класса Mixture с указанными частицами. В случае, если среди них есть заряженные частицы, считается, что смесь также содержит электроны.
//
//                **Параметры**:
//
//                    * ``const std::string particle_names`` - строка с названиями частиц, разделенными запятыми
//
//                    * ``const std::string &interactions_filename`` - путь к файлу базы данных взаимодействий (значение по умолчанию ``"interaction.yaml"``)
//
//                    * ``const std::string &particles_filename`` - путь к файлу базы данных частиц, необходимо для загрузки данных об электроне (значение по умолчанию ``"particles.yaml"``)
//
//                    * ``bool anharmonic`` - имеют ли молекулы ангармонический колебательный спектр (значение по умолчанию ``true``)
//
//                    * ``bool rigid_rotatora`` - являются ли молекулы жесткими ротаторами (значение по умолчанию ``true``)
//            </doc>
//        </cdoc>
//        */
//        // </constructor_docs>
//
//        // <method_docs>
//        std::string get_names();
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает строку с названиями частиц в смеси
//            </doc>
//        </fdoc>
//        */
//
//        int get_n_particles();
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает число частиц в смеси
//            </doc>
//        </fdoc>
//        */
//
//        int get_n_vibr_levels();
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает суммарное число колебательных уровней в смеси
//            </doc>
//        </fdoc>
//        */
//
//        std::vector<int> get_n_vibr_levels_array();
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает массив числа колебательных уровней молекул в смеси
//            </doc>
//        </fdoc>
//        */
//
//        arma::vec convert_molar_frac_to_mass(const arma::vec &x);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Переводит массив молярных долей в массив массовых долей
//
//                **Параметры**:
//
//                    * ``const arma::vec &x`` - массив молярных долей компонент смеси
//
//            </doc>
//        </fdoc>
//        */
//
//        arma::vec convert_mass_frac_to_molar(const arma::vec &y);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Переводит массив массовых долей в массив молярных долей
//
//                **Параметры**:
//
//                    * ``const arma::vec &y`` - массив массовых долей компонент смеси
//
//            </doc>
//        </fdoc>
//        */
//        
//        kappa::Molecule molecule(const std::string &name);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает объект, соответствующей молекуле, имеющейся в смеси
//
//                **Параметры**:
//
//                    * ``const std::string &name`` - название молекулы
//
//            </doc>
//        </fdoc>
//        */
//
//        kappa::Atom atom(const std::string &name);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает объект, соответствующей атому, имеющемуся в смеси
//
//                **Параметры**:
//
//                    * ``const std::string &name`` - название молекулы
//
//            </doc>
//        </fdoc>
//        */
//        
//        kappa::Interaction interaction(const kappa::Molecule &molecule1, const kappa::Molecule &molecule2);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает объект взаимодействия двух молекул, имеющихся в смеси
//
//                **Параметры**:
//
//                    * ``const kappa::Molecule &molecule1`` - первая молекула
//
//                    * ``const kappa::Molecule &molecule2`` - вторая молекула
//
//            </doc>
//        </fdoc>
//        */
//
//        kappa::Interaction interaction(const kappa::Molecule &molecule, const kappa::Atom &atom);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает объект взаимодействия молекулы и атома, имеющихся в смеси
//
//                **Параметры**:
//
//                    * ``const kappa::Molecule &molecule`` - молекула
//
//                    * ``const kappa::Atom &atom`` - атом
//
//            </doc>
//        </fdoc>
//        */
//
//        kappa::Interaction interaction(const kappa::Atom &atom, const kappa::Molecule &molecule);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает объект взаимодействия молекулы и атома, имеющихся в смеси
//
//                **Параметры**:
//
//                    * ``const kappa::Atom &atom`` - атом
//
//                    * ``const kappa::Molecule &molecule`` - молекула
//
//            </doc>
//        </fdoc>
//        */
//
//        kappa::Interaction interaction(const kappa::Atom &atom1, const kappa::Atom &atom2);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает объект взаимодействия двух атомов, имеющихся в смеси
//
//                **Параметры**:
//
//                    * ``const kappa::Molecule &molecule1`` - первая молекула
//
//                    * ``const kappa::Molecule &molecule2`` - вторая молекула
//
//            </doc>
//        </fdoc>
//        */
//
//        double debye_length(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons);
//        /*
//        <fdoc>
//            <doc lang="rus">
//               Расчет Дебаевской длины
//
//                **Параметры**:
//
//                    * ``double T`` - температура
//
//                    * ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//                    * ``const arma::vec &n_atom`` - Вектор числовых плотностей атомов
//
//                    * ``double n_electrons`` - числовая плотность электронов
//
//            </doc>
//        </fdoc>
//        */
//
//        double debye_length(double T, const std::vector<arma::vec> &n_vl_molecule, double n_electrons);
//        /*
//        <fdoc>
//            <doc lang="rus">
//               Расчет Дебаевской длины (в смеси без атомов)
//
//                **Параметры**:
//
//                    * ``double T`` - температура
//
//                    * ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//                    * ``double n_electrons`` - числовая плотность электронов
//
//            </doc>
//        </fdoc>
//        */
//
//        double debye_length(double T, const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons);
//        /*
//        <fdoc>
//            <doc lang="rus">
//               Расчет Дебаевской длины
//
//                **Параметры**:
//
//                    * ``double T`` - температура
//
//                    * ``const arma::vec &n_molecule`` - Вектор числовых плотностей молекул
//
//                    * ``const arma::vec &n_atom`` - Вектор числовых плотностей атомов
//
//                    * ``double n_electrons`` - числовая плотность электронов
//
//            </doc>
//        </fdoc>
//        */
//
//        double debye_length(double T, const arma::vec &n);
//        /*
//        <fdoc>
//            <doc lang="rus">
//               Расчет Дебаевской длины
//
//                **Параметры**:
//
//                    * ``const arma::vec &n`` - Вектор числовых плотностей (сначала в векторе идут плотности молекул, затем атомов, затем электронов)
//					  (сначала в векторе идут плотности молекул, затем атомов, затем электронов, порядок молекул и атомов должен быть таким же, что использовался и при создании смеси)
//
//            </doc>
//        </fdoc>
//        */
//
//        double compute_n(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons=0.0);
//		/*
//		<fdoc>
//			<doc lang="rus">
//				Вычисление числовой плотности смеси
//
//				**Параметры**:
//
//					* ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//					* ``const arma::vec &n_atom`` - Вектор числовых плотностей атомов
//
//                    * ``double n_electrons`` - числовая плотность электронов (значение по умолчанию - 0)
//			</doc>
//		</fdoc>
//		*/
//
//        double compute_n(const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons=0.0);
//        /*
//		<fdoc>
//			<doc lang="rus">
//				Вычисление числовой плотности смеси
//
//				**Параметры**:
//
//					* ``const arma::vec &n_molecule`` - Вектор числовых плотностей молекул
//
//                    * ``const arma::vec &n_atom`` - Вектор числовых плотностей атомов
//
//                    * ``double n_electrons`` - числовая плотность электронов (значение по умолчанию - 0)
//		    </doc>
//        </fdoc>
//        */
//
//        double compute_n(const std::vector<arma::vec> &n_vl_molecule, double n_electrons=0.0);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление числовой плотности смеси
//
//                **Параметры**:
//
//                    * ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//                    * ``double n_electrons`` - числовая плотность электронов (значение по умолчанию - 0)
//            </doc>
//        </fdoc>
//        */
//
//        double compute_n(const arma::vec &n);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление числовой плотности смеси
//
//                **Параметры**:
//
//                    * ``const arma::vec &n`` - Вектор числовых плотностей (сначала в векторе идут плотности молекул, затем атомов, затем электронов)
//					  (сначала в векторе идут плотности молекул, затем атомов, затем электронов, порядок молекул и атомов должен быть таким же, что использовался и при создании смеси)
//            </doc>
//        </fdoc>
//        */
//
//		arma::vec compute_n_molecule(const std::vector<arma::vec> &n_vl_molecule);
//		/*
//		<fdoc>
//			<doc lang="rus">
//				Вычисление числовой плотности молекул 
//
//				**Параметры**:
//
//					* ``const std::vector<arma::vec> &n_vl_molecule`` -  Массив заселенностей колебательных уровней молекул
//			</doc>
//		</fdoc>
//		*/
//        
//        double compute_density(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons=0.0);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление плотности смеси
//                    
//                **Параметры**:
//
//                    * ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//                    * ``const arma::vec &n_atom`` - Вектор числовых плотностей атомов
//
//                    * ``double n_electrons`` - числовая плотность электронов (значение по умолчанию - 0)
//            </doc>
//        </fdoc>
//        */
//
//        double compute_density(const std::vector<arma::vec> &n_vl_molecule, double n_electrons=0.0);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление плотности смеси
//                    
//                **Параметры**:
//
//                    * ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//                    * ``double n_electrons`` - числовая плотность электронов (значение по умолчанию - 0)
//            </doc>
//        </fdoc>
//        */
//
//        double compute_density(const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons=0.0);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление плотности смеси
//                    
//                **Параметры**:
//
//                    * ``const arma::vec &n_molecule`` - Вектор числовых плотностей молекул
//
//                    * ``const arma::vec &n_atom`` - Вектор числовых плотностей атомов
//
//                    * ``double n_electrons`` - числовая плотность электронов (значение по умолчанию - 0)
//            </doc>
//        </fdoc>
//        */
//
//        double compute_density(const arma::vec &n);
//		/*
//		<fdoc>
//			<doc lang="rus">
//				Вычисление плотности смеси
//
//				**Параметры**:
//			
//					* ``const arma::vec &n`` - Вектор числовых плотностей (сначала в векторе идут плотности молекул, затем атомов, затем электронов)
//					  (сначала в векторе идут плотности молекул, затем атомов, затем электронов, порядок молекул и атомов должен быть таким же, что использовался и при создании смеси)
//			
//			</doc>
//		</fdoc>
//		*/
//
//        double c_tr(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons=0.0);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление удельной теплоемкости поступательных степеней свободы 
//
//                **Параметры**:
//
//                    * ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//                    * ``const arma::vec &n_atom`` - Вектор числовых плотностей атомов
//
//                    * ``double n_electrons`` - числовая плотность электронов (значение по умолчанию - 0)
//            </doc>
//        </fdoc>
//        */
//
//        double c_tr(const std::vector<arma::vec> &n_vl_molecule, double n_electrons=0.0);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление удельной теплоемкости поступательных степеней свободы 
//
//                **Параметры**:
//
//                    * ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//                    * ``double n_electrons`` - числовая плотность электронов (значение по умолчанию - 0)
//            </doc>
//        </fdoc>
//        */
//        
//        double c_tr(const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons=0.0);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление удельной теплоемкости поступательных степеней свободы 
//
//                **Параметры**:
//
//                    * ``const arma::vec &n_molecule`` - Вектор числовых плотностей молекул
//
//                    * ``const arma::vec &n_atom`` - Вектор числовых плотностей атомов
//
//                    * ``double n_electrons`` - числовая плотность электронов (значение по умолчанию - 0)
//            </doc>
//        </fdoc>
//        */
//
//        double c_tr(const arma::vec &n);
//		/*
//		<fdoc>
//			<doc lang="rus">
//				Вычисление удельной теплоемкости поступательных степеней свободы 
//					
//				**Параметры**:
//
//					* ``const arma::vec &n`` - Вектор числовых плотностей частиц смеси
//					  (сначала в векторе идут плотности молекул, затем атомов, затем электронов, порядок молекул и атомов должен быть таким же, что использовался и при создании смеси)
//			</doc>
//		</fdoc>
//		*/
//
//        double c_rot(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons=0.0);
//		/*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление удельной теплоемкости вращательных степеней свободы 
//
//                **Параметры**:
//
//                	* ``double T`` - Температура
//
//                    * ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//                    * ``const arma::vec &n_atom`` - Вектор числовых плотностей атомов
//
//                    * ``double n_electrons`` - числовая плотность электронов (значение по умолчанию - 0)
//            </doc>
//        </fdoc>
//        */
//
//        double c_rot(double T, const std::vector<arma::vec> &n_vl_molecule, double n_electrons=0.0);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление удельной теплоемкости вращательных степеней свободы 
//
//                **Параметры**:
//
//                	* ``double T`` - Температура
//
//                    * ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//                    * ``double n_electrons`` - числовая плотность электронов (значение по умолчанию - 0)
//            </doc>
//        </fdoc>
//        */
//
//		double c_rot(double T, const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons=0.0);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление удельной теплоемкости вращательных степеней свободы 
//
//                **Параметры**:
//
//                	* ``double T`` - Температура
//
//                    * ``const arma::vec &n_molecule`` - Вектор числовых плотностей молекул
//
//                    * ``const arma::vec &n_atom`` - Вектор числовых плотностей атомов
//
//                    * ``double n_electrons`` - числовая плотность электронов (значение по умолчанию - 0)
//            </doc>
//        </fdoc>
//        */
//
//        double c_rot(double T, const arma::vec &n);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление удельной теплоемкости вращательных степеней свободы 
//                    
//                **Параметры**:
//
//                	* ``double T`` - Температура
//
//                    * ``const arma::vec &n`` - Вектор числовых плотностей частиц смеси (сначала в векторе идут плотности молекул, затем атомов, затем электронов)
//					  (сначала в векторе идут плотности молекул, затем атомов, затем электронов, порядок молекул и атомов должен быть таким же, что использовался и при создании смеси)
//            </doc>
//        </fdoc>
//        */
//
//        void compute_transport_coefficients(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons, kappa::models_omega model=kappa::models_omega::model_omega_esa, double perturbation=1e-9);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление коэффициентов переноса
//                    
//                **Параметры**:
//
//                	* ``double T`` - Температура
//
//                    * ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//                    * ``const arma::vec &n_atom`` - Вектор числовых плотностей атомов
//
//                    * ``double n_electrons`` - числовая плотность электронов
//
//                    * ``kappa::models_omega model`` - модель для расчета Омега-интегралов (значение модели по умолчанию ``kappa::models_omega::model_omega_esa``):
//                        
//                        * ``kappa::models_omega::model_omega_rs`` - модель твердых сфер (RS)
//                        
//                        * ``kappa::models_omega::model_omega_vss`` - модель сфер переменного диаметра (VSS)
//                        
//                        * ``kappa::models_omega::model_omega_bornmayer`` - модель Борна–Майера
//                        
//                        * ``kappa::models_omega::model_omega_lennardjones`` - модель Леннарда–Джонса
//                        
//                        * ``kappa::models_omega::model_omega_esa`` - модель ESA (феноменологическая модель)
//
//                    * ``double perturbation`` - малое возмущение, используемое для избежания вырождения систем при малых (исчезающих) концентрациях компонент (значение по умолчанию - 1e-9)
//            </doc>
//        </fdoc>
//        */
//
//        void compute_transport_coefficients(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, kappa::models_omega model=kappa::models_omega::model_omega_esa, double perturbation=1e-9);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление коэффициентов переноса в смеси без электронов
//                    
//                **Параметры**:
//
//                	* ``double T`` - Температура
//
//                    * ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//                    * ``const arma::vec &n_atom`` - Вектор числовых плотностей атомов
//
//                    * ``kappa::models_omega model`` - модель для расчета Омега-интегралов (значение модели по умолчанию ``kappa::models_omega::model_omega_esa``):
//                        
//                        * ``kappa::models_omega::model_omega_rs`` - модель твердых сфер (RS)
//                        
//                        * ``kappa::models_omega::model_omega_vss`` - модель сфер переменного диаметра (VSS)
//                        
//                        * ``kappa::models_omega::model_omega_bornmayer`` - модель Борна–Майера
//                        
//                        * ``kappa::models_omega::model_omega_lennardjones`` - модель Леннарда–Джонса
//                        
//                        * ``kappa::models_omega::model_omega_esa`` - модель ESA (феноменологическая модель)
//
//                    * ``double perturbation`` - малое возмущение, используемое для избежания вырождения систем при малых (исчезающих) концентрациях компонент (значение по умолчанию - 1e-9)
//            </doc>
//        </fdoc>
//        */
//
//        void compute_transport_coefficients(double T, const std::vector<arma::vec> &n_vl_molecule, double n_electrons, kappa::models_omega model=kappa::models_omega::model_omega_esa, double perturbation=1e-9);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление коэффициентов переноса в смеси без атомов
//                    
//                **Параметры**:
//
//                	* ``double T`` - Температура
//
//                    * ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//                    * ``double n_electrons`` - числовая плотность электронов
//
//                    * ``kappa::models_omega model`` - модель для расчета Омега-интегралов (значение модели по умолчанию ``kappa::models_omega::model_omega_esa``):
//                        
//                        * ``kappa::models_omega::model_omega_rs`` - модель твердых сфер (RS)
//                        
//                        * ``kappa::models_omega::model_omega_vss`` - модель сфер переменного диаметра (VSS)
//                        
//                        * ``kappa::models_omega::model_omega_bornmayer`` - модель Борна–Майера
//                        
//                        * ``kappa::models_omega::model_omega_lennardjones`` - модель Леннарда–Джонса
//                        
//                        * ``kappa::models_omega::model_omega_esa`` - модель ESA (феноменологическая модель)
//
//                    * ``double perturbation`` - малое возмущение, используемое для избежания вырождения систем при малых (исчезающих) концентрациях компонент (значение по умолчанию - 1e-9)
//            </doc>
//        </fdoc>
//        */
//
//        void compute_transport_coefficients(double T, const std::vector<arma::vec> &n_vl_molecule, kappa::models_omega model=kappa::models_omega::model_omega_esa, double perturbation=1e-9);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление коэффициентов переноса в смеси без атомов и электронов
//                    
//                **Параметры**:
//
//                	* ``double T`` - Температура
//
//                    * ``const std::vector<arma::vec> &n_vl_molecule`` - Массив заселенностей колебательных уровней молекул
//
//                    * ``kappa::models_omega model`` - модель для расчета Омега-интегралов (значение модели по умолчанию ``kappa::models_omega::model_omega_esa``):
//                        
//                        * ``kappa::models_omega::model_omega_rs`` - модель твердых сфер (RS)
//                        
//                        * ``kappa::models_omega::model_omega_vss`` - модель сфер переменного диаметра (VSS)
//                        
//                        * ``kappa::models_omega::model_omega_bornmayer`` - модель Борна–Майера
//                        
//                        * ``kappa::models_omega::model_omega_lennardjones`` - модель Леннарда–Джонса
//                        
//                        * ``kappa::models_omega::model_omega_esa`` - модель ESA (феноменологическая модель)
//
//                    * ``double perturbation`` - малое возмущение, используемое для избежания вырождения систем при малых (исчезающих) концентрациях компонент (значение по умолчанию - 1e-9)
//            </doc>
//        </fdoc>
//        */
//
//        void compute_transport_coefficients(double T, const arma::vec &n, kappa::models_omega model=kappa::models_omega::model_omega_esa, double perturbation=1e-9);
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Вычисление коэффициентов переноса в смеси
//                    
//                **Параметры**:
//
//                	* ``double T`` - Температура
//
//					* ``const arma::vec &n`` - Вектор числовых плотностей частиц смеси
//					  (сначала в векторе идут плотности молекул, затем атомов, затем электронов, порядок молекул и атомов должен быть таким же, что использовался и при создании смеси)
//
//                    * ``kappa::models_omega model`` - модель для расчета Омега-интегралов (значение модели по умолчанию ``kappa::models_omega::model_omega_esa``):
//                        
//                        * ``kappa::models_omega::model_omega_rs`` - модель твердых сфер (RS)
//                        
//                        * ``kappa::models_omega::model_omega_vss`` - модель сфер переменного диаметра (VSS)
//                        
//                        * ``kappa::models_omega::model_omega_bornmayer`` - модель Борна–Майера
//                        
//                        * ``kappa::models_omega::model_omega_lennardjones`` - модель Леннарда–Джонса
//                        
//                        * ``kappa::models_omega::model_omega_esa`` - модель ESA (феноменологическая модель)
//
//                    * ``double perturbation`` - малое возмущение, используемое для избежания вырождения систем при малых (исчезающих) концентрациях компонент (значение по умолчанию - 1e-9)
//            </doc>
//        </fdoc>
//        */
//
//        double get_thermal_conductivity();
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает коэффициент теплопроводности (для расчета необходимо вызвать функцию ``compute_transport_coefficients``)
//            </doc>
//        </fdoc>
//        */
//
//        double get_shear_viscosity();
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает коэффициент сдиговой вязкости (для расчета необходимо вызвать функцию ``compute_transport_coefficients``)
//            </doc>
//        </fdoc>
//        */
//
//        double get_bulk_viscosity();
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает коэффициент объемной вязкости (для расчета необходимо вызвать функцию ``compute_transport_coefficients``)
//            </doc>
//        </fdoc>
//        */
//
//        arma::vec get_thermodiffusion();
//        /*
//        <fdoc>
//            <doc lang="rus">
//                Возвращает вектор коэффициентов термодиффузии (для расчета необходимо вызвать функцию ``compute_transport_coefficients``).
//                Сначала идут коэффициенты термодиффузии каждого колебательного уровня каждой молекулы в смеси, затем коэффициенты термодиффузии атомов в смеси,
//                затем (при наличии) - электронов; коэффициенты для молекул и атомов идут в том же порядке, что и при создании смеси
//            </doc>
//        </fdoc>
//        */
//        // </method_docs>
//
//        // void compute_transport_coefficients(double T, const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons, kappa::models_omega model=kappa::models_omega::model_omega_esa, double perturbation=1e-9); TODO: Write!
//        // void compute_transport_coefficients(double T, const arma::vec &n_molecule, const arma::vec &n_atom, kappa::models_omega model=kappa::models_omega::model_omega_esa, double perturbation=1e-9);
//		arma::mat get_diffusion();
//    };
//}
//
//#endif /* kappa_mixture_hpp */
/*!
\file mixture.hpp
*/

#ifndef kappa_mixture_hpp
#define kappa_mixture_hpp

#include <vector>    
#include <map>
#include <numeric>
#include <string>

#define ARMA_DONT_PRINT_ERRORS // to avoid warning message during solve() linear systems

#include <armadillo> 

#include "approximation.hpp"
#include "numeric.hpp"
#include "models.h"

namespace kappa {

	class Mixture : public Approximation {

		using Approximation::c_tr;
		using Approximation::c_rot;
		using Approximation::debye_length;

	public:

		Mixture(const std::vector<kappa::Molecule> &i_molecules, const std::vector<kappa::Atom> &i_atoms, const std::string &interactions_filename = "interaction.yaml", const std::string &particles_filename = "particles.yaml");
		Mixture(const std::vector<kappa::Molecule> &i_molecules, const std::string &interactions_filename = "interaction.yaml", const std::string &particles_filename = "particles.yaml");
		Mixture(const std::vector<kappa::Atom> &i_atoms, const std::string &interactions_filename = "interaction.yaml", const std::string &particles_filename = "particles.yaml");
		Mixture(const kappa::Molecule &molecule, const kappa::Atom &atom, const std::string &interactions_filename = "interaction.yaml", const std::string &particles_filename = "particles.yaml");
		Mixture(const kappa::Atom &atom, const std::string &interactions_filename = "interaction.yaml", const std::string &particles_filename = "particles.yaml");
		Mixture(const kappa::Molecule &molecule, const std::string &interactions_filename = "interaction.yaml", const std::string &particles_filename = "particles.yaml");
		Mixture(const std::string particle_names, const std::string &interactions_filename = "interaction.yaml", const std::string &particles_filename = "particles.yaml", bool anharmonic = true, bool rigid_rotators = true);

		// Returns a string with the names of particles in the mixture
		std::string get_names();

		// Returns the number of particles in the mixture
		int get_n_particles();

		// Returns the sum of the number of vibrational levels in the mixtures
		int get_n_vibr_levels();

		// Returns an array of the number of vibrational levels of molecules in the mixture
		std::vector<int> get_n_vibr_levels_array();

		// Translates an array of molar fractions into an array of mass fractions
		// @param const arma::vec &x - array of molar fractions of the mixture
		arma::vec convert_molar_frac_to_mass(const arma::vec &x);

		// Converts an array of mass fractions into an array of molar fractions
		// @param const arma::vec &y - array of mass fractions of the mixture
		arma::vec convert_mass_frac_to_molar(const arma::vec &y);

		// Returns an object corresponding to the molecule present in the mixture
		// @param const std::string &name - name of the molecule
		kappa::Molecule molecule(const std::string &name);

		// Returns an object corresponding to the atom present in the mixture
		// @param const std::string &name - name of the atom
		kappa::Atom atom(const std::string &name);

		kappa::Interaction interaction(const kappa::Molecule &molecule1, const kappa::Molecule &molecule2);
		kappa::Interaction interaction(const kappa::Molecule &molecule, const kappa::Atom &atom);
		kappa::Interaction interaction(const kappa::Atom &atom, const kappa::Molecule &molecule);
		kappa::Interaction interaction(const kappa::Atom &atom1, const kappa::Atom &atom2);

		double debye_length(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons);
		double debye_length(double T, const std::vector<arma::vec> &n_vl_molecule, double n_electrons);
		double debye_length(double T, const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons);
		double debye_length(double T, const arma::vec &n);

		// Calculation of the number density of the mixture
		// @param const std::vector<arma::vec> &n_vl_molecule - the population of the vibrational levels of molecules
		// @param const arma::vec &n_molecule - Vector of number densities of molecules
		// @param const arma::vec &n_atom - Vector of number densities of atoms
		// @param double n_electrons - the number density of electrons (the default value is 0)
		// @param const arma::vec &n - Vector of number densities
		double compute_n(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons = 0.0);
		double compute_n(const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons = 0.0);
		double compute_n(const std::vector<arma::vec> &n_vl_molecule, double n_electrons = 0.0);
		double compute_n(const arma::vec &n);

		// Calculation of the number density of molecules
		arma::vec compute_n_molecule(const std::vector<arma::vec> &n_vl_molecule);

		// Number density of the mixture
		arma::vec compute_density_array(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons = 0.0);
		arma::vec compute_density_array(const std::vector<arma::vec> &n_vl_molecule);
		double compute_density(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons = 0.0);
		double compute_density(const std::vector<arma::vec> &n_vl_molecule, double n_electrons = 0.0);
		double compute_density(const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons = 0.0);
		double compute_density(const arma::vec &n);

		// Compute the pressure of the mixture
		double compute_pressure(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons = 0.0);

		// Calculation of the specific heat of translational degrees of freedom
		double c_tr(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons = 0.0);
		double c_tr(const std::vector<arma::vec> &n_vl_molecule, double n_electrons = 0.0);
		double c_tr(const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons = 0.0);
		double c_tr(const arma::vec &n);

		// Calculation of the specific heat of rotational degrees of freedom
		double c_rot(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons = 0.0);
		double c_rot(double T, const std::vector<arma::vec> &n_vl_molecule, double n_electrons = 0.0);
		double c_rot(double T, const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons = 0.0);
		double c_rot(double T, const arma::vec &n);

		void compute_transport_coefficients(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons, kappa::models_omega model = kappa::models_omega::model_omega_esa, double perturbation = 1e-9);
		void compute_transport_coefficients(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, kappa::models_omega model = kappa::models_omega::model_omega_esa, double perturbation = 1e-9);
		// void compute_transport_coefficients(double T, const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons, kappa::models_omega model=kappa::models_omega::model_omega_esa, double perturbation=1e-9);	
		// void compute_transport_coefficients(double T, const arma::vec &n_molecule, const arma::vec &n_atom, kappa::models_omega model=kappa::models_omega::model_omega_esa, double perturbation=1e-9);
		void compute_transport_coefficients(double T, const std::vector<arma::vec> &n_vl_molecule, double n_electrons, kappa::models_omega model = kappa::models_omega::model_omega_esa, double perturbation = 1e-9);
		void compute_transport_coefficients(double T, const std::vector<arma::vec> &n_vl_molecule, kappa::models_omega model = kappa::models_omega::model_omega_esa, double perturbation = 1e-9);
		void compute_transport_coefficients(double T, const arma::vec &n, kappa::models_omega model = kappa::models_omega::model_omega_esa, double perturbation = 1e-9);

		double get_thermal_conductivity();
		double get_shear_viscosity();
		double get_bulk_viscosity();
		arma::vec get_thermodiffusion();

		// Returns the matrix of diffusion coefficients (for calculation it is necessary to call the function compute_transport_coefficients)
		arma::mat get_diffusion();
		arma::mat get_lite_diffusion();
		arma::vec get_binary_diffusion();
		// arma::mat binary_diffusion(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);


	protected:

		int inter_index(int i, int j);
		std::vector<std::string> split_string(std::string input);

		arma::vec molecule_charges;
		arma::vec atom_charges;
		arma::vec molecule_charges_sq;
		arma::vec atom_charges_sq; // squared charges

		arma::mat omega_11;
		arma::mat omega_12;
		arma::mat omega_13;
		arma::mat omega_22;
		arma::mat rot_rel_times;
		arma::vec c_rot_arr;
		arma::vec c_rot_rigid_rot_arr;

		arma::mat shear_viscosity_LHS;
		arma::mat thermal_conductivity_LHS;
		arma::mat thermal_conductivity_rigid_rot_LHS;

		arma::mat diffusion_LHS;
		arma::mat diffusion_rigid_rot_LHS;

		arma::mat bulk_viscosity_LHS;
		arma::mat bulk_viscosity_rigid_rot_LHS;

		arma::vec thermal_conductivity_RHS;
		arma::vec thermal_conductivity_rigid_rot_RHS;

		arma::vec diffusion_RHS;
		arma::vec diffusion_rigid_rot_RHS;

		arma::vec shear_viscosity_RHS;
		arma::vec bulk_viscosity_RHS;
		arma::vec bulk_viscosity_rigid_rot_RHS;
		arma::vec thermal_conductivity_rigid_rot_coeffs;
		arma::vec thermal_conductivity_coeffs;

		arma::vec diffusion_coeffs;
		arma::vec diffusion_rigid_rot_coeffs;

		arma::vec shear_viscosity_coeffs;
		arma::vec bulk_viscosity_coeffs;
		arma::vec bulk_viscosity_rigid_rot_coeffs;

		arma::vec empty_n_atom; // if our mixture has no atoms, we pass this as a dummy parameter to functions which expect numeric density of atomic species
		std::vector<arma::vec> empty_n_vl_molecule;

		arma::vec this_n_atom;
		std::vector<arma::vec> this_n_vl_mol;
		arma::vec this_n_molecules;

		double this_total_n;
		double this_total_dens;
		double this_ctr;
		double this_crot;
		double this_n_electrons;

		std::vector<kappa::Molecule> molecules;
		std::vector<kappa::Atom> atoms;
		kappa::Particle electron;

		int num_molecules;
		int num_atoms;
		int n_vibr_levels_total;
		int n_particles;
		bool cache_on;
		bool all_rigid_rotators = true;
		bool is_ionized;

		// this allows to easily calculate which elements of state-to-state matrix we need to fill, given the molecule indices and the vibrational levels
		// vl_offset[0] = 0; 
		// vl_offset[1] = molecules[0].num_vibr_levels[0]; 
		// vl_offset[2] = molecules[0].num_vibr_levels[0] + molecules[1].num_vibr_levels[0], etc.
		std::vector<int> vl_offset;

		std::vector<kappa::Interaction> interactions;
		std::map<std::string, int> molecule_name_map;
		std::map<std::string, int> atom_name_map;

		void init_matrices(const std::string &particles_filename);
		void add_interactions(const std::string &filename);

		void check_n_vl_molecule(const std::vector<arma::vec> &n_vl_molecule); // check sizes of arrays and test values for non-negativity
		void check_n_molecule(const arma::vec &n_molecule); // check sizes of arrays and test values for non-negativity
		void check_n_atom(const arma::vec &n_atom);
		void check_n(const arma::vec &n);

		void compute_omega11(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		void compute_omega12(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		void compute_omega13(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		void compute_omega22(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);

		// Particle + e- interactions
		void compute_omega11(double T, double debye_length, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		void compute_omega12(double T, double debye_length, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		void compute_omega13(double T, double debye_length, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		void compute_omega22(double T, double debye_length, kappa::models_omega model = kappa::models_omega::model_omega_esa);

		void compute_c_rot_rigid_rot(double T);
		void compute_c_rot(double T); // compute c_rot arrays
		void compute_full_crot_rigid_rot(double T);
		void compute_full_crot(double T); // compute c_rot of mixture
		void compute_rot_rel_times(double T, double n, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		void inplace_compute_n_molecule(const std::vector<arma::vec> &n_vl_molecule);

		double th_cond, sh_visc, b_visc;

		arma::vec th_diff;
		arma::mat diff;
		arma::mat lite_diff;
		arma::vec binary_diff;

		const arma::mat &compute_bulk_viscosity_LHS(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		const arma::vec &compute_bulk_viscosity_RHS(double T);
		const arma::vec &compute_bulk_viscosity_coeffs(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);

		const arma::mat &compute_bulk_viscosity_rigid_rot_LHS(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		const arma::vec &compute_bulk_viscosity_rigid_rot_RHS(double T);
		const arma::vec &compute_bulk_viscosity_rigid_rot_coeffs(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);

		const arma::mat &compute_thermal_conductivity_LHS(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		const arma::vec &compute_thermal_conductivity_RHS(double T);
		const arma::vec &compute_thermal_conductivity_coeffs(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);

		const arma::mat &compute_thermal_conductivity_rigid_rot_LHS(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		const arma::vec &compute_thermal_conductivity_rigid_rot_RHS(double T);
		const arma::vec &compute_thermal_conductivity_rigid_rot_coeffs(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);

		// const arma::mat &compute_diffusion_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
		// const arma::vec &compute_diffusion_LHS(double T);
		// const arma::vec &compute_diffusion_RHS(double T);
		// const arma::vec &compute_diffusion_RHS(double T, int b, int n);
		// const arma::vec &compute_diffusion_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
		const arma::mat &compute_diffusion_LHS(double T);
		const arma::mat &lite_compute_diffusion_LHS(double T);
		const arma::mat &lite_compute_diffusion_LHS(double T, int b, int n, std::string ParticleType);
		const arma::vec &compute_diffusion_RHS(double T, int b, int n);
		const arma::vec &lite_compute_diffusion_RHS(double T, int b, int n, std::string ParticleType);
		const arma::vec &compute_diffusion_coeffs(double T, int b, int n);
		const arma::vec &lite_compute_diffusion_coeffs(double T, int b, int n, std::string ParticleType);
		// const arma::mat &compute_diffusion_rigid_rot_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
		// const arma::vec &compute_diffusion_rigid_rot_RHS(double T);
		// const arma::vec &compute_diffusion_rigid_rot_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);

		const arma::mat &compute_shear_viscosity_LHS(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		const arma::vec &compute_shear_viscosity_RHS(double T);
		const arma::vec &compute_shear_viscosity_coeffs(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);

		double bulk_viscosity(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		double thermal_conductivity(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		double shear_viscosity(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);

		void diffusion(double T);
		void lite_diffusion(double T);

		void thermodiffusion(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);
		void binary_diffusion(double T, kappa::models_omega model = kappa::models_omega::model_omega_esa);

	}; // class Mixture
} // namespace kappa
#endif /* kappa_mixture_hpp */
