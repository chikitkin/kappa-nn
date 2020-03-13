//
//  models.h
//  String values of models
//

#ifndef kappa_models_h
#define kappa_models_h

namespace kappa {

	enum class models_prob_vv {model_prob_vv_fho};
	enum class models_prob_vt {model_prob_vt_fho};
	enum class models_prob_diss {model_prob_diss_thresh_cmass_vibr, model_prob_diss_thresh_vibr, model_prob_diss_thresh_cmass, model_prob_diss_thresh};
	enum class models_cs_elastic {model_cs_el_rs, model_cs_el_vss};
	enum class models_cs_vv {model_cs_vv_rs_fho, model_cs_vv_vss_fho};
	enum class models_cs_vt {model_cs_vt_rs_fho, model_cs_vt_vss_fho};
	enum class models_cs_diss {model_cs_diss_rs_thresh_cmass_vibr, model_cs_diss_rs_thresh_vibr, model_cs_diss_rs_thresh_cmass, model_cs_diss_rs_thresh,
					       model_cs_diss_vss_thresh_cmass_vibr, model_cs_diss_vss_thresh_vibr, model_cs_diss_vss_thresh_cmass, model_cs_diss_vss_thresh, model_cs_diss_ilt};
	enum class models_k_vv {model_k_vv_rs_fho, model_k_vv_vss_fho, model_k_vv_ssh, model_k_vv_billing};
	enum class models_k_vt {model_k_vt_rs_fho, model_k_vt_vss_fho, model_k_vt_ssh, model_k_vt_phys4entry, model_k_vt_billing};
    enum class models_k_exch {model_k_exch_arrh_scanlon, model_k_exch_arrh_park, model_k_exch_warnatz, model_k_exch_rf, model_k_exch_polak,
    						  model_k_exch_maliat_D6k_arrh_scanlon, model_k_exch_maliat_3T_arrh_scanlon, model_k_exch_maliat_infty_arrh_scanlon,
    						  model_k_exch_maliat_D6k_arrh_park, model_k_exch_maliat_3T_arrh_park, model_k_exch_maliat_infty_arrh_park};
	enum class models_k_diss {model_k_diss_rs_thresh_cmass_vibr, model_k_diss_rs_thresh_vibr, model_k_diss_rs_thresh_cmass, model_k_diss_rs_thresh,
					          model_k_diss_vss_thresh_cmass_vibr, model_k_diss_vss_thresh_vibr, model_k_diss_vss_thresh_cmass, model_k_diss_vss_thresh,
					          model_k_diss_arrh_scanlon, model_k_diss_arrh_park,
						      model_k_diss_tm_D6k_arrh_scanlon,  model_k_diss_tm_3T_arrh_scanlon, model_k_diss_tm_infty_arrh_scanlon,
						      model_k_diss_tm_D6k_arrh_park,  model_k_diss_tm_3T_arrh_park, model_k_diss_tm_infty_arrh_park,
					          model_k_diss_phys4entry, model_k_diss_ilt};
    enum class models_omega {model_omega_rs, model_omega_vss, model_omega_bornmayer, model_omega_lennardjones, model_omega_esa};
}

#endif /* models_h */