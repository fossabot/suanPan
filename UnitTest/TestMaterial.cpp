#include <catch/catch.hpp>
#include <Material/Material1D/Elastic/Elastic1D.h>

TEST_CASE("Elastic1D Basic Function", "[Material.Elastic1D]") {
	unique_ptr<Material> mat_obj = make_unique<Elastic1D>(0, 2E5);

	mat_obj->Material::initialize();
	mat_obj->initialize();

	for(auto I = 0; I < 10; ++I) {
		mat_obj->update_incre_status(1.2E-4);
		mat_obj->commit_status();
	}

	REQUIRE(Approx(240) == mat_obj->get_trial_stress()(0));
}
