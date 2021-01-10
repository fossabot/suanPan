#include <catch/catch.hpp>
#include <Toolbox/Quaternion.hpp>

TEST_CASE("Quaternion Basic Function", "[Utility.Quaternion]") {
	const Quaternion<double> A(2., 3., 4., 5.);
	const Quaternion<double> B(1., -2., 6., 3.);

	const auto C = A + B;

	REQUIRE(Approx(3) == C.real());

	const auto D = A * B;

	REQUIRE(Approx(-31) == D.real());
	REQUIRE(Approx(-19) == D.imag()(0));
	REQUIRE(Approx(-3) == D.imag()(1));
	REQUIRE(Approx(37) == D.imag()(2));

	const auto E = B * A;

	REQUIRE(Approx(-31) == E.real());
	REQUIRE(Approx(17) == E.imag()(0));
	REQUIRE(Approx(35) == E.imag()(1));
	REQUIRE(Approx(-15) == E.imag()(2));
}
