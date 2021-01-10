#include <catch/catch.hpp>
#include <Toolbox/utility.h>

TEST_CASE("Binomial Compute Basic Function", "[Utility.Binomial]") {
	REQUIRE(165 == suanpan::binomial(11, 3));
	REQUIRE(6435 == suanpan::binomial(15, 7));
	REQUIRE(455 == suanpan::binomial(15, 3));
	REQUIRE(126 == suanpan::binomial(9, 4));
	REQUIRE(220 == suanpan::binomial(12, 3));
}

TEST_CASE("Sign", "[Utility.Sign]") {
	REQUIRE(Approx(1) == suanpan::sign(std::numeric_limits<double>::epsilon()));
	REQUIRE(Approx(-1) == suanpan::sign(-std::numeric_limits<double>::epsilon()));
	REQUIRE(Approx(0) == suanpan::sign(0));
}
