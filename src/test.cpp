#include <hpx/hpx_init.hpp>
#include <immintrin.h>
#include <tigersph/simd.hpp>

void advance_cosmos(float &a, float &adot, float dt) {
	const float c0 = 0.6;
	const auto dadt = [=](float adot) {
		return adot * dt;
	};
	const auto dadotdt = [dt, c0](float a) {
		return -c0 * dt / (a * a);
	};
	const float da1 = dadt(adot);
	const float dadot1 = dadotdt(a);
	const float da2 = dadt(adot + 0.5 * dadot1);
	const float dadot2 = dadotdt(a + 0.5 * da1);
	const float da3 = dadt(adot + 0.5 * dadot2);
	const float dadot3 = dadotdt(a + 0.5 * da2);
	const float da4 = dadt(adot + dadot3);
	const float dadot4 = dadotdt(a + da3);
	a += (da1 + 2.0 * (da2 + da3) + da4) / 6.0;
	adot += (dadot1 + 2.0 * (dadot2 + dadot3) + dadot4) / 6.0;

}

int hpx_main(int argc, char *argv[]) {

	float a = 1.0;
	float adot = 1.0;
	float t = 0.0;
	float dt = 0.01;
	for (int i = 0; t < 100.0; i++) {
		printf("%e %e %e\n", t, a, adot);
		advance_cosmos(a, adot, dt);
		if (a <= 0.0) {
			break;
		}
		t += dt;
	}
	return hpx::finalize();
}
int main(int argc, char *argv[]) {
	std::vector < std::string > cfg = { "hpx.commandline.allow_unknown=1" };
	hpx::init(argc, argv, cfg);
}

