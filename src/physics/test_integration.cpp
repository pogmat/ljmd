#include "physics.h"
#include "gtest/gtest.h"

class IntegrationTest
    : public ::testing::TestWithParam<std::tuple<double, double>> {

      protected:
        mdsys_t *sys;
        // double force_params = GetParam();

        void SetUp() {
                sys = new mdsys_t;
                sys->natoms = 2;
                sys->mass = 1.0;
			
#if defined(MPI_ENABLED)
                arr_seg_t proc_seg;
                sys->proc_seg = &proc_seg;
#endif
					
                sys->r = new vec3_t[2]();
				sys->v = new vec3_t[2]();
                sys->f = new vec3_t[2]();

			
        }

        void TearDown() {
                delete[] sys->r;

			
				delete[] sys->v;

                delete[] sys->f;

#if defined(MPI_ENABLED)
                sys->proc_seg = nullptr;
#endif

                delete sys;
        }
};

TEST_P(IntegrationTest, testVerlet1) {
        ASSERT_NE(sys, nullptr);


#if defined(MPI_ENABLED)
        sys->proc_seg->size = 2;
        sys->proc_seg->idx = 0;
#endif
	
	sys->dt = std::get<0>(GetParam());
	
        sys->r[0].y = 1.0;
	sys->r[1].y = -1.0;
	
	sys->f[0].x = std::get<1>(GetParam());
	sys->f[1].x = std::get<1>(GetParam());
	
	verlet_1(sys);
	
	double exp_v_increm = 0.5*std::get<0>(GetParam())*std::get<1>(GetParam()) / mvsq2e;
	double exp_r_increm = exp_v_increm*std::get<0>(GetParam());
		
	ASSERT_DOUBLE_EQ(sys->r[0].x,exp_r_increm);
	ASSERT_DOUBLE_EQ(sys->r[0].y,1.0);
	ASSERT_DOUBLE_EQ(sys->r[0].z,0.0);
	ASSERT_DOUBLE_EQ(sys->v[0].x,exp_v_increm);
	ASSERT_DOUBLE_EQ(sys->v[0].y,0.0);
	ASSERT_DOUBLE_EQ(sys->v[0].z,0.0);

	ASSERT_DOUBLE_EQ(sys->r[1].x,exp_r_increm);
	ASSERT_DOUBLE_EQ(sys->r[1].y,-1.0);
	ASSERT_DOUBLE_EQ(sys->r[1].z,0.0);
	ASSERT_DOUBLE_EQ(sys->v[1].x,exp_v_increm);
	ASSERT_DOUBLE_EQ(sys->v[1].y,0.0);
	ASSERT_DOUBLE_EQ(sys->v[1].z,0.0);

}

TEST_P(IntegrationTest, testVerlet2) {
        ASSERT_NE(sys, nullptr);


#if defined(MPI_ENABLED)
        sys->proc_seg->size = 2;
        sys->proc_seg->idx = 0;
#endif
	
	sys->dt = std::get<0>(GetParam());
	
        sys->r[0].y = 1.0;
	sys->r[1].y = -1.0;
	
	sys->f[0].x = std::get<1>(GetParam());
	sys->f[1].x = std::get<1>(GetParam());
	
	verlet_2(sys);
	
	double exp_v_increm = 0.5*std::get<0>(GetParam())*std::get<1>(GetParam()) / mvsq2e;
		
	ASSERT_DOUBLE_EQ(sys->r[0].x,0.0);
	ASSERT_DOUBLE_EQ(sys->r[0].y,1.0);
	ASSERT_DOUBLE_EQ(sys->r[0].z,0.0);
	ASSERT_DOUBLE_EQ(sys->v[0].x,exp_v_increm);
	ASSERT_DOUBLE_EQ(sys->v[0].y,0.0);
	ASSERT_DOUBLE_EQ(sys->v[0].z,0.0);

	ASSERT_DOUBLE_EQ(sys->r[1].x,0.0);
	ASSERT_DOUBLE_EQ(sys->r[1].y,-1.0);
	ASSERT_DOUBLE_EQ(sys->r[1].z,0.0);
	ASSERT_DOUBLE_EQ(sys->v[1].x,exp_v_increm);
	ASSERT_DOUBLE_EQ(sys->v[1].y,0.0);
	ASSERT_DOUBLE_EQ(sys->v[1].z,0.0);

}

INSTANTIATE_TEST_SUITE_P(IntegrationTest_parametric, IntegrationTest,
                         ::testing::Combine(::testing::Values(0.0, 0.5, 1.0),
                                            ::testing::Values(0.0, 1.0, 2.0)));
