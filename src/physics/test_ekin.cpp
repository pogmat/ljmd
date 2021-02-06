#include "gtest/gtest.h"


#include "physics.h"





class Ekin_T_Test : public ::testing::Test {

protected:
        mdsys_t *sys;

        void SetUp() {
                sys = new mdsys_t;
                sys->natoms = 2;
                sys->mass = 1.0;
                sys->v = new vec3_t[2];

                sys->v[0].x = 0.0;
                sys->v[1].x = 0.0;
		sys->v[0].y = 1.0;
                sys->v[1].y = -1.0;
		sys->v[0].z = 0.0;
                sys->v[1].z = 0.0;

        }

        void TearDown() {
                delete[] sys->v;
			
                delete sys;
        }
};

TEST(Ekin_T_TestEmpty, empty) {
        mdsys_t *sys = new mdsys_t;
        sys->natoms = 0;

        ekin(sys);
	
        ASSERT_DOUBLE_EQ(sys->ekin, 0.0);
        ASSERT_DOUBLE_EQ(sys->temp, 0.0);

        delete sys;
}



TEST_F(Ekin_T_Test, test1) {
        ASSERT_NE(sys, nullptr);
        ASSERT_DOUBLE_EQ(sys->v[0].y, 1.0);
        ASSERT_DOUBLE_EQ(sys->v[1].y, -1.0);
	
        ekin(sys);
		double exp_ekin = mvsq2e;
		double exp_temp = 2*exp_ekin / (3* kboltz );
	
        ASSERT_DOUBLE_EQ(sys->ekin, exp_ekin);
        ASSERT_DOUBLE_EQ(sys->temp, exp_temp);
}
