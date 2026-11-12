/***************************************************************************
 *            test_bounder.cpp
 *
 *  Copyright  2018-20  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Integration.
 *
 *  Ariadne Integration is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Integration is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Integration. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file tests/test_bounder.cpp
 *  @brief Test bounder behavior.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "integration/bounder.hpp"
#include "function/formula.hpp"
#include "function/taylor_model.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;
using namespace std;

class TestBounder : public Ariadne::Testing::ContextFixture
{
  private:
    std::unique_ptr<BounderInterface> bounder_ptr;
    EffectiveScalarMultivariateFunction x, y;
  public:
    unsigned _saved_verbosity;
  public:
    TestBounder(Ariadne::Testing::Context& context, const BounderInterface& b, unsigned verbosity)
        : ContextFixture(context), bounder_ptr(b.clone()),
          _saved_verbosity(Logger::instance().configuration().verbosity())
    {
        Logger::instance().configuration().set_verbosity(verbosity);
        x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        y=EffectiveScalarMultivariateFunction::coordinate(2,1);
    }



    ~TestBounder() { Logger::instance().configuration().set_verbosity(_saved_verbosity); }

    Void test_print_name() {
        ARIADNE_TESTING_NOTE(*bounder_ptr);
    }

    Void test_D_mismatch() {
        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType D={ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType hsug(0.25_x);

        ARIADNE_TESTING_EXPECT_EXCEPTION(bounder_ptr->compute(f,D,suggest(hsug)));
    }

    Void test_DA_mismatch() {
        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType D={ExactIntervalType(-0.25_x,0.25_x),ExactIntervalType(-0.25_x,0.25_x)};
        ExactBoxType A={ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType hsug(0.25_x);

        ARIADNE_TESTING_EXPECT_EXCEPTION(bounder_ptr->compute(f,D,A,suggest(hsug)));
    }

    Void test_suggested_step_acceptable() {
        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType dom={ExactIntervalType(-0.25_x,0.25_x),ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType hsug(0.25_x);

        StepSizeType h;
        UpperBoxType B;
        std::tie(h,B) = bounder_ptr->compute(f,dom,suggest(hsug));

        ARIADNE_TESTING_NOTE(h);
        ARIADNE_TESTING_NOTE(B);
        ARIADNE_TESTING_EXPECT_EQUAL(h,hsug);
        ARIADNE_TESTING_PROVE(definitely(is_bounded(B)));
    }

    Void test_suggested_step_not_acceptable() {
        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType dom={ExactIntervalType(-0.25_x,0.25_x),ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType hsug(1.0_x);

        StepSizeType h;
        UpperBoxType B;
        std::tie(h,B) = bounder_ptr->compute(f,dom,suggest(hsug));

        ARIADNE_TESTING_NOTE(h);
        ARIADNE_TESTING_NOTE(B);
        ARIADNE_TESTING_EXPECT_LESS(h, hsug);
        ARIADNE_TESTING_PROVE(definitely(is_bounded(B)));
    }

    Void test_unbounded() {
        EffectiveVectorMultivariateFunction f={pow(x,3),-pow(y,5)};
        ExactBoxType dom={ExactIntervalType(-20.0_x,20.0_x),ExactIntervalType(-20.0_x,20.0_x)};
        StepSizeType hsug(10.0_x);

        ARIADNE_TESTING_EXPECT_THROWS(BoundingNotFoundException, bounder_ptr->compute(f,dom,suggest(hsug)));
    }

    Void test_overflowing_vector_field() {
        // The binary endpoints 2^500 are exactly representable on every
        // supported platform. The Jacobian (3*2^1000) is finite, while the
        // cubic field (2^1500) overflows double-precision bounds. No enclosure
        // can be found above the bounder's minimum admissible step.
        EffectiveVectorMultivariateFunction f={pow(x,3),-y};
        ExactBoxType dom={{-0x1p500_x,0x1p500_x},{-1.0_x,1.0_x}};
        StepSizeType step(0.125_x);
        ARIADNE_TESTING_EXPECT_THROWS(BoundingNotFoundException, bounder_ptr->compute(f,dom,suggest(step)));
    }

    Void test_time_invariant_with_parameter() {
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction x1=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction a0=EffectiveScalarMultivariateFunction::coordinate(3,2);

        EffectiveVectorMultivariateFunction f={x0+a0,-x1};
        ExactBoxType D={ExactIntervalType(-0.25_x,0.25_x),ExactIntervalType(-0.25_x,0.25_x)};
        ExactBoxType A={ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType td(1.0_x);
        StepSizeType hsug(1.0_x);

        StepSizeType h;
        UpperBoxType B;

        std::tie(h,B) = bounder_ptr->compute(f,D,A,suggest(hsug));

        ARIADNE_TESTING_EXPECT_LESS_EQUAL(h, hsug);
        ARIADNE_TESTING_PROVE(definitely(is_bounded(B)));
    }

    Void test_time_variant_without_parameter() {
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction x1=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(3,2);

        EffectiveVectorMultivariateFunction f={x0,-x1+2*t};
        ExactBoxType D={ExactIntervalType(-0.25_x,0.25_x),ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType td(1.0_x);
        StepSizeType hsug(1.0_x);

        StepSizeType h;
        UpperBoxType B;

        std::tie(h,B) = bounder_ptr->compute(f,D,td,suggest(hsug));

        ARIADNE_TESTING_EXPECT_LESS_EQUAL(h, hsug);
        ARIADNE_TESTING_PROVE(definitely(is_bounded(B)));
    }

    Void test_time_variant_with_parameter() {
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(4,0);
        EffectiveScalarMultivariateFunction x1=EffectiveScalarMultivariateFunction::coordinate(4,1);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(4,2);
        EffectiveScalarMultivariateFunction a0=EffectiveScalarMultivariateFunction::coordinate(4,3);

        EffectiveVectorMultivariateFunction f={x0+a0,-x1+2*t};
        ExactBoxType D={ExactIntervalType(-0.25_x,0.25_x),ExactIntervalType(-0.25_x,0.25_x)};
        ExactBoxType A={ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType td(1.0_x);
        StepSizeType hsug(1.0_x);

        StepSizeType h;
        UpperBoxType B;

        std::tie(h,B) = bounder_ptr->compute(f,D,td,A,suggest(hsug));

        ARIADNE_TESTING_EXPECT_LESS_EQUAL(h, hsug);
        ARIADNE_TESTING_PROVE(definitely(is_bounded(B)));
    }

    Void test_handle_forwarding() {
        BounderHandle original{EulerBounder()};
        BounderHandle assigned{EulerBounder()};
        assigned=original;
        BounderHandle const& self=assigned;
        assigned=self;
        std::unique_ptr<BounderHandle> cloned(assigned.clone());
        std::ostringstream printed;
        cloned->_write(printed);
        ARIADNE_TESTING_EXPECT_EQUAL(printed.str(),"EulerBounder");

        ExactBoxType D={{0.0_x,0.0_x}};
        ExactBoxType A={{0.0_x,0.25_x}};
        StepSizeType step(0.125_x);
        Suggestion<StepSizeType> suggestion=suggest(step);
        static_cast<StepSizeType&>(suggestion)=step;
        ARIADNE_TESTING_EXPECT_EQUAL(static_cast<StepSizeType const&>(suggestion),step);
        auto x0=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
        auto xp=EffectiveScalarMultivariateFunction::coordinate(2u,0u);
        auto p=EffectiveScalarMultivariateFunction::coordinate(2u,1u);
        auto xt=EffectiveScalarMultivariateFunction::coordinate(2u,0u);
        auto t=EffectiveScalarMultivariateFunction::coordinate(2u,1u);
        auto xtp=EffectiveScalarMultivariateFunction::coordinate(3u,0u);
        auto tp=EffectiveScalarMultivariateFunction::coordinate(3u,1u);
        auto ap=EffectiveScalarMultivariateFunction::coordinate(3u,2u);

        // Both logging configurations must preserve the computed enclosure.
        auto previous_verbosity=Logger::instance().configuration().verbosity();
        for(unsigned verbosity:{0u,3u}) {
            Logger::instance().configuration().set_verbosity(verbosity);
            auto autonomous=cloned->compute(EffectiveVectorMultivariateFunction({x0}),D,suggestion);
            auto parametric=cloned->compute(EffectiveVectorMultivariateFunction({xp+p}),D,A,suggestion);
            auto timed=cloned->compute(EffectiveVectorMultivariateFunction({xt+t}),D,step,suggestion);
            auto full=cloned->compute(EffectiveVectorMultivariateFunction({xtp+tp+ap}),D,step,A,suggestion);
            for(auto const& computed:{autonomous,parametric,timed,full}) {
                ARIADNE_TESTING_EXPECT_LESS_EQUAL(computed.first, step);
                ARIADNE_TESTING_PROVE(definitely(is_bounded(computed.second)));
                ARIADNE_TESTING_EXPECT_EQUAL(computed.second.dimension(),1u);
                ARIADNE_TESTING_PROVE(refines(UpperBoxType(D),computed.second));
            }
        }
        Logger::instance().configuration().set_verbosity(previous_verbosity);
    }

    Void test_time_variant_dimension_errors() {
        EulerBounder bounder;
        StepSizeType step(0.125_x);
        ExactBoxType D={{0.0_x,1.0_x}};
        auto x0=EffectiveScalarMultivariateFunction::coordinate(2u,0u);
        EffectiveVectorMultivariateFunction timed={x0};
        ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, bounder.compute(timed,ExactBoxType(0u),step,ExactBoxType(0u),suggest(step)));
        ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, bounder.compute(timed,D,step,D,suggest(step)));
    }
};


int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_bounder");
    auto bounder_v0 = Ariadne::Testing::fixture_variant<TestBounder>(ariadne_test_suite, "EulerBounder/verbosity=0", [](Ariadne::Testing::Context& ctx) { return TestBounder(ctx, EulerBounder(), 0u); });
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v0, test_print_name);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v0, test_D_mismatch);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v0, test_DA_mismatch);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v0, test_suggested_step_acceptable);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v0, test_suggested_step_not_acceptable);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v0, test_unbounded);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v0, test_overflowing_vector_field);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v0, test_time_invariant_with_parameter);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v0, test_time_variant_without_parameter);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v0, test_time_variant_with_parameter);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v0, test_handle_forwarding);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v0, test_time_variant_dimension_errors);
    auto bounder_v3 = Ariadne::Testing::fixture_variant<TestBounder>(ariadne_test_suite, "EulerBounder/verbosity=3", [](Ariadne::Testing::Context& ctx) { return TestBounder(ctx, EulerBounder(), 3u); });
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v3, test_print_name);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v3, test_D_mismatch);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v3, test_DA_mismatch);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v3, test_suggested_step_acceptable);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v3, test_suggested_step_not_acceptable);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v3, test_unbounded);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v3, test_overflowing_vector_field);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v3, test_time_invariant_with_parameter);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v3, test_time_variant_without_parameter);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v3, test_time_variant_with_parameter);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v3, test_handle_forwarding);
    ARIADNE_TESTING_ADD_VARIANT_TEST(bounder_v3, test_time_variant_dimension_errors);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
