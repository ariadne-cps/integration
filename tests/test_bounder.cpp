/***************************************************************************
 *            test_bounder.cpp
 *
 *  Copyright  2018-20  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "integration/bounder.hpp"
#include "function/formula.hpp"
#include "function/taylor_model.hpp"

#include "macro/testing.hpp"

using namespace Ariadne;
using namespace std;

class TestBounder
{
  private:
    std::unique_ptr<BounderInterface> bounder_ptr;
    EffectiveScalarMultivariateFunction x, y;
  public:
    TestBounder(const BounderInterface& b)
        : bounder_ptr(b.clone())
    {
        x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        y=EffectiveScalarMultivariateFunction::coordinate(2,1);
    }

    Int test() {
        ARIADNE_TEST_CALL(test_print_name());
        ARIADNE_TEST_CALL(test_D_mismatch());
        ARIADNE_TEST_CALL(test_DA_mismatch());
        ARIADNE_TEST_CALL(test_suggested_step_acceptable());
        ARIADNE_TEST_CALL(test_suggested_step_not_acceptable());
        ARIADNE_TEST_CALL(test_unbounded());
        ARIADNE_TEST_CALL(test_overflowing_vector_field());
        ARIADNE_TEST_CALL(test_time_invariant_with_parameter());
        ARIADNE_TEST_CALL(test_time_variant_without_parameter());
        ARIADNE_TEST_CALL(test_time_variant_with_parameter());
        ARIADNE_TEST_CALL(test_handle_forwarding());
        ARIADNE_TEST_CALL(test_time_variant_dimension_errors());
        return 0;
    }

    Void test_print_name() {
        ARIADNE_TEST_PRINT(*bounder_ptr);
    }

    Void test_D_mismatch() {
        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType D={ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType hsug(0.25_x);

        ARIADNE_TEST_FAIL(bounder_ptr->compute(f,D,suggest(hsug)));
    }

    Void test_DA_mismatch() {
        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType D={ExactIntervalType(-0.25_x,0.25_x),ExactIntervalType(-0.25_x,0.25_x)};
        ExactBoxType A={ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType hsug(0.25_x);

        ARIADNE_TEST_FAIL(bounder_ptr->compute(f,D,A,suggest(hsug)));
    }

    Void test_suggested_step_acceptable() {
        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType dom={ExactIntervalType(-0.25_x,0.25_x),ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType hsug(0.25_x);

        StepSizeType h;
        UpperBoxType B;
        std::tie(h,B) = bounder_ptr->compute(f,dom,suggest(hsug));

        ARIADNE_TEST_PRINT(h);
        ARIADNE_TEST_PRINT(B);
        ARIADNE_TEST_EQUAL(h,hsug);
        ARIADNE_TEST_ASSERT(definitely(is_bounded(B)));
    }

    Void test_suggested_step_not_acceptable() {
        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType dom={ExactIntervalType(-0.25_x,0.25_x),ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType hsug(1.0_x);

        StepSizeType h;
        UpperBoxType B;
        std::tie(h,B) = bounder_ptr->compute(f,dom,suggest(hsug));

        ARIADNE_TEST_PRINT(h);
        ARIADNE_TEST_PRINT(B);
        ARIADNE_TEST_COMPARE(h,<,hsug);
        ARIADNE_TEST_ASSERT(definitely(is_bounded(B)));
    }

    Void test_unbounded() {
        EffectiveVectorMultivariateFunction f={pow(x,3),-pow(y,5)};
        ExactBoxType dom={ExactIntervalType(-20.0_x,20.0_x),ExactIntervalType(-20.0_x,20.0_x)};
        StepSizeType hsug(10.0_x);

        ARIADNE_TEST_THROWS(bounder_ptr->compute(f,dom,suggest(hsug)),BoundingNotFoundException);
    }

    Void test_overflowing_vector_field() {
        // The binary endpoints 2^500 are exactly representable on every
        // supported platform. The Jacobian (3*2^1000) is finite, while the
        // cubic field (2^1500) overflows double-precision bounds. No enclosure
        // can be found above the bounder's minimum admissible step.
        EffectiveVectorMultivariateFunction f={pow(x,3),-y};
        ExactBoxType dom={{-0x1p500_x,0x1p500_x},{-1.0_x,1.0_x}};
        StepSizeType step(0.125_x);
        ARIADNE_TEST_THROWS(bounder_ptr->compute(f,dom,suggest(step)),BoundingNotFoundException);
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

        ARIADNE_TEST_COMPARE(h,<=,hsug);
        ARIADNE_TEST_ASSERT(definitely(is_bounded(B)));
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

        ARIADNE_TEST_COMPARE(h,<=,hsug);
        ARIADNE_TEST_ASSERT(definitely(is_bounded(B)));
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

        ARIADNE_TEST_COMPARE(h,<=,hsug);
        ARIADNE_TEST_ASSERT(definitely(is_bounded(B)));
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
        ARIADNE_TEST_EQUALS(printed.str(),"EulerBounder");

        ExactBoxType D={{0.0_x,0.0_x}};
        ExactBoxType A={{0.0_x,0.25_x}};
        StepSizeType step(0.125_x);
        Suggestion<StepSizeType> suggestion=suggest(step);
        static_cast<StepSizeType&>(suggestion)=step;
        ARIADNE_TEST_EQUALS(static_cast<StepSizeType const&>(suggestion),step);
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
                ARIADNE_TEST_COMPARE(computed.first,<=,step);
                ARIADNE_TEST_ASSERT(definitely(is_bounded(computed.second)));
                ARIADNE_TEST_EQUALS(computed.second.dimension(),1u);
                ARIADNE_TEST_ASSERT(refines(UpperBoxType(D),computed.second));
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
        ARIADNE_TEST_THROWS(bounder.compute(timed,ExactBoxType(0u),step,ExactBoxType(0u),suggest(step)),std::runtime_error);
        ARIADNE_TEST_THROWS(bounder.compute(timed,D,step,D,suggest(step)),std::runtime_error);
    }
};

Int main() {

    List<BounderHandle> bounders = { EulerBounder() };

    for(unsigned verbosity:{0u,3u}) {
        Logger::instance().configuration().set_verbosity(verbosity);
        for (BounderHandle bounder : bounders)
            TestBounder(bounder).test();
    }

    return ARIADNE_TEST_FAILURES;
}
