#include <algorithm>
#include <chrono>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "drawing/command_line_interface.hpp"
#include "logging/logging.hpp"
#include "integrator_benchmark.hpp"

using namespace Ariadne;
using namespace IntegratorBenchmark;

namespace {

using Clock=std::chrono::steady_clock;

constexpr double benchmark_step_maximum_error=1e-8;
constexpr Nat maximum_halvings=20u;
constexpr Nat refinement_iterations=6u;

struct Trial {
    Bool accepted;
    double elapsed_ms;
    std::optional<FlowStepModelType> flow;
};

struct MaximumStepResult {
    std::string display;
};

FlowStepModelType run_integrator_step(
        const IntegratorInterface& integrator,
        const BenchmarkCase& benchmark_case,
        const StepSizeType& step)
{
    Logger::instance().increase_level(3);
    try {
        FlowStepModelType flow=integrator.flow_step(
            benchmark_case.vector_field,
            benchmark_case.domain,
            step);
        Logger::instance().decrease_level(3);
        return flow;
    } catch(...) {
        Logger::instance().decrease_level(3);
        throw;
    }
}

Trial try_step(
        const NamedIntegrator& named_integrator,
        const BenchmarkCase& benchmark_case,
        const StepSizeType& step)
{
    LOGGING_PRINTLN_AT(3,
        "Trial start: integrator="<<named_integrator.key
        <<", function="<<benchmark_case.name
        <<", step="<<step)

    const auto start=Clock::now();
    try {
        auto flow=run_integrator_step(
            *named_integrator.integrator,
            benchmark_case,
            step);

        const double elapsed_ms=
            std::chrono::duration<double,std::milli>(
                Clock::now()-start).count();
        const auto remainder=flow.error();

        LOGGING_PRINTLN_AT(3,
            "Trial accepted: integrator="<<named_integrator.key
            <<", function="<<benchmark_case.name
            <<", step="<<step
            <<", elapsed="<<elapsed_ms<<" ms"
            <<", remainder="<<remainder)

        return {true,elapsed_ms,std::move(flow)};
    } catch(const FlowTimeStepException& exception) {
        LOGGING_PRINTLN_AT(3,
            "Trial rejected: integrator="<<named_integrator.key
            <<", function="<<benchmark_case.name
            <<", step="<<step
            <<", reason=FlowTimeStepException"
            <<", message="<<exception.what())
    } catch(const IncompleteFlowException& exception) {
        LOGGING_PRINTLN_AT(3,
            "Trial rejected: integrator="<<named_integrator.key
            <<", function="<<benchmark_case.name
            <<", step="<<step
            <<", reason=IncompleteFlowException"
            <<", message="<<exception.what())
    } catch(const FlowBoundsException& exception) {
        LOGGING_PRINTLN_AT(3,
            "Trial rejected: integrator="<<named_integrator.key
            <<", function="<<benchmark_case.name
            <<", step="<<step
            <<", reason=FlowBoundsException"
            <<", message="<<exception.what())
    }

    const double elapsed_ms=
        std::chrono::duration<double,std::milli>(
            Clock::now()-start).count();

    return {false,elapsed_ms,std::nullopt};
}

MaximumStepResult benchmark_pair(
        const NamedIntegrator& named_integrator,
        const BenchmarkCase& benchmark_case)
{
    LOGGING_PRINTLN_AT(3,
        "Maximum-step search: integrator="<<named_integrator.key
        <<", function="<<benchmark_case.name)

    StepSizeType lower=0.0_x;
    StepSizeType upper=0.03125_x;
    Trial upper_trial=try_step(named_integrator,benchmark_case,upper);

    if(upper_trial.accepted) {
        lower=upper;
        for(Nat i=0u; i!=maximum_halvings; ++i) {
            StepSizeType candidate=upper+upper;
            Trial trial=try_step(
                named_integrator,benchmark_case,candidate);
            if(!trial.accepted) {
                upper=candidate;
                upper_trial=std::move(trial);
                break;
            }
            lower=candidate;
            upper=candidate;
            upper_trial=std::move(trial);
        }

        if(upper_trial.accepted) {
            std::ostringstream stream;
            stream << ">=" << lower;

            LOGGING_PRINTLN_AT(2,
                "Maximum-step search completed: integrator="
                <<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", result="<<stream.str())

            return {stream.str()};
        }
    } else {
        for(Nat i=0u; i!=maximum_halvings; ++i) {
            upper=hlf(upper);
            Trial trial=try_step(
                named_integrator,benchmark_case,upper);
            if(trial.accepted) {
                lower=upper;
                upper=upper+upper;
                break;
            }
        }

        if(lower==0.0_x) {
            std::ostringstream stream;
            stream << "<" << upper;

            LOGGING_PRINTLN_AT(2,
                "Maximum-step search completed: integrator="
                <<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", result="<<stream.str())

            return {stream.str()};
        }
    }

    LOGGING_PRINTLN_AT(3,
        "Maximum-step bracket: integrator="<<named_integrator.key
        <<", function="<<benchmark_case.name
        <<", lower="<<lower
        <<", upper="<<upper)

    for(Nat i=0u; i!=refinement_iterations; ++i) {
        const StepSizeType middle=hlf(lower+upper);
        Trial trial=try_step(
            named_integrator,benchmark_case,middle);
        if(trial.accepted) {
            lower=middle;
        } else {
            upper=middle;
        }
    }

    LOGGING_PRINTLN_AT(2,
        "Maximum-step search completed: integrator="
        <<named_integrator.key
        <<", function="<<benchmark_case.name
        <<", accepted="<<lower
        <<", rejected="<<upper)

    std::ostringstream stream;
    stream << lower;
    return {stream.str()};
}

void print_table(
        const std::vector<NamedIntegrator>& integrators,
        const std::vector<BenchmarkCase>& cases,
        const std::vector<std::vector<MaximumStepResult>>& results)
{
    constexpr std::size_t column_spacing=2u;

    std::size_t function_width=std::string("function").size();
    for(const auto& benchmark_case : cases) {
        function_width=std::max(
            function_width,
            benchmark_case.name.size());
    }

    std::vector<std::size_t> integrator_widths(integrators.size());
    for(SizeType integrator_index=0u;
            integrator_index!=integrators.size();
            ++integrator_index) {
        auto& width=integrator_widths[integrator_index];
        width=integrators[integrator_index].key.size();
        for(SizeType case_index=0u;
                case_index!=cases.size();
                ++case_index) {
            width=std::max(
                width,
                results[case_index][integrator_index].display.size());
        }
    }

    std::ostringstream table;
    table
        <<"Maximum-step integrator benchmark\n"
        <<"step_maximum_error="<<benchmark_step_maximum_error<<"\n"
        <<"Values are largest accepted steps found; "
        <<">= means no failing upper bound was found; "
        <<"< means no accepted step was found\n";

    table
        <<std::left
        <<std::setw(static_cast<int>(function_width))
        <<"function";

    for(SizeType integrator_index=0u;
            integrator_index!=integrators.size();
            ++integrator_index) {
        table
            <<std::string(column_spacing,' ')
            <<std::left
            <<std::setw(static_cast<int>(
                integrator_widths[integrator_index]))
            <<integrators[integrator_index].key;
    }
    table<<"\n";

    for(SizeType case_index=0u; case_index!=cases.size(); ++case_index) {
        table
            <<std::left
            <<std::setw(static_cast<int>(function_width))
            <<cases[case_index].name;

        for(SizeType integrator_index=0u;
                integrator_index!=integrators.size();
                ++integrator_index) {
            table
                <<std::string(column_spacing,' ')
                <<std::left
                <<std::setw(static_cast<int>(
                    integrator_widths[integrator_index]))
                <<results[case_index][integrator_index].display;
        }
        table<<"\n";
    }

    LOGGING_PRINTLN(table.str())
}

} // namespace

Int main(Int argc, const char* argv[]) {
    Logger::instance().configuration().set_verbosity(1u);

    BenchmarkCommandLineSelection selection;
    const auto benchmark_parsers=
        make_benchmark_argument_parsers(selection);

    auto command_line=make_kernel_command_line_interface();
    command_line.add(CommandLineContribution(
        "integrator-benchmark",
        benchmark_parsers));

    if(not command_line.acquire(argc,argv)) {
        return -1;
    }

    const auto cases=make_cases(selection.benchmarks);
    if(selection.list_cases) {
        print_benchmark_cases(cases);
        return 0;
    }
    if(cases.empty()) {
        throw std::runtime_error(
            "Benchmark selector matched no cases.");
    }

    LOGGING_PRINTLN_AT(2,
        "Integrator selection="<<selection.integrator
        <<", step maximum error="<<benchmark_step_maximum_error
        <<", benchmark cases="<<cases.size())

    const auto integrators=select_integrators(
        selection.integrator,benchmark_step_maximum_error);

    std::vector<std::vector<MaximumStepResult>> results(
        cases.size(),
        std::vector<MaximumStepResult>(integrators.size()));

    for(SizeType integrator_index=0u;
            integrator_index!=integrators.size();
            ++integrator_index) {
        LOGGING_PRINTLN_AT(1,
            "Integrator="<<integrators[integrator_index].key)

        for(SizeType case_index=0u;
                case_index!=cases.size();
                ++case_index) {
            results[case_index][integrator_index]=benchmark_pair(
                integrators[integrator_index],
                cases[case_index]);
        }
    }

    print_table(integrators,cases,results);
    return 0;
}
