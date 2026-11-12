/***************************************************************************
 *            benchmark_integrators_single_step.cpp
 *
 *  Copyright  2026  Luca Geretti
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

/** @file benchmarks/benchmark_integrators_single_step.cpp
 *  @brief Benchmark benchmark integrators single step performance.
 */

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <string>
#include <utility>
#include <tuple>
#include <vector>

#include "ariadne-benchmarking.hpp"
#include "integrator_benchmark.hpp"

using namespace Ariadne;
using namespace IntegratorBenchmark;
namespace Bench=Ariadne::Benchmarking;

namespace {

constexpr double minimum_total_measurement_ns=1000000.0;
constexpr Nat maximum_repetitions=1024u;
constexpr double benchmark_step_maximum_error=1e-5;

struct Measurement {
    double elapsed_ns;
    FlowStepModelType flow;
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

struct BenchmarkResult {
    std::string integrator;
    std::string function;
    Nat repetitions;
    double time_us_per_step;
    Bench::Statistics timing;
    decltype(std::declval<FlowStepModelType>().error()) remainder;
    decltype(std::declval<FlowStepModelType>().errors()) component_remainders;
};

StepSizeType benchmark_step() {
    return 0.03125_x;
}

Measurement measure(
        const IntegratorInterface& integrator,
        const BenchmarkCase& benchmark_case,
        const StepSizeType& step,
        Nat repetitions)
{
    // The original pilot does one untimed invocation before measurement.
    FlowStepModelType flow=run_integrator_step(
        integrator,benchmark_case,step);
    const double elapsed_ns=Bench::measure_ns([&] {
        for(Nat i=0u;i!=repetitions;++i)
            flow=run_integrator_step(integrator,benchmark_case,step);
    });
    return {elapsed_ns,std::move(flow)};
}

Nat determine_repetitions(
        const std::vector<NamedIntegrator>& integrators,
        const std::vector<BenchmarkCase>& cases,
        const StepSizeType& step)
{
    double minimum_elapsed_ns=std::numeric_limits<double>::infinity();
    LOGGING_PRINTLN_AT(2,
        "Calibrating repetitions with step="<<step
        <<", minimum measurement time="<<minimum_total_measurement_ns/1e6<<" ms")

    for(const auto& named_integrator:integrators) {
        for(const auto& benchmark_case:cases) {
            LOGGING_PRINTLN_AT(3,
                "Calibration start: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name)
            const auto measurement=measure(
                *named_integrator.integrator,benchmark_case,step,1u);
            LOGGING_PRINTLN_AT(2,
                "Calibration: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", elapsed="<<measurement.elapsed_ns/1e6<<" ms")
            LOGGING_PRINTLN_AT(3,
                "Calibration remainder: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", remainder="<<measurement.flow.error()
                <<", components="<<measurement.flow.errors())
            minimum_elapsed_ns=std::min(minimum_elapsed_ns,measurement.elapsed_ns);
        }
    }

    if(minimum_elapsed_ns<=0.0) {
        LOGGING_PRINTLN_AT(2,
            "Calibration selected maximum repetitions="<<maximum_repetitions
            <<" because the minimum measured duration was non-positive")
        return maximum_repetitions;
    }
    const Nat repetitions=Bench::calibrated_repetitions(
        minimum_elapsed_ns,minimum_total_measurement_ns,1u,maximum_repetitions);
    LOGGING_PRINTLN_AT(2,"Calibration selected repetitions="<<repetitions)
    return repetitions;
}

std::vector<BenchmarkResult> run_benchmarks(
        const std::vector<NamedIntegrator>& integrators,
        const std::vector<BenchmarkCase>& cases,
        const StepSizeType& step,
        Nat repetitions,
        Bench::Runner const& runner,
        unsigned samples)
{
    std::vector<BenchmarkResult> results;
    results.reserve(integrators.size()*cases.size());
    LOGGING_PRINTLN_AT(2,
        "Running benchmark with step="<<step
        <<", step maximum error="<<benchmark_step_maximum_error
        <<", repetitions="<<repetitions
        <<", samples="<<samples)

    for(const auto& named_integrator:integrators) {
        LOGGING_PRINTLN_AT(1,"Integrator="<<named_integrator.key)
        for(const auto& benchmark_case:cases) {
            LOGGING_PRINTLN_AT(3,
                "Benchmark start: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name)

            // Ensure the flow object exists before the timed region,
            // matching the one-time untimed setup in the old benchmark.
            FlowStepModelType flow=run_integrator_step(
                *named_integrator.integrator,benchmark_case,step);
            auto const times=runner.run(
                named_integrator.key+"/"+benchmark_case.name,samples,
                [] {},
                [&] {
                    for(Nat i=0u;i!=repetitions;++i)
                        flow=run_integrator_step(
                            *named_integrator.integrator,benchmark_case,step);
                },
                [&] { (void)flow.error(); },
                static_cast<double>(repetitions));
            auto const statistics=Bench::summarize(times);
            double const time_us_per_step=statistics.median/1000.0;
            if(not std::isfinite(time_us_per_step) or time_us_per_step<=0.0)
                throw std::runtime_error("Invalid benchmark timing.");

            LOGGING_PRINTLN_AT(2,
                "Benchmark result: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", median="<<time_us_per_step<<" us/step"
                <<", remainder="<<flow.error())
            LOGGING_PRINTLN_AT(3,
                "Benchmark component remainders: integrator="
                <<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", components="<<flow.errors())
            results.push_back({
                named_integrator.key,benchmark_case.name,repetitions,
                time_us_per_step,statistics,flow.error(),flow.errors()
            });
        }
    }
    return results;
}

void print_results(
        const StepSizeType& step,
        Nat repetitions,
        const std::vector<BenchmarkResult>& results,
        Bench::SummaryFields const& columns)
{
    // Preserve both views and the diagnostic remainder fields. The timed
    // statistics are normalized nanoseconds per integrator step.
    std::ostringstream header;
    header<<"step="<<step
          <<", step_maximum_error="<<benchmark_step_maximum_error
          <<", repetitions="<<repetitions
          <<"; timings in ns per step";
    const auto make_report=[&](std::string title,bool by_integrator) {
        Bench::Report report(std::move(title),columns);
        report.set_parameter_header(header.str());
        report.set_dimension_columns({
            {by_integrator?"Integrator":"Function"},
            {by_integrator?"Function":"Integrator"},
            {"Median (us)",0u,Bench::TableAlignment::Right},
            {"Remainder"},
            {"Component remainders"}
        });
        auto opts=Bench::TableOptions::ruled();
        opts.group_columns={0u};
        report.set_table_options(opts);
        // Remainder values need not be assignable. Sort pointers to the
        // immutable measurements instead of copying or moving results.
        std::vector<BenchmarkResult const*> sorted;
        sorted.reserve(results.size());
        for(auto const& result:results) sorted.push_back(&result);
        std::sort(sorted.begin(),sorted.end(),
            [by_integrator](auto const* left,auto const* right) {
                if(by_integrator)
                    return std::tie(left->integrator,left->function)
                         < std::tie(right->integrator,right->function);
                return std::tie(left->function,left->integrator)
                     < std::tie(right->function,right->integrator);
            });
        for(auto const* result:sorted) {
            std::ostringstream rem,components;
            rem<<result->remainder;
            components<<result->component_remainders;
            report.add_row({
                by_integrator?result->integrator:result->function,
                by_integrator?result->function:result->integrator,
                Bench::Table::number(result->time_us_per_step,6u),
                rem.str(),components.str()
            },result->timing);
        }
        report.print();
    };
    make_report("Single-step integrator benchmark: by integrator",true);
    make_report("Single-step integrator benchmark: by function",false);
}

} // namespace

Int main(Int argc, const char* argv[]) {
    BenchmarkCommandLineSelection selection;
    const auto benchmark_parsers=
        make_benchmark_argument_parsers(selection);

    Bench::CommandLine command_line;
    command_line.add(CommandLineContribution(
        "integrator-benchmark",benchmark_parsers));
    unsigned samples=1u;
    command_line.add(CommandLineContribution("integrator-samples",{
        make_valued_argument_parser("","samples",
            "Measured rounds per integrator/case (default: 1)",
            [&samples](String const& value)->VoidFunction {
                unsigned const count=Bench::parse_unsigned<unsigned>(value);
                if(count==0u) throw std::invalid_argument("samples must be positive");
                return [&samples,count] { samples=count; };
            })
    }));

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
        <<", benchmark cases="<<cases.size())

    const auto integrators=select_integrators(
        selection.integrator,benchmark_step_maximum_error);
    const StepSizeType step=benchmark_step();
    const Nat repetitions=determine_repetitions(integrators,cases,step);
    Bench::Runner runner(command_line,integrators.size()*cases.size());
    const auto results=run_benchmarks(
        integrators,cases,step,repetitions,runner,samples);

    print_results(step,repetitions,results,command_line.summary_fields());
    return 0;
}
