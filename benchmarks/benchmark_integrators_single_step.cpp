#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "drawing/command_line_interface.hpp"
#include "logging/logging.hpp"
#include "integrator_benchmark.hpp"

using namespace Ariadne;
using namespace IntegratorBenchmark;

namespace {

using Clock = std::chrono::steady_clock;
using Duration = std::chrono::nanoseconds;

constexpr auto minimum_total_measurement_time=std::chrono::milliseconds(1);
constexpr Nat maximum_repetitions=1024u;
constexpr double benchmark_step_maximum_error=1e-5;

struct Measurement {
    Duration elapsed;
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
    FlowStepModelType flow=run_integrator_step(
        integrator,
        benchmark_case,
        step);

    const auto start=Clock::now();
    for(Nat i=0u; i!=repetitions; ++i) {
        flow=run_integrator_step(
            integrator,
            benchmark_case,
            step);
    }
    const auto stop=Clock::now();

    return {
        std::chrono::duration_cast<Duration>(stop-start),
        std::move(flow)
    };
}

Nat determine_repetitions(
        const std::vector<NamedIntegrator>& integrators,
        const std::vector<BenchmarkCase>& cases,
        const StepSizeType& step)
{
    Duration minimum_elapsed=Duration::max();
    const double minimum_measurement_time_ms=
        std::chrono::duration<double,std::milli>(
            minimum_total_measurement_time).count();

    LOGGING_PRINTLN_AT(2,
        "Calibrating repetitions with step="<<step
        <<", minimum measurement time="
        <<minimum_measurement_time_ms
        <<" ms")

    for(const auto& named_integrator : integrators) {
        for(const auto& benchmark_case : cases) {
            LOGGING_PRINTLN_AT(3,
                "Calibration start: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name)

            const auto measurement=measure(
                *named_integrator.integrator,
                benchmark_case,
                step,
                1u);

            const double elapsed_ms=
                std::chrono::duration<double,std::milli>(
                    measurement.elapsed).count();

            LOGGING_PRINTLN_AT(2,
                "Calibration: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", elapsed="<<elapsed_ms<<" ms")

            LOGGING_PRINTLN_AT(3,
                "Calibration remainder: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", remainder="<<measurement.flow.error()
                <<", components="<<measurement.flow.errors())

            if(measurement.elapsed<minimum_elapsed) {
                minimum_elapsed=measurement.elapsed;
            }
        }
    }

    if(minimum_elapsed<=Duration::zero()) {
        LOGGING_PRINTLN_AT(2,
            "Calibration selected maximum repetitions="
            <<maximum_repetitions
            <<" because the minimum measured duration was non-positive")
        return maximum_repetitions;
    }

    const double required=
        std::chrono::duration<double>(minimum_total_measurement_time).count()
        / std::chrono::duration<double>(minimum_elapsed).count();

    Nat repetitions=static_cast<Nat>(std::ceil(required));
    if(repetitions<1u) {
        repetitions=1u;
    }
    if(repetitions>maximum_repetitions) {
        repetitions=maximum_repetitions;
    }

    LOGGING_PRINTLN_AT(2,
        "Calibration selected repetitions="<<repetitions)

    return repetitions;
}

std::vector<BenchmarkResult> run_benchmarks(
        const std::vector<NamedIntegrator>& integrators,
        const std::vector<BenchmarkCase>& cases,
        const StepSizeType& step,
        Nat repetitions)
{
    std::vector<BenchmarkResult> results;
    results.reserve(integrators.size()*cases.size());

    LOGGING_PRINTLN_AT(2,
        "Running benchmark with step="<<step
        <<", step maximum error="<<benchmark_step_maximum_error
        <<", repetitions="<<repetitions)

    for(const auto& named_integrator : integrators) {
        LOGGING_PRINTLN_AT(1,
            "Integrator="<<named_integrator.key)

        for(const auto& benchmark_case : cases) {
            LOGGING_PRINTLN_AT(3,
                "Benchmark start: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name)

            const auto measurement=measure(
                *named_integrator.integrator,
                benchmark_case,
                step,
                repetitions);

            const double time_us_per_step=
                std::chrono::duration<double,std::micro>(
                    measurement.elapsed).count()
                / static_cast<double>(repetitions);

            if(!(time_us_per_step>0.0) || !std::isfinite(time_us_per_step)) {
                LOGGING_PRINTLN(
                    "Invalid timing for integrator "<<named_integrator.name
                    <<" and function "<<benchmark_case.name
                    <<": "<<time_us_per_step<<" us")
                throw std::runtime_error("Invalid benchmark timing.");
            }

            LOGGING_PRINTLN_AT(2,
                "Benchmark result: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", time="<<time_us_per_step<<" us/step"
                <<", remainder="<<measurement.flow.error())

            LOGGING_PRINTLN_AT(3,
                "Benchmark component remainders: integrator="
                <<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", components="<<measurement.flow.errors())

            results.push_back({
                named_integrator.key,
                benchmark_case.name,
                repetitions,
                time_us_per_step,
                measurement.flow.error(),
                measurement.flow.errors()
            });
        }
    }

    return results;
}

void print_results(
        const StepSizeType& step,
        Nat repetitions,
        const std::vector<BenchmarkResult>& results)
{
    constexpr std::size_t column_spacing=2u;

    struct DisplayRow {
        std::string integrator;
        std::string function;
        std::string time_us_per_step;
        std::string remainder;
        std::string component_remainders;
    };

    auto to_string=[](const auto& value) {
        std::ostringstream stream;
        stream<<std::fixed<<std::setprecision(6)<<value;
        return stream.str();
    };

    std::vector<DisplayRow> rows;
    rows.reserve(results.size());
    for(const auto& result : results) {
        rows.push_back({
            result.integrator,
            result.function,
            to_string(result.time_us_per_step),
            to_string(result.remainder),
            to_string(result.component_remainders)
        });
    }

    std::size_t integrator_width=std::string("integrator").size();
    std::size_t function_width=std::string("function").size();
    std::size_t time_width=std::string("time_us_per_step").size();
    std::size_t remainder_width=std::string("remainder").size();
    std::size_t component_remainders_width=
        std::string("component_remainders").size();

    for(const auto& row : rows) {
        integrator_width=std::max(integrator_width,row.integrator.size());
        function_width=std::max(function_width,row.function.size());
        time_width=std::max(time_width,row.time_us_per_step.size());
        remainder_width=std::max(remainder_width,row.remainder.size());
        component_remainders_width=std::max(
            component_remainders_width,
            row.component_remainders.size());
    }

    const std::size_t table_width=
        integrator_width+function_width+time_width+remainder_width
        +component_remainders_width+4u*column_spacing;
    const std::string separator(table_width,'-');

    const auto print_column=[](
            std::ostringstream& table,
            const std::string& value,
            std::size_t width)
    {
        table
            <<std::left
            <<std::setw(static_cast<int>(width))
            <<value;
    };

    const auto print_metrics=[
            &print_column,
            time_width,
            remainder_width,
            component_remainders_width](
            std::ostringstream& table,
            const DisplayRow& row)
    {
        table<<std::string(column_spacing,' ');
        print_column(table,row.time_us_per_step,time_width);
        table<<std::string(column_spacing,' ');
        print_column(table,row.remainder,remainder_width);
        table<<std::string(column_spacing,' ');
        print_column(
            table,
            row.component_remainders,
            component_remainders_width);
    };

    const auto print_metric_headers=[
            &print_column,
            time_width,
            remainder_width,
            component_remainders_width](
            std::ostringstream& table)
    {
        table<<std::string(column_spacing,' ');
        print_column(table,"time_us_per_step",time_width);
        table<<std::string(column_spacing,' ');
        print_column(table,"remainder",remainder_width);
        table<<std::string(column_spacing,' ');
        print_column(
            table,
            "component_remainders",
            component_remainders_width);
    };

    std::ostringstream table;
    table
        <<"Single-step integrator benchmark\n"
        <<"step="<<step
        <<", step_maximum_error="<<benchmark_step_maximum_error
        <<", repetitions="<<repetitions<<"\n\n";

    auto by_integrator=rows;
    std::sort(
        by_integrator.begin(),
        by_integrator.end(),
        [](const DisplayRow& lhs,const DisplayRow& rhs) {
            if(lhs.integrator!=rhs.integrator) {
                return lhs.integrator<rhs.integrator;
            }
            return lhs.function<rhs.function;
        });

    table<<"View by integrator\n";
    print_column(table,"integrator",integrator_width);
    table<<std::string(column_spacing,' ');
    print_column(table,"function",function_width);
    print_metric_headers(table);
    table<<"\n"<<separator<<"\n";

    for(std::size_t row_index=0u;
            row_index!=by_integrator.size();
            ++row_index) {
        const auto& row=by_integrator[row_index];
        print_column(table,row.integrator,integrator_width);
        table<<std::string(column_spacing,' ');
        print_column(table,row.function,function_width);
        print_metrics(table,row);
        table<<"\n";

        const bool end_of_group=
            row_index+1u==by_integrator.size()
            || by_integrator[row_index+1u].integrator!=row.integrator;
        if(end_of_group) {
            table<<separator<<"\n";
        }
    }

    auto by_function=rows;
    std::sort(
        by_function.begin(),
        by_function.end(),
        [](const DisplayRow& lhs,const DisplayRow& rhs) {
            if(lhs.function!=rhs.function) {
                return lhs.function<rhs.function;
            }
            return lhs.integrator<rhs.integrator;
        });

    table<<"\nView by function\n";
    print_column(table,"function",function_width);
    table<<std::string(column_spacing,' ');
    print_column(table,"integrator",integrator_width);
    print_metric_headers(table);
    table<<"\n"<<separator<<"\n";

    for(std::size_t row_index=0u;
            row_index!=by_function.size();
            ++row_index) {
        const auto& row=by_function[row_index];
        print_column(table,row.function,function_width);
        table<<std::string(column_spacing,' ');
        print_column(table,row.integrator,integrator_width);
        print_metrics(table,row);
        table<<"\n";

        const bool end_of_group=
            row_index+1u==by_function.size()
            || by_function[row_index+1u].function!=row.function;
        if(end_of_group) {
            table<<separator<<"\n";
        }
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
        <<", benchmark cases="<<cases.size())

    const auto integrators=select_integrators(
        selection.integrator,benchmark_step_maximum_error);
    const StepSizeType step=benchmark_step();
    const Nat repetitions=determine_repetitions(integrators,cases,step);
    const auto results=run_benchmarks(
        integrators,cases,step,repetitions);

    print_results(step,repetitions,results);
    return 0;
}
