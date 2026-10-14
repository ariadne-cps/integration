#pragma once

#include <cstdlib>
#include <iostream>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ariadne-kernel.hpp"

namespace IntegratorBenchmark {

using namespace Ariadne;

enum class AlgebraicClass {
    Affine,
    Polynomial,
    Analytic
};

enum class Sparsity {
    Sparse,
    Intermediate,
    Dense
};

enum class CouplingStrength {
    Weak,
    Moderate,
    Strong
};

enum class CouplingStructure {
    Decoupled,
    Directed,
    Bidirectional,
    Global
};

enum class FlowDegreeGrowth {
    None,
    Slow,
    Moderate,
    Fast,
    Unbounded
};

enum class ScaleSeparation {
    Low,
    Moderate,
    High
};

enum class Anisotropy {
    Low,
    Moderate,
    High
};

enum class RotationShear {
    Low,
    Moderate,
    High
};

struct DegreeRange {
    std::optional<Nat> minimum;
    std::optional<Nat> maximum;

    bool contains(Nat polynomial_degree_value) const {
        return (!minimum || polynomial_degree_value>=*minimum)
            && (!maximum || polynomial_degree_value<=*maximum);
    }
};

struct BenchmarkCharacteristics {
    Nat dimension;
    AlgebraicClass algebraic_class;
    std::optional<Nat> polynomial_degree;
    Sparsity sparsity;
    CouplingStrength coupling_strength;
    CouplingStructure coupling_structure;
    FlowDegreeGrowth flow_degree_growth;
    ScaleSeparation scale_separation;
    Anisotropy anisotropy;
    RotationShear rotation_shear;
};

struct BenchmarkCase {
    std::string name;
    ValidatedVectorMultivariateFunction vector_field;
    ExactBoxType domain;
    BenchmarkCharacteristics characteristics;
};

struct BenchmarkSelector {
    std::set<Nat> dimensions;
    std::set<AlgebraicClass> algebraic_classes;
    std::optional<DegreeRange> polynomial_degree;
    std::set<Sparsity> sparsities;
    std::set<CouplingStrength> coupling_strengths;
    std::set<CouplingStructure> coupling_structures;
    std::set<FlowDegreeGrowth> flow_degree_growths;
    std::set<ScaleSeparation> scale_separations;
    std::set<Anisotropy> anisotropies;
    std::set<RotationShear> rotation_shears;
};

struct BenchmarkCommandLineSelection {
    std::string integrator="all";
    BenchmarkSelector benchmarks;
    bool list_cases=false;
};

struct NamedIntegrator {
    std::string key;
    std::string name;
    std::unique_ptr<IntegratorInterface> integrator;
};

inline ExactBoxType box2() {
    return ExactBoxType({
        {-0.25_x,0.25_x},
        {-0.25_x,0.25_x}
    });
}

inline ExactBoxType box4() {
    return ExactBoxType({
        {-0.25_x,0.25_x},
        {-0.25_x,0.25_x},
        {-0.25_x,0.25_x},
        {-0.25_x,0.25_x}
    });
}

inline BenchmarkCase make_affine_2d_dense_rotating() {
    // A1: isolate geometric wrapping in a linear system.  The dense,
    // non-normal linear part rotates and shears sets without introducing
    // polynomial truncation error.  A stateful QR-preconditioned method
    // should benefit from an oriented affine/Taylor-model representation,
    // whereas an axis-aligned box is expected to wrap.
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[0]-2*x[1],
        2*x[0]-x[1]
    });
    return {"affine_2d_dense_rotating",f,box2(),{
        2u,
        AlgebraicClass::Affine,
        1u,
        Sparsity::Dense,
        CouplingStrength::Strong,
        CouplingStructure::Bidirectional,
        FlowDegreeGrowth::None,
        ScaleSeparation::Low,
        Anisotropy::Moderate,
        RotationShear::High
    }};
}

inline BenchmarkCase make_affine_4d_sparse_multiscale() {
    // A2: isolate time-scale separation in an axis-aligned, sparse affine
    // system with only weak one-way coupling.  This should favour methods
    // that preserve simple scaled coordinates; QR rotation is not expected
    // to help much while the principal directions remain axis-aligned.
    auto x=ValidatedScalarMultivariateFunction::coordinates(4);
    ValidatedVectorMultivariateFunction f({
        -16*x[0]+x[1],
        -4*x[1],
        -x[2]+x[3]/4,
        -x[3]/4
    });
    return {"affine_4d_sparse_multiscale",f,box4(),{
        4u,
        AlgebraicClass::Affine,
        1u,
        Sparsity::Sparse,
        CouplingStrength::Weak,
        CouplingStructure::Directed,
        FlowDegreeGrowth::None,
        ScaleSeparation::High,
        Anisotropy::High,
        RotationShear::Low
    }};
}

inline BenchmarkCase make_polynomial_low_2d_sparse_weak() {
    // P1: baseline low-degree nonlinear case.  Sparsity and weak coupling
    // keep mixed dependencies limited and the spatial degree of the flow
    // grows gradually.  A graded Taylor representation should be able to
    // exploit this structure without requiring geometric preconditioning.
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[0]+sqr(x[0]),
        -x[1]+x[0]/4
    });
    return {"polynomial_low_2d_sparse_weak",f,box2(),{
        2u,
        AlgebraicClass::Polynomial,
        2u,
        Sparsity::Sparse,
        CouplingStrength::Weak,
        CouplingStructure::Directed,
        FlowDegreeGrowth::Slow,
        ScaleSeparation::Low,
        Anisotropy::Low,
        RotationShear::Low
    }};
}

inline BenchmarkCase make_polynomial_low_2d_dense_strong() {
    // P2: keep degree low but make both components strongly coupled.  The
    // purpose is to stress the generation of mixed monomials independently
    // of high polynomial degree.  A full multivariate Taylor model, and in
    // particular a total-degree series method, should be a natural fit.
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[0]+x[1]+x[0]*x[1],
        x[0]-x[1]+sqr(x[0])-sqr(x[1])
    });
    return {"polynomial_low_2d_dense_strong",f,box2(),{
        2u,
        AlgebraicClass::Polynomial,
        2u,
        Sparsity::Dense,
        CouplingStrength::Strong,
        CouplingStructure::Bidirectional,
        FlowDegreeGrowth::Slow,
        ScaleSeparation::Low,
        Anisotropy::Moderate,
        RotationShear::Moderate
    }};
}

inline BenchmarkCase make_polynomial_low_2d_sparse_fastgrowth() {
    // P3: low vector-field degree but fast growth of spatial degree in the
    // time-series coefficients, with sparse asymmetric coupling.  This is
    // intended to expose insufficient MaximumSpacialOrder in graded methods;
    // an oriented set representation cannot compensate for missing spatial
    // polynomial order.
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[1],
        x[1]-x[0]-sqr(x[0])*x[1]
    });
    return {"polynomial_low_2d_sparse_fastgrowth",f,box2(),{
        2u,
        AlgebraicClass::Polynomial,
        3u,
        Sparsity::Sparse,
        CouplingStrength::Strong,
        CouplingStructure::Bidirectional,
        FlowDegreeGrowth::Fast,
        ScaleSeparation::Low,
        Anisotropy::Moderate,
        RotationShear::Moderate
    }};
}

inline BenchmarkCase make_polynomial_low_2d_rotating_anisotropic() {
    // G1: moderate polynomial complexity with a strongly non-normal linear
    // part.  The main stress is anisotropy plus rotation/shear rather than
    // rapid degree growth.  This is the case where stateful QR
    // preconditioning should have the clearest opportunity to outperform an
    // axis-aligned Taylor-model representation.
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[0]-8*x[1]+sqr(x[0])/4,
        x[0]-x[1]
    });
    return {"polynomial_low_2d_rotating_anisotropic",f,box2(),{
        2u,
        AlgebraicClass::Polynomial,
        2u,
        Sparsity::Dense,
        CouplingStrength::Strong,
        CouplingStructure::Bidirectional,
        FlowDegreeGrowth::Slow,
        ScaleSeparation::Low,
        Anisotropy::High,
        RotationShear::High
    }};
}

inline BenchmarkCase make_polynomial_low_4d_sparse_strong_propagation() {
    // P6: the vector field is sparse, but coupling forms a nonlinear ring so
    // successive time derivatives propagate dependencies through all four
    // variables.  This stresses flow densification rather than initial
    // density.  Sparse/graded representations may win early, while a full
    // multivariate Taylor model becomes relevant as mixed terms accumulate.
    auto x=ValidatedScalarMultivariateFunction::coordinates(4);
    ValidatedVectorMultivariateFunction f({
        x[1]+x[0]*x[1],
        x[2]+x[1]*x[2],
        x[3]+x[2]*x[3],
        -x[0]+x[3]*x[0]
    });
    return {"polynomial_low_4d_sparse_strong_propagation",f,box4(),{
        4u,
        AlgebraicClass::Polynomial,
        2u,
        Sparsity::Sparse,
        CouplingStrength::Strong,
        CouplingStructure::Global,
        FlowDegreeGrowth::Slow,
        ScaleSeparation::Low,
        Anisotropy::Moderate,
        RotationShear::Moderate
    }};
}

inline BenchmarkCase make_analytic_2d_dense() {
    // T1: dense non-polynomial analytic dynamics.  This stresses elementary
    // function expansion and rigorous remainder bounds rather than finite
    // polynomial degree.  Picard/series-bounder approaches are expected to
    // be competitive when paired with Taylor models carrying interval
    // remainders.
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        sin(x[0])+exp(x[1])-1,
        cos(x[1])+x[0]*x[1]-1
    });
    return {"analytic_2d_dense",f,box2(),{
        2u,
        AlgebraicClass::Analytic,
        std::nullopt,
        Sparsity::Dense,
        CouplingStrength::Strong,
        CouplingStructure::Bidirectional,
        FlowDegreeGrowth::Unbounded,
        ScaleSeparation::Low,
        Anisotropy::Moderate,
        RotationShear::Moderate
    }};
}

inline std::vector<BenchmarkCase> make_benchmark_library() {
    std::vector<BenchmarkCase> cases;
    cases.push_back(make_affine_2d_dense_rotating());
    cases.push_back(make_affine_4d_sparse_multiscale());
    cases.push_back(make_polynomial_low_2d_sparse_weak());
    cases.push_back(make_polynomial_low_2d_dense_strong());
    cases.push_back(make_polynomial_low_2d_sparse_fastgrowth());
    cases.push_back(make_polynomial_low_2d_rotating_anisotropic());
    cases.push_back(make_polynomial_low_4d_sparse_strong_propagation());
    cases.push_back(make_analytic_2d_dense());
    return cases;
}

template<class T>
inline bool matches_allowed_values(
        const std::set<T>& allowed,
        const T& actual)
{
    return allowed.empty() || allowed.find(actual)!=allowed.end();
}

inline bool matches(
        const BenchmarkCase& benchmark_case,
        const BenchmarkSelector& selector)
{
    const auto& c=benchmark_case.characteristics;

    if(!matches_allowed_values(selector.dimensions,c.dimension)) {
        return false;
    }
    if(!matches_allowed_values(
            selector.algebraic_classes,c.algebraic_class)) {
        return false;
    }
    if(selector.polynomial_degree) {
        if(!c.polynomial_degree
                || !selector.polynomial_degree->contains(
                    *c.polynomial_degree)) {
            return false;
        }
    }
    if(!matches_allowed_values(selector.sparsities,c.sparsity)) {
        return false;
    }
    if(!matches_allowed_values(
            selector.coupling_strengths,c.coupling_strength)) {
        return false;
    }
    if(!matches_allowed_values(
            selector.coupling_structures,c.coupling_structure)) {
        return false;
    }
    if(!matches_allowed_values(
            selector.flow_degree_growths,c.flow_degree_growth)) {
        return false;
    }
    if(!matches_allowed_values(
            selector.scale_separations,c.scale_separation)) {
        return false;
    }
    if(!matches_allowed_values(selector.anisotropies,c.anisotropy)) {
        return false;
    }
    if(!matches_allowed_values(
            selector.rotation_shears,c.rotation_shear)) {
        return false;
    }
    return true;
}

inline std::vector<BenchmarkCase> select_benchmark_cases(
        const std::vector<BenchmarkCase>& library,
        const BenchmarkSelector& selector)
{
    std::vector<BenchmarkCase> selected;
    for(const auto& benchmark_case : library) {
        if(matches(benchmark_case,selector)) {
            selected.push_back(benchmark_case);
        }
    }
    return selected;
}

inline std::vector<BenchmarkCase> make_cases(
        const BenchmarkSelector& selector)
{
    return select_benchmark_cases(make_benchmark_library(),selector);
}

inline std::vector<BenchmarkCase> make_cases() {
    return make_benchmark_library();
}

inline const char* characteristic_name(AlgebraicClass value) {
    switch(value) {
        case AlgebraicClass::Affine: return "affine";
        case AlgebraicClass::Polynomial: return "polynomial";
        case AlgebraicClass::Analytic: return "analytic";
        default: return "unknown";
    }
}

inline const char* characteristic_name(Sparsity value) {
    switch(value) {
        case Sparsity::Sparse: return "sparse";
        case Sparsity::Intermediate: return "intermediate";
        case Sparsity::Dense: return "dense";
        default: return "unknown";
    }
}

inline const char* characteristic_name(CouplingStrength value) {
    switch(value) {
        case CouplingStrength::Weak: return "weak";
        case CouplingStrength::Moderate: return "moderate";
        case CouplingStrength::Strong: return "strong";
        default: return "unknown";
    }
}

inline const char* characteristic_name(CouplingStructure value) {
    switch(value) {
        case CouplingStructure::Decoupled: return "decoupled";
        case CouplingStructure::Directed: return "directed";
        case CouplingStructure::Bidirectional: return "bidirectional";
        case CouplingStructure::Global: return "global";
        default: return "unknown";
    }
}

inline const char* characteristic_name(FlowDegreeGrowth value) {
    switch(value) {
        case FlowDegreeGrowth::None: return "none";
        case FlowDegreeGrowth::Slow: return "slow";
        case FlowDegreeGrowth::Moderate: return "moderate";
        case FlowDegreeGrowth::Fast: return "fast";
        case FlowDegreeGrowth::Unbounded: return "unbounded";
        default: return "unknown";
    }
}

inline const char* characteristic_name(ScaleSeparation value) {
    switch(value) {
        case ScaleSeparation::Low: return "low";
        case ScaleSeparation::Moderate: return "moderate";
        case ScaleSeparation::High: return "high";
        default: return "unknown";
    }
}

inline const char* characteristic_name(Anisotropy value) {
    switch(value) {
        case Anisotropy::Low: return "low";
        case Anisotropy::Moderate: return "moderate";
        case Anisotropy::High: return "high";
        default: return "unknown";
    }
}

inline const char* characteristic_name(RotationShear value) {
    switch(value) {
        case RotationShear::Low: return "low";
        case RotationShear::Moderate: return "moderate";
        case RotationShear::High: return "high";
        default: return "unknown";
    }
}

inline std::string describe(const BenchmarkCase& benchmark_case) {
    const auto& c=benchmark_case.characteristics;
    std::ostringstream stream;
    stream
        <<"dimension="<<c.dimension
        <<", class="<<characteristic_name(c.algebraic_class)
        <<", degree=";
    if(c.polynomial_degree) {
        stream<<*c.polynomial_degree;
    } else {
        stream<<"n/a";
    }
    stream
        <<", sparsity="<<characteristic_name(c.sparsity)
        <<", coupling-strength="<<characteristic_name(c.coupling_strength)
        <<", coupling-structure="<<characteristic_name(c.coupling_structure)
        <<", flow-growth="<<characteristic_name(c.flow_degree_growth)
        <<", scale-separation="<<characteristic_name(c.scale_separation)
        <<", anisotropy="<<characteristic_name(c.anisotropy)
        <<", rotation-shear="<<characteristic_name(c.rotation_shear);
    return stream.str();
}

inline void print_benchmark_cases(
        const std::vector<BenchmarkCase>& cases)
{
    std::cout<<"Benchmark cases\n";
    for(const auto& benchmark_case : cases) {
        std::cout
            <<benchmark_case.name<<"\n"
            <<"  "<<describe(benchmark_case)<<"\n";
    }
}

inline std::vector<NamedIntegrator> make_all_integrators(double step_error) {
    constexpr DegreeType order_value=8;

    std::vector<NamedIntegrator> integrators;
    integrators.push_back({
        "taylor-picard",
        "TaylorPicardIntegrator",
        std::make_unique<TaylorPicardIntegrator>(StepMaximumError(step_error))
    });
    integrators.push_back({
        "graded-taylor-picard",
        "GradedTaylorPicardIntegrator",
        std::make_unique<GradedTaylorPicardIntegrator>(
            StepMaximumError(step_error),Order(order_value))
    });
    integrators.push_back({
        "taylor-series",
        "TaylorSeriesIntegrator",
        std::make_unique<TaylorSeriesIntegrator>(
            StepMaximumError(step_error),Order(order_value))
    });
    integrators.push_back({
        "taylor-series-bounder",
        "TaylorSeriesBounderIntegrator",
        std::make_unique<TaylorSeriesBounderIntegrator>(
            StepMaximumError(step_error),Order(order_value))
    });
    integrators.push_back({
        "graded-taylor-series",
        "GradedTaylorSeriesIntegrator",
        std::make_unique<GradedTaylorSeriesIntegrator>(
            StepMaximumError(step_error))
    });
    integrators.push_back({
        "preconditioned-graded-taylor-series",
        "PreconditionedGradedTaylorSeriesIntegrator",
        std::make_unique<PreconditionedGradedTaylorSeriesIntegrator>(
            StepMaximumError(step_error))
    });
    return integrators;
}

inline std::vector<NamedIntegrator> select_integrators(
        const std::string& selection,
        double step_error)
{
    auto all=make_all_integrators(step_error);
    if(selection=="all") {
        return all;
    }

    std::vector<NamedIntegrator> selected;
    for(auto& named_integrator : all) {
        if(named_integrator.key==selection) {
            selected.push_back(std::move(named_integrator));
            return selected;
        }
    }

    std::cerr << "Unknown integrator '" << selection << "'.\n";
    std::cerr << "Available integrators:\n";
    std::cerr << "  all\n";
    for(const auto& named_integrator : all) {
        std::cerr << "  " << named_integrator.key << "\n";
    }
    std::exit(1);
}

inline bool is_integrator_selection(const std::string& value) {
    return value=="all"
        || value=="taylor-picard"
        || value=="graded-taylor-picard"
        || value=="taylor-series"
        || value=="taylor-series-bounder"
        || value=="graded-taylor-series"
        || value=="preconditioned-graded-taylor-series";
}

inline std::vector<std::string> split_selector_values(
        const std::string& values)
{
    if(values.empty()) {
        throw std::runtime_error("Empty benchmark selector value.");
    }

    std::vector<std::string> result;
    std::size_t begin=0u;
    while(begin<=values.size()) {
        const std::size_t end=values.find(',',begin);
        const std::string value=values.substr(
            begin,
            end==std::string::npos
                ? std::string::npos
                : end-begin);
        if(value.empty()) {
            throw std::runtime_error(
                "Empty value in benchmark selector list.");
        }
        result.push_back(value);
        if(end==std::string::npos) {
            break;
        }
        begin=end+1u;
    }
    return result;
}

inline Nat parse_selector_nat(const std::string& value) {
    if(value.empty() || value[0]=='-') {
        throw std::runtime_error(
            "Expected a non-negative integer, got '"+value+"'.");
    }

    std::size_t consumed=0u;
    const unsigned long parsed=std::stoul(value,&consumed);
    if(consumed!=value.size()) {
        throw std::runtime_error(
            "Expected a non-negative integer, got '"+value+"'.");
    }
    return static_cast<Nat>(parsed);
}

inline DegreeRange parse_degree_range(const std::string& value) {
    const std::size_t separator=value.find(':');
    if(separator==std::string::npos) {
        const Nat polynomial_degree_value=parse_selector_nat(value);
        return {polynomial_degree_value,polynomial_degree_value};
    }
    if(value.find(':',separator+1u)!=std::string::npos) {
        throw std::runtime_error(
            "Invalid degree range '"+value+"'.");
    }

    const std::string minimum=value.substr(0u,separator);
    const std::string maximum=value.substr(separator+1u);
    if(minimum.empty() && maximum.empty()) {
        throw std::runtime_error(
            "Degree range ':' does not restrict the degree.");
    }

    DegreeRange range;
    if(!minimum.empty()) {
        range.minimum=parse_selector_nat(minimum);
    }
    if(!maximum.empty()) {
        range.maximum=parse_selector_nat(maximum);
    }
    if(range.minimum && range.maximum
            && *range.minimum>*range.maximum) {
        throw std::runtime_error(
            "Degree range minimum exceeds maximum in '"+value+"'.");
    }
    return range;
}

inline AlgebraicClass parse_algebraic_class(const std::string& value) {
    if(value=="affine") return AlgebraicClass::Affine;
    if(value=="polynomial") return AlgebraicClass::Polynomial;
    if(value=="analytic") return AlgebraicClass::Analytic;
    throw std::runtime_error(
        "Unknown algebraic class '"+value+"'.");
}

inline Sparsity parse_sparsity(const std::string& value) {
    if(value=="sparse") return Sparsity::Sparse;
    if(value=="intermediate") return Sparsity::Intermediate;
    if(value=="dense") return Sparsity::Dense;
    throw std::runtime_error("Unknown sparsity '"+value+"'.");
}

inline CouplingStrength parse_coupling_strength(
        const std::string& value)
{
    if(value=="weak") return CouplingStrength::Weak;
    if(value=="moderate") return CouplingStrength::Moderate;
    if(value=="strong") return CouplingStrength::Strong;
    throw std::runtime_error(
        "Unknown coupling strength '"+value+"'.");
}

inline CouplingStructure parse_coupling_structure(
        const std::string& value)
{
    if(value=="decoupled") return CouplingStructure::Decoupled;
    if(value=="directed") return CouplingStructure::Directed;
    if(value=="bidirectional") return CouplingStructure::Bidirectional;
    if(value=="global") return CouplingStructure::Global;
    throw std::runtime_error(
        "Unknown coupling structure '"+value+"'.");
}

inline FlowDegreeGrowth parse_flow_degree_growth(
        const std::string& value)
{
    if(value=="none") return FlowDegreeGrowth::None;
    if(value=="slow") return FlowDegreeGrowth::Slow;
    if(value=="moderate") return FlowDegreeGrowth::Moderate;
    if(value=="fast") return FlowDegreeGrowth::Fast;
    if(value=="unbounded") return FlowDegreeGrowth::Unbounded;
    throw std::runtime_error(
        "Unknown flow-degree growth '"+value+"'.");
}

inline ScaleSeparation parse_scale_separation(
        const std::string& value)
{
    if(value=="low") return ScaleSeparation::Low;
    if(value=="moderate") return ScaleSeparation::Moderate;
    if(value=="high") return ScaleSeparation::High;
    throw std::runtime_error(
        "Unknown scale separation '"+value+"'.");
}

inline Anisotropy parse_anisotropy(const std::string& value) {
    if(value=="low") return Anisotropy::Low;
    if(value=="moderate") return Anisotropy::Moderate;
    if(value=="high") return Anisotropy::High;
    throw std::runtime_error(
        "Unknown anisotropy '"+value+"'.");
}

inline RotationShear parse_rotation_shear(const std::string& value) {
    if(value=="low") return RotationShear::Low;
    if(value=="moderate") return RotationShear::Moderate;
    if(value=="high") return RotationShear::High;
    throw std::runtime_error(
        "Unknown rotation/shear '"+value+"'.");
}

inline List<ArgumentParser> make_benchmark_argument_parsers(
        BenchmarkCommandLineSelection& selection)
{
    List<ArgumentParser> parsers;

    parsers.append(make_valued_argument_parser(
        "",
        "integrator",
        "Choose the integrator or 'all'",
        [&selection](String const& value) -> VoidFunction {
            const std::string parsed=value;
            if(!is_integrator_selection(parsed)) {
                throw std::runtime_error(
                    "Unknown integrator '"+parsed+"'.");
            }
            return [&selection,parsed] {
                selection.integrator=parsed;
            };
        }));

    parsers.append(make_unvalued_argument_parser(
        "",
        "list-cases",
        "List benchmark cases selected by the filters",
        [&selection] {
            selection.list_cases=true;
        }));

    parsers.append(make_valued_argument_parser(
        "",
        "dimension",
        "Select dimensions as a comma-separated list, e.g. 2,4",
        [&selection](String const& value) -> VoidFunction {
            std::set<Nat> parsed;
            for(const auto& item : split_selector_values(value)) {
                parsed.insert(parse_selector_nat(item));
            }
            return [&selection,parsed] {
                selection.benchmarks.dimensions.insert(
                    parsed.begin(),parsed.end());
            };
        }));

    parsers.append(make_valued_argument_parser(
        "",
        "class",
        "Select algebraic classes: affine, polynomial, analytic",
        [&selection](String const& value) -> VoidFunction {
            std::set<AlgebraicClass> parsed;
            for(const auto& item : split_selector_values(value)) {
                parsed.insert(parse_algebraic_class(item));
            }
            return [&selection,parsed] {
                selection.benchmarks.algebraic_classes.insert(
                    parsed.begin(),parsed.end());
            };
        }));

    parsers.append(make_valued_argument_parser(
        "",
        "degree",
        "Select polynomial degree range: N, N:M, N:, or :M",
        [&selection](String const& value) -> VoidFunction {
            const DegreeRange parsed=parse_degree_range(value);
            return [&selection,parsed] {
                selection.benchmarks.polynomial_degree=parsed;
            };
        }));

    parsers.append(make_valued_argument_parser(
        "",
        "sparsity",
        "Select sparsity: sparse, intermediate, dense",
        [&selection](String const& value) -> VoidFunction {
            std::set<Sparsity> parsed;
            for(const auto& item : split_selector_values(value)) {
                parsed.insert(parse_sparsity(item));
            }
            return [&selection,parsed] {
                selection.benchmarks.sparsities.insert(
                    parsed.begin(),parsed.end());
            };
        }));

    parsers.append(make_valued_argument_parser(
        "",
        "coupling-strength",
        "Select coupling strength: weak, moderate, strong",
        [&selection](String const& value) -> VoidFunction {
            std::set<CouplingStrength> parsed;
            for(const auto& item : split_selector_values(value)) {
                parsed.insert(parse_coupling_strength(item));
            }
            return [&selection,parsed] {
                selection.benchmarks.coupling_strengths.insert(
                    parsed.begin(),parsed.end());
            };
        }));

    parsers.append(make_valued_argument_parser(
        "",
        "coupling-structure",
        "Select coupling structure: decoupled, directed, bidirectional, global",
        [&selection](String const& value) -> VoidFunction {
            std::set<CouplingStructure> parsed;
            for(const auto& item : split_selector_values(value)) {
                parsed.insert(parse_coupling_structure(item));
            }
            return [&selection,parsed] {
                selection.benchmarks.coupling_structures.insert(
                    parsed.begin(),parsed.end());
            };
        }));

    parsers.append(make_valued_argument_parser(
        "",
        "flow-growth",
        "Select flow-degree growth: none, slow, moderate, fast, unbounded",
        [&selection](String const& value) -> VoidFunction {
            std::set<FlowDegreeGrowth> parsed;
            for(const auto& item : split_selector_values(value)) {
                parsed.insert(parse_flow_degree_growth(item));
            }
            return [&selection,parsed] {
                selection.benchmarks.flow_degree_growths.insert(
                    parsed.begin(),parsed.end());
            };
        }));

    parsers.append(make_valued_argument_parser(
        "",
        "scale-separation",
        "Select scale separation: low, moderate, high",
        [&selection](String const& value) -> VoidFunction {
            std::set<ScaleSeparation> parsed;
            for(const auto& item : split_selector_values(value)) {
                parsed.insert(parse_scale_separation(item));
            }
            return [&selection,parsed] {
                selection.benchmarks.scale_separations.insert(
                    parsed.begin(),parsed.end());
            };
        }));

    parsers.append(make_valued_argument_parser(
        "",
        "anisotropy",
        "Select anisotropy: low, moderate, high",
        [&selection](String const& value) -> VoidFunction {
            std::set<Anisotropy> parsed;
            for(const auto& item : split_selector_values(value)) {
                parsed.insert(parse_anisotropy(item));
            }
            return [&selection,parsed] {
                selection.benchmarks.anisotropies.insert(
                    parsed.begin(),parsed.end());
            };
        }));

    parsers.append(make_valued_argument_parser(
        "",
        "rotation-shear",
        "Select rotation/shear: low, moderate, high",
        [&selection](String const& value) -> VoidFunction {
            std::set<RotationShear> parsed;
            for(const auto& item : split_selector_values(value)) {
                parsed.insert(parse_rotation_shear(item));
            }
            return [&selection,parsed] {
                selection.benchmarks.rotation_shears.insert(
                    parsed.begin(),parsed.end());
            };
        }));

    return parsers;
}

} // namespace IntegratorBenchmark
