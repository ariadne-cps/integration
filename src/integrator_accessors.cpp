/***************************************************************************
 *            src/integrator_accessors.cpp
 *
 *  Copyright  2026  Luca Geretti
 ****************************************************************************/

/** @file src/integrator_accessors.cpp
 *  @brief Non-template integration configuration and state accessors.
 */

#include "integration/integrator.hpp"

namespace Ariadne {

const Generator<StepMaximumError> step_maximum_error{};
const Generator<StepSweepThreshold> step_sweep_threshold{};
const Generator<Order> order{};
const Generator<SpacialOrder> spacial_order{};
const Generator<TemporalOrder> temporal_order{};
const Generator<MinimumTemporalOrder> minimum_temporal_order{};
const Generator<MinimumSpacialOrder> minimum_spacial_order{};
const Generator<MaximumSpacialOrder> maximum_spacial_order{};
const Generator<MaximumTemporalOrder> maximum_temporal_order{};


StepMaximumError::StepMaximumError(ApproximateDouble x) : Attribute<ApproximateDouble>(x) {}

StepMaximumError::StepMaximumError(double x) : Attribute<ApproximateDouble>(x) {}

StepSweepThreshold::StepSweepThreshold(const ApproximateDouble& value) : Attribute<ApproximateDouble>(value) {}

Order::Order(DegreeType v) : Attribute<DegreeType>(v) {}

SpacialOrder::SpacialOrder(DegreeType v) : Attribute<DegreeType>(v) {}

TemporalOrder::TemporalOrder(DegreeType v) : Attribute<DegreeType>(v) {}

MinimumSpacialOrder::MinimumSpacialOrder(DegreeType v) : Attribute<DegreeType>(v) {}

MinimumTemporalOrder::MinimumTemporalOrder(DegreeType v) : Attribute<DegreeType>(v) {}

MaximumSpacialOrder::MaximumSpacialOrder(DegreeType v) : Attribute<DegreeType>(v) {}

MaximumTemporalOrder::MaximumTemporalOrder(DegreeType v) : Attribute<DegreeType>(v) {}

FlowStepModelType::FlowStepModelType(const ValidatedVectorMultivariateFunctionPatch& fsmt) : ValidatedVectorMultivariateFunctionPatch(fsmt) {}

DegreeType TaylorPicardIntegrator::minimum_temporal_order() const { return this->_minimum_temporal_order; }

Void TaylorPicardIntegrator::set_minimum_temporal_order(Nat m) {
        ARIADNE_PRECONDITION(m <= std::numeric_limits<DegreeType>::max());
        this->_minimum_temporal_order = static_cast<DegreeType>(m);
    }

DegreeType TaylorPicardIntegrator::maximum_temporal_order() const { return this->_maximum_temporal_order; }

Void TaylorPicardIntegrator::set_maximum_temporal_order(Nat m) {
        ARIADNE_PRECONDITION(m <= std::numeric_limits<DegreeType>::max());
        this->_maximum_temporal_order = static_cast<DegreeType>(m);
    }

const Sweeper<FloatDP>& TaylorPicardIntegrator::sweeper() const { return this->_sweeper; }

ExactDouble TaylorPicardIntegrator::step_maximum_error() const { return this->_step_maximum_error; }

Void TaylorPicardIntegrator::set_step_maximum_error(ApproximateDouble e) { _step_maximum_error = cast_exact(e); }

TaylorPicardIntegrator* TaylorPicardIntegrator::clone() const { return new TaylorPicardIntegrator(*this); }

DegreeType GradedTaylorPicardIntegrator::order() const { return this->_order; }

Void GradedTaylorPicardIntegrator::set_order(Nat m) {
        /** @brief Construct the integration object. */
        ARIADNE_PRECONDITION(m <= std::numeric_limits<DegreeType>::max());
        const DegreeType order = static_cast<DegreeType>(m);
        this->_order = order;
        this->_sweeper = GradedThresholdSweeper<FloatDP>(DoublePrecision(), order, FloatDP(cast_exact(this->_sweep_threshold), DoublePrecision()));
        this->set_function_factory(ValidatedFunctionPatchFactory(make_taylor_function_patch_factory(this->_sweeper)));
    }

ExactDouble GradedTaylorPicardIntegrator::step_maximum_error() const { return this->_step_maximum_error; }

Void GradedTaylorPicardIntegrator::set_step_maximum_error(ApproximateDouble e) { _step_maximum_error = cast_exact(e); }

ExactDouble GradedTaylorPicardIntegrator::error_refinement_minimum_improvement_percentage() const { return this->_error_refinement_minimum_improvement_percentage; }

Void GradedTaylorPicardIntegrator::set_error_refinement_minimum_improvement_percentage(ApproximateDouble e) { _error_refinement_minimum_improvement_percentage = cast_exact(e); }

Void GradedTaylorPicardIntegrator::set_maximum_error_refinement_iterations(Nat n) { _maximum_error_refinement_iterations = n; }

Nat GradedTaylorPicardIntegrator::maximum_error_refinement_iterations() const { return _maximum_error_refinement_iterations; }

Void GradedTaylorPicardIntegrator::set_diagnostics(Bool diagnostics) { _diagnostics = diagnostics; }

GradedTaylorPicardIntegrator* GradedTaylorPicardIntegrator::clone() const { return new GradedTaylorPicardIntegrator(*this); }

DegreeType TaylorSeriesIntegrator::order() const { return this->_order; }

Void TaylorSeriesIntegrator::set_order(DegreeType n) { this->_order = n; }

Bool TaylorSeriesIntegrator::has_step_maximum_error() const { return this->_has_step_maximum_error; }

ExactDouble TaylorSeriesIntegrator::step_maximum_error() const {
        /** @brief Construct the integration object. */
        ARIADNE_PRECONDITION(this->has_step_maximum_error());
        return this->_step_maximum_error;
    }

Void TaylorSeriesIntegrator::set_step_maximum_error(ApproximateDouble e) {
        this->_step_maximum_error = cast_exact(e);
        this->_has_step_maximum_error = true;
    }

const Sweeper<FloatDP>& TaylorSeriesIntegrator::sweeper() const { return this->_sweeper; }

Void TaylorSeriesIntegrator::set_sweeper(const Sweeper<FloatDP>& sweeper) {
        _sweeper = sweeper;
        this->set_function_factory(ValidatedFunctionPatchFactory(make_taylor_function_patch_factory(_sweeper)));
    }

TaylorSeriesIntegrator* TaylorSeriesIntegrator::clone() const { return new TaylorSeriesIntegrator(*this); }

TaylorSeriesBounderIntegrator* TaylorSeriesBounderIntegrator::clone() const { return new TaylorSeriesBounderIntegrator(*this); }

DegreeType GradedTaylorSeriesIntegrator::minimum_spacial_order() const { return this->_minimum_spacial_order; }

Void GradedTaylorSeriesIntegrator::set_minimum_spacial_order(DegreeType n) { this->_minimum_spacial_order = n; }

DegreeType GradedTaylorSeriesIntegrator::minimum_temporal_order() const { return this->_minimum_temporal_order; }

Void GradedTaylorSeriesIntegrator::set_minimum_temporal_order(DegreeType m) { this->_minimum_temporal_order = m; }

DegreeType GradedTaylorSeriesIntegrator::maximum_spacial_order() const { return this->_maximum_spacial_order; }

Void GradedTaylorSeriesIntegrator::set_maximum_spacial_order(DegreeType n) { this->_maximum_spacial_order = n; }

DegreeType GradedTaylorSeriesIntegrator::maximum_temporal_order() const { return this->_maximum_temporal_order; }

Void GradedTaylorSeriesIntegrator::set_maximum_temporal_order(DegreeType m) { this->_maximum_temporal_order = m; }

const Sweeper<FloatDP>& GradedTaylorSeriesIntegrator::sweeper() const { return this->_sweeper; }

Void GradedTaylorSeriesIntegrator::set_sweeper(const Sweeper<FloatDP>& sweeper) {
        _sweeper = sweeper;
        this->set_function_factory(ValidatedFunctionPatchFactory(make_taylor_function_patch_factory(_sweeper)));
    }

ExactDouble GradedTaylorSeriesIntegrator::step_maximum_error() const { return this->_step_maximum_error; }

Void GradedTaylorSeriesIntegrator::set_step_maximum_error(ApproximateDouble e) { _step_maximum_error = cast_exact(e); }

GradedTaylorSeriesIntegrator* GradedTaylorSeriesIntegrator::clone() const { return new GradedTaylorSeriesIntegrator(*this); }

Bool PreconditionedTaylorSeriesState::has_physical_mapping() const { return static_cast<Bool>(_physical_mapping != nullptr); }

const ValidatedVectorMultivariateFunctionPatch& PreconditionedTaylorSeriesState::physical_mapping() const { return *_physical_mapping; }

const Vector<FloatDP>& PreconditionedTaylorSeriesState::centre() const { return _centre; }

const Matrix<FloatDP>& PreconditionedTaylorSeriesState::linear_map() const { return _linear_map; }

const ExactBoxType& PreconditionedTaylorSeriesState::local_domain() const { return _local_domain; }

const ValidatedVectorMultivariateFunctionPatch& PreconditionedTaylorSeriesState::normalised_mapping() const { return _normalised_mapping; }

const ExactBoxType PreconditionedTaylorSeriesState::parameter_domain() const { return _normalised_mapping.domain(); }

const StepSizeType& PreconditionedTaylorSeriesStep::time_step() const { return _time_step; }

const ValidatedVectorMultivariateFunctionPatch& PreconditionedTaylorSeriesStep::flowpipe_mapping() const { return _flowpipe_mapping; }

const ValidatedVectorMultivariateFunctionPatch& PreconditionedTaylorSeriesStep::evolved_mapping() const { return _evolved_mapping; }

const PreconditionedTaylorSeriesState& PreconditionedTaylorSeriesStep::final_state() const { return _final_state; }

TaylorSeriesPreconditioning PreconditionedGradedTaylorSeriesIntegrator::preconditioning() const { return _preconditioning; }

Void PreconditionedGradedTaylorSeriesIntegrator::set_preconditioning(TaylorSeriesPreconditioning value) { _preconditioning = value; }

Bool PreconditionedGradedTaylorSeriesIntegrator::diagnostics() const { return _diagnostics; }

Void PreconditionedGradedTaylorSeriesIntegrator::set_diagnostics(Bool value) { _diagnostics = value; }

Bool PreconditionedGradedTaylorSeriesIntegrator::carried_expansion_diagnostics() const { return _carried_expansion_diagnostics; }

Void PreconditionedGradedTaylorSeriesIntegrator::set_carried_expansion_diagnostics(Bool value) { _carried_expansion_diagnostics = value; }

PreconditionedGradedTaylorSeriesIntegrator* PreconditionedGradedTaylorSeriesIntegrator::clone() const {
        return new PreconditionedGradedTaylorSeriesIntegrator(*this);
    }

DegreeType AffineIntegrator::spacial_order() const { return this->_spacial_order; }

DegreeType AffineIntegrator::temporal_order() const { return this->_temporal_order; }

AffineIntegrator* AffineIntegrator::clone() const { return new AffineIntegrator(*this); }

PreconditionedTaylorSeriesState::PreconditionedTaylorSeriesState(
        Vector<FloatDP> centre,
        Matrix<FloatDP> linear_map,
        ExactBoxType local_domain,
        ValidatedVectorMultivariateFunctionPatch normalised_mapping) : _centre(std::move(centre)),
          _linear_map(std::move(linear_map)),
          _local_domain(std::move(local_domain)),
          _normalised_mapping(std::move(normalised_mapping)) {}

PreconditionedTaylorSeriesState::PreconditionedTaylorSeriesState(
        Vector<FloatDP> centre,
        Matrix<FloatDP> linear_map,
        ExactBoxType local_domain,
        ValidatedVectorMultivariateFunctionPatch normalised_mapping,
        ValidatedVectorMultivariateFunctionPatch physical_mapping) : _centre(std::move(centre)),
          _linear_map(std::move(linear_map)),
          _local_domain(std::move(local_domain)),
          _normalised_mapping(std::move(normalised_mapping)),
          _physical_mapping(std::make_shared<const ValidatedVectorMultivariateFunctionPatch>(
              std::move(physical_mapping))) {}

PreconditionedTaylorSeriesStep::PreconditionedTaylorSeriesStep(
        StepSizeType time_step,
        ValidatedVectorMultivariateFunctionPatch flowpipe_mapping,
        ValidatedVectorMultivariateFunctionPatch evolved_mapping,
        PreconditionedTaylorSeriesState final_state) : _time_step(time_step),
          _flowpipe_mapping(std::move(flowpipe_mapping)),
          _evolved_mapping(std::move(evolved_mapping)),
          _final_state(std::move(final_state)) {}

} // namespace Ariadne
