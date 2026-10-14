/***************************************************************************
 *            integration_submodule.cpp
 *
 *  Copyright  2009-26  Pieter Collins, Luca Geretti
 *
 ****************************************************************************/

#include "pybind11.hpp"
#include "kernel-utilities.hpp"

#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "integration/integrator_interface.hpp"
#include "integration/integrator.hpp"
#include "integration/runge_kutta_integrator.hpp"

using namespace Ariadne;

namespace Ariadne {

class IntegratorWrapper
  : public pybind11::wrapper<IntegratorInterface>
{
  public:
    IntegratorInterface* clone() const {
        return this->get_override("clone")(); }
    FlowStepModelType flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D) const {
        return this->get_override("flow_step")(vf,D); }
    FlowStepModelType flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, const StepSizeType& h) const {
        return this->get_override("flow_step")(vf,D,h); }
    FlowStepModelType flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, const StepSizeType& h, const UpperBoxType& B) const {
        return this->get_override("flow_step")(vf,D,h,B); }
    Void _write(OutputStream& os) const {
        this->get_override("_write")(os); }
};

} // namespace Ariadne

Void export_integrators(pybind11::module& module)
{
    pybind11::class_<FlowStepModelType,pybind11::bases<ValidatedVectorMultivariateFunctionPatch>> flow_step_model_class(module,"FlowStepModelType");
    flow_step_model_class.def("__str__", &__cstr__<FlowStepModelType>);

    pybind11::class_<FlowModelType> flow_model_class(module,"FlowModelType");
    flow_model_class.def("size",&FlowModelType::size);
    // Use a lambda here to prevent errors when using Clang/GCC
    flow_model_class.def("__len__",[](FlowModelType const& fm){return fm.size();});
    // NOTE: The export below gives 'ValueError: vector::reserve' at runtime when compiled using Clang
    // flow_model_class.def("__len__",&FlowModelType::size);
    // NOTE: The export below produces 'internal compiler error: in fold_convert_loc' with GCC
    // flow_model_class.def("__len__",(SizeType(FlowModelType::*)()const)&FlowModelType::size);
    flow_model_class.def("__getitem__",&__getitem__<FlowModelType,SizeType>);
    flow_model_class.def("__str__", &__cstr__<FlowModelType>);

    pybind11::class_<IntegratorInterface,IntegratorWrapper> integrator_interface_class(module,"IntegratorInterface");
    integrator_interface_class.def("flow_step",(FlowStepModelType(IntegratorInterface::*)(const ValidatedVectorMultivariateFunction&, const ExactBoxType&)const)&IntegratorInterface::flow_step);
    integrator_interface_class.def("flow_step",(FlowStepModelType(IntegratorInterface::*)(const ValidatedVectorMultivariateFunction&, const ExactBoxType&, const StepSizeType&)const)&IntegratorInterface::flow_step);
    integrator_interface_class.def("flow_step",(FlowStepModelType(IntegratorInterface::*)(const ValidatedVectorMultivariateFunction&,const ExactBoxType&,const StepSizeType&,const UpperBoxType&)const)&IntegratorInterface::flow_step);
    integrator_interface_class.def("__str__", &__cstr__<IntegratorInterface>);

    pybind11::class_<TaylorPicardIntegrator,IntegratorInterface> taylor_picard_integrator_class(module,"TaylorPicardIntegrator");
    taylor_picard_integrator_class.def(pybind11::init<ApproximateDouble>());
    taylor_picard_integrator_class.def("minimum_temporal_order",&TaylorPicardIntegrator::minimum_temporal_order);
    taylor_picard_integrator_class.def("maximum_temporal_order",&TaylorPicardIntegrator::maximum_temporal_order);
    taylor_picard_integrator_class.def("step_maximum_error",&TaylorPicardIntegrator::step_maximum_error);
    taylor_picard_integrator_class.def("set_minimum_temporal_order",&TaylorPicardIntegrator::set_minimum_temporal_order);
    taylor_picard_integrator_class.def("set_maximum_temporal_order",&TaylorPicardIntegrator::set_maximum_temporal_order);
    taylor_picard_integrator_class.def("set_step_maximum_error",&TaylorPicardIntegrator::set_step_maximum_error);

    pybind11::class_<GradedTaylorPicardIntegrator,IntegratorInterface> unbounded_taylor_picard_integrator_class(module, "GradedTaylorPicardIntegrator");
    unbounded_taylor_picard_integrator_class.def(pybind11::init<ApproximateDouble,Order>());
    unbounded_taylor_picard_integrator_class.def("order",&GradedTaylorPicardIntegrator::order);
    unbounded_taylor_picard_integrator_class.def("step_maximum_error",&GradedTaylorPicardIntegrator::step_maximum_error);
    unbounded_taylor_picard_integrator_class.def("error_refinement_minimum_improvement_percentage",&GradedTaylorPicardIntegrator::error_refinement_minimum_improvement_percentage);
    unbounded_taylor_picard_integrator_class.def("set_order",&GradedTaylorPicardIntegrator::set_order);
    unbounded_taylor_picard_integrator_class.def("set_step_maximum_error",&GradedTaylorPicardIntegrator::set_step_maximum_error);
    unbounded_taylor_picard_integrator_class.def("set_error_refinement_minimum_improvement_percentage",&GradedTaylorPicardIntegrator::set_error_refinement_minimum_improvement_percentage);

    pybind11::class_<TaylorSeriesIntegrator,IntegratorInterface> taylor_series_integrator_class(module,"TaylorSeriesIntegrator");
    taylor_series_integrator_class.def(pybind11::init<ApproximateDouble,Nat>());
    taylor_series_integrator_class.def("order",&TaylorSeriesIntegrator::order);
    taylor_series_integrator_class.def("has_step_maximum_error",&TaylorSeriesIntegrator::has_step_maximum_error);
    taylor_series_integrator_class.def("step_maximum_error",&TaylorSeriesIntegrator::step_maximum_error);
    taylor_series_integrator_class.def("set_order",&TaylorSeriesIntegrator::set_order);
    taylor_series_integrator_class.def("set_step_maximum_error",&TaylorSeriesIntegrator::set_step_maximum_error);

    pybind11::class_<GradedTaylorSeriesIntegrator,IntegratorInterface> graded_taylor_series_integrator_class(module,"GradedTaylorSeriesIntegrator");
    graded_taylor_series_integrator_class.def(pybind11::init<ApproximateDouble>());
    graded_taylor_series_integrator_class.def("maximum_spacial_order",&GradedTaylorSeriesIntegrator::maximum_spacial_order);
    graded_taylor_series_integrator_class.def("maximum_temporal_order",&GradedTaylorSeriesIntegrator::maximum_temporal_order);
    graded_taylor_series_integrator_class.def("step_maximum_error",&GradedTaylorSeriesIntegrator::step_maximum_error);
    graded_taylor_series_integrator_class.def("set_maximum_spacial_order",&GradedTaylorSeriesIntegrator::set_maximum_spacial_order);
    graded_taylor_series_integrator_class.def("set_maximum_temporal_order",&GradedTaylorSeriesIntegrator::set_maximum_temporal_order);
    graded_taylor_series_integrator_class.def("set_step_maximum_error",&GradedTaylorSeriesIntegrator::set_step_maximum_error);

    pybind11::class_<RungeKutta4Integrator> runge_kutta_4_integrator_class(module,"RungeKutta4Integrator");
    runge_kutta_4_integrator_class.def(pybind11::init<ApproximateDouble>());
    runge_kutta_4_integrator_class.def("step", &RungeKutta4Integrator::step);
    runge_kutta_4_integrator_class.def("evolve", &RungeKutta4Integrator::evolve);
}


Void integration_submodule(pybind11::module& module)
{
    export_integrators(module);
}
