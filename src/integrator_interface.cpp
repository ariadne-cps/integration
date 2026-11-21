/***************************************************************************
 *            src/integrator_interface.cpp
 *
 *  Copyright  2026  Luca Geretti
 ****************************************************************************/

/** @file src/integrator_interface.cpp
 *  @brief Non-template integration interface and error operations.
 */

#include "integration/integrator.hpp"

namespace Ariadne {

FlowBoundsException::FlowBoundsException(const StringType& what) : std::runtime_error(what) {}

FlowTimeStepException::FlowTimeStepException(const StringType& what) : std::runtime_error(what) {}

const FlowStepModelType& IncompleteFlowException::computed_model() const { return *this->_computed_model; }

OutputStream& operator<<(OutputStream& os, const IntegratorInterface& integrator) {
    integrator._write(os);
    return os;
}

} // namespace Ariadne
