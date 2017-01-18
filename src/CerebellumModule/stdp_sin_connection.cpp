/*
 *  stdp_sin_connection.cpp
 *
 */

#include "network.h"
#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "stdp_sin_connection.h"
#include "event.h"
#include "nestmodule.h"

namespace nest
{
//
// Implementation of class STDPSinCommonProperties.
//

STDPSinCommonProperties::STDPSinCommonProperties()
  : CommonSynapseProperties()
  , A_plus_( 1.0 )
  , A_minus_( 1.0 )
  , Peak_( 100.0 )
  , Exponent_( 2 )
  , Wmin_( 0.0 )
  , Wmax_( 200.0 )
{
  inv_tau_ = atan((float) Exponent_)/Peak_;
  factor_ = 1.0f/(exp(-atan((float)Exponent_))*pow(sin(atan((float)Exponent_)),(int) Exponent_));
}

void
STDPSinCommonProperties::get_status( DictionaryDatum& d ) const
{
  CommonSynapseProperties::get_status( d );

  def< double_t >( d, "A_plus", A_plus_ );
  def< double_t >( d, "A_minus", A_minus_ );
  def< double_t >( d, "peak", Peak_ );
  def< int_t >( d, "exponent", Exponent_ );
  def< double_t >( d, "Wmin", Wmin_ );
  def< double_t >( d, "Wmax", Wmax_ );
}

void
STDPSinCommonProperties::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  CommonSynapseProperties::set_status( d, cm );

  updateValue< double_t >( d, "A_plus", A_plus_ );
  updateValue< double_t >( d, "A_minus", A_minus_ );
  updateValue< int_t >( d, "exponent", Exponent_ );
  updateValue< double_t >( d, "peak", Peak_ );
  updateValue< double_t >( d, "Wmin", Wmin_ );
  updateValue< double_t >( d, "Wmax", Wmax_ );

  inv_tau_ = atan((float) Exponent_)/Peak_;
  factor_ = 1.0f/(exp(-atan((float)Exponent_))*pow(sin(atan((float)Exponent_)),(int) Exponent_)); 
}

} // of namespace nest
