/*
 *  stdp_sin_connection.cpp
 *
 */

#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "stdp_sin_connection.h"
#include "event.h"
#include "nestmodule.h"

namespace mynest
{
//
// Implementation of class STDPSinCommonProperties.
//
float STDPSinCommonProperties::terms[11][11] = {{1,0,0,0,0,0,0,0,0,0,0},
    {A,-A,0,0,0,0,0,0,0,0,0},
    {3.0f/2.0f*pow(A,2),-4.0f/2.0f*pow(A,2),1.0f/2.0f*pow(A,2),0,0,0,0,0,0,0,0},
    {10.0f/4.0f*pow(A,3),-15.0f/4.0f*pow(A,3),6.0f/4.0f*pow(A,3),-1.0f/4.0f*pow(A,3),0,0,0,0,0,0,0},
    {35.0f/8.0f*pow(A,4),-56.0f/8.0f*pow(A,4),28.0f/8.0f*pow(A,4),-8.0f/8.0f*pow(A,4),1.0f/8.0f*pow(A,4),0,0,0,0,0,0},
    {126.0f/16.0f*pow(A,5),-210.0f/16.0f*pow(A,5),120.0f/16.0f*pow(A,5),-45.0f/16.0f*pow(A,5),10.0f/16.0f*pow(A,5),-1.0f/16.0f*pow(A,5),0,0,0,0,0},
    {231.0f/16.0f*pow(A,6),-99.0f/4.0f*pow(A,6),495.0f/32.0f*pow(A,6),-55.0f/8.0f*pow(A,6),66.0f/32.0f*pow(A,6),-3.0f/8.0f*pow(A,6),1.0f/32.0f*pow(A,6),0,0,0,0},
    {429.0f/16.0f*pow(A,7),-3003.0f/64.0f*pow(A,7),1001.0f/32.0f*pow(A,7),-1001.0f/64.0f*pow(A,7),91.0f/16.0f*pow(A,7),-91.0f/64.0f*pow(A,7),7.0f/32.0f*pow(A,7),-1.0f/64.0f*pow(A,7),0,0,0},
    {6435.0f/128.0f*pow(A,8),-715.0f/8.0f*pow(A,8),1001.0f/16.0f*pow(A,8),-273.0f/8.0f*pow(A,8),455.0f/32.0f*pow(A,8),-35.0f/8.0f*pow(A,8),15.0f/16.0f*pow(A,8),-1.0f/8.0f*pow(A,8),1.0f/128.0f*pow(A,8),0,0},
    {12155.0f/128.0f*pow(A,9),-21879.0f/128.0f*pow(A,9),1989.0f/16.0f*pow(A,9),-4641.0f/64.0f*pow(A,9),1071.0f/32.0f*pow(A,9),-765.0f/64.0f*pow(A,9),51.0f/16.0f*pow(A,9),-153.0f/256.0f*pow(A,9),9.0f/128.0f*pow(A,9),-1.0f/256.0f*pow(A,9),0},
    {46189.0f/256.0f*pow(A,10),-20995.0f/64.0f*pow(A,10),62985.0f/256.0f*pow(A,10),-4845.0f/32.0f*pow(A,10),4845.0f/64.0f*pow(A,10),-969.0f/32.0f*pow(A,10),4845.0f/512.0f*pow(A,10),-285.0f/128.0f*pow(A,10),95.0f/256.0f*pow(A,10),-5.0f/128.0f*pow(A,10),1.0f/512.0f*pow(A,10)}};


STDPSinCommonProperties::STDPSinCommonProperties()
  : CommonSynapseProperties()
  , vt_( 0 )
  , A_plus_( 1.0 )
  , A_minus_( 1.0 )
  , Wmin_( 0.0 )
  , Wmax_( 200.0 )
  
{
  
}

void
STDPSinCommonProperties::get_status( DictionaryDatum& d ) const
{
  CommonSynapseProperties::get_status( d );

  if ( vt_ != 0 )
    def< long >( d, "vt", vt_->get_gid() );
  else
    def< long >( d, "vt", -1 );

  def< float >( d, "A_plus", A_plus_ );
  def< float >( d, "A_minus", A_minus_ );
  def< float >( d, "Wmin", Wmin_ );
  def< float >( d, "Wmax", Wmax_ );
}

void
STDPSinCommonProperties::set_status( const DictionaryDatum& d, nest::ConnectorModel& cm )
{
  CommonSynapseProperties::set_status( d, cm );

  long vtgid;
  if ( updateValue< long >( d, "vt", vtgid ) )
  {
    vt_ = dynamic_cast< volume_transmitter* >( nest::kernel().node_manager.get_node(vtgid, nest::kernel().vp_manager.get_thread_id() ) );

    if ( vt_ == 0 )
      throw nest::BadProperty( "STDP sin source must be volume transmitter" );
  }

  updateValue< float >( d, "A_plus", A_plus_ );
  updateValue< float >( d, "A_minus", A_minus_ );
  
  updateValue< float >( d, "Wmin", Wmin_ );
  updateValue< float >( d, "Wmax", Wmax_ );
}

nest::Node*
STDPSinCommonProperties::get_node()
{
  if ( vt_ == 0 )
    throw nest::BadProperty( "No volume transmitter has been assigned to the STDP Sin synapse." );
  else
    return vt_;
}

} // of namespace nest
