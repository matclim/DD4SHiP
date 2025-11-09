//==========================================================================
//  AIDA Detector description implementation 
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================
//
//
// Display using:
// $> geoDisplay examples/ClientTests/compact/SHiP_HPL_Fibre_Trackers.xml
//
//==========================================================================
#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Printout.h>

using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
  double       tol     = 1e-5 * dd4hep::mm;
  xml_det_t    x_det   = e;
  xml_dim_t    x_box   = x_det.child(_U(box));
  xml_dim_t    x_rot   = x_det.child(_U(rotation));
  xml_dim_t    x_pos   = x_det.child(_U(position));
  xml_det_t    x_fibre = x_det.child(_Unicode(fibre));
  xml_det_t    x_core   = x_det.child(_Unicode(core));
  std::string  nam     = x_det.nameStr();
  const double thick   = x_fibre.thickness();
  const double delta   = 2e0*x_fibre.rmax();
  const int    num_x   = int(2e0*x_box.x() / delta);
  const int    num_x_small = num_x - 1;
//  const int    num_z   = int(2e0*x_box.z() / (delta+2*tol));
  const double num_z   =  x_det.attr<int>("n_fibre_layers");
  Tube   fibre(0., x_fibre.rmax()-tol, x_fibre.y()-tol);
  Volume fibre_vol("fibre", fibre, description.material(x_fibre.materialStr()));
  fibre_vol.setAttributes(description, x_fibre.regionStr(), x_fibre.limitsStr(), x_fibre.visStr());
  
  Tube   fibre_core(0., fibre.rMax()-thick, fibre.dZ());
  Volume fibre_core_vol("core", fibre_core, description.material(x_core.materialStr()));
  fibre_core_vol.setAttributes(description, x_core.regionStr(), x_core.limitsStr(), x_core.visStr());

  fibre_vol.placeVolume(fibre_core_vol);

  printout(INFO, "SHiP_HPL_Fibre_Trackers", "%s: Straw: rmax: %7.3f y: %7.3f mat: %s vis: %s solid: %s",
           nam.c_str(), x_fibre.rmax(), x_fibre.y(), x_fibre.materialStr().c_str(),
           x_fibre.visStr().c_str(), fibre.type());
  if( x_core.hasChild(_U(sensitive)) )  {
    sens.setType("tracker");
    fibre_core_vol.setSensitiveDetector(sens);
  }

  // Envelope: make envelope box 'tol' bigger on each side, so that the fibres
  Box    box(x_box.x()+tol, x_box.y()+tol, x_box.z()+tol);
  Volume box_vol(nam, box, description.air());
  box_vol.setAttributes(description, x_box.regionStr(), x_box.limitsStr(), x_box.visStr());

  Box    big_layer(x_box.x(), x_box.y(), x_fibre.rmax());
  Volume big_layer_vol("big_layer", big_layer, description.air());
  big_layer_vol.setVisAttributes(description.visAttributes("VisibleGray"));
  
  Box    small_layer(x_box.x(), x_box.y(), x_fibre.rmax());
  Volume small_layer_vol("small_layer", small_layer, description.air());
  small_layer_vol.setVisAttributes(description.visAttributes("VisibleGray"));
  
  printout(INFO, "SHiP_HPL_Fibre_Trackers", "%s: Layer:   nx: %7d nz: %7d delta: %7.3f", nam.c_str(), num_x, num_z, delta);
 //Big layer creation
  Rotation3D rot(RotationZYX(0e0, 0e0, M_PI/2e0));
  for( int ix=0; ix < num_x; ++ix )  {
    double x = -box.x() + (double(ix)+0.5) * (delta + 2e0*tol);
    PlacedVolume pv = big_layer_vol.placeVolume(fibre_vol, Transform3D(rot,Position(x, 0e0, 0e0)));
    pv.addPhysVolID("fibre", ix);
  }
  
  for( int ix=0; ix < num_x_small; ++ix )  {
    double x = -box.x() + (double(ix)+0.5) * (delta + 2e0*tol) + x_fibre.rmax();
    PlacedVolume pv = small_layer_vol.placeVolume(fibre_vol, Transform3D(rot,Position(x, 0e0, 0e0)));
    pv.addPhysVolID("fibre", ix);
  }


  for( int iz=0; iz < num_z; ++iz )  {
    // leave 'tol' space between the layers
    if(iz%2 == 0){
    	double z = -box.z() + (double(iz)+0.5) * (2.0*tol + delta);
    	PlacedVolume pv = box_vol.placeVolume(big_layer_vol, Position(0e0, 0e0, z));
    	pv.addPhysVolID("big_layer", iz);
    }
    else{
    	double z = -box.z() + (double(iz)+0.5) * (2.0*tol + delta);
    	PlacedVolume pv = box_vol.placeVolume(small_layer_vol, Position(0e0, 0e0, z));
    	pv.addPhysVolID("small_layer", iz);
    }
  }
  printout(INFO, "SHiP_HPL_Fibre_Trackers", "%s: Created %d layers of %d fibres each.", nam.c_str(), num_z, num_x);
  
  DetElement   sdet  (nam, x_det.id());
  Volume       mother(description.pickMotherVolume(sdet));
  Rotation3D   rot3D (RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0)));
  Transform3D  trafo (rot3D, Position(x_pos.x(0), x_pos.y(0), x_pos.z(0)));
  PlacedVolume pv = mother.placeVolume(box_vol, trafo);
  pv.addPhysVolID("system", x_det.id());
  sdet.setPlacement(pv);  // associate the placed volume to the detector element
  printout(INFO, "SHiP_HPL_Fibre_Trackers", "%s: Detector construction finished.", nam.c_str());
  return sdet;
}
DECLARE_DETELEMENT(DD4hep_SHiP_HPL_Fibre_Tracker,create_detector)
