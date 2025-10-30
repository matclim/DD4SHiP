//==========================================================================
//  AIDA Detector description implementation 
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Climescu
//
//==========================================================================
#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Printout.h>
#include <iostream>
using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
  double       tol     = 1e-5 * dd4hep::mm;
  xml_det_t    x_det   = e;
  xml_dim_t    x_detbox   = x_det.child(_U(box));
  xml_dim_t    x_rot   = x_det.child(_U(rotation));
  xml_dim_t    x_pos   = x_det.child(_U(position));
  xml_det_t    x_bar = x_det.child(_Unicode(bar));
  std::string  nam     = x_det.nameStr();
  //vertical bars by default
  const double spacing   =  x_bar.attr<int>("spacing");
  //const int    num_x   =static_cast<int>(x_bar.y()/x_bar.x());
//  const int    num_x   = x_bar.num_x();
//  const int    num_x   = x_det.child(_Unicode(num_x))
  const int num_x = x_bar.attr<int>("num_x");
  Box   bar(x_bar.x()-tol, x_bar.y()-tol,x_bar.z()-tol);
  Volume bar_vol("bar", bar, description.material(x_bar.materialStr()));
  bar_vol.setAttributes(description, x_bar.regionStr(), x_bar.limitsStr(), x_bar.visStr());
  

  printout(INFO, "LayerOfBars", "%s: Bars: x: %7.3f y: %7.3f z: %7.3f mat: %s vis: %s solid: %s",
           nam.c_str(), x_bar.x(), x_bar.y(), x_bar.z(), x_bar.materialStr().c_str(),
           x_bar.visStr().c_str(), bar.type());
  sens.setType("calorimeter");
  bar_vol.setSensitiveDetector(sens);

  // Envelope: make envelope box 'tol' bigger on each side
  Box    detbox(x_detbox.x()+tol, x_detbox.y()+tol, x_detbox.z()+tol);
  Volume detbox_vol(nam, detbox, description.air());
  detbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
//  box_vol.setVisAttributes(description.visAttributes(""));

  Box    layerbox(x_detbox.x()+tol, x_detbox.y()+tol, x_bar.z()+tol);
  Volume layerbox_vol("layerbox", layerbox, description.air());
  layerbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  layerbox_vol.setVisAttributes(description.visAttributes(x_detbox.visStr()));
  
  printout(INFO, "LayerOfBars", "%s: Layer:   nx: %7d spacing: %7.3f", nam.c_str(), num_x, spacing);
  Rotation3D rot(RotationZYX(0e0, 0e0, M_PI/2e0));
  for( int ix=0; ix < num_x; ++ix )  {
//    double x = -box.x() + (double(ix)+0.5) * (spacing + 2e0*tol);
    double x = x_bar.x()/2. + static_cast<double>(ix) *  spacing; 
  //  printout(INFO, "LayerOfBars", "%s: PLACED Layer:   x: %7.3f ", x);
    std::cout << "X value of placed bar " << x << std::endl;
    PlacedVolume pv = layerbox_vol.placeVolume(bar_vol, Transform3D(rot,Position(x, 0e0, 0e0)));
    pv.addPhysVolID("bar", ix);
  }
  //for( int iz=0; iz < num_z; ++iz )  {
  //  // leave 'tol' space between the layers
  //  double z = -box.z() + (double(iz)+0.5) * (2.0*tol + delta);
  //  PlacedVolume pv = box_vol.placeVolume(layer_vol, Position(0e0, 0e0, z));
  //  pv.addPhysVolID("layer", iz);
  //}
//  printout(INFO, "LayerOfBars", "%s: Created %d layers of %d bars each.", nam.c_str(), num_z, num_x);
  PlacedVolume pv2 = detbox_vol.placeVolume(layerbox_vol, Transform3D(rot,Position(0e0, 0e0, 0e0)));
  pv2.addPhysVolID("layerbox", 0e0);
  
  DetElement   sdet  (nam, x_det.id());
  Volume       mother(description.pickMotherVolume(sdet));
  Rotation3D   rot3D (RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0)));
  Transform3D  trafo (rot3D, Position(x_pos.x(0), x_pos.y(0), x_pos.z(0)));
  PlacedVolume pv = mother.placeVolume(detbox_vol, trafo);
  pv.addPhysVolID("system", x_det.id());
  sdet.setPlacement(pv);  // associate the placed volume to the detector element
  printout(INFO, "LayerOfBars", "%s: Detector construction finished.", nam.c_str());
  return sdet;
}

DECLARE_DETELEMENT(DD4hep_LayerOfBars,create_detector)
