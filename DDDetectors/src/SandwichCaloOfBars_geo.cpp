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
  double       tol     = 0 * dd4hep::mm;
  xml_det_t    x_det   = e;
  xml_dim_t    x_detbox   = x_det.child(_U(box));
  xml_dim_t    x_rot   = x_det.child(_U(rotation));
  xml_dim_t    x_pos   = x_det.child(_U(position));
  xml_det_t    x_bar = x_det.child(_Unicode(bar));
  xml_det_t    x_passive_layer = x_det.child(_Unicode(passive_layer));
  xml_det_t    x_split = x_det.child(_Unicode(split));
  std::string  nam     = x_det.nameStr();
  //vertical bars by default
  const double x_spacing   =  x_bar.attr<double>("x_spacing");
  const double z_spacing   = 2*(x_bar.z()+x_passive_layer.z());  
  const double splitlayer   =  x_bar.attr<int>("splitlayer");
  const double extrazgap   =  x_bar.attr<double>("extrazgap");
  const int num_z   =  x_bar.attr<int>("num_z");
  const int num_x   =  x_bar.attr<int>("num_x");
  Box   bar(x_bar.x()-tol, x_bar.y()-tol,(x_bar.z()-tol)/2.);
  Box   passive_layer_box(x_passive_layer.x()-tol, x_passive_layer.y()-tol,(x_passive_layer.z()-tol)/2.);
  Volume bar_vol("bar", bar, description.material(x_bar.materialStr()));
  Volume passive_layer_vol("passive_layer", passive_layer_box, description.material(x_passive_layer.materialStr()));
  bar_vol.setAttributes(description, x_bar.regionStr(), x_bar.limitsStr(), x_bar.visStr());
  passive_layer_vol.setAttributes(description, x_passive_layer.regionStr(), x_passive_layer.limitsStr(), x_passive_layer.visStr());
  

  printout(INFO, "SandwichCalo", "%s: Bars: x: %7.3f y: %7.3f z: %7.3f mat: %s vis: %s solid: %s",
           nam.c_str(), x_bar.x(), x_bar.y(), x_bar.z(), x_bar.materialStr().c_str(),
           x_bar.visStr().c_str(), bar.type());
  sens.setType("calorimeter");
  bar_vol.setSensitiveDetector(sens);

  // Envelope: make envelope box 'tol' bigger on each side
  Box    detbox(x_detbox.x()+tol, x_detbox.y()+tol, (x_detbox.z()+tol)/2.);
  Volume detbox_vol(nam, detbox, description.air());
  detbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  
//  box_vol.setVisAttributes(description.visAttributes(""));

  Box    det_layerbox(x_detbox.x()+tol, x_detbox.y()+tol, (x_bar.z()+tol)/2.);
  Volume det_layerbox_vol("det_layerbox", det_layerbox, description.air());
  det_layerbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  det_layerbox_vol.setVisAttributes(description.visAttributes(x_detbox.visStr()));
  
  printout(INFO, "SandwichCalo", "%s: Layer:   nx: %7d x spacing: %7.3f", nam.c_str(), num_x, x_spacing);
//  Rotation3D rot(RotationZYX(0e0, 0e0, M_PI/2e0));
  Rotation3D rot(RotationZYX(0e0, 0e0, 0e0));
  
  if( x_bar.hasChild(_U(sensitive)) )  {
    sens.setType("calorimeter");
    bar_vol.setSensitiveDetector(sens);
  }

  
  //Loop for x-wise placement -> build the sensitive bar layer 
  for( int ix=0; ix < num_x; ++ix )  {
    double x = x_bar.x()/2. + static_cast<double>(ix) *  x_spacing; 
    std::cout << "X value of placed bar " << x << std::endl;
    PlacedVolume pv = det_layerbox_vol.placeVolume(bar_vol, Transform3D(rot,Position(x, 0e0, 0e0)));
    pv.addPhysVolID("bar", ix);
  }
  //Loop for z-wide placement -> build the calorimeter sandwich
  double z_layer =0.;
  Rotation3D rot_layers;
  for( int iz=0; iz < num_z; ++iz )  {
    // leave 'tol' space between the layers
    //z_layer += x_passive_layer.z()+x_bar.z();
    if(iz%2==1) {
    	    z_layer += x_bar.z()/2.;
	    rot_layers = RotationZYX(M_PI/2e0,0e0,0e0);
    	    PlacedVolume pv_det = detbox_vol.placeVolume(det_layerbox_vol, Transform3D(rot_layers,Position(x_bar.y(),-x_bar.y() , z_layer)));
    	    pv_det.addPhysVolID("layer", iz*2);
    	    z_layer += x_bar.z()/2.;
    
    
    }
    else{ 
    	z_layer += x_bar.z()/2.;
	rot_layers = RotationZYX(0e0, 0e0, 0e0);
    	PlacedVolume pv_det = detbox_vol.placeVolume(det_layerbox_vol, Transform3D(rot_layers,Position(0e0, 0e0, z_layer)));
        pv_det.addPhysVolID("layer", iz*2);
    	z_layer += x_bar.z()/2.;

    }
    std::cout << "Zlayer Det " << z_layer << " DET Z " << x_bar.z() << std::endl;
    //PlacedVolume pv_passive = detbox_vol.placeVolume(passive_layer_vol,Transform3D(rot_layers,Position(0e0, 0e0, z_passive)) );
    //z_layer += x_bar.z()+x_passive_layer.z();
    z_layer += x_passive_layer.z()/2.;
    PlacedVolume pv_passive = detbox_vol.placeVolume(passive_layer_vol,Transform3D(rot_layers,Position(x_passive_layer.x(), 0e0, z_layer)));
    std::cout << "Zlayer passive " << z_layer << " passive Z "<< x_passive_layer.z() << std::endl;
    pv_passive.addPhysVolID("passivelayer", iz*2+1);
    z_layer += x_passive_layer.z()/2.;
    z_layer += extrazgap;
  }
//  printout(INFO, "SandwichCalo", "%s: Created %d layers of %d bars each.", nam.c_str(), num_z, num_x);
  //PlacedVolume pv2 = detbox_vol.placeVolume(det_layerbox_vol, Transform3D(rot,Position(0e0, 0e0, 0e0)));
  //pv2.addPhysVolID("det_layerbox", 0e0);
  
  DetElement   sdet  (nam, x_det.id());
  Volume       mother(description.pickMotherVolume(sdet));
  Rotation3D   rot3D (RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0)));
  Transform3D  trafo (rot3D, Position(x_pos.x(0), x_pos.y(0), x_pos.z(0)));
// PlacedVolume pv2 = mother.placeVolume(passive_layer_vol, trafo);
  PlacedVolume pv = mother.placeVolume(detbox_vol, trafo);
  //pv2.addPhysVolID("system", x_det.id());
  pv.addPhysVolID("system", x_det.id());
  //sdet.setPlacement(pv2);  // associate the placed volume to the detector element
  sdet.setPlacement(pv);  // associate the placed volume to the detector element
  printout(INFO, "SandwichCalo", "%s: Detector construction finished.", nam.c_str());
  return sdet;
}

DECLARE_DETELEMENT(DD4hep_SandwichCalo,create_detector)
