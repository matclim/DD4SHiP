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
// Date       : 09.11.2025
//==========================================================================
#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Printout.h>
#include <iostream>
using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
  //Calo scintillator bars' feature extraction
  double       tol     = 0 * dd4hep::mm;
  xml_det_t    x_det   = e;
  xml_dim_t    x_detbox   = x_det.child(_U(box));
  xml_dim_t    x_rot   = x_det.child(_U(rotation));
  xml_dim_t    x_pos   = x_det.child(_U(position));
  xml_det_t    x_widebar = x_det.child(_Unicode(bar));
  xml_det_t    x_passive_layer = x_det.child(_Unicode(passive_layer));
  std::string  nam     = x_det.nameStr();
  //vertical bars by default
  const double widebar_x_spacing   =  x_widebar.attr<double>("x_extra_spacing");
  const double extrazgap   =  x_widebar.attr<double>("extrazgap");
  const std::string calo_layer_codes = x_det.attr<std::string>("layer_codes");
  const int num_z   =  static_cast<unsigned>(calo_layer_codes.size()); 
  const int widebar_num_x   =  x_widebar.attr<unsigned>("num_x");

  
  //Bar definition
  Box   widebar((x_widebar.x()-tol)/2., (x_widebar.y()-tol)/2.,(x_widebar.z()-tol)/2.);
  Box   passive_layer_box((x_passive_layer.x()-tol)/2., (x_passive_layer.y()-tol)/2.,(x_passive_layer.z()-tol)/2.);
  Volume widebar_vol("hcal_widebar", widebar, description.material(x_widebar.materialStr()));
  Volume passive_layer_vol("passive_layer", passive_layer_box, description.material(x_passive_layer.materialStr()));
  widebar_vol.setAttributes(description, x_widebar.regionStr(), x_widebar.limitsStr(), x_widebar.visStr());
  passive_layer_vol.setAttributes(description, x_passive_layer.regionStr(), x_passive_layer.limitsStr(), x_passive_layer.visStr());

  printout(INFO, "SandwichCalo", "%s: Bars: x: %7.3f y: %7.3f z: %7.3f mat: %s vis: %s solid: %s",
           nam.c_str(), x_widebar.x(), x_widebar.y(), x_widebar.z(), x_widebar.materialStr().c_str(),
           x_widebar.visStr().c_str(), widebar.type());
  sens.setType("calorimeter");

  // Envelope: make envelope box 'tol' bigger on each side
  Box    detbox((x_detbox.x()+tol)/2., (x_detbox.y()+tol)/2., (x_detbox.z()+tol)/2.);
  Volume detbox_vol(nam, detbox, description.air());
  detbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  
//  box_vol.setVisAttributes(description.visAttributes(""));

  Box    det_wide_layerbox((x_detbox.x()+tol)/2., (x_detbox.y()+tol)/2., (x_widebar.z()+tol)/2.);
  Volume det_wide_layerbox_vol("hcal_det_wide_layerbox", det_wide_layerbox, description.air());
  det_wide_layerbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  det_wide_layerbox_vol.setVisAttributes(description.visAttributes(x_detbox.visStr()));
  
  printout(INFO, "HCAL", "%s: Layer:   nx: %7d x spacing: %7.3f", nam.c_str(), widebar_num_x, widebar_x_spacing);
//  Rotation3D rot(RotationZYX(0e0, 0e0, M_PI/2e0));
  Rotation3D rot(RotationZYX(0e0, 0e0, 0e0));
  
   // sens.setType("calorimeter");
    widebar_vol.setSensitiveDetector(sens);
  //}

  //Loop for x-wise placement -> build the sensitive bar layer 
  //
  //Build Wide bar layers

 // long DetectorCode = 9 * 1e15;
 // long HCALCode = 2 * 1e12;

  //int DetectorCode = 9 * 1e8;
  //int HCALCode = 2 * 1e7;
  
  double wideboxwidth = x_widebar.x()*static_cast<double>(widebar_num_x);
  double xpos = -(wideboxwidth+tol)/2.;
   
  for( int ix=0; ix < widebar_num_x; ++ix )  {
    xpos += x_widebar.x()/2.;  
    PlacedVolume pv = det_wide_layerbox_vol.placeVolume(widebar_vol, Transform3D(rot,Position(xpos, 0e0, 0e0)));
    pv.addPhysVolID("widebar", ix);
    xpos += x_widebar.x()/2. + widebar_x_spacing; 
  }



  //Loop for z-wide placement -> build the calorimeter sandwich
  double z_layer = -x_detbox.z()/2.;
  Rotation3D rot_layers;



  for( int iz=0; iz < num_z; ++iz )  {
    // leave 'tol' space between the layers

    
    std::cout << static_cast<int>(calo_layer_codes[iz]) - '0' << std::endl;
    switch(static_cast<int>(calo_layer_codes[iz]) - '0'){
    
	    case 1:{
		//Place wide layer vertically

    		z_layer += x_widebar.z()/2.;
	    	rot_layers = RotationZYX(M_PI/2e0,0e0,0e0);
    	    	PlacedVolume pv_det = detbox_vol.placeVolume(det_wide_layerbox_vol, Transform3D(rot_layers,Position(0., 0., z_layer)));
    	    	pv_det.addPhysVolID("hcal_layer", iz);
    		//z_layer += x_widebar.z()+x_passive_layer.z();
    		z_layer += x_widebar.z()/2.;
   		break;
	    }
	    case 2:{
		//Place wide layer horizontally	
		
		z_layer += x_widebar.z()/2.;
		rot_layers = RotationZYX(0e0, 0e0, 0e0);
    		PlacedVolume pv_det = detbox_vol.placeVolume(det_wide_layerbox_vol, Transform3D(rot_layers,Position(0., 0., z_layer)));
    	    	pv_det.addPhysVolID("hcal_layer", iz);
    		//z_layer += x_widebar.z()+x_passive_layer.z();
    		z_layer += x_widebar.z()/2.;
	  	break; 
	    }
	    case 7:{
    		z_layer += x_passive_layer.z()/2.;
    		//PlacedVolume pv_passive = detbox_vol.placeVolume(passive_layer_vol,Transform3D(rot_layers,Position(0e0, 0e0, z_passive)) );
    		PlacedVolume pv_passive = detbox_vol.placeVolume(passive_layer_vol,Transform3D(rot_layers,Position(0., 0. , z_layer)));
    		z_layer += x_passive_layer.z()/2.;
    		pv_passive.addPhysVolID("hcal_passivelayer", iz);
    		z_layer += extrazgap;
	  	break; 
	    }
    }
    std::cout << "Zlayer Det " << z_layer << std::endl;

    if(static_cast<int>(calo_layer_codes[iz]) - '0' != 5 && static_cast<int>(calo_layer_codes[iz]) - '0' != 6){
    }
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
  printout(INFO, "SHiP HCAL", "%s: Detector construction finished.", nam.c_str());
  return sdet;
}

DECLARE_DETELEMENT(DD4hep_SHiPHCAL,create_detector)
