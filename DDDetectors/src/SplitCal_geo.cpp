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
  double       tol     = 0 * dd4hep::mm;
  xml_det_t    x_det   = e;
  xml_dim_t    x_detbox   = x_det.child(_U(box));
  xml_dim_t    x_rot   = x_det.child(_U(rotation));
  xml_dim_t    x_pos   = x_det.child(_U(position));
  xml_det_t    x_widebar = x_det.child(_Unicode(widebar));
  xml_det_t    x_thinbar = x_det.child(_Unicode(thinbar));
  xml_det_t    x_passive_layer = x_det.child(_Unicode(passive_layer));
  xml_det_t    x_split = x_det.child(_Unicode(split));
  std::string  nam     = x_det.nameStr();
  //vertical bars by default
  const double splitlayer   =  x_det.attr<int>("splitlayer");
  const double widebar_x_spacing   =  x_widebar.attr<double>("x_spacing");
  const double thinbar_x_spacing   =  x_thinbar.attr<double>("x_spacing");
  const double extrazgap   =  x_widebar.attr<double>("extrazgap");
  const std::string calo_layer_codes = x_det.attr<std::string>("layer_codes");
  const int num_z   =  static_cast<unsigned>(calo_layer_codes.size()); 
  const int widebar_num_x   =  x_widebar.attr<unsigned>("num_x");
  const int thinbar_num_x   =  x_thinbar.attr<unsigned>("num_x");
  Box   widebar(x_widebar.x()-tol, x_widebar.y()-tol,(x_widebar.z()-tol)/2.);
  Box   thinbar(x_thinbar.x()-tol, x_thinbar.y()-tol,(x_thinbar.z()-tol)/2.);
  Box   passive_layer_box(x_passive_layer.x()-tol, x_passive_layer.y()-tol,(x_passive_layer.z()-tol)/2.);
  Box   split_box(x_split.x()-tol, x_split.y()-tol,(x_split.z()-tol)/2.);
  Volume widebar_vol("widebar", widebar, description.material(x_widebar.materialStr()));
  Volume thinbar_vol("thinbar", thinbar, description.material(x_thinbar.materialStr()));
  Volume passive_layer_vol("passive_layer", passive_layer_box, description.material(x_passive_layer.materialStr()));
  Volume split_vol("split", split_box, description.material(x_split.materialStr()));
  widebar_vol.setAttributes(description, x_widebar.regionStr(), x_widebar.limitsStr(), x_widebar.visStr());
  thinbar_vol.setAttributes(description, x_thinbar.regionStr(), x_thinbar.limitsStr(), x_thinbar.visStr());
  passive_layer_vol.setAttributes(description, x_passive_layer.regionStr(), x_passive_layer.limitsStr(), x_passive_layer.visStr());
  split_vol.setAttributes(description, x_split.regionStr(), x_split.limitsStr(), x_split.visStr());
 

  printout(INFO, "SandwichCalo", "%s: Bars: x: %7.3f y: %7.3f z: %7.3f mat: %s vis: %s solid: %s",
           nam.c_str(), x_widebar.x(), x_widebar.y(), x_widebar.z(), x_widebar.materialStr().c_str(),
           x_widebar.visStr().c_str(), widebar.type());
  sens.setType("calorimeter");
  widebar_vol.setSensitiveDetector(sens);
  thinbar_vol.setSensitiveDetector(sens);

  // Envelope: make envelope box 'tol' bigger on each side
  Box    detbox(x_detbox.x()+tol, x_detbox.y()+tol, x_detbox.z()+tol);
  Volume detbox_vol(nam, detbox, description.air());
  detbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  
//  box_vol.setVisAttributes(description.visAttributes(""));

  Box    det_wide_layerbox(x_detbox.x()+tol, x_detbox.y()+tol, x_widebar.z()+tol);
  Volume det_wide_layerbox_vol("det_wide_layerbox", det_wide_layerbox, description.air());
  det_wide_layerbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  det_wide_layerbox_vol.setVisAttributes(description.visAttributes(x_detbox.visStr()));
  
  Box    det_thin_layerbox(x_detbox.x()+tol, x_detbox.y()+tol, x_widebar.z()+tol);
  Volume det_thin_layerbox_vol("det_thin_layerbox", det_thin_layerbox, description.air());
  det_thin_layerbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  det_thin_layerbox_vol.setVisAttributes(description.visAttributes(x_detbox.visStr()));
  
  printout(INFO, "SandwichCalo", "%s: Layer:   nx: %7d x spacing: %7.3f", nam.c_str(), widebar_num_x, widebar_x_spacing);
//  Rotation3D rot(RotationZYX(0e0, 0e0, M_PI/2e0));
  Rotation3D rot(RotationZYX(0e0, 0e0, 0e0));
  
  if( x_widebar.hasChild(_U(sensitive)) )  {
    sens.setType("calorimeter");
    widebar_vol.setSensitiveDetector(sens);
  }

  if( x_thinbar.hasChild(_U(sensitive)) )  {
    sens.setType("calorimeter");
    thinbar_vol.setSensitiveDetector(sens);
  }
  
  //Loop for x-wise placement -> build the sensitive bar layer 
  //
  //Wide bar layers
  for( int ix=0; ix < widebar_num_x; ++ix )  {
    double x = x_widebar.x()/2. + static_cast<double>(ix) *  widebar_x_spacing; 
    std::cout << "X value of placed wide bar " << x << std::endl;
    PlacedVolume pv = det_wide_layerbox_vol.placeVolume(widebar_vol, Transform3D(rot,Position(x, 0e0, 0e0)));
    pv.addPhysVolID("widebar", ix);
  }
  //Thin bar layers
  for( int ix=0; ix < thinbar_num_x; ++ix )  {
    double x = x_thinbar.x()/2. + static_cast<double>(ix) *  thinbar_x_spacing; 
    std::cout << "X value of placed thin bar " << x << std::endl;
    PlacedVolume pv = det_thin_layerbox_vol.placeVolume(thinbar_vol, Transform3D(rot,Position(x, 0e0, 0e0)));
    pv.addPhysVolID("thinbar", ix);
  }
  //Loop for z-wide placement -> build the calorimeter sandwich
  double z_layer =0.;
  Rotation3D rot_layers;
  for( int iz=0; iz < num_z; ++iz )  {
    // leave 'tol' space between the layers

    
    std::cout << static_cast<int>(calo_layer_codes[iz]) - '0' << std::endl;
    switch(static_cast<int>(calo_layer_codes[iz]) - '0'){
    
	    case 1:{
		//Place wide layer vertically
    		z_layer += x_widebar.z()/2.;
	    	rot_layers = RotationZYX(M_PI/2e0,0e0,0e0);
    	    	PlacedVolume pv_det = detbox_vol.placeVolume(det_wide_layerbox_vol, Transform3D(rot_layers,Position(x_widebar.y(),-x_widebar.y() , z_layer)));
    	    	pv_det.addPhysVolID("layer", iz*2);
    		//z_layer += x_widebar.z()+x_passive_layer.z();
    		z_layer += x_widebar.z()/2.;
   		break;
	    }
	    case 2:{
		//Place wide layer horizontally	
    		z_layer += x_widebar.z()/2.;
		rot_layers = RotationZYX(0e0, 0e0, 0e0);
    		PlacedVolume pv_det = detbox_vol.placeVolume(det_wide_layerbox_vol, Transform3D(rot_layers,Position(0e0, 0e0, z_layer)));
        	pv_det.addPhysVolID("layer", iz*2);
    		//z_layer += x_widebar.z()+x_passive_layer.z();
    		z_layer += x_widebar.z()/2.;
	  	break; 
	    }
	    case 3:{
		//Place thin layer vertically
    		z_layer += x_thinbar.z()/2.;
	    	rot_layers = RotationZYX(M_PI/2e0,0e0,0e0);
    	    	PlacedVolume pv_det = detbox_vol.placeVolume(det_thin_layerbox_vol, Transform3D(rot_layers,Position(x_thinbar.y(),-x_thinbar.y() , z_layer)));
    	    	pv_det.addPhysVolID("layer", iz*2);
    		z_layer += x_thinbar.z()/2.;
    		//z_layer += x_thinbar.z()+x_passive_layer.z();
   		break;
            }
	    case 4:{
		//Place thin layer horizontally	
    		z_layer += x_thinbar.z()/2.;
		rot_layers = RotationZYX(0e0, 0e0, 0e0);
    		PlacedVolume pv_det = detbox_vol.placeVolume(det_thin_layerbox_vol, Transform3D(rot_layers,Position(0e0, 0e0, z_layer)));
        	pv_det.addPhysVolID("layer", iz*2);
    		z_layer += x_thinbar.z()/2.;
    		//z_layer += x_thinbar.z()+x_passive_layer.z();
		break;
            }
    }
    std::cout << "Zlayer Det " << z_layer << std::endl;

    z_layer += x_passive_layer.z()/2.;
    //PlacedVolume pv_passive = detbox_vol.placeVolume(passive_layer_vol,Transform3D(rot_layers,Position(0e0, 0e0, z_passive)) );
    PlacedVolume pv_passive = detbox_vol.placeVolume(passive_layer_vol,Transform3D(rot_layers,Position(x_passive_layer.x(), 0e0, z_layer)));
    z_layer += x_passive_layer.z()/2.;
    std::cout << "Zlayer passive " << z_layer << std::endl;
    pv_passive.addPhysVolID("passivelayer", iz*2+1);
    z_layer += extrazgap;
    if(iz==splitlayer){
    	z_layer += x_split.z()+x_widebar.z()+x_passive_layer.z();
    	PlacedVolume pv_split = detbox_vol.placeVolume(split_vol,Transform3D(rot,Position(x_split.x(), 0e0, z_layer)));
    	pv_split.addPhysVolID("split_layer", 9000);
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
  printout(INFO, "SplitCal", "%s: Detector construction finished.", nam.c_str());
  return sdet;
}

DECLARE_DETELEMENT(DD4hep_SplitCal,create_detector)
